//! Host-owned, bounded command admission above the synchronous shell executor.
//!
//! Workers own detached engines. Observation and cancellation never acquire the
//! project lock. Only a successful, still-current result is committed; failed
//! or cancelled work can leave external resource/file effects but no live engine
//! delta. Receipts survive consumption of the result until bounded eviction.

use crate::{
    digest_utils::sha256_prefixed_bytes,
    engine::{GentleEngine, OperationProgress},
    engine_shell::{
        ShellCommand, ShellExecutionOptions, ShellRunResult, execute_shell_command_with_options,
        parse_shell_line,
    },
    runtime_status::{
        RuntimeParentScope, RuntimeStatusFrameKind, RuntimeStatusFrameState,
        runtime_status_registry,
    },
};
use serde::Serialize;
use std::sync::atomic::{AtomicBool, Ordering};
use std::{
    collections::BTreeMap,
    sync::{Arc, Mutex, RwLock},
    thread,
};

const MAX_ACTIVE: usize = 4;
const MAX_RETAINED: usize = 32;
const MAX_RESULT_BYTES: usize = 16 * 1024 * 1024;

/// Cached receipt. Counts describe operations, never a runtime estimate or a
/// biological acceptance verdict. Approval remains the shell's responsibility.
#[derive(Debug, Clone, Serialize)]
pub struct CommandReceipt {
    pub schema: &'static str,
    pub job_id: u64,
    pub owner_instance: u64,
    pub result_instance: Option<u64>,
    /// Revision at commit, so UI consumers can detect later edits without
    /// confusing the command's own structural changes with stale results.
    pub result_structural_revision: Option<u64>,
    pub runtime_frame_id: Option<String>,
    pub command_sha256: String,
    pub state: RuntimeStatusFrameState,
    pub phase: String,
    pub cancel_requested: bool,
    /// Absent when failure/rollback leaves no reliable completed journal count.
    pub computed_operations: Option<usize>,
    pub committed_operations: usize,
    pub completed_steps: Option<usize>,
    pub total_steps: Option<usize>,
    pub result_sha256: Option<String>,
    pub error: Option<String>,
}

struct Job {
    cancel: Arc<AtomicBool>,
    published: bool,
    receipt: CommandReceipt,
    result: Option<Result<ShellRunResult, String>>,
}

#[derive(Default)]
struct Store {
    next: u64,
    jobs: BTreeMap<u64, Arc<Mutex<Job>>>,
}

/// One local host's admission/control service; no tool probing on submission.
#[derive(Clone, Default)]
pub struct CommandExecutionService(Arc<Mutex<Store>>);

fn terminal(state: RuntimeStatusFrameState) -> bool {
    matches!(
        state,
        RuntimeStatusFrameState::Completed
            | RuntimeStatusFrameState::Failed
            | RuntimeStatusFrameState::Cancelled
    )
}

fn fork_admitted_command(
    live: &GentleEngine,
    owner: u64,
    structural_revision: u64,
    journal_len: usize,
) -> Result<crate::engine::DetachedEngineExecution, String> {
    if live.instance_id() != owner {
        return Err("Command belongs to a replaced project instance".into());
    }
    if live.structural_revision() != structural_revision || live.journal_len() != journal_len {
        return Err("Command became stale after admission; execution did not start".into());
    }
    Ok(live.fork_detached_execution())
}

impl CommandExecutionService {
    /// First audited migration group. Other shell routes retain their existing
    /// host policy until their completion/file/job semantics have been audited.
    pub fn manages(command: &ShellCommand) -> bool {
        matches!(
            command,
            ShellCommand::ReferencePrepare { .. }
                | ShellCommand::Workflow { .. }
                | ShellCommand::Op { .. }
                | ShellCommand::MacrosRun { .. }
                | ShellCommand::MacrosTemplateRun { .. }
                | ShellCommand::CandidatesMacro { .. }
                | ShellCommand::CandidatesTemplateRun { .. }
                | ShellCommand::PrimersExecuteGeneIsoformAssayStudyWorkflow { .. }
        )
    }
    /// Admit the exact parsed command. Capacity or engine-lock contention is
    /// reported immediately rather than turning the caller into a waiter.
    pub fn submit(
        &self,
        engine: Arc<RwLock<GentleEngine>>,
        command_text: String,
        options: ShellExecutionOptions,
    ) -> Result<u64, String> {
        if command_text.len() > 1024 * 1024 {
            return Err("Interactive command exceeds 1 MiB; use a bound @file request".into());
        }
        let cancel = Arc::new(AtomicBool::new(false));
        let caller_cancel = cancel.clone();
        self.submit_work_with_cancel(
            engine,
            sha256_prefixed_bytes(command_text.as_bytes()),
            cancel,
            move |engine, callback| {
                let command = parse_shell_line(&command_text)?;
                let caller_callback = options.progress_callback.clone();
                let options = ShellExecutionOptions {
                    progress_callback: Some(Arc::new(Mutex::new(Box::new(move |event| {
                        if !callback
                            .lock()
                            .map(|mut f| f(event.clone()))
                            .unwrap_or(false)
                        {
                            return false;
                        }
                        let keep_going = caller_callback
                            .as_ref()
                            .map(|cb| cb.lock().map(|mut f| f(event)).unwrap_or(false))
                            .unwrap_or(true);
                        if !keep_going {
                            caller_cancel.store(true, Ordering::Release);
                        }
                        keep_going
                    })))),
                    ..options
                };
                execute_shell_command_with_options(engine, &command, &options)
            },
        )
    }

    /// Admit trusted in-process work with the same snapshot, cancellation and
    /// commit rules as shell commands. Callers own input/approval validation.
    pub(crate) fn submit_work<F>(
        &self,
        engine: Arc<RwLock<GentleEngine>>,
        command_sha256: String,
        work: F,
    ) -> Result<u64, String>
    where
        F: FnOnce(
                &mut GentleEngine,
                crate::engine_shell::ShellProgressCallback,
            ) -> Result<ShellRunResult, String>
            + Send
            + 'static,
    {
        self.submit_work_with_cancel(
            engine,
            command_sha256,
            Arc::new(AtomicBool::new(false)),
            work,
        )
    }

    fn submit_work_with_cancel<F>(
        &self,
        engine: Arc<RwLock<GentleEngine>>,
        command_sha256: String,
        cancel: Arc<AtomicBool>,
        work: F,
    ) -> Result<u64, String>
    where
        F: FnOnce(
                &mut GentleEngine,
                crate::engine_shell::ShellProgressCallback,
            ) -> Result<ShellRunResult, String>
            + Send
            + 'static,
    {
        let (owner, admitted_revision, admitted_journal_len) = {
            let live = engine.try_read().map_err(|_| {
                "Project is busy; command was not admitted. Retry shortly.".to_string()
            })?;
            (
                live.instance_id(),
                live.structural_revision(),
                live.journal_len(),
            )
        };
        let (id, job) = {
            let mut store = self.0.lock().map_err(|_| "Command registry unavailable")?;
            let active = store
                .jobs
                .values()
                .filter(|job| !terminal(job.lock().unwrap().receipt.state))
                .count();
            if active >= MAX_ACTIVE {
                return Err("Command capacity reached; command was not admitted".into());
            }
            while store.jobs.len() >= MAX_RETAINED {
                // Unconsumed outputs are never evicted on the submitting thread.
                let removable = store
                    .jobs
                    .iter()
                    .find(|(_, job)| {
                        let job = job.lock().unwrap();
                        terminal(job.receipt.state) && job.result.is_none()
                    })
                    .map(|(id, _)| *id);
                let Some(id) = removable else {
                    return Err(
                        "Consume completed command results before admitting more work".into(),
                    );
                };
                store.jobs.remove(&id);
            }
            store.next += 1;
            let id = store.next;
            let job = Arc::new(Mutex::new(Job {
                cancel: cancel.clone(),
                published: false,
                receipt: CommandReceipt {
                    schema: "gentle.command_execution.v1",
                    job_id: id,
                    owner_instance: owner,
                    result_instance: None,
                    result_structural_revision: None,
                    runtime_frame_id: None,
                    command_sha256,
                    state: RuntimeStatusFrameState::Waiting,
                    phase: "admitted".into(),
                    cancel_requested: false,
                    computed_operations: None,
                    committed_operations: 0,
                    completed_steps: None,
                    total_steps: None,
                    result_sha256: None,
                    error: None,
                },
                result: None,
            }));
            store.jobs.insert(id, job.clone());
            (id, job)
        };
        let parent = RuntimeParentScope::capture();
        let worker_job = job.clone();
        #[cfg(test)]
        let tool_overrides = crate::tool_overrides::scoped_tool_overrides_snapshot();
        let spawn = thread::Builder::new().name(format!("gentle-command-{id}")).stack_size(16 * 1024 * 1024).spawn(move || {
            #[cfg(test)]
            let _tool_overrides = crate::tool_overrides::ScopedToolOverridesSnapshotGuard::install(tool_overrides);
            let _parent = RuntimeParentScope::enter(parent);
            let frame = Arc::new(runtime_status_registry().push(RuntimeStatusFrameKind::BackgroundJob, format!("Command {id}")));
            worker_job.lock().unwrap().receipt.runtime_frame_id = Some(frame.frame_id().to_string());
            let run = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| -> Result<ShellRunResult, String> {
                if cancel.load(Ordering::Acquire) { return Err("Command cancelled before snapshot".into()); }
                frame.update_phase("snapshot");
                let mut detached = {
                    let live = engine.read().map_err(|_| "Project lock unavailable")?;
                    fork_admitted_command(&live, owner, admitted_revision, admitted_journal_len)?
                };
                let baseline = detached.engine().journal_len();
                if cancel.load(Ordering::Acquire) { return Err("Command cancelled after snapshot; execution did not start".into()); }
                {
                    let mut job = worker_job.lock().unwrap();
                    job.receipt.state = RuntimeStatusFrameState::Running;
                    job.receipt.phase = "executing".into();
                }
                frame.update_phase("executing");
                let callback_job = worker_job.clone();
                let callback_frame = frame.clone();
                let callback: crate::engine_shell::ShellProgressCallback = Arc::new(Mutex::new(Box::new(move |event| {
                    let mut job = callback_job.lock().unwrap();
                    if job.cancel.load(Ordering::Acquire) { return false; }
                    match event {
                        OperationProgress::Workflow { completed, total } => {
                            job.receipt.completed_steps = Some(completed);
                            job.receipt.total_steps = Some(total);
                            callback_frame.update_detail(format!("{completed}/{total} workflow steps computed; awaiting commit"));
                        }
                        OperationProgress::GenomePrepare(p) => {
                            job.receipt.phase = p.phase.clone();
                            callback_frame.update_from_progress(p.phase, p.item, p.bytes_done, p.bytes_total, p.percent);
                        }
                        OperationProgress::PrimerDesign(p) => {
                            job.receipt.phase = format!("{}: {}", p.design_kind, p.stage);
                            callback_frame.update_detail(p.detail);
                        }
                        _ => {}
                    }
                    true
                })));
                let result = work(detached.engine_mut(), callback);
                let retained_operations = detached.engine().journal_len().saturating_sub(baseline);
                worker_job.lock().unwrap().receipt.computed_operations = if result.is_err() && retained_operations == 0 { None } else { Some(retained_operations) };
                let result = result?;
                if detached.engine().instance_id() != owner {
                    return Err("Managed command replaced its project instance; use the explicit project load/import route instead".into());
                }
                let bytes = serde_json::to_vec(&result.output).map_err(|e| e.to_string())?;
                if bytes.len() > MAX_RESULT_BYTES { return Err("Command result exceeds the 16 MiB interactive retention limit; use a file-export operation".into()); }
                let digest = sha256_prefixed_bytes(&bytes);
                drop(bytes);
                frame.update_phase("committing");
                let old = {
                    let mut live = engine.write().map_err(|_| "Project lock unavailable")?;
                    let mut job = worker_job.lock().unwrap();
                    if cancel.load(Ordering::Acquire) { return Err("Command cancelled before commit; computed engine changes discarded".into()); }
                    let old = live.commit_detached_execution(&mut detached).map_err(|e| e.to_string())?;
                    // Cancellation and commit are serialized by the same small control lock.
                    job.published = true;
                    job.receipt.phase = "committed".into();
                    job.receipt.result_instance = Some(live.instance_id());
                    job.receipt.result_structural_revision = Some(live.structural_revision());
                    job.receipt.committed_operations = job.receipt.computed_operations.unwrap_or(0);
                    job.receipt.result_sha256 = Some(digest);
                    old
                };
                drop(old);
                Ok(result)
            })).unwrap_or_else(|_| Err("Command worker panicked; no success receipt was issued".into()));
            let mut job = worker_job.lock().unwrap();
            job.receipt.cancel_requested = cancel.load(Ordering::Acquire);
            if let Err(error) = &run {
                job.receipt.state = if job.receipt.cancel_requested { RuntimeStatusFrameState::Cancelled } else { RuntimeStatusFrameState::Failed };
                job.receipt.phase = if job.receipt.cancel_requested { "cancelled" } else { "failed" }.into();
                job.receipt.error = Some(error.clone());
                if job.receipt.cancel_requested { frame.cancel(error); } else { frame.fail(error); }
            } else {
                job.receipt.state = RuntimeStatusFrameState::Completed;
                job.receipt.phase = "completed".into();
                frame.update_state(RuntimeStatusFrameState::Completed);
            }
            job.result = Some(run);
        });
        if let Err(error) = spawn {
            self.0.lock().unwrap().jobs.remove(&id);
            return Err(format!("Could not start command worker: {error}"));
        }
        Ok(id)
    }

    /// Cached observation only: no engine lock, I/O, probing or dispatch.
    pub fn status(&self, id: u64) -> Option<CommandReceipt> {
        let job = self.0.lock().ok()?.jobs.get(&id)?.clone();
        let job = job.lock().ok()?;
        let mut receipt = job.receipt.clone();
        receipt.cancel_requested = job.cancel.load(Ordering::Acquire);
        if receipt.cancel_requested && !terminal(receipt.state) {
            receipt.phase = "cancelling".into();
        }
        Some(receipt)
    }

    /// Cooperative request, not a terminal cancellation or proof of child exit.
    pub fn cancel(&self, id: u64) -> bool {
        let Some(job) = self.0.lock().ok().and_then(|s| s.jobs.get(&id).cloned()) else {
            return false;
        };
        let mut job = job.lock().unwrap();
        if terminal(job.receipt.state) || job.published {
            return false;
        }
        job.receipt.cancel_requested = true;
        job.cancel.store(true, Ordering::Release);
        job.receipt.phase = "cancelling".into();
        let frame_id = job.receipt.runtime_frame_id.clone();
        drop(job);
        if let Some(frame_id) = frame_id {
            runtime_status_registry().request_frame_cancel(&frame_id);
        }
        true
    }

    /// Move, rather than clone, a completed output. The small receipt remains.
    pub fn take_result(&self, id: u64) -> Option<Result<ShellRunResult, String>> {
        let job = self.0.lock().ok()?.jobs.get(&id)?.clone();
        let result = job.lock().ok()?.result.take();
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::{Engine, Operation, Workflow};
    use std::{
        sync::mpsc,
        time::{Duration, Instant},
    };

    fn create(id: &str) -> Operation {
        Operation::CreateSequenceFromText {
            sequence_text: "ATGC".into(),
            output_id: Some(id.into()),
            name: None,
            circular: false,
        }
    }

    fn success() -> Result<ShellRunResult, String> {
        Ok(ShellRunResult {
            state_changed: false,
            output: serde_json::json!({"ok": true}),
        })
    }

    fn wait(service: &CommandExecutionService, id: u64) -> Result<ShellRunResult, String> {
        let deadline = Instant::now() + Duration::from_secs(10);
        loop {
            if let Some(result) = service.take_result(id) {
                return result;
            }
            assert!(
                Instant::now() < deadline,
                "command did not finish: {:?}",
                service.status(id)
            );
            thread::sleep(Duration::from_millis(2));
        }
    }

    #[test]
    fn admitted_baseline_rejects_structural_edits_before_snapshot_but_allows_metadata() {
        let mut live = GentleEngine::new();
        let owner = live.instance_id();
        let revision = live.structural_revision();
        let journal_len = live.journal_len();
        live.auxiliary_metadata_mut()
            .insert("live-display-context".into(), serde_json::json!(true));
        assert!(fork_admitted_command(&live, owner, revision, journal_len).is_ok());
        live.apply(create("later-edit")).unwrap();
        assert!(
            fork_admitted_command(&live, owner, revision, journal_len)
                .unwrap_err()
                .contains("stale after admission")
        );
    }

    #[test]
    fn held_command_leaves_engine_and_control_available_and_cancel_discards_work() {
        let service = CommandExecutionService::default();
        let engine = Arc::new(RwLock::new(GentleEngine::new()));
        let (entered_tx, entered_rx) = mpsc::channel();
        let (release_tx, release_rx) = mpsc::channel();
        let id = service
            .submit_work(engine.clone(), "synthetic".into(), move |e, _| {
                e.apply(create("discarded")).unwrap();
                entered_tx.send(()).unwrap();
                release_rx.recv_timeout(Duration::from_secs(5)).unwrap();
                success()
            })
            .unwrap();
        entered_rx.recv_timeout(Duration::from_secs(5)).unwrap();
        let live = engine
            .try_write()
            .expect("work must not hold live engine lock");
        assert_eq!(
            service.status(id).unwrap().state,
            RuntimeStatusFrameState::Running
        );
        assert!(
            service.cancel(id),
            "control must work even while another caller holds the engine lock"
        );
        let frame_id = service.status(id).unwrap().runtime_frame_id.unwrap();
        let snapshot =
            runtime_status_registry().snapshot(crate::runtime_status::RuntimeStatusTrigger::Shell);
        let frame = snapshot
            .frames
            .iter()
            .find(|frame| frame.frame_id == frame_id)
            .unwrap();
        assert!(frame.cancel_requested);
        assert_eq!(frame.phase.as_deref(), Some("cancelling"));
        assert!(live.state().sequences.is_empty());
        drop(live);
        release_tx.send(()).unwrap();
        assert!(wait(&service, id).is_err());
        let receipt = service.status(id).unwrap();
        assert_eq!(receipt.state, RuntimeStatusFrameState::Cancelled);
        assert_eq!(receipt.computed_operations, Some(1));
        assert_eq!(receipt.committed_operations, 0);
        assert_eq!(engine.read().unwrap().journal_len(), 0);
        assert!(!service.cancel(id));
    }

    #[test]
    fn replaced_project_rejects_held_command_and_retains_failure() {
        let service = CommandExecutionService::default();
        let engine = Arc::new(RwLock::new(GentleEngine::new()));
        let (tx, rx) = mpsc::channel();
        let (release_tx, release_rx) = mpsc::channel();
        let id = service
            .submit_work(engine.clone(), "synthetic".into(), move |e, _| {
                tx.send(()).unwrap();
                release_rx.recv_timeout(Duration::from_secs(5)).unwrap();
                e.apply(create("old-project")).unwrap();
                success()
            })
            .unwrap();
        rx.recv_timeout(Duration::from_secs(5)).unwrap();
        *engine.write().unwrap() = GentleEngine::new();
        release_tx.send(()).unwrap();
        assert!(wait(&service, id).unwrap_err().contains("project instance"));
        assert_eq!(service.status(id).unwrap().committed_operations, 0);
        assert!(engine.read().unwrap().state().sequences.is_empty());
    }

    #[test]
    fn shell_workflow_retains_exact_request_and_result_hashes_and_journals_once() {
        let service = CommandExecutionService::default();
        let engine = Arc::new(RwLock::new(GentleEngine::new()));
        let wf = Workflow {
            run_id: "synthetic-run".into(),
            ops: vec![create("one"), create("two")],
        };
        let text = format!("workflow '{}'", serde_json::to_string(&wf).unwrap());
        let id = service
            .submit(
                engine.clone(),
                text.clone(),
                ShellExecutionOptions::default(),
            )
            .unwrap();
        let result = wait(&service, id).unwrap();
        let receipt = service.status(id).unwrap();
        assert_eq!(
            receipt.command_sha256,
            sha256_prefixed_bytes(text.as_bytes())
        );
        assert_eq!(
            receipt.result_sha256,
            Some(sha256_prefixed_bytes(
                &serde_json::to_vec(&result.output).unwrap()
            ))
        );
        assert_eq!(receipt.committed_operations, 2);
        assert_eq!(receipt.completed_steps, Some(2));
        assert_eq!(engine.read().unwrap().journal_len(), 2);
        assert!(
            !service.cancel(id),
            "completion must not be relabelled cancelled"
        );
        assert!(service.take_result(id).is_none());
    }

    #[test]
    fn local_genome_prepare_streams_progress_and_journals_once_even_with_warm_cache() {
        let _makeblastdb = crate::tool_overrides::ScopedToolOverrideGuard::set(
            crate::genomes::MAKEBLASTDB_ENV_BIN,
            "__missing_test_makeblastdb__",
        );
        let _blastn = crate::tool_overrides::ScopedToolOverrideGuard::set(
            crate::genomes::BLASTN_ENV_BIN,
            "__missing_test_blastn__",
        );
        let _blastdbcmd = crate::tool_overrides::ScopedToolOverrideGuard::set(
            crate::genomes::BLASTDBCMD_ENV_BIN,
            "__missing_test_blastdbcmd__",
        );
        let td = tempfile::tempdir().unwrap();
        let fasta = td.path().join("synthetic.fa");
        let annotation = td.path().join("synthetic.gtf");
        let catalog = td.path().join("catalog.json");
        std::fs::write(&fasta, ">chr1\nACGTACGT\n").unwrap();
        std::fs::write(
            &annotation,
            "chr1\tsynthetic\tgene\t1\t8\t.\t+\t.\tgene_id \"G1\";\n",
        )
        .unwrap();
        std::fs::write(&catalog, serde_json::to_vec(&serde_json::json!({"Synthetic": {
            "sequence_local": fasta, "annotations_local": annotation, "cache_dir": td.path().join("cache")
        }})).unwrap()).unwrap();
        let service = CommandExecutionService::default();
        let engine = Arc::new(RwLock::new(GentleEngine::new()));
        for expected_journal in 1..=2 {
            let progress = Arc::new(std::sync::atomic::AtomicUsize::new(0));
            let observed = progress.clone();
            let options = ShellExecutionOptions {
                progress_callback: Some(Arc::new(Mutex::new(Box::new(move |event| {
                    if matches!(event, OperationProgress::GenomePrepare(_)) {
                        observed.fetch_add(1, Ordering::Relaxed);
                    }
                    true
                })))),
                ..ShellExecutionOptions::default()
            };
            let id = service
                .submit(
                    engine.clone(),
                    format!(
                        "genomes prepare Synthetic --catalog '{}'",
                        catalog.display()
                    ),
                    options,
                )
                .unwrap();
            let result = wait(&service, id).unwrap();
            assert!(result.output.get("binary_preflight").is_some());
            assert!(progress.load(Ordering::Relaxed) > 0);
            assert_eq!(engine.read().unwrap().journal_len(), expected_journal);
            assert_eq!(service.status(id).unwrap().committed_operations, 1);
        }
    }

    #[test]
    fn failure_discards_detached_prefix_and_retains_exact_error() {
        let service = CommandExecutionService::default();
        let engine = Arc::new(RwLock::new(GentleEngine::new()));
        let id = service
            .submit_work(engine.clone(), "synthetic".into(), |e, _| {
                e.apply(create("discarded")).unwrap();
                Err("synthetic failure receipt".into())
            })
            .unwrap();
        assert_eq!(wait(&service, id).unwrap_err(), "synthetic failure receipt");
        assert_eq!(service.status(id).unwrap().computed_operations, Some(1));
        assert_eq!(service.status(id).unwrap().committed_operations, 0);
        assert!(engine.read().unwrap().state().sequences.is_empty());

        let rolled_back = service
            .submit_work(engine.clone(), "rollback".into(), |e, callback| {
                e.with_rollback_on_error(|e| -> Result<ShellRunResult, String> {
                    e.apply(create("rolled-back")).unwrap();
                    assert!(callback.lock().unwrap()(OperationProgress::Workflow {
                        completed: 1,
                        total: 2
                    }));
                    Err("synthetic rollback".into())
                })
            })
            .unwrap();
        assert!(wait(&service, rolled_back).is_err());
        let receipt = service.status(rolled_back).unwrap();
        assert_eq!(
            receipt.computed_operations, None,
            "rolled-back journal is not proof of zero work"
        );
        assert_eq!(receipt.completed_steps, Some(1));
        assert_eq!(receipt.committed_operations, 0);
    }

    #[test]
    fn command_cannot_merge_a_replacement_engine_as_an_ordinary_delta() {
        let service = CommandExecutionService::default();
        let engine = Arc::new(RwLock::new(GentleEngine::new()));
        let owner = engine.read().unwrap().instance_id();
        let id = service
            .submit_work(engine.clone(), "replace".into(), |e, _| {
                *e = GentleEngine::new();
                success()
            })
            .unwrap();
        assert!(
            wait(&service, id)
                .unwrap_err()
                .contains("replaced its project instance")
        );
        assert_eq!(engine.read().unwrap().instance_id(), owner);
    }

    #[test]
    fn caller_progress_cancellation_is_not_reported_as_failure_or_partial_success() {
        let service = CommandExecutionService::default();
        let engine = Arc::new(RwLock::new(GentleEngine::new()));
        let wf = Workflow {
            run_id: "cancel-test".into(),
            ops: vec![create("one"), create("two")],
        };
        let options = ShellExecutionOptions {
            progress_callback: Some(Arc::new(Mutex::new(Box::new(|event| {
                !matches!(event, OperationProgress::Workflow { completed: 1, .. })
            })))),
            ..ShellExecutionOptions::default()
        };
        let id = service
            .submit(
                engine.clone(),
                format!("workflow '{}'", serde_json::to_string(&wf).unwrap()),
                options,
            )
            .unwrap();
        assert!(wait(&service, id).is_err());
        let receipt = service.status(id).unwrap();
        assert_eq!(receipt.state, RuntimeStatusFrameState::Cancelled);
        assert_eq!(receipt.completed_steps, Some(1));
        assert_eq!(receipt.computed_operations, Some(1));
        assert_eq!(receipt.committed_operations, 0);
        assert!(engine.read().unwrap().state().sequences.is_empty());
    }

    #[test]
    fn admission_under_write_lock_is_nonblocking_and_does_not_start_work() {
        let service = CommandExecutionService::default();
        let engine = Arc::new(RwLock::new(GentleEngine::new()));
        let _guard = engine.write().unwrap();
        assert!(
            service
                .submit_work(engine.clone(), "synthetic".into(), |_, _| panic!(
                    "must not run"
                ))
                .unwrap_err()
                .contains("not admitted")
        );
    }

    #[test]
    fn capacity_is_bounded_and_terminal_receipts_survive_result_consumption() {
        let service = CommandExecutionService::default();
        let engine = Arc::new(RwLock::new(GentleEngine::new()));
        let mut releases = Vec::new();
        for _ in 0..MAX_ACTIVE {
            let (tx, rx) = mpsc::channel();
            let id = service
                .submit_work(engine.clone(), "synthetic".into(), move |_, _| {
                    rx.recv_timeout(Duration::from_secs(5)).unwrap();
                    success()
                })
                .unwrap();
            releases.push((id, tx));
        }
        assert!(
            service
                .submit_work(engine, "synthetic".into(), |_, _| success())
                .is_err()
        );
        for (id, tx) in releases {
            tx.send(()).unwrap();
            wait(&service, id).unwrap();
            assert!(service.status(id).is_some());
        }
    }

    #[test]
    fn shell_parent_frame_is_inherited_across_expanded_stack_worker() {
        let service = CommandExecutionService::default();
        let engine = Arc::new(RwLock::new(GentleEngine::new()));
        let (tx, rx) = mpsc::channel();
        let wf = Workflow {
            run_id: "parent-test".into(),
            ops: vec![create("one")],
        };
        let callback = Arc::new(Mutex::new(Box::new(move |_| {
            let snapshot = runtime_status_registry()
                .snapshot(crate::runtime_status::RuntimeStatusTrigger::Shell);
            tx.send(snapshot).unwrap();
            true
        })
            as Box<dyn FnMut(OperationProgress) -> bool + Send>));
        let id = service
            .submit(
                engine,
                format!("workflow '{}'", serde_json::to_string(&wf).unwrap()),
                ShellExecutionOptions {
                    progress_callback: Some(callback),
                    ..ShellExecutionOptions::default()
                },
            )
            .unwrap();
        let snapshot = rx.recv_timeout(Duration::from_secs(5)).unwrap();
        assert!(snapshot.frames.iter().any(|shell| {
            shell.label == "shared shell command"
                && snapshot.frames.iter().any(|parent| {
                    shell.parent_frame_id.as_ref() == Some(&parent.frame_id)
                        && parent.label == format!("Command {id}")
                })
        }));
        wait(&service, id).unwrap();
    }
}
