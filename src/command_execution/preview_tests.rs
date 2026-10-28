//! Inline synthetic service fixtures reuse tutorial 06.07's invented CDS.
//! Channels and receipt facts establish stages; no model, native UI or private data.

use super::*;
use crate::engine::{Engine, Operation, Workflow, lifecycle_slice3_tests as fixture};
use std::{
    sync::mpsc,
    time::{Duration, Instant},
};

fn command(op: &Operation, named: bool) -> String {
    if named {
        let Operation::PlanDnaSequenceDesign {
            request,
            path: None,
        } = op
        else {
            panic!("no-path plan required");
        };
        format!(
            "sequence-design plan {}",
            crate::engine_shell::quote_shell_arg(&serde_json::to_string(request).unwrap())
        )
    } else {
        format!(
            "op {}",
            crate::engine_shell::quote_shell_arg(&serde_json::to_string(op).unwrap())
        )
    }
}

fn wait(service: &CommandExecutionService, id: u64) -> Result<ShellRunResult, String> {
    let deadline = Instant::now() + Duration::from_secs(20);
    loop {
        if let Some(result) = service.take_result(id) {
            return result;
        }
        assert!(
            Instant::now() < deadline,
            "preview worker timed out: {:?}",
            service.status(id)
        );
        std::thread::yield_now();
    }
}

fn held(
    service: &CommandExecutionService,
    engine: Arc<RwLock<GentleEngine>>,
) -> (u64, mpsc::Sender<()>) {
    let (ready_tx, ready_rx) = mpsc::channel();
    let (go_tx, go_rx) = mpsc::channel();
    let mut first = true;
    let options = ShellExecutionOptions {
        progress_callback: Some(Arc::new(Mutex::new(Box::new(move |event| {
            if first && matches!(event, OperationProgress::SequenceDesign { .. }) {
                first = false;
                ready_tx.send(()).unwrap();
                go_rx.recv_timeout(Duration::from_secs(20)).unwrap();
            }
            true
        })))),
        ..Default::default()
    };
    let text = command(&fixture::plan(None, true), false);
    let id = service
        .submit_parsed(engine, parse_shell_line(&text).unwrap(), text, options)
        .unwrap();
    ready_rx.recv_timeout(Duration::from_secs(20)).unwrap();
    (id, go_tx)
}

#[test]
fn lifecycle_slice5_preview_matches_direct_without_live_journal_or_history() {
    for loaded in [false, true] {
        for named in [false, true] {
            for parsed in [false, true] {
                let service = CommandExecutionService::default();
                let engine = Arc::new(RwLock::new(fixture::live_with_redo()));
                let baseline = engine.read().unwrap().clone(); // Test oracle only.
                let op = fixture::plan(None, loaded);
                let mut direct = baseline.clone();
                let expected = direct
                    .apply(op.clone())
                    .unwrap()
                    .dna_sequence_design
                    .unwrap();
                let text = command(&op, named);
                let id = if parsed {
                    service.submit_parsed(
                        engine.clone(),
                        parse_shell_line(&text).unwrap(),
                        text,
                        Default::default(),
                    )
                } else {
                    service.submit(engine.clone(), text, Default::default())
                }
                .unwrap();
                let run = wait(&service, id).unwrap();
                assert!(!run.state_changed);
                assert_eq!(
                    run.output["result"]["dna_sequence_design"],
                    serde_json::to_value(&expected).unwrap()
                );
                let receipt = service.status(id).unwrap();
                assert_eq!(receipt.publication, "preview");
                assert!(!receipt.requires_project_commit);
                assert_eq!(receipt.computed_operations, Some(1));
                assert_eq!(receipt.committed_operations, 0);
                assert_eq!(receipt.result_instance, Some(baseline.instance_id()));
                assert_eq!(
                    receipt.result_structural_revision,
                    Some(baseline.structural_revision())
                );
                assert_eq!(
                    receipt.result_sha256,
                    Some(sha256_prefixed_bytes(
                        &serde_json::to_vec(&run.output).unwrap()
                    ))
                );
                assert!(!service.cancel(id));
                assert!(service.take_result(id).is_none());
                assert!(service.receipts().iter().any(|r| r.job_id == id));
                let live = engine.read().unwrap();
                assert_eq!(
                    serde_json::to_value(live.state()).unwrap(),
                    serde_json::to_value(baseline.state()).unwrap()
                );
                assert_eq!(
                    (
                        live.execution_revision(),
                        live.mutation_revision(),
                        live.structural_revision()
                    ),
                    (
                        baseline.execution_revision(),
                        baseline.mutation_revision(),
                        baseline.structural_revision()
                    )
                );
                assert_eq!(
                    serde_json::to_value(live.history_summary()).unwrap(),
                    serde_json::to_value(baseline.history_summary()).unwrap()
                );
            }
        }
    }
}

#[test]
fn lifecycle_slice5_held_mutation_survives_preview_publication() {
    let service = CommandExecutionService::default();
    let engine = Arc::new(RwLock::new(fixture::live_with_redo()));
    let (ready_tx, ready_rx) = mpsc::channel();
    let (go_tx, go_rx) = mpsc::channel();
    let mutation = service
        .submit_work(
            engine.clone(),
            "synthetic-mutation".into(),
            move |worker, _| {
                ready_tx.send(()).unwrap();
                go_rx.recv_timeout(Duration::from_secs(20)).unwrap();
                let result = worker.apply(fixture::create()).map_err(|e| e.to_string())?;
                Ok(ShellRunResult {
                    state_changed: true,
                    output: serde_json::to_value(result).unwrap(),
                })
            },
        )
        .unwrap();
    ready_rx.recv_timeout(Duration::from_secs(20)).unwrap();
    let preview = service
        .submit(
            engine.clone(),
            command(&fixture::plan(None, true), false),
            Default::default(),
        )
        .unwrap();
    wait(&service, preview).unwrap();
    go_tx.send(()).unwrap();
    assert!(wait(&service, mutation).unwrap().state_changed);
    assert_eq!(service.status(mutation).unwrap().committed_operations, 1);
    assert!(
        engine
            .read()
            .unwrap()
            .state()
            .sequences
            .contains_key("synthetic_temporary")
    );
}

#[test]
fn lifecycle_slice5_two_previews_do_not_invalidate_each_other() {
    let service = CommandExecutionService::default();
    let engine = Arc::new(RwLock::new(fixture::live_with_redo()));
    let journal = engine.read().unwrap().journal_len();
    let (a, release_a) = held(&service, engine.clone());
    let (b, release_b) = held(&service, engine.clone());
    release_a.send(()).unwrap();
    release_b.send(()).unwrap();
    assert_eq!(
        wait(&service, a).unwrap().output,
        wait(&service, b).unwrap().output
    );
    assert_eq!(engine.read().unwrap().journal_len(), journal);
}

#[test]
fn lifecycle_slice5_preview_tolerates_a_read_only_journal_append_after_snapshot() {
    let service = CommandExecutionService::default();
    let engine = Arc::new(RwLock::new(fixture::live_with_redo()));
    let (id, release) = held(&service, engine.clone());
    {
        let mut live = engine.write().unwrap();
        let structural = live.structural_revision();
        let mutation = live.mutation_revision();
        let journal = live.journal_len();
        live.apply(fixture::snv()).unwrap();
        assert_eq!(live.structural_revision(), structural);
        assert_eq!(live.mutation_revision(), mutation);
        assert_eq!(live.journal_len(), journal + 1);
    }
    let expected = engine.read().unwrap().clone(); // Post-S1 oracle; preview must not replace it.
    release.send(()).unwrap();
    wait(&service, id).unwrap();
    let receipt = service.status(id).unwrap();
    assert_eq!(receipt.publication, "preview");
    assert_eq!(receipt.committed_operations, 0);
    let live = engine.read().unwrap();
    assert_eq!(
        serde_json::to_value(&*live).unwrap(),
        serde_json::to_value(&expected).unwrap()
    );
    assert_eq!(live.execution_revision(), expected.execution_revision());
    assert_eq!(live.undo_available(), expected.undo_available());
    assert_eq!(live.redo_available(), expected.redo_available());
}

#[test]
fn lifecycle_slice5_published_preview_apply_revalidates_only_its_actual_constraints() {
    use crate::dna_sequence::DNAsequence;
    use gentle_protocol::sequence_design::DnaSequenceDesignReport;

    for change in ["unrelated", "source_dna", "source_features", "output_id"] {
        let service = CommandExecutionService::default();
        let engine = Arc::new(RwLock::new(GentleEngine::from_state(fixture::state())));
        let id = service
            .submit(
                engine.clone(),
                command(&fixture::plan(None, true), false),
                Default::default(),
            )
            .unwrap();
        let run = wait(&service, id).unwrap();
        let report: DnaSequenceDesignReport =
            serde_json::from_value(run.output["result"]["dna_sequence_design"].clone()).unwrap();
        let approval = report.approval_digest.clone().unwrap();
        let output_id = report.request.output_seq_id.clone();
        let published_revision = service
            .status(id)
            .unwrap()
            .result_structural_revision
            .unwrap();
        let mut live = engine.write().unwrap();
        match change {
            "unrelated" => {
                live.apply(fixture::create()).unwrap();
            }
            "source_dna" => {
                // Synthetic synonymous change, still MFK: this is different source DNA.
                live.state_mut().sequences.insert(
                    "synthetic_source".into(),
                    DNAsequence::from_sequence("ATGTTCAAATAA").unwrap(),
                );
            }
            "source_features" => {
                live.state_mut()
                    .sequences
                    .get_mut("synthetic_source")
                    .unwrap()
                    .features_mut()
                    .push(gb_io::seq::Feature {
                        kind: "misc_feature".into(),
                        location: gb_io::seq::Location::simple_range(0, 3),
                        qualifiers: vec![(
                            "note".into(),
                            Some("Synthetic annotation added after preview".into()),
                        )],
                    });
            }
            "output_id" => {
                live.apply(Operation::CreateSequenceFromText {
                    sequence_text: "ATGC".into(),
                    output_id: Some(output_id.clone()),
                    name: None,
                    circular: false,
                })
                .unwrap();
            }
            _ => unreachable!(),
        }
        assert_ne!(live.structural_revision(), published_revision);
        let source_before_apply = live.state().sequences["synthetic_source"].get_forward_string();
        let journal_before_apply = live.journal_len();
        // Submit the original, unmodified proposal/digest explicitly. No automatic Apply.
        let applied = live.apply(Operation::ApplyDnaSequenceDesign {
            proposal: Box::new(report),
            approval_digest: approval,
        });
        if change == "unrelated" {
            let applied = applied.unwrap();
            let receipt = applied.dna_sequence_design_receipt.unwrap();
            assert!(receipt.output_constraints_verified);
            assert!(!receipt.search_claims_verified);
            assert_eq!(
                live.state().sequences[&output_id].get_forward_string(),
                "ATGTTTAAGTAA"
            );
            assert_eq!(live.journal_len(), journal_before_apply + 1);
        } else {
            let error = applied.unwrap_err();
            assert!(
                error.message.contains(if change == "output_id" {
                    "never overwrites"
                } else {
                    "prepare and review a new preview"
                }),
                "{error}"
            );
            if change == "output_id" {
                assert_eq!(
                    live.state().sequences[&output_id].get_forward_string(),
                    "ATGC"
                );
            } else {
                assert!(!live.state().sequences.contains_key(&output_id));
            }
            assert_eq!(live.journal_len(), journal_before_apply);
        }
        assert_eq!(
            live.state().sequences["synthetic_source"].get_forward_string(),
            source_before_apply
        );
    }
}

#[test]
fn lifecycle_slice5_rejects_structural_edit_and_replacement_before_publication() {
    for replace in [false, true] {
        let service = CommandExecutionService::default();
        let engine = Arc::new(RwLock::new(fixture::live_with_redo()));
        let (id, release) = held(&service, engine.clone());
        {
            let mut live = engine.write().unwrap();
            if replace {
                *live = GentleEngine::from_state(fixture::state());
            } else {
                live.apply(fixture::create()).unwrap();
            }
        }
        let expected = serde_json::to_value(engine.read().unwrap().state()).unwrap();
        release.send(()).unwrap();
        assert!(wait(&service, id).unwrap_err().contains("stale"));
        let receipt = service.status(id).unwrap();
        assert!(receipt.result_instance.is_none());
        assert!(receipt.result_sha256.is_none());
        assert_eq!(receipt.committed_operations, 0);
        assert_eq!(
            serde_json::to_value(engine.read().unwrap().state()).unwrap(),
            expected
        );
    }
}

#[test]
fn lifecycle_slice5_keeps_live_display_metadata_and_save_as() {
    let service = CommandExecutionService::default();
    let engine = Arc::new(RwLock::new(fixture::live_with_redo()));
    let td = tempfile::tempdir().unwrap();
    let path = td.path().join("new-host.json");
    let (id, release) = held(&service, engine.clone());
    {
        let mut live = engine.write().unwrap();
        live.display_state_mut().linear_view_start_bp = 3;
        live.auxiliary_metadata_mut()
            .insert("synthetic-live-edit".into(), serde_json::json!(true));
        live.set_host_project_path(Some(&path)).unwrap();
    }
    let expected = engine.read().unwrap().clone();
    release.send(()).unwrap();
    wait(&service, id).unwrap();
    let live = engine.read().unwrap();
    assert_eq!(
        serde_json::to_value(live.state()).unwrap(),
        serde_json::to_value(expected.state()).unwrap()
    );
    assert_eq!(live.mutation_revision(), expected.mutation_revision());
    assert_eq!(live.redo_available(), expected.redo_available());
    assert_eq!(live.host_project_path(), Some(path.as_path()));
}

#[test]
fn lifecycle_slice5_cancellation_during_search_and_after_computation_refuses_result() {
    for after_computation in [false, true] {
        let service = CommandExecutionService::default();
        let engine = Arc::new(RwLock::new(fixture::live_with_redo()));
        let journal = engine.read().unwrap().journal_len();
        let (id, release) = held(&service, engine.clone());
        if after_computation {
            let lock = engine.write().unwrap();
            release.send(()).unwrap();
            let deadline = Instant::now() + Duration::from_secs(20);
            // This count is set after execute_command returns, before publication.
            while service.status(id).unwrap().computed_operations != Some(1) {
                assert!(Instant::now() < deadline, "computation did not finish");
                std::thread::yield_now();
            }
            assert!(service.cancel(id));
            drop(lock);
        } else {
            assert!(service.cancel(id));
            release.send(()).unwrap();
        }
        assert!(
            wait(&service, id)
                .unwrap_err()
                .contains("cancelled before preview publication")
        );
        let receipt = service.status(id).unwrap();
        assert_eq!(receipt.state, RuntimeStatusFrameState::Cancelled);
        assert!(receipt.result_sha256.is_none());
        assert_eq!(receipt.committed_operations, 0);
        assert_eq!(engine.read().unwrap().journal_len(), journal);
    }
}

#[test]
fn lifecycle_slice5_existing_worker_file_input_and_noncoverage_controls() {
    let td = tempfile::tempdir().unwrap();
    let source = td.path().join("synthetic's operation.json");
    let op = fixture::plan(None, true);
    std::fs::write(
        &source,
        serde_json::to_string_pretty(&op)
            .unwrap()
            .replace('\n', "\r\n"),
    )
    .unwrap();
    let service = CommandExecutionService::default();
    let engine = Arc::new(RwLock::new(fixture::live_with_redo()));
    let text = format!(
        "op {}",
        crate::engine_shell::quote_shell_arg(&format!("@{}", source.display()))
    );
    let id = service
        .submit_parsed(
            engine.clone(),
            parse_shell_line(&text).unwrap(),
            text,
            Default::default(),
        )
        .unwrap();
    wait(&service, id).unwrap();
    assert_eq!(service.status(id).unwrap().publication, "preview");
    let report = td.path().join("report.json");
    let controls = [
        command(
            &fixture::plan(Some(report.display().to_string()), true),
            false,
        ),
        command(&fixture::snv(), false),
        format!(
            "workflow {}",
            crate::engine_shell::quote_shell_arg(
                &serde_json::to_string(&Workflow {
                    run_id: "synthetic-control".into(),
                    ops: vec![op.clone()]
                })
                .unwrap()
            )
        ),
    ];
    for text in controls {
        let id = service
            .submit(engine.clone(), text, Default::default())
            .unwrap();
        wait(&service, id).unwrap();
        assert_eq!(service.status(id).unwrap().publication, "project");
        assert_eq!(service.status(id).unwrap().committed_operations, 1);
    }
    let id = service
        .submit_work(
            engine,
            "trusted-preview-control".into(),
            move |worker, _| {
                let result = worker.apply(op).map_err(|e| e.to_string())?;
                Ok(ShellRunResult {
                    state_changed: false,
                    output: serde_json::to_value(result).unwrap(),
                })
            },
        )
        .unwrap();
    wait(&service, id).unwrap();
    assert_eq!(service.status(id).unwrap().publication, "project");
    assert_eq!(service.status(id).unwrap().committed_operations, 1);
}

#[test]
fn lifecycle_slice5_unresolved_scientific_outcomes_are_successfully_published_not_feasible() {
    for frozen in [false, true] {
        let mut op = fixture::plan(None, true);
        if let Operation::PlanDnaSequenceDesign { request, .. } = &mut op {
            if frozen {
                request.protected_intervals = vec![request.cds.clone()];
            } else {
                request.max_evaluations = 1;
            }
        }
        let service = CommandExecutionService::default();
        let engine = Arc::new(RwLock::new(fixture::live_with_redo()));
        let id = service
            .submit(engine, command(&op, false), Default::default())
            .unwrap();
        let output = wait(&service, id).unwrap().output;
        let report = &output["result"]["dna_sequence_design"];
        assert_eq!(
            report["status"],
            if frozen {
                "proven_infeasible"
            } else {
                "search_exhausted"
            }
        );
        assert!(report["approval_digest"].is_null());
        assert!(report["output_sequence"].is_null());
        assert_eq!(
            service.status(id).unwrap().state,
            RuntimeStatusFrameState::Completed
        );
        assert_eq!(service.status(id).unwrap().publication, "preview");
    }
}

#[test]
fn lifecycle_slice5_preview_has_no_compact_fallback_for_oversized_output() {
    assert!(
        bounded_output(&serde_json::json!({"synthetic": "x".repeat(MAX_RESULT_BYTES)})).is_err()
    );
    // Publication infrastructure, not a reachable large solver fixture.
    let service = CommandExecutionService::default();
    let engine = Arc::new(RwLock::new(fixture::live_with_redo()));
    let id = service.submit_work_with_cancel(engine.clone(), "synthetic-output-boundary".into(),
        Arc::new(AtomicBool::new(false)), Lane::Execution, true, |worker, _| {
            let result = worker.apply(fixture::plan(None, true)).map_err(|e| e.to_string())?;
            Ok(WorkResult { run: ShellRunResult { state_changed: false, output: serde_json::json!({"result":result,"synthetic_padding":"x".repeat(MAX_RESULT_BYTES)}) }, publication: Publication::Preview })
        }).unwrap();
    assert!(wait(&service, id).unwrap_err().contains("16 MiB"));
    assert_eq!(service.status(id).unwrap().publication, "pending");
    assert!(service.status(id).unwrap().result_sha256.is_none());
}

#[test]
fn lifecycle_slice5_cancel_before_snapshot_and_invalid_request_never_publish() {
    let service = CommandExecutionService::default();
    let engine = Arc::new(RwLock::new(fixture::live_with_redo()));
    let journal = engine.read().unwrap().journal_len();
    let id = service
        .submit_work_with_cancel(
            engine.clone(),
            "synthetic-early-cancel".into(),
            Arc::new(AtomicBool::new(true)),
            Lane::Execution,
            true,
            |_, _| panic!("work must not start after pre-snapshot cancellation"),
        )
        .unwrap();
    assert!(
        wait(&service, id)
            .unwrap_err()
            .contains("cancelled before snapshot")
    );
    assert_eq!(
        service.status(id).unwrap().state,
        RuntimeStatusFrameState::Cancelled
    );
    assert_eq!(service.status(id).unwrap().computed_operations, None);
    let mut op = fixture::plan(None, true);
    if let Operation::PlanDnaSequenceDesign { request, .. } = &mut op {
        request.output_seq_id = "synthetic_source".into(); // Existing ID is refused at preview, too.
    }
    let id = service
        .submit(engine.clone(), command(&op, false), Default::default())
        .unwrap();
    assert!(wait(&service, id).unwrap_err().contains("never overwrites"));
    assert!(service.status(id).unwrap().result_instance.is_none());
    assert!(service.status(id).unwrap().result_sha256.is_none());
    let id = service
        .submit(
            engine.clone(),
            command(&fixture::plan(Some(String::new()), true), false),
            Default::default(),
        )
        .unwrap();
    assert!(wait(&service, id).is_err());
    let receipt = service.status(id).unwrap();
    assert!(
        receipt.requires_project_commit,
        "an invalid path is still on the legacy project path"
    );
    assert_ne!(receipt.publication, "preview");
    assert_eq!(engine.read().unwrap().journal_len(), journal);
}
