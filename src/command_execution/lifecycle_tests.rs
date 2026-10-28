//! Synthetic managed-publication controls; no model, native GUI or private sequence.

use super::*;
use crate::engine::{Engine, Workflow, lifecycle_slice3_tests as fixture};
use std::{
    sync::mpsc,
    time::{Duration, Instant},
};

fn wait(service: &CommandExecutionService, id: u64) -> Result<ShellRunResult, String> {
    let deadline = Instant::now() + Duration::from_secs(10);
    loop {
        if let Some(result) = service.take_result(id) {
            return result;
        }
        assert!(
            Instant::now() < deadline,
            "managed lifecycle worker timed out"
        );
        std::thread::sleep(Duration::from_millis(2));
    }
}

#[test]
fn observed_transcripts_managed_publication_keeps_live_project_and_journal_unchanged() {
    let td = tempfile::tempdir().unwrap();
    let (_, request) = GentleEngine::observed_transcripts_test_fixture(td.path());
    let service = CommandExecutionService::default();
    let engine = Arc::new(RwLock::new(fixture::live_with_redo()));
    let before = engine.read().unwrap();
    let state = serde_json::to_value(before.state()).unwrap();
    let journal = before.journal_len();
    let history = before.history_summary();
    drop(before);
    let command = format!(
        "rna-reads transcript-structures {}",
        crate::engine_shell::quote_shell_arg(&serde_json::to_string(&request).unwrap())
    );
    let id = service
        .submit(engine.clone(), command, Default::default())
        .unwrap();
    let result = wait(&service, id).unwrap();
    assert!(!result.state_changed);
    assert_eq!(
        result.output["result"]["observed_transcript_structures"]["status"],
        "completed"
    );
    let receipt = service.status(id).unwrap();
    assert_eq!(receipt.publication, "artifact_only");
    assert!(!receipt.requires_project_commit);
    assert_eq!(receipt.committed_operations, 0);
    let report: crate::engine::ObservedTranscriptStructuresReport =
        serde_json::from_value(result.output["result"]["observed_transcript_structures"].clone())
            .unwrap();
    let binding =
        GentleEngine::observed_transcript_structures_result_binding(&request, &report).unwrap();
    assert_eq!(
        serde_json::to_value(receipt.observed_structures_result.as_ref().unwrap()).unwrap(),
        serde_json::to_value(binding).unwrap()
    );
    assert_eq!(
        crate::agent_feedback::AgentExecutionStatus::from_shell_output(&result.output),
        crate::agent_feedback::AgentExecutionStatus::Completed
    );
    assert!(!service.cancel(id));
    let after = engine.read().unwrap();
    assert_eq!(serde_json::to_value(after.state()).unwrap(), state);
    assert_eq!(after.journal_len(), journal);
    assert_eq!(
        (
            after.history_summary().operation_log_count,
            after.undo_available(),
            after.redo_available()
        ),
        (
            history.operation_log_count,
            history.undo_count,
            history.redo_count
        )
    );
}

#[test]
fn observed_annotation_managed_publication_does_not_need_the_live_project_lock() {
    let root = tempfile::tempdir().unwrap();
    let request = GentleEngine::observed_annotation_test_fixture(root.path());
    let service = CommandExecutionService::default();
    let engine = Arc::new(RwLock::new(fixture::live_with_redo()));
    let held = engine.write().unwrap();
    let before = serde_json::to_value(held.state()).unwrap();
    let journal = held.journal_len();
    let redo = held.redo_available();
    let command = format!(
        "rna-reads transcript-annotation-prepare {}",
        crate::engine_shell::quote_shell_arg(&serde_json::to_string(&request).unwrap())
    );
    let parsed = crate::engine_shell::parse_shell_line(&command).unwrap();
    let id = service
        .submit_parsed(engine.clone(), parsed, command, Default::default())
        .unwrap();
    let result = wait(&service, id).unwrap();
    assert!(!result.state_changed);
    let receipt = service.status(id).unwrap();
    assert_eq!(receipt.publication, "artifact_only");
    assert!(!receipt.requires_project_commit);
    assert_eq!(receipt.committed_operations, 0);
    assert_eq!(
        receipt.observed_annotation_result.as_ref().unwrap().outcome,
        "completed"
    );
    assert_eq!(before, serde_json::to_value(held.state()).unwrap());
    assert_eq!((journal, redo), (held.journal_len(), held.redo_available()));
    assert!(!service.cancel(id));
}

#[test]
fn lifecycle_slice3_managed_workflow_publishes_journals_not_project_changes() {
    {
        let service = CommandExecutionService::default();
        let engine = Arc::new(RwLock::new(fixture::live_with_redo()));
        let before = engine.read().unwrap();
        let state = serde_json::to_value(before.state()).unwrap();
        let revision = before.mutation_revision();
        let structural = before.structural_revision();
        let history = before.history_summary();
        let owner = before.instance_id();
        drop(before);
        let command = format!(
            "workflow {}",
            crate::engine_shell::shell_quote(
                &serde_json::to_string(&Workflow {
                    run_id: "managed-synthetic".into(),
                    ops: vec![fixture::snv(), fixture::plan(None, true)]
                })
                .unwrap()
            )
        );
        let id = service
            .submit(engine.clone(), command, Default::default())
            .unwrap();
        assert!(!wait(&service, id).unwrap().state_changed);
        let receipt = service.status(id).unwrap();
        assert_eq!(receipt.publication, "project");
        assert!(receipt.requires_project_commit);
        let after = engine.read().unwrap();
        assert_eq!(serde_json::to_value(after.state()).unwrap(), state);
        assert_eq!(
            (
                after.instance_id(),
                after.mutation_revision(),
                after.structural_revision()
            ),
            (owner, revision, structural)
        );
        assert_eq!(
            (after.undo_available(), after.redo_available()),
            (history.undo_count, history.redo_count)
        );
        assert_eq!(receipt.committed_operations, 2);
        assert_eq!(
            after.history_summary().operation_log_count,
            history.operation_log_count + receipt.committed_operations
        );
        drop(after);
        engine.write().unwrap().redo_last_operation().unwrap();
        assert!(
            engine
                .read()
                .unwrap()
                .state()
                .sequences
                .contains_key("synthetic_temporary")
        );
    }
}

#[test]
fn lifecycle_slice3_managed_inspection_keeps_concurrent_journal_stale_guard() {
    let service = CommandExecutionService::default();
    let engine = Arc::new(RwLock::new(fixture::live_with_redo()));
    let (ready_tx, ready_rx) = mpsc::channel();
    let (go_tx, go_rx) = mpsc::channel();
    let held = service
        .submit_work(
            engine.clone(),
            "synthetic-held-snapshot".into(),
            move |_, _| {
                ready_tx.send(()).unwrap();
                go_rx.recv_timeout(Duration::from_secs(10)).unwrap();
                Ok(ShellRunResult {
                    state_changed: false,
                    output: serde_json::json!({"synthetic":true}),
                })
            },
        )
        .unwrap();
    ready_rx.recv_timeout(Duration::from_secs(10)).unwrap();
    let id = service
        .submit(
            engine.clone(),
            format!(
                "op {}",
                crate::engine_shell::shell_quote(&serde_json::to_string(&fixture::snv()).unwrap())
            ),
            Default::default(),
        )
        .unwrap();
    assert!(!wait(&service, id).unwrap().state_changed);
    go_tx.send(()).unwrap();
    assert!(
        wait(&service, held)
            .unwrap_err()
            .to_lowercase()
            .contains("stale")
    );
}

#[test]
fn lifecycle_slice3_managed_cancel_discards_journal_not_already_written_report() {
    let service = CommandExecutionService::default();
    let engine = Arc::new(RwLock::new(fixture::live_with_redo()));
    let baseline = engine.read().unwrap().history_summary().operation_log_count;
    let td = tempfile::tempdir().unwrap();
    let report = td.path().join("synthetic-preview.json");
    let op = fixture::plan(Some(report.display().to_string()), true);
    let (ready_tx, ready_rx) = mpsc::channel();
    let (go_tx, go_rx) = mpsc::channel();
    let id = service
        .submit_work(
            engine.clone(),
            "synthetic-cancel-window".into(),
            move |engine, _| {
                let result = engine.apply(op).map_err(|e| e.to_string())?;
                ready_tx.send(()).unwrap();
                go_rx.recv_timeout(Duration::from_secs(10)).unwrap();
                Ok(ShellRunResult {
                    state_changed: false,
                    output: serde_json::to_value(result).unwrap(),
                })
            },
        )
        .unwrap();
    ready_rx.recv_timeout(Duration::from_secs(10)).unwrap();
    assert!(report.is_file());
    assert!(service.cancel(id));
    go_tx.send(()).unwrap();
    assert!(wait(&service, id).is_err());
    assert!(report.is_file());
    assert_eq!(
        engine.read().unwrap().history_summary().operation_log_count,
        baseline
    );
    assert_eq!(service.status(id).unwrap().committed_operations, 0);
}

#[test]
fn lifecycle_slice4_managed_worker_copies_binding_and_commit_keeps_live_save_as() {
    let service = CommandExecutionService::default();
    let mut live = fixture::live_with_redo();
    let td = tempfile::tempdir().unwrap();
    let old = td.path().join("old.json");
    let new = td.path().join("new.json");
    let report = td.path().join("worker-report.json");
    live.set_host_project_path(Some(&old)).unwrap();
    let engine = Arc::new(RwLock::new(live));
    let (ready_tx, ready_rx) = mpsc::channel();
    let (go_tx, go_rx) = mpsc::channel();
    let operation = fixture::plan(Some(report.display().to_string()), true);
    let id = service
        .submit_work(
            engine.clone(),
            "synthetic-host-binding".into(),
            move |worker, _| {
                assert_eq!(worker.host_project_path(), Some(old.as_path()));
                let output = worker.apply(operation).map_err(|e| e.to_string())?;
                ready_tx.send(()).unwrap();
                go_rx.recv_timeout(Duration::from_secs(10)).unwrap();
                Ok(ShellRunResult {
                    state_changed: false,
                    output: serde_json::to_value(output).unwrap(),
                })
            },
        )
        .unwrap();
    ready_rx.recv_timeout(Duration::from_secs(10)).unwrap();
    engine
        .write()
        .unwrap()
        .set_host_project_path(Some(&new))
        .unwrap();
    go_tx.send(()).unwrap();
    assert!(!wait(&service, id).unwrap().state_changed);
    assert!(report.is_file());
    assert_eq!(
        engine.read().unwrap().host_project_path(),
        Some(new.as_path())
    );
    assert_eq!(service.status(id).unwrap().committed_operations, 1);
}
