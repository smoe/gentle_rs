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
fn lifecycle_slice3_managed_op_and_workflow_publish_journals_not_project_changes() {
    for workflow in [false, true] {
        let service = CommandExecutionService::default();
        let engine = Arc::new(RwLock::new(fixture::live_with_redo()));
        let before = engine.read().unwrap();
        let state = serde_json::to_value(before.state()).unwrap();
        let revision = before.mutation_revision();
        let structural = before.structural_revision();
        let history = before.history_summary();
        let owner = before.instance_id();
        drop(before);
        let command = if workflow {
            format!(
                "workflow {}",
                crate::engine_shell::shell_quote(
                    &serde_json::to_string(&Workflow {
                        run_id: "managed-synthetic".into(),
                        ops: vec![fixture::snv(), fixture::plan(None, true)]
                    })
                    .unwrap()
                )
            )
        } else {
            format!(
                "op {}",
                crate::engine_shell::shell_quote(
                    &serde_json::to_string(&fixture::plan(None, true)).unwrap()
                )
            )
        };
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
        assert_eq!(receipt.committed_operations, if workflow { 2 } else { 1 });
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
