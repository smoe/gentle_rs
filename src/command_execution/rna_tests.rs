//! Synthetic shared-demo execution regressions; channels control race windows, never private reads.

use super::*;
use crate::engine::{Engine, Operation, ProjectState, ScreenRnaReadPanelRequest};
use std::{
    fs,
    path::Path,
    sync::mpsc,
    time::{Duration, Instant},
};

fn demo(root: &Path) -> (Arc<RwLock<GentleEngine>>, ScreenRnaReadPanelRequest) {
    let directory = root.join("demo");
    GentleEngine::new()
        .prepare_rna_read_panel_demo(&directory.display().to_string())
        .unwrap();
    let state: ProjectState =
        serde_json::from_slice(&fs::read(directory.join("artifacts/starter.json")).unwrap())
            .unwrap();
    let engine = GentleEngine::from_state(state);
    let build =
        serde_json::from_slice(&fs::read(directory.join("artifacts/build.json")).unwrap()).unwrap();
    let panel = engine.build_rna_read_panel(&build).unwrap();
    let mut request: ScreenRnaReadPanelRequest =
        serde_json::from_slice(&fs::read(directory.join("artifacts/screen.json")).unwrap())
            .unwrap();
    request.expected_panel_sha256 = Some(panel.content_sha256);
    (Arc::new(RwLock::new(engine)), request)
}

fn command(request: &ScreenRnaReadPanelRequest) -> ShellCommand {
    ShellCommand::Op {
        payload: serde_json::to_string(&Operation::ScreenRnaReadPanel {
            request: Box::new(request.clone()),
        })
        .unwrap(),
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
            "job did not finish: {:?}",
            service.status(id)
        );
        thread::sleep(Duration::from_millis(2));
    }
}

#[test]
fn rna_panel_service_late_cancel_and_live_edit_keep_published_completion() {
    let td = tempfile::tempdir().unwrap();
    let (live, request) = demo(td.path());
    let baseline = live.read().unwrap().journal_len();
    let service = CommandExecutionService::default();
    let (ready_tx, ready_rx) = mpsc::channel();
    let (release_tx, release_rx) = mpsc::channel();
    let cancel = Arc::new(AtomicBool::new(false));
    let caller_cancel = cancel.clone();
    let command = command(&request);
    let id = service
        .submit_work_with_cancel(
            live.clone(),
            "synthetic-command".into(),
            cancel,
            Lane::Execution,
            false,
            move |engine, callback| {
                let result = execute_command(
                    engine,
                    command,
                    ShellExecutionOptions::default(),
                    callback,
                    caller_cancel,
                )?;
                ready_tx.send(()).unwrap();
                release_rx.recv_timeout(Duration::from_secs(20)).unwrap();
                Ok(result)
            },
        )
        .unwrap();
    ready_rx.recv_timeout(Duration::from_secs(20)).unwrap();
    assert!(Path::new(&request.output_dir).join("report.json").exists());
    live.write()
        .unwrap()
        .apply(Operation::CreateSequenceFromText {
            sequence_text: "ACGT".into(),
            output_id: Some("concurrent_edit".into()),
            name: None,
            circular: false,
        })
        .unwrap();
    assert!(
        service.cancel(id),
        "cancel arrives after file publication but before host retention"
    );
    release_tx.send(()).unwrap();
    wait(&service, id).unwrap();
    let receipt = service.status(id).unwrap();
    assert_eq!(receipt.state, RuntimeStatusFrameState::Completed);
    assert!(receipt.cancel_requested);
    assert_eq!(
        receipt.artifact_result.as_ref().unwrap().outcome,
        "completed"
    );
    assert_eq!(receipt.committed_operations, 0);
    assert!(!receipt.requires_project_commit);
    assert_eq!(live.read().unwrap().journal_len(), baseline + 1);
    assert!(
        live.read()
            .unwrap()
            .state()
            .sequences
            .contains_key("concurrent_edit")
    );
    assert!(!service.cancel(id));
}

#[test]
fn rna_panel_service_retention_does_not_stale_another_detached_job() {
    let td = tempfile::tempdir().unwrap();
    let (live, request) = demo(td.path());
    let service = CommandExecutionService::default();
    let (ready_tx, ready_rx) = mpsc::channel();
    let (release_tx, release_rx) = mpsc::channel();
    let mutation = service
        .submit_work(
            live.clone(),
            "synthetic-mutation".into(),
            move |engine, _| {
                engine
                    .apply(Operation::CreateSequenceFromText {
                        sequence_text: "ACGT".into(),
                        output_id: Some("kept".into()),
                        name: None,
                        circular: false,
                    })
                    .unwrap();
                ready_tx.send(()).unwrap();
                release_rx.recv_timeout(Duration::from_secs(20)).unwrap();
                Ok(ShellRunResult {
                    state_changed: true,
                    output: serde_json::json!({"ok":true}),
                })
            },
        )
        .unwrap();
    ready_rx.recv_timeout(Duration::from_secs(20)).unwrap();
    let before = live.read().unwrap().journal_len();
    let id = service
        .submit_parsed(
            live.clone(),
            command(&request),
            "synthetic-screen".into(),
            ShellExecutionOptions::default(),
        )
        .unwrap();
    wait(&service, id).unwrap();
    assert_eq!(live.read().unwrap().journal_len(), before);
    release_tx.send(()).unwrap();
    wait(&service, mutation).unwrap();
    assert!(live.read().unwrap().state().sequences.contains_key("kept"));
}

#[test]
fn rna_panel_service_partial_cancel_and_control_error_keep_bound_reports() {
    for fault in [false, true] {
        let td = tempfile::tempdir().unwrap();
        let (live, request) = demo(td.path());
        let marker = request.cancel_path.clone().unwrap();
        let options = ShellExecutionOptions {
            progress_callback: Some(Arc::new(Mutex::new(Box::new(move |event| {
                if matches!(event, OperationProgress::RnaReadPanel(ref p) if p.phase == "file_checksum_verified" && p.file_ordinal == 1)
                {
                    if fault {
                        fs::write(&marker, "foreign marker").unwrap();
                    } else {
                        return false;
                    }
                }
                true
            })))),
            ..Default::default()
        };
        let service = CommandExecutionService::default();
        let id = service
            .submit_parsed(
                live.clone(),
                command(&request),
                "synthetic-partial".into(),
                options,
            )
            .unwrap();
        let output = wait(&service, id).unwrap();
        let receipt = service.status(id).unwrap();
        let binding = receipt.artifact_result.unwrap();
        assert_eq!(
            binding.outcome,
            if fault { "control_error" } else { "cancelled" }
        );
        assert_eq!(
            receipt.state,
            if fault {
                RuntimeStatusFrameState::Failed
            } else {
                RuntimeStatusFrameState::Cancelled
            }
        );
        assert_eq!(
            output.output["result"]["rna_read_panel_screen"]["files_completed"]
                .as_array()
                .unwrap()
                .len(),
            1
        );
        let inspection = GentleEngine::new()
            .inspect_rna_read_panel(&binding.report_path, None, 0)
            .unwrap();
        GentleEngine::validate_rna_read_panel_result_binding(
            &binding,
            &binding.report_path,
            inspection.1.as_ref(),
            None,
        )
        .unwrap();
        assert_eq!(live.read().unwrap().journal_len(), 0);
    }
}

#[test]
fn rna_panel_control_and_metadata_do_not_compete_with_execution_slots_or_project_lock() {
    let td = tempfile::tempdir().unwrap();
    let (live, request) = demo(td.path());
    let service = CommandExecutionService::default();
    let mut releases = vec![];
    let mut jobs = vec![];
    for _ in 0..MAX_ACTIVE {
        let (ready_tx, ready_rx) = mpsc::channel();
        let (release_tx, release_rx) = mpsc::channel();
        jobs.push(
            service
                .submit_work(live.clone(), "held".into(), move |_, _| {
                    ready_tx.send(()).unwrap();
                    release_rx.recv_timeout(Duration::from_secs(20)).unwrap();
                    Ok(ShellRunResult {
                        state_changed: false,
                        output: serde_json::json!({"ok":true}),
                    })
                })
                .unwrap(),
        );
        ready_rx.recv_timeout(Duration::from_secs(20)).unwrap();
        releases.push(release_tx);
    }
    let guard = live.write().unwrap();
    let cancel = ShellCommand::Op {
        payload: serde_json::to_string(&Operation::CancelRnaReadPanel {
            request: Box::new(request.clone()),
        })
        .unwrap(),
    };
    let id = service
        .submit_parsed(
            live.clone(),
            cancel,
            "cancel".into(),
            ShellExecutionOptions::default(),
        )
        .unwrap();
    assert_eq!(wait(&service, id).unwrap().output["cancel_requested"], true);
    let inspect = Operation::InspectRnaReadPanel {
        path: request.panel_path.clone(),
        expected_result: None,
        verify_artifacts: false,
        file_id: None,
        rows: 0,
    };
    let id = service
        .submit_parsed(
            live.clone(),
            ShellCommand::Op {
                payload: serde_json::to_string(&inspect).unwrap(),
            },
            "inspect".into(),
            ShellExecutionOptions::default(),
        )
        .unwrap();
    wait(&service, id).unwrap();
    assert!(service.status(id).unwrap().owner_instance.is_none());
    let expensive = Operation::InspectRnaReadPanel {
        path: request.panel_path,
        expected_result: None,
        verify_artifacts: true,
        file_id: None,
        rows: 0,
    };
    assert!(
        service
            .submit_parsed(
                live.clone(),
                ShellCommand::Op {
                    payload: serde_json::to_string(&expensive).unwrap()
                },
                "verify".into(),
                ShellExecutionOptions::default()
            )
            .unwrap_err()
            .contains("capacity")
    );
    drop(guard);
    for release in releases {
        release.send(()).unwrap();
    }
    for id in jobs {
        wait(&service, id).unwrap();
    }
}

#[test]
fn rna_panel_service_allowlist_and_large_result_retention_are_explicit() {
    let mutator = ShellCommand::Op {
        payload: serde_json::to_string(&Operation::CreateSequenceFromText {
            sequence_text: "ACGT".into(),
            output_id: None,
            name: None,
            circular: false,
        })
        .unwrap(),
    };
    assert!(rna_artifacts::admission(&mutator).needs_project());
    let td = tempfile::tempdir().unwrap();
    let (live, request) = demo(td.path());
    let command = command(&request);
    assert!(!rna_artifacts::admission(&command).needs_project());
    let service = CommandExecutionService::default();
    let cancel = Arc::new(AtomicBool::new(false));
    let caller_cancel = cancel.clone();
    let id = service
        .submit_work_with_cancel(
            live.clone(),
            "large-synthetic".into(),
            cancel,
            Lane::Execution,
            false,
            move |engine, callback| {
                let mut result = execute_command(
                    engine,
                    command,
                    ShellExecutionOptions::default(),
                    callback,
                    caller_cancel,
                )?;
                result.run.output["padding"] = serde_json::json!("x".repeat(MAX_RESULT_BYTES));
                Ok(result)
            },
        )
        .unwrap();
    let output = wait(&service, id).unwrap();
    assert_eq!(output.output["schema"], "gentle.command_artifact_result.v1");
    assert_eq!(
        crate::agent_feedback::AgentExecutionStatus::from_shell_output(&output.output),
        crate::agent_feedback::AgentExecutionStatus::Completed
    );
    let binding = service.status(id).unwrap().artifact_result.unwrap();
    assert!(Path::new(&binding.report_path).is_file());
    assert_eq!(binding.outcome, "completed");
    assert_eq!(live.read().unwrap().journal_len(), 0);
}

#[test]
fn rna_panel_registered_build_keeps_project_commit_and_mutation_feedback() {
    let td = tempfile::tempdir().unwrap();
    let (live, _) = demo(td.path());
    let mut build: crate::engine::BuildRnaReadPanelRequest =
        serde_json::from_slice(&fs::read(td.path().join("demo/artifacts/build.json")).unwrap())
            .unwrap();
    build.output_path = td.path().join("registered.json").display().to_string();
    build.register_in_project = Some(true);
    let command = ShellCommand::Op {
        payload: serde_json::to_string(&Operation::BuildRnaReadPanel {
            request: Box::new(build),
        })
        .unwrap(),
    };
    assert!(rna_artifacts::admission(&command).needs_project());
    let service = CommandExecutionService::default();
    let id = service
        .submit_parsed(
            live.clone(),
            command,
            "registered-build".into(),
            ShellExecutionOptions::default(),
        )
        .unwrap();
    assert!(wait(&service, id).unwrap().state_changed);
    let receipt = service.status(id).unwrap();
    assert!(receipt.requires_project_commit);
    assert_eq!(receipt.publication, "project");
    assert!(receipt.artifact_result.is_none());
    assert_eq!(
        GentleEngine::rna_read_panel_project_entries(live.read().unwrap().state())
            .unwrap()
            .len(),
        1
    );
}
