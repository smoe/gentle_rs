//! Synthetic FASTQ and caller-defined oligos test managed publication, never private libraries.
//! All inputs are created here in disposable directories; channels pin publication races.

use super::*;
use crate::agent_feedback::AgentExecutionStatus;
use crate::engine::{Engine, Operation, RnaPrimerPanel, TracePrimerReadsRequest};
use std::{
    fs,
    path::Path,
    sync::mpsc,
    time::{Duration, Instant},
};

fn request(root: &Path) -> TracePrimerReadsRequest {
    let mut panel: RnaPrimerPanel = serde_json::from_value(serde_json::json!({
        "schema":"gentle.rna_primer_panel.v1","content_sha256":"",
        "primers":[{"primer_id":"synthetic","assay_label":"Managed trace test",
            "sequence":"ACGTCAGTACCGTTAGCCTAGG",
            "components":[{"component_id":"core","role":"target_core","start":0,"end":22}],
            "core_orientation":"sense_5prime","expected_targets":[],"variant":"synthetic",
            "source":"Hand-written managed-worker regression, not a physical oligo"}],"policy":{}}))
    .unwrap();
    GentleEngine::seal_rna_primer_panel(&mut panel).unwrap();
    let sequence = format!(
        "GG{}GATCCATGAGCTAGTCGACATCGTAGCTACGATGCTAGCATCG",
        panel.primers[0].sequence
    );
    let files = (0..2).map(|i| {
        let bytes = format!("@invented_{i}\n{sequence}\n+\n{}\n", "I".repeat(sequence.len()));
        let path = root.join(format!("input's {i}.fastq"));
        fs::write(&path, bytes.as_bytes()).unwrap();
        serde_json::json!({"sample_id":"invented","barcode":"synthetic","priming":"synthetic",
            "path":path,"declared_complete":true,"expected_sha256":sha256_prefixed_bytes(bytes.as_bytes())})
    }).collect::<Vec<_>>();
    serde_json::from_value(
        serde_json::json!({"schema":"gentle.rna_primer_trace_request.v1",
        "primer_panel":{"kind":"inline","panel":panel},"files":files,
        "output_dir":root.join("fresh output"),"cancel_path":root.join("cancel marker")}),
    )
    .unwrap()
}

fn command(request: &TracePrimerReadsRequest) -> ShellCommand {
    ShellCommand::Op {
        payload: serde_json::to_string(&Operation::TracePrimerReads {
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
            "trace job did not finish: {:?}",
            service.status(id)
        );
        thread::sleep(Duration::from_millis(2));
    }
}

#[test]
fn primer_trace_managed_partial_reports_keep_outcomes_without_project_lock() {
    for fault in [false, true] {
        let td = tempfile::tempdir().unwrap();
        let request = request(td.path());
        let marker = request.cancel_path.clone().unwrap();
        let live = Arc::new(RwLock::new(GentleEngine::new()));
        let guard = live.write().unwrap();
        let before = serde_json::to_value(guard.state()).unwrap();
        let journal = guard.journal_len();
        let options = ShellExecutionOptions {
            progress_callback: Some(Arc::new(Mutex::new(Box::new(move |event| {
                if matches!(event, OperationProgress::RnaReadPanel(ref p) if p.file_ordinal == 2 && p.processed_reads == 0)
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
                "synthetic-trace".into(),
                options,
            )
            .unwrap();
        let output = wait(&service, id).unwrap();
        let receipt = service.status(id).unwrap();
        let binding = receipt.primer_trace_result.as_ref().unwrap();
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
            AgentExecutionStatus::from_shell_output(&output.output),
            if fault {
                AgentExecutionStatus::Failed
            } else {
                AgentExecutionStatus::Cancelled
            }
        );
        assert_eq!(receipt.publication, "artifact_only");
        assert!(!receipt.requires_project_commit);
        assert!(receipt.owner_instance.is_none());
        assert!(receipt.artifact_result.is_none());
        assert_eq!(receipt.committed_operations, 0);
        let inspection = GentleEngine::inspect_primer_trace(
            &binding.report_path,
            Some(binding),
            Some("file_0000"),
            1,
        )
        .unwrap();
        assert_eq!(inspection.report.files_completed.len(), 1);
        assert_eq!(inspection.report.inputs_not_completed.len(), 1);
        assert_eq!(inspection.rows.len(), 1);
        assert!(!inspection.artifact_content_verified);
        assert_eq!(serde_json::to_value(guard.state()).unwrap(), before);
        assert_eq!(guard.journal_len(), journal);
    }
}

#[test]
fn primer_trace_managed_compact_result_survives_late_cancel_and_project_replacement() {
    let td = tempfile::tempdir().unwrap();
    let request = request(td.path());
    let live = Arc::new(RwLock::new(GentleEngine::new()));
    let service = CommandExecutionService::default();
    let (ready_tx, ready_rx) = mpsc::channel();
    let (release_tx, release_rx) = mpsc::channel();
    let cancel = Arc::new(AtomicBool::new(false));
    let caller_cancel = cancel.clone();
    let command = command(&request);
    let id = service
        .submit_work_with_cancel(
            live.clone(),
            "synthetic-trace".into(),
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
                ready_tx.send(()).unwrap();
                release_rx.recv_timeout(Duration::from_secs(20)).unwrap();
                Ok(result)
            },
        )
        .unwrap();
    ready_rx.recv_timeout(Duration::from_secs(20)).unwrap();
    assert!(Path::new(&request.output_dir).join("report.json").is_file());
    let replacement = GentleEngine::new();
    let replacement_instance = replacement.instance_id();
    *live.write().unwrap() = replacement;
    assert!(service.cancel(id));
    release_tx.send(()).unwrap();
    let output = wait(&service, id).unwrap();
    let receipt = service.status(id).unwrap();
    assert_eq!(receipt.state, RuntimeStatusFrameState::Completed);
    assert!(receipt.cancel_requested);
    assert_eq!(
        output.output["schema"],
        "gentle.command_primer_trace_result.v1"
    );
    assert_eq!(
        AgentExecutionStatus::from_shell_output(&output.output),
        AgentExecutionStatus::Completed
    );
    let binding = receipt.primer_trace_result.as_ref().unwrap();
    assert_eq!(binding.outcome, "completed");
    assert_eq!(receipt.publication, "artifact_only");
    assert_eq!(receipt.committed_operations, 0);
    assert!(!receipt.requires_project_commit);
    assert!(receipt.artifact_result.is_none());
    let inspection =
        GentleEngine::inspect_primer_trace(&binding.report_path, Some(binding), None, 0).unwrap();
    assert_eq!(inspection.report.files_completed.len(), 2);
    assert_eq!(live.read().unwrap().instance_id(), replacement_instance);
    assert_eq!(live.read().unwrap().journal_len(), 0);
    assert!(!service.cancel(id));
}
