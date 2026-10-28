//! Synthetic S2.1c shared-route, approval and lifecycle tests using the inline GTF fixture.

use super::*;
use crate::engine::{Engine, LifecycleRequest, Operation, Workflow, lifecycle_slice3_tests};
use crate::engine_shell::{execute_shell_command, parse_shell_line, quote_shell_arg};

#[test]
fn observed_annotation_forged_agent_flags_cannot_bypass_write_confirmation() {
    use crate::agent_planner::{AgentPlanCandidate, AgentPlanCandidateKind, AgentPlanResult};
    let root = tempfile::tempdir().unwrap();
    let mut r = fixture(root.path());
    r.cancel_path = Some(root.path().join("unwritten marker").display().to_string());
    for cancel in [false, true] {
        let op = if cancel {
            Operation::CancelObservedTranscriptAnnotation {
                request: Box::new(r.clone()),
            }
        } else {
            Operation::PrepareObservedTranscriptAnnotation {
                request: Box::new(r.clone()),
            }
        };
        let candidate = AgentPlanCandidate {
            candidate_id: "annotation".into(),
            title: "Synthetic write approval".into(),
            kind: AgentPlanCandidateKind::Op,
            execution_mode: crate::agent_planner::AgentPlanExecutionMode::Auto,
            operation: Some(serde_json::to_value(op).unwrap()),
            mutating: false,
            requires_confirmation: false,
            ..Default::default()
        };
        let plan = AgentPlanResult {
            schema: "gentle.agent_plan_result.v1".into(),
            assistant_message: String::new(),
            questions: vec![],
            candidates: vec![candidate],
        };
        let mut engine = lifecycle_slice3_tests::live_with_redo();
        let before = serde_json::to_value(engine.state()).unwrap();
        let journal = engine.journal_len();
        let attempt = crate::agent_execution::execute_agent_plan_candidate_with_feedback(
            &mut engine,
            &plan,
            "annotation",
            false,
            &Default::default(),
        );
        assert!(attempt.result.unwrap_err().contains("explicit --confirm"));
        assert_eq!(
            attempt.feedback.status,
            crate::agent_feedback::AgentExecutionStatus::Blocked
        );
        assert!(!Path::new(&r.output_dir).exists());
        assert!(!Path::new(r.cancel_path.as_ref().unwrap()).exists());
        assert_eq!(before, serde_json::to_value(engine.state()).unwrap());
        assert_eq!(journal, engine.journal_len());
    }
}

#[test]
fn observed_annotation_operation_preserves_constructor_content_redo_and_revisions() {
    let root = tempfile::tempdir().unwrap();
    let mut r = fixture(root.path());
    r.cancel_path = Some(root.path().join("cancel marker").display().to_string());
    let op = Operation::PrepareObservedTranscriptAnnotation {
        request: Box::new(r.clone()),
    };
    let initial = lifecycle_slice3_tests::state();
    let constructed =
        GentleEngine::from_state_for(initial.clone(), LifecycleRequest::Operation(&op));
    assert_eq!(
        serde_json::to_value(constructed.state()).unwrap(),
        serde_json::to_value(initial).unwrap()
    );
    let mut engine = lifecycle_slice3_tests::live_with_redo();
    let before = serde_json::to_value(engine.state()).unwrap();
    let revisions = (engine.mutation_revision(), engine.structural_revision());
    let history = engine.history_summary();
    let journal = engine.journal_len();
    let result = engine.apply(op.clone()).unwrap();
    assert!(result.observed_transcript_annotation.is_some());
    assert_eq!(before, serde_json::to_value(engine.state()).unwrap());
    assert_eq!(
        revisions,
        (engine.mutation_revision(), engine.structural_revision())
    );
    assert_eq!(
        (engine.undo_available(), engine.redo_available()),
        (history.undo_count, history.redo_count)
    );
    assert_eq!(engine.journal_len(), journal + 1);
    let cancel = Operation::CancelObservedTranscriptAnnotation {
        request: Box::new(r),
    };
    assert!(LifecycleRequest::Operation(&cancel).preserves_project_content());
    let workflow = Workflow {
        run_id: "synthetic annotation".into(),
        ops: vec![op, cancel],
    };
    assert!(LifecycleRequest::Workflow(&workflow).preserves_project_content());
    let mixed = Workflow {
        run_id: "mixed".into(),
        ops: vec![workflow.ops[0].clone(), lifecycle_slice3_tests::create()],
    };
    assert!(!LifecycleRequest::Workflow(&mixed).preserves_project_content());
}

#[test]
fn observed_annotation_shared_parser_quotes_lf_crlf_requests_and_discovery_requires_approval() {
    let root = tempfile::tempdir().unwrap();
    let original = fixture(root.path());
    let mut engine = GentleEngine::new();
    for crlf in [false, true] {
        let mut r = original.clone();
        r.output_dir = root
            .path()
            .join(format!("output {crlf}"))
            .display()
            .to_string();
        r.cancel_path = Some(
            root.path()
                .join(format!("marker {crlf}"))
                .display()
                .to_string(),
        );
        let path = root.path().join("request's file.json");
        let raw = serde_json::to_string_pretty(&r).unwrap();
        fs::write(&path, if crlf { raw.replace('\n', "\r\n") } else { raw }).unwrap();
        let payload = quote_shell_arg(&format!("@{}", path.display()));
        let command = parse_shell_line(&format!(
            "rna-reads transcript-annotation-prepare {payload}"
        ))
        .unwrap();
        assert!(command.is_state_mutating());
        let before = serde_json::to_value(engine.state()).unwrap();
        let run = execute_shell_command(&mut engine, &command).unwrap();
        assert!(!run.state_changed);
        assert_eq!(
            run.output["result"]["observed_transcript_annotation"]["report"]["analysis"]["outcome"],
            "completed"
        );
        assert_eq!(before, serde_json::to_value(engine.state()).unwrap());
        let journal = engine.journal_len();
        let cancel =
            parse_shell_line(&format!("rna-reads transcript-annotation-cancel {payload}")).unwrap();
        assert!(cancel.is_state_mutating());
        let cancelled = execute_shell_command(&mut engine, &cancel).unwrap();
        assert!(!cancelled.state_changed);
        assert_eq!(cancelled.output["termination_confirmed"], false);
        assert_eq!(engine.journal_len(), journal);
        assert_eq!(before, serde_json::to_value(engine.state()).unwrap());
    }
    let discovered = execute_shell_command(
        &mut engine,
        &parse_shell_line("introspect capabilities --kind operation").unwrap(),
    )
    .unwrap();
    for id in [
        "PrepareObservedTranscriptAnnotation",
        "CancelObservedTranscriptAnnotation",
    ] {
        let row = discovered.output["capabilities"]
            .as_array()
            .unwrap()
            .iter()
            .find(|d| d["id"] == id)
            .unwrap();
        assert_eq!(row["requires_confirmation"], true);
        assert_eq!(row["mutating"], "external");
        let op: Operation = serde_json::from_value(row["operation_template"].clone()).unwrap();
        assert!(LifecycleRequest::Operation(&op).preserves_project_content());
        assert!(
            crate::engine::operation_approval_policy(&op)
                .unwrap()
                .requires_confirmation()
        );
    }
}

#[test]
fn observed_annotation_cancel_is_idempotent_request_bound_and_never_clobbers_sources() {
    let root = tempfile::tempdir().unwrap();
    let mut r = fixture(root.path());
    let marker = root.path().join("cancel marker");
    r.cancel_path = Some(marker.display().to_string());
    let original = fs::read(&r.gtf.path).unwrap();
    let expected = semantic_hash(&r, r.budgets.max_retained_bytes).unwrap();
    GentleEngine::cancel_observed_transcript_annotation(&r).unwrap();
    assert_eq!(fs::read_to_string(&marker).unwrap(), expected);
    assert!(
        GentleEngine::cancel_observed_transcript_annotation(&r)
            .unwrap()
            .contains("already")
    );
    let mut wrong = r.clone();
    wrong.selector.gene_id = "another_gene".into();
    assert!(GentleEngine::cancel_observed_transcript_annotation(&wrong).is_err());
    let mut alias = r.clone();
    alias.cancel_path = Some(alias.gtf.path.clone());
    assert!(GentleEngine::cancel_observed_transcript_annotation(&alias).is_err());
    assert_eq!(fs::read(&r.gtf.path).unwrap(), original);
    let cancelled = run(r.clone());
    assert!(cancelled.report.analysis.is_none());
    assert_eq!(
        cancelled.report.status,
        ObservedAnnotationExecutionStatus::Cancelled
    );
    // Cancellation does not replay source data or require source bytes after publication.
    fs::remove_file(&r.gtf.path).unwrap();
    fs::remove_file(&r.dictionary.path).unwrap();
    assert!(GentleEngine::cancel_observed_transcript_annotation(&r).is_ok());
}

#[test]
fn observed_annotation_result_binding_checks_request_receipt_path_and_raw_bytes() {
    let root = tempfile::tempdir().unwrap();
    let r = fixture(root.path());
    let p = run(r.clone());
    let binding = GentleEngine::observed_transcript_annotation_result_binding(&r, &p).unwrap();
    assert_eq!(binding.outcome, "completed");
    assert_eq!(
        binding.gene_outcome,
        Some(ObservedAnnotationOutcome::Completed)
    );
    for i in 0..5 {
        let mut changed = p.clone();
        match i {
            0 => changed.report.request.selector.gene_id = "other".into(),
            1 => changed.report.content_sha256 = format!("sha256:{}", "0".repeat(64)),
            2 => {
                changed.report_file.path =
                    root.path().join("wrong report.json").display().to_string()
            }
            3 => changed.report_file.sha256 = format!("sha256:{}", "0".repeat(64)),
            _ => changed.report_file.bytes += 1,
        }
        assert!(GentleEngine::observed_transcript_annotation_result_binding(&r, &changed).is_err());
    }
    let mut forged = p;
    forged.report.analysis = None;
    forged.report.annotation = None;
    forged.report.content_sha256.clear();
    forged.report.content_sha256 =
        semantic_hash(&forged.report, r.budgets.max_retained_bytes).unwrap();
    assert!(GentleEngine::observed_transcript_annotation_result_binding(&r, &forged).is_err());
}

#[test]
fn observed_annotation_agent_outcome_keeps_refusal_and_control_separate() {
    use crate::agent_feedback::AgentExecutionStatus as Status;
    let root = tempfile::tempdir().unwrap();
    let original = fixture(root.path());
    for (i, expected) in [
        Status::Completed,
        Status::Partial,
        Status::Cancelled,
        Status::Failed,
    ]
    .into_iter()
    .enumerate()
    {
        let mut r = original.clone();
        r.output_dir = root
            .path()
            .join(format!("status {i}"))
            .display()
            .to_string();
        if i == 1 {
            r.selector.gene_id = "missing".into();
        }
        if i == 3 {
            let marker = root.path().join("wrong marker");
            fs::write(&marker, b"wrong identity").unwrap();
            r.cancel_path = Some(marker.display().to_string());
        }
        let p = GentleEngine::new()
            .prepare_observed_transcript_annotation(r.clone(), &mut |_| i != 2)
            .unwrap();
        let binding = GentleEngine::observed_transcript_annotation_result_binding(&r, &p).unwrap();
        let compact = json!({"schema":"gentle.command_observed_annotation_result.v1","observed_annotation_result":binding});
        let full = json!({"result":{"observed_transcript_annotation":p}});
        assert_eq!(Status::from_shell_output(&compact), expected);
        assert_eq!(Status::from_shell_output(&full), expected);
    }
}
