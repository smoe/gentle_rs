//! Temporary hand-crafted SAM/GTF file-boundary fixtures, not real RNA or aligner evidence.
//! Recreate through `observed_sam_` tests; all inputs/outputs stay in private tempdirs.

use super::*;

fn bound(path: &Path) -> PrimerPayloadFile {
    let raw = fs::read(path).unwrap();
    PrimerPayloadFile {
        path: path.display().to_string(),
        bytes: raw.len() as u64,
        sha256: sha256_prefixed_bytes(&raw),
    }
}

pub(crate) fn fixture() -> (tempfile::TempDir, Request) {
    let (root, input, _) = kernel::sam_intake_test_fixture();
    let source = &input.annotation_sources[0];
    let annotation = root.path().join("normalized annotation.json");
    fs::write(
        &annotation,
        serde_json::to_vec_pretty(&source.annotation).unwrap(),
    )
    .unwrap();
    let r = Request {
        schema: OBSERVED_SAM_REQUEST_SCHEMA.into(),
        sam: input.sam,
        manifest: input.manifest,
        dictionary: input.dictionary,
        conversion: input.conversion,
        reference_id: input.reference_id,
        reference: input.reference,
        reference_gene: input.reference_gene.clone(),
        comparison_genes: input.comparison_genes,
        annotation_sources: vec![ObservedSamAnnotationBinding::Normalized {
            selector: input.reference_gene,
            provider: source.provider.clone(),
            annotation: bound(&annotation),
        }],
        version_policy: input.version_policy,
        unpadded_reference: input.unpadded_reference,
        required_flank_bases: input.required_flank_bases,
        text_budgets: input.limits,
        kernel_budgets: input.kernel_budgets,
        budgets: ObservedSamBudgets {
            max_input_bytes: 1_000_000,
            max_annotation_bytes: 65536,
            max_output_bytes: 2 * 1024 * 1024,
            max_output_row_bytes: 65536,
            max_units: 100,
            max_retained_bytes: 8 * 1024 * 1024,
            max_projection_work: 1_000_000,
        },
        output_dir: root
            .path()
            .join("SAM gene's evidence")
            .display()
            .to_string(),
        cancel_path: None,
    };
    (root, r)
}

fn run(r: &Request) -> ObservedSamReport {
    GentleEngine::new()
        .summarize_observed_transcript_structures_from_sam(r.clone(), &mut |_| true)
        .unwrap()
}

fn data(report: &ObservedSamReport) -> ObservedSamData {
    serde_json::from_slice(
        &fs::read(
            Path::new(&report.resolved_output_dir)
                .join("artifacts")
                .join(DATA),
        )
        .unwrap(),
    )
    .unwrap()
}

#[test]
fn observed_sam_report_binding_and_saved_inspection_work_when_all_original_sources_are_offline() {
    let (_root, r) = fixture();
    let engine = GentleEngine::new();
    let before = encode(engine.state()).unwrap();
    let journal = engine.journal_len();
    let report = engine
        .summarize_observed_transcript_structures_from_sam(r.clone(), &mut |_| true)
        .unwrap();
    assert_eq!(report.status, ObservedSamExecutionStatus::Completed);
    assert_eq!(
        report.coordinate_authority,
        ObservedSamCoordinateAuthority::Declared
    );
    assert!(report.reference_verification.is_none());
    assert_eq!(report.counts.as_ref().unwrap().query_units, 4);
    assert_eq!(report.artifacts.len(), 7);
    assert_eq!(encode(engine.state()).unwrap(), before);
    assert_eq!(engine.journal_len(), journal);
    let binding = GentleEngine::observed_sam_result_binding(&r, &report).unwrap();
    let all = data(&report);
    assert_eq!(all.schema, OBSERVED_SAM_DATA_SCHEMA);
    assert_eq!(all.analysis.as_ref().unwrap().units.len(), 4);
    for file in input_files(&r) {
        fs::remove_file(&file.path).unwrap();
    }
    let inspected = GentleEngine::inspect_observed_transcript_structures_from_sam(
        &binding.report_path,
        Some(&binding),
        TABLES[0],
        4,
        true,
    )
    .unwrap();
    assert!(inspected.artifacts_verified);
    assert_eq!(inspected.rows.len(), 4);
    assert_eq!(inspected.rows[0]["unit"]["source_sha256"], r.sam.sha256);
    assert!(inspected.rows[0].get("record_index").is_none());
    assert_eq!(
        serde_json::to_value(&inspected.report).unwrap(),
        serde_json::to_value(&report).unwrap()
    );
}

#[test]
fn observed_sam_corrupt_artifacts_do_not_pass_explicit_verification_or_preview() {
    let (_root, r) = fixture();
    let report = run(&r);
    let binding = GentleEngine::observed_sam_result_binding(&r, &report).unwrap();
    fs::write(
        Path::new(&report.resolved_output_dir)
            .join("artifacts")
            .join(TABLES[0]),
        "row_json\n{}\n",
    )
    .unwrap();
    let metadata = GentleEngine::inspect_observed_transcript_structures_from_sam(
        &binding.report_path,
        Some(&binding),
        TABLES[0],
        0,
        false,
    )
    .unwrap();
    assert!(!metadata.artifacts_verified);
    assert!(metadata.rows.is_empty());
    assert!(
        GentleEngine::inspect_observed_transcript_structures_from_sam(
            &binding.report_path,
            Some(&binding),
            TABLES[0],
            0,
            true
        )
        .is_err()
    );
    assert!(
        GentleEngine::inspect_observed_transcript_structures_from_sam(
            &binding.report_path,
            Some(&binding),
            TABLES[0],
            1,
            false
        )
        .is_err()
    );
}

#[test]
fn observed_sam_forged_authority_and_mismatched_bindings_are_rejected() {
    let (_root, r) = fixture();
    let mut report = run(&r);
    let binding = GentleEngine::observed_sam_result_binding(&r, &report).unwrap();
    let mut wrong = binding.clone();
    wrong.request_sha256 = format!("sha256:{}", "f".repeat(64));
    assert!(
        GentleEngine::inspect_observed_transcript_structures_from_sam(
            &binding.report_path,
            Some(&wrong),
            TABLES[0],
            0,
            false
        )
        .is_err()
    );
    report.coordinate_authority = ObservedSamCoordinateAuthority::Verified;
    report.content_sha256.clear();
    report.content_sha256 = digest(&report).unwrap();
    fs::write(
        &binding.report_path,
        serde_json::to_vec_pretty(&report).unwrap(),
    )
    .unwrap();
    assert_eq!(
        GentleEngine::inspect_observed_transcript_structures_from_sam(
            &binding.report_path,
            None,
            TABLES[0],
            0,
            false
        )
        .unwrap_err()
        .message,
        "observed_sam_report_binding"
    );
}

fn prepared(
    r: &mut Request,
    root: &Path,
    refused: bool,
) -> ObservedTranscriptAnnotationPreparation {
    let mut prep = GentleEngine::observed_annotation_test_fixture(root);
    prep.reference_id = r.reference_id.clone();
    prep.reference = r.reference.clone();
    prep.selector = r.reference_gene.clone();
    prep.dictionary = r.dictionary.clone();
    let next_start = if refused { 111 } else { 121 };
    let raw = format!(
        "ctg_a\ttest\ttranscript\t101\t140\t.\t+\t.\tgene_id \"synthetic_gene\"; transcript_id \"tx.7\"; transcript_version \"7\";\n\
ctg_a\ttest\texon\t101\t110\t.\t+\t.\tgene_id \"synthetic_gene\"; transcript_id \"tx.7\"; transcript_version \"7\";\n\
ctg_a\ttest\texon\t{next_start}\t140\t.\t+\t.\tgene_id \"synthetic_gene\"; transcript_id \"tx.7\"; transcript_version \"7\";\n"
    );
    fs::write(&prep.gtf.path, raw).unwrap();
    prep.gtf = bound(Path::new(&prep.gtf.path));
    let result = GentleEngine::new()
        .prepare_observed_transcript_annotation(prep, &mut |_| true)
        .unwrap();
    let annotation = result
        .report
        .annotation
        .as_ref()
        .map(|a| PrimerPayloadFile {
            path: Path::new(&result.report.request.output_dir)
                .join(&a.path)
                .display()
                .to_string(),
            bytes: a.bytes,
            sha256: a.sha256.clone(),
        });
    r.annotation_sources = vec![ObservedSamAnnotationBinding::Prepared {
        selector: r.reference_gene.clone(),
        report: result.report_file.clone(),
        annotation,
    }];
    result
}

#[test]
fn observed_sam_joins_prepared_selectors_exactly_without_version_or_symbol_guessing() {
    let (root, mut r) = fixture();
    let prep = prepared(&mut r, root.path(), false);
    let mut mismatch = r.clone();
    mismatch.reference_gene.gene_id.push_str(".1");
    if let ObservedSamAnnotationBinding::Prepared { selector, .. } =
        &mut mismatch.annotation_sources[0]
    {
        *selector = mismatch.reference_gene.clone();
    }
    assert_eq!(
        GentleEngine::new()
            .summarize_observed_transcript_structures_from_sam(mismatch, &mut |_| true)
            .unwrap_err()
            .message,
        "observed_sam_selector_binding_mismatch"
    );
    assert!(!Path::new(&r.output_dir).exists());
    // A portable copy must not depend on the old preparation output/GTF paths.
    let copied_report = root.path().join("copied receipt.json");
    fs::copy(&prep.report_file.path, &copied_report).unwrap();
    if let ObservedSamAnnotationBinding::Prepared {
        report,
        annotation: Some(annotation),
        ..
    } = &mut r.annotation_sources[0]
    {
        let copied = root.path().join("copied annotation.json");
        fs::copy(&annotation.path, &copied).unwrap();
        fs::remove_file(&annotation.path).unwrap();
        *annotation = bound(&copied);
        *report = bound(&copied_report);
    } else {
        panic!("completed preparation must include annotation");
    }
    fs::remove_file(&prep.report_file.path).unwrap();
    fs::remove_dir(Path::new(&prep.report.request.output_dir).join("artifacts")).unwrap();
    fs::remove_dir(&prep.report.request.output_dir).unwrap();
    fs::remove_file(&prep.report.request.gtf.path).unwrap();
    let report = run(&r);
    let saved = data(&report);
    assert_eq!(
        saved.annotations[0].availability,
        ObservedSamAnnotationAvailability::Available
    );
    assert_eq!(
        saved.annotations[0]
            .preparation
            .as_ref()
            .unwrap()
            .analysis
            .as_ref()
            .unwrap()
            .transcript_versions["tx.7"],
        ObservedAnnotationVersionState::Declared("7".into())
    );
    assert_eq!(
        report.gene_resolutions[0]
            .comparison
            .as_ref()
            .unwrap()
            .resolution,
        kernel::GeneResolution::Resolved
    );
}

#[test]
fn observed_sam_refused_preparation_remains_typed_unavailable_not_empty_or_missing_geometry() {
    let (root, mut r) = fixture();
    prepared(&mut r, root.path(), true);
    let report = run(&r);
    let saved = data(&report);
    assert_eq!(
        saved.annotations[0].availability,
        ObservedSamAnnotationAvailability::Refused
    );
    assert_eq!(
        report.gene_resolutions[0].availability,
        ObservedSamAnnotationAvailability::Refused
    );
    assert!(report.gene_resolutions[0].comparison.is_none());
    assert!(
        saved
            .analysis
            .as_ref()
            .unwrap()
            .annotation_bindings
            .is_empty()
    );
    assert!(!saved.analysis.as_ref().unwrap().observations.is_empty());
}

#[test]
fn observed_sam_cancellation_is_request_bound_and_does_not_require_original_inputs() {
    let (root, mut r) = fixture();
    r.cancel_path = Some(root.path().join("cancel.flag").display().to_string());
    for file in input_files(&r) {
        fs::remove_file(&file.path).unwrap();
    }
    GentleEngine::cancel_observed_transcript_structures_from_sam(&r).unwrap();
    let raw = fs::read(r.cancel_path.as_ref().unwrap()).unwrap();
    GentleEngine::cancel_observed_transcript_structures_from_sam(&r).unwrap();
    assert_eq!(fs::read(r.cancel_path.as_ref().unwrap()).unwrap(), raw);
    let report = run(&r);
    assert_eq!(report.status, ObservedSamExecutionStatus::Cancelled);
    assert!(report.counts.is_none());
    assert!(data(&report).analysis.is_none());
    let mut wrong = r.clone();
    wrong.required_flank_bases += 1;
    assert!(GentleEngine::cancel_observed_transcript_structures_from_sam(&wrong).is_err());
    assert_eq!(fs::read(r.cancel_path.as_ref().unwrap()).unwrap(), raw);
}

#[test]
fn observed_sam_control_error_is_distinct_and_late_cancellation_cannot_retract_publication() {
    let (root, mut r) = fixture();
    r.cancel_path = Some(root.path().join("cancel.flag").display().to_string());
    let marker = r.cancel_path.clone().unwrap();
    let report = GentleEngine::new()
        .summarize_observed_transcript_structures_from_sam(r.clone(), &mut |_| {
            fs::write(&marker, "wrong identity").unwrap();
            true
        })
        .unwrap();
    assert_eq!(report.status, ObservedSamExecutionStatus::ControlError);
    assert!(report.control_error.is_some());
    assert!(report.counts.is_none());
    let (root, mut r) = fixture();
    r.cancel_path = Some(root.path().join("cancel.flag").display().to_string());
    let copy = r.clone();
    let report = GentleEngine::new()
        .summarize_observed_transcript_structures_from_sam(r, &mut |p| {
            if matches!(
                p,
                OperationProgress::Workflow {
                    completed: 5,
                    total: 5
                }
            ) {
                GentleEngine::cancel_observed_transcript_structures_from_sam(&copy).unwrap();
                false
            } else {
                true
            }
        })
        .unwrap();
    assert_eq!(report.status, ObservedSamExecutionStatus::Completed);
    assert!(
        Path::new(&report.resolved_output_dir)
            .join("report.json")
            .exists()
    );
}

#[test]
fn observed_sam_output_and_control_collisions_are_refused_without_source_mutation() {
    let (_root, mut r) = fixture();
    r.cancel_path = Some(r.sam.path.clone());
    let raw = fs::read(&r.sam.path).unwrap();
    assert_eq!(
        GentleEngine::cancel_observed_transcript_structures_from_sam(&r)
            .unwrap_err()
            .message,
        "observed_sam_control_collision"
    );
    assert_eq!(fs::read(&r.sam.path).unwrap(), raw);
    r.cancel_path = None;
    r.budgets.max_output_bytes = 10;
    assert!(
        GentleEngine::new()
            .summarize_observed_transcript_structures_from_sam(r.clone(), &mut |_| true)
            .is_err()
    );
    assert!(!Path::new(&r.output_dir).exists());
}

#[test]
fn observed_sam_forged_agent_flags_cannot_bypass_write_approval() {
    use crate::agent_planner::{
        AgentPlanCandidate, AgentPlanCandidateKind, AgentPlanExecutionMode, AgentPlanResult,
    };
    use crate::engine::{Operation, lifecycle_slice3_tests};
    let (root, mut r) = fixture();
    let marker = root.path().join("approval marker");
    r.cancel_path = Some(marker.display().to_string());
    for op in [
        Operation::SummarizeObservedTranscriptStructuresFromSam {
            request: Box::new(r.clone()),
        },
        Operation::CancelObservedTranscriptStructuresFromSam {
            request: Box::new(r.clone()),
        },
    ] {
        let plan = AgentPlanResult {
            schema: "gentle.agent_plan_result.v1".into(),
            assistant_message: String::new(),
            questions: vec![],
            candidates: vec![AgentPlanCandidate {
                candidate_id: "sam".into(),
                title: "Synthetic approval regression".into(),
                kind: AgentPlanCandidateKind::Op,
                execution_mode: AgentPlanExecutionMode::Auto,
                operation: Some(serde_json::to_value(op).unwrap()),
                mutating: false,
                requires_confirmation: false,
                ..Default::default()
            }],
        };
        let mut engine = lifecycle_slice3_tests::live_with_redo();
        let state = encode(engine.state()).unwrap();
        let before = engine.history_summary();
        let result = crate::agent_execution::execute_agent_plan_candidate_with_feedback(
            &mut engine,
            &plan,
            "sam",
            false,
            &Default::default(),
        );
        assert!(result.result.unwrap_err().contains("explicit --confirm"));
        assert_eq!(encode(engine.state()).unwrap(), state);
        let after = engine.history_summary();
        assert_eq!(
            (
                before.undo_count,
                before.redo_count,
                before.operation_log_count
            ),
            (
                after.undo_count,
                after.redo_count,
                after.operation_log_count
            )
        );
        assert!(!Path::new(&r.output_dir).exists());
        assert!(!marker.exists());
    }
}

#[test]
fn observed_sam_exact_output_row_and_total_byte_limits_are_atomic() {
    let (root, r) = fixture();
    let baseline = run(&r);
    let largest_row = TABLES
        .iter()
        .flat_map(|name| {
            let mut reader = csv::ReaderBuilder::new()
                .delimiter(b'\t')
                .from_path(
                    Path::new(&baseline.resolved_output_dir)
                        .join("artifacts")
                        .join(name),
                )
                .unwrap();
            reader
                .records()
                .map(|row| row.unwrap()[0].len())
                .collect::<Vec<_>>()
        })
        .max()
        .unwrap();
    assert!(largest_row > 1);
    let mut exact = r.clone();
    exact.budgets.max_output_row_bytes = largest_row;
    exact.output_dir = root.path().join("exact rows").display().to_string();
    run(&exact);
    exact.budgets.max_output_row_bytes -= 1;
    exact.output_dir = root.path().join("refused rows").display().to_string();
    assert!(
        GentleEngine::new()
            .summarize_observed_transcript_structures_from_sam(exact.clone(), &mut |_| true)
            .is_err()
    );
    assert!(!Path::new(&exact.output_dir).exists());

    // Only the budget integer's digit count changes this self-describing report's size.
    // Keep output-path byte lengths equal, rather than guessing a platform's temp prefix.
    let total = baseline.artifacts.iter().map(|a| a.bytes).sum::<u64>()
        + fs::metadata(Path::new(&baseline.resolved_output_dir).join("report.json"))
            .unwrap()
            .len();
    let old_digits = r.budgets.max_output_bytes.to_string().len() as u64;
    let bound = |margin| {
        let mut n = total;
        loop {
            let next = total - old_digits + n.to_string().len() as u64 - margin;
            if next == n {
                break n;
            }
            n = next;
        }
    };
    let name_len = Path::new(&r.output_dir)
        .file_name()
        .unwrap()
        .to_str()
        .unwrap()
        .len();
    let mut exact = r.clone();
    exact.output_dir = root.path().join("e".repeat(name_len)).display().to_string();
    exact.budgets.max_output_bytes = bound(0);
    let report = run(&exact);
    assert_eq!(
        report.artifacts.iter().map(|a| a.bytes).sum::<u64>()
            + fs::metadata(Path::new(&report.resolved_output_dir).join("report.json"))
                .unwrap()
                .len(),
        exact.budgets.max_output_bytes
    );
    exact.output_dir = root.path().join("f".repeat(name_len)).display().to_string();
    exact.budgets.max_output_bytes = bound(1);
    assert!(
        GentleEngine::new()
            .summarize_observed_transcript_structures_from_sam(exact.clone(), &mut |_| true)
            .is_err()
    );
    assert!(!Path::new(&exact.output_dir).exists());
}

#[test]
fn observed_sam_publisher_preserves_replaced_artifacts_and_reports_orphans() {
    let (root, r) = fixture();
    let report = run(&r);
    let mut copy = report.clone();
    copy.resolved_output_dir = root.path().join("fault output").display().to_string();
    let stage = tempfile::tempdir_in(root.path()).unwrap();
    for artifact in &copy.artifacts {
        fs::copy(
            Path::new(&report.resolved_output_dir)
                .join("artifacts")
                .join(&artifact.path),
            stage.path().join(&artifact.path),
        )
        .unwrap();
    }
    let moved = root.path().join("moved originals");
    let mut rename_blocked = false;
    let result = publication::publish(stage, &copy, &mut |phase| {
        if phase == publication::Phase::ArtifactsVisible {
            let path = Path::new(&copy.resolved_output_dir).join("artifacts");
            if let Err(error) = fs::rename(&path, &moved) {
                // The locked-directory refusal is itself the observable Windows guard.
                assert!(
                    cfg!(windows) && error.kind() == std::io::ErrorKind::PermissionDenied,
                    "unexpected artifact rename failure: {error}"
                );
                rename_blocked = true;
                return Err(bad("synthetic blocked replacement fault"));
            }
            fs::create_dir(&path).unwrap();
            fs::write(path.join("foreign"), "preserve").unwrap();
        }
        Ok(())
    })
    .unwrap_err();
    if rename_blocked {
        assert!(
            result
                .message
                .contains("synthetic blocked replacement fault")
        );
        assert!(!moved.exists());
        assert!(
            !Path::new(&copy.resolved_output_dir).exists(),
            "owned output must be cleaned after refused replacement"
        );
        return;
    }
    assert!(result.message.contains("observed_sam_orphaned_output"));
    assert_eq!(
        fs::read_to_string(Path::new(&copy.resolved_output_dir).join("artifacts/foreign")).unwrap(),
        "preserve"
    );
    assert!(
        !Path::new(&copy.resolved_output_dir)
            .join("report.json")
            .exists()
    );
    assert_eq!(fs::read_dir(moved).unwrap().count(), 7);
}

#[test]
fn observed_sam_publisher_preserves_unowned_artifact_after_fault() {
    let (root, request) = fixture();
    let report = run(&request);
    let mut copy = report.clone();
    copy.resolved_output_dir = root
        .path()
        .join("unowned artifact fault")
        .display()
        .to_string();
    let stage = tempfile::tempdir_in(root.path()).unwrap();
    for artifact in &copy.artifacts {
        fs::copy(
            Path::new(&report.resolved_output_dir)
                .join("artifacts")
                .join(&artifact.path),
            stage.path().join(&artifact.path),
        )
        .unwrap();
    }
    let out = Path::new(&copy.resolved_output_dir);
    let error = publication::publish(stage, &copy, &mut |phase| {
        if phase == publication::Phase::ArtifactsVisible {
            fs::write(out.join("artifacts/foreign"), b"unowned evidence").unwrap();
            return Err(bad("synthetic unowned artifact fault"));
        }
        Ok(())
    })
    .unwrap_err();
    assert!(error.message.contains("observed_sam_orphaned_output"));
    assert_eq!(
        fs::read(out.join("artifacts/foreign")).unwrap(),
        b"unowned evidence"
    );
    assert!(!out.join("report.json").exists());
    assert_eq!(fs::read_dir(out.join("artifacts")).unwrap().count(), 1);
}
