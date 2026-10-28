//! Inline synthetic genomic Phase 2 producer/consumer tests; no real reads, genes or aligner.

use super::tests::{PAYLOAD, bound_file, limits, prepare, source};
use super::*;
use crate::observed_transcript_structures as kernel;
use serde_json::json;

#[path = "observed_preservation_tests.rs"]
mod preservation;

fn fixture(root: &Path) -> (GentleEngine, SummarizeObservedTranscriptStructuresRequest) {
    fixture_mode(root, true)
}

fn fixture_mode(
    root: &Path,
    source_continue: bool,
) -> (GentleEngine, SummarizeObservedTranscriptStructuresRequest) {
    let (engine, trace) = source(root, true, PAYLOAD);
    let prepared = engine
        .prepare_primer_payload_alignment_queries(prepare(trace.clone(), root), &mut |_| true)
        .unwrap();
    let mut req = super::alignment_tests::request(root, trace, &prepared);
    let n = prepared.queries[0].query_length;
    let half = n / 2;
    let reference = &mut req.references[0];
    reference.scope = PrimerPayloadReferenceScope::Genome;
    let mut annotation: PrimerPayloadAnnotation =
        read_json(&reference.annotation, &limits()).unwrap();
    for (i, t) in annotation.transcripts.iter_mut().enumerate() {
        t.genomic_contig = Some(t.target_name.clone());
        t.genomic_strand = Some(if i == 0 { '+' } else { '-' });
        t.genomic_exons = vec![(0, half), (half + 5, n + 5)];
    }
    fs::write(
        &reference.annotation.path,
        serde_json::to_vec(&annotation).unwrap(),
    )
    .unwrap();
    reference.annotation = bound_file(Path::new(&reference.annotation.path));
    // Invented spliced alignment: exact coordinate consumption, not claimed base/search accuracy.
    let paf=prepared.queries.iter().map(|q|format!("{}\t{n}\t0\t{n}\t+\ttoy_foreign\t{}\t0\t{}\t{n}\t{n}\t60\tcg:Z:{half}M5N{}M\ttp:A:P\n",q.query_id,n+10,n+5,n-half)).collect::<String>();
    fs::write(&req.batches[0].paf.path, paf).unwrap();
    req.batches[0].paf = bound_file(Path::new(&req.batches[0].paf.path));
    let reference = req.references[0].clone();
    let report = engine
        .summarize_primer_payload_alignments(req, &mut |_| source_continue)
        .unwrap();
    let binding =
        GentleEngine::primer_payload_alignment_result_binding(&report.request, &report).unwrap();
    let request = SummarizeObservedTranscriptStructuresRequest {
        schema: "gentle.observed_transcript_structures_request.v1".into(),
        source_report: bound_file(Path::new(&binding.report_path)),
        source: binding,
        batch_id: "synthetic_batch".into(),
        reference_gene: kernel::ReferenceGene {
            source_id: "synthetic_annotation".into(),
            gene_id: "foreign_gene".into(),
        },
        comparison_genes: Vec::new(),
        annotation_sources: vec![ObservedTranscriptAnnotation {
            source_id: "synthetic_annotation".into(),
            provider: Some("synthetic".into()),
            reference,
        }],
        required_flank_bases: 2,
        budgets: limits(),
        kernel_budgets: kernel::Budgets {
            max_queries: 100,
            max_records: 100,
            max_segments_per_query: 100,
            max_operations: 10000,
            max_skips_per_chain: 100,
            max_chain_families: 100,
            max_alternatives_per_record: 100,
            max_annotation_rows: 100,
            max_model_comparisons_per_family: 100,
            max_total_work: 1_000_000,
        },
        output_dir: root.join("observed gene's evidence").display().to_string(),
        cancel_path: None,
    };
    (engine, request)
}

impl GentleEngine {
    /// Shared inline synthetic producer for adapter tests, absent from production builds.
    pub(crate) fn observed_transcripts_test_fixture(
        root: &Path,
    ) -> (Self, SummarizeObservedTranscriptStructuresRequest) {
        fixture(root)
    }
}

#[test]
fn observed_transcripts_raw_lifecycle_preserves_unreconciled_content_and_history() {
    let td = tempfile::tempdir().unwrap();
    let (_, request) = fixture(td.path());
    let op = Operation::SummarizeObservedTranscriptStructures {
        request: Box::new(request),
    };
    let initial = crate::engine::lifecycle_slice3_tests::state();
    let mut engine =
        GentleEngine::from_state_for(initial.clone(), LifecycleRequest::Operation(&op));
    assert_eq!(encode(&initial).unwrap(), encode(engine.state()).unwrap());
    let mutation = engine.mutation_revision();
    let structural = engine.structural_revision();
    let history = engine.history_summary();
    let result = engine.apply(op).unwrap();
    assert!(result.observed_transcript_structures.is_some());
    assert_eq!(encode(&initial).unwrap(), encode(engine.state()).unwrap());
    assert_eq!(
        (mutation, structural),
        (engine.mutation_revision(), engine.structural_revision())
    );
    let after = engine.history_summary();
    assert_eq!(
        (history.undo_count, history.redo_count),
        (after.undo_count, after.redo_count)
    );
    assert_eq!(after.operation_log_count, history.operation_log_count + 1);
}

#[test]
fn observed_transcripts_incomplete_batch_counts_are_absent_not_zero() {
    let td = tempfile::tempdir().unwrap();
    let (engine, r) = fixture_mode(td.path(), false);
    let report = engine
        .summarize_observed_transcript_structures(r, &mut |_| true)
        .unwrap();
    assert_eq!(report.status, "unavailable");
    assert_eq!(report.source_status, "cancelled");
    assert!(report.counts.is_none());
    let saved = data(&report);
    assert!(saved.analysis.is_none());
    assert!(!saved.unavailable.is_empty());
}

#[test]
fn observed_transcripts_late_cancel_and_control_fault_do_not_become_success() {
    for fault in [false, true] {
        let td = tempfile::tempdir().unwrap();
        let (engine, mut r) = fixture(td.path());
        let marker = td.path().join("observed marker");
        r.cancel_path = Some(marker.display().to_string());
        let report = engine
            .summarize_observed_transcript_structures(r, &mut |p| {
                if fault {
                    fs::write(&marker, "wrong identity").unwrap();
                    true
                } else {
                    !matches!(p, OperationProgress::Workflow { completed: 1, .. })
                }
            })
            .unwrap();
        assert_eq!(
            report.status,
            if fault { "control_error" } else { "cancelled" }
        );
        assert_eq!(report.counts.is_none(), fault);
        assert_eq!(report.control_error.is_some(), fault);
    }
}

#[test]
fn observed_transcripts_resource_and_annotation_errors_publish_no_report() {
    for kind in ["annotation", "input", "output", "kernel"] {
        let td = tempfile::tempdir().unwrap();
        let (engine, mut r) = fixture(td.path());
        match kind {
            "annotation" => {
                r.annotation_sources[0].reference.annotation.sha256 =
                    format!("sha256:{}", "0".repeat(64))
            }
            "input" => r.budgets.max_input_bytes = 1,
            "output" => r.budgets.max_output_bytes = 1,
            _ => r.kernel_budgets.max_total_work = 1,
        }
        let output = PathBuf::from(&r.output_dir);
        assert!(
            engine
                .summarize_observed_transcript_structures(r, &mut |_| true)
                .is_err(),
            "{kind}"
        );
        assert!(!output.join("report.json").exists());
    }
}

#[test]
fn observed_transcripts_rehashed_duplicate_missing_or_foreign_rows_are_refused() {
    for kind in ["duplicate", "missing", "foreign", "geometry"] {
        let td = tempfile::tempdir().unwrap();
        let (engine, mut r) = fixture(td.path());
        let mut source: PrimerPayloadAlignmentReport =
            read_json(&r.source_report, &limits()).unwrap();
        let path = Path::new(&source.request.output_dir).join("artifacts/read_anchor_evidence.tsv");
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(&path)
            .unwrap();
        let row = reader.records().next().unwrap().unwrap();
        let mut e: PrimerPayloadReadEvidence = serde_json::from_str(&row[0]).unwrap();
        drop(reader);
        match kind {
            "foreign" => e.read_id = "foreign_read".into(),
            "geometry" => e.segments[0].target_interval.1 += 1,
            _ => {}
        }
        let count = if kind == "missing" {
            0
        } else if kind == "duplicate" {
            2
        } else {
            1
        };
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_path(&path)
            .unwrap();
        writer.write_record(["evidence_json"]).unwrap();
        for _ in 0..count {
            writer
                .write_record([serde_json::to_string(&e).unwrap()])
                .unwrap();
        }
        writer.flush().unwrap();
        drop(writer);
        let file = bound_file(&path);
        let a = source
            .artifacts
            .iter_mut()
            .find(|a| a.path == "read_anchor_evidence.tsv")
            .unwrap();
        a.sha256 = file.sha256;
        a.bytes = file.bytes;
        source.evidence_rows = count;
        source.content_sha256.clear();
        source.content_sha256 = hash(&source).unwrap();
        fs::write(&r.source_report.path, serde_json::to_vec(&source).unwrap()).unwrap();
        r.source_report = bound_file(Path::new(&r.source_report.path));
        r.source = GentleEngine::primer_payload_alignment_result_binding(&source.request, &source)
            .unwrap();
        assert!(
            engine
                .summarize_observed_transcript_structures(r, &mut |_| true)
                .is_err(),
            "{kind}"
        );
    }
}

fn data(report: &ObservedTranscriptStructuresReport) -> ObservedTranscriptStructuresData {
    serde_json::from_slice(
        &fs::read(Path::new(&report.request.output_dir).join("artifacts/structures.json")).unwrap(),
    )
    .unwrap()
}

#[test]
fn observed_transcripts_unknown_target_length_withholds_the_whole_view() {
    let td = tempfile::tempdir().unwrap();
    let (engine, mut r) = fixture(td.path());
    let original = &r.annotation_sources[0].reference.annotation;
    let mut annotation: PrimerPayloadAnnotation = read_json(original, &limits()).unwrap();
    annotation
        .transcripts
        .retain(|t| t.target_name == "toy_expected");
    let path = td.path().join("independent annotation.json");
    fs::write(&path, serde_json::to_vec(&annotation).unwrap()).unwrap();
    r.annotation_sources[0].reference.annotation = bound_file(&path);
    r.reference_gene.gene_id = "expected_gene".into();
    let report = engine
        .summarize_observed_transcript_structures(r, &mut |_| true)
        .unwrap();
    let data = data(&report);
    let analysis = data.analysis.unwrap();
    assert!(analysis.observations.is_empty());
    assert!(analysis.families.is_empty());
    assert_eq!(analysis.records.len(), 1);
    assert!(!analysis.records[0].all_query_views_assessed);
    assert!(analysis.records[0].definite_chains.is_empty());
    assert_eq!(data.unavailable[0].reason, "target_length_unavailable");
}

#[test]
fn observed_transcripts_agent_feedback_does_not_change_phase2_outcome_rules() {
    use crate::agent_feedback::AgentExecutionStatus as Status;
    for (outcome, status) in [
        ("completed", Status::Completed),
        ("cancelled", Status::Cancelled),
        ("control_error", Status::Failed),
        ("unavailable", Status::Partial),
    ] {
        let compact = json!({"schema":"gentle.command_observed_structures_result.v1","observed_structures_result":{"schema":"gentle.observed_transcript_structures_result_binding.v1","outcome":outcome}});
        let raw = json!({"result":{"observed_transcript_structures":{"schema":"gentle.observed_transcript_structures_report.v1","status":outcome}}});
        assert_eq!(Status::from_shell_output(&compact), status);
        assert_eq!(Status::from_shell_output(&raw), status);
    }
    assert_eq!(
        Status::from_shell_output(
            &json!({"result":{"primer_payload_alignments":{"schema":"gentle.primer_payload_alignment_report.v1","status":"unavailable"}}})
        ),
        Status::Failed
    );
}

#[test]
fn observed_transcripts_bound_intake_keeps_phase2_bytes_and_configurable_gene() {
    let td = tempfile::tempdir().unwrap();
    let (engine, r) = fixture(td.path());
    let source: PrimerPayloadAlignmentReport = read_json(&r.source_report, &limits()).unwrap();
    let saved = std::iter::once(PathBuf::from(&r.source_report.path))
        .chain(source.artifacts.iter().map(|a| {
            Path::new(&source.request.output_dir)
                .join("artifacts")
                .join(&a.path)
        }))
        .chain([
            PathBuf::from(&source.request.queries.path),
            PathBuf::from(&source.request.trace.report_path),
        ])
        .map(|p| {
            let bytes = fs::read(&p).unwrap();
            (p, bytes)
        })
        .collect::<Vec<_>>();
    for p in [
        source.request.references[0].fasta.path.clone(),
        source.request.batches[0].paf.path.clone(),
        td.path()
            .join("invented input's reads.fastq")
            .display()
            .to_string(),
    ] {
        fs::remove_file(p).unwrap();
    }
    let state = encode(engine.state()).unwrap();
    let journal = engine.journal_len();
    let a = engine
        .summarize_observed_transcript_structures(r.clone(), &mut |_| true)
        .unwrap();
    assert_eq!(a.status, "completed");
    assert_eq!(a.counts.as_ref().unwrap().chain_families, 1);
    assert_eq!(state, encode(engine.state()).unwrap());
    assert_eq!(journal, engine.journal_len());
    let a_data = data(&a);
    let s = a_data.analysis.unwrap();
    assert_eq!(s.reference_gene.gene_id, "foreign_gene");
    assert_eq!(s.observations[0].comparisons[0].models[0].model_strand, '-');
    let chain = s.families[0].chain_id.clone();
    let mut other = r.clone();
    other.reference_gene.gene_id = "expected_gene".into();
    other.output_dir = td.path().join("other gene").display().to_string();
    let b = engine
        .summarize_observed_transcript_structures(other, &mut |_| true)
        .unwrap();
    assert_ne!(a.request_sha256, b.request_sha256);
    assert_eq!(data(&b).analysis.unwrap().families[0].chain_id, chain);
    let mut absent = r.clone();
    absent.reference_gene.gene_id = "absent_gene".into();
    absent.output_dir = td.path().join("missing gene").display().to_string();
    let missing = engine
        .summarize_observed_transcript_structures(absent, &mut |_| true)
        .unwrap();
    assert_eq!(
        data(&missing).analysis.unwrap().gene_resolutions[0].resolution,
        kernel::GeneResolution::MissingGene
    );
    let binding = GentleEngine::observed_transcript_structures_result_binding(&r, &a).unwrap();
    let inspection = GentleEngine::inspect_observed_transcript_structures(
        &binding.report_path,
        Some(&binding),
        "read_structures.tsv",
        1,
        true,
    )
    .unwrap();
    assert_eq!(inspection.rows.len(), 1);
    assert!(inspection.artifacts_verified);
    for (p, bytes) in saved {
        assert_eq!(fs::read(p).unwrap(), bytes);
    }
    let mut outer = GentleEngine::new()
        .apply(Operation::InspectPrimerPayload {
            path: r.source_report.path.clone(),
            expected_result: Some(Box::new(r.source.clone())),
            rows: 0,
            verify_artifacts: false,
        })
        .unwrap();
    outer.primer_payload_alignments = Some(Box::new(source));
    let serialized = serde_json::to_value(outer).unwrap();
    assert!(serialized.get("observed_transcript_structures").is_none());
    assert!(
        serialized
            .get("observed_transcript_structures_inspection")
            .is_none()
    );
    assert!(
        engine
            .summarize_observed_transcript_structures(r, &mut |_| true)
            .is_err()
    );
}

#[test]
fn observed_transcripts_source_tampering_and_request_bound_cancel_fail_closed() {
    let td = tempfile::tempdir().unwrap();
    let (engine, mut r) = fixture(td.path());
    let mut wrong = r.clone();
    wrong.source.report_content_sha256 = "sha256:wrong".into();
    assert!(
        engine
            .summarize_observed_transcript_structures(wrong, &mut |_| true)
            .is_err()
    );
    assert!(!Path::new(&r.output_dir).exists());
    let marker = td.path().join("cancel marker");
    r.cancel_path = Some(marker.display().to_string());
    GentleEngine::cancel_observed_transcript_structures(&r).unwrap();
    let mut wrong = r.clone();
    wrong.reference_gene.gene_id = "expected_gene".into();
    assert!(
        engine
            .summarize_observed_transcript_structures(wrong, &mut |_| true)
            .is_err()
    );
    let result = engine
        .summarize_observed_transcript_structures(r, &mut |_| true)
        .unwrap();
    assert_eq!(result.status, "cancelled");
    assert!(result.counts.is_none());
    assert!(data(&result).analysis.is_none());
}

#[test]
fn observed_transcripts_hashes_cover_full_evidence_and_retained_rows() {
    for crlf in [false, true] {
        let td = tempfile::tempdir().unwrap();
        let (engine, mut r) = fixture(td.path());
        let mut report: PrimerPayloadAlignmentReport =
            read_json(&r.source_report, &limits()).unwrap();
        let a = report
            .artifacts
            .iter_mut()
            .find(|a| a.path == "read_anchor_evidence.tsv")
            .unwrap();
        let path = Path::new(&report.request.output_dir)
            .join("artifacts")
            .join(&a.path);
        let original = fs::read_to_string(&path).unwrap();
        let text = if crlf {
            original.replace('\n', "\r\n")
        } else {
            original
        };
        fs::write(&path, text.as_bytes()).unwrap();
        a.sha256 = sha256_prefixed_bytes(text.as_bytes());
        a.bytes = text.len() as u64;
        report.content_sha256.clear();
        report.content_sha256 = hash(&report).unwrap();
        fs::write(
            &r.source_report.path,
            serde_json::to_vec_pretty(&report).unwrap(),
        )
        .unwrap();
        r.source_report = bound_file(Path::new(&r.source_report.path));
        r.source = GentleEngine::primer_payload_alignment_result_binding(&report.request, &report)
            .unwrap();
        let result = engine
            .summarize_observed_transcript_structures(r.clone(), &mut |_| true)
            .unwrap();
        let bind =
            GentleEngine::observed_transcript_structures_result_binding(&r, &result).unwrap();
        let path = Path::new(&r.output_dir).join("artifacts/read_structures.tsv");
        let mut bytes = fs::read(&path).unwrap();
        bytes.push(b' ');
        fs::write(path, bytes).unwrap();
        assert!(
            GentleEngine::inspect_observed_transcript_structures(
                &bind.report_path,
                Some(&bind),
                "read_structures.tsv",
                1,
                false
            )
            .is_err()
        );
    }
}

#[test]
fn observed_transcripts_shared_parser_discovery_and_approval_are_explicit() {
    let td = tempfile::tempdir().unwrap();
    let (mut engine, r) = fixture(td.path());
    let path = td.path().join("request's file.json");
    fs::write(&path, serde_json::to_vec(&r).unwrap()).unwrap();
    let line = format!(
        "rna-reads transcript-structures {}",
        crate::engine_shell::quote_shell_arg(&format!("@{}", path.display()))
    );
    let command = crate::engine_shell::parse_shell_line(&line).unwrap();
    assert!(command.is_state_mutating());
    let result = crate::engine_shell::execute_shell_command(&mut engine, &command).unwrap();
    assert!(!result.state_changed);
    assert_eq!(
        result.output["result"]["observed_transcript_structures"]["status"],
        "completed"
    );
    let capability = crate::engine_shell::execute_shell_command(
        &mut engine,
        &crate::engine_shell::parse_shell_line("introspect capabilities --kind operation").unwrap(),
    )
    .unwrap();
    for id in [
        "SummarizeObservedTranscriptStructures",
        "InspectObservedTranscriptStructures",
        "CancelObservedTranscriptStructures",
    ] {
        let row = capability.output["capabilities"]
            .as_array()
            .unwrap()
            .iter()
            .find(|c| c["id"] == id)
            .unwrap();
        let op: Operation = serde_json::from_value(row["operation_template"].clone()).unwrap();
        assert!(crate::engine::LifecycleRequest::Operation(&op).preserves_project_content());
        assert_eq!(
            row["requires_confirmation"],
            json!(id != "InspectObservedTranscriptStructures")
        );
    }
}
