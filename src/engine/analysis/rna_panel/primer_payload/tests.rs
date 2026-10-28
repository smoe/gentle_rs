//! Inline invented oligos, FASTQ and annotation geometry; deterministic inputs, no assay evidence.

use super::*;

#[test]
fn primer_payload_annotation_preparation_preserves_legacy_catalog_touching_policy() {
    // Invented geometry, recreated in a private temp directory through the real
    // bound JSON reader. This is catalog admission, not a FASTA/alignment replay.
    let td = tempfile::tempdir().unwrap();
    let digest = format!("sha256:{}", "1".repeat(64));
    let annotation = PrimerPayloadAnnotation {
        schema: PRIMER_PAYLOAD_ANNOTATION_SCHEMA.into(),
        reference_id: "invented".into(),
        assembly: "invented".into(),
        release: "invented".into(),
        reference_fasta_sha256: digest.clone(),
        coordinate_convention: "zero_based_half_open".into(),
        transcripts: vec![PrimerPayloadTranscript {
            transcript_id: "invented_tx".into(),
            gene_ids: vec!["invented_gene".into()],
            target_name: "invented_contig".into(),
            target_length: 1000,
            genomic_contig: Some("invented_contig".into()),
            genomic_strand: Some('+'),
            genomic_exons: vec![(100, 110), (110, 120)],
        }],
    };
    let path = td.path().join("invented touching annotation.json");
    fs::write(&path, serde_json::to_vec_pretty(&annotation).unwrap()).unwrap();
    let r = PrimerPayloadReference {
        reference_id: "invented".into(),
        scope: PrimerPayloadReferenceScope::Genome,
        assembly: "invented".into(),
        release: "invented".into(),
        fasta: PrimerPayloadFile {
            path: "not opened by catalog".into(),
            sha256: digest,
            bytes: 0,
        },
        annotation: bound_file(&path),
    };
    let catalog = reference::catalog(&r, &limits()).unwrap();
    assert_eq!(
        catalog.annotation.transcripts[0].genomic_exons,
        [(100, 110), (110, 120)]
    );
}

pub(super) fn limits() -> PrimerPayloadBudgets {
    PrimerPayloadBudgets {
        max_input_bytes: 20_000_000,
        max_output_bytes: 20_000_000,
        max_records: 100,
        max_queries: 100,
        max_query_bases: 10_000,
        max_line_bytes: 100_000,
        max_alignments_per_query: 100,
        max_cigar_operations: 1000,
        max_annotation_rows: 100,
        max_reference_context_bases: 100_000,
        max_chain_comparisons: 100_000,
    }
}

pub(super) fn source(
    root: &Path,
    fasta: bool,
    sequence: &str,
) -> (GentleEngine, RnaPrimerTraceResultBinding) {
    source_reads(root, fasta, &[format!("AGTCCGATCA{sequence}")])
}

pub(super) fn source_reads(
    root: &Path,
    fasta: bool,
    reads: &[String],
) -> (GentleEngine, RnaPrimerTraceResultBinding) {
    let engine = GentleEngine::new();
    let mut panel: RnaPrimerPanel = serde_json::from_value(serde_json::json!({
        "schema":RNA_PRIMER_PANEL_SCHEMA,"content_sha256":"", "policy":{"max_full_edits":0,"max_core_edits_for_attribution":0},
        "primers":[{"primer_id":"invented_primer","assay_label":"Synthetic only",
            "sequence":"AGTCCGATCA", "components":[{"component_id":"core","role":"target_core","start":0,"end":10}],
            "core_orientation":"sense_5prime","expected_targets":[{"locus_id":"invented_expected_locus","transcript_feature_id":null}],"variant":"invented","source":"Inline Phase 2 fixture"}]
    })).unwrap();
    GentleEngine::seal_rna_primer_panel(&mut panel).unwrap();
    let bytes = reads
        .iter()
        .map(|read| format!("@reused_name\n{read}\n+\n{}\n", "I".repeat(read.len())))
        .collect::<String>();
    let input = root.join("invented input's reads.fastq");
    fs::write(&input, bytes.as_bytes()).unwrap();
    let request: TracePrimerReadsRequest = serde_json::from_value(serde_json::json!({
        "schema":RNA_PRIMER_TRACE_REQUEST_SCHEMA,"primer_panel":{"kind":"inline","panel":panel},
        "files":[{"sample_id":"toy_sample","barcode":"toy_barcode","priming":"declared_toy",
            "path":input,"declared_complete":true,"expected_sha256":sha256_prefixed_bytes(bytes.as_bytes())}],
        "export_payload_fasta":fasta,"output_dir":root.join("phase one")
    })).unwrap();
    let report = engine
        .trace_primer_reads(request.clone(), &mut |_| true)
        .unwrap();
    let binding = GentleEngine::primer_trace_result_binding(&request, &report).unwrap();
    (engine, binding)
}

pub(super) fn prepare(
    binding: RnaPrimerTraceResultBinding,
    root: &Path,
) -> PreparePrimerPayloadAlignmentQueriesRequest {
    PreparePrimerPayloadAlignmentQueriesRequest {
        schema: PRIMER_PAYLOAD_PREPARE_SCHEMA.into(),
        trace: binding,
        file_ids: vec!["file_0000".into()],
        cohort: PrimerPayloadCohort {
            cohort_id: "invented_cohort".into(),
            label: "Synthetic cohort".into(),
            all_verified_units: true,
            selection_provenance: "All verified units of this invented trace".into(),
        },
        source_mode: PrimerPayloadSourceMode::RetainedPrimaryFasta,
        anchor_selection: PrimerPayloadAnchorSelection::Primary,
        internal_primer_policy: PrimerPayloadInternalPrimerPolicy::Mask,
        min_payload_bases: 12,
        budgets: limits(),
        output_dir: root.join("fresh queries").display().to_string(),
        cancel_path: None,
    }
}

pub(super) const PAYLOAD: &str = "GATCCATGAGCTAGTCGACATCGTAGCTACGATGCTAGCATCG";

pub(super) fn bound_file(path: &Path) -> PrimerPayloadFile {
    let bytes = fs::read(path).unwrap();
    PrimerPayloadFile {
        path: path.display().to_string(),
        sha256: sha256_prefixed_bytes(&bytes),
        bytes: bytes.len() as u64,
    }
}

fn interval_hit(start: usize, end: usize, strand: PrimerReadStrand) -> RnaPrimerHit {
    RnaPrimerHit {
        primer_id: Some("abstract_interval_only".into()),
        anchor_class: "synthetic_geometry".into(),
        strand,
        start,
        end,
        edits: 0,
        core_start: Some(start),
        core_end: Some(end),
        core_edits: Some(0),
        terminal: false,
    }
}

#[test]
fn primer_payload_masks_exhaustive_positions_match_basewise_oracle() {
    // Abstract intervals, not biologically admissible oligos or another primer matcher.
    let length = 12usize;
    let mut cases = 0;
    for start in 0..length {
        for end in start + 1..=length {
            for other_start in 0..length {
                for other_end in other_start + 1..=length {
                    for strand in [
                        PrimerReadStrand::Forward,
                        PrimerReadStrand::ReverseComplement,
                    ] {
                        let anchor = interval_hit(start, end, strand);
                        let other = interval_hit(other_start, other_end, PrimerReadStrand::Forward);
                        let hits = vec![anchor.clone(), other];
                        let got = prepare::payload_geometry(length, &anchor, &hits).unwrap();
                        let positions = match strand {
                            PrimerReadStrand::Forward => (end..length).collect::<Vec<_>>(),
                            PrimerReadStrand::ReverseComplement => {
                                (0..start).rev().collect::<Vec<_>>()
                            }
                        };
                        let mut expected = Vec::<(usize, usize)>::new();
                        let mut usable = 0;
                        for (index, raw) in positions.iter().enumerate() {
                            if hits.iter().any(|h| (h.start..h.end).contains(raw)) {
                                if let Some(last) = expected.last_mut()
                                    && last.1 == index
                                {
                                    last.1 += 1;
                                } else {
                                    expected.push((index, index + 1));
                                }
                            } else {
                                usable += 1;
                            }
                        }
                        assert_eq!(
                            got.excluded, expected,
                            "{start}..{end}, {other_start}..{other_end}, {strand:?}"
                        );
                        assert_eq!(got.usable, usable);
                        assert_eq!(got.interval.1 - got.interval.0, positions.len());
                        cases += 1;
                    }
                }
            }
        }
    }
    assert_eq!(cases, 12_168);
}

#[test]
fn primer_payload_masks_reject_invalid_anchor_and_hit_geometry() {
    let anchor = interval_hit(4, 7, PrimerReadStrand::Forward);
    assert!(prepare::payload_geometry(12, &anchor, &[]).is_err());
    for invalid in [
        interval_hit(5, 5, anchor.strand),
        interval_hit(7, 3, anchor.strand),
        interval_hit(8, 13, anchor.strand),
    ] {
        assert!(prepare::payload_geometry(12, &invalid, std::slice::from_ref(&invalid)).is_err());
        assert!(prepare::payload_geometry(12, &anchor, &[anchor.clone(), invalid]).is_err());
    }
}

#[test]
fn primer_payload_masks_handle_empty_fully_masked_and_extreme_coordinates() {
    for anchor in [
        interval_hit(4, 12, PrimerReadStrand::Forward),
        interval_hit(0, 8, PrimerReadStrand::ReverseComplement),
    ] {
        let got = prepare::payload_geometry(12, &anchor, std::slice::from_ref(&anchor)).unwrap();
        assert_eq!(got.usable, 0);
        assert!(got.excluded.is_empty());
    }
    for (anchor, other) in [
        (
            interval_hit(0, 4, PrimerReadStrand::Forward),
            interval_hit(4, 12, PrimerReadStrand::Forward),
        ),
        (
            interval_hit(8, 12, PrimerReadStrand::ReverseComplement),
            interval_hit(0, 8, PrimerReadStrand::Forward),
        ),
    ] {
        let got = prepare::payload_geometry(12, &anchor, &[anchor.clone(), other]).unwrap();
        assert_eq!(got.usable, 0);
        assert_eq!(got.excluded, vec![(0, 8)]);
    }
    let anchor = interval_hit(usize::MAX - 20, usize::MAX - 10, PrimerReadStrand::Forward);
    let got = prepare::payload_geometry(
        usize::MAX,
        &anchor,
        &[anchor.clone(), interval_hit(0, 5, anchor.strand)],
    )
    .unwrap();
    assert_eq!(got.usable, 10);
    assert!(got.excluded.is_empty());
}

#[test]
fn primer_payload_fraction_zero_or_inconsistent_denominator_is_unavailable() {
    assert_eq!(evidence::target_fraction(Some(0), 0), None);
    assert_eq!(evidence::target_fraction(Some(1), 0), None);
    assert_eq!(evidence::target_fraction(Some(11), 10), None);
    assert_eq!(evidence::target_fraction(None, 10), None);
    assert_eq!(evidence::target_fraction(Some(0), 10), Some(0));
    assert_eq!(evidence::target_fraction(Some(5), 10), Some(500));
    assert_eq!(
        evidence::target_fraction(Some(usize::MAX), usize::MAX),
        Some(1000)
    );
}

#[test]
fn primer_payload_masked_queries_round_trip_both_orientations_without_reopening_original() {
    let td = tempfile::tempdir().unwrap();
    let forward = format!("AGTCCGATCA{PAYLOAD}AGTCCGATCA{PAYLOAD}");
    let reverse =
        String::from_utf8(GentleEngine::reverse_complement_bytes(forward.as_bytes())).unwrap();
    let (engine, trace_binding) = source_reads(td.path(), false, &[forward, reverse]);
    let trace = *GentleEngine::inspect_primer_trace(
        &trace_binding.report_path,
        Some(&trace_binding),
        None,
        0,
    )
    .unwrap()
    .report;
    let before = fs::read(&trace.request.files[0].path).unwrap();
    let mut request = prepare(trace_binding.clone(), td.path());
    request.source_mode = PrimerPayloadSourceMode::ReplayVerifiedFastq;
    request.anchor_selection = PrimerPayloadAnchorSelection::All;
    let report = engine
        .prepare_primer_payload_alignment_queries(request.clone(), &mut |_| true)
        .unwrap();
    assert_eq!(report.queries.len(), 4);
    let summary = |record| {
        let mut rows = report
            .queries
            .iter()
            .filter(|q| q.record_index == record)
            .map(|q| {
                (
                    q.query_length,
                    q.usable_bases,
                    q.excluded_query_intervals.clone(),
                    q.sequence_sha256.clone(),
                )
            })
            .collect::<Vec<_>>();
        rows.sort();
        rows
    };
    assert_eq!(summary(0), summary(1));
    assert_eq!(
        report
            .queries
            .iter()
            .filter(|q| !q.excluded_query_intervals.is_empty())
            .count(),
        2
    );
    assert_eq!(fs::read(&trace.request.files[0].path).unwrap(), before);
    fs::remove_file(&trace.request.files[0].path).unwrap();
    let file = bound_file(&Path::new(&request.output_dir).join("report.json"));
    let accepted = binding::queries(&file, &trace, &trace_binding, &request.budgets, 0).unwrap();
    assert_eq!(accepted.queries, report.queries);
}

#[test]
fn primer_payload_prepare_retained_fasta_is_bound_without_reopening_fastq() {
    let td = tempfile::tempdir().unwrap();
    let (engine, binding) = source(td.path(), true, PAYLOAD);
    let request = prepare(binding.clone(), td.path());
    let original: RnaPrimerTraceReport =
        *GentleEngine::inspect_primer_trace(&binding.report_path, Some(&binding), None, 0)
            .unwrap()
            .report;
    fs::remove_file(&original.request.files[0].path).unwrap();
    let report = engine
        .prepare_primer_payload_alignment_queries(request.clone(), &mut |_| true)
        .unwrap();
    assert_eq!(report.status, "completed");
    assert!(!report.original_fastq_reverified);
    assert_eq!(report.queries.len(), 1);
    assert_eq!(
        report.queries[0].raw_payload_interval,
        (10, 10 + PAYLOAD.len())
    );
    assert_eq!(
        report.queries[0].sequence_sha256,
        sha256_prefixed_str(PAYLOAD)
    );
    assert_eq!(report.queries[0].primer_ids, vec!["invented_primer"]);
    assert!(report.queries[0].query_id.starts_with("q_"));
    assert!(!report.queries[0].query_id.contains("reused_name"));
    assert!(
        engine
            .prepare_primer_payload_alignment_queries(request, &mut |_| true)
            .is_err()
    );
}

#[test]
fn primer_payload_prepare_replay_requires_explicit_mode_and_preserves_input() {
    let td = tempfile::tempdir().unwrap();
    let (engine, binding) = source(td.path(), false, PAYLOAD);
    let mut request = prepare(binding.clone(), td.path());
    assert!(
        engine
            .prepare_primer_payload_alignment_queries(request.clone(), &mut |_| true)
            .unwrap_err()
            .message
            .contains("retained_fasta_missing")
    );
    assert!(!Path::new(&request.output_dir).exists());
    let original =
        GentleEngine::inspect_primer_trace(&binding.report_path, Some(&binding), None, 0).unwrap();
    let path = &original.report.request.files[0].path;
    let before = fs::read(path).unwrap();
    request.source_mode = PrimerPayloadSourceMode::ReplayVerifiedFastq;
    request.anchor_selection = PrimerPayloadAnchorSelection::All;
    let report = engine
        .prepare_primer_payload_alignment_queries(request, &mut |_| true)
        .unwrap();
    assert!(report.original_fastq_reverified);
    assert_eq!(report.queries.len(), 1);
    assert_eq!(fs::read(path).unwrap(), before);
}

#[test]
fn primer_payload_prepare_rejects_corrupted_trace_rows_and_wrong_selection() {
    let td = tempfile::tempdir().unwrap();
    let (engine, binding) = source(td.path(), true, PAYLOAD);
    let mut request = prepare(binding, td.path());
    request.file_ids.push("foreign_file".into());
    assert!(
        engine
            .prepare_primer_payload_alignment_queries(request.clone(), &mut |_| true)
            .unwrap_err()
            .message
            .contains("source_selection")
    );
    request.file_ids.pop();
    let rows = Path::new(&request.trace.report_path)
        .parent()
        .unwrap()
        .join("artifacts/file_0000/reads.tsv");
    let mut bytes = fs::read(&rows).unwrap();
    bytes.push(b'\n');
    fs::write(rows, bytes).unwrap();
    assert!(
        engine
            .prepare_primer_payload_alignment_queries(request.clone(), &mut |_| true)
            .is_err()
    );
    assert!(!Path::new(&request.output_dir).exists());
}

#[test]
fn primer_payload_prepare_all_anchors_masks_other_primers_without_counting_them_as_payload() {
    let td = tempfile::tempdir().unwrap();
    let sequence = format!("{PAYLOAD}AGTCCGATCA{PAYLOAD}");
    let (engine, binding) = source(td.path(), false, &sequence);
    let mut request = prepare(binding, td.path());
    request.source_mode = PrimerPayloadSourceMode::ReplayVerifiedFastq;
    request.anchor_selection = PrimerPayloadAnchorSelection::All;
    let report = engine
        .prepare_primer_payload_alignment_queries(request, &mut |_| true)
        .unwrap();
    assert_eq!(report.queries.len(), 2);
    assert_ne!(report.queries[0].anchor_id, report.queries[1].anchor_id);
    let masked = report
        .queries
        .iter()
        .find(|q| !q.excluded_query_intervals.is_empty())
        .unwrap();
    assert_eq!(masked.query_length - masked.usable_bases, 10);
    let fasta =
        fs::read_to_string(Path::new(&report.request.output_dir).join("artifacts/queries.fasta"))
            .unwrap();
    assert!(fasta.contains("NNNNNNNNNN"));
}

#[test]
fn primer_payload_prepare_control_errors_do_not_publish_unverified_units() {
    let td = tempfile::tempdir().unwrap();
    let (engine, binding) = source(td.path(), true, PAYLOAD);
    let mut request = prepare(binding, td.path());
    let marker = td.path().join("control marker");
    request.cancel_path = Some(marker.display().to_string());
    let report = engine
        .prepare_primer_payload_alignment_queries(request, &mut |_| {
            fs::write(&marker, "not the request digest").unwrap();
            true
        })
        .unwrap();
    assert_eq!(report.status, "control_error");
    assert!(report.files_completed.is_empty());
    assert!(report.queries.is_empty());
    assert_eq!(report.files_not_completed, vec!["file_0000"]);
}

#[test]
fn primer_payload_prepare_schema_rejects_unknown_fields_and_implicit_replay() {
    let td = tempfile::tempdir().unwrap();
    let (_, binding) = source(td.path(), true, PAYLOAD);
    let request = prepare(binding, td.path());
    let mut value = serde_json::to_value(&request).unwrap();
    value["run_aligner"] = serde_json::json!(true);
    assert!(serde_json::from_value::<PreparePrimerPayloadAlignmentQueriesRequest>(value).is_err());
    let mut value = serde_json::to_value(request).unwrap();
    value.as_object_mut().unwrap().remove("source_mode");
    assert!(serde_json::from_value::<PreparePrimerPayloadAlignmentQueriesRequest>(value).is_err());
}
