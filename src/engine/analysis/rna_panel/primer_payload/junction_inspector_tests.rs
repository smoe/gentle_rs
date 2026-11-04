//! Hand-crafted synthetic reads, oligos, references and PAF rows recreated in private temp dirs.
//! Tests use the real trace/query/payload producers and verified saved-file boundary; no real RNA.

use super::super::{alignment_tests, tests as fixture};
use super::*;
use crate::engine_shell::{execute_shell_command, parse_shell_line, quote_shell_arg};

const PRIMER: &str = "AGTCCGATCA";

struct Fixture {
    engine: GentleEngine,
    prepared: PrimerPayloadQueries,
    report: PrimerPayloadAlignmentReport,
    inspect: InspectPrimerPayloadJunctionsRequest,
}

fn bound(path: &Path) -> PrimerPayloadFile {
    fixture::bound_file(path)
}

fn inspection(report: &PrimerPayloadAlignmentReport) -> InspectPrimerPayloadJunctionsRequest {
    let source =
        GentleEngine::primer_payload_alignment_result_binding(&report.request, report).unwrap();
    InspectPrimerPayloadJunctionsRequest {
        schema: PRIMER_JUNCTION_INSPECTION_REQUEST_SCHEMA.into(),
        source_report: bound(Path::new(&source.report_path)),
        source,
        file_id: None,
        primer_id: None,
        offset: 0,
        rows: 100,
        budgets: primer_junction_inspection_default_budgets(),
    }
}

fn saved(
    root: &Path,
    reads: &[String],
    all_anchors: bool,
    ambiguous: bool,
    outcome: &str,
) -> Fixture {
    let (engine, trace) = fixture::source_reads(root, true, reads);
    from_trace(root, engine, trace, all_anchors, ambiguous, outcome)
}

fn from_trace(
    root: &Path,
    engine: GentleEngine,
    trace: RnaPrimerTraceResultBinding,
    all_anchors: bool,
    ambiguous: bool,
    outcome: &str,
) -> Fixture {
    let mut prepare = fixture::prepare(trace.clone(), root);
    if all_anchors {
        prepare.source_mode = PrimerPayloadSourceMode::ReplayVerifiedFastq;
        prepare.anchor_selection = PrimerPayloadAnchorSelection::All;
    }
    let prepared = engine
        .prepare_primer_payload_alignment_queries(prepare, &mut |_| true)
        .unwrap();
    let mut request = alignment_tests::request(root, trace, &prepared);
    for primer in prepared
        .queries
        .iter()
        .flat_map(|q| q.primer_ids.iter())
        .collect::<BTreeSet<_>>()
    {
        if !request
            .expected_targets
            .iter()
            .any(|e| &e.primer_id == primer)
        {
            let mut expected = request.expected_targets[0].clone();
            expected.primer_id = primer.to_string();
            request.expected_targets.push(expected);
        }
    }

    // Extend the invented reference for long/masked multi-anchor queries while preserving
    // the fixed ten-base primer-adjacent seed socket. PAF rows are synthetic observations.
    let length = prepared
        .queries
        .iter()
        .map(|q| q.query_length)
        .max()
        .unwrap_or(fixture::PAYLOAD.len())
        + 10;
    let reference = &mut request.references[0];
    fs::write(
        &reference.fasta.path,
        format!(
            ">toy_expected\n{}\n>toy_foreign\n{PRIMER}{}\n",
            "C".repeat(length),
            "G".repeat(length - 10)
        ),
    )
    .unwrap();
    reference.fasta = bound(Path::new(&reference.fasta.path));
    let mut annotation: PrimerPayloadAnnotation =
        serde_json::from_slice(&fs::read(&reference.annotation.path).unwrap()).unwrap();
    annotation.reference_fasta_sha256 = reference.fasta.sha256.clone();
    for t in &mut annotation.transcripts {
        t.target_length = length;
    }
    fs::write(
        &reference.annotation.path,
        serde_json::to_vec(&annotation).unwrap(),
    )
    .unwrap();
    reference.annotation = bound(Path::new(&reference.annotation.path));
    let paf = fs::read_to_string(&request.batches[0].paf.path).unwrap();
    let mut rows = String::new();
    for line in paf.lines() {
        let mut fields = line.split('\t').map(str::to_owned).collect::<Vec<_>>();
        fields[6] = length.to_string();
        rows.push_str(&fields.join("\t"));
        rows.push('\n');
        if ambiguous {
            fields[5] = "toy_expected".into();
            rows.push_str(&fields.join("\t"));
            rows.push('\n');
        }
    }
    fs::write(&request.batches[0].paf.path, rows).unwrap();
    request.batches[0].paf = bound(Path::new(&request.batches[0].paf.path));
    let report = match outcome {
        "cancelled" => engine
            .summarize_primer_payload_alignments(request, &mut |_| false)
            .unwrap(),
        "control_error" => {
            let marker = root.join("synthetic-control-marker");
            request.cancel_path = Some(marker.display().to_string());
            engine
                .summarize_primer_payload_alignments(request, &mut |_| {
                    fs::write(&marker, "wrong-request-identity").unwrap();
                    true
                })
                .unwrap()
        }
        _ => engine
            .summarize_primer_payload_alignments(request, &mut |_| true)
            .unwrap(),
    };
    let inspect = inspection(&report);
    Fixture {
        engine,
        prepared,
        report,
        inspect,
    }
}

fn one(root: &Path) -> Fixture {
    saved(
        root,
        &[format!("{PRIMER}{}", fixture::PAYLOAD)],
        false,
        false,
        "completed",
    )
}

fn run(f: &Fixture) -> PrimerPayloadJunctionInspection {
    f.engine
        .inspect_primer_payload_junctions(f.inspect.clone(), &mut |_| true)
        .unwrap()
}

fn reseal(f: &mut Fixture) {
    f.report.content_sha256.clear();
    f.report.content_sha256 = hash(&f.report).unwrap();
    let path = Path::new(&f.inspect.source.report_path);
    fs::write(path, serde_json::to_vec_pretty(&f.report).unwrap()).unwrap();
    f.inspect = inspection(&f.report);
}

fn table<T: Serialize>(f: &mut Fixture, name: &str, header: &str, rows: &[T]) {
    let mut csv = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(Vec::new());
    csv.write_record([header]).unwrap();
    for row in rows {
        csv.write_record([serde_json::to_string(row).unwrap()])
            .unwrap();
    }
    let bytes = csv.into_inner().unwrap();
    let path = Path::new(&f.inspect.source.report_path)
        .parent()
        .unwrap()
        .join("artifacts")
        .join(name);
    fs::write(&path, &bytes).unwrap();
    let artifact = f
        .report
        .artifacts
        .iter_mut()
        .find(|a| a.path == name)
        .unwrap();
    artifact.bytes = bytes.len() as u64;
    artifact.sha256 = sha256_prefixed_bytes(&bytes);
}

fn custom_trace(
    root: &Path,
    modify: impl FnOnce(&mut RnaPrimerPanel),
) -> (GentleEngine, RnaPrimerTraceResultBinding) {
    let (engine, binding) =
        fixture::source_reads(root, true, &[format!("{PRIMER}{}", fixture::PAYLOAD)]);
    let old: RnaPrimerTraceReport =
        serde_json::from_slice(&fs::read(&binding.report_path).unwrap()).unwrap();
    let mut request = old.request;
    let RnaPrimerPanelBinding::Inline { panel } = &mut request.primer_panel else {
        panic!("synthetic inline panel")
    };
    modify(panel);
    panel.content_sha256.clear();
    GentleEngine::seal_rna_primer_panel(panel).unwrap();
    request.output_dir = root
        .join("alternate synthetic phase one")
        .display()
        .to_string();
    let report = engine
        .trace_primer_reads(request.clone(), &mut |_| true)
        .unwrap();
    let binding = GentleEngine::primer_trace_result_binding(&request, &report).unwrap();
    (engine, binding)
}

#[test]
fn primer_junction_inspector_same_geometry_primer_alternatives_keep_their_own_edit_evidence() {
    let td = tempfile::tempdir().unwrap();
    let (engine, trace) = custom_trace(td.path(), |panel| {
        panel.policy.max_full_edits = 1;
        panel.policy.max_core_edits_for_attribution = 1;
        let mut second = panel.primers[0].clone();
        second.primer_id = "alternative_primer".into();
        second.sequence = "AGTCTGATCA".into();
        panel.primers.push(second);
    });
    let f = from_trace(td.path(), engine, trace, false, false, "completed");
    let inspected = run(&f);
    let row = &inspected.rows[0];
    assert_eq!(
        row.queries[0].primer_ids,
        vec!["alternative_primer", "invented_primer"]
    );
    let hit = |id: &str| {
        row.trace
            .hits
            .iter()
            .find(|h| {
                h.primer_id.as_deref() == Some(id)
                    && h.start == row.queries[0].anchor.start
                    && h.end == row.queries[0].anchor.end
            })
            .unwrap()
    };
    assert_eq!(hit("invented_primer").edits, 0);
    assert_eq!(hit("alternative_primer").edits, 1);
    assert_eq!(hit("invented_primer").anchor_class, "full_exact");
    assert_eq!(hit("alternative_primer").anchor_class, "full_approximate");
    assert_eq!(
        inspected.profiles[0].primer_ids,
        vec!["alternative_primer", "invented_primer"]
    );
    assert_eq!(inspected.profiles[0].supporting_reads, 1);
}

#[test]
fn primer_junction_inspector_core_only_observation_never_invents_complete_handle_trim() {
    let td = tempfile::tempdir().unwrap();
    let handle = "TTGGTTAA";
    let (engine, trace) = custom_trace(td.path(), |panel| {
        let primer = &mut panel.primers[0];
        primer.sequence = format!("{handle}{PRIMER}");
        primer.components = vec![
            RnaPrimerComponent {
                component_id: "handle".into(),
                role: "common_handle".into(),
                start: 0,
                end: handle.len(),
            },
            RnaPrimerComponent {
                component_id: "core".into(),
                role: "target_core".into(),
                start: handle.len(),
                end: handle.len() + PRIMER.len(),
            },
        ];
    });
    let f = from_trace(td.path(), engine, trace, false, false, "completed");
    let inspected = run(&f);
    let row = &inspected.rows[0];
    assert_eq!(row.queries[0].anchor.anchor_class, "core_exact");
    assert_eq!(
        row.queries[0].raw_payload_interval,
        (PRIMER.len(), PRIMER.len() + fixture::PAYLOAD.len())
    );
    assert_eq!(
        row.queries[0].anchor.end - row.queries[0].anchor.start,
        PRIMER.len()
    );
    assert_eq!(
        inspected.primer_panel.primers[0].sequence.len(),
        handle.len() + PRIMER.len()
    );
}

#[test]
fn primer_junction_inspector_orients_saved_queries_without_reorienting_original_trace() {
    let td = tempfile::tempdir().unwrap();
    let forward = format!("GGG{PRIMER}{}", fixture::PAYLOAD);
    let reverse = GentleEngine::reverse_complement(&forward);
    let f = saved(td.path(), &[forward, reverse], false, false, "completed");
    let inspected = run(&f);
    assert_eq!(inspected.rows.len(), 2);
    assert_eq!(inspected.cohort_distinct_reads, 2);
    let (a, b) = (&inspected.rows[0], &inspected.rows[1]);
    assert_eq!(a.queries[0].anchor.strand, PrimerReadStrand::Forward);
    assert_eq!(
        b.queries[0].anchor.strand,
        PrimerReadStrand::ReverseComplement
    );
    assert_ne!(
        a.queries[0].raw_payload_interval,
        b.queries[0].raw_payload_interval
    );
    assert_eq!(a.queries[0].sequence_sha256, b.queries[0].sequence_sha256);
    assert_eq!(a.queries[0].query_length, b.queries[0].query_length);
    assert_eq!(a.payload_previews[&a.queries[0].query_id], fixture::PAYLOAD);
    assert_eq!(b.payload_previews[&b.queries[0].query_id], fixture::PAYLOAD);

    assert_eq!(
        a.trace.hits,
        a.queries
            .iter()
            .map(|q| q.anchor.clone())
            .collect::<Vec<_>>()
    );
    assert!(!a.evidence[0].terminal_seeds.is_empty());
    assert!(inspected.saved_artifacts_verified);
    assert!(!inspected.original_reads_reverified);
    assert!(!inspected.reference_bases_reverified);
    assert_eq!(
        serde_json::to_value(inspected.profiles).unwrap(),
        serde_json::to_value(&f.report.profiles).unwrap()
    );
}

#[test]
fn primer_junction_inspector_distinct_records_are_not_duplicate_read_names_or_anchor_counts() {
    let td = tempfile::tempdir().unwrap();
    let read = format!("{PRIMER}{}", fixture::PAYLOAD);
    let f = saved(td.path(), &[read.clone(), read], false, false, "completed");
    let inspected = run(&f);
    assert_eq!(
        inspected.rows[0].trace.read_id,
        inspected.rows[1].trace.read_id
    );
    assert_ne!(
        inspected.rows[0].trace.record_index,
        inspected.rows[1].trace.record_index
    );
    assert_eq!(inspected.cohort_distinct_reads, 2);
    let profile = inspected
        .profiles
        .iter()
        .find(|p| !p.primer_ids.is_empty())
        .unwrap();
    assert_eq!(
        (
            profile.supporting_reads,
            profile.eligible_reads,
            profile.anchor_observations
        ),
        (2, 2, 2)
    );

    let td = tempfile::tempdir().unwrap();
    let f = saved(
        td.path(),
        &[format!(
            "{PRIMER}{}{PRIMER}{}",
            fixture::PAYLOAD,
            fixture::PAYLOAD
        )],
        true,
        false,
        "completed",
    );
    let inspected = run(&f);
    assert_eq!(inspected.cohort_distinct_reads, 1);
    assert_eq!(inspected.rows[0].queries.len(), 2);
    assert!(inspected.rows[0].trace.hits.len() >= 2);
    let profile = inspected
        .profiles
        .iter()
        .find(|p| !p.primer_ids.is_empty())
        .unwrap();
    assert_eq!(
        (
            profile.supporting_reads,
            profile.eligible_reads,
            profile.anchor_observations
        ),
        (1, 1, 2)
    );
}

#[test]
fn primer_junction_inspector_ambiguous_partner_equivalence_is_not_split_into_gene_counts() {
    let td = tempfile::tempdir().unwrap();
    let f = saved(
        td.path(),
        &[format!("{PRIMER}{}", fixture::PAYLOAD)],
        false,
        true,
        "completed",
    );
    let inspected = run(&f);
    assert_eq!(inspected.rows[0].evidence[0].classification, "ambiguous");
    assert_eq!(
        inspected.rows[0].evidence[0].primer_edge_alternatives.len(),
        2
    );
    let profile = &inspected.profiles[0];
    assert_eq!(
        profile.partner_gene_ids,
        vec!["expected_gene", "foreign_gene"]
    );
    assert_eq!(
        (
            profile.supporting_reads,
            profile.ambiguous_reads,
            profile.anchor_observations
        ),
        (1, 1, 1)
    );
    assert_eq!(inspected.profiles.len(), 1);
    assert!(inspected.breakpoint_clusters.is_empty());
}

#[test]
fn primer_junction_inspector_filters_and_pagination_keep_global_cohort_and_profile_denominators() {
    let td = tempfile::tempdir().unwrap();
    let mut f = saved(
        td.path(),
        &[format!("{PRIMER}{}", fixture::PAYLOAD), "A".repeat(80)],
        false,
        false,
        "completed",
    );
    let whole = run(&f);
    assert_eq!(whole.cohort_distinct_reads, 2);
    assert_eq!(whole.rows.len(), 2);
    let unanchored = whole.rows.iter().find(|r| r.queries.is_empty()).unwrap();
    assert!(
        unanchored
            .trace
            .anchor
            .as_ref()
            .is_none_or(|a| a.primer_id.is_none())
    );
    assert_eq!(unanchored.unassessed[0].reason, "no_attributed_anchor");
    f.inspect.primer_id = Some("invented_primer".into());
    f.inspect.rows = 0;
    let selected = run(&f);
    assert_eq!(selected.cohort_distinct_reads, 2);
    assert_eq!(selected.matched_trace_rows, 1);
    assert!(selected.rows.is_empty());
    assert!(selected.rows_truncated);
    assert_eq!(
        serde_json::to_value(selected.profiles).unwrap(),
        serde_json::to_value(&whole.profiles).unwrap()
    );
    f.inspect.primer_id = None;
    f.inspect.offset = 1;
    f.inspect.rows = 1;
    let page = run(&f);
    assert_eq!(
        (
            page.cohort_distinct_reads,
            page.matched_trace_rows,
            page.rows.len()
        ),
        (2, 2, 1)
    );
    assert_eq!(page.rows[0].trace.record_index, 1);
}

#[test]
fn primer_junction_inspector_never_reopens_original_fastq_references_annotations_or_paf() {
    let td = tempfile::tempdir().unwrap();
    let f = one(td.path());
    let trace: RnaPrimerTraceReport =
        serde_json::from_slice(&fs::read(&f.report.request.trace.report_path).unwrap()).unwrap();
    for unit in &trace.files_completed {
        fs::remove_file(&unit.input.path).unwrap();
    }
    for r in &f.report.request.references {
        fs::remove_file(&r.fasta.path).unwrap();
        fs::remove_file(&r.annotation.path).unwrap();
    }
    for b in &f.report.request.batches {
        fs::remove_file(&b.paf.path).unwrap();
    }
    let inspected = run(&f);
    assert_eq!(inspected.rows.len(), 1);
    assert_eq!(
        inspected.rows[0].queries[0].query_id,
        f.prepared.queries[0].query_id
    );
    assert!(!inspected.original_reads_reverified);
    assert!(!inspected.reference_bases_reverified);
}

#[test]
fn primer_junction_inspector_missing_retained_artifacts_refuse_without_implicit_replay() {
    let td = tempfile::tempdir().unwrap();
    let f = one(td.path());
    let trace: RnaPrimerTraceReport =
        serde_json::from_slice(&fs::read(&f.report.request.trace.report_path).unwrap()).unwrap();
    let unit = &trace.files_completed[0];
    let payload = unit
        .artifacts
        .iter()
        .find(|a| a.path.ends_with("payload.fasta"))
        .unwrap();
    let path = Path::new(&f.report.request.trace.report_path)
        .parent()
        .unwrap()
        .join("artifacts")
        .join(&payload.path);
    fs::remove_file(path).unwrap();
    assert!(
        f.engine
            .inspect_primer_payload_junctions(f.inspect, &mut |_| true)
            .is_err()
    );
}

#[test]
fn primer_junction_inspector_raw_hash_binding_and_nonselected_artifacts_are_verified() {
    for which in 0..3 {
        let td = tempfile::tempdir().unwrap();
        let mut f = one(td.path());
        f.inspect.rows = 0;
        match which {
            0 => f.inspect.source_report.sha256 = format!("sha256:{}", "0".repeat(64)),
            1 => {
                let path = Path::new(&f.inspect.source.report_path)
                    .parent()
                    .unwrap()
                    .join("artifacts/partner_profiles.tsv");
                let mut bytes = fs::read(&path).unwrap();
                bytes.push(b' ');
                fs::write(path, bytes).unwrap();
            }
            _ => {
                let path = Path::new(&f.inspect.source.report_path)
                    .parent()
                    .unwrap()
                    .join("artifacts/read_anchor_evidence.tsv");
                let mut bytes = fs::read(&path).unwrap();
                *bytes.last_mut().unwrap() = b' ';
                fs::write(path, bytes).unwrap();
            }
        }
        assert!(
            f.engine
                .inspect_primer_payload_junctions(f.inspect, &mut |_| true)
                .is_err(),
            "case {which}"
        );
    }
}

#[test]
fn primer_junction_inspector_rehashed_missing_rows_or_profile_arithmetic_cannot_pass() {
    for which in 0..2 {
        let td = tempfile::tempdir().unwrap();
        let mut f = one(td.path());
        if which == 0 {
            table::<PrimerPayloadReadEvidence>(
                &mut f,
                "read_anchor_evidence.tsv",
                "evidence_json",
                &[],
            );
            f.report.evidence_rows = 0;
        } else {
            f.report.profiles[0].supporting_reads += 1;
            let profiles = f.report.profiles.clone();
            table(
                &mut f,
                "partner_profiles.tsv",
                "partner_profile_json",
                &profiles,
            );
        }
        reseal(&mut f);
        assert!(
            f.engine
                .inspect_primer_payload_junctions(f.inspect, &mut |_| true)
                .is_err(),
            "case {which}"
        );
    }
}

fn evidence_rows(f: &Fixture) -> Vec<PrimerPayloadReadEvidence> {
    let path = Path::new(&f.inspect.source.report_path)
        .parent()
        .unwrap()
        .join("artifacts/read_anchor_evidence.tsv");
    let mut csv = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .unwrap();
    csv.records()
        .map(|r| serde_json::from_str(&r.unwrap()[0]).unwrap())
        .collect()
}

fn synchronize_evidence(f: &mut Fixture, rows: &[PrimerPayloadReadEvidence]) {
    let trace: RnaPrimerTraceReport =
        serde_json::from_slice(&fs::read(&f.report.request.trace.report_path).unwrap()).unwrap();
    let references = f
        .report
        .request
        .batches
        .iter()
        .filter(|b| f.report.batches_completed.contains(&b.batch_id))
        .map(|b| b.reference_id.clone())
        .collect::<Vec<_>>();
    let (profiles, clusters) = evidence::aggregate(
        rows,
        &f.prepared,
        &trace,
        &references,
        f.report.request.policy.breakpoint_tolerance_bases,
    );
    table(f, "read_anchor_evidence.tsv", "evidence_json", rows);
    table(f, "partner_profiles.tsv", "partner_profile_json", &profiles);
    table(
        f,
        "breakpoint_clusters.tsv",
        "breakpoint_cluster_json",
        &clusters,
    );
    f.report.evidence_rows = rows.len();
    f.report.profiles = profiles;
    f.report.breakpoint_clusters = clusters;
    reseal(f);
}

#[test]
fn primer_junction_inspector_rehashed_mapped_rows_cannot_claim_producer_unmapped() {
    let td = tempfile::tempdir().unwrap();
    let mut f = one(td.path());
    let mut rows = evidence_rows(&f);
    assert!(!rows[0].segments.is_empty());
    assert_eq!(
        f.report.request.batches[0].accounting[&rows[0].query_id],
        PrimerPayloadQueryAccounting::Mapped
    );
    rows[0].classification = "producer_reported_unmapped".into();
    synchronize_evidence(&mut f, &rows);
    assert!(
        f.engine
            .inspect_primer_payload_junctions(f.inspect, &mut |_| true)
            .is_err()
    );
}

#[test]
fn primer_junction_inspector_rehashed_invalid_internal_boundary_strand_is_refused() {
    let td = tempfile::tempdir().unwrap();
    let mut f = one(td.path());
    let mut rows = evidence_rows(&f);
    let mut left = rows[0].primer_edge_alternatives[0].clone();
    let right = left.clone();
    left.strand = '?';
    rows[0]
        .internal_query_boundaries
        .push(PrimerPayloadBoundary {
            left_query_end: 1,
            right_query_start: 2,
            left,
            right,
        });
    synchronize_evidence(&mut f, &rows);
    assert!(
        f.engine
            .inspect_primer_payload_junctions(f.inspect, &mut |_| true)
            .is_err()
    );
}

#[test]
fn primer_junction_inspector_saved_partial_outcomes_and_missing_batches_are_not_completed_or_zero_evidence()
 {
    for outcome in ["cancelled", "control_error"] {
        let td = tempfile::tempdir().unwrap();
        let f = saved(
            td.path(),
            &[format!("{PRIMER}{}", fixture::PAYLOAD)],
            false,
            false,
            outcome,
        );
        assert_eq!(f.report.status, outcome);
        let inspected = run(&f);
        assert_eq!(inspected.source_status, outcome);
        assert_eq!(
            inspected.source_control_error.is_some(),
            outcome == "control_error"
        );
        assert_eq!(inspected.trace_status, "completed");
        assert_eq!(inspected.query_status, "completed");
        assert_eq!(inspected.batches_not_completed.len(), 1);
        assert_eq!(inspected.cohort_distinct_reads, 1);
        assert_eq!(inspected.rows[0].queries.len(), 1);
        assert!(inspected.rows[0].evidence.is_empty());
    }
}

#[test]
fn primer_junction_inspector_shared_shell_typed_route_identity_and_read_only_effects() {
    let td = tempfile::tempdir().unwrap();
    let mut f = one(td.path());
    let before = encode(f.engine.state()).unwrap();
    let revision = f.engine.mutation_revision();
    let undo = f.engine.undo_available();
    let result = f
        .engine
        .apply(Operation::InspectPrimerPayloadJunctions {
            request: Box::new(f.inspect.clone()),
        })
        .unwrap();
    let expected = serde_json::to_value(result.primer_payload_junctions.unwrap()).unwrap();
    assert_eq!(encode(f.engine.state()).unwrap(), before);
    assert_eq!(f.engine.mutation_revision(), revision);
    assert_eq!(f.engine.undo_available(), undo);
    let journal = f.engine.journal_len();
    let command = parse_shell_line(&format!(
        "rna-reads primer-junctions {}",
        quote_shell_arg(&serde_json::to_string(&f.inspect).unwrap())
    ))
    .unwrap();
    assert!(!command.is_state_mutating());
    let actual = execute_shell_command(&mut f.engine, &command).unwrap();
    assert!(!actual.state_changed);
    assert_eq!(
        actual.output["result"]["primer_payload_junctions"],
        expected
    );
    assert_eq!(f.engine.journal_len(), journal);
    assert_eq!(encode(f.engine.state()).unwrap(), before);
    let path = td.path().join("junction request's metadata.json");
    fs::write(&path, serde_json::to_vec(&f.inspect).unwrap()).unwrap();
    let command = parse_shell_line(&format!(
        "rna-reads primer-junctions {}",
        quote_shell_arg(&format!("@{}", path.display()))
    ))
    .unwrap();
    assert_eq!(
        execute_shell_command(&mut f.engine, &command)
            .unwrap()
            .output["result"]["primer_payload_junctions"],
        expected
    );
    let inspected = run(&f);
    let mut basis = inspected.clone();
    basis.content_sha256.clear();
    assert_eq!(inspected.content_sha256, hash(&basis).unwrap());
}

#[test]
fn primer_junction_inspector_budget_cancellation_and_unknown_fields_refuse_without_output() {
    let td = tempfile::tempdir().unwrap();
    let f = one(td.path());
    assert!(
        f.engine
            .inspect_primer_payload_junctions(f.inspect.clone(), &mut |_| false)
            .is_err()
    );
    for which in 0..3 {
        let mut req = f.inspect.clone();
        match which {
            0 => req.budgets.max_input_bytes = 1,
            1 => req.rows = 1001,
            _ => req.budgets.max_records = 1,
        }
        if which == 2 {
            // One-record fixture is valid at this cap; a response smaller than one row is not.
            req.budgets.max_output_bytes = 1;
        }
        assert!(
            f.engine
                .inspect_primer_payload_junctions(req, &mut |_| true)
                .is_err(),
            "case {which}"
        );
    }
    let mut json = serde_json::to_value(&f.inspect).unwrap();
    json["implicit_fastq_replay"] = serde_json::json!(true);
    assert!(serde_json::from_value::<InspectPrimerPayloadJunctionsRequest>(json).is_err());
}
