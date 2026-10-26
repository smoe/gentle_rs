//! Invented source units and deliberately rehashed query tampering; no private data or trusted-producer claim.

use super::tests::{PAYLOAD, bound_file, limits, prepare, source_reads};
use super::*;

fn masked_read() -> String {
    format!("AGTCCGATCA{PAYLOAD}AGTCCGATCA{PAYLOAD}")
}

fn prepared(
    root: &Path,
    reads: &[String],
    policy: PrimerPayloadInternalPrimerPolicy,
    minimum: usize,
) -> (
    GentleEngine,
    RnaPrimerTraceResultBinding,
    RnaPrimerTraceReport,
    PrimerPayloadQueries,
) {
    let (engine, binding) = source_reads(root, true, reads);
    let trace = *GentleEngine::inspect_primer_trace(&binding.report_path, Some(&binding), None, 0)
        .unwrap()
        .report;
    let mut request = prepare(binding.clone(), root);
    request.internal_primer_policy = policy;
    request.min_payload_bases = minimum;
    let report = engine
        .prepare_primer_payload_alignment_queries(request, &mut |_| true)
        .unwrap();
    (engine, binding, trace, report)
}

fn persisted(report: &mut PrimerPayloadQueries) -> PrimerPayloadFile {
    report.request_sha256 = hash(&report.request).unwrap();
    report.content_sha256.clear();
    report.content_sha256 = hash(report).unwrap();
    let path = Path::new(&report.request.output_dir).join("report.json");
    fs::write(&path, serde_json::to_vec_pretty(report).unwrap()).unwrap();
    bound_file(&path)
}

fn accept(
    report: &PrimerPayloadQueries,
    trace: &RnaPrimerTraceReport,
    binding: &RnaPrimerTraceResultBinding,
) -> Result<PrimerPayloadQueries, EngineError> {
    binding::queries(
        &bound_file(&Path::new(&report.request.output_dir).join("report.json")),
        trace,
        binding,
        &limits(),
        0,
    )
}

fn source_bytes(
    trace: &RnaPrimerTraceReport,
    binding: &RnaPrimerTraceResultBinding,
) -> BTreeMap<PathBuf, Vec<u8>> {
    let root = Path::new(&binding.report_path)
        .parent()
        .unwrap()
        .join("artifacts");
    let mut paths = vec![PathBuf::from(&binding.report_path)];
    paths.extend(trace.request.files.iter().map(|f| PathBuf::from(&f.path)));
    for artifact in trace
        .artifacts
        .iter()
        .chain(trace.files_completed.iter().flat_map(|f| &f.artifacts))
    {
        paths.push(PathBuf::from(artifact_file(&root, artifact).unwrap().path));
    }
    paths
        .into_iter()
        .map(|p| {
            let bytes = fs::read(&p).unwrap();
            (p, bytes)
        })
        .collect()
}

#[test]
fn primer_payload_bound_exclusions_accept_valid_dispositions_in_any_order() {
    for (read, policy, minimum, reason) in [
        (
            format!("AGTCCGATCA{PAYLOAD}"),
            PrimerPayloadInternalPrimerPolicy::Mask,
            12,
            None,
        ),
        (
            masked_read(),
            PrimerPayloadInternalPrimerPolicy::Mask,
            12,
            None,
        ),
        (
            masked_read(),
            PrimerPayloadInternalPrimerPolicy::Unassessed,
            12,
            Some("additional_primer_interval_unassessed"),
        ),
        (
            format!("AGTCCGATCA{PAYLOAD}"),
            PrimerPayloadInternalPrimerPolicy::Mask,
            100,
            Some("insufficient_primer_free_payload"),
        ),
        (
            masked_read(),
            PrimerPayloadInternalPrimerPolicy::Unassessed,
            1000,
            Some("insufficient_primer_free_payload"),
        ),
        (
            PAYLOAD.into(),
            PrimerPayloadInternalPrimerPolicy::Mask,
            12,
            Some("no_attributed_anchor"),
        ),
    ] {
        let td = tempfile::tempdir().unwrap();
        let (_engine, binding, trace, mut report) =
            prepared(td.path(), &[read.clone(), read], policy, minimum);
        let before = source_bytes(&trace, &binding);
        let accepted = accept(&report, &trace, &binding).unwrap();
        assert_eq!(encode(&accepted).unwrap(), encode(&report).unwrap());
        if let Some(reason) = reason {
            assert_eq!(report.unassessed.len(), 2);
            assert!(report.unassessed.iter().all(|r| r.reason == reason));
        } else {
            assert_eq!(report.queries.len(), 2);
            assert!(report.unassessed.is_empty());
        }
        report.queries.reverse();
        report.unassessed.reverse();
        persisted(&mut report);
        assert!(accept(&report, &trace, &binding).is_ok());
        assert_eq!(source_bytes(&trace, &binding), before);
    }
}

#[test]
fn primer_payload_rehashed_eligible_anchor_cannot_be_hidden_as_excluded() {
    let td = tempfile::tempdir().unwrap();
    let (engine, binding, trace, mut report) = prepared(
        td.path(),
        &[format!("AGTCCGATCA{PAYLOAD}")],
        PrimerPayloadInternalPrimerPolicy::Mask,
        12,
    );
    let before = source_bytes(&trace, &binding);
    let baseline = alignment_tests::request(td.path(), binding.clone(), &report);
    engine
        .summarize_primer_payload_alignments(baseline.clone(), &mut |_| true)
        .unwrap();
    let query = report.queries.remove(0);
    report.unassessed.push(PrimerPayloadUnassessed {
        file_id: query.file_id,
        record_index: query.record_index,
        read_id: query.read_id,
        anchor_id: Some(query.anchor_id),
        primer_ids: query.primer_ids,
        reason: "insufficient_primer_free_payload".into(),
    });
    let fasta = Path::new(&report.request.output_dir).join("artifacts/queries.fasta");
    fs::write(&fasta, b"").unwrap();
    let raw = bound_file(&fasta);
    report.artifacts[0].sha256 = raw.sha256;
    report.artifacts[0].bytes = raw.bytes;
    let manifest = persisted(&mut report);
    let mut request = baseline;
    request.queries = manifest;
    request.query_manifest_sha256 = report.content_sha256.clone();
    request.batches[0].query_manifest_sha256 = report.content_sha256.clone();
    request.batches[0].accounting.clear();
    let paf = td.path().join("empty forged batch.paf");
    fs::write(&paf, b"").unwrap();
    request.batches[0].paf = bound_file(&paf);
    request.output_dir = td.path().join("refused import").display().to_string();
    let error = engine
        .summarize_primer_payload_alignments(request.clone(), &mut |_| true)
        .unwrap_err();
    assert!(
        error.message.contains("exclusion_policy_binding"),
        "{error:?}"
    );
    assert!(!Path::new(&request.output_dir).exists());
    assert_eq!(source_bytes(&trace, &binding), before);
}

#[test]
fn primer_payload_rehashed_exclusion_cannot_move_between_repeated_read_records() {
    let td = tempfile::tempdir().unwrap();
    let (_engine, binding, trace, mut report) = prepared(
        td.path(),
        &[masked_read(), masked_read()],
        PrimerPayloadInternalPrimerPolicy::Unassessed,
        12,
    );
    assert_eq!(report.unassessed.len(), 2);
    assert_eq!(report.unassessed[0].read_id, report.unassessed[1].read_id);
    assert_ne!(
        report.unassessed[0].record_index,
        report.unassessed[1].record_index
    );
    assert!(accept(&report, &trace, &binding).is_ok());
    let id = report.unassessed[0].anchor_id.clone();
    report.unassessed[0].anchor_id = report.unassessed[1].anchor_id.clone();
    report.unassessed[1].anchor_id = id;
    persisted(&mut report);
    let error = accept(&report, &trace, &binding).unwrap_err();
    assert!(error.message.contains("exclusion_binding"), "{error:?}");
}

#[test]
fn primer_payload_rehashed_exclusion_fields_and_reason_are_source_bound() {
    let td = tempfile::tempdir().unwrap();
    let (_engine, binding, trace, baseline) = prepared(
        td.path(),
        &[masked_read()],
        PrimerPayloadInternalPrimerPolicy::Unassessed,
        12,
    );
    let before = source_bytes(&trace, &binding);
    assert!(accept(&baseline, &trace, &binding).is_ok());
    for case in [
        "empty_primers",
        "duplicate_primers",
        "wrong_reason",
        "empty_reason",
        "wrong_read",
        "wrong_record",
        "wrong_unit",
    ] {
        let mut report = baseline.clone();
        let row = &mut report.unassessed[0];
        match case {
            "empty_primers" => row.primer_ids.clear(),
            "duplicate_primers" => row.primer_ids.push(row.primer_ids[0].clone()),
            "wrong_reason" => row.reason = "insufficient_primer_free_payload".into(),
            "empty_reason" => row.reason.clear(),
            "wrong_read" => row.read_id = "foreign_read".into(),
            "wrong_record" => row.record_index = u64::MAX,
            _ => row.file_id = "foreign_unit".into(),
        }
        persisted(&mut report);
        assert!(accept(&report, &trace, &binding).is_err(), "{case}");
        assert_eq!(source_bytes(&trace, &binding), before);
    }
}

#[test]
fn primer_payload_rehashed_policy_flip_cannot_change_a_retained_disposition() {
    for policy in [
        PrimerPayloadInternalPrimerPolicy::Mask,
        PrimerPayloadInternalPrimerPolicy::Unassessed,
    ] {
        let td = tempfile::tempdir().unwrap();
        let (_engine, binding, trace, mut report) =
            prepared(td.path(), &[masked_read()], policy, 12);
        assert!(accept(&report, &trace, &binding).is_ok());
        report.request.internal_primer_policy = if policy == PrimerPayloadInternalPrimerPolicy::Mask
        {
            PrimerPayloadInternalPrimerPolicy::Unassessed
        } else {
            PrimerPayloadInternalPrimerPolicy::Mask
        };
        persisted(&mut report);
        let error = accept(&report, &trace, &binding).unwrap_err();
        assert!(
            error.message.contains("exclusion_policy_binding"),
            "{error:?}"
        );
    }
}

#[test]
fn primer_payload_rehashed_missing_duplicate_and_foreign_exclusions_are_refused() {
    let td = tempfile::tempdir().unwrap();
    let (_engine, binding, trace, baseline) = prepared(
        td.path(),
        &[masked_read(), masked_read()],
        PrimerPayloadInternalPrimerPolicy::Unassessed,
        12,
    );
    assert!(accept(&baseline, &trace, &binding).is_ok());
    for case in ["missing", "duplicate", "foreign"] {
        let mut report = baseline.clone();
        if case == "missing" {
            report.unassessed.pop();
        } else {
            let mut row = report.unassessed[0].clone();
            if case == "foreign" {
                row.anchor_id = Some(format!("sha256:{}", "9".repeat(64)));
            }
            report.unassessed.push(row);
        }
        persisted(&mut report);
        assert!(accept(&report, &trace, &binding).is_err(), "{case}");
    }
}

#[test]
fn primer_payload_bound_import_record_budget_is_global_across_units() {
    let td = tempfile::tempdir().unwrap();
    let (engine, original) = source_reads(td.path(), true, &[format!("AGTCCGATCA{PAYLOAD}")]);
    let first =
        *GentleEngine::inspect_primer_trace(&original.report_path, Some(&original), None, 0)
            .unwrap()
            .report;
    let mut trace_request = first.request.clone();
    let mut second = trace_request.files[0].clone();
    let input = td.path().join("second input.fastq");
    let read = format!("AGTCCGATCA{PAYLOAD}T");
    fs::write(
        &input,
        format!("@reused_name\n{read}\n+\n{}\n", "I".repeat(read.len())),
    )
    .unwrap();
    second.path = input.display().to_string();
    second.expected_sha256 = Some(bound_file(&input).sha256);
    second.expected_md5 = None;
    trace_request.files.push(second);
    trace_request.output_dir = td.path().join("two-unit trace").display().to_string();
    let trace = engine
        .trace_primer_reads(trace_request.clone(), &mut |_| true)
        .unwrap();
    let binding = GentleEngine::primer_trace_result_binding(&trace_request, &trace).unwrap();
    let mut request = prepare(binding.clone(), td.path());
    request.file_ids.push("file_0001".into());
    let report = engine
        .prepare_primer_payload_alignment_queries(request, &mut |_| true)
        .unwrap();
    assert_eq!(report.queries.len(), 2);
    let file = bound_file(&Path::new(&report.request.output_dir).join("report.json"));
    assert!(binding::queries(&file, &trace, &binding, &limits(), 0).is_ok());
    let mut small = limits();
    small.max_records = 1;
    let error = binding::queries(&file, &trace, &binding, &small, 0).unwrap_err();
    assert!(error.message.contains("trace_record_budget"), "{error:?}");
}

#[test]
fn primer_payload_anchor_identity_preserves_published_encoding() {
    let td = tempfile::tempdir().unwrap();
    let (_engine, binding, trace, report) = prepared(
        td.path(),
        &[format!("AGTCCGATCA{PAYLOAD}")],
        PrimerPayloadInternalPrimerPolicy::Mask,
        12,
    );
    let root = Path::new(&binding.report_path)
        .parent()
        .unwrap()
        .join("artifacts");
    let mut row = prepare::source_rows(
        &root,
        &trace.files_completed[0],
        &trace.primer_panel,
        &limits(),
    )
    .unwrap()
    .remove(0);
    row.record_index = 5;
    row.read_id = "abstract_read".into();
    let mut hit = row.anchor.clone().unwrap();
    // Pin the existing fingerprint encoding; the synthetic strand change is not source evidence.
    for (strand, expected) in [
        (
            PrimerReadStrand::Forward,
            "sha256:409ca3317ce4c39c5fecaef99971395e4b3cd12ed1559e39d37a711f6e6eec81",
        ),
        (
            PrimerReadStrand::ReverseComplement,
            "sha256:1e90921487798d28b27616b7448e84bc946e7fafc1bb2cd670ba5d9ec553a5fc",
        ),
    ] {
        hit.strand = strand;
        let got = prepare::anchor_identity(
            &format!("sha256:{}", "0".repeat(64)),
            "arbitrary_file",
            &row,
            &hit,
            &report.queries[0].primer_ids,
        )
        .unwrap();
        assert_eq!(got, expected);
    }
}

#[test]
fn primer_payload_rehashed_anchor_out_of_range_is_refused_before_arithmetic() {
    let td = tempfile::tempdir().unwrap();
    let (_engine, binding, trace, baseline) = prepared(
        td.path(),
        &[format!("AGTCCGATCA{PAYLOAD}")],
        PrimerPayloadInternalPrimerPolicy::Mask,
        12,
    );
    for strand in [
        PrimerReadStrand::Forward,
        PrimerReadStrand::ReverseComplement,
    ] {
        let mut report = baseline.clone();
        report.queries[0].anchor.strand = strand;
        report.queries[0].anchor.end = usize::MAX;
        persisted(&mut report);
        let error = accept(&report, &trace, &binding).unwrap_err();
        assert!(error.message.contains("anchor_geometry"), "{error:?}");
    }
}
