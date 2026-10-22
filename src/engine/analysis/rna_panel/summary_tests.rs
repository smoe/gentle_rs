//! Hand-crafted fixed-LCG loci and literal FASTQ only; recreated in temporary directories.
//! These checks exercise report aggregation, not biological abundance or real-data acceptance.

use super::tests::{bases, fastq, gzip, request, starter};
use super::*;
use md5::Digest;

fn fixture(
    directory: &Path,
    sample: &str,
) -> (
    GentleEngine,
    RnaReadPanelScreenReport,
    ExportRnaReadPanelSummaryRequest,
) {
    fs::create_dir_all(directory).unwrap();
    let (engine, build) = starter(2, directory);
    let panel = engine.build_rna_read_panel(&build).unwrap();
    let sequence = bases(100, 150);
    let other = bases(101, 150);
    let mut screen = request(&panel, directory, &sequence[..40]);
    let mut bytes = vec![];
    let unmatched = "N".repeat(40);
    for (i, read) in [
        &sequence[..40],
        &sequence[60..100],
        &other[60..100],
        &unmatched,
    ]
    .into_iter()
    .enumerate()
    {
        let record = String::from_utf8(fastq(read, b'>')).unwrap();
        bytes.extend_from_slice(
            record
                .replace("@read same-header-may-recur", &format!("@{sample}_{i}"))
                .as_bytes(),
        );
    }
    let bytes = gzip(&bytes);
    fs::write(&screen.files[0].path, &bytes).unwrap();
    screen.files[0].sample_id = sample.into();
    screen.files[0].expected_md5 = Some(input::hex(&md5::Md5::digest(&bytes)));
    screen.policy.kmer_len = 15;
    screen.policy.min_chain_bases = 30;
    screen.policy.min_junction_flank_bp = 5;
    let report = engine.screen_rna_read_panel(screen, &mut |_| true).unwrap();
    let summary = ExportRnaReadPanelSummaryRequest {
        schema: RNA_PANEL_SUMMARY_REQUEST_SCHEMA.into(),
        report_paths: vec![directory.join("run/report.json").display().to_string()],
        output_dir: directory.join("summary").display().to_string(),
    };
    (engine, report, summary)
}

fn save_report(path: &str, report: &mut RnaReadPanelScreenReport) {
    report.content_sha256.clear();
    report.content_sha256 = sha256_prefixed_bytes(&encode(report).unwrap());
    fs::write(path, serde_json::to_vec_pretty(report).unwrap()).unwrap();
}

#[test]
fn rna_panel_summary_counts_reads_without_resolving_ambiguous_isoforms_or_reopening_inputs() {
    let td = tempfile::tempdir().unwrap();
    let (engine, report, request) = fixture(td.path(), "sample_a");
    let before = encode(engine.state()).unwrap();
    fs::remove_file(&report.files_completed[0].input.path).unwrap();
    fs::remove_file(td.path().join("run/artifacts/file_0000/nominations.tsv")).unwrap();
    let summary = engine.export_rna_read_panel_summary(&request).unwrap();
    let sample = &summary.samples[0];
    assert_eq!(
        (
            sample.reads,
            sample.single_candidate_reads,
            sample.ambiguous_reads,
            sample.no_candidate_reads
        ),
        (4, 2, 1, 1)
    );
    assert_eq!(summary.completeness, "complete");
    assert!(!summary.original_inputs_reverified);
    assert!(!summary.nomination_sidecars_reverified);
    assert_eq!(summary.transcripts.len(), 4);
    let row = |member: &str| {
        summary
            .transcripts
            .iter()
            .find(|t| t.member_id == member)
            .unwrap()
    };
    assert_eq!(
        (
            row("locus_0:0").single_candidate_reads,
            row("locus_0:0").ambiguous_compatible_reads
        ),
        (1, 1)
    );
    assert_eq!(
        (
            row("locus_0:1").single_candidate_reads,
            row("locus_0:1").ambiguous_compatible_reads
        ),
        (0, 1)
    );
    assert_eq!(row("locus_1:0").single_candidate_reads, 1);
    assert_eq!(row("locus_1:1").ambiguous_compatible_reads, 0);
    assert_eq!(summary.ambiguous_sets.len(), 1);
    assert_eq!(
        summary.ambiguous_sets[0].member_ids,
        vec!["locus_0:0", "locus_0:1"]
    );
    assert_eq!(summary.ambiguous_sets[0].reads, 1);
    assert_eq!(encode(engine.state()).unwrap(), before);
    for artifact in &summary.artifacts {
        let bytes = fs::read(Path::new(&request.output_dir).join(&artifact.path)).unwrap();
        assert_eq!(sha256_prefixed_bytes(&bytes), artifact.sha256);
        assert_eq!(bytes.len() as u64, artifact.bytes);
        assert!(!bytes.contains(&b'\r'));
    }
    let stored = Path::new(&request.output_dir).join("report.json");
    let shown = engine
        .inspect_rna_read_panel(&stored.to_string_lossy(), None, 0)
        .unwrap()
        .3
        .unwrap();
    assert_eq!(encode(&shown).unwrap(), encode(&summary).unwrap());
    assert!(
        engine
            .inspect_rna_read_panel(&stored.to_string_lossy(), None, 1)
            .is_err()
    );
    assert!(engine.export_rna_read_panel_summary(&request).is_err());
    let report_path = &request.report_paths[0];
    let source = String::from_utf8(fs::read(report_path).unwrap()).unwrap();
    for (label, bytes) in [
        ("lf", source.clone()),
        ("crlf", source.replace('\n', "\r\n")),
    ] {
        fs::write(report_path, bytes).unwrap();
        let mut repeated = request.clone();
        repeated.output_dir = td.path().join(label).display().to_string();
        let same = engine.export_rna_read_panel_summary(&repeated).unwrap();
        assert_eq!(encode(&same).unwrap(), encode(&summary).unwrap());
    }
    let mut tampered = serde_json::to_value(&summary).unwrap();
    tampered["samples"][0]["reads"] = 999.into();
    fs::write(&stored, serde_json::to_vec(&tampered).unwrap()).unwrap();
    assert!(
        engine
            .inspect_rna_read_panel(&stored.to_string_lossy(), None, 0)
            .is_err()
    );
}

#[test]
fn rna_panel_summary_combines_samples_deterministically_and_requires_identical_policy() {
    let td = tempfile::tempdir().unwrap();
    let (engine, _, mut request) = fixture(&td.path().join("first"), "sample_a");
    let (_, mut other, second) = fixture(&td.path().join("second"), "sample_b");
    request.report_paths.extend(second.report_paths.clone());
    let summary = engine.export_rna_read_panel_summary(&request).unwrap();
    assert_eq!(summary.samples.len(), 2);
    assert_eq!(summary.transcripts.len(), 8);
    assert_eq!(summary.ambiguous_sets.len(), 2);
    let mut reversed = request.clone();
    reversed.report_paths.reverse();
    reversed.output_dir = td.path().join("reverse order").display().to_string();
    let repeated = engine.export_rna_read_panel_summary(&reversed).unwrap();
    assert_eq!(encode(&summary).unwrap(), encode(&repeated).unwrap());
    other.request.policy.min_chain_bases += 1;
    save_report(&second.report_paths[0], &mut other);
    request.output_dir = td.path().join("incompatible policy").display().to_string();
    let error = engine.export_rna_read_panel_summary(&request).unwrap_err();
    assert!(error.message.contains("incompatible_panel_or_policy"));
    assert!(!Path::new(&request.output_dir).exists());
}

#[test]
fn rna_panel_summary_aggregates_chunks_but_refuses_duplicate_inputs_and_conflicting_sample_bindings()
 {
    let td = tempfile::tempdir().unwrap();
    let (engine, first, mut request) = fixture(&td.path().join("a"), "sample_a");
    let (_, mut second, second_request) = fixture(&td.path().join("b"), "sample_b");
    second.request.files[0].sample_id = "sample_a".into();
    second.files_completed[0].input.sample_id = "sample_a".into();
    save_report(&second_request.report_paths[0], &mut second);
    request
        .report_paths
        .extend(second_request.report_paths.clone());
    let summary = engine.export_rna_read_panel_summary(&request).unwrap();
    assert_eq!(summary.samples.len(), 1);
    assert_eq!(summary.samples[0].reads, 8);
    assert_eq!(summary.ambiguous_sets[0].reads, 2);
    request.output_dir = td.path().join("duplicate input").display().to_string();
    second.files_completed[0].observed_sha256 = first.files_completed[0].observed_sha256.clone();
    second.files_completed[0].observed_md5 = first.files_completed[0].observed_md5.clone();
    second.files_completed[0].input.expected_md5 =
        first.files_completed[0].input.expected_md5.clone();
    second.request.files[0] = second.files_completed[0].input.clone();
    save_report(&second_request.report_paths[0], &mut second);
    assert!(
        engine
            .export_rna_read_panel_summary(&request)
            .unwrap_err()
            .message
            .contains("duplicate_input")
    );
    assert!(!Path::new(&request.output_dir).exists());
    second.files_completed[0].observed_sha256 = format!("sha256:{}", "1".repeat(64));
    second.files_completed[0].input.barcode = "different_barcode".into();
    second.request.files[0] = second.files_completed[0].input.clone();
    save_report(&second_request.report_paths[0], &mut second);
    assert!(
        engine
            .export_rna_read_panel_summary(&request)
            .unwrap_err()
            .message
            .contains("sample_binding_conflict")
    );
}

#[test]
fn rna_panel_summary_cancelled_sources_stay_partial_including_zero_verified_units() {
    let td = tempfile::tempdir().unwrap();
    let (engine, mut report, mut request) = fixture(td.path(), "sample_a");
    report.status = "cancelled".into();
    let mut incomplete = report.request.files[0].clone();
    incomplete.path.push_str(".not-completed");
    report.request.files.push(incomplete.clone());
    report.files_not_completed = vec![incomplete.path];
    save_report(&request.report_paths[0], &mut report);
    let partial = engine.export_rna_read_panel_summary(&request).unwrap();
    assert_eq!(partial.completeness, "partial");
    assert_eq!(partial.samples[0].reads, 4);
    report.files_completed.clear();
    report.files_not_completed = report
        .request
        .files
        .iter()
        .map(|f| f.path.clone())
        .collect();
    save_report(&request.report_paths[0], &mut report);
    request.output_dir = td.path().join("zero verified").display().to_string();
    let empty = engine.export_rna_read_panel_summary(&request).unwrap();
    assert_eq!(empty.completeness, "partial");
    assert!(empty.samples.is_empty());
    assert!(
        empty.transcripts.is_empty(),
        "unassessed samples must not acquire zero-count rows"
    );
    assert_eq!(empty.source_runs[0].files_not_completed.len(), 2);
}

#[test]
fn rna_panel_summary_rejects_tampering_unknown_members_and_inconsistent_counts_atomically() {
    let td = tempfile::tempdir().unwrap();
    let (engine, report, request) = fixture(td.path(), "sample_a");
    let mut invalid = report.clone();
    invalid.files_completed[0].reads += 1;
    fs::write(
        &request.report_paths[0],
        serde_json::to_vec(&invalid).unwrap(),
    )
    .unwrap();
    assert!(
        engine
            .export_rna_read_panel_summary(&request)
            .unwrap_err()
            .message
            .contains("panel_report_binding")
    );
    save_report(&request.report_paths[0], &mut invalid);
    assert!(
        engine
            .export_rna_read_panel_summary(&request)
            .unwrap_err()
            .message
            .contains("class_total")
    );
    for key in [
        "[\"unknown:0\"]",
        "[\"locus_0:1\",\"locus_0:0\"]",
        "[\"locus_0:0\",\"locus_0:0\"]",
        "[ \"locus_0:0\" ]",
    ] {
        let mut invalid = report.clone();
        let count = invalid.files_completed[0]
            .equivalence_classes
            .pop_first()
            .unwrap()
            .1;
        invalid.files_completed[0]
            .equivalence_classes
            .insert(key.into(), count);
        save_report(&request.report_paths[0], &mut invalid);
        assert!(
            engine.export_rna_read_panel_summary(&request).is_err(),
            "{key}"
        );
        assert!(!Path::new(&request.output_dir).exists());
    }
    let mut overflow = report.clone();
    overflow.files_completed[0].classes =
        BTreeMap::from([("partial".into(), u64::MAX), ("unmatched".into(), 1)]);
    save_report(&request.report_paths[0], &mut overflow);
    assert!(
        engine
            .export_rna_read_panel_summary(&request)
            .unwrap_err()
            .message
            .contains("count_overflow")
    );
}

#[test]
fn rna_panel_summary_rejects_duplicate_reports_tampered_panels_and_invalid_requests() {
    let td = tempfile::tempdir().unwrap();
    let (engine, _, request) = fixture(&td.path().join("input"), "sample_a");
    for schema in ["wrong", RNA_PANEL_SUMMARY_REQUEST_SCHEMA] {
        let mut invalid = request.clone();
        invalid.schema = schema.into();
        invalid.report_paths.clear();
        assert!(engine.export_rna_read_panel_summary(&invalid).is_err());
    }
    let mut invalid = request.clone();
    invalid.report_paths = vec![request.report_paths[0].clone(); 129];
    assert!(engine.export_rna_read_panel_summary(&invalid).is_err());
    invalid.report_paths = vec![request.report_paths[0].clone(); 2];
    assert!(
        engine
            .export_rna_read_panel_summary(&invalid)
            .unwrap_err()
            .message
            .contains("duplicate_report")
    );
    invalid = request.clone();
    invalid.output_dir = td.path().join("input/../escape").display().to_string();
    assert!(engine.export_rna_read_panel_summary(&invalid).is_err());
    let path = td.path().join("input/run/artifacts/panel.json");
    let mut panel: RnaReadPanel = load(&path.to_string_lossy()).unwrap();
    panel.templates[0].transcript_id = "tampered".into();
    fs::write(path, serde_json::to_vec(&panel).unwrap()).unwrap();
    assert!(engine.export_rna_read_panel_summary(&request).is_err());
    assert!(!Path::new(&request.output_dir).exists());
}

#[test]
fn rna_panel_summary_shared_shell_reopens_summary_without_mutating_project_or_history_on_show() {
    let td = tempfile::tempdir().unwrap();
    let (mut engine, _, request) = fixture(&td.path().join("with space and 'quote"), "sample_a");
    let before = encode(engine.state()).unwrap();
    let command = crate::engine_shell::parse_shell_line(&format!(
        "rna-reads panel-export-summary {}",
        crate::engine_shell::quote_shell_arg(&serde_json::to_string(&request).unwrap())
    ))
    .unwrap();
    assert!(
        command.is_state_mutating(),
        "external writes must keep confirmation"
    );
    let result = crate::engine_shell::execute_shell_command(&mut engine, &command).unwrap();
    assert!(!result.state_changed);
    assert_eq!(
        result.output["result"]["rna_read_panel_summary"]["samples"][0]["reads"],
        4
    );
    let journal = engine.journal().len();
    let show = crate::engine_shell::parse_shell_line(&format!(
        "rna-reads panel-show {}",
        crate::engine_shell::quote_shell_arg(
            &Path::new(&request.output_dir)
                .join("report.json")
                .to_string_lossy()
        )
    ))
    .unwrap();
    assert!(!show.is_state_mutating());
    let shown = crate::engine_shell::execute_shell_command(&mut engine, &show).unwrap();
    assert_eq!(
        shown.output["result"]["rna_read_panel_summary"],
        result.output["result"]["rna_read_panel_summary"]
    );
    assert_eq!(engine.journal().len(), journal);
    assert_eq!(encode(engine.state()).unwrap(), before);
}
