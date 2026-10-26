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
    // Synthetic count mutations below have no measured work history; preserve absence,
    // not stale counter values from the original screen fixture.
    report.index_work_profile = None;
    for file in &mut report.files_completed {
        file.work_profile = None;
    }
    if report.schema == RNA_PANEL_SCREEN_SCHEMA {
        report.inputs_not_completed =
            Some(report.request.files[report.files_completed.len()..].to_vec());
    }
    report.content_sha256.clear();
    report.content_sha256 = sha256_prefixed_bytes(&encode(report).unwrap());
    fs::write(path, serde_json::to_vec_pretty(report).unwrap()).unwrap();
}

fn expected_summary_for_request(
    summary: &RnaReadPanelSummary,
    request: &ExportRnaReadPanelSummaryRequest,
) -> RnaReadPanelSummary {
    let mut expected = summary.clone();
    expected.submitted_request_sha256 = Some(sha256_prefixed_bytes(&encode(request).unwrap()));
    expected.content_sha256.clear();
    expected.content_sha256 = sha256_prefixed_bytes(&encode(&expected).unwrap());
    expected
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
            sample.co_nominated_reads,
            sample.no_candidate_reads
        ),
        (4, 2, Some(1), 1)
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
            row("locus_0:0").co_nominated_reads
        ),
        (1, Some(1))
    );
    assert_eq!(
        (
            row("locus_0:1").single_candidate_reads,
            row("locus_0:1").co_nominated_reads
        ),
        (0, Some(1))
    );
    assert_eq!(row("locus_1:0").single_candidate_reads, 1);
    assert_eq!(row("locus_1:1").co_nominated_reads, Some(0));
    assert_eq!(summary.co_nomination_sets.as_ref().unwrap().len(), 1);
    assert_eq!(
        summary.co_nomination_sets.as_ref().unwrap()[0].member_ids,
        vec!["locus_0:0", "locus_0:1"]
    );
    assert_eq!(summary.co_nomination_sets.as_ref().unwrap()[0].reads, 1);
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
        let expected = expected_summary_for_request(&summary, &repeated);
        assert_eq!(encode(&same).unwrap(), encode(&expected).unwrap());
        assert_ne!(
            same.submitted_request_sha256,
            summary.submitted_request_sha256
        );
        GentleEngine::rna_read_panel_summary_result_binding(&repeated, &same).unwrap();
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
    assert_eq!(summary.co_nomination_sets.as_ref().unwrap().len(), 2);
    let mut reversed = request.clone();
    reversed.report_paths.reverse();
    reversed.output_dir = td.path().join("reverse order").display().to_string();
    let repeated = engine.export_rna_read_panel_summary(&reversed).unwrap();
    let expected = expected_summary_for_request(&summary, &reversed);
    assert_eq!(encode(&expected).unwrap(), encode(&repeated).unwrap());
    assert_ne!(
        repeated.submitted_request_sha256,
        summary.submitted_request_sha256
    );
    GentleEngine::rna_read_panel_summary_result_binding(&reversed, &repeated).unwrap();
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
    assert_eq!(summary.co_nomination_sets.as_ref().unwrap()[0].reads, 2);
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
    let coverage = &partial.sample_coverage.as_ref().unwrap()[0];
    assert_eq!(coverage.status, "partial");
    assert_eq!(coverage.reads, Some(4));
    assert_eq!(coverage.inputs_not_completed[0].sample_id, "sample_a");
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
    let coverage = &empty.sample_coverage.as_ref().unwrap()[0];
    assert_eq!(coverage.status, "unassessed");
    assert_eq!(coverage.reads, None);
    let table =
        fs::read_to_string(Path::new(&request.output_dir).join("artifacts/samples.tsv")).unwrap();
    assert!(table.contains("\t0\t2\tunassessed\t\t\t\t\t\n"), "{table}");
}

#[test]
fn rna_panel_summary_names_an_unscreened_second_library_without_zero_reads() {
    let td = tempfile::tempdir().unwrap();
    let (engine, mut report, request) = fixture(td.path(), "sample_a");
    report.status = "cancelled".into();
    let mut missing = report.request.files[0].clone();
    missing.sample_id = "sample_b".into();
    missing.barcode = "barcode_b".into();
    missing.path.push_str(".not-yet-screened");
    report.request.files.push(missing.clone());
    report.files_not_completed = vec![missing.path.clone()];
    save_report(&request.report_paths[0], &mut report);
    let summary = engine.export_rna_read_panel_summary(&request).unwrap();
    let coverage = summary.sample_coverage.as_ref().unwrap();
    assert_eq!(coverage.len(), 2);
    assert_eq!(coverage[0].status, "complete");
    assert_eq!(coverage[0].reads, Some(4));
    assert_eq!(coverage[1].sample_id, "sample_b");
    assert_eq!(coverage[1].barcode, "barcode_b");
    assert_eq!(coverage[1].status, "unassessed");
    assert_eq!(coverage[1].reads, None);
    assert_eq!(coverage[1].inputs_not_completed[0].path, missing.path);
    assert_eq!(summary.completeness, "partial");
    assert!(
        summary
            .transcripts
            .iter()
            .all(|t| t.sample_id == "sample_a")
    );
}

#[test]
fn rna_panel_summary_preserves_joint_class_member_counts_and_legacy_uncertainty() {
    let td = tempfile::tempdir().unwrap();
    let (engine, mut report, mut request) = fixture(td.path(), "sample_a");
    let ids = vec!["locus_0:0".to_string(), "locus_1:0".to_string()];
    let file = &mut report.files_completed[0];
    file.classes = BTreeMap::from([
        ("multi_locus_ambiguous".into(), 2),
        ("multi_segment_candidate".into(), 2),
    ]);
    file.equivalence_classes = file
        .classes
        .iter()
        .map(|(class, count)| {
            (
                serde_json::to_string(&RnaReadPanelNominationSetKey {
                    class: class.clone(),
                    member_ids: ids.clone(),
                })
                .unwrap(),
                *count,
            )
        })
        .collect();
    save_report(&request.report_paths[0], &mut report);
    let summary = engine.export_rna_read_panel_summary(&request).unwrap();
    let sets = summary.co_nomination_sets.as_ref().unwrap();
    assert_eq!(sets.len(), 2);
    assert_eq!(sets[0].member_ids, sets[1].member_ids);
    assert_ne!(sets[0].class, sets[1].class);
    assert_eq!(sets.iter().map(|s| s.reads).sum::<u64>(), 4);
    assert_eq!(summary.class_member_counts_complete, Some(true));
    // V1 deliberately lost the association. Converting its stored counts must never reconstruct it.
    report.schema = RNA_PANEL_SCREEN_LEGACY_SCHEMA.into();
    report.inputs_not_completed = None;
    report.files_completed[0].equivalence_classes =
        BTreeMap::from([(serde_json::to_string(&ids).unwrap(), 4)]);
    save_report(&request.report_paths[0], &mut report);
    request.output_dir = td.path().join("legacy").display().to_string();
    let legacy = engine.export_rna_read_panel_summary(&request).unwrap();
    assert_eq!(legacy.class_member_counts_complete, Some(false));
    let set = &legacy.co_nomination_sets.as_ref().unwrap()[0];
    assert_eq!(set.class, "legacy_class_unresolved");
    assert_eq!(set.reads, 4);
    assert_eq!(set.member_ids, ids);
    assert_eq!(legacy.samples[0].classes["multi_segment_candidate"], 2);
}

#[test]
fn rna_panel_summary_legacy_json_keeps_digest_and_never_acquires_new_claims() {
    let td = tempfile::tempdir().unwrap();
    let (engine, _, request) = fixture(td.path(), "sample_a");
    let summary = engine.export_rna_read_panel_summary(&request).unwrap();
    let mut old = serde_json::to_value(summary).unwrap();
    old["schema"] = RNA_PANEL_SUMMARY_LEGACY_SCHEMA.into();
    let sets = old
        .as_object_mut()
        .unwrap()
        .remove("co_nomination_sets")
        .unwrap();
    old["ambiguous_sets"] = sets;
    for set in old["ambiguous_sets"].as_array_mut().unwrap() {
        set.as_object_mut().unwrap().remove("class");
        set.as_object_mut().unwrap().remove("sequence_sharing");
    }
    old.as_object_mut()
        .unwrap()
        .remove("shared_transcript_sequences");
    old.as_object_mut()
        .unwrap()
        .remove("source_overlap_assessment");
    for name in ["sample_coverage", "class_member_counts_complete"] {
        old.as_object_mut().unwrap().remove(name);
    }
    for row in old["samples"].as_array_mut().unwrap() {
        let count = row
            .as_object_mut()
            .unwrap()
            .remove("co_nominated_reads")
            .unwrap();
        row["ambiguous_reads"] = count;
    }
    for row in old["transcripts"].as_array_mut().unwrap() {
        let count = row
            .as_object_mut()
            .unwrap()
            .remove("co_nominated_reads")
            .unwrap();
        row["ambiguous_compatible_reads"] = count;
    }
    for run in old["source_runs"].as_array_mut().unwrap() {
        run.as_object_mut().unwrap().remove("inputs_not_completed");
    }
    old["content_sha256"] = "".into();
    old["content_sha256"] = sha256_prefixed_bytes(&encode(&old).unwrap()).into();
    let path = td.path().join("legacy-summary.json");
    for newline in ["\n", "\r\n"] {
        fs::write(
            &path,
            serde_json::to_string_pretty(&old)
                .unwrap()
                .replace('\n', newline),
        )
        .unwrap();
        let shown = engine
            .inspect_rna_read_panel(&path.to_string_lossy(), None, 0)
            .unwrap()
            .3
            .unwrap();
        assert_eq!(encode(&shown).unwrap(), encode(&old).unwrap());
        assert!(shown.co_nomination_sets.is_none());
        assert!(shown.class_member_counts_complete.is_none());
        assert!(shown.transcripts[0].co_nominated_reads.is_none());
    }
}

#[test]
fn rna_panel_summary_verification_is_explicit_bounded_and_hash_bound() {
    let td = tempfile::tempdir().unwrap();
    let (engine, _, request) = fixture(&td.path().join("with space and 'quote"), "sample_a");
    let summary = engine.export_rna_read_panel_summary(&request).unwrap();
    let path = Path::new(&request.output_dir).join("report.json");
    let path = path.to_str().unwrap();
    let verified = engine
        .verify_rna_read_panel_summary_artifacts(path, &summary)
        .unwrap();
    assert_eq!(verified.summary_sha256, summary.content_sha256);
    assert_eq!(
        verified.verified_bytes,
        summary.artifacts.iter().map(|a| a.bytes).sum::<u64>()
    );
    assert!(!verified.original_inputs_reverified);
    let target = Path::new(&request.output_dir).join("artifacts/transcripts.tsv");
    let bytes = fs::read(&target).unwrap();
    let mut edited = bytes.clone();
    edited[0] = b'X';
    fs::write(&target, &edited).unwrap();
    assert!(engine.inspect_rna_read_panel(path, None, 0).is_ok());
    assert!(
        engine
            .verify_rna_read_panel_summary_artifacts(path, &summary)
            .unwrap_err()
            .message
            .contains("artifact_hash_or_size")
    );
    fs::write(
        &target,
        String::from_utf8(bytes.clone())
            .unwrap()
            .replace('\n', "\r\n"),
    )
    .unwrap();
    assert!(
        engine.inspect_rna_read_panel(path, None, 0).is_ok(),
        "default show remains metadata-only"
    );
    assert!(
        engine
            .verify_rna_read_panel_summary_artifacts(path, &summary)
            .is_err()
    );
    fs::write(&target, &bytes).unwrap();
    assert!(
        engine
            .verify_rna_read_panel_summary_artifacts(path, &summary)
            .is_ok()
    );
    let mut traversal = summary.clone();
    traversal.artifacts[0].path = "artifacts/../report.json".into();
    traversal.content_sha256.clear();
    traversal.content_sha256 = sha256_prefixed_bytes(&encode(&traversal).unwrap());
    assert!(
        engine
            .verify_rna_read_panel_summary_artifacts(path, &traversal)
            .unwrap_err()
            .message
            .contains("allowlist")
    );
    let mut excessive = summary.clone();
    excessive.artifacts[0].bytes = 256_000_001;
    excessive.content_sha256.clear();
    excessive.content_sha256 = sha256_prefixed_bytes(&encode(&excessive).unwrap());
    assert!(
        engine
            .verify_rna_read_panel_summary_artifacts(path, &excessive)
            .unwrap_err()
            .message
            .contains("byte_budget")
    );
    #[cfg(unix)]
    {
        fs::remove_file(&target).unwrap();
        std::os::unix::fs::symlink(td.path().join("missing"), &target).unwrap();
        assert!(
            engine
                .verify_rna_read_panel_summary_artifacts(path, &summary)
                .unwrap_err()
                .message
                .contains("file_binding")
        );
    }
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
    let journal = engine.journal_len();
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
    assert_eq!(engine.journal_len(), journal);
    assert_eq!(encode(engine.state()).unwrap(), before);
}
