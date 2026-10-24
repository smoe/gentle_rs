//! Hand-crafted synthetic JSON and seeded toy loci/reads, never private or biological evidence.

use super::tests::{request, starter};
use super::*;

#[test]
fn rna_panel_legacy_request_bytes_and_cancel_identity_survive_optional_fields() {
    let bytes = br#"{"cancel_path":null,"files":[],"md5_manifest_path":null,"output_dir":"run","panel_path":"panel.json","policy":{"kmer_len":15,"max_anchors_per_read":100000,"max_equivalence_classes":4096,"max_index_postings":2000000,"min_chain_bases":45,"min_junction_flank_bp":5,"write_locus_fastas":false},"schema":"gentle.rna_read_panel_screen_request.v1"}"#;
    for text in [
        String::from_utf8(bytes.to_vec()).unwrap(),
        String::from_utf8(bytes.to_vec())
            .unwrap()
            .replace(',', ",\r\n"),
    ] {
        let request: ScreenRnaReadPanelRequest = serde_json::from_str(&text).unwrap();
        assert!(request.expected_panel_sha256.is_none());
        assert_eq!(encode(&request).unwrap(), bytes);
        assert_eq!(
            control::cancel_identity(&request).unwrap(),
            sha256_prefixed_bytes(bytes)
        );
    }
}

#[test]
fn rna_panel_reviewed_digest_rejects_another_valid_panel_before_streaming() {
    let td = tempfile::tempdir().unwrap();
    let (engine, mut build) = starter(2, td.path());
    let second = engine.build_rna_read_panel(&build).unwrap();
    build.sources.truncate(1);
    build.output_path = td.path().join("first.json").display().to_string();
    let first = engine.build_rna_read_panel(&build).unwrap();
    let mut request = request(&second, td.path(), &second.templates[0].sequence);
    request.expected_panel_sha256 = Some(first.content_sha256);
    let mut progress = 0;
    let error = engine
        .screen_rna_read_panel(request.clone(), &mut |_| {
            progress += 1;
            true
        })
        .unwrap_err();
    assert!(error.message.contains("panel_expected_digest_mismatch"));
    assert_eq!(progress, 0);
    assert!(!Path::new(&request.output_dir).exists());
}

#[test]
fn rna_panel_effective_binding_survives_manifest_changes_and_rejects_tampering() {
    let td = tempfile::tempdir().unwrap();
    let (mut engine, build) = starter(1, td.path());
    let panel = engine.build_rna_read_panel(&build).unwrap();
    let mut submitted = request(&panel, td.path(), &panel.templates[0].sequence);
    submitted.expected_panel_sha256 = Some(panel.content_sha256.clone());
    submitted.output_dir = td.path().join("run with spaces").display().to_string();
    let md5 = submitted.files[0].expected_md5.take().unwrap();
    let manifest = td.path().join("checksums.txt");
    fs::write(
        &manifest,
        format!("{}  reads.fastq.gz\n", md5.to_uppercase()),
    )
    .unwrap();
    submitted.md5_manifest_path = Some(manifest.display().to_string());
    let report = engine
        .screen_rna_read_panel(submitted.clone(), &mut |_| true)
        .unwrap();
    assert_eq!(
        report.request.files[0].expected_md5.as_deref(),
        Some(md5.as_str())
    );
    let binding = GentleEngine::rna_read_panel_screen_result_binding(&submitted, &report).unwrap();
    assert_ne!(
        binding.submitted_request_sha256,
        binding.effective_request_sha256
    );
    fs::write(&manifest, "replaced after publication").unwrap();
    fs::remove_file(&submitted.files[0].path).unwrap();
    let before = serde_json::to_value(engine.state()).unwrap();
    let line = format!(
        "rna-reads panel-show {} --binding {}",
        crate::engine_shell::quote_shell_arg(&binding.report_path),
        crate::engine_shell::quote_shell_arg(&serde_json::to_string(&binding).unwrap())
    );
    let command = crate::engine_shell::parse_shell_line(&line).unwrap();
    assert!(!command.is_state_mutating());
    let output = crate::engine_shell::execute_shell_command(&mut engine, &command).unwrap();
    assert_eq!(
        output.output["result"]["rna_read_panel_screen"]["content_sha256"],
        report.content_sha256
    );
    assert_eq!(serde_json::to_value(engine.state()).unwrap(), before);
    for field in [
        "submitted_request_sha256",
        "effective_request_sha256",
        "panel_sha256",
        "report_content_sha256",
        "report_path",
        "outcome",
    ] {
        let mut value = serde_json::to_value(&binding).unwrap();
        value[field] = serde_json::json!("different");
        let wrong: RnaReadPanelResultBinding = serde_json::from_value(value).unwrap();
        assert!(
            GentleEngine::validate_rna_read_panel_result_binding(
                &wrong,
                &binding.report_path,
                Some(&report),
                None
            )
            .is_err(),
            "{field}"
        );
    }
    let mut forged = report.clone();
    forged.request.policy.min_chain_bases += 1;
    forged.content_sha256.clear();
    forged.content_sha256 = sha256_prefixed_bytes(&encode(&forged).unwrap());
    assert!(GentleEngine::rna_read_panel_screen_result_binding(&submitted, &forged).is_err());
}

#[test]
fn rna_panel_summary_binding_and_legacy_report_hashes_are_explicit() {
    let td = tempfile::tempdir().unwrap();
    let (engine, build) = starter(1, td.path());
    let panel = engine.build_rna_read_panel(&build).unwrap();
    let submitted = request(&panel, td.path(), &panel.templates[0].sequence);
    let report = engine
        .screen_rna_read_panel(submitted.clone(), &mut |_| true)
        .unwrap();
    let request = ExportRnaReadPanelSummaryRequest {
        schema: RNA_PANEL_SUMMARY_REQUEST_SCHEMA.into(),
        report_paths: vec![
            Path::new(&submitted.output_dir)
                .join("report.json")
                .display()
                .to_string(),
        ],
        output_dir: td.path().join("summary").display().to_string(),
    };
    let summary = engine.export_rna_read_panel_summary(&request).unwrap();
    let binding = GentleEngine::rna_read_panel_summary_result_binding(&request, &summary).unwrap();
    GentleEngine::validate_rna_read_panel_result_binding(
        &binding,
        &binding.report_path,
        None,
        Some(&summary),
    )
    .unwrap();
    let mut wrong = request.clone();
    wrong.report_paths[0].push_str("-other");
    assert!(GentleEngine::rna_read_panel_summary_result_binding(&wrong, &summary).is_err());
    let mut legacy = report;
    legacy.submitted_request_sha256 = None;
    legacy.content_sha256.clear();
    legacy.content_sha256 = sha256_prefixed_bytes(&encode(&legacy).unwrap());
    let value = serde_json::to_value(&legacy).unwrap();
    assert!(value.get("submitted_request_sha256").is_none());
    let decoded: RnaReadPanelScreenReport = serde_json::from_value(value).unwrap();
    validate_screen_report_binding(&decoded).unwrap();
    assert_eq!(encode(&legacy).unwrap(), encode(&decoded).unwrap());
    assert!(
        GentleEngine::rna_read_panel_screen_result_binding(&submitted, &decoded).is_err(),
        "legacy provenance is not invented"
    );
    let mut legacy = summary;
    legacy.submitted_request_sha256 = None;
    legacy.content_sha256.clear();
    legacy.content_sha256 = sha256_prefixed_bytes(&encode(&legacy).unwrap());
    assert!(
        serde_json::to_value(&legacy)
            .unwrap()
            .get("submitted_request_sha256")
            .is_none()
    );
}
