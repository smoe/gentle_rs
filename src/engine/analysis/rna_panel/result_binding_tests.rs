//! Hand-crafted synthetic JSON and seeded toy loci/reads, never private or biological evidence.

use super::tests::{request, starter};
use super::*;

#[test]
fn rna_panel_screen_draft_roundtrips_without_opening_inputs_or_overwriting_files() {
    // Hand-crafted invalid-for-execution draft; none of these input paths exists.
    let mut engine = GentleEngine::default();
    let request: ScreenRnaReadPanelRequest = serde_json::from_value(serde_json::json!({
        "schema": RNA_PANEL_SCREEN_REQUEST_SCHEMA, "panel_path": "absent-panel.json",
        "files": [], "md5_manifest_path": "absent-manifest", "output_dir": "absent-output"
    }))
    .unwrap();
    let before = serde_json::to_value(engine.state()).unwrap();
    let journal = engine.journal_len();
    let td = tempfile::tempdir().unwrap();
    let path = td.path().join("draft with 'quotes'.json");
    let line = format!(
        "rna-reads panel-screen --request {} --save {}",
        crate::engine_shell::quote_shell_arg(&serde_json::to_string(&request).unwrap()),
        crate::engine_shell::quote_shell_arg(&path.display().to_string())
    );
    let command = crate::engine_shell::parse_shell_line(&line).unwrap();
    assert!(
        command.is_state_mutating(),
        "external write needs agent approval"
    );
    let saved = crate::engine_shell::execute_shell_command(&mut engine, &command).unwrap();
    assert!(!saved.state_changed);
    let bytes = fs::read(&path).unwrap();
    assert_eq!(
        saved.output["artifact"]["sha256"],
        sha256_prefixed_bytes(&bytes)
    );
    assert!(crate::engine_shell::execute_shell_command(&mut engine, &command).is_err());
    assert_eq!(fs::read(&path).unwrap(), bytes);
    for text in [
        String::from_utf8(bytes).unwrap(),
        serde_json::to_string_pretty(&request)
            .unwrap()
            .replace('\n', "\r\n"),
        format!(
            "#!/usr/bin/env gentle\r\n{}",
            serde_json::to_string_pretty(&request)
                .unwrap()
                .replace('\n', "\r\n")
        ),
    ] {
        fs::write(&path, text).unwrap();
        let load = format!(
            "rna-reads panel-screen --request {}",
            crate::engine_shell::quote_shell_arg(&format!("@{}", path.display()))
        );
        let command = crate::engine_shell::parse_shell_line(&load).unwrap();
        assert!(!command.is_state_mutating());
        let result = crate::engine_shell::execute_shell_command(&mut engine, &command).unwrap();
        assert_eq!(
            result.output["request"],
            serde_json::to_value(&request).unwrap()
        );
        assert!(result.output["artifact"].is_null());
    }
    assert_eq!(serde_json::to_value(engine.state()).unwrap(), before);
    assert_eq!(engine.journal_len(), journal);
    assert!(!td.path().join("absent-output").exists());
    assert!(
        engine
            .rna_read_panel_screen_draft(
                &request,
                Some(&td.path().join("../escape.json").display().to_string())
            )
            .is_err()
    );
}

#[test]
fn rna_panel_screen_draft_file_is_lazy_bounded_and_does_not_open_referenced_inputs() {
    // Synthetic metadata, including a sparse oversized file: no biological reads or large allocation.
    let td = tempfile::tempdir().unwrap();
    let path = td.path().join("lazy 'draft' \u{00fc}.json");
    let line = format!(
        "rna-reads panel-screen --request {}",
        crate::engine_shell::quote_shell_arg(&format!("@{}", path.display()))
    );
    let command = crate::engine_shell::parse_shell_line(&line).unwrap();
    assert!(matches!(
        command,
        crate::engine_shell::ShellCommand::RnaReadPanelLoadDraft { .. }
    ));
    assert!(!command.is_state_mutating());
    let mut engine = GentleEngine::default();
    let before = serde_json::to_value(engine.state()).unwrap();
    assert!(
        crate::engine_shell::execute_shell_command(&mut engine, &command)
            .unwrap_err()
            .contains("panel_draft_open")
    );
    let request = serde_json::json!({
        "schema": RNA_PANEL_SCREEN_REQUEST_SCHEMA, "panel_path": "absent-panel.json",
        "files": [], "md5_manifest_path": "absent-manifest", "output_dir": "absent-output"
    });
    fs::write(&path, serde_json::to_vec(&request).unwrap()).unwrap();
    let captured = crate::engine_shell::execute_shell_command(&mut engine, &command).unwrap();
    assert_eq!(
        captured.output["request"]["panel_path"],
        request["panel_path"]
    );
    assert!(captured.output["artifact"].is_null());
    fs::write(&path, b"{\"unknown_field\":true}").unwrap();
    assert!(
        crate::engine_shell::execute_shell_command(&mut engine, &command)
            .unwrap_err()
            .contains("panel_draft_json")
    );
    let file = fs::File::create(&path).unwrap();
    file.set_len(128_000_001).unwrap();
    drop(file);
    assert!(
        crate::engine_shell::execute_shell_command(&mut engine, &command)
            .unwrap_err()
            .contains("panel_draft_json_budget")
    );
    assert_eq!(
        captured.output["request"]["panel_path"],
        request["panel_path"]
    );
    assert_eq!(serde_json::to_value(engine.state()).unwrap(), before);
    assert_eq!(engine.journal_len(), 0);
}

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
