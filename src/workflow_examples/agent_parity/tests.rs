//! Synthetic, temporary tutorial contracts; no models or scientific operations.

use super::*;
use crate::workflow_examples::{TUTORIAL_AGENT_PARITY_SCHEMA, TutorialAgentParserPayload};
use serde_json::json;

fn unit(cases: serde_json::Value) -> TutorialSourceUnit {
    serde_json::from_value(json!({
        "id": "synthetic", "title": "Synthetic parser contract",
        "catalog": {"order": 1, "path": "guide.md", "type": "gui_cli_walkthrough",
                    "status": "manual/hybrid", "source": "hand_written_markdown"},
        "agent_parity": {"schema": TUTORIAL_AGENT_PARITY_SCHEMA, "cases": cases}
    }))
    .unwrap()
}

fn case(id: &str, command: &str, execution: &str, mutating: bool) -> serde_json::Value {
    json!({"id": id, "command": command, "execution": execution, "mutating": mutating})
}

fn write_guide(root: &Path, unit: &TutorialSourceUnit, extra: &str) {
    let mut text = unit
        .agent_parity
        .as_ref()
        .unwrap()
        .cases
        .iter()
        .map(|case| format!("```text\n{}\n```\n", case.command))
        .collect::<String>();
    text.push_str(extra);
    fs::write(root.join("guide.md"), text).unwrap();
}

#[test]
fn discovers_all_contracts_and_keeps_semantic_and_parser_counts_distinct() {
    let temp = tempfile::tempdir().unwrap();
    let first = unit(json!([
        case("inspect", "promoters tss-collection sample", "ask", false),
        case("ui", "ui open sequence-window sample", "auto", false),
        case("forget", "promoters tss-forget sample", "ask", true)
    ]));
    write_guide(temp.path(), &first, "");
    let mut another = first.clone();
    another.id = "new_tutorial_not_in_a_whitelist".into();
    let mut legacy = first.clone();
    legacy.agent_parity = None;
    let report = check_tutorial_agent_parity(&[first, another, legacy], temp.path());
    assert!(report.findings.is_empty(), "{:?}", report.findings);
    assert_eq!(
        report.summary,
        TutorialAgentParitySummary {
            tutorials: 2,
            cases: 6,
            declared_mutating: 2,
            parser_state_mutating: 4
        }
    );
}

#[test]
fn source_loader_still_rejects_author_declared_auto_mutation() {
    let temp = tempfile::tempdir().unwrap();
    let unit = unit(json!([case(
        "unsafe_auto",
        "promoters tss-forget sample",
        "auto",
        true
    )]));
    fs::write(
        temp.path().join("source.json"),
        serde_json::to_vec(&unit).unwrap(),
    )
    .unwrap();
    let error = crate::workflow_examples::load_tutorial_source_units(temp.path()).unwrap_err();
    assert!(error.contains("unsafe_auto"), "{error}");
    assert!(error.contains("must use execution='ask'"), "{error}");
}

#[test]
fn collects_stable_findings_in_source_order_including_history_auto() {
    let temp = tempfile::tempdir().unwrap();
    let unit = unit(json!([
        case("bad", "not-a-shell-command", "ask", false),
        case("misdeclared", "state-summary", "ask", true),
        case(
            "auto_inspect",
            "promoters tss-collection sample",
            "auto",
            false
        ),
        case("undo", "history undo", "auto", false),
        case("redo", "history redo", "auto", false)
    ]));
    fs::write(temp.path().join("guide.md"), "unrelated prose").unwrap();
    let report = check_tutorial_agent_parity(&[unit], temp.path());
    let codes: Vec<_> = report
        .findings
        .iter()
        .map(|finding| finding.code.as_str())
        .collect();
    assert_eq!(
        codes,
        [
            "guide_missing_command",
            "parse_failed",
            "guide_missing_command",
            "declared_mutation_not_parser_mutating",
            "guide_missing_command",
            "auto_blocked_by_runtime",
            "guide_missing_command",
            "auto_blocked_by_runtime",
            "guide_missing_command",
            "auto_blocked_by_runtime"
        ]
    );
    let message = report.failure_message().unwrap();
    assert!(message.starts_with("tutorial agent parity failed (10 findings):\n"));
    for line in message.lines().skip(1) {
        let finding: serde_json::Value = serde_json::from_str(line).unwrap();
        assert_eq!(finding["tutorial_id"], "synthetic");
        assert!(finding["case_id"].is_string());
    }
}

#[test]
fn payload_paths_use_repo_root_and_shared_quoting_not_cwd() {
    let temp = tempfile::Builder::new()
        .prefix("tutorial root's ")
        .tempdir()
        .unwrap();
    let unit = unit(json!([case(
        "inspect",
        "promoters tss-inventory '@request file.json'",
        "ask",
        false
    )]));
    write_guide(temp.path(), &unit, "");
    fs::write(
        temp.path().join("request file.json"),
        r#"{"seq_id":"sample's quoted id", "gene_query":"SYNTHETIC", "collection_id":"synthetic_tss"}"#,
    )
    .unwrap();
    let report = check_tutorial_agent_parity(&[unit], temp.path());
    assert!(report.findings.is_empty(), "{:?}", report.findings);
    assert_eq!(fs::read_dir(temp.path()).unwrap().count(), 2);
}

#[test]
fn learner_template_matches_semantic_json_and_keeps_unapproved_placeholders() {
    for newline in ["\n", "\r\n"] {
        let temp = tempfile::tempdir().unwrap();
        let mut unit = unit(json!([case(
            "materialize",
            "promoters tss-materialize @learner.json",
            "ask",
            true
        )]));
        let case = &mut unit.agent_parity.as_mut().unwrap().cases[0];
        case.parser_payload = Some(TutorialAgentParserPayload {
            file: "learner.json".into(),
            template: "template.json".into(),
        });
        let value = json!({
            "inventory": {"seq_id": "synthetic", "gene_query": "SYNTHETIC", "collection_id": "synthetic_tss"},
            "expected_approval_sha256": "sha256:PASTE_FROM_YOUR_PREVIEW",
            "selected_tss_ids": ["tss_PASTE_FROM_YOUR_PREVIEW"]
        });
        fs::write(
            temp.path().join("template.json"),
            serde_json::to_string_pretty(&value)
                .unwrap()
                .replace('\n', newline),
        )
        .unwrap();
        write_guide(
            temp.path(),
            &unit,
            &format!("```json\n{value}\n```\n").replace('\n', newline),
        );
        let report = check_tutorial_agent_parity(&[unit.clone()], temp.path());
        assert!(report.findings.is_empty(), "{:?}", report.findings);
        let guide = fs::read_to_string(temp.path().join("guide.md")).unwrap();
        let command = parser_command(
            &unit.agent_parity.as_ref().unwrap().cases[0],
            &temp.path().canonicalize().unwrap(),
            &guide,
        )
        .unwrap();
        let ShellCommand::Op { payload } = parse_shell_line(&command).unwrap() else {
            panic!("expected operation");
        };
        assert!(payload.contains("sha256:PASTE_FROM_YOUR_PREVIEW"));
        assert!(payload.contains("tss_PASTE_FROM_YOUR_PREVIEW"));
        assert!(!temp.path().join("learner.json").exists());
    }
}

#[test]
fn rejects_template_drift_bad_json_and_unused_mapping() {
    let temp = tempfile::tempdir().unwrap();
    let mut unit = unit(json!([case(
        "preview",
        "promoters tss-inventory @learner.json",
        "ask",
        false
    )]));
    unit.agent_parity.as_mut().unwrap().cases[0].parser_payload =
        Some(TutorialAgentParserPayload {
            file: "learner.json".into(),
            template: "template.json".into(),
        });
    write_guide(
        temp.path(),
        &unit,
        "```json\n{\"seq_id\":\"original\"}\n```\n",
    );
    for (template, expected) in [
        ("{\"seq_id\":\"changed\"}", "payload_not_in_guide"),
        ("{broken", "payload_invalid_json"),
    ] {
        fs::write(temp.path().join("template.json"), template).unwrap();
        let report = check_tutorial_agent_parity(&[unit.clone()], temp.path());
        assert_eq!(report.findings[0].code, expected);
    }
    unit.agent_parity.as_mut().unwrap().cases[0]
        .parser_payload
        .as_mut()
        .unwrap()
        .file = "unused.json".into();
    let report = check_tutorial_agent_parity(&[unit], temp.path());
    assert_eq!(report.findings[0].code, "payload_binding_invalid");
}

#[test]
fn invalid_typed_json_is_not_made_valid_by_a_matching_template() {
    let temp = tempfile::tempdir().unwrap();
    let mut unit = unit(json!([case(
        "preview",
        "promoters tss-inventory @learner.json",
        "ask",
        false
    )]));
    unit.agent_parity.as_mut().unwrap().cases[0].parser_payload =
        Some(TutorialAgentParserPayload {
            file: "learner.json".into(),
            template: "template.json".into(),
        });
    fs::write(temp.path().join("template.json"), "{}").unwrap();
    write_guide(temp.path(), &unit, "```json\n{}\n```\n");
    let report = check_tutorial_agent_parity(&[unit], temp.path());
    assert_eq!(report.findings[0].code, "parse_failed");
}

#[test]
fn rejects_raw_absolute_traversal_and_windows_paths_before_normalizing() {
    let temp = tempfile::tempdir().unwrap();
    fs::write(temp.path().join("request.json"), "{}").unwrap();
    fs::create_dir(temp.path().join("child")).unwrap();
    for path in [
        "../request.json",
        "child/../request.json",
        "/request.json",
        "/tmp/request.json",
        "//server/share/request.json",
        "C:request.json",
        "C:/request.json",
        r"C:\request.json",
        r"\request.json",
        r"\\server\request.json",
    ] {
        let command = format!(
            "promoters tss-inventory {}",
            shell_quote(&format!("@{path}"))
        );
        let unit = unit(json!([case("unsafe", &command, "ask", false)]));
        write_guide(temp.path(), &unit, "");
        let report = check_tutorial_agent_parity(&[unit], temp.path());
        assert_eq!(
            report.findings[0].code, "unsafe_path",
            "{path}: {:?}",
            report.findings
        );
    }
}

#[test]
fn rejects_rooted_guide_template_and_payload_paths_before_file_lookup() {
    let temp = tempfile::tempdir().unwrap();
    let mut original = unit(json!([case(
        "preview",
        "promoters tss-inventory @learner.json",
        "ask",
        false
    )]));
    original.agent_parity.as_mut().unwrap().cases[0].parser_payload =
        Some(TutorialAgentParserPayload {
            file: "learner.json".into(),
            template: "template.json".into(),
        });
    let payload =
        json!({"seq_id": "synthetic", "gene_query": "SYNTHETIC", "collection_id": "synthetic_tss"});
    fs::write(temp.path().join("template.json"), payload.to_string()).unwrap();
    let extra = format!("```json\n{payload}\n```\n");
    for raw in [
        "/missing.json",
        "//server/share/missing.json",
        r"\missing.json",
        "C:missing.json",
    ] {
        for surface in ["guide", "template", "payload"] {
            let mut invalid = original.clone();
            match surface {
                "guide" => invalid.catalog.as_mut().unwrap().path = raw.into(),
                "template" => {
                    invalid.agent_parity.as_mut().unwrap().cases[0]
                        .parser_payload
                        .as_mut()
                        .unwrap()
                        .template = raw.into();
                }
                "payload" => {
                    let case = &mut invalid.agent_parity.as_mut().unwrap().cases[0];
                    case.parser_payload.as_mut().unwrap().file = raw.into();
                    case.command = format!(
                        "promoters tss-inventory {}",
                        shell_quote(&format!("@{raw}"))
                    );
                }
                _ => unreachable!(),
            }
            write_guide(temp.path(), &invalid, &extra);
            let report = check_tutorial_agent_parity(&[invalid], temp.path());
            assert_eq!(report.findings.len(), 1, "{surface}: {raw}: {report:?}");
            assert_eq!(
                report.findings[0].code, "unsafe_path",
                "{surface}: {raw}: {report:?}"
            );
        }
    }
}

#[test]
fn guide_template_and_missing_payload_use_the_same_bound_file_policy() {
    let temp = tempfile::tempdir().unwrap();
    let mut unit = unit(json!([case(
        "preview",
        "promoters tss-inventory @missing.json",
        "ask",
        false
    )]));
    write_guide(temp.path(), &unit, "");
    let report = check_tutorial_agent_parity(&[unit.clone()], temp.path());
    assert_eq!(report.findings[0].code, "file_unavailable");
    unit.catalog.as_mut().unwrap().path = "../guide.md".into();
    let report = check_tutorial_agent_parity(&[unit.clone()], temp.path());
    assert_eq!(report.findings[0].code, "unsafe_path");
    unit.catalog.as_mut().unwrap().path = "guide.md".into();
    unit.agent_parity.as_mut().unwrap().cases[0].parser_payload =
        Some(TutorialAgentParserPayload {
            file: "missing.json".into(),
            template: "../template.json".into(),
        });
    let report = check_tutorial_agent_parity(&[unit], temp.path());
    assert_eq!(report.findings[0].code, "unsafe_path");
}

#[cfg(unix)]
#[test]
fn rejects_symlink_escape_but_allows_in_repository_symlinks() {
    let temp = tempfile::tempdir().unwrap();
    let outside = tempfile::tempdir().unwrap();
    fs::write(
        outside.path().join("request.json"),
        r#"{"seq_id":"synthetic"}"#,
    )
    .unwrap();
    std::os::unix::fs::symlink(
        outside.path().join("request.json"),
        temp.path().join("escape.json"),
    )
    .unwrap();
    let mut unit = unit(json!([case(
        "preview",
        "promoters tss-inventory @escape.json",
        "ask",
        false
    )]));
    write_guide(temp.path(), &unit, "");
    let report = check_tutorial_agent_parity(&[unit.clone()], temp.path());
    assert_eq!(report.findings[0].code, "unsafe_path");
    fs::write(
        temp.path().join("local.json"),
        r#"{"seq_id":"synthetic", "gene_query":"SYNTHETIC", "collection_id":"synthetic_tss"}"#,
    )
    .unwrap();
    std::os::unix::fs::symlink("local.json", temp.path().join("link.json")).unwrap();
    unit.agent_parity.as_mut().unwrap().cases[0].command =
        "promoters tss-inventory @link.json".into();
    write_guide(temp.path(), &unit, "");
    assert!(
        check_tutorial_agent_parity(&[unit], temp.path())
            .findings
            .is_empty()
    );
}
