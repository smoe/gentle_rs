//! Real CLI admission tests using hand-authored temporary source units.
//! See tutorial_agent_parity_README.md for provenance and recreation.

use serde_json::json;
use std::{fs, process::Command};

#[test]
fn tutorial_check_reports_parity_findings_before_workflow_loading_or_execution() {
    for newline in ["\n", "\r\n"] {
        let root = tempfile::Builder::new()
            .prefix("tutorial repo's ")
            .tempdir()
            .unwrap();
        let cwd = tempfile::tempdir().unwrap();
        let sources = root.path().join("docs/tutorial/sources");
        fs::create_dir_all(&sources).unwrap();
        let source = json!({
            "schema": "gentle.tutorial_source.v4", "id": "boundary", "title": "Synthetic boundary",
            "catalog": {"order": 1, "path": "guide.md", "type": "gui_cli_walkthrough",
                        "status": "manual/hybrid", "source": "hand_written_markdown"},
            "agent_parity": {"schema": "gentle.tutorial_agent_parity.v1", "cases": [
                {"id": "bad_grammar", "command": "not-a-shell-command", "execution": "ask", "mutating": false},
                {"id": "auto_forget", "command": "promoters tss-forget never_execute", "execution": "auto", "mutating": false}
            ]}
        });
        fs::write(
            sources.join("boundary.json"),
            serde_json::to_string_pretty(&source)
                .unwrap()
                .replace('\n', newline),
        )
        .unwrap();
        fs::write(
            root.path().join("guide.md"),
            "```text\nnot-a-shell-command\npromoters tss-forget never_execute\n```\n"
                .replace('\n', newline),
        )
        .unwrap();
        // Missing replay inputs deliberately prove that admission runs first.
        let before_source = fs::read(sources.join("boundary.json")).unwrap();
        let before_guide = fs::read(root.path().join("guide.md")).unwrap();
        let output = Command::new(env!("CARGO_BIN_EXE_gentle_examples_docs"))
            .args(["tutorial-check", "--repo-root"])
            .arg(root.path())
            .args([
                "--manifest",
                "docs/tutorial/manifest.json",
                "--source",
                "missing-workflows",
                "--tutorial-output",
                "never-generated",
            ])
            .current_dir(cwd.path())
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(1), "{output:?}");
        assert!(output.stdout.is_empty());
        let stderr = String::from_utf8(output.stderr).unwrap();
        assert!(
            stderr.contains("tutorial agent parity failed (2 findings):"),
            "{stderr}"
        );
        let findings: Vec<serde_json::Value> = stderr
            .lines()
            .filter_map(|line| serde_json::from_str(line).ok())
            .collect();
        assert_eq!(findings.len(), 2, "{stderr}");
        assert_eq!(findings[0]["code"], "parse_failed");
        assert_eq!(findings[0]["tutorial_id"], "boundary");
        assert_eq!(findings[0]["case_id"], "bad_grammar");
        assert_eq!(findings[1]["code"], "auto_blocked_by_runtime");
        assert_eq!(findings[1]["case_id"], "auto_forget");
        assert_eq!(fs::read_dir(cwd.path()).unwrap().count(), 0);
        assert!(!root.path().join("never-generated").exists());
        assert_eq!(fs::read_dir(root.path()).unwrap().count(), 2);
        assert_eq!(
            fs::read(sources.join("boundary.json")).unwrap(),
            before_source
        );
        assert_eq!(
            fs::read(root.path().join("guide.md")).unwrap(),
            before_guide
        );
    }
}
