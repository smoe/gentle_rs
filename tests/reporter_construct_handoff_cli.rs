use serde_json::Value;
use std::process::Command;

#[test]
fn reporter_plan_handoff_cli_writes_output_schema_and_macro_command() {
    let temp = tempfile::tempdir().expect("tempdir");
    let output = temp.path().join("handoff.json");
    let output_result = Command::new(env!("CARGO_BIN_EXE_gentle_cli"))
        .args([
            "reporters",
            "plan-handoff",
            "docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/promoter_reporter_candidates.json",
            "--output",
            output.to_str().expect("utf-8 output path"),
        ])
        .output()
        .expect("run gentle_cli reporters plan-handoff");
    assert!(
        output_result.status.success(),
        "gentle_cli reporters plan-handoff failed: {}",
        String::from_utf8_lossy(&output_result.stderr)
    );

    let payload: Value = serde_json::from_str(
        &std::fs::read_to_string(&output).expect("read reporter handoff output"),
    )
    .expect("parse handoff JSON");
    assert_eq!(
        payload["schema"].as_str(),
        Some("gentle.reporter_construct_handoff.v1")
    );
    let commands = payload["commands"]
        .as_array()
        .expect("commands array in handoff output");
    assert!(commands.iter().any(|command| {
        command["command"]
            .as_str()
            .unwrap_or_default()
            .contains("allele_paired_promoter_luciferase_reporter")
    }));
}
