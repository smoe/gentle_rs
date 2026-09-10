//! CLI-only presentation of exported TSS profiles; engine/MCP results stay complete.

use gentle::engine::{OpResult, Operation};
use gentle::engine_shell::ShellCommand;
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::path::Path;

pub(super) fn export_directory(op: &Operation) -> Option<String> {
    match op {
        Operation::ComputeTssTfbsProfiles { export, .. } => {
            export.as_ref().map(|request| request.output_dir.clone())
        }
        Operation::ExportTssTfbsProfiles { request, .. } => Some(request.output_dir.clone()),
        _ => None,
    }
}

// Read only the destination from the externally tagged operation. In particular,
// do not deserialize a second copy of a report-only export's score vectors.
#[derive(Deserialize)]
struct ExportProjection {
    #[serde(rename = "ComputeTssTfbsProfiles")]
    compute: Option<ComputeProjection>,
    #[serde(rename = "ExportTssTfbsProfiles")]
    export: Option<ExportOnlyProjection>,
}

#[derive(Deserialize)]
struct ComputeProjection {
    export: Option<Directory>,
}

#[derive(Deserialize)]
struct ExportOnlyProjection {
    request: Directory,
}

#[derive(Deserialize)]
struct Directory {
    output_dir: String,
}

pub(super) fn shell_export_directory(command: &ShellCommand) -> Option<String> {
    let ShellCommand::Op { payload } = command else {
        return None;
    };
    let text = super::load_json_arg(payload).ok()?;
    let projected: ExportProjection = serde_json::from_str(&text).ok()?;
    projected
        .compute
        .and_then(|op| op.export)
        .or_else(|| projected.export.map(|op| op.request))
        .map(|directory| directory.output_dir)
}

fn summarize_result(result: &mut Value, output_dir: &str) {
    let Some(receipt) = result
        .get("tss_tfbs_profile_receipt")
        .filter(|v| v.is_object())
    else {
        // Never discard a compute-only result: there is no saved report to read.
        return;
    };
    let report = &result["tss_tfbs_profiles"];
    let file = |name: &str| {
        json!({
            "path": Path::new(output_dir).join(name),
            "sha256": receipt["outputs"][name],
        })
    };
    let summary = json!({
        "schema": "gentle.tss_tfbs_profile_cli_summary.v1",
        "output_dir": output_dir,
        "report": file("report.json"),
        "index": file("index.json"),
        // The receipt deliberately does not include its own hash.
        "receipt_path": Path::new(output_dir).join("receipt.json"),
        "tss_count": receipt["tss_count"],
        "page_count": receipt["page_count"],
        "hashed_output_file_count": receipt["outputs"].as_object().map(|files| files.len()),
        "input_manifest_sha256": receipt["input_manifest_sha256"],
        "source_revision": receipt["source_revision"],
        "producer_revision": receipt["producer_revision"],
        "exporter_revision": receipt["exporter_revision"],
        "verification": report.get("verification").cloned()
            .unwrap_or(json!("not_reassessed_report_only_export")),
        "reference": report["reference"],
        "report_warnings": report.get("warnings"),
        "non_claims": receipt["non_claims"],
    });
    let object = result
        .as_object_mut()
        .expect("a result containing a receipt is an object");
    object.remove("tss_tfbs_profiles");
    object.remove("tss_tfbs_profile_receipt");
    object.insert("tss_tfbs_profile_summary".into(), summary);
}

pub(super) fn shell_output(
    mut output: Value,
    output_dir: Option<&str>,
    full_report: bool,
) -> Value {
    if !full_report {
        if let (Some(result), Some(directory)) = (output.get_mut("result"), output_dir) {
            summarize_result(result, directory);
        }
    }
    output
}

#[derive(Serialize)]
#[serde(untagged)]
pub(super) enum OperationOutput {
    Full(Box<OpResult>),
    Summary(Value),
}

pub(super) fn operation_output(
    mut result: OpResult,
    output_dir: Option<&str>,
    full_report: bool,
) -> Result<OperationOutput, String> {
    let Some(directory) =
        output_dir.filter(|_| !full_report && result.tss_tfbs_profile_receipt.is_some())
    else {
        return Ok(OperationOutput::Full(Box::new(result)));
    };
    // Take ownership so we can omit millions of scores before serialization,
    // without cloning them or maintaining a second list of base result fields.
    let report = result.tss_tfbs_profiles.take();
    let mut projected = serde_json::to_value(result).map_err(|e| e.to_string())?;
    if let Some(report) = report {
        projected["tss_tfbs_profiles"] = json!({
            "verification": report.verification,
            "reference": report.reference,
            "warnings": report.warnings,
        });
    }
    summarize_result(&mut projected, directory);
    Ok(OperationOutput::Summary(projected))
}

#[cfg(test)]
mod tests {
    use super::*;
    use gentle::engine::{Engine, GentleEngine, ProjectState};
    use gentle::engine_shell::{
        ShellExecutionOptions, execute_shell_command_with_options, parse_shell_tokens,
    };

    fn operation_json(result: &OpResult, directory: Option<&str>, full: bool) -> Value {
        serde_json::to_value(operation_output(result.clone(), directory, full).unwrap()).unwrap()
    }

    fn compute_command(output: Option<&Path>) -> ShellCommand {
        let fixture =
            Path::new(env!("CARGO_MANIFEST_DIR")).join("test_files/fixtures/tss_profiles");
        let tokens = [
            "features",
            "tss-tfbs-profiles",
            "--manifest",
            fixture.join("manifest.json").to_str().unwrap(),
            "--panel",
            fixture.join("panel.json").to_str().unwrap(),
            "--expected-genome-id",
            "synthetic-genome-v1",
            "--output-dir",
            output.unwrap_or(Path::new("unused")).to_str().unwrap(),
        ]
        .map(str::to_string);
        let mut command = parse_shell_tokens(&tokens).unwrap();
        if output.is_none() {
            let ShellCommand::Op { payload } = &mut command else {
                panic!("typed operation")
            };
            let mut op: Operation = serde_json::from_str(payload).unwrap();
            let Operation::ComputeTssTfbsProfiles { export, .. } = &mut op else {
                panic!("compute")
            };
            *export = None;
            *payload = serde_json::to_string(&op).unwrap();
        }
        command
    }

    #[test]
    fn tss_export_stdout_is_compact_and_bound_without_changing_shared_results() {
        let temp = tempfile::tempdir().unwrap();
        let root = temp.path().canonicalize().unwrap();
        let destination = root.join("export");
        let command = compute_command(Some(&destination));
        let directory = shell_export_directory(&command).unwrap();
        let mut engine = GentleEngine::from_state(ProjectState::default());
        let run = execute_shell_command_with_options(
            &mut engine,
            &command,
            &ShellExecutionOptions::default(),
        )
        .unwrap();
        assert!(!run.state_changed);
        let original = run.output;
        let full = shell_output(original.clone(), Some(&directory), true);
        assert_eq!(full, original);
        let report = &original["result"]["tss_tfbs_profiles"];
        assert!(!report["windows"].as_array().unwrap().is_empty());
        let saved: Value =
            serde_json::from_slice(&std::fs::read(destination.join("report.json")).unwrap())
                .unwrap();
        assert_eq!(&saved, report);

        let concise = shell_output(original.clone(), Some(&directory), false);
        let summary = &concise["result"]["tss_tfbs_profile_summary"];
        let receipt = &original["result"]["tss_tfbs_profile_receipt"];
        assert_eq!(summary["report"]["sha256"], receipt["report_sha256"]);
        assert_eq!(summary["verification"], report["verification"]);
        assert!(
            summary["verification"]
                .as_str()
                .unwrap()
                .contains("prepared_reference_not_assessed")
        );
        assert_eq!(summary["reference"], report["reference"]);
        assert_eq!(summary["report_warnings"], report["warnings"]);
        assert_eq!(summary["non_claims"], receipt["non_claims"]);
        assert_eq!(
            summary["receipt_path"],
            destination.join("receipt.json").to_str().unwrap()
        );
        assert!(concise["result"].get("tss_tfbs_profiles").is_none());
        let compact_len = serde_json::to_vec_pretty(&concise).unwrap().len();
        assert!(compact_len < 8_000, "stdout is {compact_len} bytes");
        let mut large = original.clone();
        large["result"]["tss_tfbs_profiles"]["windows"] =
            json!(vec![json!({"scores": vec![0.25; 10_000]}); 10]);
        assert_eq!(shell_output(large, Some(&directory), false), concise);

        let result: OpResult = serde_json::from_value(original["result"].clone()).unwrap();
        assert_eq!(
            operation_json(&result, Some(&directory), false),
            concise["result"]
        );
        assert_eq!(
            operation_json(&result, Some(&directory), true),
            original["result"]
        );

        let reexport = Operation::ExportTssTfbsProfiles {
            report: result.tss_tfbs_profiles.unwrap(),
            request: gentle_protocol::tss_profiles::ExportTssProfilesRequest {
                output_dir: root.join("reexport").to_str().unwrap().into(),
                rendering: Default::default(),
                formats: vec![gentle_protocol::tss_profiles::TssExportFormat::Svg],
            },
        };
        let directory = export_directory(&reexport).unwrap();
        assert_eq!(
            shell_export_directory(&ShellCommand::Op {
                payload: serde_json::to_string(&reexport).unwrap()
            }),
            Some(directory.clone())
        );
        let result = engine.apply(reexport).unwrap();
        let summary = operation_json(&result, Some(&directory), false);
        assert_eq!(
            summary["tss_tfbs_profile_summary"]["verification"],
            "not_reassessed_report_only_export"
        );
        assert!(summary["tss_tfbs_profile_summary"]["report_warnings"].is_null());
    }

    #[test]
    fn tss_compute_without_export_keeps_the_only_full_report() {
        let command = compute_command(None);
        assert!(shell_export_directory(&command).is_none());
        let ShellCommand::Op { payload } = command else {
            panic!("typed operation")
        };
        let op: Operation = serde_json::from_str(&payload).unwrap();
        assert!(export_directory(&op).is_none());
        let result = GentleEngine::from_state(ProjectState::default())
            .apply(op)
            .unwrap();
        let expected = serde_json::to_value(&result).unwrap();
        assert_eq!(operation_json(&result, None, false), expected);
        // Even a stale destination hint must not hide the only copy of the data.
        assert_eq!(operation_json(&result, Some("unused"), false), expected);
        let wrapped = json!({"result": expected});
        assert_eq!(
            shell_output(wrapped.clone(), Some("unused"), false),
            wrapped
        );
    }

    #[test]
    fn tss_presentation_leaves_unrelated_results_unchanged() {
        let original =
            json!({"result": {"warnings": ["keep me"], "other_report": {"scores": [1, 2]}}});
        assert_eq!(
            shell_output(original.clone(), Some("unused"), false),
            original
        );
        assert!(
            shell_export_directory(&ShellCommand::Op {
                payload: "{\"StateSummary\":{}}".into()
            })
            .is_none()
        );
    }
}
