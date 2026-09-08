//! CLI regressions for publication request dispatch.
//!
//! Test inputs are synthetic JSON written into a temporary directory; no
//! external or biological fixture data is used.

use std::{fs, process::Command};

#[test]
fn publication_help_succeeds_without_inputs_or_outputs() {
    let temp = tempfile::tempdir().expect("tempdir");
    let binary = env!("CARGO_BIN_EXE_gentle_publication_report");
    for flag in ["--help", "-h"] {
        let output = Command::new(binary)
            .arg(flag)
            .current_dir(temp.path())
            .output()
            .expect("publication help");
        assert!(output.status.success(), "{output:?}");
        assert!(
            String::from_utf8_lossy(&output.stdout).contains("Usage: gentle_publication_report")
        );
        assert!(output.stderr.is_empty());
    }
    let missing_inputs = Command::new(binary)
        .current_dir(temp.path())
        .output()
        .expect("publication without inputs");
    assert_eq!(missing_inputs.status.code(), Some(2));
    assert_eq!(fs::read_dir(temp.path()).expect("read tempdir").count(), 0);
}

#[test]
fn legacy_publication_requests_reject_isoform_projection_options() {
    let temp = tempfile::tempdir().expect("tempdir");
    let request_path = temp.path().join("legacy-request.json");
    fs::write(
        &request_path,
        serde_json::to_vec_pretty(&serde_json::json!({
            "schema": "gentle.gene_set_publication_request.v1"
        }))
        .unwrap(),
    )
    .expect("legacy request");
    let binary = env!("CARGO_BIN_EXE_gentle_publication_report");

    for options in [
        vec!["--profile", "review"],
        vec!["--blocks", "gene.gene1.overview"],
    ] {
        let output = Command::new(binary)
            .arg(&request_path)
            .arg(temp.path().join("out"))
            .args(options)
            .output()
            .expect("run publication CLI");
        assert!(!output.status.success());
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(
            stderr.contains(
                "--profile and --blocks require gentle.gene_isoform_assay_publication_request.v1"
            ),
            "{stderr}"
        );
    }
}
