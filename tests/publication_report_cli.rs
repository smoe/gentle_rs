//! CLI regressions for publication request dispatch.
//!
//! Test inputs are synthetic JSON written into a temporary directory; no
//! external or biological fixture data is used.

use std::{fs, process::Command};

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
