//! Real CLI/shared-shell expression drafts with synthetic temporary inputs.
//! See protein_expression_handoff_README.md for provenance and recreation.

use gentle::{dna_sequence::DNAsequence, engine::ProjectState, engine_shell::quote_shell_arg};
use serde_json::{Value, json};
use std::{
    collections::BTreeMap,
    fs,
    path::{Path, PathBuf},
    process::Command,
};

fn snapshot(root: &Path) -> BTreeMap<PathBuf, Option<Vec<u8>>> {
    fn visit(root: &Path, path: &Path, files: &mut BTreeMap<PathBuf, Option<Vec<u8>>>) {
        for entry in fs::read_dir(path).unwrap() {
            let path = entry.unwrap().path();
            let relative = path.strip_prefix(root).unwrap().to_path_buf();
            if path.is_dir() {
                files.insert(relative, None);
                visit(root, &path, files);
            } else {
                files.insert(relative, Some(fs::read(path).unwrap()));
            }
        }
    }
    let mut files = BTreeMap::new();
    visit(root, root, &mut files);
    files
}

#[test]
fn protein_expression_handoff_cli_is_read_only_with_explicit_provider_and_catalog_gates() {
    for newline in ["\n", "\r\n"] {
        for case in [
            "neutral",
            "selected",
            "protein",
            "in_house",
            "excluded",
            "missing_product",
            "missing_permission",
            "partial",
            "unsupported_capability",
        ] {
            let temp = tempfile::Builder::new()
                .prefix("expression repo's ")
                .tempdir()
                .unwrap();
            let root = temp.path();
            let mut sequence = DNAsequence::from_sequence(if case == "protein" {
                "MHHHHHHA"
            } else {
                "ATGCACCACCACCACCACCACGCTTAA"
            })
            .unwrap();
            if case == "protein" {
                sequence.set_molecule_type("protein");
            }
            let mut state = ProjectState::default();
            state
                .sequences
                .insert("synthetic_product".to_string(), sequence);
            let state_path = root.join("state's project.gentle.json");
            state.save_to_path(&state_path.to_string_lossy()).unwrap();
            let mut requirements: Value = serde_json::from_str(include_str!(
                "../docs/examples/plans/protein_expression_requirements_in_house.json"
            ))
            .unwrap();
            requirements["outsourcing"]["allowed"] = json!(case != "in_house");
            requirements["outsourcing"]["approved_providers"] = json!([if case == "excluded" {
                "metabion"
            } else {
                "geneart"
            }]);
            if case == "missing_permission" {
                requirements.as_object_mut().unwrap().remove("outsourcing");
            }
            if case == "partial" {
                requirements = json!({"schema":"gentle.protein_expression_requirements.v1", "yield_goal":{"metric":"active_protein"}, "outsourcing":{"allowed":true}});
            }
            let requirements_path = root.join("reviewed requirements.json");
            fs::write(
                &requirements_path,
                serde_json::to_string_pretty(&requirements)
                    .unwrap()
                    .replace('\n', newline),
            )
            .unwrap();
            for directory in ["home", "config", "system", "runtime", "tmp"] {
                fs::create_dir(root.join(directory)).unwrap();
            }
            if case == "unsupported_capability" {
                // An overlay replaces the provider, not just its display name.
                let mut catalog: Value =
                    serde_json::from_str(include_str!("../assets/external_service_providers.json"))
                        .unwrap();
                let mut provider = catalog["providers"][0].clone();
                provider["capabilities"]
                    .as_array_mut()
                    .unwrap()
                    .retain(|capability| capability["service_kind"] != "protein_expression");
                catalog["providers"] = json!([provider]);
                let directory = root.join(".gentle/catalogs");
                fs::create_dir_all(&directory).unwrap();
                fs::write(
                    directory.join("external_service_providers.json"),
                    serde_json::to_string_pretty(&catalog)
                        .unwrap()
                        .replace('\n', newline),
                )
                .unwrap();
            }
            let mut tokens = vec![
                "planning".to_string(),
                "protein-expression-handoff".to_string(),
                "--seq-id".to_string(),
                if case == "missing_product" {
                    "absent"
                } else {
                    "synthetic_product"
                }
                .to_string(),
                "--requirements".to_string(),
                format!("@{}", requirements_path.display()),
                "--format".to_string(),
                "json".to_string(),
            ];
            if case != "neutral" {
                tokens.extend(["--provider".to_string(), "geneart".to_string()]);
            }
            let before = snapshot(root);
            let mut outputs = Vec::new();
            for shell in [false, true] {
                let mut command = Command::new(env!("CARGO_BIN_EXE_gentle_cli"));
                command
                    .env_remove("RUST_MIN_STACK")
                    .env("GENTLE_ASSET_ROOT", env!("CARGO_MANIFEST_DIR"))
                    .env("GENTLE_SYSTEM_CONFIG_ROOT", root.join("system"))
                    .env("GENTLE_PROJECT_ROOT", root)
                    .env("HOME", root.join("home"))
                    .env("XDG_CONFIG_HOME", root.join("config"))
                    .env("XDG_RUNTIME_DIR", root.join("runtime"))
                    .env("TMPDIR", root.join("tmp"))
                    .env("TMP", root.join("tmp"))
                    .env("TEMP", root.join("tmp"))
                    .current_dir(root)
                    .arg("--state")
                    .arg(&state_path);
                if shell {
                    command.arg("shell").arg(
                        tokens
                            .iter()
                            .map(|token| quote_shell_arg(token))
                            .collect::<Vec<_>>()
                            .join(" "),
                    );
                } else {
                    command.args(&tokens);
                }
                let output = command.output().unwrap();
                assert!(
                    output.status.success(),
                    "{case}: {}",
                    String::from_utf8_lossy(&output.stderr)
                );
                let mut report: Value = serde_json::from_slice(&output.stdout).unwrap();
                assert_eq!(report["schema"], "gentle.protein_expression_handoff.v2");
                assert_eq!(
                    snapshot(root),
                    before,
                    "{case}: no input rewrites or new files/directories"
                );
                let expected = match case {
                    "neutral" => "not_selected",
                    "selected" | "protein" | "partial" => "projected",
                    _ => "withheld",
                };
                assert_eq!(report["provider_binding"]["status"], expected, "{case}");
                if expected != "projected" {
                    assert!(report["provider_projection"].is_null());
                    assert!(
                        report["service_handoff_candidates"]
                            .as_array()
                            .unwrap()
                            .is_empty()
                    );
                    assert!(
                        report["suggested_next_actions"]
                            .as_array()
                            .unwrap()
                            .iter()
                            .all(|action| !action["shell_line"]
                                .as_str()
                                .unwrap()
                                .starts_with("services "))
                    );
                } else {
                    let identity = if case == "protein" {
                        "protein_seq_id"
                    } else {
                        "seq_id"
                    };
                    assert_eq!(
                        report["provider_projection"]["source_target"][identity],
                        "synthetic_product"
                    );
                }
                if case == "partial" {
                    assert!(report["neutral_service_draft"]["optimization_target"].is_null());
                    assert!(!report["missing_questions"].as_array().unwrap().is_empty());
                }
                // Time of report generation is the sole expected adapter difference.
                report
                    .as_object_mut()
                    .unwrap()
                    .remove("generated_at_unix_ms");
                outputs.push(report);
            }
            assert_eq!(
                outputs[0], outputs[1],
                "{case}: direct CLI/shared shell parity"
            );
        }
    }
}
