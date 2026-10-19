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
            "quote_unavailable",
            "blocked_preflight",
            "unknown_provider",
            "invalid_product",
        ] {
            let temp = tempfile::Builder::new()
                .prefix("expression repo's ")
                .tempdir()
                .unwrap();
            let root = temp.path();
            let mut sequence = DNAsequence::from_sequence(if case == "protein" {
                "MHHHHHHA"
            } else if case == "invalid_product" {
                "ACGNNNT"
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
            if case == "unknown_provider" {
                requirements["outsourcing"]["approved_providers"] = json!([]);
            }
            let requirements_path = root.join("reviewed requirements.json");
            fs::write(
                &requirements_path,
                serde_json::to_string_pretty(&requirements)
                    .unwrap()
                    .replace('\n', newline),
            )
            .unwrap();
            for directory in [
                "home", "config", "cache", "data", "system", "runtime", "tmp",
            ] {
                fs::create_dir(root.join(directory)).unwrap();
            }
            if matches!(
                case,
                "unsupported_capability" | "quote_unavailable" | "blocked_preflight"
            ) {
                // An overlay replaces the provider, not just its display name.
                let mut catalog: Value =
                    serde_json::from_str(include_str!("../assets/external_service_providers.json"))
                        .unwrap();
                let mut provider = catalog["providers"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .find(|provider| provider["provider"] == "geneart")
                    .unwrap()
                    .clone();
                if case == "unsupported_capability" {
                    provider["capabilities"]
                        .as_array_mut()
                        .unwrap()
                        .retain(|capability| capability["service_kind"] != "protein_expression");
                } else if case == "quote_unavailable" {
                    let capability = provider["capabilities"]
                        .as_array_mut()
                        .unwrap()
                        .iter_mut()
                        .find(|capability| capability["service_kind"] == "protein_expression")
                        .unwrap();
                    capability["quote_handoff_supported"] = json!(false);
                    capability["direct_api_implemented"] = json!(false);
                } else {
                    let rule = provider["validation_rules"]
                        .as_array_mut()
                        .unwrap()
                        .iter_mut()
                        .find(|rule| rule["service_kind"] == "protein_expression")
                        .unwrap();
                    rule["required_source_fields"] = json!(["synthetic_missing_source_field"]);
                }
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
                if matches!(case, "selected" | "quote_unavailable" | "blocked_preflight") {
                    "text"
                } else {
                    "json"
                }
                .to_string(),
            ];
            if case != "neutral" {
                tokens.extend([
                    "--provider".to_string(),
                    if case == "unknown_provider" {
                        "synthetic_unknown_expression_provider"
                    } else {
                        "geneart"
                    }
                    .to_string(),
                ]);
            }
            let before = snapshot(root);
            let mut outputs = Vec::new();
            let make_command = || {
                let mut command = Command::new(env!("CARGO_BIN_EXE_gentle_cli"));
                for (key, _) in std::env::vars_os() {
                    if key.to_str().is_some_and(|key| key.starts_with("GENTLE_")) {
                        command.env_remove(key);
                    }
                }
                command
                    .env_remove("RUST_MIN_STACK")
                    .env("GENTLE_ASSET_ROOT", env!("CARGO_MANIFEST_DIR"))
                    .env("GENTLE_SYSTEM_CONFIG_ROOT", root.join("system"))
                    .env("GENTLE_PROJECT_ROOT", root)
                    .env("HOME", root.join("home"))
                    .env("USERPROFILE", root.join("home"))
                    .env("XDG_CONFIG_HOME", root.join("config"))
                    .env("XDG_CACHE_HOME", root.join("cache"))
                    .env("XDG_DATA_HOME", root.join("data"))
                    .env("XDG_RUNTIME_DIR", root.join("runtime"))
                    .env("TMPDIR", root.join("tmp"))
                    .env("TMP", root.join("tmp"))
                    .env("TEMP", root.join("tmp"))
                    .current_dir(root)
                    .arg("--state")
                    .arg(&state_path);
                command
            };
            for shell in [false, true] {
                let mut command = make_command();
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
                    "selected" | "protein" | "partial" | "quote_unavailable"
                    | "blocked_preflight" => "projected",
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
                                .starts_with("services project-"))
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
                    let request = serde_json::to_string(&report["provider_projection"]).unwrap();
                    let mut preflight = make_command();
                    if shell {
                        preflight.arg("shell").arg(format!(
                            "services project-preflight {}",
                            quote_shell_arg(&request)
                        ));
                    } else {
                        preflight.args(["services", "project-preflight", &request]);
                    }
                    let result = preflight.output().unwrap();
                    assert!(
                        result.status.success(),
                        "{case}: {}",
                        String::from_utf8_lossy(&result.stderr)
                    );
                    let readiness: Value = serde_json::from_slice(&result.stdout).unwrap();
                    let expected_eligible =
                        !matches!(case, "quote_unavailable" | "blocked_preflight");
                    assert_eq!(
                        readiness["eligible"], expected_eligible,
                        "{case}, shell={shell}, newline={newline:?}: {readiness}"
                    );
                    assert_eq!(
                        readiness["quote_handoff_available"], expected_eligible,
                        "{case}, shell={shell}, newline={newline:?}: {readiness}"
                    );
                    let blocking_issues = readiness["blocking_issues"].as_array().unwrap();
                    assert_eq!(
                        blocking_issues.is_empty(),
                        expected_eligible,
                        "{case}, shell={shell}, newline={newline:?}: {readiness}"
                    );
                    if !expected_eligible {
                        assert_eq!(readiness["capability_status"], "blocked", "{readiness}");
                        let reason = if case == "quote_unavailable" {
                            "no implemented quote handoff or submission adapter"
                        } else {
                            "synthetic_missing_source_field"
                        };
                        assert!(
                            blocking_issues
                                .iter()
                                .any(|issue| issue.as_str().unwrap().contains(reason)),
                            "{case}, shell={shell}, newline={newline:?}: {readiness}"
                        );
                    }
                    assert_eq!(readiness["direct_submission_available"], false);
                    assert_eq!(
                        snapshot(root),
                        before,
                        "{case}: local preflight must not write"
                    );
                }
                if case == "partial" {
                    assert_eq!(
                        report["neutral_service_draft"]["optimization_target"]["yield_goal"]["metric"],
                        "active_protein"
                    );
                    assert!(
                        report["neutral_service_draft"]["optimization_target"]
                            .get("chassis")
                            .is_none()
                    );
                    assert!(!report["missing_questions"].as_array().unwrap().is_empty());
                }
                if matches!(case, "missing_product" | "invalid_product") {
                    assert!(report["neutral_service_draft"].is_null());
                }
                if case == "protein" {
                    assert_eq!(
                        report["neutral_service_draft"]["source_target"]["coding_dna_status"],
                        "unresolved"
                    );
                    assert_eq!(
                        report["neutral_service_draft"]["source_target"]["reverse_translation_status"],
                        "not_performed"
                    );
                }
                if matches!(
                    case,
                    "quote_unavailable" | "blocked_preflight" | "partial" | "protein"
                ) {
                    assert!(
                        report["suggested_next_actions"]
                            .as_array()
                            .unwrap()
                            .iter()
                            .all(|action| !action["shell_line"]
                                .as_str()
                                .unwrap()
                                .starts_with("services project-quote "))
                    );
                }
                if let Some(text) = report["text_report"].as_str() {
                    for action in report["suggested_next_actions"].as_array().unwrap() {
                        let line = action["shell_line"].as_str().unwrap();
                        if line.starts_with("services project-") {
                            assert!(
                                text.contains(line),
                                "{case}: text omitted an available action"
                            );
                        }
                    }
                    assert_eq!(
                        text.contains("services project-quote "),
                        case == "selected",
                        "{case}: text must not invent quote readiness"
                    );
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
