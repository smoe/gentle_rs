//! Bounded MCP adapters for stored primer reports and panel-specificity handoffs.
//!
//! All execution uses fixed shared shell routes. Confirmation precedes state or
//! file access for exports, handoff preparation, and specificity finalization.

use super::*;

const REPORT_PATHS: &[&str] = &[
    "primers list-reports",
    "primers show-report",
    "primers export-report",
    "primers list-qpcr-reports",
    "primers show-qpcr-report",
    "primers export-qpcr-report",
    "primers list-transcript-assay-panels",
    "primers show-transcript-assay-panel",
    "primers export-transcript-assay-panel",
    "primers list-transcript-assay-fallbacks",
    "primers show-transcript-assay-fallback",
    "primers export-transcript-assay-fallback",
];

// Adapter field names and scalar types only; scientific policy validation stays
// in the shared parser/engine, including the defaults for omitted fields.
const PLAN_FLAGS: &[(&str, &str, &str)] = &[
    (
        "max_target_amplicon_bp",
        "--max-target-amplicon-bp",
        "integer",
    ),
    (
        "readiness_max_amplicon_bp",
        "--readiness-max-amplicon-bp",
        "integer",
    ),
    (
        "exploratory_max_amplicon_bp",
        "--exploratory-max-amplicon-bp",
        "integer",
    ),
    ("report_detail", "--report-detail", "string"),
    ("full_alignment", "--full-alignment", "string"),
    (
        "min_primer_coverage_fraction",
        "--min-primer-coverage-fraction",
        "number",
    ),
    (
        "max_3prime_mismatches",
        "--max-3prime-mismatches",
        "integer",
    ),
    (
        "three_prime_window_bp",
        "--three-prime-window-bp",
        "integer",
    ),
    (
        "min_total_mismatches_to_unintended_target",
        "--min-total-mismatches-to-unintended-target",
        "integer",
    ),
    ("max_hits_per_primer", "--max-hits-per-primer", "integer"),
    (
        "allow_same_gene_splice_variants",
        "--allow-same-gene-splice-variants",
        "boolean",
    ),
    ("avoid_known_variants", "--avoid-known-variants", "boolean"),
    ("avoid_rmsk_repeats", "--avoid-rmsk-repeats", "boolean"),
    ("avoid_low_complexity", "--avoid-low-complexity", "boolean"),
    ("catalog_path", "--catalog", "string"),
    ("cache_dir", "--cache-dir", "string"),
];

pub(super) fn command_paths(name: &str) -> Option<&'static [&'static str]> {
    match name {
        "primer_reports" => Some(REPORT_PATHS),
        "transcript_assay_specificity_plan" => Some(&["primers transcript-assay-specificity-plan"]),
        "transcript_assay_specificity_finalize" => {
            Some(&["primers transcript-assay-specificity-finalize"])
        }
        _ => None,
    }
}

pub(super) fn descriptors() -> Vec<Value> {
    let mut plan_properties = json!({
        "panel_report_id": {"type":"string"},
        "target_genome_id": {"type":"string"},
        "output_dir": {"type":"string"},
        "state_path": {"type":"string"},
        "confirm": {"type":"boolean", "const":true},
        "reviewed_off_target_allowlist": {"type":"array", "items":{"type":"object"}}
    });
    for (field, _, kind) in PLAN_FLAGS {
        plan_properties[*field] = json!({"type":kind});
    }
    plan_properties["report_detail"]["enum"] = json!(["compact", "full"]);
    plan_properties["full_alignment"]["enum"] = json!(["disabled", "best-effort", "required"]);
    vec![
        json!({
            "name":"primer_reports",
            "inputSchema": {
                "type":"object", "additionalProperties":false,
                "required":["family","action"],
                "properties": {
                    "family":{"type":"string","enum":["primer","qpcr","transcript_assay_panel","transcript_assay_fallback"]},
                    "action":{"type":"string","enum":["list","show","export"]},
                    "report_id":{"type":"string","description":"Required for show/export; use the execution_id for a fallback."},
                    "path":{"type":"string","description":"Required only for export."},
                    "confirm":{"type":"boolean","description":"Export requires explicit true; list/show need no approval."},
                    "state_path":{"type":"string"}
                }
            }
        }),
        json!({
            "name":"transcript_assay_specificity_plan",
            "inputSchema":{"type":"object","additionalProperties":false,
                "required":["panel_report_id","target_genome_id","output_dir","confirm"],
                "properties":plan_properties}
        }),
        json!({
            "name":"transcript_assay_specificity_finalize",
            "inputSchema":{"type":"object","additionalProperties":false,
                "required":["handoff_path","execution_manifest","confirm"],
                "properties":{
                    "handoff_path":{"type":"string"},
                    "execution_manifest":{"type":"object","description":"gentle.transcript_assay_panel_specificity_execution_manifest.v1 process evidence, not a biological decision."},
                    "path":{"type":"string"},
                    "state_path":{"type":"string"},
                    "confirm":{"type":"boolean","const":true}
                }}
        }),
    ]
}

fn validate_args<'a>(name: &str, arguments: &'a Value) -> Result<&'a Map<String, Value>, String> {
    let args = arguments
        .as_object()
        .ok_or("Expected an arguments object")?;
    let descriptor = descriptors()
        .into_iter()
        .find(|item| item["name"] == name)
        .ok_or("Unknown primer tool")?;
    let schema = &descriptor["inputSchema"];
    for required in schema["required"].as_array().unwrap() {
        let field = required.as_str().unwrap();
        if !args.contains_key(field) {
            return Err(format!("Missing required argument '{field}'"));
        }
    }
    for (field, value) in args {
        let property = schema["properties"]
            .get(field)
            .ok_or_else(|| format!("Unknown argument '{field}' for '{name}'"))?;
        let valid = match property["type"].as_str() {
            Some("string") => value.as_str().is_some_and(|s| !s.trim().is_empty()),
            Some("integer") => value.as_u64().is_some(),
            Some("number") => value.as_f64().is_some_and(f64::is_finite),
            Some("boolean") => value.is_boolean(),
            Some("array") => value
                .as_array()
                .is_some_and(|a| a.iter().all(Value::is_object)),
            Some("object") => value.is_object(),
            _ => false,
        };
        if !valid
            || property
                .get("enum")
                .is_some_and(|choices| !choices.as_array().unwrap().contains(value))
        {
            return Err(format!("Invalid value for '{field}'"));
        }
    }
    Ok(args)
}

fn tokens(name: &str, args: &Map<String, Value>) -> Result<Vec<String>, String> {
    let string = |field| required_string_arg(args, field);
    if name == "primer_reports" {
        let family = match string("family")?.as_str() {
            "primer" => 0,
            "qpcr" => 1,
            "transcript_assay_panel" => 2,
            "transcript_assay_fallback" => 3,
            _ => return Err("Unknown report family".into()),
        };
        let action = string("action")?;
        let offset = match action.as_str() {
            "list" => 0,
            "show" => 1,
            "export" => 2,
            _ => return Err("Unknown report action".into()),
        };
        if (offset == 0 && args.contains_key("report_id"))
            || (offset != 2 && args.contains_key("path"))
        {
            return Err("report_id is only for show/export; path is only for export".into());
        }
        let mut tokens: Vec<String> = REPORT_PATHS[family * 3 + offset]
            .split_whitespace()
            .map(str::to_owned)
            .collect();
        if offset > 0 {
            tokens.push(string("report_id")?);
        }
        if offset == 2 {
            tokens.push(string("path")?);
        }
        return Ok(tokens);
    }
    if name == "transcript_assay_specificity_finalize" {
        let mut tokens = vec![
            "primers".into(),
            "transcript-assay-specificity-finalize".into(),
            string("handoff_path")?,
            args["execution_manifest"].to_string(),
        ];
        append_string_flag(&mut tokens, "--path", optional_string_arg(args, "path")?);
        return Ok(tokens);
    }
    if args.contains_key("max_target_amplicon_bp")
        && (args.contains_key("readiness_max_amplicon_bp")
            || args.contains_key("exploratory_max_amplicon_bp"))
    {
        return Err("Use either max_target_amplicon_bp or separate readiness/exploratory ceilings, not both".into());
    }
    let mut tokens = vec![
        "primers".into(),
        "transcript-assay-specificity-plan".into(),
        string("panel_report_id")?,
        "--target-genome".into(),
        string("target_genome_id")?,
        "--output-dir".into(),
        string("output_dir")?,
    ];
    for (field, flag, kind) in PLAN_FLAGS {
        if let Some(value) = args.get(*field) {
            if *kind == "boolean" {
                if value == true {
                    tokens.push((*flag).into());
                }
            } else {
                tokens.push((*flag).into());
                tokens.push(
                    value
                        .as_str()
                        .map(str::to_owned)
                        .unwrap_or_else(|| value.to_string()),
                );
            }
        }
    }
    if let Some(allowlist) = args.get("reviewed_off_target_allowlist") {
        tokens.extend([
            "--reviewed-off-target-allowlist".into(),
            allowlist.to_string(),
        ]);
    }
    Ok(tokens)
}

pub(super) fn call(default_state_path: &str, name: &str, arguments: &Value) -> Value {
    let run = || -> Result<Value, String> {
        let args = validate_args(name, arguments)?;
        let writes = name != "primer_reports" || args["action"] == "export";
        if writes {
            require_confirm_true(args, name)?;
        }
        let tokens = tokens(name, args)?;
        if name == "transcript_assay_specificity_finalize" {
            run_shell_tool_with_optional_persist(default_state_path, args, tokens, name)
        } else {
            run_non_mutating_shell_tool(default_state_path, args, tokens, name)
        }
    };
    match run() {
        Ok(output) => tool_result_json(output, false),
        Err(error) => tool_result_text(error, "text", true),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn report_routes_match_shared_shell_and_do_not_create_state() {
        let dir = tempfile::tempdir().unwrap();
        let state_path = dir.path().join("absent.json");
        for family in [
            "primer",
            "qpcr",
            "transcript_assay_panel",
            "transcript_assay_fallback",
        ] {
            for action in ["list", "show", "export"] {
                let mut args = json!({"family":family,"action":action});
                if action != "list" {
                    args["report_id"] = json!("missing-report");
                }
                if action == "export" {
                    args["path"] = json!(dir.path().join("missing-output.json"));
                    args["confirm"] = json!(true);
                }
                let shell_tokens = tokens("primer_reports", args.as_object().unwrap()).unwrap();
                let expected = run_non_mutating_shell_tool(
                    state_path.to_str().unwrap(),
                    args.as_object().unwrap(),
                    shell_tokens,
                    "primer_reports",
                );
                let result = call(state_path.to_str().unwrap(), "primer_reports", &args);
                match expected {
                    Ok(payload) => assert_eq!(result, tool_result_json(payload, false)),
                    Err(error) => assert_eq!(result, tool_result_text(error, "text", true)),
                }
            }
        }
        assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 0);
    }

    #[test]
    fn writing_routes_require_confirmation_before_state_or_file_access() {
        let dir = tempfile::tempdir().unwrap();
        let state = dir.path().join("invalid.json");
        std::fs::write(&state, "not a project").unwrap();
        let mut cases = vec![
            (
                "transcript_assay_specificity_plan",
                json!({"panel_report_id":"panel","target_genome_id":"genome","output_dir":dir.path().join("output"),"confirm":false}),
            ),
            (
                "transcript_assay_specificity_finalize",
                json!({"handoff_path":"missing.json","execution_manifest":{},"confirm":false}),
            ),
        ];
        for family in [
            "primer",
            "qpcr",
            "transcript_assay_panel",
            "transcript_assay_fallback",
        ] {
            cases.push(("primer_reports", json!({"family":family,"action":"export","report_id":"report","path":dir.path().join("export.json"),"confirm":false})));
        }
        for (name, mut args) in cases {
            for confirm in [Some(json!(false)), None] {
                if let Some(value) = confirm {
                    args["confirm"] = value;
                } else {
                    args.as_object_mut().unwrap().remove("confirm");
                }
                let result = call(state.to_str().unwrap(), name, &args);
                assert_eq!(result["isError"], true, "{name}: {result}");
                assert!(result.to_string().contains("confirm"), "{result}");
                assert!(!result.to_string().contains("Could not load state"));
            }
        }
        assert_eq!(std::fs::read_to_string(state).unwrap(), "not a project");
        assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 1);
    }

    #[test]
    fn invalid_arguments_cannot_inject_routes_or_silently_ignore_fields() {
        for args in [
            json!({"family":"primer","action":"agents ask"}),
            json!({"family":"primer","action":"list","path":"ignored"}),
            json!({"family":"primer","action":"list","command":"agents ask"}),
            json!({"family":"primer","action":"show"}),
            json!({"family":"primer","action":"list","report_id":"ignored"}),
            json!({"family":"primer","action":"list","state_path":3}),
            json!([]),
        ] {
            assert_eq!(call("unused", "primer_reports", &args)["isError"], true);
        }
    }

    #[test]
    fn plan_fields_use_shared_parser_policy_and_finalize_preserves_manifest() {
        let args = json!({"panel_report_id":"panel","target_genome_id":"human","output_dir":"out with spaces","confirm":true,
            "readiness_max_amplicon_bp":800,"exploratory_max_amplicon_bp":2000,
            "min_primer_coverage_fraction":0.9,"max_3prime_mismatches":1,
            "three_prime_window_bp":5,"min_total_mismatches_to_unintended_target":3,
            "max_hits_per_primer":123,"report_detail":"compact","full_alignment":"required",
            "allow_same_gene_splice_variants":true,"avoid_known_variants":true,
            "avoid_rmsk_repeats":true,"avoid_low_complexity":false,
            "catalog_path":"catalog.json","cache_dir":"cache","reviewed_off_target_allowlist":[]});
        let args = validate_args("transcript_assay_specificity_plan", &args).unwrap();
        let command =
            parse_shell_tokens(&tokens("transcript_assay_specificity_plan", args).unwrap())
                .unwrap();
        let crate::engine_shell::ShellCommand::PrimersTranscriptAssaySpecificityPlan {
            panel_report_id,
            target_genome_id,
            output_dir,
            policy,
            catalog_path,
            cache_dir,
        } = command
        else {
            panic!("Wrong shared route");
        };
        assert_eq!(
            (
                panel_report_id.as_str(),
                target_genome_id.as_str(),
                output_dir.as_str()
            ),
            ("panel", "human", "out with spaces")
        );
        assert_eq!(policy.readiness_max_target_amplicon_bp, Some(800));
        assert_eq!(policy.max_target_amplicon_bp, 2000);
        assert_eq!(policy.min_primer_coverage_fraction, 0.9);
        assert_eq!(policy.max_3prime_mismatches, 1);
        assert_eq!(policy.three_prime_window_bp, 5);
        assert_eq!(policy.min_total_mismatches_to_unintended_target, 3);
        assert_eq!(policy.max_hits_per_primer, 123);
        assert!(
            policy.allow_same_gene_splice_variants
                && policy.avoid_known_variants
                && policy.avoid_rmsk_repeats
        );
        assert!(!policy.avoid_low_complexity);
        assert_eq!(catalog_path.as_deref(), Some("catalog.json"));
        assert_eq!(cache_dir.as_deref(), Some("cache"));
        let args = json!({"handoff_path":"handoff.json","execution_manifest":{"schema":"gentle.transcript_assay_panel_specificity_execution_manifest.v1","commands":[]},"path":"out.json","confirm":true});
        let args = validate_args("transcript_assay_specificity_finalize", &args).unwrap();
        let command =
            parse_shell_tokens(&tokens("transcript_assay_specificity_finalize", args).unwrap())
                .unwrap();
        let crate::engine_shell::ShellCommand::PrimersTranscriptAssaySpecificityFinalize {
            handoff_path,
            execution_manifest_json,
            path,
        } = command
        else {
            panic!("Wrong finalize route");
        };
        assert_eq!(handoff_path, "handoff.json");
        assert_eq!(
            serde_json::from_str::<Value>(&execution_manifest_json).unwrap(),
            args["execution_manifest"]
        );
        assert_eq!(path.as_deref(), Some("out.json"));
    }
}
