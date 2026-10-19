//! Typed, read-only MCP exposure of the existing expression handoff command.
//!
//! Transport validation admits inline JSON objects, not file-path operands.
//! Product readiness, provider eligibility and all interpretation stay shared.

use super::*;

#[derive(Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
struct Arguments {
    state_path: Option<String>,
    seq_id: Option<String>,
    objective: Option<Value>,
    requirements: Option<Value>,
    provider: Option<String>,
    profile_scope: Option<String>,
    format: Option<String>,
}

pub(super) fn descriptor() -> Value {
    json!({
        "name": "protein_expression_handoff",
        "title": "Protein Expression Handoff",
        "description": "Inspect the shared read-only protein-expression handoff.",
        "inputSchema": {
            "type": "object",
            "properties": {
                "state_path": {"type": "string", "description": "Saved project to inspect; defaults to server startup state path. Unsaved GUI edits are not visible."},
                "seq_id": {"type": "string", "description": "Explicit product id; omission never selects a loaded product automatically."},
                "objective": {"type": "object", "description": "Inline gentle.planning_objective.v1; omit to use the saved objective."},
                "requirements": {"type": "object", "description": "Inline gentle.protein_expression_requirements.v1. Omitted topics stay unanswered; no example or provider default is loaded."},
                "provider": {"type": "string", "description": "Explicit optional catalog provider id. Selection is not outsourcing permission or a request to submit."},
                "profile_scope": {"type": "string", "enum": ["effective"]},
                "format": {"type": "string", "enum": ["json", "text"], "description": "Default json; text adds the shared text_report inside the same structured report."}
            },
            "additionalProperties": false
        }
    })
}

fn command(args: &Arguments) -> Result<Vec<String>, String> {
    let mut tokens = vec!["planning".into(), "protein-expression-handoff".into()];
    for (flag, value) in [
        ("--seq-id", &args.seq_id),
        ("--provider", &args.provider),
        ("--profile-scope", &args.profile_scope),
        ("--format", &args.format),
    ] {
        if let Some(value) = value {
            tokens.extend([flag.into(), value.clone()]);
        }
    }
    for (flag, value) in [
        ("--objective", &args.objective),
        ("--requirements", &args.requirements),
    ] {
        if let Some(value) = value {
            if !value.is_object() {
                return Err(format!(
                    "{flag} must be an inline JSON object, not a path or string"
                ));
            }
            tokens.extend([flag.into(), value.to_string()]);
        }
    }
    Ok(tokens)
}

pub(super) fn call(default_state_path: &str, arguments: &Value) -> Value {
    let result = (|| {
        if let Some((field, _)) = arguments
            .as_object()
            .and_then(|args| args.iter().find(|(_, value)| value.is_null()))
        {
            return Err(format!(
                "MCP argument '{field}' must not be null; omit unanswered inputs instead"
            ));
        }
        let args: Arguments = serde_json::from_value(if arguments.is_null() {
            json!({})
        } else {
            arguments.clone()
        })
        .map_err(|err| format!("Invalid protein_expression_handoff arguments: {err}"))?;
        let mut state_args = Map::new();
        if let Some(path) = &args.state_path {
            state_args.insert("state_path".into(), json!(path));
        }
        run_non_mutating_shell_tool(
            default_state_path,
            &state_args,
            command(&args)?,
            "protein_expression_handoff",
        )
    })();
    match result {
        Ok(report) => tool_result_json(report, false),
        Err(error) => tool_result_text(error, "text", true),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine_shell::{execute_shell_command, quote_shell_arg};

    fn canonical(mut report: Value) -> Value {
        report
            .as_object_mut()
            .unwrap()
            .remove("generated_at_unix_ms");
        report
    }

    #[test]
    fn protein_expression_handoff_mcp_discovery_is_read_only_and_registered() {
        let tools = tool_list();
        let tool = tools
            .as_array()
            .unwrap()
            .iter()
            .find(|tool| tool["name"] == "protein_expression_handoff")
            .unwrap();
        assert_eq!(tool["mutating"], false);
        assert_eq!(
            tool["commandPaths"],
            json!(["planning protein-expression-handoff"])
        );
        assert_eq!(tool["inputSchema"]["additionalProperties"], false);
        assert!(tool["inputSchema"].get("required").is_none());
        let descriptor = gentle_protocol::capability_descriptor(
            CapabilitySource::McpTool,
            "protein_expression_handoff",
        )
        .unwrap();
        assert_eq!(
            descriptor.mutating,
            gentle_protocol::CapabilityMutation::ReadOnly
        );
    }

    #[test]
    fn protein_expression_handoff_mcp_and_inner_agent_match_shared_report_without_writes() {
        // Hand-crafted synthetic CDS and reviewed fixtures from plans/README.md,
        // never a natural product, provider call or live model invocation.
        let temp = tempfile::tempdir().unwrap();
        let state_path = temp.path().join("expression review.json");
        let mut state = ProjectState::default();
        state.sequences.insert(
            "expression_review_product".into(),
            crate::dna_sequence::DNAsequence::from_sequence("ATGCACCACCACCACCACCACGCTTAA").unwrap(),
        );
        state.save_to_path(&state_path.to_string_lossy()).unwrap();
        let before = std::fs::read(&state_path).unwrap();
        for (product, requirements, provider, expected_status) in [
            (None, None, None, "not_selected"),
            (
                Some("expression_review_product"),
                None,
                None,
                "not_selected",
            ),
            (
                Some("expression_review_product"),
                Some(include_str!(
                    "../../docs/examples/plans/protein_expression_requirements_in_house.json"
                )),
                None,
                "not_selected",
            ),
            (
                Some("expression_review_product"),
                Some(include_str!(
                    "../../docs/examples/plans/protein_expression_requirements_in_house.json"
                )),
                Some("geneart"),
                "withheld",
            ),
            (
                Some("expression_review_product"),
                Some(include_str!(
                    "../../docs/examples/plans/protein_expression_requirements_review_only.json"
                )),
                Some("geneart"),
                "projected",
            ),
        ] {
            for format in ["json", "text"] {
                let mut arguments = json!({"state_path": state_path, "format": format});
                if let Some(product) = product {
                    arguments["seq_id"] = json!(product);
                }
                if let Some(requirements) = requirements {
                    arguments["requirements"] = serde_json::from_str(requirements).unwrap();
                }
                if let Some(provider) = provider {
                    arguments["provider"] = json!(provider);
                }
                let args: Arguments = serde_json::from_value(arguments.clone()).unwrap();
                let tokens = command(&args).unwrap();
                let mut engine = GentleEngine::from_state(state.clone());
                let expected =
                    execute_shell_command(&mut engine, &parse_shell_tokens(&tokens).unwrap())
                        .unwrap();
                assert!(!expected.state_changed);
                assert_eq!(
                    expected.output["provider_binding"]["status"],
                    expected_status
                );
                let mcp = mcp_tool_call_for_capability_surface_tests(
                    DEFAULT_MCP_STATE_PATH,
                    "protein_expression_handoff",
                    arguments.clone(),
                );
                assert_eq!(mcp["isError"], false, "{mcp}");
                assert_eq!(
                    canonical(mcp["structuredContent"].clone()),
                    canonical(expected.output.clone())
                );
                // Same preparation path used by the inner assistant, with no
                // model and no provider-dependent prose treated as an oracle.
                let mut agent_engine = GentleEngine::from_state(state.clone());
                let line = tokens
                    .iter()
                    .map(|token| quote_shell_arg(token))
                    .collect::<Vec<_>>()
                    .join(" ");
                let response: crate::agent_bridge::AgentResponse = serde_json::from_value(json!({
                    "schema": "gentle.agent_response.v1",
                    "assistant_message": "Synthetic local review, not expression evidence",
                    "suggested_commands": [{"command": line, "execution": "ask"}]
                }))
                .unwrap();
                let parsed =
                    crate::engine_shell::parse_shell_line(&response.suggested_commands[0].command)
                        .unwrap();
                let agent = execute_shell_command_with_options(
                    &mut agent_engine,
                    &parsed,
                    &ShellExecutionOptions {
                        allow_agent_commands: false,
                        ..ShellExecutionOptions::default()
                    },
                )
                .unwrap();
                assert!(!agent.state_changed);
                assert_eq!(canonical(agent.output), canonical(expected.output));
                assert_eq!(std::fs::read(&state_path).unwrap(), before);
                assert_eq!(std::fs::read_dir(temp.path()).unwrap().count(), 1);
            }
        }
    }

    #[test]
    fn protein_expression_handoff_mcp_refuses_unknown_fields_paths_and_shared_invalid_inputs() {
        for arguments in [
            json!({"requirements": "@private.json"}),
            json!({"requirements": null}),
            json!({"format": null}),
            json!({"objective": "@private.json"}),
            json!({"provider": true}),
            json!({"confirm": true}),
            json!({"format": "yaml"}),
            json!({"profile_scope": "global"}),
            json!({"requirements": {"schema": "gentle.protein_expression_requirements.v1", "invented": true}}),
        ] {
            let temp = tempfile::tempdir().unwrap();
            let state_path = temp.path().join("review.json");
            ProjectState::default()
                .save_to_path(&state_path.to_string_lossy())
                .unwrap();
            let before = std::fs::read(&state_path).unwrap();
            let result = mcp_tool_call_for_capability_surface_tests(
                &state_path.to_string_lossy(),
                "protein_expression_handoff",
                arguments.clone(),
            );
            assert_eq!(result["isError"], true, "{arguments}: {result}");
            assert_eq!(std::fs::read(&state_path).unwrap(), before);
            assert_eq!(std::fs::read_dir(temp.path()).unwrap().count(), 1);
        }
    }
}
