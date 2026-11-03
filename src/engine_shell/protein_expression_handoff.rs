//! Inline-only admission for named, read-only expression-handoff adapters.
//!
//! CLI file operands retain their existing grammar. These transport bindings
//! cannot turn an objective/requirements string into implicit file ingress.

#[cfg(any(test, feature = "js-interface", feature = "lua-interface"))]
use super::{
    GentleEngine, ProjectState, ShellCommand, ShellRunResult, execute_shell_command,
    parse_shell_tokens,
};
#[cfg(any(test, feature = "js-interface", feature = "lua-interface"))]
use gentle_protocol::{EngineError, ErrorCode};
use serde::Deserialize;
use serde_json::Value;

#[derive(Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
struct Arguments {
    seq_id: Option<String>,
    objective: Option<Value>,
    requirements: Option<Value>,
    provider: Option<String>,
    profile_scope: Option<String>,
    format: Option<String>,
}

/// Check before removing transport-specific fields or normalizing Lua values.
pub(crate) fn reject_null_fields(options: &Value, label: &str) -> Result<(), String> {
    if let Some((field, _)) = options
        .as_object()
        .and_then(|args| args.iter().find(|(_, value)| value.is_null()))
    {
        return Err(format!(
            "{label} '{field}' must not be null; omit unanswered inputs instead"
        ));
    }
    Ok(())
}

/// Preserve omission and empty strings; the existing parser owns their meaning.
/// MCP keeps its existing state-load/parser error boundary around these tokens.
pub(crate) fn inline_tokens(options: &Value) -> Result<Vec<String>, String> {
    reject_null_fields(options, "Protein-expression handoff argument")?;
    if !options.is_object() {
        return Err(
            "Invalid protein_expression_handoff arguments: expected an options object".into(),
        );
    }
    let args: Arguments = serde_json::from_value(options.clone())
        .map_err(|err| format!("Invalid protein_expression_handoff arguments: {err}"))?;
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

#[cfg(any(test, feature = "js-interface", feature = "lua-interface"))]
pub(crate) fn prepare(options: &Value) -> Result<ShellCommand, String> {
    let tokens = inline_tokens(options)?;
    let command = parse_shell_tokens(&tokens)?;
    if !matches!(
        command,
        ShellCommand::PlanningProteinExpressionHandoff { .. }
    ) {
        return Err("Protein-expression handoff parsed to an unexpected command".into());
    }
    Ok(command)
}

#[cfg(any(test, feature = "js-interface", feature = "lua-interface"))]
fn read_only_output(run: ShellRunResult) -> Result<Value, EngineError> {
    if run.state_changed {
        return Err(EngineError::new(
            ErrorCode::Internal,
            "Protein-expression handoff unexpectedly changed project state",
        ));
    }
    Ok(run.output)
}

/// Inspect a private in-memory snapshot; never load/save the caller's project.
#[cfg(any(test, feature = "js-interface", feature = "lua-interface"))]
pub(crate) fn inspect(state: ProjectState, options: &Value) -> Result<Value, EngineError> {
    let command = prepare(options).map_err(EngineError::invalid_input)?;
    let mut engine = GentleEngine::from_state(state);
    let run = execute_shell_command(&mut engine, &command).map_err(EngineError::invalid_input)?;
    read_only_output(run)
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    #[test]
    fn protein_expression_handoff_inline_admission_preserves_operands() {
        // Synthetic identifiers exercise transport spelling, not real products.
        let options = json!({"seq_id":"", "objective":{}, "format":"text"});
        let command = prepare(&options).unwrap();
        let tokens = inline_tokens(&options).unwrap();
        assert!(matches!(
            command,
            ShellCommand::PlanningProteinExpressionHandoff { .. }
        ));
        assert!(tokens.windows(2).any(|pair| pair == ["--seq-id", ""]));
        assert!(tokens.windows(2).any(|pair| pair == ["--objective", "{}"]));
        assert!(
            inline_tokens(&json!({}))
                .unwrap()
                .iter()
                .all(|token| token != "--objective")
        );
        let report = inspect(ProjectState::default(), &json!({})).unwrap();
        assert_eq!(report["schema"], "gentle.protein_expression_handoff.v2");
        assert_eq!(report["status"], "needs_product_definition");
        assert_eq!(
            prepare(&json!({"provider":""})).err().unwrap(),
            parse_shell_tokens(&[
                "planning".into(),
                "protein-expression-handoff".into(),
                "--provider".into(),
                "".into()
            ])
            .unwrap_err()
        );
    }

    #[test]
    fn protein_expression_handoff_inline_admission_refuses_paths_nulls_and_unknown_fields() {
        for options in [
            json!(null),
            json!([]),
            json!("@private.json"),
            json!({"objective":null}),
            json!({"requirements":null}),
            json!({"objective":"@private.json"}),
            json!({"requirements":"C:\\private.json"}),
            json!({"requirements":"existing.json"}),
            json!({"objective":"a/b"}),
            json!({"confirm":true}),
            json!({"state_path":"other.json"}),
            json!({"format":"yaml"}),
            json!({"profile_scope":"global"}),
            json!({"provider":false}),
        ] {
            assert!(prepare(&options).is_err(), "{options}");
        }
        let temp = tempfile::tempdir().unwrap();
        let file = temp.path().join("existing review's \u{3bc}.json");
        std::fs::write(
            &file,
            include_str!("../../docs/examples/plans/protein_expression_requirements_in_house.json"),
        )
        .unwrap();
        let before = std::fs::read(&file).unwrap();
        assert!(prepare(&json!({"requirements":file.to_str().unwrap()})).is_err());
        assert_eq!(std::fs::read(file).unwrap(), before);
    }

    #[test]
    fn protein_expression_handoff_report_only_guard_rejects_state_changes() {
        let error = read_only_output(ShellRunResult {
            state_changed: true,
            output: json!({}),
        })
        .unwrap_err();
        assert_eq!(error.code, ErrorCode::Internal);
        assert_eq!(
            read_only_output(ShellRunResult {
                state_changed: false,
                output: json!({"status":"needs_product_definition"})
            })
            .unwrap()["status"],
            "needs_product_definition"
        );
    }
}
