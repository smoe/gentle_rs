//! Deterministic execution for stored agent plans.
//!
//! This module consumes previously generated `gentle.agent_plan_result.v1`
//! payloads and executes one selected candidate without re-planning. That
//! keeps the execution step auditable, replayable, and adapter-equivalent
//! across CLI, MCP, and external orchestrators.

use crate::{
    agent_planner::{AgentPlanCandidate, AgentPlanCandidateKind, AgentPlanResult},
    engine::GentleEngine,
    engine_shell::{
        ShellCommand, ShellExecutionOptions, ShellRunResult, execute_shell_command_with_options,
        parse_shell_line,
    },
};
use serde::{Deserialize, Serialize};
use serde_json::Value;

/// Stored execution result for one selected plan candidate.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct AgentExecutionResult {
    pub schema: String,
    pub candidate_id: String,
    pub title: String,
    pub kind: AgentPlanCandidateKind,
    pub mutating: bool,
    pub requires_confirmation: bool,
    pub confirmed: bool,
    pub state_changed: bool,
    pub output: Value,
}

fn ui_intent_to_shell_command(value: &Value) -> Result<String, String> {
    let obj = value
        .as_object()
        .ok_or_else(|| "ui_intent payload must be a JSON object".to_string())?;
    let action = obj
        .get("action")
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .ok_or_else(|| "ui_intent.action must be a non-empty string".to_string())?;
    let target = obj
        .get("target")
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .ok_or_else(|| "ui_intent.target must be a non-empty string".to_string())?;
    let mut tokens = vec!["ui".to_string(), action.to_string(), target.to_string()];
    if let Some(genome_id) = obj
        .get("genome_id")
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|value| !value.is_empty())
    {
        tokens.push("--genome-id".to_string());
        tokens.push(genome_id.to_string());
    }
    if obj.get("helpers").and_then(Value::as_bool).unwrap_or(false) {
        tokens.push("--helpers".to_string());
    }
    for (flag, key) in [
        ("--catalog", "catalog_path"),
        ("--cache-dir", "cache_dir"),
        ("--filter", "filter"),
        ("--species", "species"),
    ] {
        if let Some(value) = obj
            .get(key)
            .and_then(Value::as_str)
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            tokens.push(flag.to_string());
            tokens.push(value.to_string());
        }
    }
    if obj.get("latest").and_then(Value::as_bool).unwrap_or(false) {
        tokens.push("--latest".to_string());
    }
    Ok(tokens.join(" "))
}

fn shell_command_for_candidate(candidate: &AgentPlanCandidate) -> Result<Option<String>, String> {
    match candidate.kind {
        AgentPlanCandidateKind::Shell => Ok(candidate
            .shell_command
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(ToString::to_string)),
        AgentPlanCandidateKind::UiIntent => candidate
            .ui_intent
            .as_ref()
            .map(ui_intent_to_shell_command)
            .transpose(),
        AgentPlanCandidateKind::Op | AgentPlanCandidateKind::Workflow => Ok(None),
    }
}

fn command_is_blocked(parsed: &ShellCommand) -> bool {
    matches!(
        parsed,
        ShellCommand::AgentsAsk { .. }
            | ShellCommand::AgentsPlan { .. }
            | ShellCommand::AgentsExecutePlan { .. }
    )
}

fn execute_shell_like_candidate(
    engine: &mut GentleEngine,
    candidate: &AgentPlanCandidate,
    command_text: &str,
    options: &ShellExecutionOptions,
) -> Result<ShellRunResult, String> {
    let parsed = parse_shell_line(command_text).map_err(|err| {
        format!(
            "Could not parse shell command for planner candidate '{}': {err}",
            candidate.candidate_id
        )
    })?;
    if command_is_blocked(&parsed) {
        return Err(
            "Agent plan execution blocks nested agent assistant/planner commands".to_string(),
        );
    }
    let nested_options = ShellExecutionOptions {
        allow_screenshots: options.allow_screenshots,
        allow_agent_commands: false,
        progress_callback: options.progress_callback.clone(),
    };
    execute_shell_command_with_options(engine, &parsed, &nested_options)
}

fn execute_structured_candidate(
    engine: &mut GentleEngine,
    candidate: &AgentPlanCandidate,
    confirm: bool,
    options: &ShellExecutionOptions,
) -> Result<ShellRunResult, String> {
    match candidate.kind {
        AgentPlanCandidateKind::Shell | AgentPlanCandidateKind::UiIntent => {
            let command_text = shell_command_for_candidate(candidate)?
                .ok_or_else(|| "planner candidate is missing shell command".to_string())?;
            execute_shell_like_candidate(engine, candidate, &command_text, options)
        }
        AgentPlanCandidateKind::Op => {
            if (candidate.mutating || candidate.requires_confirmation) && !confirm {
                return Err(format!(
                    "Planner candidate '{}' requires explicit --confirm before operation execution",
                    candidate.candidate_id
                ));
            }
            let payload = serde_json::to_string(
                candidate
                    .operation
                    .as_ref()
                    .ok_or_else(|| "planner candidate missing operation payload".to_string())?,
            )
            .map_err(|err| format!("Could not encode operation payload: {err}"))?;
            execute_shell_command_with_options(
                engine,
                &ShellCommand::Op { payload },
                &ShellExecutionOptions {
                    allow_screenshots: options.allow_screenshots,
                    allow_agent_commands: false,
                    progress_callback: options.progress_callback.clone(),
                },
            )
        }
        AgentPlanCandidateKind::Workflow => {
            if (candidate.mutating || candidate.requires_confirmation) && !confirm {
                return Err(format!(
                    "Planner candidate '{}' requires explicit --confirm before workflow execution",
                    candidate.candidate_id
                ));
            }
            let payload = serde_json::to_string(
                candidate
                    .workflow
                    .as_ref()
                    .ok_or_else(|| "planner candidate missing workflow payload".to_string())?,
            )
            .map_err(|err| format!("Could not encode workflow payload: {err}"))?;
            execute_shell_command_with_options(
                engine,
                &ShellCommand::Workflow { payload },
                &ShellExecutionOptions {
                    allow_screenshots: options.allow_screenshots,
                    allow_agent_commands: false,
                    progress_callback: options.progress_callback.clone(),
                },
            )
        }
    }
}

/// Execute one selected candidate from a stored plan without re-planning.
pub fn execute_agent_plan_candidate(
    engine: &mut GentleEngine,
    plan: &AgentPlanResult,
    candidate_id: &str,
    confirm: bool,
    options: &ShellExecutionOptions,
) -> Result<AgentExecutionResult, String> {
    plan.validate()?;
    let selected = plan
        .candidates
        .iter()
        .find(|candidate| candidate.candidate_id == candidate_id)
        .ok_or_else(|| format!("planner candidate '{}' was not found", candidate_id))?;
    selected.validate()?;
    let run = execute_structured_candidate(engine, selected, confirm, options)?;
    Ok(AgentExecutionResult {
        schema: "gentle.agent_execution_result.v1".to_string(),
        candidate_id: selected.candidate_id.clone(),
        title: selected.title.clone(),
        kind: selected.kind,
        mutating: selected.mutating,
        requires_confirmation: selected.requires_confirmation,
        confirmed: confirm,
        state_changed: run.state_changed,
        output: run.output,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        agent_planner::{AgentPlanExecutionMode, AgentPlanResult},
        engine::ProjectState,
    };

    #[test]
    fn ui_intent_candidate_is_translated_to_ui_shell_command() {
        let value = serde_json::json!({
            "action": "open",
            "target": "agent-assistant",
            "helpers": true,
            "species": "human"
        });
        let shell = ui_intent_to_shell_command(&value).expect("ui intent to shell");
        assert_eq!(
            shell,
            "ui open agent-assistant --helpers --species human".to_string()
        );
    }

    #[test]
    fn operation_candidate_requires_confirm() {
        let plan = AgentPlanResult {
            schema: "gentle.agent_plan_result.v1".to_string(),
            assistant_message: "ready".to_string(),
            questions: vec![],
            candidates: vec![AgentPlanCandidate {
                candidate_id: "candidate-1".to_string(),
                title: "set".to_string(),
                rationale: String::new(),
                kind: AgentPlanCandidateKind::Op,
                mutating: true,
                requires_confirmation: true,
                execution_mode: AgentPlanExecutionMode::Ask,
                shell_command: None,
                operation: Some(serde_json::json!({
                    "SetParameter": {"name":"max_fragments_per_container","value": 42}
                })),
                workflow: None,
                ui_intent: None,
            }],
        };
        let mut engine = GentleEngine::from_state(ProjectState::default());
        let err = execute_agent_plan_candidate(
            &mut engine,
            &plan,
            "candidate-1",
            false,
            &ShellExecutionOptions::default(),
        )
        .expect_err("confirm required");
        assert!(err.contains("--confirm"));
    }
}
