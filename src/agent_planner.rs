//! Machine-facing prose compiler for GENtle agent systems.
//!
//! This module accepts free prose plus optional state context and produces a
//! typed, auditable execution plan. The plan is intentionally transport-agnostic
//! at the boundary: it can be generated from the local assistant bridge today
//! and consumed later by CLI, MCP, or external orchestrators without re-planning.

use crate::{
    agent_bridge::{
        AGENT_BASE_URL_ENV, AGENT_CONNECT_TIMEOUT_SECS_ENV, AGENT_MAX_RESPONSE_BYTES_ENV,
        AGENT_MAX_RETRIES_ENV, AGENT_MODEL_ENV, AGENT_READ_TIMEOUT_SECS_ENV,
        AGENT_TIMEOUT_SECS_ENV, AgentExecutionIntent, AgentSystemTransport,
        invoke_agent_support_with_env_overrides, load_agent_system_catalog,
    },
    engine::EngineStateSummary,
    engine_shell::parse_shell_line,
};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::{collections::HashMap, fs};

const AGENT_PLAN_REQUEST_SCHEMA: &str = "gentle.agent_plan_request.v1";
const AGENT_PLAN_RESULT_SCHEMA: &str = "gentle.agent_plan_result.v1";

/// Planner request accepted by the machine-facing prose compiler.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct AgentPlanRequest {
    pub schema: String,
    pub system_id: String,
    pub prompt: String,
    pub state_summary: Option<EngineStateSummary>,
    pub max_candidates: Option<usize>,
    pub allow_mutating_candidates: Option<bool>,
}

/// Candidate execution kind returned by the planner.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum AgentPlanCandidateKind {
    #[default]
    Shell,
    Op,
    Workflow,
    UiIntent,
}

/// Candidate execution mode returned by the planner.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum AgentPlanExecutionMode {
    #[default]
    Ask,
    Auto,
}

impl AgentPlanExecutionMode {
    pub fn from_agent_execution(intent: AgentExecutionIntent) -> Self {
        match intent {
            AgentExecutionIntent::Auto => Self::Auto,
            AgentExecutionIntent::Chat | AgentExecutionIntent::Ask => Self::Ask,
        }
    }
}

/// One decision-complete planner candidate.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct AgentPlanCandidate {
    pub candidate_id: String,
    pub title: String,
    pub rationale: String,
    pub kind: AgentPlanCandidateKind,
    pub mutating: bool,
    pub requires_confirmation: bool,
    pub execution_mode: AgentPlanExecutionMode,
    pub shell_command: Option<String>,
    pub operation: Option<Value>,
    pub workflow: Option<Value>,
    pub ui_intent: Option<Value>,
}

impl AgentPlanCandidate {
    fn payload_count(&self) -> usize {
        [
            self.shell_command.is_some(),
            self.operation.is_some(),
            self.workflow.is_some(),
            self.ui_intent.is_some(),
        ]
        .into_iter()
        .filter(|value| *value)
        .count()
    }

    /// Validate the candidate envelope and exactly-one payload constraint.
    pub fn validate(&self) -> Result<(), String> {
        if self.candidate_id.trim().is_empty() {
            return Err("planner candidate_id cannot be empty".to_string());
        }
        if self.title.trim().is_empty() {
            return Err(format!(
                "planner candidate '{}' title cannot be empty",
                self.candidate_id
            ));
        }
        if self.payload_count() != 1 {
            return Err(format!(
                "planner candidate '{}' must set exactly one payload field",
                self.candidate_id
            ));
        }
        match self.kind {
            AgentPlanCandidateKind::Shell => {
                let Some(command) = self
                    .shell_command
                    .as_deref()
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                else {
                    return Err(format!(
                        "planner candidate '{}' kind=shell requires non-empty shell_command",
                        self.candidate_id
                    ));
                };
                parse_shell_line(command).map_err(|err| {
                    format!(
                        "planner candidate '{}' shell_command is invalid: {err}",
                        self.candidate_id
                    )
                })?;
            }
            AgentPlanCandidateKind::Op => {
                if !matches!(self.operation, Some(Value::Object(_))) {
                    return Err(format!(
                        "planner candidate '{}' kind=op requires object field 'operation'",
                        self.candidate_id
                    ));
                }
            }
            AgentPlanCandidateKind::Workflow => {
                if !matches!(self.workflow, Some(Value::Object(_))) {
                    return Err(format!(
                        "planner candidate '{}' kind=workflow requires object field 'workflow'",
                        self.candidate_id
                    ));
                }
            }
            AgentPlanCandidateKind::UiIntent => {
                if !matches!(self.ui_intent, Some(Value::Object(_))) {
                    return Err(format!(
                        "planner candidate '{}' kind=ui_intent requires object field 'ui_intent'",
                        self.candidate_id
                    ));
                }
            }
        }
        Ok(())
    }
}

/// Structured plan result returned by the prose compiler.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct AgentPlanResult {
    pub schema: String,
    pub assistant_message: String,
    pub questions: Vec<String>,
    pub candidates: Vec<AgentPlanCandidate>,
}

impl AgentPlanRequest {
    /// Validate request shape and normalize defaults.
    pub fn validate(&self) -> Result<(), String> {
        if self.schema.trim() != AGENT_PLAN_REQUEST_SCHEMA {
            return Err(format!(
                "unsupported planner request schema '{}' (expected '{}')",
                self.schema, AGENT_PLAN_REQUEST_SCHEMA
            ));
        }
        if self.system_id.trim().is_empty() {
            return Err("planner request system_id cannot be empty".to_string());
        }
        if self.prompt.trim().is_empty() {
            return Err("planner request prompt cannot be empty".to_string());
        }
        Ok(())
    }
}

impl AgentPlanResult {
    /// Validate the result envelope and all nested candidates.
    pub fn validate(&self) -> Result<(), String> {
        if self.schema.trim() != AGENT_PLAN_RESULT_SCHEMA {
            return Err(format!(
                "unsupported planner result schema '{}' (expected '{}')",
                self.schema, AGENT_PLAN_RESULT_SCHEMA
            ));
        }
        if self.assistant_message.trim().is_empty()
            && self.questions.is_empty()
            && self.candidates.is_empty()
        {
            return Err(
                "planner result must include assistant_message, questions, or candidates"
                    .to_string(),
            );
        }
        for candidate in &self.candidates {
            candidate.validate()?;
        }
        Ok(())
    }
}

fn candidate_id(index_1based: usize) -> String {
    format!("candidate-{index_1based}")
}

fn shell_candidate_from_command(
    index_1based: usize,
    title: Option<&str>,
    rationale: Option<&str>,
    command: &str,
    execution_mode: AgentPlanExecutionMode,
) -> Result<AgentPlanCandidate, String> {
    let parsed = parse_shell_line(command).map_err(|err| {
        format!("planner candidate shell command '{command}' could not be parsed: {err}")
    })?;
    Ok(AgentPlanCandidate {
        candidate_id: candidate_id(index_1based),
        title: title
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .unwrap_or("Shell command")
            .to_string(),
        rationale: rationale.unwrap_or_default().trim().to_string(),
        kind: AgentPlanCandidateKind::Shell,
        mutating: parsed.is_state_mutating(),
        requires_confirmation: false,
        execution_mode,
        shell_command: Some(command.trim().to_string()),
        operation: None,
        workflow: None,
        ui_intent: None,
    })
}

fn builtin_echo_plan_result(
    prompt: &str,
    max_candidates: Option<usize>,
    allow_mutating_candidates: bool,
) -> Result<AgentPlanResult, String> {
    let mut candidates = vec![];
    for line in prompt.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let index_1based = candidates.len() + 1;
        if let Some(command) = trimmed.strip_prefix("auto:") {
            let candidate = shell_candidate_from_command(
                index_1based,
                Some("Auto shell command"),
                Some("Extracted from 'auto:' line in builtin echo planner"),
                command.trim(),
                AgentPlanExecutionMode::Auto,
            )?;
            candidates.push(candidate);
            continue;
        }
        if let Some(command) = trimmed.strip_prefix("ask:") {
            let candidate = shell_candidate_from_command(
                index_1based,
                Some("Shell command"),
                Some("Extracted from 'ask:' line in builtin echo planner"),
                command.trim(),
                AgentPlanExecutionMode::Ask,
            )?;
            candidates.push(candidate);
            continue;
        }
        if let Some(payload) = trimmed.strip_prefix("op:") {
            let operation = serde_json::from_str::<Value>(payload.trim())
                .map_err(|err| format!("builtin planner op payload is not valid JSON: {err}"))?;
            candidates.push(AgentPlanCandidate {
                candidate_id: candidate_id(index_1based),
                title: "Operation payload".to_string(),
                rationale: "Extracted from 'op:' line in builtin echo planner".to_string(),
                kind: AgentPlanCandidateKind::Op,
                mutating: true,
                requires_confirmation: true,
                execution_mode: AgentPlanExecutionMode::Ask,
                shell_command: None,
                operation: Some(operation),
                workflow: None,
                ui_intent: None,
            });
            continue;
        }
        if let Some(payload) = trimmed.strip_prefix("workflow:") {
            let workflow = serde_json::from_str::<Value>(payload.trim()).map_err(|err| {
                format!("builtin planner workflow payload is not valid JSON: {err}")
            })?;
            candidates.push(AgentPlanCandidate {
                candidate_id: candidate_id(index_1based),
                title: "Workflow payload".to_string(),
                rationale: "Extracted from 'workflow:' line in builtin echo planner".to_string(),
                kind: AgentPlanCandidateKind::Workflow,
                mutating: true,
                requires_confirmation: true,
                execution_mode: AgentPlanExecutionMode::Ask,
                shell_command: None,
                operation: None,
                workflow: Some(workflow),
                ui_intent: None,
            });
            continue;
        }
        if let Some(payload) = trimmed.strip_prefix("ui:") {
            let ui_intent = serde_json::from_str::<Value>(payload.trim())
                .map_err(|err| format!("builtin planner ui payload is not valid JSON: {err}"))?;
            candidates.push(AgentPlanCandidate {
                candidate_id: candidate_id(index_1based),
                title: "UI intent".to_string(),
                rationale: "Extracted from 'ui:' line in builtin echo planner".to_string(),
                kind: AgentPlanCandidateKind::UiIntent,
                mutating: false,
                requires_confirmation: false,
                execution_mode: AgentPlanExecutionMode::Ask,
                shell_command: None,
                operation: None,
                workflow: None,
                ui_intent: Some(ui_intent),
            });
        }
    }

    let mut filtered = vec![];
    for candidate in candidates {
        if !allow_mutating_candidates && candidate.mutating {
            continue;
        }
        filtered.push(candidate);
        if let Some(limit) = max_candidates {
            if filtered.len() >= limit {
                break;
            }
        }
    }

    let result = AgentPlanResult {
        schema: AGENT_PLAN_RESULT_SCHEMA.to_string(),
        assistant_message:
            "Builtin echo planner compiled your prose into deterministic plan candidates."
                .to_string(),
        questions: vec![],
        candidates: filtered,
    };
    result.validate()?;
    Ok(result)
}

fn env_overrides_from_options(
    base_url_override: Option<&str>,
    model_override: Option<&str>,
    timeout_seconds: Option<u64>,
    connect_timeout_seconds: Option<u64>,
    read_timeout_seconds: Option<u64>,
    max_retries: Option<usize>,
    max_response_bytes: Option<usize>,
) -> HashMap<String, String> {
    let mut env_overrides = HashMap::new();
    if let Some(base_url) = base_url_override
        .map(str::trim)
        .filter(|value| !value.is_empty())
    {
        env_overrides.insert(AGENT_BASE_URL_ENV.to_string(), base_url.to_string());
    }
    if let Some(model) = model_override
        .map(str::trim)
        .filter(|value| !value.is_empty())
    {
        env_overrides.insert(AGENT_MODEL_ENV.to_string(), model.to_string());
    }
    if let Some(timeout_seconds) = timeout_seconds.filter(|value| *value > 0) {
        env_overrides.insert(
            AGENT_TIMEOUT_SECS_ENV.to_string(),
            timeout_seconds.to_string(),
        );
    }
    if let Some(connect_timeout_seconds) = connect_timeout_seconds.filter(|value| *value > 0) {
        env_overrides.insert(
            AGENT_CONNECT_TIMEOUT_SECS_ENV.to_string(),
            connect_timeout_seconds.to_string(),
        );
    }
    if let Some(read_timeout_seconds) = read_timeout_seconds.filter(|value| *value > 0) {
        env_overrides.insert(
            AGENT_READ_TIMEOUT_SECS_ENV.to_string(),
            read_timeout_seconds.to_string(),
        );
    }
    if let Some(max_retries) = max_retries {
        env_overrides.insert(AGENT_MAX_RETRIES_ENV.to_string(), max_retries.to_string());
    }
    if let Some(max_response_bytes) = max_response_bytes.filter(|value| *value > 0) {
        env_overrides.insert(
            AGENT_MAX_RESPONSE_BYTES_ENV.to_string(),
            max_response_bytes.to_string(),
        );
    }
    env_overrides
}

/// Compile prose into a typed machine-facing plan.
pub fn invoke_agent_plan_with_env_overrides(
    catalog_path: Option<&str>,
    system_id: &str,
    prompt: &str,
    state_summary: Option<&EngineStateSummary>,
    env_overrides: Option<&HashMap<String, String>>,
    max_candidates: Option<usize>,
    allow_mutating_candidates: bool,
) -> Result<AgentPlanResult, String> {
    let (resolved_catalog_path, catalog) = load_agent_system_catalog(catalog_path)?;
    let mut system = catalog.resolve_system(system_id)?;
    if let Some(overrides) = env_overrides {
        for (key, value) in overrides {
            let key = key.trim();
            if key.is_empty() {
                continue;
            }
            system.env.insert(key.to_string(), value.to_string());
        }
    }

    let mut plan = if matches!(system.transport, AgentSystemTransport::BuiltinEcho) {
        builtin_echo_plan_result(prompt, max_candidates, allow_mutating_candidates)?
    } else {
        let invocation = invoke_agent_support_with_env_overrides(
            Some(resolved_catalog_path.as_str()),
            system_id,
            prompt,
            state_summary,
            if system.env.is_empty() {
                None
            } else {
                Some(&system.env)
            },
        )?;
        let mut candidates = vec![];
        for suggestion in invocation.response.suggested_commands.iter() {
            let index_1based = candidates.len() + 1;
            let candidate = shell_candidate_from_command(
                index_1based,
                suggestion.title.as_deref(),
                suggestion.rationale.as_deref(),
                &suggestion.command,
                AgentPlanExecutionMode::from_agent_execution(suggestion.execution),
            )?;
            if !allow_mutating_candidates && candidate.mutating {
                continue;
            }
            candidates.push(candidate);
            if let Some(limit) = max_candidates {
                if candidates.len() >= limit {
                    break;
                }
            }
        }
        AgentPlanResult {
            schema: AGENT_PLAN_RESULT_SCHEMA.to_string(),
            assistant_message: invocation.response.assistant_message,
            questions: invocation.response.questions,
            candidates,
        }
    };
    if let Some(limit) = max_candidates {
        plan.candidates.truncate(limit);
    }
    plan.validate()?;
    Ok(plan)
}

/// Build a plan request envelope from shell/adapter inputs.
pub fn build_plan_request(
    system_id: &str,
    prompt: &str,
    state_summary: Option<&EngineStateSummary>,
    max_candidates: Option<usize>,
    allow_mutating_candidates: bool,
) -> Result<AgentPlanRequest, String> {
    let request = AgentPlanRequest {
        schema: AGENT_PLAN_REQUEST_SCHEMA.to_string(),
        system_id: system_id.trim().to_string(),
        prompt: prompt.to_string(),
        state_summary: state_summary.cloned(),
        max_candidates,
        allow_mutating_candidates: Some(allow_mutating_candidates),
    };
    request.validate()?;
    Ok(request)
}

/// Convenience helper used by CLI/shared-shell routes.
pub fn plan_from_shell_options(
    catalog_path: Option<&str>,
    system_id: &str,
    prompt: &str,
    state_summary: Option<&EngineStateSummary>,
    base_url_override: Option<&str>,
    model_override: Option<&str>,
    timeout_seconds: Option<u64>,
    connect_timeout_seconds: Option<u64>,
    read_timeout_seconds: Option<u64>,
    max_retries: Option<usize>,
    max_response_bytes: Option<usize>,
    max_candidates: Option<usize>,
    allow_mutating_candidates: bool,
) -> Result<AgentPlanResult, String> {
    let env_overrides = env_overrides_from_options(
        base_url_override,
        model_override,
        timeout_seconds,
        connect_timeout_seconds,
        read_timeout_seconds,
        max_retries,
        max_response_bytes,
    );
    invoke_agent_plan_with_env_overrides(
        catalog_path,
        system_id,
        prompt,
        state_summary,
        if env_overrides.is_empty() {
            None
        } else {
            Some(&env_overrides)
        },
        max_candidates,
        allow_mutating_candidates,
    )
}

/// Parse a stored plan JSON string.
pub fn parse_agent_plan_json(raw_json: &str) -> Result<AgentPlanResult, String> {
    let plan = serde_json::from_str::<AgentPlanResult>(raw_json)
        .map_err(|err| format!("Invalid planner result JSON: {err}"))?;
    plan.validate()?;
    Ok(plan)
}

/// Load a stored plan from raw JSON or an `@path` file reference.
pub fn load_agent_plan_from_argument(raw: &str) -> Result<AgentPlanResult, String> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return Err("planner result JSON or @path cannot be empty".to_string());
    }
    if let Some(path) = trimmed.strip_prefix('@') {
        let text = fs::read_to_string(path)
            .map_err(|err| format!("Could not read planner result file '{path}': {err}"))?;
        return parse_agent_plan_json(&text);
    }
    parse_agent_plan_json(trimmed)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builtin_echo_planner_supports_multiple_candidate_kinds() {
        let plan = builtin_echo_plan_result(
            "auto: state-summary\nop: {\"SetParameter\":{\"name\":\"x\",\"value\":1}}\nui: {\"action\":\"open\",\"target\":\"agent-assistant\"}",
            None,
            true,
        )
        .expect("builtin plan");
        assert_eq!(plan.candidates.len(), 3);
        assert_eq!(plan.candidates[0].kind, AgentPlanCandidateKind::Shell);
        assert_eq!(plan.candidates[1].kind, AgentPlanCandidateKind::Op);
        assert_eq!(plan.candidates[2].kind, AgentPlanCandidateKind::UiIntent);
    }

    #[test]
    fn planner_candidate_validate_rejects_multiple_payloads() {
        let candidate = AgentPlanCandidate {
            candidate_id: "candidate-1".to_string(),
            title: "bad".to_string(),
            kind: AgentPlanCandidateKind::Shell,
            shell_command: Some("state-summary".to_string()),
            operation: Some(serde_json::json!({"SetParameter":{"name":"x","value":1}})),
            ..Default::default()
        };
        let err = candidate.validate().expect_err("multiple payloads");
        assert!(err.contains("exactly one payload field"));
    }
}
