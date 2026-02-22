use crate::engine::EngineStateSummary;
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::{
    collections::{HashMap, HashSet},
    fs,
    io::{ErrorKind, Write},
    path::{Path, PathBuf},
    process::{Command, Output, Stdio},
    thread,
    time::{Duration, SystemTime, UNIX_EPOCH},
};

pub const DEFAULT_AGENT_SYSTEM_CATALOG_PATH: &str = "assets/agent_systems.json";
const AGENT_SYSTEMS_SCHEMA: &str = "gentle.agent_systems.v1";
const AGENT_REQUEST_SCHEMA: &str = "gentle.agent_request.v1";
const AGENT_RESPONSE_SCHEMA: &str = "gentle.agent_response.v1";
const AGENT_SYSTEMS_SCHEMA_PREFIX: &str = "gentle.agent_systems.v";
const AGENT_REQUEST_SCHEMA_PREFIX: &str = "gentle.agent_request.v";
const AGENT_RESPONSE_SCHEMA_PREFIX: &str = "gentle.agent_response.v";
const AGENT_SCHEMA_SUPPORTED_MAJOR: u32 = 1;
const AGENT_INVOKE_RETRY_ATTEMPTS: usize = 3;
const AGENT_INVOKE_RETRY_BASE_DELAY_MS: u64 = 250;
const OPENAI_DEFAULT_MODEL: &str = "gpt-5";
const OPENAI_DEFAULT_BASE_URL: &str = "https://api.openai.com/v1";
const OPENAI_COMPAT_DEFAULT_MODEL: &str = "llama3.1";
const OPENAI_COMPAT_DEFAULT_BASE_URL: &str = "http://127.0.0.1:11434/v1";
pub const OPENAI_API_KEY_ENV: &str = "OPENAI_API_KEY";
pub const AGENT_BASE_URL_ENV: &str = "GENTLE_AGENT_BASE_URL";
pub const AGENT_MODEL_ENV: &str = "GENTLE_AGENT_MODEL";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AgentBridgeErrorCode {
    InvalidInput,
    CatalogRead,
    CatalogParse,
    CatalogValidation,
    SchemaValidation,
    SchemaUnsupported,
    AdapterUnavailable,
    AdapterTransient,
    AdapterFailed,
    ResponseParse,
    ResponseValidation,
}

impl AgentBridgeErrorCode {
    fn as_str(self) -> &'static str {
        match self {
            Self::InvalidInput => "INVALID_INPUT",
            Self::CatalogRead => "CATALOG_READ",
            Self::CatalogParse => "CATALOG_PARSE",
            Self::CatalogValidation => "CATALOG_VALIDATION",
            Self::SchemaValidation => "SCHEMA_VALIDATION",
            Self::SchemaUnsupported => "SCHEMA_UNSUPPORTED",
            Self::AdapterUnavailable => "ADAPTER_UNAVAILABLE",
            Self::AdapterTransient => "ADAPTER_TRANSIENT",
            Self::AdapterFailed => "ADAPTER_FAILED",
            Self::ResponseParse => "RESPONSE_PARSE",
            Self::ResponseValidation => "RESPONSE_VALIDATION",
        }
    }
}

fn agent_err(code: AgentBridgeErrorCode, message: impl Into<String>) -> String {
    format!("AGENT_{}: {}", code.as_str(), message.into())
}

fn parse_schema_major(schema: &str, prefix: &str) -> Option<u32> {
    schema
        .trim()
        .strip_prefix(prefix)
        .and_then(|raw| raw.parse::<u32>().ok())
}

fn require_supported_schema(schema: &str, prefix: &str, context: &str) -> Result<u32, String> {
    let Some(major) = parse_schema_major(schema, prefix) else {
        return Err(agent_err(
            AgentBridgeErrorCode::SchemaValidation,
            format!(
                "{context} schema '{}' is invalid (expected '{}{}')",
                schema, prefix, AGENT_SCHEMA_SUPPORTED_MAJOR
            ),
        ));
    };
    if major != AGENT_SCHEMA_SUPPORTED_MAJOR {
        return Err(agent_err(
            AgentBridgeErrorCode::SchemaUnsupported,
            format!(
                "{context} schema major version {} is not supported (supported: {})",
                major, AGENT_SCHEMA_SUPPORTED_MAJOR
            ),
        ));
    }
    Ok(major)
}

fn is_extension_key(key: &str) -> bool {
    key.starts_with("x_") || key.starts_with("x-")
}

fn ensure_known_keys(
    object: &serde_json::Map<String, Value>,
    known: &[&str],
    context: &str,
    code: AgentBridgeErrorCode,
) -> Result<(), String> {
    let known_set = known.iter().copied().collect::<HashSet<_>>();
    let mut allowed_fields = known
        .iter()
        .map(|value| value.to_string())
        .collect::<Vec<_>>();
    allowed_fields.sort();
    for key in object.keys() {
        if !known_set.contains(key.as_str()) && !is_extension_key(key) {
            return Err(agent_err(
                code,
                format!(
                    "{context} contains unsupported field '{}' (allowed: {}, x_/x- extensions)",
                    key,
                    allowed_fields.join(", ")
                ),
            ));
        }
    }
    Ok(())
}

fn has_path_separator(value: &str) -> bool {
    value.contains(std::path::MAIN_SEPARATOR) || value.contains('/')
}

fn is_executable_file(path: &Path) -> bool {
    if !path.is_file() {
        return false;
    }
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        return fs::metadata(path)
            .map(|meta| meta.permissions().mode() & 0o111 != 0)
            .unwrap_or(false);
    }
    #[cfg(not(unix))]
    {
        true
    }
}

fn resolve_executable_path(program: &str) -> Option<PathBuf> {
    let trimmed = program.trim();
    if trimmed.is_empty() {
        return None;
    }
    if has_path_separator(trimmed) {
        let candidate = PathBuf::from(trimmed);
        return is_executable_file(&candidate).then_some(candidate);
    }
    let path_var = std::env::var_os("PATH")?;
    for entry in std::env::split_paths(&path_var) {
        let candidate = entry.join(trimmed);
        if is_executable_file(&candidate) {
            return Some(candidate);
        }
    }
    None
}

fn resolve_api_key(system: &AgentSystemSpec) -> Option<String> {
    system
        .env
        .get(OPENAI_API_KEY_ENV)
        .map(|value| value.trim().to_string())
        .filter(|value| !value.is_empty())
        .or_else(|| {
            std::env::var(OPENAI_API_KEY_ENV)
                .ok()
                .map(|value| value.trim().to_string())
                .filter(|value| !value.is_empty())
        })
}

fn resolve_model(system: &AgentSystemSpec, default: &str) -> String {
    system
        .env
        .get(AGENT_MODEL_ENV)
        .map(|value| value.trim().to_string())
        .filter(|value| !value.is_empty())
        .or_else(|| {
            system
                .model
                .as_deref()
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(str::to_string)
        })
        .unwrap_or_else(|| default.to_string())
}

fn normalize_base_url(raw: &str) -> Option<String> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return None;
    }
    let with_scheme = if trimmed.contains("://") {
        trimmed.to_string()
    } else {
        format!("http://{trimmed}")
    };
    Some(with_scheme.trim_end_matches('/').to_string())
}

fn resolve_base_url(system: &AgentSystemSpec, default: &str) -> String {
    system
        .env
        .get(AGENT_BASE_URL_ENV)
        .and_then(|value| normalize_base_url(value))
        .or_else(|| {
            system
                .base_url
                .as_deref()
                .and_then(|value| normalize_base_url(value))
        })
        .or_else(|| normalize_base_url(default))
        .unwrap_or_else(|| default.trim_end_matches('/').to_string())
}

fn openai_compat_endpoint_candidates(base_url: &str) -> Vec<String> {
    let normalized = normalize_base_url(base_url)
        .or_else(|| normalize_base_url(OPENAI_COMPAT_DEFAULT_BASE_URL))
        .unwrap_or_else(|| {
            OPENAI_COMPAT_DEFAULT_BASE_URL
                .trim_end_matches('/')
                .to_string()
        });
    let mut endpoints = vec![format!("{normalized}/chat/completions")];
    if !normalized.ends_with("/v1") {
        let v1 = format!("{normalized}/v1/chat/completions");
        if !endpoints.contains(&v1) {
            endpoints.push(v1);
        }
    }
    endpoints
}

fn openai_compat_endpoint_candidates_for_system(system: &AgentSystemSpec) -> Vec<String> {
    let primary_base = resolve_base_url(system, OPENAI_COMPAT_DEFAULT_BASE_URL);
    let has_explicit_override = system
        .env
        .get(AGENT_BASE_URL_ENV)
        .map(|value| !value.trim().is_empty())
        .unwrap_or(false);
    let mut endpoints = openai_compat_endpoint_candidates(&primary_base);
    if has_explicit_override {
        return endpoints;
    }
    // For local OpenAI-compatible services, probe common host/port defaults
    // when catalog defaults do not match the running daemon.
    let local_base_candidates = [
        "http://localhost:11964",
        "http://127.0.0.1:11964",
        "http://localhost:11434",
        "http://127.0.0.1:11434",
    ];
    for base in local_base_candidates {
        for endpoint in openai_compat_endpoint_candidates(base) {
            if !endpoints.contains(&endpoint) {
                endpoints.push(endpoint);
            }
        }
    }
    endpoints
}

fn openai_model_list_endpoint_candidates(base_url: &str) -> Vec<String> {
    let normalized = normalize_base_url(base_url)
        .or_else(|| normalize_base_url(OPENAI_COMPAT_DEFAULT_BASE_URL))
        .unwrap_or_else(|| {
            OPENAI_COMPAT_DEFAULT_BASE_URL
                .trim_end_matches('/')
                .to_string()
        });
    let mut endpoints = vec![format!("{normalized}/models")];
    if !normalized.ends_with("/v1") {
        let v1 = format!("{normalized}/v1/models");
        if !endpoints.contains(&v1) {
            endpoints.push(v1);
        }
    }
    endpoints
}

fn extract_models_from_openai_models_payload(value: &Value) -> Vec<String> {
    let mut out: Vec<String> = vec![];
    let mut seen: HashSet<String> = HashSet::new();
    let mut push = |raw: &str| {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return;
        }
        if seen.insert(trimmed.to_string()) {
            out.push(trimmed.to_string());
        }
    };
    if let Some(data) = value.get("data").and_then(Value::as_array) {
        for item in data {
            if let Some(id) = item.get("id").and_then(Value::as_str) {
                push(id);
            }
        }
    }
    if let Some(models) = value.get("models").and_then(Value::as_array) {
        for item in models {
            if let Some(id) = item.get("id").and_then(Value::as_str) {
                push(id);
                continue;
            }
            if let Some(name) = item.get("name").and_then(Value::as_str) {
                push(name);
            }
        }
    }
    if let Some(items) = value.as_array() {
        for item in items {
            if let Some(id) = item.as_str() {
                push(id);
            }
        }
    }
    out
}

pub fn discover_openai_models(
    base_url: &str,
    api_key: Option<&str>,
) -> Result<Vec<String>, String> {
    let endpoints = openai_model_list_endpoint_candidates(base_url);
    let client = reqwest::blocking::Client::builder()
        .timeout(Duration::from_secs(20))
        .build()
        .map_err(|e| format!("could not build OpenAI model-discovery client: {e}"))?;
    let mut first_error: Option<String> = None;
    for endpoint in endpoints {
        let mut request = client.get(&endpoint);
        if let Some(key) = api_key.map(str::trim).filter(|value| !value.is_empty()) {
            request = request.bearer_auth(key);
        }
        match request.send() {
            Ok(response) => {
                let status = response.status();
                let body = response
                    .text()
                    .map_err(|e| format!("could not read model discovery response body: {e}"))?;
                if !status.is_success() {
                    if (status.as_u16() == 404 || status.as_u16() == 405) && first_error.is_none() {
                        first_error = Some(format!(
                            "model discovery endpoint not supported at {endpoint} (status={status})"
                        ));
                        continue;
                    }
                    return Err(format!(
                        "model discovery failed at {endpoint} (status={status}): {}",
                        body.trim()
                    ));
                }
                let value = serde_json::from_str::<Value>(&body).map_err(|e| {
                    format!("model discovery endpoint returned invalid JSON at {endpoint}: {e}")
                })?;
                let models = extract_models_from_openai_models_payload(&value);
                if models.is_empty() {
                    return Err(format!(
                        "model discovery at {endpoint} succeeded but returned no model ids"
                    ));
                }
                return Ok(models);
            }
            Err(e) => {
                if first_error.is_none() {
                    first_error = Some(format!("request failed at {endpoint}: {e}"));
                }
            }
        }
    }
    Err(first_error
        .unwrap_or_else(|| "model discovery failed: no endpoint candidate succeeded".to_string()))
}

pub fn agent_system_availability(system: &AgentSystemSpec) -> AgentSystemAvailability {
    match system.transport {
        AgentSystemTransport::BuiltinEcho => AgentSystemAvailability {
            available: true,
            reason: None,
        },
        AgentSystemTransport::ExternalJsonStdio => {
            if system.command.is_empty() {
                return AgentSystemAvailability {
                    available: false,
                    reason: Some("no command configured".to_string()),
                };
            }
            let program = system.command[0].trim();
            if program.is_empty() {
                return AgentSystemAvailability {
                    available: false,
                    reason: Some("empty command executable".to_string()),
                };
            }
            match resolve_executable_path(program) {
                Some(path) => AgentSystemAvailability {
                    available: true,
                    reason: Some(format!("executable found at {}", path.display())),
                },
                None => AgentSystemAvailability {
                    available: false,
                    reason: Some(format!("executable '{}' not found on PATH", program)),
                },
            }
        }
        AgentSystemTransport::NativeOpenai => {
            if resolve_api_key(system).is_none() {
                return AgentSystemAvailability {
                    available: false,
                    reason: Some(format!("{OPENAI_API_KEY_ENV} is not set")),
                };
            }
            AgentSystemAvailability {
                available: true,
                reason: None,
            }
        }
        AgentSystemTransport::NativeOpenaiCompat => {
            let base_url = resolve_base_url(system, OPENAI_COMPAT_DEFAULT_BASE_URL);
            let endpoints = openai_compat_endpoint_candidates_for_system(system);
            let model = resolve_model(system, OPENAI_COMPAT_DEFAULT_MODEL);
            let preview = endpoints
                .iter()
                .take(2)
                .cloned()
                .collect::<Vec<_>>()
                .join(" or ");
            let extra = endpoints.len().saturating_sub(2);
            AgentSystemAvailability {
                available: true,
                reason: Some(format!(
                    "expects OpenAI-compatible local server near {} (model '{}'; tries {}{}; key optional)",
                    base_url,
                    model,
                    preview,
                    if extra > 0 {
                        format!(" + {} fallback endpoint(s)", extra)
                    } else {
                        String::new()
                    }
                )),
            }
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum AgentSystemTransport {
    #[default]
    ExternalJsonStdio,
    NativeOpenai,
    NativeOpenaiCompat,
    BuiltinEcho,
}

impl AgentSystemTransport {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::ExternalJsonStdio => "external_json_stdio",
            Self::NativeOpenai => "native_openai",
            Self::NativeOpenaiCompat => "native_openai_compat",
            Self::BuiltinEcho => "builtin_echo",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct AgentSystemSpec {
    pub id: String,
    pub label: String,
    pub description: Option<String>,
    pub transport: AgentSystemTransport,
    pub command: Vec<String>,
    pub model: Option<String>,
    pub base_url: Option<String>,
    pub env: HashMap<String, String>,
    pub working_dir: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct AgentSystemAvailability {
    pub available: bool,
    pub reason: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct AgentSystemCatalog {
    pub schema: String,
    pub systems: Vec<AgentSystemSpec>,
}

impl AgentSystemCatalog {
    pub fn effective_catalog_path(catalog_path: Option<&str>) -> String {
        catalog_path
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .unwrap_or(DEFAULT_AGENT_SYSTEM_CATALOG_PATH)
            .to_string()
    }

    pub fn from_json_file(path: &str) -> Result<Self, String> {
        let text = fs::read_to_string(path).map_err(|e| {
            agent_err(
                AgentBridgeErrorCode::CatalogRead,
                format!("could not read agent catalog '{}': {}", path, e),
            )
        })?;
        let mut catalog = serde_json::from_str::<Self>(&text).map_err(|e| {
            agent_err(
                AgentBridgeErrorCode::CatalogParse,
                format!("could not parse agent catalog '{}': {}", path, e),
            )
        })?;
        if catalog.schema.trim().is_empty() {
            catalog.schema = AGENT_SYSTEMS_SCHEMA.to_string();
        }
        require_supported_schema(
            &catalog.schema,
            AGENT_SYSTEMS_SCHEMA_PREFIX,
            "agent catalog",
        )?;
        let mut normalized: Vec<AgentSystemSpec> = vec![];
        let mut seen_ids = HashSet::new();
        for mut system in catalog.systems {
            let id = system.id.trim().to_string();
            if id.is_empty() {
                return Err(agent_err(
                    AgentBridgeErrorCode::CatalogValidation,
                    "agent catalog contains a system with empty 'id'",
                ));
            }
            if !seen_ids.insert(id.clone()) {
                return Err(agent_err(
                    AgentBridgeErrorCode::CatalogValidation,
                    format!("agent catalog has duplicate system id '{}'", id),
                ));
            }
            system.id = id;
            if system.label.trim().is_empty() {
                system.label = system.id.clone();
            }
            if matches!(system.transport, AgentSystemTransport::ExternalJsonStdio) {
                if system.command.is_empty() {
                    return Err(agent_err(
                        AgentBridgeErrorCode::CatalogValidation,
                        format!(
                            "agent system '{}' uses external_json_stdio but has no command",
                            system.id
                        ),
                    ));
                }
                if system.command.iter().any(|part| part.trim().is_empty()) {
                    return Err(agent_err(
                        AgentBridgeErrorCode::CatalogValidation,
                        format!(
                            "agent system '{}' command contains an empty token",
                            system.id
                        ),
                    ));
                }
            }
            normalized.push(system);
        }
        normalized.sort_by(|left, right| left.id.cmp(&right.id));
        catalog.systems = normalized;
        Ok(catalog)
    }

    pub fn resolve_system(&self, system_id: &str) -> Result<AgentSystemSpec, String> {
        let id = system_id.trim();
        if id.is_empty() {
            return Err(agent_err(
                AgentBridgeErrorCode::InvalidInput,
                "agent system ID cannot be empty",
            ));
        }
        self.systems
            .iter()
            .find(|system| system.id == id)
            .cloned()
            .ok_or_else(|| {
                let mut known = self
                    .systems
                    .iter()
                    .map(|system| system.id.clone())
                    .collect::<Vec<_>>();
                known.sort();
                agent_err(
                    AgentBridgeErrorCode::CatalogValidation,
                    format!(
                        "agent system '{}' not found in catalog (available: {})",
                        id,
                        if known.is_empty() {
                            "none".to_string()
                        } else {
                            known.join(", ")
                        }
                    ),
                )
            })
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
struct AgentRequestPayload {
    schema: String,
    system_id: String,
    prompt: String,
    sent_at_unix_ms: u128,
    state_summary: Option<EngineStateSummary>,
}

impl Default for AgentRequestPayload {
    fn default() -> Self {
        Self {
            schema: AGENT_REQUEST_SCHEMA.to_string(),
            system_id: String::new(),
            prompt: String::new(),
            sent_at_unix_ms: 0,
            state_summary: None,
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum AgentExecutionIntent {
    Chat,
    #[default]
    Ask,
    Auto,
}

impl AgentExecutionIntent {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Chat => "chat",
            Self::Ask => "ask",
            Self::Auto => "auto",
        }
    }

    pub fn parse(raw: &str) -> Option<Self> {
        match raw.trim().to_ascii_lowercase().as_str() {
            "chat" | "question" | "none" => Some(Self::Chat),
            "ask" | "confirm" | "confirmation" => Some(Self::Ask),
            "auto" | "autorun" | "run" => Some(Self::Auto),
            _ => None,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct AgentSuggestedCommand {
    pub title: Option<String>,
    pub rationale: Option<String>,
    pub command: String,
    pub execution: AgentExecutionIntent,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct AgentResponse {
    pub schema: String,
    pub assistant_message: String,
    pub questions: Vec<String>,
    pub suggested_commands: Vec<AgentSuggestedCommand>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AgentInvocationOutcome {
    pub catalog_path: String,
    pub system_id: String,
    pub system_label: String,
    pub transport: String,
    pub command: Vec<String>,
    pub request: Value,
    pub response: AgentResponse,
    pub raw_stdout: String,
    pub raw_stderr: String,
    pub exit_code: Option<i32>,
    pub elapsed_ms: u128,
}

fn now_unix_ms() -> u128 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_millis())
        .unwrap_or(0)
}

fn build_agent_request(
    system_id: &str,
    prompt: &str,
    state_summary: Option<&EngineStateSummary>,
) -> Result<(AgentRequestPayload, Value, String), String> {
    let payload = AgentRequestPayload {
        schema: AGENT_REQUEST_SCHEMA.to_string(),
        system_id: system_id.trim().to_string(),
        prompt: prompt.to_string(),
        sent_at_unix_ms: now_unix_ms(),
        state_summary: state_summary.cloned(),
    };
    let request_value = serde_json::to_value(&payload).map_err(|e| {
        agent_err(
            AgentBridgeErrorCode::SchemaValidation,
            format!("could not serialize agent request payload: {e}"),
        )
    })?;
    validate_agent_request_value(&request_value)?;
    let request_json = serde_json::to_string(&payload).map_err(|e| {
        agent_err(
            AgentBridgeErrorCode::SchemaValidation,
            format!("could not encode agent request JSON: {e}"),
        )
    })?;
    Ok((payload, request_value, request_json))
}

fn validate_agent_request_value(value: &Value) -> Result<(), String> {
    let Some(object) = value.as_object() else {
        return Err(agent_err(
            AgentBridgeErrorCode::SchemaValidation,
            "agent request must be a JSON object",
        ));
    };
    ensure_known_keys(
        object,
        &[
            "schema",
            "system_id",
            "prompt",
            "sent_at_unix_ms",
            "state_summary",
        ],
        "agent request",
        AgentBridgeErrorCode::SchemaValidation,
    )?;
    let schema = object
        .get("schema")
        .and_then(Value::as_str)
        .ok_or_else(|| {
            agent_err(
                AgentBridgeErrorCode::SchemaValidation,
                "agent request 'schema' must be a string",
            )
        })?;
    require_supported_schema(schema, AGENT_REQUEST_SCHEMA_PREFIX, "agent request")?;
    let system_id = object
        .get("system_id")
        .and_then(Value::as_str)
        .map(str::trim)
        .ok_or_else(|| {
            agent_err(
                AgentBridgeErrorCode::SchemaValidation,
                "agent request 'system_id' must be a string",
            )
        })?;
    if system_id.is_empty() {
        return Err(agent_err(
            AgentBridgeErrorCode::SchemaValidation,
            "agent request 'system_id' cannot be empty",
        ));
    }
    let prompt = object
        .get("prompt")
        .and_then(Value::as_str)
        .ok_or_else(|| {
            agent_err(
                AgentBridgeErrorCode::SchemaValidation,
                "agent request 'prompt' must be a string",
            )
        })?;
    if prompt.trim().is_empty() {
        return Err(agent_err(
            AgentBridgeErrorCode::SchemaValidation,
            "agent request 'prompt' cannot be empty",
        ));
    }
    if object
        .get("sent_at_unix_ms")
        .and_then(Value::as_u64)
        .is_none()
    {
        return Err(agent_err(
            AgentBridgeErrorCode::SchemaValidation,
            "agent request 'sent_at_unix_ms' must be an unsigned integer",
        ));
    }
    if let Some(state_summary) = object.get("state_summary") {
        if !(state_summary.is_null() || state_summary.is_object()) {
            return Err(agent_err(
                AgentBridgeErrorCode::SchemaValidation,
                "agent request 'state_summary' must be object or null",
            ));
        }
    }
    Ok(())
}

fn replace_placeholders(
    token: &str,
    prompt: &str,
    request_json: &str,
    state_summary_json: &str,
) -> String {
    token
        .replace("{prompt}", prompt)
        .replace("{request_json}", request_json)
        .replace("{state_summary_json}", state_summary_json)
}

fn parse_questions_field(value: Option<&Value>) -> Result<Vec<String>, String> {
    let Some(value) = value else {
        return Ok(vec![]);
    };
    let Value::Array(items) = value else {
        return Err(agent_err(
            AgentBridgeErrorCode::ResponseValidation,
            "agent response 'questions' must be an array of strings",
        ));
    };
    let mut parsed = Vec::with_capacity(items.len());
    for (idx, item) in items.iter().enumerate() {
        let question = item.as_str().ok_or_else(|| {
            agent_err(
                AgentBridgeErrorCode::ResponseValidation,
                format!("agent response 'questions[{idx}]' must be a string"),
            )
        })?;
        let trimmed = question.trim();
        if !trimmed.is_empty() {
            parsed.push(trimmed.to_string());
        }
    }
    Ok(parsed)
}

fn parse_suggested_command_value(
    value: &Value,
    idx: usize,
) -> Result<AgentSuggestedCommand, String> {
    let Value::Object(obj) = value else {
        return Err(agent_err(
            AgentBridgeErrorCode::ResponseValidation,
            format!("agent response 'suggested_commands[{idx}]' must be an object"),
        ));
    };
    ensure_known_keys(
        obj,
        &["title", "rationale", "command", "execution"],
        "agent suggested command",
        AgentBridgeErrorCode::ResponseValidation,
    )?;

    let command = obj
        .get("command")
        .and_then(Value::as_str)
        .map(str::trim)
        .ok_or_else(|| {
            agent_err(
                AgentBridgeErrorCode::ResponseValidation,
                format!("agent response 'suggested_commands[{idx}].command' must be a string"),
            )
        })?;
    if command.is_empty() {
        return Err(agent_err(
            AgentBridgeErrorCode::ResponseValidation,
            format!("agent response 'suggested_commands[{idx}].command' cannot be empty"),
        ));
    }
    let execution = if let Some(value) = obj.get("execution") {
        let raw = value.as_str().ok_or_else(|| {
            agent_err(
                AgentBridgeErrorCode::ResponseValidation,
                format!("agent response 'suggested_commands[{idx}].execution' must be a string"),
            )
        })?;
        AgentExecutionIntent::parse(raw).ok_or_else(|| {
            agent_err(
                AgentBridgeErrorCode::ResponseValidation,
                format!(
                    "agent response 'suggested_commands[{idx}].execution' has unsupported value '{}'",
                    raw
                ),
            )
        })?
    } else {
        AgentExecutionIntent::Ask
    };
    let title = if let Some(value) = obj.get("title") {
        Some(
            value
                .as_str()
                .ok_or_else(|| {
                    agent_err(
                        AgentBridgeErrorCode::ResponseValidation,
                        format!(
                            "agent response 'suggested_commands[{idx}].title' must be a string"
                        ),
                    )
                })?
                .trim()
                .to_string(),
        )
        .filter(|value| !value.is_empty())
    } else {
        None
    };
    let rationale = if let Some(value) = obj.get("rationale") {
        Some(
            value
                .as_str()
                .ok_or_else(|| {
                    agent_err(
                        AgentBridgeErrorCode::ResponseValidation,
                        format!(
                            "agent response 'suggested_commands[{idx}].rationale' must be a string"
                        ),
                    )
                })?
                .trim()
                .to_string(),
        )
        .filter(|value| !value.is_empty())
    } else {
        None
    };

    Ok(AgentSuggestedCommand {
        title,
        rationale,
        command: command.to_string(),
        execution,
    })
}

fn parse_suggested_commands_field(
    value: Option<&Value>,
) -> Result<Vec<AgentSuggestedCommand>, String> {
    let Some(value) = value else {
        return Ok(vec![]);
    };
    let Value::Array(items) = value else {
        return Err(agent_err(
            AgentBridgeErrorCode::ResponseValidation,
            "agent response 'suggested_commands' must be an array of command objects",
        ));
    };
    let mut parsed = Vec::with_capacity(items.len());
    for (idx, item) in items.iter().enumerate() {
        parsed.push(parse_suggested_command_value(item, idx)?);
    }
    Ok(parsed)
}

fn parse_agent_response_value(value: Value) -> Result<AgentResponse, String> {
    let Some(obj) = value.as_object() else {
        return Err(agent_err(
            AgentBridgeErrorCode::ResponseValidation,
            "agent response must be a JSON object",
        ));
    };
    ensure_known_keys(
        obj,
        &[
            "schema",
            "assistant_message",
            "questions",
            "suggested_commands",
        ],
        "agent response",
        AgentBridgeErrorCode::ResponseValidation,
    )?;

    let schema = obj
        .get("schema")
        .and_then(Value::as_str)
        .ok_or_else(|| {
            agent_err(
                AgentBridgeErrorCode::ResponseValidation,
                "agent response 'schema' must be a string",
            )
        })?
        .trim()
        .to_string();
    require_supported_schema(&schema, AGENT_RESPONSE_SCHEMA_PREFIX, "agent response")?;

    let assistant_message = obj
        .get("assistant_message")
        .and_then(Value::as_str)
        .ok_or_else(|| {
            agent_err(
                AgentBridgeErrorCode::ResponseValidation,
                "agent response 'assistant_message' must be a string",
            )
        })?
        .trim()
        .to_string();
    let questions = parse_questions_field(obj.get("questions"))?;
    let suggested_commands = parse_suggested_commands_field(obj.get("suggested_commands"))?;

    if assistant_message.is_empty() && questions.is_empty() && suggested_commands.is_empty() {
        return Err(agent_err(
            AgentBridgeErrorCode::ResponseValidation,
            "agent response must include assistant_message, questions, or suggested_commands",
        ));
    }

    Ok(AgentResponse {
        schema,
        assistant_message,
        questions,
        suggested_commands,
    })
}

fn parse_agent_response(stdout: &str) -> Result<AgentResponse, String> {
    let trimmed = stdout.trim();
    if trimmed.is_empty() {
        return Err(agent_err(
            AgentBridgeErrorCode::ResponseParse,
            "agent produced empty stdout",
        ));
    }
    let value = serde_json::from_str::<Value>(trimmed).map_err(|e| {
        agent_err(
            AgentBridgeErrorCode::ResponseParse,
            format!("agent stdout was not valid JSON: {e}"),
        )
    })?;
    parse_agent_response_value(value)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ExternalAttemptErrorKind {
    Unavailable,
    Transient,
    Fatal,
}

#[derive(Debug, Clone)]
struct ExternalAttemptError {
    kind: ExternalAttemptErrorKind,
    message: String,
}

fn classify_io_error_kind(kind: ErrorKind) -> ExternalAttemptErrorKind {
    match kind {
        ErrorKind::NotFound | ErrorKind::PermissionDenied => ExternalAttemptErrorKind::Unavailable,
        ErrorKind::Interrupted
        | ErrorKind::WouldBlock
        | ErrorKind::TimedOut
        | ErrorKind::ConnectionAborted
        | ErrorKind::ConnectionReset
        | ErrorKind::NotConnected
        | ErrorKind::BrokenPipe => ExternalAttemptErrorKind::Transient,
        _ => ExternalAttemptErrorKind::Fatal,
    }
}

fn is_transient_exit_code(code: Option<i32>) -> bool {
    matches!(code, Some(75 | 104 | 110 | 111))
}

fn retry_backoff_duration(attempt: usize) -> Duration {
    let shift = (attempt.saturating_sub(1)).min(6) as u32;
    let multiplier = 1_u64 << shift;
    Duration::from_millis(AGENT_INVOKE_RETRY_BASE_DELAY_MS.saturating_mul(multiplier))
}

fn invoke_external_json_stdio_once(
    system: &AgentSystemSpec,
    rendered_command: &[String],
    request_json: &str,
) -> Result<Output, ExternalAttemptError> {
    let mut cmd = Command::new(&rendered_command[0]);
    if rendered_command.len() > 1 {
        cmd.args(&rendered_command[1..]);
    }
    if let Some(working_dir) = system
        .working_dir
        .as_deref()
        .map(str::trim)
        .filter(|value| !value.is_empty())
    {
        cmd.current_dir(working_dir);
    }
    for (key, value) in &system.env {
        cmd.env(key, value);
    }
    cmd.stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());

    let mut child = cmd.spawn().map_err(|e| ExternalAttemptError {
        kind: classify_io_error_kind(e.kind()),
        message: format!(
            "could not spawn agent command '{}': {}",
            rendered_command.join(" "),
            e
        ),
    })?;
    if let Some(mut stdin) = child.stdin.take() {
        stdin
            .write_all(request_json.as_bytes())
            .map_err(|e| ExternalAttemptError {
                kind: classify_io_error_kind(e.kind()),
                message: format!("could not write request JSON to agent stdin: {e}"),
            })?;
    }
    let output = child.wait_with_output().map_err(|e| ExternalAttemptError {
        kind: classify_io_error_kind(e.kind()),
        message: format!("could not wait for agent process: {e}"),
    })?;
    if !output.status.success() {
        let stdout = String::from_utf8_lossy(&output.stdout).to_string();
        let stderr = String::from_utf8_lossy(&output.stderr).to_string();
        let kind = if is_transient_exit_code(output.status.code()) {
            ExternalAttemptErrorKind::Transient
        } else {
            ExternalAttemptErrorKind::Fatal
        };
        return Err(ExternalAttemptError {
            kind,
            message: format!(
                "agent command failed (exit={:?}): {}{}",
                output.status.code(),
                stderr.trim(),
                if stdout.trim().is_empty() {
                    String::new()
                } else {
                    format!(" | stdout: {}", stdout.trim())
                }
            ),
        });
    }
    Ok(output)
}

fn extract_openai_output_text(response_json: &Value) -> Option<String> {
    if let Some(text) = response_json.get("output_text").and_then(Value::as_str) {
        let trimmed = text.trim();
        if !trimmed.is_empty() {
            return Some(trimmed.to_string());
        }
    }
    let output_items = response_json.get("output").and_then(Value::as_array)?;
    let mut collected = String::new();
    for item in output_items {
        let Some(content) = item.get("content").and_then(Value::as_array) else {
            continue;
        };
        for block in content {
            let text = block
                .get("text")
                .and_then(Value::as_str)
                .or_else(|| block.get("output_text").and_then(Value::as_str));
            if let Some(text) = text {
                let trimmed = text.trim();
                if !trimmed.is_empty() {
                    if !collected.is_empty() {
                        collected.push('\n');
                    }
                    collected.push_str(trimmed);
                }
            }
        }
    }
    if collected.trim().is_empty() {
        None
    } else {
        Some(collected)
    }
}

fn extract_openai_chat_completions_text(response_json: &Value) -> Option<String> {
    let choices = response_json.get("choices")?.as_array()?;
    let first = choices.first()?;
    let message = first.get("message")?;
    if let Some(content) = message.get("content").and_then(Value::as_str) {
        let trimmed = content.trim();
        if !trimmed.is_empty() {
            return Some(trimmed.to_string());
        }
    }
    let content_blocks = message.get("content").and_then(Value::as_array)?;
    let mut collected = String::new();
    for block in content_blocks {
        if let Some(text) = block.get("text").and_then(Value::as_str) {
            let trimmed = text.trim();
            if !trimmed.is_empty() {
                if !collected.is_empty() {
                    collected.push('\n');
                }
                collected.push_str(trimmed);
            }
        }
    }
    if collected.trim().is_empty() {
        None
    } else {
        Some(collected)
    }
}

fn extract_openai_error_code(body: &str) -> Option<String> {
    let value = serde_json::from_str::<Value>(body).ok()?;
    let error = value.get("error")?;
    error
        .get("code")
        .and_then(Value::as_str)
        .or_else(|| error.get("type").and_then(Value::as_str))
        .map(|value| value.trim().to_string())
        .filter(|value| !value.is_empty())
}

fn classify_openai_http_error(status: reqwest::StatusCode, body: &str) -> ExternalAttemptError {
    let message_prefix = format!("OpenAI API error (status={}): {}", status, body.trim());
    let error_code = extract_openai_error_code(body).unwrap_or_default();
    if status.as_u16() == 429 && error_code.eq_ignore_ascii_case("insufficient_quota") {
        return ExternalAttemptError {
            kind: ExternalAttemptErrorKind::Unavailable,
            message: format!(
                "{message_prefix}\nHint: OpenAI reported insufficient quota. Check API project billing/usage at https://platform.openai.com/usage and https://platform.openai.com/settings/organization/billing/overview ."
            ),
        };
    }
    let kind = if status.is_server_error() || status.as_u16() == 429 {
        ExternalAttemptErrorKind::Transient
    } else if status.as_u16() == 401 || status.as_u16() == 403 {
        ExternalAttemptErrorKind::Unavailable
    } else {
        ExternalAttemptErrorKind::Fatal
    };
    ExternalAttemptError {
        kind,
        message: message_prefix,
    }
}

fn invoke_native_openai_once(
    system: &AgentSystemSpec,
    request_json: &str,
) -> Result<(String, String), ExternalAttemptError> {
    let api_key = resolve_api_key(system).ok_or_else(|| ExternalAttemptError {
        kind: ExternalAttemptErrorKind::Unavailable,
        message: format!("{OPENAI_API_KEY_ENV} is not set"),
    })?;
    let model = resolve_model(system, OPENAI_DEFAULT_MODEL);
    let base_url = resolve_base_url(system, OPENAI_DEFAULT_BASE_URL);
    let endpoint = format!("{base_url}/responses");
    let system_prompt = "You are a GENtle agent bridge.\nReturn STRICT JSON only with this schema:\n{\"schema\":\"gentle.agent_response.v1\",\"assistant_message\":\"string\",\"questions\":[\"string\"],\"suggested_commands\":[{\"title\":\"string\",\"rationale\":\"string\",\"command\":\"string\",\"execution\":\"chat|ask|auto\"}]}\nUse only keys from the schema. Extensions may use x_ prefix. Do not include markdown fences.";
    let payload = json!({
        "model": model,
        "input": [
            {
                "role": "system",
                "content": [
                    { "type": "input_text", "text": system_prompt }
                ]
            },
            {
                "role": "user",
                "content": [
                    {
                        "type": "input_text",
                        "text": format!("GENtle agent request JSON:\n{request_json}")
                    }
                ]
            }
        ]
    });
    let client = reqwest::blocking::Client::builder()
        .timeout(Duration::from_secs(60))
        .build()
        .map_err(|e| ExternalAttemptError {
            kind: ExternalAttemptErrorKind::Fatal,
            message: format!("could not build OpenAI client: {e}"),
        })?;
    let response = client
        .post(&endpoint)
        .bearer_auth(api_key)
        .header("Content-Type", "application/json")
        .json(&payload)
        .send()
        .map_err(|e| ExternalAttemptError {
            kind: if e.is_timeout() || e.is_connect() || e.is_request() {
                ExternalAttemptErrorKind::Transient
            } else {
                ExternalAttemptErrorKind::Fatal
            },
            message: format!("OpenAI request failed: {e}"),
        })?;
    let status = response.status();
    let body = response.text().map_err(|e| ExternalAttemptError {
        kind: ExternalAttemptErrorKind::Fatal,
        message: format!("could not read OpenAI response body: {e}"),
    })?;
    if !status.is_success() {
        return Err(classify_openai_http_error(status, &body));
    }
    let response_json = serde_json::from_str::<Value>(&body).map_err(|e| ExternalAttemptError {
        kind: ExternalAttemptErrorKind::Fatal,
        message: format!("OpenAI API returned invalid JSON: {e}"),
    })?;
    let text = extract_openai_output_text(&response_json).ok_or_else(|| ExternalAttemptError {
        kind: ExternalAttemptErrorKind::Fatal,
        message: "OpenAI API response did not contain output_text".to_string(),
    })?;
    Ok((text, body))
}

fn invoke_native_openai_compat_once(
    system: &AgentSystemSpec,
    request_json: &str,
) -> Result<(String, String), ExternalAttemptError> {
    let api_key = resolve_api_key(system);
    let model = resolve_model(system, OPENAI_COMPAT_DEFAULT_MODEL);
    let base_url = resolve_base_url(system, OPENAI_COMPAT_DEFAULT_BASE_URL);
    let endpoints = openai_compat_endpoint_candidates_for_system(system);
    let system_prompt = "You are a GENtle agent bridge.\nReturn STRICT JSON only with this schema:\n{\"schema\":\"gentle.agent_response.v1\",\"assistant_message\":\"string\",\"questions\":[\"string\"],\"suggested_commands\":[{\"title\":\"string\",\"rationale\":\"string\",\"command\":\"string\",\"execution\":\"chat|ask|auto\"}]}\nUse only keys from the schema. Extensions may use x_ prefix. Do not include markdown fences.";
    let payload = json!({
        "model": model,
        "messages": [
            {
                "role": "system",
                "content": system_prompt
            },
            {
                "role": "user",
                "content": format!("GENtle agent request JSON:\n{request_json}")
            }
        ],
        "temperature": 0.2
    });
    let client = reqwest::blocking::Client::builder()
        .timeout(Duration::from_secs(60))
        .build()
        .map_err(|e| ExternalAttemptError {
            kind: ExternalAttemptErrorKind::Fatal,
            message: format!("could not build OpenAI-compatible client: {e}"),
        })?;
    let mut first_path_error: Option<ExternalAttemptError> = None;
    for (idx, endpoint) in endpoints.iter().enumerate() {
        let mut request = client
            .post(endpoint)
            .header("Content-Type", "application/json");
        if let Some(key) = api_key.as_deref() {
            request = request.bearer_auth(key);
        }
        let response = match request.json(&payload).send() {
            Ok(response) => response,
            Err(e) => {
                let kind = if e.is_timeout() || e.is_connect() || e.is_request() {
                    ExternalAttemptErrorKind::Transient
                } else {
                    ExternalAttemptErrorKind::Fatal
                };
                let error = ExternalAttemptError {
                    kind,
                    message: format!("OpenAI-compatible request failed ({endpoint}): {e}"),
                };
                let can_try_fallback = matches!(kind, ExternalAttemptErrorKind::Transient)
                    && idx + 1 < endpoints.len();
                if can_try_fallback {
                    if first_path_error.is_none() {
                        first_path_error = Some(error);
                    }
                    continue;
                }
                return Err(error);
            }
        };
        let status = response.status();
        let body = response.text().map_err(|e| ExternalAttemptError {
            kind: ExternalAttemptErrorKind::Fatal,
            message: format!("could not read OpenAI-compatible response body: {e}"),
        })?;
        if !status.is_success() {
            let classified = classify_openai_http_error(status, &body);
            let can_try_fallback = (status.as_u16() == 404 || status.as_u16() == 405)
                && idx + 1 < endpoints.len()
                && !matches!(classified.kind, ExternalAttemptErrorKind::Unavailable);
            if can_try_fallback {
                if first_path_error.is_none() {
                    first_path_error = Some(classified);
                }
                continue;
            }
            return Err(classified);
        }
        let response_json =
            serde_json::from_str::<Value>(&body).map_err(|e| ExternalAttemptError {
                kind: ExternalAttemptErrorKind::Fatal,
                message: format!("OpenAI-compatible API returned invalid JSON: {e}"),
            })?;
        let text = extract_openai_chat_completions_text(&response_json).ok_or_else(|| {
            ExternalAttemptError {
                kind: ExternalAttemptErrorKind::Fatal,
                message: "OpenAI-compatible response did not contain choices[0].message.content"
                    .to_string(),
            }
        })?;
        return Ok((text, body));
    }
    Err(first_path_error.unwrap_or_else(|| ExternalAttemptError {
        kind: ExternalAttemptErrorKind::Fatal,
        message: format!(
            "OpenAI-compatible call failed before receiving a response body (base_url={base_url})"
        ),
    }))
}

fn builtin_echo_response(prompt: &str) -> AgentResponse {
    let mut suggested: Vec<AgentSuggestedCommand> = vec![];
    for line in prompt.lines() {
        let trimmed = line.trim();
        if let Some(command) = trimmed.strip_prefix("auto:") {
            let command = command.trim();
            if !command.is_empty() {
                suggested.push(AgentSuggestedCommand {
                    title: Some("Auto suggestion (demo)".to_string()),
                    rationale: Some("Extracted from 'auto:' line in prompt".to_string()),
                    command: command.to_string(),
                    execution: AgentExecutionIntent::Auto,
                });
            }
            continue;
        }
        if let Some(command) = trimmed.strip_prefix("ask:") {
            let command = command.trim();
            if !command.is_empty() {
                suggested.push(AgentSuggestedCommand {
                    title: Some("Confirm suggestion (demo)".to_string()),
                    rationale: Some("Extracted from 'ask:' line in prompt".to_string()),
                    command: command.to_string(),
                    execution: AgentExecutionIntent::Ask,
                });
            }
            continue;
        }
    }

    AgentResponse {
        schema: AGENT_RESPONSE_SCHEMA.to_string(),
        assistant_message: format!(
            "Builtin echo agent received your prompt. Configure an external agent system in '{}' for real project assistance.",
            DEFAULT_AGENT_SYSTEM_CATALOG_PATH
        ),
        questions: vec![],
        suggested_commands: suggested,
    }
}

pub fn load_agent_system_catalog(
    catalog_path: Option<&str>,
) -> Result<(String, AgentSystemCatalog), String> {
    let resolved_path = AgentSystemCatalog::effective_catalog_path(catalog_path);
    let catalog = AgentSystemCatalog::from_json_file(&resolved_path)?;
    Ok((resolved_path, catalog))
}

pub fn invoke_agent_support_with_env_overrides(
    catalog_path: Option<&str>,
    system_id: &str,
    prompt: &str,
    state_summary: Option<&EngineStateSummary>,
    env_overrides: Option<&HashMap<String, String>>,
) -> Result<AgentInvocationOutcome, String> {
    if prompt.trim().is_empty() {
        return Err(agent_err(
            AgentBridgeErrorCode::InvalidInput,
            "agent prompt cannot be empty",
        ));
    }
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
    let availability = agent_system_availability(&system);
    if !availability.available {
        return Err(agent_err(
            AgentBridgeErrorCode::AdapterUnavailable,
            availability
                .reason
                .unwrap_or_else(|| "agent system is not available".to_string()),
        ));
    }
    let (_payload, request_value, request_json) =
        build_agent_request(&system.id, prompt, state_summary)?;
    let start = std::time::Instant::now();

    let outcome = match system.transport {
        AgentSystemTransport::BuiltinEcho => {
            let response = builtin_echo_response(prompt);
            AgentInvocationOutcome {
                catalog_path: resolved_catalog_path,
                system_id: system.id,
                system_label: system.label,
                transport: system.transport.as_str().to_string(),
                command: vec![],
                request: request_value,
                response,
                raw_stdout: String::new(),
                raw_stderr: String::new(),
                exit_code: Some(0),
                elapsed_ms: start.elapsed().as_millis(),
            }
        }
        AgentSystemTransport::ExternalJsonStdio => {
            if system.command.is_empty() {
                return Err(agent_err(
                    AgentBridgeErrorCode::CatalogValidation,
                    format!("agent system '{}' has no command configured", system.id),
                ));
            }
            let state_summary_json = request_value
                .get("state_summary")
                .cloned()
                .unwrap_or(Value::Null)
                .to_string();
            let rendered_command = system
                .command
                .iter()
                .map(|token| {
                    replace_placeholders(token, prompt, &request_json, &state_summary_json)
                })
                .collect::<Vec<_>>();
            let mut output: Option<Output> = None;
            let mut last_transient_message: Option<String> = None;
            for attempt in 1..=AGENT_INVOKE_RETRY_ATTEMPTS {
                match invoke_external_json_stdio_once(&system, &rendered_command, &request_json) {
                    Ok(result) => {
                        output = Some(result);
                        break;
                    }
                    Err(error) => match error.kind {
                        ExternalAttemptErrorKind::Unavailable => {
                            return Err(agent_err(
                                AgentBridgeErrorCode::AdapterUnavailable,
                                error.message,
                            ));
                        }
                        ExternalAttemptErrorKind::Fatal => {
                            return Err(agent_err(
                                AgentBridgeErrorCode::AdapterFailed,
                                error.message,
                            ));
                        }
                        ExternalAttemptErrorKind::Transient => {
                            last_transient_message = Some(error.message.clone());
                            if attempt >= AGENT_INVOKE_RETRY_ATTEMPTS {
                                return Err(agent_err(
                                    AgentBridgeErrorCode::AdapterTransient,
                                    format!(
                                        "agent command remained transiently unavailable after {} attempts: {}",
                                        AGENT_INVOKE_RETRY_ATTEMPTS, error.message
                                    ),
                                ));
                            }
                            thread::sleep(retry_backoff_duration(attempt));
                        }
                    },
                }
            }
            if output.is_none() {
                return Err(agent_err(
                    AgentBridgeErrorCode::AdapterTransient,
                    format!(
                        "agent command did not complete after {} attempts{}",
                        AGENT_INVOKE_RETRY_ATTEMPTS,
                        last_transient_message
                            .as_ref()
                            .map(|value| format!(": {value}"))
                            .unwrap_or_default()
                    ),
                ));
            }
            let output = output.expect("checked is_some above");
            let stdout = String::from_utf8_lossy(&output.stdout).to_string();
            let stderr = String::from_utf8_lossy(&output.stderr).to_string();
            let response = parse_agent_response(&stdout)?;
            AgentInvocationOutcome {
                catalog_path: resolved_catalog_path,
                system_id: system.id,
                system_label: system.label,
                transport: system.transport.as_str().to_string(),
                command: rendered_command,
                request: request_value,
                response,
                raw_stdout: stdout,
                raw_stderr: stderr,
                exit_code: output.status.code(),
                elapsed_ms: start.elapsed().as_millis(),
            }
        }
        AgentSystemTransport::NativeOpenai => {
            let mut stdout: Option<String> = None;
            let mut raw_body: Option<String> = None;
            let mut last_transient_message: Option<String> = None;
            for attempt in 1..=AGENT_INVOKE_RETRY_ATTEMPTS {
                match invoke_native_openai_once(&system, &request_json) {
                    Ok((result_text, raw_json_body)) => {
                        stdout = Some(result_text);
                        raw_body = Some(raw_json_body);
                        break;
                    }
                    Err(error) => match error.kind {
                        ExternalAttemptErrorKind::Unavailable => {
                            return Err(agent_err(
                                AgentBridgeErrorCode::AdapterUnavailable,
                                error.message,
                            ));
                        }
                        ExternalAttemptErrorKind::Fatal => {
                            return Err(agent_err(
                                AgentBridgeErrorCode::AdapterFailed,
                                error.message,
                            ));
                        }
                        ExternalAttemptErrorKind::Transient => {
                            last_transient_message = Some(error.message.clone());
                            if attempt >= AGENT_INVOKE_RETRY_ATTEMPTS {
                                return Err(agent_err(
                                    AgentBridgeErrorCode::AdapterTransient,
                                    format!(
                                        "native OpenAI call remained transiently unavailable after {} attempts: {}",
                                        AGENT_INVOKE_RETRY_ATTEMPTS, error.message
                                    ),
                                ));
                            }
                            thread::sleep(retry_backoff_duration(attempt));
                        }
                    },
                }
            }
            if stdout.is_none() {
                return Err(agent_err(
                    AgentBridgeErrorCode::AdapterTransient,
                    format!(
                        "native OpenAI call did not complete after {} attempts{}",
                        AGENT_INVOKE_RETRY_ATTEMPTS,
                        last_transient_message
                            .as_ref()
                            .map(|value| format!(": {value}"))
                            .unwrap_or_default()
                    ),
                ));
            }
            let stdout = stdout.expect("checked is_some above");
            let response = parse_agent_response(&stdout)?;
            AgentInvocationOutcome {
                catalog_path: resolved_catalog_path,
                system_id: system.id,
                system_label: system.label,
                transport: system.transport.as_str().to_string(),
                command: vec![],
                request: request_value,
                response,
                raw_stdout: stdout,
                raw_stderr: raw_body.unwrap_or_default(),
                exit_code: Some(0),
                elapsed_ms: start.elapsed().as_millis(),
            }
        }
        AgentSystemTransport::NativeOpenaiCompat => {
            let mut stdout: Option<String> = None;
            let mut raw_body: Option<String> = None;
            let mut last_transient_message: Option<String> = None;
            for attempt in 1..=AGENT_INVOKE_RETRY_ATTEMPTS {
                match invoke_native_openai_compat_once(&system, &request_json) {
                    Ok((result_text, raw_json_body)) => {
                        stdout = Some(result_text);
                        raw_body = Some(raw_json_body);
                        break;
                    }
                    Err(error) => match error.kind {
                        ExternalAttemptErrorKind::Unavailable => {
                            return Err(agent_err(
                                AgentBridgeErrorCode::AdapterUnavailable,
                                error.message,
                            ));
                        }
                        ExternalAttemptErrorKind::Fatal => {
                            return Err(agent_err(
                                AgentBridgeErrorCode::AdapterFailed,
                                error.message,
                            ));
                        }
                        ExternalAttemptErrorKind::Transient => {
                            last_transient_message = Some(error.message.clone());
                            if attempt >= AGENT_INVOKE_RETRY_ATTEMPTS {
                                return Err(agent_err(
                                    AgentBridgeErrorCode::AdapterTransient,
                                    format!(
                                        "native OpenAI-compatible call remained transiently unavailable after {} attempts: {}",
                                        AGENT_INVOKE_RETRY_ATTEMPTS, error.message
                                    ),
                                ));
                            }
                            thread::sleep(retry_backoff_duration(attempt));
                        }
                    },
                }
            }
            if stdout.is_none() {
                return Err(agent_err(
                    AgentBridgeErrorCode::AdapterTransient,
                    format!(
                        "native OpenAI-compatible call did not complete after {} attempts{}",
                        AGENT_INVOKE_RETRY_ATTEMPTS,
                        last_transient_message
                            .as_ref()
                            .map(|value| format!(": {value}"))
                            .unwrap_or_default()
                    ),
                ));
            }
            let stdout = stdout.expect("checked is_some above");
            let response = parse_agent_response(&stdout)?;
            AgentInvocationOutcome {
                catalog_path: resolved_catalog_path,
                system_id: system.id,
                system_label: system.label,
                transport: system.transport.as_str().to_string(),
                command: vec![],
                request: request_value,
                response,
                raw_stdout: stdout,
                raw_stderr: raw_body.unwrap_or_default(),
                exit_code: Some(0),
                elapsed_ms: start.elapsed().as_millis(),
            }
        }
    };
    Ok(outcome)
}

pub fn invoke_agent_support(
    catalog_path: Option<&str>,
    system_id: &str,
    prompt: &str,
    state_summary: Option<&EngineStateSummary>,
) -> Result<AgentInvocationOutcome, String> {
    invoke_agent_support_with_env_overrides(catalog_path, system_id, prompt, state_summary, None)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builtin_echo_extracts_demo_suggestions() {
        let response = builtin_echo_response("auto: state-summary\nask: capabilities");
        assert_eq!(response.suggested_commands.len(), 2);
        assert_eq!(
            response.suggested_commands[0].execution,
            AgentExecutionIntent::Auto
        );
        assert_eq!(
            response.suggested_commands[1].execution,
            AgentExecutionIntent::Ask
        );
    }

    #[test]
    fn parse_agent_response_rejects_plain_text() {
        let err = parse_agent_response("hello world").expect_err("plain text should fail");
        assert!(err.starts_with("AGENT_RESPONSE_PARSE:"));
    }

    #[test]
    fn parse_agent_response_parses_command_objects() {
        let response = parse_agent_response(
            r#"{
  "schema": "gentle.agent_response.v1",
  "assistant_message": "ready",
  "suggested_commands": [
    {"command":"state-summary","execution":"auto"},
    {"command":"capabilities","execution":"ask"}
  ]
}"#,
        )
        .expect("json response parse");
        assert_eq!(response.suggested_commands.len(), 2);
        assert_eq!(
            response.suggested_commands[0].execution,
            AgentExecutionIntent::Auto
        );
    }

    #[test]
    fn parse_agent_response_rejects_unknown_critical_fields() {
        let err = parse_agent_response(
            r#"{
  "schema": "gentle.agent_response.v1",
  "assistant_message": "ready",
  "commands": []
}"#,
        )
        .expect_err("unknown canonical-critical field should fail");
        assert!(err.starts_with("AGENT_RESPONSE_VALIDATION:"));
    }

    #[test]
    fn parse_agent_response_rejects_future_schema_major() {
        let err = parse_agent_response(
            r#"{
  "schema": "gentle.agent_response.v2",
  "assistant_message": "ready",
  "questions": [],
  "suggested_commands": []
}"#,
        )
        .expect_err("future schema major should fail");
        assert!(err.starts_with("AGENT_SCHEMA_UNSUPPORTED:"));
    }

    #[test]
    fn availability_native_openai_uses_system_env_override() {
        let mut system = AgentSystemSpec {
            id: "openai_native".to_string(),
            label: "OpenAI Native".to_string(),
            description: None,
            transport: AgentSystemTransport::NativeOpenai,
            command: vec![],
            model: Some("gpt-5".to_string()),
            base_url: None,
            env: HashMap::new(),
            working_dir: None,
        };
        let unavailable = agent_system_availability(&system);
        assert!(!unavailable.available);
        system
            .env
            .insert(OPENAI_API_KEY_ENV.to_string(), "sk-test".to_string());
        let available = agent_system_availability(&system);
        assert!(available.available);
    }

    #[test]
    fn availability_external_stdio_requires_command() {
        let system = AgentSystemSpec {
            id: "stdio".to_string(),
            label: "stdio".to_string(),
            description: None,
            transport: AgentSystemTransport::ExternalJsonStdio,
            command: vec![],
            model: None,
            base_url: None,
            env: HashMap::new(),
            working_dir: None,
        };
        let availability = agent_system_availability(&system);
        assert!(!availability.available);
    }

    #[test]
    fn availability_native_openai_compat_allows_missing_key() {
        let system = AgentSystemSpec {
            id: "local-compat".to_string(),
            label: "local-compat".to_string(),
            description: None,
            transport: AgentSystemTransport::NativeOpenaiCompat,
            command: vec![],
            model: Some("llama3.1".to_string()),
            base_url: Some("http://127.0.0.1:11434/v1".to_string()),
            env: HashMap::new(),
            working_dir: None,
        };
        let availability = agent_system_availability(&system);
        assert!(availability.available);
    }

    #[test]
    fn availability_native_openai_compat_uses_env_base_url_override() {
        let mut system = AgentSystemSpec {
            id: "local-compat".to_string(),
            label: "local-compat".to_string(),
            description: None,
            transport: AgentSystemTransport::NativeOpenaiCompat,
            command: vec![],
            model: Some("llama3.1".to_string()),
            base_url: Some("http://127.0.0.1:11434/v1".to_string()),
            env: HashMap::new(),
            working_dir: None,
        };
        system.env.insert(
            AGENT_BASE_URL_ENV.to_string(),
            "http://localhost:11964".to_string(),
        );
        system
            .env
            .insert(AGENT_MODEL_ENV.to_string(), "deepseek-r1:8b".to_string());
        let availability = agent_system_availability(&system);
        let reason = availability.reason.unwrap_or_default();
        assert!(availability.available);
        assert!(reason.contains("localhost:11964"));
        assert!(reason.contains("deepseek-r1:8b"));
        assert!(reason.contains("/v1/chat/completions"));
    }

    #[test]
    fn openai_compat_endpoint_candidates_include_v1_fallback_for_root_base() {
        let endpoints = openai_compat_endpoint_candidates("http://localhost:11964");
        assert_eq!(
            endpoints,
            vec![
                "http://localhost:11964/chat/completions".to_string(),
                "http://localhost:11964/v1/chat/completions".to_string()
            ]
        );
    }

    #[test]
    fn openai_compat_endpoint_candidates_for_system_include_local_fallbacks_without_override() {
        let system = AgentSystemSpec {
            id: "local-compat".to_string(),
            label: "local-compat".to_string(),
            description: None,
            transport: AgentSystemTransport::NativeOpenaiCompat,
            command: vec![],
            model: Some("llama3.1".to_string()),
            base_url: Some("http://127.0.0.1:10000/v1".to_string()),
            env: HashMap::new(),
            working_dir: None,
        };
        let endpoints = openai_compat_endpoint_candidates_for_system(&system);
        assert!(
            endpoints
                .iter()
                .any(|value| value == "http://127.0.0.1:10000/v1/chat/completions")
        );
        assert!(
            endpoints
                .iter()
                .any(|value| value == "http://localhost:11964/chat/completions")
        );
    }

    #[test]
    fn openai_compat_endpoint_candidates_for_system_respects_explicit_override() {
        let mut system = AgentSystemSpec {
            id: "local-compat".to_string(),
            label: "local-compat".to_string(),
            description: None,
            transport: AgentSystemTransport::NativeOpenaiCompat,
            command: vec![],
            model: Some("llama3.1".to_string()),
            base_url: Some("http://127.0.0.1:10000/v1".to_string()),
            env: HashMap::new(),
            working_dir: None,
        };
        system.env.insert(
            AGENT_BASE_URL_ENV.to_string(),
            "http://localhost:11964".to_string(),
        );
        let endpoints = openai_compat_endpoint_candidates_for_system(&system);
        assert_eq!(
            endpoints,
            vec![
                "http://localhost:11964/chat/completions".to_string(),
                "http://localhost:11964/v1/chat/completions".to_string()
            ]
        );
    }

    #[test]
    fn extract_models_from_openai_models_payload_reads_standard_data_shape() {
        let value = serde_json::json!({
            "object": "list",
            "data": [
                {"id":"deepseek-r1:8b"},
                {"id":"qwen3:0.6b"}
            ]
        });
        let models = extract_models_from_openai_models_payload(&value);
        assert_eq!(
            models,
            vec!["deepseek-r1:8b".to_string(), "qwen3:0.6b".to_string()]
        );
    }

    #[test]
    fn classify_openai_insufficient_quota_is_non_transient_with_hint() {
        let body = r#"{
  "error": {
    "message": "You exceeded your current quota",
    "type": "insufficient_quota",
    "code": "insufficient_quota"
  }
}"#;
        let err = classify_openai_http_error(reqwest::StatusCode::TOO_MANY_REQUESTS, body);
        assert_eq!(err.kind, ExternalAttemptErrorKind::Unavailable);
        assert!(err.message.contains("OpenAI API error (status=429"));
        assert!(err.message.contains("insufficient quota"));
    }
}
