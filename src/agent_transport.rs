//! Shared transport/runtime helpers for GENtle agent systems.
//!
//! This module exposes the machine-facing transport metadata needed by both the
//! local assistant UX and external orchestrators. It intentionally keeps the
//! public surface read-only and deterministic: catalog loading, availability,
//! preflight summaries, and OpenAI-compatible model discovery.

use crate::agent_bridge::{
    AGENT_BASE_URL_ENV, AGENT_CONNECT_TIMEOUT_SECS_ENV, AGENT_MAX_RESPONSE_BYTES_ENV,
    AGENT_MAX_RETRIES_ENV, AGENT_MODEL_ENV, AGENT_READ_TIMEOUT_SECS_ENV, AGENT_TIMEOUT_SECS_ENV,
    ANTHROPIC_API_KEY_ENV, OPENAI_API_KEY_ENV, extract_anthropic_error_code,
    extract_models_from_openai_models_payload, extract_openai_error_code, redact_sensitive_text,
};
use serde::{Deserialize, Serialize};
use std::time::Duration;

pub use crate::agent_bridge::{
    AgentInvocationRuntime, AgentSystemAvailability, AgentSystemCatalog, AgentSystemSpec,
    AgentSystemTransport, DEFAULT_AGENT_SYSTEM_CATALOG_PATH, agent_system_availability,
    discover_anthropic_models, discover_openai_models, load_agent_system_catalog,
};

const AGENT_REQUEST_TIMEOUT_SECS_DEFAULT: u64 = 180;
const AGENT_CONNECT_TIMEOUT_SECS_DEFAULT: u64 = 10;
const AGENT_MAX_RETRIES_DEFAULT: usize = 2;
const AGENT_MAX_RESPONSE_BYTES_DEFAULT: usize = 1_048_576;
const OPENAI_DEFAULT_MODEL: &str = "gpt-5";
const OPENAI_DEFAULT_BASE_URL: &str = "https://api.openai.com/v1";
const ANTHROPIC_DEFAULT_MODEL: &str = "claude-sonnet-4-6";
const ANTHROPIC_DEFAULT_BASE_URL: &str = "https://api.anthropic.com/v1";
const ANTHROPIC_API_VERSION: &str = "2023-06-01";
const OPENAI_COMPAT_DEFAULT_MODEL: &str = "unspecified";
const OPENAI_COMPAT_DEFAULT_BASE_URL: &str = "http://127.0.0.1:11434/v1";

/// Stable read-only preflight result for one configured agent system.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct AgentSystemPreflight {
    pub schema: String,
    pub catalog_path: String,
    pub system_id: String,
    pub system_label: String,
    pub transport: String,
    pub available: bool,
    pub availability_reason: Option<String>,
    pub base_url: Option<String>,
    pub model: Option<String>,
    pub endpoint_candidates: Vec<String>,
    pub model_endpoint_candidates: Vec<String>,
    pub timeout_secs: u64,
    pub connect_timeout_secs: u64,
    pub read_timeout_secs: u64,
    pub max_retries: usize,
    pub max_response_bytes: usize,
    pub warnings: Vec<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub live_probe: Option<AgentSystemLiveProbe>,
}

/// Live, non-generating model-discovery probe for one configured agent system.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct AgentSystemLiveProbe {
    pub enabled: bool,
    pub attempted_endpoints: Vec<String>,
    pub selected_endpoint: Option<String>,
    pub reachable: bool,
    pub auth_ok: bool,
    pub model_list_ok: bool,
    pub selected_model_seen: bool,
    pub status_class: AgentLiveProbeStatusClass,
    pub message: String,
    pub provider_error_code: Option<String>,
}

/// Coarse, user-facing live-probe classification.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum AgentLiveProbeStatusClass {
    Ok,
    MissingKey,
    AuthFailed,
    QuotaOrBilling,
    ModelMissing,
    EndpointUnreachable,
    UnsupportedTransport,
    #[default]
    ProviderError,
}

impl AgentLiveProbeStatusClass {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Ok => "ok",
            Self::MissingKey => "missing_key",
            Self::AuthFailed => "auth_failed",
            Self::QuotaOrBilling => "quota_or_billing",
            Self::ModelMissing => "model_missing",
            Self::EndpointUnreachable => "endpoint_unreachable",
            Self::UnsupportedTransport => "unsupported_transport",
            Self::ProviderError => "provider_error",
        }
    }
}

fn parse_positive_u64(raw: &str) -> Option<u64> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return None;
    }
    trimmed.parse::<u64>().ok().filter(|value| *value > 0)
}

fn parse_positive_usize(raw: &str) -> Option<usize> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return None;
    }
    trimmed.parse::<usize>().ok().filter(|value| *value > 0)
}

fn parse_nonnegative_usize(raw: &str) -> Option<usize> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return None;
    }
    trimmed.parse::<usize>().ok()
}

fn resolve_u64_override(system: &AgentSystemSpec, key: &str) -> Option<u64> {
    system
        .env
        .get(key)
        .and_then(|value| parse_positive_u64(value))
        .or_else(|| {
            std::env::var(key)
                .ok()
                .and_then(|value| parse_positive_u64(&value))
        })
}

fn resolve_usize_override(system: &AgentSystemSpec, key: &str) -> Option<usize> {
    system
        .env
        .get(key)
        .and_then(|value| parse_positive_usize(value))
        .or_else(|| {
            std::env::var(key)
                .ok()
                .and_then(|value| parse_positive_usize(&value))
        })
}

fn resolve_nonnegative_usize_override(system: &AgentSystemSpec, key: &str) -> Option<usize> {
    system
        .env
        .get(key)
        .and_then(|value| parse_nonnegative_usize(value))
        .or_else(|| {
            std::env::var(key)
                .ok()
                .and_then(|value| parse_nonnegative_usize(&value))
        })
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
                .map(ToString::to_string)
                .filter(|value| !value.is_empty())
        })
        .unwrap_or_else(|| default.to_string())
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

fn resolve_env_key(system: &AgentSystemSpec, key: &str) -> Option<String> {
    system
        .env
        .get(key)
        .map(|value| value.trim().to_string())
        .filter(|value| !value.is_empty())
        .or_else(|| {
            std::env::var(key)
                .ok()
                .map(|value| value.trim().to_string())
                .filter(|value| !value.is_empty())
        })
}

fn resolve_openai_api_key(system: &AgentSystemSpec) -> Option<String> {
    resolve_env_key(system, OPENAI_API_KEY_ENV)
}

fn resolve_anthropic_api_key(system: &AgentSystemSpec) -> Option<String> {
    resolve_env_key(system, ANTHROPIC_API_KEY_ENV)
}

fn openai_compat_endpoint_candidates(base_url: &str) -> Vec<String> {
    let normalized =
        normalize_base_url(base_url).unwrap_or_else(|| base_url.trim_end_matches('/').to_string());
    let mut endpoints = vec![format!("{normalized}/chat/completions")];
    if !normalized.ends_with("/v1") {
        let v1 = format!("{normalized}/v1/chat/completions");
        if !endpoints.contains(&v1) {
            endpoints.push(v1);
        }
    }
    endpoints
}

fn model_endpoint_candidates(base_url: &str) -> Vec<String> {
    let normalized =
        normalize_base_url(base_url).unwrap_or_else(|| base_url.trim_end_matches('/').to_string());
    let mut endpoints = vec![format!("{normalized}/models")];
    if !normalized.ends_with("/v1") {
        let v1 = format!("{normalized}/v1/models");
        if !endpoints.contains(&v1) {
            endpoints.push(v1);
        }
    }
    endpoints
}

fn agent_runtime(system: &AgentSystemSpec) -> AgentInvocationRuntime {
    let timeout_secs = resolve_u64_override(system, AGENT_TIMEOUT_SECS_ENV)
        .unwrap_or(AGENT_REQUEST_TIMEOUT_SECS_DEFAULT);
    let connect_timeout_secs = resolve_u64_override(system, AGENT_CONNECT_TIMEOUT_SECS_ENV)
        .unwrap_or(AGENT_CONNECT_TIMEOUT_SECS_DEFAULT);
    let read_timeout_secs =
        resolve_u64_override(system, AGENT_READ_TIMEOUT_SECS_ENV).unwrap_or(timeout_secs);
    let max_retries = resolve_nonnegative_usize_override(system, AGENT_MAX_RETRIES_ENV)
        .unwrap_or(AGENT_MAX_RETRIES_DEFAULT);
    let max_response_bytes = resolve_usize_override(system, AGENT_MAX_RESPONSE_BYTES_ENV)
        .unwrap_or(AGENT_MAX_RESPONSE_BYTES_DEFAULT);
    AgentInvocationRuntime {
        timeout_secs,
        connect_timeout_secs: Some(connect_timeout_secs),
        read_timeout_secs: Some(read_timeout_secs),
        max_retries,
        max_response_bytes,
        endpoint_candidates: vec![],
        attempted_endpoints: vec![],
        selected_endpoint: None,
    }
}

fn is_unspecified_model_name(raw: &str) -> bool {
    let trimmed = raw.trim();
    trimmed.is_empty() || trimmed.eq_ignore_ascii_case(OPENAI_COMPAT_DEFAULT_MODEL)
}

fn clipped_message(raw: &str) -> String {
    let trimmed = raw.trim();
    let mut out = trimmed.chars().take(700).collect::<String>();
    if trimmed.chars().count() > 700 {
        out.push_str(" ...[truncated]");
    }
    out
}

fn live_probe_with_status(
    status_class: AgentLiveProbeStatusClass,
    message: impl Into<String>,
) -> AgentSystemLiveProbe {
    AgentSystemLiveProbe {
        enabled: true,
        status_class,
        message: redact_sensitive_text(&message.into()),
        ..Default::default()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ModelListAuth {
    OpenaiBearer,
    Anthropic,
}

fn read_live_probe_body(
    mut response: reqwest::blocking::Response,
    max_response_bytes: usize,
) -> Result<String, String> {
    use std::io::Read;
    let mut body = Vec::new();
    response
        .by_ref()
        .take(max_response_bytes as u64 + 1)
        .read_to_end(&mut body)
        .map_err(|err| format!("could not read model-list response body: {err}"))?;
    if body.len() > max_response_bytes {
        return Err(format!(
            "model-list response body exceeded {} bytes (set {} to increase)",
            max_response_bytes, AGENT_MAX_RESPONSE_BYTES_ENV
        ));
    }
    String::from_utf8(body).map_err(|err| format!("model-list response body is not UTF-8: {err}"))
}

fn classify_model_list_http_error(
    status: reqwest::StatusCode,
    endpoint: &str,
    body: &str,
    provider: &str,
) -> (AgentLiveProbeStatusClass, bool, String, Option<String>) {
    let provider_error_code = if provider == "Anthropic" {
        extract_anthropic_error_code(body)
    } else {
        extract_openai_error_code(body)
    };
    let body_preview = clipped_message(body);
    let lower_body = body_preview.to_ascii_lowercase();
    let lower_code = provider_error_code
        .as_deref()
        .unwrap_or_default()
        .to_ascii_lowercase();
    let status_class = if status.as_u16() == 401 || status.as_u16() == 403 {
        AgentLiveProbeStatusClass::AuthFailed
    } else if status.as_u16() == 429
        && (lower_code.contains("insufficient_quota")
            || lower_code.contains("quota")
            || lower_code.contains("billing")
            || lower_body.contains("insufficient_quota")
            || lower_body.contains("quota")
            || lower_body.contains("billing"))
    {
        AgentLiveProbeStatusClass::QuotaOrBilling
    } else {
        AgentLiveProbeStatusClass::ProviderError
    };
    let auth_ok = !matches!(status_class, AgentLiveProbeStatusClass::AuthFailed);
    let message = match status_class {
        AgentLiveProbeStatusClass::AuthFailed => format!(
            "{provider} model-list endpoint rejected authentication at {endpoint} (status={status}). Check the API key/token for this provider."
        ),
        AgentLiveProbeStatusClass::QuotaOrBilling => format!(
            "{provider} reported quota or billing trouble during model-list probe at {endpoint} (status={status}). The probe did not send a generation request. Response: {body_preview}"
        ),
        _ => format!(
            "{provider} model-list probe failed at {endpoint} (status={status}). Response: {body_preview}"
        ),
    };
    (
        status_class,
        auth_ok,
        redact_sensitive_text(&message),
        provider_error_code,
    )
}

fn build_model_list_live_probe(
    base_url: &str,
    api_key: Option<&str>,
    selected_model: Option<&str>,
    runtime: &AgentInvocationRuntime,
    auth: ModelListAuth,
) -> AgentSystemLiveProbe {
    let endpoints = model_endpoint_candidates(base_url);
    let client = match reqwest::blocking::Client::builder()
        .connect_timeout(Duration::from_secs(
            runtime.connect_timeout_secs.unwrap_or_default().max(1),
        ))
        .timeout(Duration::from_secs(
            runtime
                .read_timeout_secs
                .unwrap_or(runtime.timeout_secs)
                .max(1),
        ))
        .build()
    {
        Ok(client) => client,
        Err(err) => {
            return live_probe_with_status(
                AgentLiveProbeStatusClass::ProviderError,
                format!("Could not build model-list probe HTTP client: {err}"),
            );
        }
    };

    let mut probe = AgentSystemLiveProbe {
        enabled: true,
        status_class: AgentLiveProbeStatusClass::EndpointUnreachable,
        message: "Model-list probe did not reach any endpoint".to_string(),
        ..Default::default()
    };
    let mut first_request_error: Option<String> = None;

    for (idx, endpoint) in endpoints.iter().enumerate() {
        probe.attempted_endpoints.push(endpoint.clone());
        let mut request = client.get(endpoint);
        if let Some(key) = api_key.map(str::trim).filter(|value| !value.is_empty()) {
            request = match auth {
                ModelListAuth::OpenaiBearer => request.bearer_auth(key),
                ModelListAuth::Anthropic => request
                    .header("x-api-key", key)
                    .header("anthropic-version", ANTHROPIC_API_VERSION),
            };
        }
        let response = match request.send() {
            Ok(response) => response,
            Err(err) => {
                if first_request_error.is_none() {
                    first_request_error = Some(redact_sensitive_text(&format!(
                        "request failed at {endpoint}: {err}"
                    )));
                }
                continue;
            }
        };

        probe.reachable = true;
        let status = response.status();
        let body = match read_live_probe_body(response, runtime.max_response_bytes.max(1)) {
            Ok(body) => body,
            Err(err) => {
                probe.auth_ok = !(status.as_u16() == 401 || status.as_u16() == 403);
                probe.status_class = AgentLiveProbeStatusClass::ProviderError;
                probe.message = redact_sensitive_text(&format!(
                    "Model-list endpoint responded at {endpoint}, but {err}"
                ));
                return probe;
            }
        };

        if !status.is_success() {
            if matches!(status.as_u16(), 404 | 405) && idx + 1 < endpoints.len() {
                continue;
            }
            let (status_class, auth_ok, message, provider_error_code) =
                classify_model_list_http_error(
                    status,
                    endpoint,
                    &body,
                    match auth {
                        ModelListAuth::OpenaiBearer => "Provider",
                        ModelListAuth::Anthropic => "Anthropic",
                    },
                );
            probe.auth_ok = auth_ok;
            probe.status_class = status_class;
            probe.message = message;
            probe.provider_error_code = provider_error_code;
            return probe;
        }

        probe.auth_ok = true;
        probe.selected_endpoint = Some(endpoint.clone());
        let value = match serde_json::from_str::<serde_json::Value>(&body) {
            Ok(value) => value,
            Err(err) => {
                probe.status_class = AgentLiveProbeStatusClass::ProviderError;
                probe.message =
                    format!("Model-list endpoint returned malformed JSON at {endpoint}: {err}");
                return probe;
            }
        };
        let models = extract_models_from_openai_models_payload(&value);
        if models.is_empty() {
            probe.status_class = AgentLiveProbeStatusClass::ProviderError;
            probe.message =
                format!("Model-list endpoint succeeded at {endpoint}, but no model ids were found");
            return probe;
        }

        probe.model_list_ok = true;
        let selected_model = selected_model
            .map(str::trim)
            .filter(|model| !model.is_empty())
            .unwrap_or(OPENAI_COMPAT_DEFAULT_MODEL);
        if is_unspecified_model_name(selected_model) {
            probe.status_class = AgentLiveProbeStatusClass::ModelMissing;
            probe.message = "Model-list probe succeeded, but no concrete model is selected; choose a discovered model or set a model override".to_string();
            return probe;
        }
        if models.iter().any(|model| model == selected_model) {
            probe.selected_model_seen = true;
            probe.status_class = AgentLiveProbeStatusClass::Ok;
            probe.message = format!(
                "Model-list probe succeeded and selected model '{selected_model}' is present"
            );
        } else {
            probe.status_class = AgentLiveProbeStatusClass::ModelMissing;
            probe.message = format!(
                "Model-list probe succeeded, but selected model '{selected_model}' was not returned by this endpoint"
            );
        }
        return probe;
    }

    probe.status_class = AgentLiveProbeStatusClass::EndpointUnreachable;
    probe.message = first_request_error
        .unwrap_or_else(|| "Model-list probe could not reach any endpoint candidate".to_string());
    probe
}

fn build_agent_live_probe(
    system: &AgentSystemSpec,
    availability: &AgentSystemAvailability,
    preflight: &AgentSystemPreflight,
) -> AgentSystemLiveProbe {
    match system.transport {
        AgentSystemTransport::BuiltinEcho | AgentSystemTransport::ExternalJsonStdio => {
            live_probe_with_status(
                AgentLiveProbeStatusClass::UnsupportedTransport,
                format!(
                    "Live endpoint probing is only supported for native_openai, native_anthropic, and native_openai_compat transports, not '{}'",
                    system.transport.as_str()
                ),
            )
        }
        AgentSystemTransport::NativeOpenai => {
            let Some(api_key) = resolve_openai_api_key(system) else {
                return live_probe_with_status(
                    AgentLiveProbeStatusClass::MissingKey,
                    format!(
                        "{OPENAI_API_KEY_ENV} is required for native OpenAI live setup probing"
                    ),
                );
            };
            let Some(base_url) = preflight.base_url.as_deref() else {
                return live_probe_with_status(
                    AgentLiveProbeStatusClass::ProviderError,
                    "No base URL resolved for native OpenAI live setup probing",
                );
            };
            build_model_list_live_probe(
                base_url,
                Some(api_key.as_str()),
                preflight.model.as_deref(),
                &agent_runtime(system),
                ModelListAuth::OpenaiBearer,
            )
        }
        AgentSystemTransport::NativeAnthropic => {
            let Some(api_key) = resolve_anthropic_api_key(system) else {
                return live_probe_with_status(
                    AgentLiveProbeStatusClass::MissingKey,
                    format!(
                        "{ANTHROPIC_API_KEY_ENV} is required for native Claude live setup probing"
                    ),
                );
            };
            let Some(base_url) = preflight.base_url.as_deref() else {
                return live_probe_with_status(
                    AgentLiveProbeStatusClass::ProviderError,
                    "No base URL resolved for native Claude live setup probing",
                );
            };
            build_model_list_live_probe(
                base_url,
                Some(api_key.as_str()),
                preflight.model.as_deref(),
                &agent_runtime(system),
                ModelListAuth::Anthropic,
            )
        }
        AgentSystemTransport::NativeOpenaiCompat => {
            if !availability.available {
                let reason = availability.reason.as_deref().unwrap_or_default();
                if !reason.contains("model is unspecified") {
                    return live_probe_with_status(
                        AgentLiveProbeStatusClass::ProviderError,
                        format!("Config preflight blocked live setup probing: {reason}"),
                    );
                }
            }
            let Some(base_url) = preflight.base_url.as_deref() else {
                return live_probe_with_status(
                    AgentLiveProbeStatusClass::ProviderError,
                    "No base URL resolved for OpenAI-compatible live setup probing",
                );
            };
            let api_key = resolve_openai_api_key(system);
            build_model_list_live_probe(
                base_url,
                api_key.as_deref(),
                preflight.model.as_deref(),
                &agent_runtime(system),
                ModelListAuth::OpenaiBearer,
            )
        }
    }
}

/// Build a deterministic preflight summary for one agent system after applying
/// optional environment-style overrides.
pub fn build_agent_system_preflight(
    catalog_path: Option<&str>,
    system_id: &str,
    env_overrides: Option<&std::collections::HashMap<String, String>>,
) -> Result<AgentSystemPreflight, String> {
    build_agent_system_preflight_with_live(catalog_path, system_id, env_overrides, false)
}

/// Build a deterministic preflight summary and optionally add a live,
/// non-generating model-discovery probe for native HTTP transports.
pub fn build_agent_system_preflight_with_live(
    catalog_path: Option<&str>,
    system_id: &str,
    env_overrides: Option<&std::collections::HashMap<String, String>>,
    live: bool,
) -> Result<AgentSystemPreflight, String> {
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
    let runtime = agent_runtime(&system);
    let mut warnings = vec![];
    let mut base_url: Option<String> = None;
    let mut model: Option<String> = None;
    let mut endpoint_candidates = vec![];
    let mut discovery_endpoint_candidates = vec![];

    match system.transport {
        AgentSystemTransport::BuiltinEcho | AgentSystemTransport::ExternalJsonStdio => {}
        AgentSystemTransport::NativeOpenai => {
            let resolved = resolve_base_url(&system, OPENAI_DEFAULT_BASE_URL);
            base_url = Some(resolved.clone());
            model = Some(resolve_model(&system, OPENAI_DEFAULT_MODEL));
            endpoint_candidates = vec![format!("{resolved}/responses")];
            discovery_endpoint_candidates = model_endpoint_candidates(&resolved);
            if resolve_openai_api_key(&system).is_none() {
                warnings.push(format!("{OPENAI_API_KEY_ENV} is not set"));
            }
        }
        AgentSystemTransport::NativeAnthropic => {
            let resolved = resolve_base_url(&system, ANTHROPIC_DEFAULT_BASE_URL);
            base_url = Some(resolved.clone());
            model = Some(resolve_model(&system, ANTHROPIC_DEFAULT_MODEL));
            endpoint_candidates = vec![format!("{resolved}/messages")];
            discovery_endpoint_candidates = model_endpoint_candidates(&resolved);
            if resolve_anthropic_api_key(&system).is_none() {
                warnings.push(format!("{ANTHROPIC_API_KEY_ENV} is not set"));
            }
        }
        AgentSystemTransport::NativeOpenaiCompat => {
            let resolved = resolve_base_url(&system, OPENAI_COMPAT_DEFAULT_BASE_URL);
            base_url = Some(resolved.clone());
            model = Some(resolve_model(&system, OPENAI_COMPAT_DEFAULT_MODEL));
            endpoint_candidates = openai_compat_endpoint_candidates(&resolved);
            discovery_endpoint_candidates = model_endpoint_candidates(&resolved);
            let host = reqwest::Url::parse(&resolved)
                .ok()
                .and_then(|url| url.host_str().map(|value| value.to_string()))
                .unwrap_or_default();
            if !host.is_empty()
                && !matches!(host.as_str(), "localhost" | "127.0.0.1" | "::1" | "0.0.0.0")
            {
                warnings.push(
                    "remote OpenAI-compatible hosts require explicit base-url override approval"
                        .to_string(),
                );
            }
        }
    }

    let mut preflight = AgentSystemPreflight {
        schema: "gentle.agent_preflight.v1".to_string(),
        catalog_path: resolved_catalog_path,
        system_id: system.id.clone(),
        system_label: system.label.clone(),
        transport: system.transport.as_str().to_string(),
        available: availability.available,
        availability_reason: availability.reason.clone(),
        base_url,
        model,
        endpoint_candidates,
        model_endpoint_candidates: discovery_endpoint_candidates,
        timeout_secs: runtime.timeout_secs,
        connect_timeout_secs: runtime.connect_timeout_secs.unwrap_or_default(),
        read_timeout_secs: runtime.read_timeout_secs.unwrap_or_default(),
        max_retries: runtime.max_retries,
        max_response_bytes: runtime.max_response_bytes,
        warnings,
        live_probe: None,
    };
    if live {
        preflight.live_probe = Some(build_agent_live_probe(&system, &availability, &preflight));
    }
    Ok(preflight)
}

/// Discover model ids for a configured agent system that uses a native HTTP
/// transport with model-list support.
pub fn discover_models_for_agent_system(
    catalog_path: Option<&str>,
    system_id: &str,
    env_overrides: Option<&std::collections::HashMap<String, String>>,
) -> Result<Vec<String>, String> {
    let (_resolved_catalog_path, catalog) = load_agent_system_catalog(catalog_path)?;
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
    let base_url = match system.transport {
        AgentSystemTransport::NativeOpenai => resolve_base_url(&system, OPENAI_DEFAULT_BASE_URL),
        AgentSystemTransport::NativeAnthropic => {
            resolve_base_url(&system, ANTHROPIC_DEFAULT_BASE_URL)
        }
        AgentSystemTransport::NativeOpenaiCompat => {
            resolve_base_url(&system, OPENAI_COMPAT_DEFAULT_BASE_URL)
        }
        other => {
            return Err(format!(
                "agent system '{}' uses transport '{}' which does not support remote model discovery",
                system.id,
                other.as_str()
            ));
        }
    };
    match system.transport {
        AgentSystemTransport::NativeAnthropic => {
            let api_key = resolve_anthropic_api_key(&system);
            discover_anthropic_models(&base_url, api_key.as_deref())
        }
        _ => {
            let api_key = resolve_openai_api_key(&system);
            discover_openai_models(&base_url, api_key.as_deref())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::{
        io::{ErrorKind, Read, Write},
        net::TcpListener,
        thread,
    };

    fn test_runtime() -> AgentInvocationRuntime {
        AgentInvocationRuntime {
            timeout_secs: 2,
            connect_timeout_secs: Some(2),
            read_timeout_secs: Some(2),
            max_retries: 0,
            max_response_bytes: 16 * 1024,
            endpoint_candidates: vec![],
            attempted_endpoints: vec![],
            selected_endpoint: None,
        }
    }

    fn spawn_model_list_server(routes: Vec<(&str, u16, &str)>) -> Option<String> {
        let listener = match TcpListener::bind("127.0.0.1:0") {
            Ok(listener) => listener,
            Err(err) if err.kind() == ErrorKind::PermissionDenied => {
                eprintln!(
                    "skipping local live-probe server test because this environment rejects localhost binds: {err}"
                );
                return None;
            }
            Err(err) => panic!("bind test model server: {err}"),
        };
        let addr = listener.local_addr().expect("test server addr");
        let routes = routes
            .into_iter()
            .map(|(path, status, body)| (path.to_string(), status, body.to_string()))
            .collect::<Vec<_>>();
        thread::spawn(move || {
            for stream in listener.incoming().take(routes.len()) {
                let Ok(mut stream) = stream else {
                    continue;
                };
                let mut buf = [0_u8; 2048];
                let read = stream.read(&mut buf).unwrap_or(0);
                let request = String::from_utf8_lossy(&buf[..read]);
                let path = request.split_whitespace().nth(1).unwrap_or("/");
                let (status, body) = routes
                    .iter()
                    .find(|(route_path, _, _)| route_path == path)
                    .map(|(_, status, body)| (*status, body.as_str()))
                    .unwrap_or((404, r#"{"error":{"code":"not_found"}}"#));
                let reason = match status {
                    200 => "OK",
                    401 => "Unauthorized",
                    403 => "Forbidden",
                    404 => "Not Found",
                    429 => "Too Many Requests",
                    _ => "Error",
                };
                let response = format!(
                    "HTTP/1.1 {status} {reason}\r\nContent-Type: application/json\r\nContent-Length: {}\r\nConnection: close\r\n\r\n{}",
                    body.len(),
                    body
                );
                let _ = stream.write_all(response.as_bytes());
            }
        });
        Some(format!("http://{addr}"))
    }

    #[test]
    fn compat_endpoint_candidates_include_v1_fallback() {
        let candidates = openai_compat_endpoint_candidates("http://127.0.0.1:11434");
        assert_eq!(
            candidates,
            vec![
                "http://127.0.0.1:11434/chat/completions".to_string(),
                "http://127.0.0.1:11434/v1/chat/completions".to_string()
            ]
        );
    }

    #[test]
    fn live_probe_classifies_missing_openai_key_without_endpoint_attempt() {
        let system = AgentSystemSpec {
            id: "openai".to_string(),
            label: "openai".to_string(),
            transport: AgentSystemTransport::NativeOpenai,
            model: Some("gpt-5".to_string()),
            base_url: Some("https://api.openai.com/v1".to_string()),
            ..Default::default()
        };
        let availability = AgentSystemAvailability {
            available: false,
            reason: Some(format!("{OPENAI_API_KEY_ENV} is not set")),
        };
        let preflight = AgentSystemPreflight {
            base_url: Some("https://api.openai.com/v1".to_string()),
            model: Some("gpt-5".to_string()),
            ..Default::default()
        };
        let previous_key = std::env::var(OPENAI_API_KEY_ENV).ok();
        unsafe {
            std::env::remove_var(OPENAI_API_KEY_ENV);
        }
        let probe = build_agent_live_probe(&system, &availability, &preflight);
        if let Some(value) = previous_key {
            unsafe {
                std::env::set_var(OPENAI_API_KEY_ENV, value);
            }
        }
        assert_eq!(probe.status_class, AgentLiveProbeStatusClass::MissingKey);
        assert!(probe.attempted_endpoints.is_empty());
        assert!(!probe.auth_ok);
    }

    #[test]
    fn live_probe_classifies_missing_anthropic_key_without_endpoint_attempt() {
        let system = AgentSystemSpec {
            id: "claude".to_string(),
            label: "Claude".to_string(),
            transport: AgentSystemTransport::NativeAnthropic,
            model: Some("claude-sonnet-4-6".to_string()),
            base_url: Some("https://api.anthropic.com/v1".to_string()),
            ..Default::default()
        };
        let availability = AgentSystemAvailability {
            available: false,
            reason: Some(format!("{ANTHROPIC_API_KEY_ENV} is not set")),
        };
        let preflight = AgentSystemPreflight {
            base_url: Some("https://api.anthropic.com/v1".to_string()),
            model: Some("claude-sonnet-4-6".to_string()),
            ..Default::default()
        };
        let previous_key = std::env::var(ANTHROPIC_API_KEY_ENV).ok();
        unsafe {
            std::env::remove_var(ANTHROPIC_API_KEY_ENV);
        }
        let probe = build_agent_live_probe(&system, &availability, &preflight);
        if let Some(value) = previous_key {
            unsafe {
                std::env::set_var(ANTHROPIC_API_KEY_ENV, value);
            }
        }
        assert_eq!(probe.status_class, AgentLiveProbeStatusClass::MissingKey);
        assert!(probe.attempted_endpoints.is_empty());
        assert!(!probe.auth_ok);
    }

    #[test]
    fn live_probe_classifies_auth_failure() {
        let Some(base_url) = spawn_model_list_server(vec![(
            "/models",
            401,
            r#"{"error":{"code":"invalid_api_key","message":"bad key"}}"#,
        )]) else {
            return;
        };
        let probe = build_model_list_live_probe(
            &base_url,
            Some("sk-test"),
            Some("gpt-5"),
            &test_runtime(),
            ModelListAuth::OpenaiBearer,
        );
        assert_eq!(probe.status_class, AgentLiveProbeStatusClass::AuthFailed);
        assert!(probe.reachable);
        assert!(!probe.auth_ok);
        assert_eq!(
            probe.provider_error_code.as_deref(),
            Some("invalid_api_key")
        );
    }

    #[test]
    fn live_probe_classifies_insufficient_quota() {
        let Some(base_url) = spawn_model_list_server(vec![(
            "/models",
            429,
            r#"{"error":{"type":"insufficient_quota","code":"insufficient_quota"}}"#,
        )]) else {
            return;
        };
        let probe = build_model_list_live_probe(
            &base_url,
            Some("sk-test"),
            Some("gpt-5"),
            &test_runtime(),
            ModelListAuth::OpenaiBearer,
        );
        assert_eq!(
            probe.status_class,
            AgentLiveProbeStatusClass::QuotaOrBilling
        );
        assert!(probe.auth_ok);
        assert_eq!(
            probe.provider_error_code.as_deref(),
            Some("insufficient_quota")
        );
    }

    #[test]
    fn live_probe_classifies_malformed_model_json() {
        let Some(base_url) = spawn_model_list_server(vec![("/models", 200, "not json")]) else {
            return;
        };
        let probe = build_model_list_live_probe(
            &base_url,
            None,
            Some("llama3"),
            &test_runtime(),
            ModelListAuth::OpenaiBearer,
        );
        assert_eq!(probe.status_class, AgentLiveProbeStatusClass::ProviderError);
        assert!(probe.reachable);
        assert!(probe.auth_ok);
        assert!(!probe.model_list_ok);
        assert!(probe.message.contains("malformed JSON"));
    }

    #[test]
    fn live_probe_classifies_selected_model_absent() {
        let Some(base_url) = spawn_model_list_server(vec![(
            "/models",
            200,
            r#"{"data":[{"id":"llama3"},{"id":"qwen3"}]}"#,
        )]) else {
            return;
        };
        let probe = build_model_list_live_probe(
            &base_url,
            None,
            Some("missing-model"),
            &test_runtime(),
            ModelListAuth::OpenaiBearer,
        );
        assert_eq!(probe.status_class, AgentLiveProbeStatusClass::ModelMissing);
        assert!(probe.model_list_ok);
        assert!(!probe.selected_model_seen);
    }

    #[test]
    fn live_probe_falls_back_to_v1_models_for_compat_roots() {
        let Some(base_url) = spawn_model_list_server(vec![
            ("/models", 404, r#"{"error":{"code":"not_found"}}"#),
            ("/v1/models", 200, r#"{"data":[{"id":"llama3"}]}"#),
        ]) else {
            return;
        };
        let probe = build_model_list_live_probe(
            &base_url,
            None,
            Some("llama3"),
            &test_runtime(),
            ModelListAuth::OpenaiBearer,
        );
        assert_eq!(probe.status_class, AgentLiveProbeStatusClass::Ok);
        assert_eq!(probe.attempted_endpoints.len(), 2);
        assert!(
            probe
                .selected_endpoint
                .as_deref()
                .unwrap_or_default()
                .ends_with("/v1/models")
        );
        assert!(probe.selected_model_seen);
    }

    #[test]
    fn live_probe_supports_anthropic_model_list_shape() {
        let Some(base_url) = spawn_model_list_server(vec![(
            "/models",
            200,
            r#"{"data":[{"id":"claude-sonnet-4-6","type":"model"}]}"#,
        )]) else {
            return;
        };
        let probe = build_model_list_live_probe(
            &base_url,
            Some("sk-ant-test"),
            Some("claude-sonnet-4-6"),
            &test_runtime(),
            ModelListAuth::Anthropic,
        );
        assert_eq!(probe.status_class, AgentLiveProbeStatusClass::Ok);
        assert!(probe.model_list_ok);
        assert!(probe.selected_model_seen);
    }
}
