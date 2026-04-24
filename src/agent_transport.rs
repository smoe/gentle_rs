//! Shared transport/runtime helpers for GENtle agent systems.
//!
//! This module exposes the machine-facing transport metadata needed by both the
//! local assistant UX and external orchestrators. It intentionally keeps the
//! public surface read-only and deterministic: catalog loading, availability,
//! preflight summaries, and OpenAI-compatible model discovery.

use crate::agent_bridge::{
    AGENT_BASE_URL_ENV, AGENT_CONNECT_TIMEOUT_SECS_ENV, AGENT_MAX_RESPONSE_BYTES_ENV,
    AGENT_MAX_RETRIES_ENV, AGENT_MODEL_ENV, AGENT_READ_TIMEOUT_SECS_ENV, AGENT_TIMEOUT_SECS_ENV,
    OPENAI_API_KEY_ENV,
};
use serde::{Deserialize, Serialize};

pub use crate::agent_bridge::{
    AgentInvocationRuntime, AgentSystemAvailability, AgentSystemCatalog, AgentSystemSpec,
    AgentSystemTransport, DEFAULT_AGENT_SYSTEM_CATALOG_PATH, agent_system_availability,
    discover_openai_models, load_agent_system_catalog,
};

const AGENT_REQUEST_TIMEOUT_SECS_DEFAULT: u64 = 180;
const AGENT_CONNECT_TIMEOUT_SECS_DEFAULT: u64 = 10;
const AGENT_MAX_RETRIES_DEFAULT: usize = 2;
const AGENT_MAX_RESPONSE_BYTES_DEFAULT: usize = 1_048_576;
const OPENAI_DEFAULT_MODEL: &str = "gpt-5";
const OPENAI_DEFAULT_BASE_URL: &str = "https://api.openai.com/v1";
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

/// Build a deterministic preflight summary for one agent system after applying
/// optional environment-style overrides.
pub fn build_agent_system_preflight(
    catalog_path: Option<&str>,
    system_id: &str,
    env_overrides: Option<&std::collections::HashMap<String, String>>,
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
            if resolve_api_key(&system).is_none() {
                warnings.push(format!("{OPENAI_API_KEY_ENV} is not set"));
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

    Ok(AgentSystemPreflight {
        schema: "gentle.agent_preflight.v1".to_string(),
        catalog_path: resolved_catalog_path,
        system_id: system.id.clone(),
        system_label: system.label.clone(),
        transport: system.transport.as_str().to_string(),
        available: availability.available,
        availability_reason: availability.reason,
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
    })
}

/// Discover model ids for a configured agent system that uses a native
/// OpenAI/OpenAI-compatible HTTP transport.
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
    let api_key = resolve_api_key(&system);
    discover_openai_models(&base_url, api_key.as_deref())
}

#[cfg(test)]
mod tests {
    use super::*;

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
}
