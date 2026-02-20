use std::{
    collections::HashMap,
    sync::{LazyLock, RwLock},
};

static TOOL_OVERRIDES: LazyLock<RwLock<HashMap<String, String>>> =
    LazyLock::new(|| RwLock::new(HashMap::new()));

fn normalized_non_empty(value: &str) -> Option<String> {
    let trimmed = value.trim();
    if trimmed.is_empty() {
        None
    } else {
        Some(trimmed.to_string())
    }
}

pub fn set_tool_override(env_var: &str, configured: &str) {
    let mut guard = TOOL_OVERRIDES
        .write()
        .expect("Tool override lock poisoned for write");
    if let Some(value) = normalized_non_empty(configured) {
        guard.insert(env_var.to_string(), value);
    } else {
        guard.remove(env_var);
    }
}

pub fn get_tool_override(env_var: &str) -> Option<String> {
    TOOL_OVERRIDES
        .read()
        .expect("Tool override lock poisoned for read")
        .get(env_var)
        .cloned()
}

pub fn configured_or_env(env_var: &str) -> String {
    get_tool_override(env_var)
        .or_else(|| {
            std::env::var(env_var)
                .ok()
                .and_then(|v| normalized_non_empty(&v))
        })
        .unwrap_or_default()
}

pub fn resolve_tool_executable(env_var: &str, default_bin: &str) -> String {
    get_tool_override(env_var)
        .or_else(|| {
            std::env::var(env_var)
                .ok()
                .and_then(|v| normalized_non_empty(&v))
        })
        .unwrap_or_else(|| default_bin.to_string())
}

pub fn active_resolution_label(env_var: &str, default_bin: &str) -> String {
    get_tool_override(env_var)
        .or_else(|| {
            std::env::var(env_var)
                .ok()
                .and_then(|v| normalized_non_empty(&v))
        })
        .unwrap_or_else(|| format!("PATH lookup: {default_bin}"))
}
