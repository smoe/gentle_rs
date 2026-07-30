//! Process-local tool-path override registry.

use std::{
    collections::HashMap,
    sync::{LazyLock, RwLock},
};

#[cfg(test)]
use std::{cell::RefCell, marker::PhantomData, rc::Rc};

static TOOL_OVERRIDES: LazyLock<RwLock<HashMap<String, String>>> =
    LazyLock::new(|| RwLock::new(HashMap::new()));

#[cfg(test)]
thread_local! {
    static SCOPED_TOOL_OVERRIDES: RefCell<HashMap<String, String>> =
        RefCell::new(HashMap::new());
}

#[cfg(test)]
pub(crate) struct ScopedToolOverrideGuard {
    env_var: String,
    previous: Option<String>,
    _not_send: PhantomData<Rc<()>>,
}

#[cfg(test)]
impl ScopedToolOverrideGuard {
    pub(crate) fn set(env_var: &str, configured: &str) -> Self {
        let configured =
            normalized_non_empty(configured).expect("scoped test tool override must not be empty");
        let previous = SCOPED_TOOL_OVERRIDES.with(|overrides| {
            overrides
                .borrow_mut()
                .insert(env_var.to_string(), configured)
        });
        Self {
            env_var: env_var.to_string(),
            previous,
            _not_send: PhantomData,
        }
    }
}

#[cfg(test)]
impl Drop for ScopedToolOverrideGuard {
    fn drop(&mut self) {
        SCOPED_TOOL_OVERRIDES.with(|overrides| {
            let mut overrides = overrides.borrow_mut();
            if let Some(previous) = self.previous.take() {
                overrides.insert(self.env_var.clone(), previous);
            } else {
                overrides.remove(&self.env_var);
            }
        });
    }
}

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
        .unwrap_or_else(|poisoned| poisoned.into_inner());
    if let Some(value) = normalized_non_empty(configured) {
        guard.insert(env_var.to_string(), value);
    } else {
        guard.remove(env_var);
    }
}

pub fn get_tool_override(env_var: &str) -> Option<String> {
    #[cfg(test)]
    if let Some(value) =
        SCOPED_TOOL_OVERRIDES.with(|overrides| overrides.borrow().get(env_var).cloned())
    {
        return Some(value);
    }
    TOOL_OVERRIDES
        .read()
        .unwrap_or_else(|poisoned| poisoned.into_inner())
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scoped_test_override_takes_precedence_and_restores_global_override() {
        const ENV_VAR: &str = "GENTLE_SCOPED_TOOL_OVERRIDE_TEST";
        set_tool_override(ENV_VAR, "global-tool");
        assert_eq!(get_tool_override(ENV_VAR).as_deref(), Some("global-tool"));

        {
            let _scoped = ScopedToolOverrideGuard::set(ENV_VAR, "scoped-tool");
            assert_eq!(get_tool_override(ENV_VAR).as_deref(), Some("scoped-tool"));
        }

        assert_eq!(get_tool_override(ENV_VAR).as_deref(), Some("global-tool"));
        set_tool_override(ENV_VAR, "");
    }
}
