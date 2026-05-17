//! Routine Assistant and Agent Assistant GUI helpers.
//!
//! This module is a move-only extraction from `app.rs`: it keeps the
//! intertwined Routine Assistant and Agent Assistant dialog/rendering helpers
//! close to `GENtleApp` while reducing the top-level app monolith.

use super::*;

impl GENtleApp {
    pub(super) fn open_routine_assistant_dialog(&mut self) {
        if self.show_routine_assistant_dialog {
            self.mark_window_open_or_focus(Self::routine_assistant_viewport_id(), true);
            return;
        }
        self.show_routine_assistant_dialog = true;
        self.mark_window_open_or_focus(Self::routine_assistant_viewport_id(), false);
        if self.routine_assistant_candidates.is_empty() {
            self.refresh_routine_assistant_candidates();
        }
        self.ensure_routine_assistant_decision_trace_started();
        let bindings_snapshot = self.routine_assistant_bindings_snapshot();
        self.update_routine_assistant_decision_trace(|trace| {
            trace.status = "draft".to_string();
            trace.bindings_snapshot = bindings_snapshot;
        });
    }

    pub(super) fn open_agent_assistant_dialog(&mut self) {
        self.refresh_agent_system_catalog();
        let was_open = self.show_agent_assistant_dialog;
        self.show_agent_assistant_dialog = true;
        self.mark_window_open_or_focus(Self::agent_assistant_viewport_id(), was_open);
    }
    pub(super) fn refresh_agent_system_catalog(&mut self) {
        let catalog_path = self.agent_catalog_path.trim().to_string();
        if !self.agent_systems.is_empty()
            && self.agent_catalog_loaded_path == catalog_path
            && self.agent_catalog_error.is_empty()
        {
            return;
        }
        self.agent_catalog_loaded_path = catalog_path.clone();
        match load_agent_system_catalog(Some(&catalog_path)) {
            Ok((_resolved, catalog)) => {
                self.agent_systems = catalog.systems;
                self.agent_catalog_error.clear();
                if self.agent_system_id.trim().is_empty()
                    || !self
                        .agent_systems
                        .iter()
                        .any(|system| system.id == self.agent_system_id)
                {
                    self.agent_system_id = self
                        .agent_systems
                        .first()
                        .map(|system| system.id.clone())
                        .unwrap_or_default();
                }
            }
            Err(err) => {
                self.agent_catalog_error = err;
                self.agent_systems.clear();
                self.agent_system_id.clear();
            }
        }
    }

    pub(super) fn selected_agent_system(&self) -> Option<AgentSystemSpec> {
        self.agent_systems
            .iter()
            .find(|system| system.id == self.agent_system_id)
            .cloned()
    }

    pub(super) fn clear_agent_preflight_output(&mut self) {
        self.agent_preflight_output = None;
    }

    pub(super) fn invalidate_agent_preflight_after_setup_input_change(&mut self) {
        self.clear_agent_preflight_output();
    }

    pub(super) fn clear_agent_model_discovery_snapshot(&mut self) {
        self.agent_model_discovery_task = None;
        self.agent_discovered_models.clear();
        self.agent_discovered_model_pick.clear();
        self.agent_model_discovery_status.clear();
        self.agent_model_discovery_source_key.clear();
        self.agent_model_discovery_failed_source_key.clear();
    }

    pub(super) fn select_agent_system_and_reset_setup(&mut self, system_id: &str) {
        if self.agent_system_id == system_id {
            return;
        }
        self.agent_system_id = system_id.to_string();
        self.clear_agent_preflight_output();
        self.clear_agent_model_discovery_snapshot();
    }

    pub(super) fn selected_agent_session_env_overrides(
        &self,
        system: &AgentSystemSpec,
    ) -> Result<HashMap<String, String>, String> {
        let mut overrides = HashMap::new();
        let session_api_key = self.agent_openai_api_key.trim();
        if !session_api_key.is_empty() {
            let key_env = match system.transport {
                AgentSystemTransport::NativeAnthropic => ANTHROPIC_API_KEY_ENV,
                _ => OPENAI_API_KEY_ENV,
            };
            overrides.insert(key_env.to_string(), session_api_key.to_string());
        }
        let override_base_url = self.agent_base_url_override.trim();
        if !override_base_url.is_empty()
            && matches!(
                system.transport,
                AgentSystemTransport::NativeOpenai
                    | AgentSystemTransport::NativeAnthropic
                    | AgentSystemTransport::NativeOpenaiCompat
            )
        {
            overrides.insert(
                AGENT_BASE_URL_ENV.to_string(),
                override_base_url.to_string(),
            );
        }
        let selected_discovered_model =
            normalize_agent_model_name(&self.agent_discovered_model_pick).filter(|picked| {
                self.agent_discovered_models
                    .iter()
                    .any(|item| item == picked)
            });
        let override_model = normalize_agent_model_name(self.agent_model_override.trim())
            .or(selected_discovered_model);
        if let Some(override_model) = override_model
            && matches!(
                system.transport,
                AgentSystemTransport::NativeOpenai
                    | AgentSystemTransport::NativeAnthropic
                    | AgentSystemTransport::NativeOpenaiCompat
            )
        {
            overrides.insert(AGENT_MODEL_ENV.to_string(), override_model);
        }
        if let Some(timeout_override) = self.parse_agent_timeout_seconds()?
            && matches!(
                system.transport,
                AgentSystemTransport::ExternalJsonStdio
                    | AgentSystemTransport::NativeOpenai
                    | AgentSystemTransport::NativeAnthropic
                    | AgentSystemTransport::NativeOpenaiCompat
            )
        {
            overrides.insert(
                AGENT_TIMEOUT_SECS_ENV.to_string(),
                timeout_override.to_string(),
            );
        }
        if let Some(connect_timeout_override) = self.parse_agent_connect_timeout_seconds()?
            && matches!(
                system.transport,
                AgentSystemTransport::NativeOpenai
                    | AgentSystemTransport::NativeAnthropic
                    | AgentSystemTransport::NativeOpenaiCompat
            )
        {
            overrides.insert(
                AGENT_CONNECT_TIMEOUT_SECS_ENV.to_string(),
                connect_timeout_override.to_string(),
            );
        }
        if let Some(read_timeout_override) = self.parse_agent_read_timeout_seconds()?
            && matches!(
                system.transport,
                AgentSystemTransport::ExternalJsonStdio
                    | AgentSystemTransport::NativeOpenai
                    | AgentSystemTransport::NativeAnthropic
                    | AgentSystemTransport::NativeOpenaiCompat
            )
        {
            overrides.insert(
                AGENT_READ_TIMEOUT_SECS_ENV.to_string(),
                read_timeout_override.to_string(),
            );
        }
        if let Some(max_retries_override) = self.parse_agent_max_retries()?
            && matches!(
                system.transport,
                AgentSystemTransport::ExternalJsonStdio
                    | AgentSystemTransport::NativeOpenai
                    | AgentSystemTransport::NativeAnthropic
                    | AgentSystemTransport::NativeOpenaiCompat
            )
        {
            overrides.insert(
                AGENT_MAX_RETRIES_ENV.to_string(),
                max_retries_override.to_string(),
            );
        }
        if let Some(max_response_bytes_override) = self.parse_agent_max_response_bytes()?
            && matches!(
                system.transport,
                AgentSystemTransport::ExternalJsonStdio
                    | AgentSystemTransport::NativeOpenai
                    | AgentSystemTransport::NativeAnthropic
                    | AgentSystemTransport::NativeOpenaiCompat
            )
        {
            overrides.insert(
                AGENT_MAX_RESPONSE_BYTES_ENV.to_string(),
                max_response_bytes_override.to_string(),
            );
        }
        Ok(overrides)
    }

    pub(super) fn selected_agent_system_with_session_overrides(
        &self,
        system: &AgentSystemSpec,
    ) -> Result<AgentSystemSpec, String> {
        let mut resolved = system.clone();
        for (key, value) in self.selected_agent_session_env_overrides(system)? {
            resolved.env.insert(key, value);
        }
        Ok(resolved)
    }

    pub(super) fn selected_agent_runtime_base_url(
        &self,
        system: &AgentSystemSpec,
    ) -> Option<String> {
        if !matches!(
            system.transport,
            AgentSystemTransport::NativeOpenai
                | AgentSystemTransport::NativeAnthropic
                | AgentSystemTransport::NativeOpenaiCompat
        ) {
            return None;
        }
        let override_base_url = self.agent_base_url_override.trim();
        if !override_base_url.is_empty() {
            return Some(override_base_url.to_string());
        }
        if let Some(catalog_base_url) = system
            .base_url
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            return Some(catalog_base_url.to_string());
        }
        Some(match system.transport {
            AgentSystemTransport::NativeOpenai => GUI_OPENAI_DEFAULT_BASE_URL.to_string(),
            AgentSystemTransport::NativeAnthropic => GUI_ANTHROPIC_DEFAULT_BASE_URL.to_string(),
            AgentSystemTransport::NativeOpenaiCompat => {
                GUI_OPENAI_COMPAT_DEFAULT_BASE_URL.to_string()
            }
            _ => return None,
        })
    }

    pub(super) fn selected_agent_model_discovery_source_key(
        &self,
        system: &AgentSystemSpec,
    ) -> Option<String> {
        let base_url = self.selected_agent_runtime_base_url(system)?;
        let key_state = self.selected_agent_model_discovery_key_label(system);
        Some(format!(
            "{}|{}|{}|{}",
            system.id,
            system.transport.as_str(),
            base_url,
            key_state
        ))
    }

    pub(super) fn selected_agent_model_discovery_key_label(
        &self,
        system: &AgentSystemSpec,
    ) -> &'static str {
        if !self.agent_openai_api_key.trim().is_empty() {
            "session-key"
        } else {
            let (env_key, label) = match system.transport {
                AgentSystemTransport::NativeAnthropic => {
                    (ANTHROPIC_API_KEY_ENV, "env-anthropic-api-key")
                }
                _ => (OPENAI_API_KEY_ENV, "env-openai-api-key"),
            };
            if std::env::var(env_key)
                .ok()
                .map(|value| !value.trim().is_empty())
                .unwrap_or(false)
            {
                label
            } else {
                "no-key"
            }
        }
    }

    pub(super) fn agent_model_discovery_failure_hint(error: &str) -> Option<&'static str> {
        let lower = error.to_ascii_lowercase();
        let auth_failed = lower.contains("401")
            || lower.contains("403")
            || lower.contains("unauthorized")
            || lower.contains("invalid_api_key")
            || lower.contains("incorrect api key")
            || lower.contains("authentication_error");
        if auth_failed
            && (lower.contains("anthropic")
                || lower.contains("x-api-key")
                || lower.contains("claude code")
                || lower.contains("claude.ai"))
        {
            if lower.contains("claude code/claude.ai") {
                return None;
            }
            return Some(ANTHROPIC_API_KEY_AUTH_HINT);
        }
        if auth_failed {
            return Some("Authentication failed. Use an OpenAI Platform API key for OPENAI_API_KEY; ChatGPT/Codex subscription tokens are not OpenAI API keys.");
        }
        if lower.contains("timed out") || lower.contains("timeout") {
            return Some(
                "The endpoint did not answer before the model-list timeout; check the Base URL or local server.",
            );
        }
        if lower.contains("connection refused")
            || lower.contains("could not connect")
            || lower.contains("dns")
        {
            return Some(
                "The model-list endpoint could not be reached; check the Base URL or start the local OpenAI-compatible server.",
            );
        }
        None
    }

    pub(super) fn agent_test_setup_uses_live_probe(system: &AgentSystemSpec) -> bool {
        matches!(
            system.transport,
            AgentSystemTransport::NativeOpenai
                | AgentSystemTransport::NativeAnthropic
                | AgentSystemTransport::NativeOpenaiCompat
        )
    }

    pub(super) fn shell_quote_command_arg(raw: &str) -> String {
        if raw
            .chars()
            .all(|ch| ch.is_ascii_alphanumeric() || matches!(ch, '/' | '.' | '_' | '-'))
        {
            return raw.to_string();
        }
        format!("'{}'", raw.replace('\'', "'\\''"))
    }

    pub(super) fn external_agent_mcp_state_path(&self) -> String {
        self.current_project_path
            .as_deref()
            .map(str::trim)
            .filter(|path| !path.is_empty())
            .unwrap_or(DEFAULT_MCP_STATE_PATH)
            .to_string()
    }

    pub(super) fn external_agent_mcp_command_snippet_for_state_path(state_path: &str) -> String {
        format!(
            "gentle_mcp --state {}",
            Self::shell_quote_command_arg(state_path)
        )
    }

    pub(super) fn external_agent_mcp_command_snippet(&self) -> String {
        Self::external_agent_mcp_command_snippet_for_state_path(
            &self.external_agent_mcp_state_path(),
        )
    }

    pub(super) fn agent_preflight_next_actions(preflight: &AgentSystemPreflight) -> Vec<String> {
        let key_hint = if preflight.transport == AgentSystemTransport::NativeAnthropic.as_str() {
            format!("Paste an Anthropic API key or set {ANTHROPIC_API_KEY_ENV}.")
        } else {
            format!(
                "Paste a session key or set {OPENAI_API_KEY_ENV}; ChatGPT/Codex subscriptions are not OpenAI API keys."
            )
        };
        if let Some(live) = &preflight.live_probe {
            let model_is_unspecified = preflight
                .model
                .as_deref()
                .map(str::trim)
                .map(|model| {
                    model.is_empty() || model.eq_ignore_ascii_case(OPENAI_COMPAT_UNSPECIFIED_MODEL)
                })
                .unwrap_or(true);
            return match live.status_class {
                AgentLiveProbeStatusClass::Ok => vec![],
                AgentLiveProbeStatusClass::MissingKey => vec![key_hint],
                AgentLiveProbeStatusClass::AuthFailed => {
                    if preflight.transport == AgentSystemTransport::NativeAnthropic.as_str() {
                        vec![ANTHROPIC_API_KEY_AUTH_HINT.to_string()]
                    } else {
                        vec![
                            "Check the API key/token for this endpoint, then run Test Setup again."
                                .to_string(),
                        ]
                    }
                }
                AgentLiveProbeStatusClass::QuotaOrBilling => vec![
                    "Check provider billing/quota; this setup probe did not intentionally generate tokens."
                        .to_string(),
                ],
                AgentLiveProbeStatusClass::ModelMissing => {
                    if model_is_unspecified {
                        vec![
                            "Pick a discovered model or set Model override before asking the assistant."
                                .to_string(),
                        ]
                    } else {
                        vec![
                            "Choose a model returned by this endpoint, or correct Base URL if the model list came from the wrong server."
                                .to_string(),
                        ]
                    }
                }
                AgentLiveProbeStatusClass::EndpointUnreachable => vec![
                    "Start the local server or correct Base URL override, then run Test Setup again."
                        .to_string(),
                ],
                AgentLiveProbeStatusClass::UnsupportedTransport => vec![
                    "Use this setup check as config-only validation for this non-HTTP transport."
                        .to_string(),
                ],
                AgentLiveProbeStatusClass::ProviderError => vec![
                    "Inspect the provider response; model discovery must return JSON with model ids."
                        .to_string(),
                ],
            };
        }

        let mut actions = Vec::new();
        if preflight.warnings.iter().any(|warning| {
            warning.contains(OPENAI_API_KEY_ENV) || warning.contains(ANTHROPIC_API_KEY_ENV)
        }) {
            actions.push(key_hint);
        }
        if preflight
            .availability_reason
            .as_deref()
            .unwrap_or_default()
            .contains("model is unspecified")
        {
            actions.push("Pick a discovered model or set Model override.".to_string());
        }
        actions
    }

    pub(super) fn parse_agent_timeout_seconds(&self) -> Result<Option<u64>, String> {
        let raw = self.agent_timeout_secs.trim();
        if raw.is_empty() {
            return Ok(None);
        }
        let parsed = raw
            .parse::<u64>()
            .map_err(|e| format!("Invalid timeout_sec '{}': {}", raw, e))?;
        if parsed == 0 {
            return Ok(None);
        }
        Ok(Some(parsed))
    }

    pub(super) fn parse_agent_connect_timeout_seconds(&self) -> Result<Option<u64>, String> {
        let raw = self.agent_connect_timeout_secs.trim();
        if raw.is_empty() {
            return Ok(None);
        }
        let parsed = raw
            .parse::<u64>()
            .map_err(|e| format!("Invalid connect_timeout_sec '{}': {}", raw, e))?;
        if parsed == 0 {
            return Ok(None);
        }
        Ok(Some(parsed))
    }

    pub(super) fn parse_agent_read_timeout_seconds(&self) -> Result<Option<u64>, String> {
        let raw = self.agent_read_timeout_secs.trim();
        if raw.is_empty() {
            return Ok(None);
        }
        let parsed = raw
            .parse::<u64>()
            .map_err(|e| format!("Invalid read_timeout_sec '{}': {}", raw, e))?;
        if parsed == 0 {
            return Ok(None);
        }
        Ok(Some(parsed))
    }

    pub(super) fn parse_agent_max_retries(&self) -> Result<Option<usize>, String> {
        let raw = self.agent_max_retries.trim();
        if raw.is_empty() {
            return Ok(None);
        }
        let parsed = raw
            .parse::<usize>()
            .map_err(|e| format!("Invalid max_retries '{}': {}", raw, e))?;
        Ok(Some(parsed))
    }

    pub(super) fn parse_agent_max_response_bytes(&self) -> Result<Option<usize>, String> {
        let raw = self.agent_max_response_bytes.trim();
        if raw.is_empty() {
            return Ok(None);
        }
        let parsed = raw
            .parse::<usize>()
            .map_err(|e| format!("Invalid max_response_bytes '{}': {}", raw, e))?;
        if parsed == 0 {
            return Ok(None);
        }
        Ok(Some(parsed))
    }

    pub(super) fn start_agent_model_discovery_task(
        &mut self,
        system: &AgentSystemSpec,
        force: bool,
    ) {
        if !matches!(
            system.transport,
            AgentSystemTransport::NativeOpenai
                | AgentSystemTransport::NativeAnthropic
                | AgentSystemTransport::NativeOpenaiCompat
        ) {
            return;
        }
        let Some(base_url) = self.selected_agent_runtime_base_url(system) else {
            return;
        };
        let Some(source_key) = self.selected_agent_model_discovery_source_key(system) else {
            return;
        };
        if !force {
            if let Some(task) = &self.agent_model_discovery_task {
                if task.source_key == source_key {
                    return;
                }
            }
            if self.agent_model_discovery_source_key == source_key
                && !self.agent_discovered_models.is_empty()
            {
                return;
            }
            if self.agent_model_discovery_failed_source_key == source_key {
                return;
            }
        }
        self.agent_model_discovery_failed_source_key.clear();
        self.agent_model_discovery_source_key = source_key.clone();
        let key_label = self.selected_agent_model_discovery_key_label(system);
        self.agent_model_discovery_status = format!(
            "Discovering models at {base_url} (auth={key_label}; timeout about 20s per endpoint) ..."
        );
        self.agent_model_discovery_task = None;
        let env_overrides = match self.selected_agent_session_env_overrides(system) {
            Ok(overrides) => overrides,
            Err(err) => {
                self.agent_model_discovery_status = err;
                return;
            }
        };
        let catalog_path = self.agent_catalog_path.trim().to_string();
        let system_id = system.id.clone();
        let (tx, rx) = mpsc::channel::<AgentModelDiscoveryTaskMessage>();
        self.agent_model_discovery_task = Some(AgentModelDiscoveryTask {
            started: Instant::now(),
            source_key: source_key.clone(),
            receiver: rx,
        });
        std::thread::spawn(move || {
            let result = discover_models_for_agent_system(
                Some(catalog_path.as_str()),
                &system_id,
                if env_overrides.is_empty() {
                    None
                } else {
                    Some(&env_overrides)
                },
            );
            let _ = tx.send(AgentModelDiscoveryTaskMessage::Done { source_key, result });
        });
    }

    pub(super) fn selected_agent_system_availability(
        &self,
        system: &AgentSystemSpec,
    ) -> (bool, Option<String>) {
        let resolved = match self.selected_agent_system_with_session_overrides(system) {
            Ok(resolved) => resolved,
            Err(err) => return (false, Some(err)),
        };
        let availability = agent_system_availability(&resolved);
        (availability.available, availability.reason)
    }

    pub(super) fn run_agent_preflight_probe(&mut self) {
        self.refresh_agent_system_catalog();
        self.clear_agent_preflight_output();
        if !self.agent_catalog_error.is_empty() {
            self.agent_status = format!("Agent catalog error: {}", self.agent_catalog_error);
            return;
        }
        let Some(selected_system) = self.selected_agent_system() else {
            self.agent_status = "Select an agent system first".to_string();
            return;
        };
        let env_overrides = match self.selected_agent_session_env_overrides(&selected_system) {
            Ok(overrides) => overrides,
            Err(err) => {
                self.agent_status = err;
                return;
            }
        };
        let live_probe = Self::agent_test_setup_uses_live_probe(&selected_system);
        match build_agent_system_preflight_with_live(
            Some(self.agent_catalog_path.trim()),
            selected_system.id.as_str(),
            if env_overrides.is_empty() {
                None
            } else {
                Some(&env_overrides)
            },
            live_probe,
        ) {
            Ok(preflight) => {
                let availability = if preflight.available {
                    "available"
                } else {
                    "unavailable"
                };
                let live_status = preflight
                    .live_probe
                    .as_ref()
                    .map(|probe| format!(", live={}", probe.status_class.as_str()))
                    .unwrap_or_default();
                self.agent_status = format!(
                    "Agent setup preflight: {} ({}{})",
                    selected_system.id, availability, live_status
                );
                self.agent_preflight_output = Some(preflight);
            }
            Err(err) => {
                self.agent_status = format!("Agent setup preflight failed: {err}");
            }
        }
    }

    pub(super) fn start_agent_assistant_request(&mut self) {
        if self.agent_task.is_some() {
            self.agent_status = "Agent request is already running".to_string();
            return;
        }
        self.refresh_agent_system_catalog();
        if !self.agent_catalog_error.is_empty() {
            self.agent_status = format!("Agent catalog error: {}", self.agent_catalog_error);
            return;
        }
        let system_id = self.agent_system_id.trim().to_string();
        if system_id.is_empty() {
            self.agent_status = "Select an agent system first".to_string();
            return;
        }
        let Some(selected_system) = self.selected_agent_system() else {
            self.agent_status = "Selected agent system is not available in catalog".to_string();
            return;
        };
        let (available, reason) = self.selected_agent_system_availability(&selected_system);
        if !available {
            self.agent_status = format!(
                "Selected agent system is unavailable: {}",
                reason.unwrap_or_else(|| "unknown reason".to_string())
            );
            return;
        }
        let prompt = self.agent_prompt.trim().to_string();
        if prompt.is_empty() {
            self.agent_status = "Agent prompt cannot be empty".to_string();
            return;
        }
        let env_overrides = match self.selected_agent_session_env_overrides(&selected_system) {
            Ok(overrides) => overrides,
            Err(err) => {
                self.agent_status = err;
                return;
            }
        };
        let timeout_seconds = self.parse_agent_timeout_seconds().ok().flatten();
        let max_retries = self.parse_agent_max_retries().ok().flatten();
        let resolved_runtime_model = env_overrides
            .get(AGENT_MODEL_ENV)
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty());
        if matches!(
            selected_system.transport,
            AgentSystemTransport::NativeOpenaiCompat
        ) && resolved_runtime_model.is_none()
        {
            let catalog_model =
                normalize_agent_model_name(selected_system.model.as_deref().unwrap_or_default());
            if let Some(catalog_model) = catalog_model {
                if !self.agent_discovered_models.is_empty()
                    && !self
                        .agent_discovered_models
                        .iter()
                        .any(|value| value == &catalog_model)
                {
                    self.agent_status = format!(
                        "Catalog model '{catalog_model}' is not available on current endpoint. Select a discovered model or set Model override."
                    );
                    return;
                }
            } else {
                self.agent_status =
                    "Model is unspecified. Discover models and select one, or set Model override."
                        .to_string();
                return;
            }
        }

        let state_summary = if self.agent_include_state_summary {
            Some(self.engine.read().unwrap().summarize_state())
        } else {
            None
        };
        let catalog_path = self.agent_catalog_path.trim().to_string();
        let job_id = self.alloc_background_job_id();
        let (tx, rx) = mpsc::channel::<AgentAskTaskMessage>();
        self.agent_status = if let Some(timeout) = timeout_seconds {
            format!(
                "Asking agent '{}' in background (timeout={}s, retries={})",
                system_id,
                timeout,
                max_retries.unwrap_or(2)
            )
        } else {
            format!("Asking agent '{}' in background", system_id)
        };
        self.push_job_event(
            BackgroundJobKind::AgentAssist,
            BackgroundJobEventPhase::Started,
            Some(job_id),
            format!("Agent request started for system '{}'", system_id),
        );
        self.agent_task = Some(AgentAskTask {
            job_id,
            started: Instant::now(),
            receiver: rx,
        });
        std::thread::spawn(move || {
            let result = invoke_agent_support_with_env_overrides(
                Some(catalog_path.as_str()),
                &system_id,
                &prompt,
                state_summary.as_ref(),
                if env_overrides.is_empty() {
                    None
                } else {
                    Some(&env_overrides)
                },
            );
            let _ = tx.send(AgentAskTaskMessage::Done { job_id, result });
        });
    }

    pub(super) fn execute_agent_suggested_command(
        &mut self,
        index_1based: usize,
        command_text: &str,
        trigger: &str,
    ) {
        let trimmed = command_text.trim();
        if trimmed.is_empty() {
            self.agent_status = format!("Suggestion #{index_1based} is empty");
            return;
        }
        let command = match parse_shell_line(trimmed) {
            Ok(command) => command,
            Err(err) => {
                self.agent_status = format!("Suggestion #{index_1based} parse error: {err}");
                self.agent_execution_log.push(AgentCommandExecutionRecord {
                    index_1based,
                    command: trimmed.to_string(),
                    trigger: trigger.to_string(),
                    ok: false,
                    state_changed: false,
                    summary: format!("parse error: {err}"),
                    executed_at_unix_ms: Self::now_unix_ms(),
                });
                return;
            }
        };
        if matches!(
            command,
            ShellCommand::AgentsAsk { .. }
                | ShellCommand::AgentsPlan { .. }
                | ShellCommand::AgentsExecutePlan { .. }
        ) {
            self.agent_status = format!(
                "Suggestion #{index_1based} rejected: agent-to-agent 'agents ...' commands are blocked"
            );
            self.agent_execution_log.push(AgentCommandExecutionRecord {
                index_1based,
                command: trimmed.to_string(),
                trigger: trigger.to_string(),
                ok: false,
                state_changed: false,
                summary: "agent-to-agent agents command blocked".to_string(),
                executed_at_unix_ms: Self::now_unix_ms(),
            });
            return;
        }
        if let Some(summary) = self.try_apply_shell_ui_intent(&command) {
            self.agent_status = format!("Suggestion #{index_1based}: {summary}");
            self.agent_execution_log.push(AgentCommandExecutionRecord {
                index_1based,
                command: trimmed.to_string(),
                trigger: trigger.to_string(),
                ok: true,
                state_changed: false,
                summary,
                executed_at_unix_ms: Self::now_unix_ms(),
            });
            if self.agent_execution_log.len() > 100 {
                let drain = self.agent_execution_log.len() - 100;
                self.agent_execution_log.drain(0..drain);
            }
            return;
        }
        let options = ShellExecutionOptions {
            allow_screenshots: false,
            allow_agent_commands: false,
            progress_callback: None,
        };
        let run = {
            let mut guard = self.engine.write().unwrap();
            execute_shell_command_with_options(&mut guard, &command, &options)
        };
        match run {
            Ok(run) => {
                if run.state_changed {
                    self.lineage_cache_valid = false;
                }
                let summary = if run.state_changed {
                    "executed (state changed)".to_string()
                } else {
                    "executed".to_string()
                };
                self.agent_status = format!("Suggestion #{index_1based}: {summary}");
                self.agent_execution_log.push(AgentCommandExecutionRecord {
                    index_1based,
                    command: trimmed.to_string(),
                    trigger: trigger.to_string(),
                    ok: true,
                    state_changed: run.state_changed,
                    summary,
                    executed_at_unix_ms: Self::now_unix_ms(),
                });
            }
            Err(err) => {
                self.agent_status = format!("Suggestion #{index_1based} failed: {err}");
                self.agent_execution_log.push(AgentCommandExecutionRecord {
                    index_1based,
                    command: trimmed.to_string(),
                    trigger: trigger.to_string(),
                    ok: false,
                    state_changed: false,
                    summary: err,
                    executed_at_unix_ms: Self::now_unix_ms(),
                });
            }
        }
        if self.agent_execution_log.len() > 100 {
            let drain = self.agent_execution_log.len() - 100;
            self.agent_execution_log.drain(0..drain);
        }
    }

    pub(super) fn try_apply_shell_ui_intent(&mut self, command: &ShellCommand) -> Option<String> {
        let ShellCommand::UiIntent {
            action,
            target,
            genome_id,
            helper_mode,
            catalog_path,
            cache_dir,
            filter,
            species,
            latest,
        } = command
        else {
            return None;
        };
        let mut selected_genome_id = genome_id
            .as_deref()
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .map(str::to_string);
        if matches!(target, UiIntentTarget::PreparedReferences) {
            self.apply_prepared_reference_intent_scope(*helper_mode, catalog_path, cache_dir);
            if selected_genome_id.is_none() {
                match self.resolve_prepared_reference_intent_selection(
                    *helper_mode,
                    catalog_path.clone(),
                    cache_dir.clone(),
                    filter.clone(),
                    species.clone(),
                    *latest,
                ) {
                    Ok(Some(resolved)) => {
                        selected_genome_id = Some(resolved);
                    }
                    Ok(None) => {}
                    Err(err) => {
                        self.app_status = format!(
                            "Could not resolve prepared-reference selection for ui intent: {err}"
                        );
                    }
                }
            }
        }
        if let Some(genome_id) = selected_genome_id
            .as_deref()
            .map(str::trim)
            .filter(|v| !v.is_empty())
        {
            self.genome_id = genome_id.to_string();
            self.invalidate_genome_genes();
        }
        match target {
            UiIntentTarget::PreparedReferences => self.open_reference_genome_inspector_dialog(),
            UiIntentTarget::PrepareReferenceGenome => self.open_reference_genome_prepare_dialog(),
            UiIntentTarget::RetrieveGenomeSequence => self.open_reference_genome_retrieve_dialog(),
            UiIntentTarget::BlastGenomeSequence => self.open_reference_genome_blast_dialog(),
            UiIntentTarget::ImportGenomeTrack => self.open_genome_bed_track_dialog(),
            UiIntentTarget::PcrDesign => self.open_pcr_design_dialog(),
            UiIntentTarget::SequencingConfirmation => self.open_sequencing_confirmation_dialog(),
            UiIntentTarget::AgentAssistant => self.open_agent_assistant_dialog(),
            UiIntentTarget::PrepareHelperGenome => self.open_helper_genome_prepare_dialog(),
            UiIntentTarget::RetrieveHelperSequence => self.open_helper_genome_retrieve_dialog(),
            UiIntentTarget::BlastHelperSequence => self.open_helper_genome_blast_dialog(),
        }
        let mut summary = format!("ui intent {} '{}'", action.as_str(), target.as_str());
        if let Some(genome_id) = selected_genome_id {
            summary.push_str(&format!(" (selected_genome_id={genome_id})"));
        }
        Some(summary)
    }

    pub(super) fn apply_prepared_reference_intent_scope(
        &mut self,
        helper_mode: bool,
        catalog_path: &Option<String>,
        cache_dir: &Option<String>,
    ) {
        let scope = if helper_mode {
            GenomeDialogScope::Helper
        } else {
            GenomeDialogScope::Reference
        };
        let normalized_catalog = catalog_path
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(str::to_string);
        let normalized_cache = cache_dir
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(str::to_string);
        self.sync_active_genome_scope_paths_from_fields();
        self.genome_dialog_scope = scope;
        let (scope_catalog, scope_cache) = self.scope_genome_paths_resolved(scope);
        let next_catalog = normalized_catalog.unwrap_or(scope_catalog);
        let next_cache = normalized_cache.unwrap_or(scope_cache);
        let catalog_changed = self.genome_catalog_path != next_catalog;
        let cache_changed = self.genome_cache_dir != next_cache;
        self.genome_catalog_path = next_catalog.clone();
        self.genome_cache_dir = next_cache.clone();
        self.set_scope_genome_paths(scope, next_catalog, next_cache);
        if catalog_changed || cache_changed {
            self.invalidate_genome_genes();
        }
    }

    pub(super) fn resolve_prepared_reference_intent_selection(
        &self,
        helper_mode: bool,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
        filter: Option<String>,
        species: Option<String>,
        latest: bool,
    ) -> Result<Option<String>, String> {
        let mut engine = self.engine.write().unwrap();
        let run = execute_shell_command_with_options(
            &mut engine,
            &ShellCommand::UiPreparedGenomes {
                helper_mode,
                catalog_path,
                cache_dir,
                filter,
                species,
                latest,
            },
            &ShellExecutionOptions::default(),
        )?;
        Ok(run
            .output
            .get("selected_genome_id")
            .and_then(|value| value.as_str())
            .map(str::to_string))
    }

    pub(super) fn execute_agent_auto_suggestions(&mut self, response: &AgentResponse) {
        for (idx, suggestion) in response.suggested_commands.iter().enumerate() {
            if suggestion.execution == AgentExecutionIntent::Auto {
                self.execute_agent_suggested_command(idx + 1, &suggestion.command, "auto");
            }
        }
    }

    pub(super) fn poll_agent_assistant_task(&mut self, ctx: &egui::Context) {
        if self.agent_task.is_none() {
            return;
        }
        ctx.request_repaint_after(Duration::from_millis(100));
        let mut done: Option<(u64, Result<AgentInvocationOutcome, String>)> = None;
        if let Some(task) = &self.agent_task {
            match task.receiver.try_recv() {
                Ok(AgentAskTaskMessage::Done { job_id, result }) => {
                    done = Some((job_id, result));
                }
                Err(mpsc::TryRecvError::Empty) => {}
                Err(mpsc::TryRecvError::Disconnected) => {
                    done = Some((task.job_id, Err("Agent worker disconnected".to_string())));
                }
            }
        }
        if let Some((job_id, outcome)) = done {
            let elapsed = self
                .agent_task
                .as_ref()
                .map(|task| task.started.elapsed().as_secs_f64())
                .unwrap_or(0.0);
            self.agent_task = None;
            match outcome {
                Ok(invocation) => {
                    let suggestion_count = invocation.response.suggested_commands.len();
                    self.agent_status = format!(
                        "Agent response received in {:.1}s (suggestions={})",
                        elapsed, suggestion_count
                    );
                    self.push_job_event(
                        BackgroundJobKind::AgentAssist,
                        BackgroundJobEventPhase::Completed,
                        Some(job_id),
                        format!(
                            "Agent '{}' completed in {:.1}s (suggestions={})",
                            invocation.system_id, elapsed, suggestion_count
                        ),
                    );
                    let response = invocation.response.clone();
                    self.agent_last_invocation = Some(invocation);
                    if self.agent_allow_auto_exec {
                        self.execute_agent_auto_suggestions(&response);
                    }
                }
                Err(err) => {
                    self.agent_status =
                        format!("Agent request failed after {:.1}s: {}", elapsed, err);
                    self.push_job_event(
                        BackgroundJobKind::AgentAssist,
                        BackgroundJobEventPhase::Failed,
                        Some(job_id),
                        format!("Agent request failed in {:.1}s: {}", elapsed, err),
                    );
                }
            }
        }
    }

    pub(super) fn poll_agent_model_discovery_task(&mut self, ctx: &egui::Context) {
        if self.agent_model_discovery_task.is_none() {
            return;
        }
        ctx.request_repaint_after(Duration::from_millis(100));
        let mut done: Option<(String, Result<Vec<String>, String>)> = None;
        if let Some(task) = &self.agent_model_discovery_task {
            match task.receiver.try_recv() {
                Ok(AgentModelDiscoveryTaskMessage::Done { source_key, result }) => {
                    done = Some((source_key, result));
                }
                Err(mpsc::TryRecvError::Empty) => {}
                Err(mpsc::TryRecvError::Disconnected) => {
                    done = Some((
                        task.source_key.clone(),
                        Err("Model discovery worker disconnected".to_string()),
                    ));
                }
            }
        }
        if let Some((source_key, result)) = done {
            let elapsed = self
                .agent_model_discovery_task
                .as_ref()
                .map(|task| task.started.elapsed().as_secs_f64())
                .unwrap_or(0.0);
            self.agent_model_discovery_task = None;
            if source_key != self.agent_model_discovery_source_key {
                return;
            }
            match result {
                Ok(models) => {
                    self.agent_model_discovery_failed_source_key.clear();
                    self.agent_discovered_models = models;
                    if self.agent_discovered_models.is_empty() {
                        self.agent_model_discovery_status =
                            format!("Model discovery returned no models ({:.1}s)", elapsed);
                        self.agent_discovered_model_pick.clear();
                    } else {
                        if !self
                            .agent_discovered_models
                            .iter()
                            .any(|item| item == &self.agent_discovered_model_pick)
                        {
                            self.agent_discovered_model_pick.clear();
                        }
                        self.agent_model_discovery_status = format!(
                            "Discovered {} model(s) in {:.1}s",
                            self.agent_discovered_models.len(),
                            elapsed
                        );
                    }
                }
                Err(err) => {
                    self.agent_discovered_models.clear();
                    self.agent_discovered_model_pick.clear();
                    self.agent_model_discovery_failed_source_key = source_key;
                    let hint = Self::agent_model_discovery_failure_hint(&err)
                        .map(|hint| format!(" {hint}"))
                        .unwrap_or_default();
                    self.agent_model_discovery_status = format!(
                        "Model discovery failed after {:.1}s: {}{}",
                        elapsed, err, hint
                    );
                }
            }
        }
    }
    pub(super) fn render_routine_assistant_contents(&mut self, ui: &mut Ui) -> bool {
        let mut close_requested = false;
        let close_hover = Self::specialist_window_close_hover_text("Routine Assistant");
        if self.render_specialist_window_nav_with_close(ui, Some(("Close", close_hover.as_str()))) {
            close_requested = true;
        }
        ui.label(
            "Apply cloning routines through one staged flow driven by shared engine commands.",
        );
        ui.small(
            "Flow: goal -> candidate routines -> compare alternatives -> parameter bindings -> preflight -> transactional run -> run-bundle export",
        );
        ui.separator();

        let selected_routine = self.routine_assistant_selected_routine();
        let has_selected = selected_routine.is_some();
        let has_preflight = self.routine_assistant_preflight_output.is_some();
        let has_execute = self.routine_assistant_execute_output.is_some();

        let stage_order = [
            RoutineAssistantStage::GoalAndCandidates,
            RoutineAssistantStage::Compare,
            RoutineAssistantStage::Parameters,
            RoutineAssistantStage::Preflight,
            RoutineAssistantStage::ExecuteAndExport,
        ];
        ui.horizontal_wrapped(|ui| {
            for stage in stage_order {
                let enabled = match stage {
                    RoutineAssistantStage::GoalAndCandidates => true,
                    RoutineAssistantStage::Compare => has_selected,
                    RoutineAssistantStage::Parameters => has_selected,
                    RoutineAssistantStage::Preflight => has_preflight,
                    RoutineAssistantStage::ExecuteAndExport => has_execute,
                };
                let resp = ui.add_enabled(
                    enabled,
                    egui::Button::new(stage.label())
                        .selected(self.routine_assistant_stage == stage),
                );
                if resp.clicked() {
                    self.routine_assistant_stage = stage;
                }
            }
        });
        if !self.routine_assistant_status.trim().is_empty() {
            ui.separator();
            ui.monospace(self.routine_assistant_status.trim());
        }
        ui.separator();

        match self.routine_assistant_stage {
            RoutineAssistantStage::GoalAndCandidates => {
                ui.horizontal(|ui| {
                    ui.label("goal");
                    ui.text_edit_singleline(&mut self.routine_assistant_goal);
                });
                ui.horizontal(|ui| {
                    ui.label("query");
                    ui.text_edit_singleline(&mut self.routine_assistant_query);
                    if ui
                        .button("Find Candidates")
                        .on_hover_text(
                            "Query routine catalog by goal/query text and load candidate routines",
                        )
                        .clicked()
                    {
                        self.refresh_routine_assistant_candidates();
                    }
                    if ui
                        .button("Reset")
                        .on_hover_text("Clear selected routine and staged assistant outputs")
                        .clicked()
                    {
                        self.maybe_mark_routine_assistant_trace_aborted();
                        self.routine_assistant_selected_routine_id.clear();
                        self.routine_assistant_compare_routine_id.clear();
                        self.routine_assistant_bindings.clear();
                        self.routine_assistant_disambiguation_answers.clear();
                        self.routine_assistant_explain_output = None;
                        self.routine_assistant_compare_output = None;
                        self.routine_assistant_preflight_output = None;
                        self.routine_assistant_execute_output = None;
                        self.routine_assistant_stage = RoutineAssistantStage::GoalAndCandidates;
                        self.routine_assistant_status =
                            "Routine Assistant: reset staged state".to_string();
                        self.routine_assistant_decision_trace = None;
                        self.ensure_routine_assistant_decision_trace_started();
                    }
                });
                ui.separator();
                self.render_routine_assistant_planning_context_strip(ui);
                if self.routine_assistant_preference_context.is_some() {
                    ui.separator();
                }
                if self.routine_assistant_candidates.is_empty() {
                    ui.small("No routine candidates loaded. Use 'Find Candidates'.");
                } else {
                    let mut choose_routine: Option<String> = None;
                    egui::ScrollArea::vertical()
                        .max_height(360.0)
                        .show(ui, |ui| {
                            for routine in &self.routine_assistant_candidates {
                                ui.group(|ui| {
                                    ui.horizontal(|ui| {
                                        ui.strong(format!(
                                            "{} ({})",
                                            routine.title, routine.routine_id
                                        ));
                                        ui.label(format!(
                                            "[family: {}, status: {}]",
                                            routine.family, routine.status
                                        ));
                                        if ui
                                            .button("Select")
                                            .on_hover_text(
                                                "Choose this routine as primary candidate and move to alternative comparison",
                                            )
                                            .clicked()
                                        {
                                            choose_routine = Some(routine.routine_id.clone());
                                        }
                                    });
                                    if let Some(summary) = routine.summary.as_deref() {
                                        ui.small(summary);
                                    }
                                    if let Some(purpose) = routine.purpose.as_deref() {
                                        ui.small(format!("purpose: {purpose}"));
                                    }
                                    if let Some(score) = routine.composite_meta_score {
                                        ui.small(format!(
                                            "planning score: {:.3} | fit: {:.3} | time: {:.2} h | cost: {:.2}",
                                            score,
                                            routine.local_fit_score.unwrap_or_default(),
                                            routine.estimated_time_hours.unwrap_or_default(),
                                            routine.estimated_cost.unwrap_or_default()
                                        ));
                                    }
                                    if let Some(estimate) = routine.planning_estimate.as_ref() {
                                        let bonus = estimate
                                            .explanation
                                            .get("routine_family_alignment_bonus")
                                            .and_then(|value| value.as_f64())
                                            .unwrap_or(0.0);
                                        if bonus > 0.0 {
                                            let sources = estimate
                                                .explanation
                                                .get("routine_family_alignment_sources")
                                                .and_then(|value| value.as_array())
                                                .map(|rows| {
                                                    rows.iter()
                                                        .filter_map(|row| row.as_str())
                                                        .collect::<Vec<_>>()
                                                        .join(", ")
                                                })
                                                .unwrap_or_default();
                                            ui.small(format!(
                                                "family-alignment bonus: +{:.2}{}",
                                                bonus,
                                                if sources.is_empty() {
                                                    String::new()
                                                } else {
                                                    format!(" ({sources})")
                                                }
                                            ));
                                        }
                                    }
                                });
                            }
                        });
                    if let Some(routine_id) = choose_routine {
                        self.routine_assistant_selected_routine_id = routine_id;
                        self.routine_assistant_compare_routine_id.clear();
                        self.routine_assistant_disambiguation_answers.clear();
                        self.routine_assistant_explain_output = None;
                        self.routine_assistant_compare_output = None;
                        self.routine_assistant_preflight_output = None;
                        self.routine_assistant_execute_output = None;
                        self.sync_routine_assistant_bindings_for_selected();
                        self.load_routine_assistant_explain();
                        self.routine_assistant_stage = RoutineAssistantStage::Compare;
                    }
                }
            }
            RoutineAssistantStage::Compare => {
                let Some(routine) = selected_routine else {
                    ui.small("Select a primary routine first in stage 1.");
                    return close_requested;
                };
                ui.strong(format!(
                    "Primary routine: {} ({})",
                    routine.title, routine.routine_id
                ));
                if let Some(summary) = routine.summary.as_deref() {
                    ui.small(summary);
                }
                if let Some(planning) = self
                    .routine_assistant_explain_output
                    .as_ref()
                    .and_then(|value| value.get("planning"))
                {
                    if let Some(estimate) = planning.get("estimate") {
                        let composite = estimate
                            .get("composite_meta_score")
                            .and_then(|value| value.as_f64());
                        let local_fit = estimate
                            .get("local_fit_score")
                            .and_then(|value| value.as_f64());
                        let time_hours = estimate
                            .get("estimated_time_hours")
                            .and_then(|value| value.as_f64());
                        let cost = estimate
                            .get("estimated_cost")
                            .and_then(|value| value.as_f64());
                        if composite.is_some()
                            || local_fit.is_some()
                            || time_hours.is_some()
                            || cost.is_some()
                        {
                            ui.small(format!(
                                "sequence-aware planning: score {} | fit {} | time {} h | cost {}",
                                composite
                                    .map(|value| format!("{value:.3}"))
                                    .unwrap_or_else(|| "-".to_string()),
                                local_fit
                                    .map(|value| format!("{value:.3}"))
                                    .unwrap_or_else(|| "-".to_string()),
                                time_hours
                                    .map(|value| format!("{value:.2}"))
                                    .unwrap_or_else(|| "-".to_string()),
                                cost.map(|value| format!("{value:.2}"))
                                    .unwrap_or_else(|| "-".to_string())
                            ));
                        }
                        let bonus = estimate
                            .get("explanation")
                            .and_then(|value| value.get("routine_family_alignment_bonus"))
                            .and_then(|value| value.as_f64())
                            .unwrap_or(0.0);
                        if bonus > 0.0 {
                            let sources = estimate
                                .get("explanation")
                                .and_then(|value| value.get("routine_family_alignment_sources"))
                                .and_then(|value| value.as_array())
                                .map(|rows| {
                                    rows.iter()
                                        .filter_map(|row| row.as_str())
                                        .collect::<Vec<_>>()
                                        .join(", ")
                                })
                                .unwrap_or_default();
                            ui.small(format!(
                                "alignment bonus: +{bonus:.2}{}",
                                if sources.is_empty() {
                                    String::new()
                                } else {
                                    format!(" ({sources})")
                                }
                            ));
                        }
                    }
                }
                self.render_routine_assistant_macro_suggestions(ui);
                ui.horizontal(|ui| {
                    if ui
                        .button("Reload Explanation")
                        .on_hover_text(
                            "Fetch routine explainability payload from shared routines explain command",
                        )
                        .clicked()
                    {
                        self.load_routine_assistant_explain();
                    }
                    if ui
                        .button("Continue to Parameters")
                        .on_hover_text("Proceed with typed parameter binding form")
                        .clicked()
                    {
                        self.routine_assistant_stage = RoutineAssistantStage::Parameters;
                    }
                });
                ui.separator();
                let alternatives = self
                    .routine_assistant_explain_output
                    .as_ref()
                    .and_then(|value| value.get("alternatives"))
                    .and_then(|value| value.as_array())
                    .cloned()
                    .unwrap_or_default();
                if alternatives.is_empty() {
                    ui.small("No explicit alternatives listed for this routine.");
                } else {
                    ui.horizontal(|ui| {
                        ui.label("compare against");
                        egui::ComboBox::from_id_salt("routine_assistant_compare_combo")
                            .selected_text(if self.routine_assistant_compare_routine_id.is_empty() {
                                "(select alternative)"
                            } else {
                                self.routine_assistant_compare_routine_id.as_str()
                            })
                            .show_ui(ui, |ui| {
                                for row in &alternatives {
                                    if let Some(routine_id) =
                                        row.get("routine_id").and_then(|value| value.as_str())
                                    {
                                        let title = row
                                            .get("title")
                                            .and_then(|value| value.as_str())
                                            .unwrap_or(routine_id);
                                        let label = format!("{title} ({routine_id})");
                                        if ui
                                            .selectable_label(
                                                self.routine_assistant_compare_routine_id
                                                    .eq_ignore_ascii_case(routine_id),
                                                label,
                                            )
                                            .clicked()
                                        {
                                            self.routine_assistant_compare_routine_id =
                                                routine_id.to_string();
                                        }
                                    }
                                }
                            });
                        if ui
                            .button("Compare")
                            .on_hover_text(
                                "Run shared routines compare command and display deterministic difference matrix",
                            )
                            .clicked()
                        {
                            self.load_routine_assistant_compare();
                        }
                    });
                }
                if let Some(compare) = &self.routine_assistant_compare_output {
                    ui.separator();
                    if let Some(rows) = compare
                        .get("comparison")
                        .and_then(|value| value.get("difference_matrix"))
                        .and_then(|value| value.as_array())
                    {
                        ui.strong("Difference matrix");
                        egui::Grid::new("routine_assistant_compare_grid")
                            .striped(true)
                            .show(ui, |ui| {
                                ui.strong("axis");
                                ui.strong("primary");
                                ui.strong("alternative");
                                ui.end_row();
                                for row in rows {
                                    let axis = row
                                        .get("axis")
                                        .and_then(|value| value.as_str())
                                        .unwrap_or("-");
                                    let left = row
                                        .get("left")
                                        .and_then(|value| value.as_str())
                                        .unwrap_or("-");
                                    let right = row
                                        .get("right")
                                        .and_then(|value| value.as_str())
                                        .unwrap_or("-");
                                    ui.monospace(axis);
                                    ui.label(left);
                                    ui.label(right);
                                    ui.end_row();
                                }
                            });
                    }
                }
                let disambiguation_questions =
                    self.routine_assistant_effective_disambiguation_questions();
                if disambiguation_questions.is_empty() {
                    ui.separator();
                    ui.small("No disambiguation questions provided for this routine pair.");
                } else {
                    self.sync_routine_assistant_disambiguation_answers_for_questions(
                        &disambiguation_questions,
                        &[],
                    );
                    ui.separator();
                    ui.strong("Disambiguation answers");
                    let mut answers_changed = false;
                    egui::Grid::new("routine_assistant_disambiguation_answers_grid")
                        .num_columns(2)
                        .striped(true)
                        .show(ui, |ui| {
                            ui.strong("question");
                            ui.strong("answer");
                            ui.end_row();
                            for row in &disambiguation_questions {
                                ui.label(row.question_text.as_str());
                                let mut answer = self
                                    .routine_assistant_disambiguation_answers
                                    .get(&row.question_id)
                                    .cloned()
                                    .unwrap_or_default();
                                let answer_resp = ui.text_edit_singleline(&mut answer);
                                if answer_resp.changed() {
                                    self.routine_assistant_disambiguation_answers
                                        .insert(row.question_id.clone(), answer);
                                    answers_changed = true;
                                }
                                ui.end_row();
                            }
                        });
                    if answers_changed {
                        let selected = self.routine_assistant_selected_routine();
                        let disambiguation_answers = self
                            .routine_assistant_disambiguation_answers_snapshot(
                                &disambiguation_questions,
                            );
                        self.update_routine_assistant_decision_trace(|trace| {
                            trace.status = "draft".to_string();
                            Self::routine_assistant_capture_selected_routine(
                                trace,
                                selected.as_ref(),
                            );
                            Self::merge_routine_assistant_disambiguation_questions(
                                &mut trace.disambiguation_questions_presented,
                                disambiguation_questions.clone(),
                            );
                            trace.disambiguation_answers = disambiguation_answers;
                        });
                    }
                }
            }
            RoutineAssistantStage::Parameters => {
                let Some(routine) = selected_routine else {
                    ui.small("Select a primary routine first in stage 1.");
                    return close_requested;
                };
                self.sync_routine_assistant_bindings_for_selected();
                ui.strong(format!(
                    "Template: {} (routine: {})",
                    routine.template_name, routine.routine_id
                ));
                if !routine.requires.is_empty() {
                    ui.label("requires");
                    for req in &routine.requires {
                        ui.small(format!("- {req}"));
                    }
                }
                self.render_routine_assistant_macro_suggestions(ui);
                self.render_routine_assistant_gibson_linearization_notice(ui, &routine);
                ui.separator();
                egui::Grid::new("routine_assistant_bindings_grid")
                    .striped(true)
                    .show(ui, |ui| {
                        ui.strong("input");
                        ui.strong("value");
                        ui.strong("type");
                        ui.end_row();
                        for port in &routine.input_ports {
                            let port_id = port
                                .get("port_id")
                                .and_then(|value| value.as_str())
                                .map(str::trim)
                                .unwrap_or("");
                            if port_id.is_empty() {
                                continue;
                            }
                            let kind = port
                                .get("kind")
                                .and_then(|value| value.as_str())
                                .unwrap_or("-");
                            let required = port
                                .get("required")
                                .and_then(|value| value.as_bool())
                                .unwrap_or(false);
                            let description = port
                                .get("description")
                                .and_then(|value| value.as_str())
                                .unwrap_or("");
                            let label = if required {
                                format!("{port_id} *")
                            } else {
                                port_id.to_string()
                            };
                            ui.label(label).on_hover_text(description);
                            let mut entry = self
                                .routine_assistant_bindings
                                .get(port_id)
                                .cloned()
                                .unwrap_or_default();
                            let edit_resp = ui.text_edit_singleline(&mut entry);
                            if edit_resp.changed() {
                                self.routine_assistant_bindings
                                    .insert(port_id.to_string(), entry.clone());
                                self.routine_assistant_preflight_output = None;
                                self.routine_assistant_execute_output = None;
                                let selected = self.routine_assistant_selected_routine();
                                let bindings_snapshot = self.routine_assistant_bindings_snapshot();
                                self.update_routine_assistant_decision_trace(|trace| {
                                    trace.status = "draft".to_string();
                                    trace.bindings_snapshot = bindings_snapshot;
                                    Self::routine_assistant_capture_selected_routine(
                                        trace,
                                        selected.as_ref(),
                                    );
                                });
                            }
                            if kind.eq_ignore_ascii_case("sequence") {
                                let compact = entry.trim();
                                let (type_text, text_color) = if compact.is_empty() {
                                    ("sequence".to_string(), egui::Color32::GRAY)
                                } else if let Some((circular, length_bp)) =
                                    self.routine_assistant_sequence_topology_for_seq_id(compact)
                                {
                                    if circular {
                                        (
                                            format!("sequence | circular | {length_bp} bp"),
                                            egui::Color32::from_rgb(190, 70, 70),
                                        )
                                    } else {
                                        (
                                            format!("sequence | linear | {length_bp} bp"),
                                            egui::Color32::from_rgb(70, 130, 80),
                                        )
                                    }
                                } else {
                                    (
                                        "sequence | missing".to_string(),
                                        egui::Color32::from_rgb(190, 70, 70),
                                    )
                                };
                                ui.label(
                                    egui::RichText::new(type_text).monospace().color(text_color),
                                );
                            } else {
                                ui.monospace(kind);
                            }
                            ui.end_row();
                        }
                    });
                ui.separator();
                ui.horizontal(|ui| {
                    if ui
                        .button("Back to Compare")
                        .on_hover_text("Return to routine alternative comparison")
                        .clicked()
                    {
                        self.routine_assistant_stage = RoutineAssistantStage::Compare;
                    }
                    if ui
                        .button("Run Preflight")
                        .on_hover_text(
                            "Run macros template-run --validate-only through shared shell executor",
                        )
                        .clicked()
                    {
                        self.run_routine_assistant_preflight();
                    }
                });
            }
            RoutineAssistantStage::Preflight => {
                if let Some(routine) = selected_routine {
                    self.render_routine_assistant_gibson_linearization_notice(ui, &routine);
                    ui.separator();
                }
                if let Some(output) = &self.routine_assistant_preflight_output {
                    let can_execute = output
                        .get("can_execute")
                        .and_then(|value| value.as_bool())
                        .unwrap_or(false);
                    ui.strong(format!(
                        "Preflight status: {}",
                        if can_execute {
                            "can execute"
                        } else {
                            "blocking errors"
                        }
                    ));
                    if let Some(preflight) = output.get("preflight") {
                        if let Some(errors) =
                            preflight.get("errors").and_then(|value| value.as_array())
                        {
                            if !errors.is_empty() {
                                ui.label("errors");
                                for err in errors {
                                    if let Some(text) = err.as_str() {
                                        ui.colored_label(
                                            egui::Color32::from_rgb(190, 70, 70),
                                            format!("- {text}"),
                                        );
                                    }
                                }
                            }
                        }
                        if let Some(warnings) =
                            preflight.get("warnings").and_then(|value| value.as_array())
                        {
                            if !warnings.is_empty() {
                                ui.label("warnings");
                                for warning in warnings {
                                    if let Some(text) = warning.as_str() {
                                        ui.small(format!("- {text}"));
                                    }
                                }
                            }
                        }
                    }
                } else {
                    ui.small("No preflight output yet. Run preflight in stage 3.");
                }
                ui.separator();
                ui.horizontal(|ui| {
                    if ui
                        .button("Back to Parameters")
                        .on_hover_text("Edit typed bindings before running again")
                        .clicked()
                    {
                        self.routine_assistant_stage = RoutineAssistantStage::Parameters;
                    }
                    let exec_resp = ui.add_enabled(
                        self.routine_assistant_can_execute(),
                        egui::Button::new("Run Transactional"),
                    );
                    if exec_resp
                        .on_hover_text(
                            "Execute macros template-run --transactional using current bindings",
                        )
                        .clicked()
                    {
                        self.run_routine_assistant_execute();
                    }
                });
            }
            RoutineAssistantStage::ExecuteAndExport => {
                if let Some(output) = &self.routine_assistant_execute_output {
                    let macro_instance = output
                        .get("macro_instance_id")
                        .and_then(|value| value.as_str())
                        .unwrap_or("-");
                    ui.strong("Transactional run completed");
                    ui.monospace(format!("macro_instance_id: {macro_instance}"));
                    if let Some(run) = output.get("run") {
                        let created = run
                            .get("created")
                            .and_then(|value| value.as_array())
                            .map(|rows| rows.len())
                            .unwrap_or(0);
                        let changed = run
                            .get("changed")
                            .and_then(|value| value.as_array())
                            .map(|rows| rows.len())
                            .unwrap_or(0);
                        ui.small(format!("created: {created}, changed: {changed}"));
                    }
                } else {
                    ui.small("No transactional run output yet.");
                }
                ui.separator();
                ui.horizontal(|ui| {
                    if ui
                        .button("Back to Preflight")
                        .on_hover_text("Inspect or re-run preflight checks")
                        .clicked()
                    {
                        self.routine_assistant_stage = RoutineAssistantStage::Preflight;
                    }
                    if ui
                        .button("Export Run Bundle")
                        .on_hover_text(
                            "Export deterministic process run bundle (inputs, parameter changes, operation log, outputs)",
                        )
                        .clicked()
                    {
                        self.export_routine_assistant_run_bundle();
                    }
                });
            }
        }
        close_requested
    }

    pub(super) fn render_routine_assistant_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_routine_assistant_dialog {
            return;
        }
        let was_open = self.show_routine_assistant_dialog;
        let mut open = self.show_routine_assistant_dialog;
        let viewport_id = Self::routine_assistant_viewport_id();
        let spec = self.hosted_window_spec_for_viewport(
            "Routine Assistant",
            Self::hosted_routine_assistant_window_id(),
            viewport_id,
            Vec2::new(980.0, 720.0),
            Vec2::new(720.0, 480.0),
        );
        if ctx.embed_viewports() {
            let mut close_requested = false;
            crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
                close_requested = self.render_routine_assistant_contents(ui);
            });
            if close_requested {
                open = false;
            }
            if ctx.input(|i| i.key_pressed(Key::Escape)) {
                open = false;
            }
            if was_open && !open {
                self.maybe_mark_routine_assistant_trace_aborted();
            }
            self.show_routine_assistant_dialog = open;
            self.finalize_viewport_open_probe(viewport_id, "Routine Assistant");
            return;
        }
        let builder = crate::egui_compat::viewport_builder_for_hosted_window(&spec);
        ctx.show_viewport_immediate(viewport_id, builder, |ctx, class| {
            self.note_viewport_focus_if_active(ctx, viewport_id);
            if class == egui::ViewportClass::EmbeddedWindow {
                let mut close_requested = false;
                crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
                    close_requested = self.render_routine_assistant_contents(ui);
                });
                if close_requested {
                    open = false;
                }
            } else {
                let mut close_requested = false;
                crate::egui_compat::show_central_panel(ctx, egui::CentralPanel::default(), |ui| {
                    close_requested = self.render_routine_assistant_contents(ui);
                });
                if close_requested || Self::viewport_close_requested_or_shortcut(ctx) {
                    open = false;
                }
            }
        });
        if ctx.input(|i| i.key_pressed(Key::Escape)) {
            open = false;
        }
        if was_open && !open {
            self.maybe_mark_routine_assistant_trace_aborted();
        }
        self.show_routine_assistant_dialog = open;
    }

    pub(super) fn render_agent_assistant_contents(&mut self, ui: &mut Ui) -> bool {
        self.refresh_agent_system_catalog();
        let mut close_requested = false;
        let close_hover = Self::specialist_window_close_hover_text("Agent Assistant");
        if self.render_specialist_window_nav_with_close(ui, Some(("Close", close_hover.as_str()))) {
            close_requested = true;
        }
        ui.label("Ask an agent system for project support, then execute suggested GENtle shell commands per reply.");
        ui.horizontal(|ui| {
            ui.label("catalog");
            ui.text_edit_singleline(&mut self.agent_catalog_path);
            if ui
                .button("Browse...")
                .on_hover_text("Browse filesystem and fill this path")
                .clicked()
            {
                if let Some(path) = rfd::FileDialog::new()
                    .add_filter("JSON", &["json"])
                    .pick_file()
                {
                    self.agent_catalog_path = path.display().to_string();
                    self.agent_catalog_loaded_path.clear();
                    self.refresh_agent_system_catalog();
                }
            }
        });
        if !self.agent_catalog_error.is_empty() {
            ui.colored_label(
                egui::Color32::from_rgb(190, 70, 70),
                format!("Catalog error: {}", self.agent_catalog_error),
            );
        }
        let mut preflight_inputs_changed = false;
        let mut requested_agent_system_id: Option<String> = None;
        ui.horizontal(|ui| {
            ui.label("system");
            egui::ComboBox::from_id_salt("agent_system_combo")
                .selected_text(if self.agent_system_id.trim().is_empty() {
                    "(choose system)"
                } else {
                    self.agent_system_id.as_str()
                })
                .show_ui(ui, |ui| {
                    for system in &self.agent_systems {
                        let (available, reason) = self.selected_agent_system_availability(system);
                        let label = if available {
                            format!("{} ({})", system.label, system.id)
                        } else {
                            format!("{} ({}) [unavailable]", system.label, system.id)
                        };
                        let mut response = ui.add(
                            egui::Button::new(label).selected(self.agent_system_id == system.id),
                        );
                        if !available {
                            response = response.on_hover_text(
                                reason.unwrap_or_else(|| "agent system unavailable".to_string()),
                            );
                        }
                        if response.clicked() {
                            requested_agent_system_id = Some(system.id.clone());
                        }
                    }
                });
        });
        if let Some(system_id) = requested_agent_system_id {
            self.select_agent_system_and_reset_setup(&system_id);
        }
        if !self.agent_systems.is_empty() {
            ui.group(|ui| {
                ui.strong("Quick start");
                ui.small(
                    "Choose whether GENtle should talk to OpenAI, Claude, a local OpenAI-compatible model, or the offline demo.",
                );
                ui.horizontal_wrapped(|ui| {
                    if let Some(openai_system_id) =
                        preferred_openai_agent_system_id(&self.agent_systems)
                    {
                        if ui
                            .button("Use OpenAI API")
                            .on_hover_text(
                                "Select the native OpenAI agent profile and use OPENAI_API_KEY for requests",
                            )
                            .clicked()
                        {
                            self.select_agent_system_and_reset_setup(&openai_system_id);
                            self.agent_base_url_override.clear();
                            self.agent_model_override.clear();
                            self.agent_discovered_model_pick.clear();
                            self.agent_status = "Selected OpenAI API quick start. Add OPENAI_API_KEY or paste a session key, then run Test Setup.".to_string();
                        }
                    }
                    if let Some(anthropic_system_id) =
                        preferred_anthropic_agent_system_id(&self.agent_systems)
                    {
                        if ui
                            .button("Use Claude API")
                            .on_hover_text(
                                "Select the native Anthropic Claude profile and use ANTHROPIC_API_KEY for requests",
                            )
                            .clicked()
                        {
                            self.select_agent_system_and_reset_setup(&anthropic_system_id);
                            self.agent_base_url_override.clear();
                            self.agent_model_override.clear();
                            self.agent_discovered_model_pick.clear();
                            self.agent_status = "Selected Claude API quick start. Add ANTHROPIC_API_KEY or paste an Anthropic API key, then run Test Setup.".to_string();
                        }
                    }
                    if let Some(local_system_id) =
                        preferred_local_agent_system_id(&self.agent_systems)
                    {
                        if ui
                            .button("Use Local Model (no OpenAI API billing)")
                            .on_hover_text(
                                "Select a local OpenAI-compatible endpoint such as Ollama, Jan, or Msty",
                            )
                            .clicked()
                        {
                            self.select_agent_system_and_reset_setup(&local_system_id);
                            self.agent_base_url_override.clear();
                            self.agent_model_override.clear();
                            self.agent_discovered_model_pick.clear();
                            self.agent_status = "Selected local-model quick start. Start your local endpoint, discover models, then run Test Setup.".to_string();
                        }
                    }
                    if self.agent_systems.iter().any(|system| system.id == "builtin_echo")
                        && ui
                            .button("Use Demo Echo")
                            .on_hover_text(
                                "Select the offline demo assistant that never contacts a remote service",
                            )
                            .clicked()
                    {
                        self.select_agent_system_and_reset_setup("builtin_echo");
                        self.agent_status = "Selected the built-in demo assistant.".to_string();
                    }
                });
                ui.small(
                    "Cloud API modes use provider API keys and talk to the provider directly. If you want a no-extra-OpenAI-bill path inside GENtle, prefer a local OpenAI-compatible endpoint.",
                );
            });
            ui.group(|ui| {
                ui.strong("External Agent / MCP");
                ui.small(
                    "Use this route when your external agent already supports MCP and has its own account/runtime.",
                );
                ui.small(
                    "ChatGPT/Codex subscriptions are not OpenAI API keys; MCP exposes GENtle tools over stdio instead of making GENtle call the OpenAI API.",
                );
                let state_path = self.external_agent_mcp_state_path();
                ui.small(format!(
                    "state path: {}{}",
                    state_path,
                    if self.current_project_path.is_some() {
                        " (active saved project)"
                    } else {
                        " (default gentle_mcp state path)"
                    }
                ));
                ui.monospace(self.external_agent_mcp_command_snippet());
            });
        }
        let mut selected_available = false;
        if let Some(system) = self.selected_agent_system() {
            let (available, reason) = self.selected_agent_system_availability(&system);
            selected_available = available;
            if let Some(description) = system.description.as_deref() {
                let trimmed = description.trim();
                if !trimmed.is_empty() {
                    ui.small(trimmed);
                }
            }
            ui.small(format!("transport: {}", system.transport.as_str()));
            if !system.command.is_empty() {
                ui.small(format!("command: {}", system.command.join(" ")));
            }
            if matches!(
                system.transport,
                AgentSystemTransport::NativeOpenai
                    | AgentSystemTransport::NativeAnthropic
                    | AgentSystemTransport::NativeOpenaiCompat
            ) {
                if let Some(source_key) = self.selected_agent_model_discovery_source_key(&system) {
                    if self.agent_model_discovery_source_key != source_key {
                        self.agent_model_discovery_source_key = source_key;
                        self.clear_agent_model_discovery_snapshot();
                        self.clear_agent_preflight_output();
                    }
                    if normalize_agent_model_name(self.agent_model_override.trim()).is_none()
                        && self.agent_model_discovery_task.is_none()
                        && self.agent_discovered_models.is_empty()
                    {
                        self.start_agent_model_discovery_task(&system, false);
                    }
                }
                let catalog_base_url = system
                    .base_url
                    .as_deref()
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .unwrap_or("(transport default)");
                if self.agent_base_url_override.trim().is_empty() {
                    ui.small(format!("base URL: {catalog_base_url}"));
                } else {
                    ui.small(format!(
                        "base URL: {} (session override)",
                        self.agent_base_url_override.trim()
                    ));
                }
                let catalog_model = system
                    .model
                    .as_deref()
                    .map(str::trim)
                    .and_then(normalize_agent_model_name)
                    .unwrap_or_else(|| match system.transport {
                        AgentSystemTransport::NativeAnthropic => {
                            GUI_ANTHROPIC_DEFAULT_MODEL.to_string()
                        }
                        _ => OPENAI_COMPAT_UNSPECIFIED_MODEL.to_string(),
                    });
                let model_override = normalize_agent_model_name(self.agent_model_override.trim());
                if let Some(model_override) = model_override {
                    ui.small(format!("model: {model_override} (session override)"));
                } else if normalize_agent_model_name(&self.agent_discovered_model_pick).is_some() {
                    ui.small(format!(
                        "model: {} (selected discovered model)",
                        self.agent_discovered_model_pick.trim()
                    ));
                } else {
                    ui.small(format!("model: {catalog_model}"));
                }
            } else {
                self.clear_agent_model_discovery_snapshot();
            }
            if !available {
                ui.colored_label(
                    egui::Color32::from_rgb(190, 70, 70),
                    format!(
                        "Unavailable: {}",
                        reason.unwrap_or_else(|| "unknown reason".to_string())
                    ),
                );
            }
            if system.id == "builtin_echo" {
                ui.small(
                    "Built-in Echo is only a demo bridge. Use one of the quick-start buttons above for a real model-backed assistant.",
                );
            }
        } else if self.agent_systems.is_empty() {
            ui.small("No systems loaded from this catalog.");
        }
        let selected_transport = self
            .selected_agent_system()
            .map(|system| system.transport)
            .unwrap_or_default();
        let (key_label, key_hint) = match selected_transport {
            AgentSystemTransport::NativeAnthropic => ("Anthropic API key", "sk-ant-..."),
            _ => ("OpenAI API key", "sk-..."),
        };
        ui.horizontal(|ui| {
            ui.label(key_label);
            let response = ui.add(
                egui::TextEdit::singleline(&mut self.agent_openai_api_key)
                    .password(true)
                    .hint_text(key_hint),
            );
            preflight_inputs_changed |= response.changed();
            if ui
                .button("Clear Key")
                .on_hover_text("Clear session-only API key override")
                .clicked()
            {
                self.agent_openai_api_key.clear();
                preflight_inputs_changed = true;
            }
        });
        ui.horizontal(|ui| {
            ui.label("Base URL override");
            let response = ui.add(
                egui::TextEdit::singleline(&mut self.agent_base_url_override)
                    .hint_text("http://localhost:11964/v1"),
            );
            preflight_inputs_changed |= response.changed();
            if ui
                .button("Clear URL")
                .on_hover_text("Clear session-only base URL override")
                .clicked()
            {
                self.agent_base_url_override.clear();
                preflight_inputs_changed = true;
            }
        });
        ui.horizontal(|ui| {
            ui.label("Model override");
            let response = ui.add(
                egui::TextEdit::singleline(&mut self.agent_model_override).hint_text("unspecified"),
            );
            preflight_inputs_changed |= response.changed();
            if ui
                .button("Clear Model")
                .on_hover_text("Clear session-only model override")
                .clicked()
            {
                self.agent_model_override.clear();
                preflight_inputs_changed = true;
            }
        });
        ui.horizontal(|ui| {
            ui.label("timeout_sec");
            let response = ui.add(
                egui::TextEdit::singleline(&mut self.agent_timeout_secs)
                    .desired_width(100.0)
                    .hint_text("default"),
            );
            preflight_inputs_changed |= response.changed();
            if ui
                .button("Clear Timeout")
                .on_hover_text("Use default timeout")
                .clicked()
            {
                self.agent_timeout_secs.clear();
                preflight_inputs_changed = true;
            }
        });
        ui.horizontal(|ui| {
            ui.label("connect_timeout_sec");
            let connect_response = ui.add(
                egui::TextEdit::singleline(&mut self.agent_connect_timeout_secs)
                    .desired_width(90.0)
                    .hint_text("default"),
            );
            preflight_inputs_changed |= connect_response.changed();
            ui.label("read_timeout_sec");
            let read_response = ui.add(
                egui::TextEdit::singleline(&mut self.agent_read_timeout_secs)
                    .desired_width(90.0)
                    .hint_text("default"),
            );
            preflight_inputs_changed |= read_response.changed();
            if ui
                .button("Clear HTTP Timeouts")
                .on_hover_text("Use default connect/read timeouts")
                .clicked()
            {
                self.agent_connect_timeout_secs.clear();
                self.agent_read_timeout_secs.clear();
                preflight_inputs_changed = true;
            }
        });
        ui.horizontal(|ui| {
            ui.label("max_retries");
            let retries_response = ui.add(
                egui::TextEdit::singleline(&mut self.agent_max_retries)
                    .desired_width(90.0)
                    .hint_text("default"),
            );
            preflight_inputs_changed |= retries_response.changed();
            ui.label("max_response_bytes");
            let bytes_response = ui.add(
                egui::TextEdit::singleline(&mut self.agent_max_response_bytes)
                    .desired_width(120.0)
                    .hint_text("default"),
            );
            preflight_inputs_changed |= bytes_response.changed();
            if ui
                .button("Clear Limits")
                .on_hover_text("Use default retry/response-size limits")
                .clicked()
            {
                self.agent_max_retries.clear();
                self.agent_max_response_bytes.clear();
                preflight_inputs_changed = true;
            }
        });
        if let Some(system) = self.selected_agent_system() {
            ui.horizontal(|ui| {
                if ui
                    .button("Test Setup")
                    .on_hover_text(
                        "Validate the current system, key, endpoint, model, and runtime settings without sending a prompt",
                    )
                    .clicked()
                {
                    self.run_agent_preflight_probe();
                }
                if self.agent_preflight_output.is_some()
                    && ui
                        .button("Clear Test")
                        .on_hover_text("Clear the latest setup-preflight snapshot")
                        .clicked()
                {
                    self.clear_agent_preflight_output();
                }
                if matches!(
                    system.transport,
                    AgentSystemTransport::NativeOpenai
                        | AgentSystemTransport::NativeAnthropic
                        | AgentSystemTransport::NativeOpenaiCompat
                ) {
                    if ui
                        .button("Discover Models")
                        .on_hover_text("Query local/server model list from current base URL")
                        .clicked()
                    {
                        self.start_agent_model_discovery_task(&system, true);
                    }
                    if let Some(task) = &self.agent_model_discovery_task {
                        ui.add(egui::Spinner::new());
                        let status = if self.agent_model_discovery_status.trim().is_empty() {
                            "Discovering models".to_string()
                        } else {
                            self.agent_model_discovery_status.clone()
                        };
                        ui.small(format!("{status} ({:.1}s)", task.started.elapsed().as_secs_f32()));
                    }
                }
            });
            if !matches!(
                system.transport,
                AgentSystemTransport::NativeOpenai
                    | AgentSystemTransport::NativeAnthropic
                    | AgentSystemTransport::NativeOpenaiCompat
            ) {
                self.clear_agent_model_discovery_snapshot();
            } else {
                if !self.agent_discovered_models.is_empty() {
                    let previous_pick = self.agent_discovered_model_pick.clone();
                    ui.horizontal(|ui| {
                        ui.label("Discovered model");
                        egui::ComboBox::from_id_salt("agent_discovered_model_combo")
                            .selected_text(if self.agent_discovered_model_pick.trim().is_empty() {
                                "(choose model)"
                            } else {
                                self.agent_discovered_model_pick.as_str()
                            })
                            .show_ui(ui, |ui| {
                                for model in &self.agent_discovered_models {
                                    ui.selectable_value(
                                        &mut self.agent_discovered_model_pick,
                                        model.clone(),
                                        model,
                                    );
                                }
                            });
                    });
                    preflight_inputs_changed |= previous_pick != self.agent_discovered_model_pick;
                    ui.small(
                        "If Model override is unspecified, the selected discovered model is used.",
                    );
                }
                if !self.agent_model_discovery_status.trim().is_empty() {
                    ui.small(self.agent_model_discovery_status.clone());
                }
            }
        }
        if preflight_inputs_changed {
            self.invalidate_agent_preflight_after_setup_input_change();
        }
        if let Some(preflight) = &self.agent_preflight_output {
            ui.group(|ui| {
                ui.strong("Setup preflight");
                let color = if preflight.available {
                    egui::Color32::from_rgb(60, 140, 80)
                } else {
                    egui::Color32::from_rgb(190, 70, 70)
                };
                ui.colored_label(
                    color,
                    format!(
                        "{} ({})",
                        if preflight.available {
                            "Available"
                        } else {
                            "Unavailable"
                        },
                        preflight.transport
                    ),
                );
                if let Some(reason) = preflight.availability_reason.as_deref() {
                    if !reason.trim().is_empty() {
                        ui.small(format!("detail: {}", reason.trim()));
                    }
                }
                if let Some(base_url) = preflight.base_url.as_deref() {
                    ui.small(format!("base URL: {base_url}"));
                }
                if let Some(model) = preflight.model.as_deref() {
                    ui.small(format!("model: {model}"));
                }
                ui.small(format!(
                    "runtime: timeout={}s | connect={}s | read={}s | retries={} | max_response_bytes={}",
                    preflight.timeout_secs,
                    preflight.connect_timeout_secs,
                    preflight.read_timeout_secs,
                    preflight.max_retries,
                    preflight.max_response_bytes
                ));
                if !preflight.endpoint_candidates.is_empty() {
                    ui.small(format!(
                        "request endpoints: {}",
                        preflight.endpoint_candidates.join(" | ")
                    ));
                }
                if !preflight.model_endpoint_candidates.is_empty() {
                    ui.small(format!(
                        "model discovery endpoints: {}",
                        preflight.model_endpoint_candidates.join(" | ")
                    ));
                }
                if !preflight.warnings.is_empty() {
                    ui.colored_label(
                        egui::Color32::from_rgb(180, 120, 50),
                        format!("warnings: {}", preflight.warnings.join(" | ")),
                    );
                }
                if let Some(live) = &preflight.live_probe {
                    ui.separator();
                    ui.strong("Live probe");
                    let color = match live.status_class {
                        AgentLiveProbeStatusClass::Ok => egui::Color32::from_rgb(60, 140, 80),
                        AgentLiveProbeStatusClass::MissingKey
                        | AgentLiveProbeStatusClass::AuthFailed
                        | AgentLiveProbeStatusClass::QuotaOrBilling
                        | AgentLiveProbeStatusClass::ModelMissing
                        | AgentLiveProbeStatusClass::EndpointUnreachable => {
                            egui::Color32::from_rgb(190, 70, 70)
                        }
                        AgentLiveProbeStatusClass::UnsupportedTransport
                        | AgentLiveProbeStatusClass::ProviderError => {
                            egui::Color32::from_rgb(180, 120, 50)
                        }
                    };
                    ui.colored_label(color, live.status_class.as_str());
                    if !live.message.trim().is_empty() {
                        ui.small(format!("detail: {}", live.message.trim()));
                    }
                    ui.small(format!(
                        "reachable={} | auth_ok={} | model_list_ok={} | selected_model_seen={}",
                        live.reachable, live.auth_ok, live.model_list_ok, live.selected_model_seen
                    ));
                    if !live.attempted_endpoints.is_empty() {
                        ui.small(format!(
                            "attempted endpoints: {}",
                            live.attempted_endpoints.join(" | ")
                        ));
                    }
                    if let Some(endpoint) = live.selected_endpoint.as_deref() {
                        ui.small(format!("selected endpoint: {endpoint}"));
                    }
                    if let Some(code) = live.provider_error_code.as_deref() {
                        ui.small(format!("provider error code: {code}"));
                    }
                }
                let next_actions = Self::agent_preflight_next_actions(preflight);
                if !next_actions.is_empty() {
                    ui.separator();
                    ui.strong("Next action");
                    for action in next_actions {
                        ui.small(action);
                    }
                }
            });
        }
        ui.small(
            "Session only: if set, this key overrides OPENAI_API_KEY for agent requests started from this GUI window.",
        );
        ui.small(
            "Session only: Base URL override applies to native_openai/native_openai_compat. For local roots (e.g. http://localhost:11964), GENtle tries /chat/completions and /v1/chat/completions on that same base URL.",
        );
        ui.small(
            "Session only: Model override applies to native_openai/native_openai_compat and maps to GENTLE_AGENT_MODEL. Value 'unspecified' means no override.",
        );
        ui.small(
            "Session only: timeout_sec maps to GENTLE_AGENT_TIMEOUT_SECS and applies to agent requests (stdio and native transports).",
        );
        ui.small(
            "Session only: connect_timeout_sec/read_timeout_sec map to GENTLE_AGENT_CONNECT_TIMEOUT_SECS/GENTLE_AGENT_READ_TIMEOUT_SECS.",
        );
        ui.small(
            "Session only: max_retries/max_response_bytes map to GENTLE_AGENT_MAX_RETRIES/GENTLE_AGENT_MAX_RESPONSE_BYTES.",
        );
        ui.horizontal(|ui| {
            ui.checkbox(
                &mut self.agent_include_state_summary,
                "Include project state summary in request",
            );
            ui.checkbox(
                &mut self.agent_allow_auto_exec,
                "Auto-run suggestions marked as 'auto'",
            );
        });
        if !agent_prompt_template_options()
            .iter()
            .any(|(id, _)| *id == self.agent_prompt_template_id)
        {
            self.agent_prompt_template_id = AGENT_PROMPT_TEMPLATE_DEFAULT_ID.to_string();
        }
        ui.horizontal(|ui| {
            ui.label("Prompt template");
            egui::ComboBox::from_id_salt("agent_prompt_template_combo")
                .selected_text(agent_prompt_template_label(&self.agent_prompt_template_id))
                .show_ui(ui, |ui| {
                    for (id, label) in agent_prompt_template_options() {
                        ui.selectable_value(
                            &mut self.agent_prompt_template_id,
                            (*id).to_string(),
                            *label,
                        );
                    }
                });
            if ui
                .button("Insert")
                .on_hover_text("Replace current prompt with selected template")
                .clicked()
            {
                self.agent_prompt =
                    agent_prompt_template_text(&self.agent_prompt_template_id).to_string();
            }
            if ui
                .button("Append")
                .on_hover_text("Append selected template below current prompt")
                .clicked()
            {
                let template_text = agent_prompt_template_text(&self.agent_prompt_template_id);
                if self.agent_prompt.trim().is_empty() {
                    self.agent_prompt = template_text.to_string();
                } else {
                    if !self.agent_prompt.ends_with('\n') {
                        self.agent_prompt.push('\n');
                    }
                    self.agent_prompt.push('\n');
                    self.agent_prompt.push_str(template_text);
                }
            }
        });
        ui.label("Prompt");
        ui.add(
            egui::TextEdit::multiline(&mut self.agent_prompt)
                .desired_rows(6)
                .desired_width(f32::INFINITY),
        );
        let running = self.agent_task.is_some();
        ui.horizontal(|ui| {
            if ui
                .add_enabled(
                    !running && selected_available,
                    egui::Button::new("Ask Agent"),
                )
                .on_hover_text("Send prompt to selected agent system")
                .clicked()
            {
                self.start_agent_assistant_request();
            }
            if ui
                .button("Clear Response")
                .on_hover_text("Clear latest agent response and status")
                .clicked()
            {
                self.agent_last_invocation = None;
                self.agent_status.clear();
            }
            if ui
                .button("Clear Execution Log")
                .on_hover_text("Clear local execution history for agent suggestions")
                .clicked()
            {
                self.agent_execution_log.clear();
            }
        });
        if let Some(task) = &self.agent_task {
            ui.horizontal(|ui| {
                ui.add(egui::Spinner::new());
                ui.label(format!(
                    "Agent request running ({:.1}s)",
                    task.started.elapsed().as_secs_f32()
                ));
            });
        }
        if !self.agent_status.is_empty() {
            ui.separator();
            ui.monospace(self.agent_status.clone());
        }

        if let Some(invocation) = self.agent_last_invocation.clone() {
            ui.separator();
            ui.label(format!(
                "Latest response from {} ({})",
                invocation.system_label, invocation.system_id
            ));
            ui.small(format!(
                "elapsed={} ms | transport={} | exit_code={:?}",
                invocation.elapsed_ms, invocation.transport, invocation.exit_code
            ));
            ui.small(format!(
                "runtime: timeout={}s | connect={:?}s | read={:?}s | max_retries={} | max_response_bytes={}",
                invocation.runtime.timeout_secs,
                invocation.runtime.connect_timeout_secs,
                invocation.runtime.read_timeout_secs,
                invocation.runtime.max_retries,
                invocation.runtime.max_response_bytes
            ));
            if !invocation.runtime.endpoint_candidates.is_empty() {
                ui.small(format!(
                    "endpoint candidates: {}",
                    invocation.runtime.endpoint_candidates.join(" | ")
                ));
            }
            if !invocation.runtime.attempted_endpoints.is_empty() {
                ui.small(format!(
                    "attempted endpoints: {}",
                    invocation.runtime.attempted_endpoints.join(" | ")
                ));
            }
            if let Some(selected_endpoint) = invocation.runtime.selected_endpoint.as_deref() {
                ui.small(format!("selected endpoint: {}", selected_endpoint));
            }
            if !invocation.response.assistant_message.trim().is_empty() {
                ui.group(|ui| {
                    ui.strong("Agent message");
                    ui.label(invocation.response.assistant_message.trim());
                });
            }
            if !invocation.response.questions.is_empty() {
                ui.group(|ui| {
                    ui.strong("Agent questions");
                    for question in &invocation.response.questions {
                        ui.label(format!("- {}", question));
                    }
                });
            }
            if invocation.response.suggested_commands.is_empty() {
                ui.small("No executable suggestions in this reply.");
            } else {
                ui.separator();
                ui.strong("Suggested commands");
                let mut run_request: Option<(usize, String)> = None;
                egui::Grid::new("agent_suggested_commands_grid")
                    .striped(true)
                    .show(ui, |ui| {
                        ui.strong("#");
                        ui.strong("Intent");
                        ui.strong("Command");
                        ui.strong("Rationale");
                        ui.strong("Action");
                        ui.end_row();
                        for (idx, suggestion) in invocation
                            .response
                            .suggested_commands
                            .iter()
                            .enumerate()
                        {
                            let index_1based = idx + 1;
                            ui.label(index_1based.to_string());
                            ui.label(suggestion.execution.as_str());
                            ui.monospace(suggestion.command.as_str());
                            ui.label(suggestion.rationale.clone().unwrap_or_default());
                            let can_run = !suggestion.command.trim().is_empty()
                                && suggestion.execution != AgentExecutionIntent::Chat;
                            let run_resp = ui
                                .add_enabled(can_run, egui::Button::new("Run"))
                                .on_hover_text(
                                    "Execute this suggested command using GENtle shared shell parser/executor",
                                );
                            if run_resp.clicked() {
                                run_request = Some((index_1based, suggestion.command.clone()));
                            }
                            ui.end_row();
                        }
                    });
                if let Some((index_1based, command)) = run_request {
                    self.execute_agent_suggested_command(index_1based, &command, "manual");
                }
            }
            if !invocation.raw_stderr.trim().is_empty() {
                ui.separator();
                ui.strong("Agent stderr");
                let mut stderr = invocation.raw_stderr.clone();
                ui.add(
                    egui::TextEdit::multiline(&mut stderr)
                        .desired_rows(4)
                        .desired_width(f32::INFINITY),
                );
            }
        }

        if !self.agent_execution_log.is_empty() {
            ui.separator();
            ui.strong("Execution log");
            egui::ScrollArea::vertical()
                .max_height(180.0)
                .show(ui, |ui| {
                    scroll_input_policy::apply_scrollarea_keyboard_navigation(
                        ui,
                        scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                    );
                    for entry in self.agent_execution_log.iter().rev() {
                        ui.label(format!(
                            "#{} [{}] {} | {} | changed={} | t={}",
                            entry.index_1based,
                            entry.trigger,
                            if entry.ok { "ok" } else { "error" },
                            entry.command,
                            entry.state_changed,
                            entry.executed_at_unix_ms
                        ));
                        ui.small(entry.summary.clone());
                    }
                });
        }
        close_requested
    }

    pub(super) fn render_agent_assistant_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_agent_assistant_dialog {
            return;
        }
        let mut open = self.show_agent_assistant_dialog;
        let viewport_id = Self::agent_assistant_viewport_id();
        let spec = self
            .hosted_window_spec_for_viewport(
                "Agent Assistant",
                Self::hosted_agent_assistant_window_id(),
                viewport_id,
                Vec2::new(980.0, 720.0),
                Vec2::new(640.0, 420.0),
            )
            .legacy_layer_id(crate::egui_compat::hosted_window_title_layer_id(
                "Agent Assistant",
            ))
            .legacy_layer_id(egui::LayerId::new(
                egui::Order::Middle,
                egui::Id::new(viewport_id),
            ));
        if ctx.embed_viewports() {
            let mut close_requested = false;
            crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
                close_requested = self.render_agent_assistant_contents(ui);
            });
            if close_requested {
                open = false;
            }
            if ctx.input(|i| i.key_pressed(Key::Escape)) {
                open = false;
            }
            self.show_agent_assistant_dialog = open;
            self.finalize_viewport_open_probe(viewport_id, "Agent Assistant");
            return;
        }
        let viewport_spec = self.hosted_window_spec_for_viewport(
            "Agent Assistant",
            Self::hosted_agent_assistant_window_id(),
            viewport_id,
            Vec2::new(980.0, 720.0),
            Vec2::new(640.0, 420.0),
        );
        let builder = crate::egui_compat::viewport_builder_for_hosted_window(&viewport_spec);
        ctx.show_viewport_immediate(viewport_id, builder, |ctx, class| {
            self.note_viewport_focus_if_active(ctx, viewport_id);
            if class == egui::ViewportClass::EmbeddedWindow {
                let mut close_requested = false;
                crate::egui_compat::show_hosted_window(ctx, &viewport_spec, &mut open, |ui| {
                    close_requested = self.render_agent_assistant_contents(ui);
                });
                if close_requested {
                    open = false;
                }
            } else {
                let mut close_requested = false;
                crate::egui_compat::show_central_panel(ctx, egui::CentralPanel::default(), |ui| {
                    close_requested = self.render_agent_assistant_contents(ui);
                });

                if close_requested || Self::viewport_close_requested_or_shortcut(ctx) {
                    open = false;
                }
            }
        });
        if ctx.input(|i| i.key_pressed(Key::Escape)) {
            open = false;
        }
        self.show_agent_assistant_dialog = open;
    }
    pub(super) fn push_unique_trace_token(values: &mut Vec<String>, token: &str) {
        let compact = token.trim();
        if compact.is_empty() || values.iter().any(|existing| existing == compact) {
            return;
        }
        values.push(compact.to_string());
    }

    pub(super) fn normalize_routine_assistant_preflight_snapshot(
        snapshot: &mut RoutineDecisionTracePreflightSnapshot,
    ) {
        let mut warnings = vec![];
        for warning in std::mem::take(&mut snapshot.warnings) {
            Self::push_unique_trace_token(&mut warnings, &warning);
        }
        snapshot.warnings = warnings;

        let mut errors = vec![];
        for error in std::mem::take(&mut snapshot.errors) {
            Self::push_unique_trace_token(&mut errors, &error);
        }
        snapshot.errors = errors;

        snapshot.contract_source = snapshot
            .contract_source
            .take()
            .map(|v| v.trim().to_string())
            .filter(|v| !v.is_empty());
    }

    pub(super) fn routine_assistant_disambiguation_question_id_from_text(text: &str) -> String {
        let mut out = String::new();
        let mut last_was_sep = false;
        for ch in text.trim().chars().flat_map(|c| c.to_lowercase()) {
            if ch.is_ascii_alphanumeric() {
                out.push(ch);
                last_was_sep = false;
            } else if !last_was_sep {
                out.push('_');
                last_was_sep = true;
            }
            if out.len() >= 48 {
                break;
            }
        }
        let compact = out.trim_matches('_').to_string();
        if compact.is_empty() {
            "question".to_string()
        } else {
            compact
        }
    }

    pub(super) fn routine_assistant_disambiguation_questions_from_output(
        output: &serde_json::Value,
    ) -> Vec<RoutineDecisionTraceDisambiguationQuestion> {
        let mut question_texts: Vec<String> = vec![];
        if let Some(rows) = output
            .get("explanation")
            .and_then(|value| value.get("disambiguation_questions"))
            .and_then(|value| value.as_array())
        {
            for row in rows {
                if let Some(text) = row.as_str() {
                    Self::push_unique_trace_token(&mut question_texts, text);
                }
            }
        }
        if let Some(rows) = output
            .get("comparison")
            .and_then(|value| value.get("disambiguation_questions"))
            .and_then(|value| value.as_array())
        {
            for row in rows {
                if let Some(text) = row.as_str() {
                    Self::push_unique_trace_token(&mut question_texts, text);
                }
            }
        }

        let mut out: Vec<RoutineDecisionTraceDisambiguationQuestion> = vec![];
        let mut used_ids: HashMap<String, usize> = HashMap::new();
        for text in question_texts {
            let question_text = text.trim().to_string();
            if question_text.is_empty() {
                continue;
            }
            let base_id =
                Self::routine_assistant_disambiguation_question_id_from_text(&question_text);
            let count = used_ids.entry(base_id.clone()).or_insert(0);
            *count += 1;
            let question_id = if *count == 1 {
                base_id
            } else {
                format!("{}_{}", base_id, *count)
            };
            out.push(RoutineDecisionTraceDisambiguationQuestion {
                question_id,
                question_text,
            });
        }
        out
    }

    pub(super) fn merge_routine_assistant_disambiguation_questions(
        existing: &mut Vec<RoutineDecisionTraceDisambiguationQuestion>,
        incoming: Vec<RoutineDecisionTraceDisambiguationQuestion>,
    ) {
        for mut row in incoming {
            row.question_id = row.question_id.trim().to_string();
            row.question_text = row.question_text.trim().to_string();
            if row.question_text.is_empty() {
                continue;
            }
            if row.question_id.is_empty() {
                row.question_id = Self::routine_assistant_disambiguation_question_id_from_text(
                    &row.question_text,
                );
            }
            if row.question_id.is_empty() {
                continue;
            }
            if existing.iter().any(|present| {
                present.question_id.eq_ignore_ascii_case(&row.question_id)
                    || present
                        .question_text
                        .eq_ignore_ascii_case(&row.question_text)
            }) {
                continue;
            }
            existing.push(row);
        }
    }

    pub(super) fn routine_assistant_commit_preflight_snapshot(
        trace: &mut RoutineDecisionTrace,
        mut snapshot: Option<RoutineDecisionTracePreflightSnapshot>,
    ) {
        if let Some(snapshot) = snapshot.as_mut() {
            Self::normalize_routine_assistant_preflight_snapshot(snapshot);
            trace.preflight_history.push(snapshot.clone());
            trace.preflight_snapshot = Some(snapshot.clone());
        } else {
            trace.preflight_snapshot = None;
        }
    }

    pub(super) fn next_routine_assistant_trace_id(&mut self) -> String {
        let counter = self.routine_assistant_trace_counter.max(1);
        self.routine_assistant_trace_counter = counter.saturating_add(1);
        format!("routine_assistant_{}_{}", Self::now_unix_ms(), counter)
    }

    pub(super) fn routine_assistant_candidate_ids_snapshot(&self) -> Vec<String> {
        let mut out: Vec<String> = vec![];
        for row in &self.routine_assistant_candidates {
            Self::push_unique_trace_token(&mut out, &row.routine_id);
        }
        out
    }

    pub(super) fn routine_assistant_construct_reasoning_seq_id(&self) -> Option<String> {
        self.active_dna_window_context()
            .map(|(seq_id, _)| seq_id)
            .or_else(|| {
                self.routine_assistant_preference_context
                    .as_ref()
                    .and_then(|context| context.construct_reasoning_seq_id.clone())
            })
            .or_else(|| {
                self.routine_assistant_decision_trace
                    .as_ref()
                    .and_then(|trace| trace.routine_preference_context.as_ref())
                    .and_then(|context| context.construct_reasoning_seq_id.clone())
            })
    }

    pub(super) fn routine_assistant_candidate_planning_scores_snapshot(
        &self,
    ) -> Vec<RoutineDecisionTraceCandidateScore> {
        let mut out = self
            .routine_assistant_candidates
            .iter()
            .map(|row| {
                let (routine_family_alignment_bonus, routine_family_alignment_sources) = row
                    .planning_estimate
                    .as_ref()
                    .map(|estimate| {
                        let bonus = estimate
                            .explanation
                            .get("routine_family_alignment_bonus")
                            .and_then(|value| value.as_f64());
                        let sources = estimate
                            .explanation
                            .get("routine_family_alignment_sources")
                            .and_then(|value| value.as_array())
                            .map(|rows| {
                                rows.iter()
                                    .filter_map(|row| row.as_str())
                                    .map(str::trim)
                                    .filter(|value| !value.is_empty())
                                    .map(|value| value.to_string())
                                    .collect::<Vec<_>>()
                            })
                            .unwrap_or_default();
                        (bonus, sources)
                    })
                    .unwrap_or((None, vec![]));
                RoutineDecisionTraceCandidateScore {
                    routine_id: row.routine_id.clone(),
                    routine_title: Some(row.title.clone()).filter(|value| !value.is_empty()),
                    routine_family: row.family.clone(),
                    passes_guardrails: row
                        .planning_estimate
                        .as_ref()
                        .map(|estimate| estimate.passes_guardrails)
                        .unwrap_or(false),
                    estimated_time_hours: row.estimated_time_hours,
                    estimated_cost: row.estimated_cost,
                    local_fit_score: row.local_fit_score,
                    composite_meta_score: row.composite_meta_score,
                    routine_family_alignment_bonus,
                    routine_family_alignment_sources,
                }
            })
            .collect::<Vec<_>>();
        out.sort_by(|left, right| {
            right
                .passes_guardrails
                .cmp(&left.passes_guardrails)
                .then_with(|| {
                    right
                        .composite_meta_score
                        .unwrap_or(f64::NEG_INFINITY)
                        .total_cmp(&left.composite_meta_score.unwrap_or(f64::NEG_INFINITY))
                })
                .then_with(|| left.routine_family.cmp(&right.routine_family))
                .then_with(|| left.routine_id.cmp(&right.routine_id))
        });
        out
    }

    pub(super) fn routine_assistant_planning_trace_artifacts(
        &self,
        selected_routine: Option<&CloningRoutineCatalogRow>,
    ) -> (
        Option<RoutinePreferenceContextRecord>,
        Vec<RoutineDecisionTraceCandidateScore>,
        Vec<MacroTemplateSuggestion>,
    ) {
        let selected_routine_id = selected_routine
            .map(|row| row.routine_id.trim().to_string())
            .or_else(|| {
                let value = self.routine_assistant_selected_routine_id.trim();
                (!value.is_empty()).then_some(value.to_string())
            });
        let selected_routine_family = selected_routine
            .map(|row| row.family.trim().to_string())
            .or_else(|| {
                self.routine_assistant_decision_trace
                    .as_ref()
                    .and_then(|trace| trace.selected_routine_family.as_deref())
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .map(|value| value.to_string())
            });
        let construct_reasoning_seq_id = self.routine_assistant_construct_reasoning_seq_id();
        let (preference_context, macro_suggestions) = self
            .engine
            .write()
            .ok()
            .map(|mut engine| {
                let context = engine.planning_routine_preference_context_record_for_sequence(
                    construct_reasoning_seq_id.as_deref(),
                );
                let suggestions = engine.suggest_macro_templates_for_routine_for_sequence(
                    selected_routine_id.as_deref(),
                    selected_routine_family.as_deref(),
                    construct_reasoning_seq_id.as_deref(),
                    6,
                );
                (Some(context), suggestions)
            })
            .unwrap_or_else(|| (None, vec![]));
        let candidate_planning_scores = self.routine_assistant_candidate_planning_scores_snapshot();
        (
            preference_context,
            candidate_planning_scores,
            macro_suggestions,
        )
    }

    pub(super) fn routine_assistant_bindings_snapshot(&self) -> BTreeMap<String, String> {
        let mut out: BTreeMap<String, String> = BTreeMap::new();
        for (key, value) in &self.routine_assistant_bindings {
            let key = key.trim();
            let value = value.trim();
            if key.is_empty() || value.is_empty() {
                continue;
            }
            out.insert(key.to_string(), value.to_string());
        }
        out
    }

    pub(super) fn routine_assistant_effective_disambiguation_questions(
        &self,
    ) -> Vec<RoutineDecisionTraceDisambiguationQuestion> {
        let mut questions: Vec<RoutineDecisionTraceDisambiguationQuestion> = vec![];
        if let Some(output) = self.routine_assistant_explain_output.as_ref() {
            Self::merge_routine_assistant_disambiguation_questions(
                &mut questions,
                Self::routine_assistant_disambiguation_questions_from_output(output),
            );
        }
        if let Some(output) = self.routine_assistant_compare_output.as_ref() {
            Self::merge_routine_assistant_disambiguation_questions(
                &mut questions,
                Self::routine_assistant_disambiguation_questions_from_output(output),
            );
        }
        if let Some(trace) = self.routine_assistant_decision_trace.as_ref() {
            Self::merge_routine_assistant_disambiguation_questions(
                &mut questions,
                trace.disambiguation_questions_presented.clone(),
            );
        }
        questions
    }

    pub(super) fn sync_routine_assistant_disambiguation_answers_for_questions(
        &mut self,
        questions: &[RoutineDecisionTraceDisambiguationQuestion],
        fallback_answers: &[RoutineDecisionTraceDisambiguationAnswer],
    ) {
        let mut fallback_by_question_id: HashMap<String, String> = HashMap::new();
        for row in fallback_answers {
            let question_id = row.question_id.trim();
            if question_id.is_empty() {
                continue;
            }
            let answer_text = row.answer_text.trim();
            if answer_text.is_empty() {
                continue;
            }
            fallback_by_question_id
                .entry(question_id.to_ascii_lowercase())
                .or_insert_with(|| answer_text.to_string());
        }

        let mut next: BTreeMap<String, String> = BTreeMap::new();
        for row in questions {
            let question_id = row.question_id.trim();
            if question_id.is_empty() {
                continue;
            }
            let answer_text = self
                .routine_assistant_disambiguation_answers
                .get(question_id)
                .map(String::as_str)
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(|value| value.to_string())
                .or_else(|| {
                    fallback_by_question_id
                        .get(&question_id.to_ascii_lowercase())
                        .cloned()
                })
                .unwrap_or_default();
            next.insert(question_id.to_string(), answer_text);
        }
        self.routine_assistant_disambiguation_answers = next;
    }

    pub(super) fn routine_assistant_disambiguation_answers_snapshot(
        &self,
        questions: &[RoutineDecisionTraceDisambiguationQuestion],
    ) -> Vec<RoutineDecisionTraceDisambiguationAnswer> {
        let mut answers_by_question_id: BTreeMap<String, String> = BTreeMap::new();
        let mut seen_question_ids: HashSet<String> = HashSet::new();
        for row in questions {
            let question_id = row.question_id.trim();
            if question_id.is_empty() {
                continue;
            }
            if !seen_question_ids.insert(question_id.to_ascii_lowercase()) {
                continue;
            }
            let answer_text = self
                .routine_assistant_disambiguation_answers
                .get(question_id)
                .map(String::as_str)
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(|value| value.to_string());
            if let Some(answer_text) = answer_text {
                answers_by_question_id.insert(question_id.to_string(), answer_text);
            }
        }
        answers_by_question_id
            .into_iter()
            .map(
                |(question_id, answer_text)| RoutineDecisionTraceDisambiguationAnswer {
                    question_id,
                    answer_text,
                },
            )
            .collect::<Vec<_>>()
    }

    pub(super) fn normalize_routine_preference_context_for_gui(
        mut context: RoutinePreferenceContextRecord,
    ) -> RoutinePreferenceContextRecord {
        context.helper_profile_id = context
            .helper_profile_id
            .take()
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty());
        context.construct_reasoning_seq_id = context
            .construct_reasoning_seq_id
            .take()
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty());
        context.helper_resolution_status = context.helper_resolution_status.trim().to_string();
        if context.helper_resolution_status.is_empty() {
            context.helper_resolution_status = "not_requested".to_string();
        }
        let normalize_vec = |values: &mut Vec<String>| {
            let mut normalized = vec![];
            for value in std::mem::take(values) {
                Self::push_unique_trace_token(&mut normalized, &value);
            }
            *values = normalized;
        };
        normalize_vec(&mut context.explicit_preferred_routine_families);
        normalize_vec(&mut context.helper_derived_preferred_routine_families);
        normalize_vec(&mut context.variant_derived_preferred_routine_families);
        normalize_vec(&mut context.effective_preferred_routine_families);
        normalize_vec(&mut context.helper_offered_functions);
        normalize_vec(&mut context.helper_component_labels);
        normalize_vec(&mut context.variant_effect_tags);
        normalize_vec(&mut context.variant_suggested_assay_ids);
        normalize_vec(&mut context.rationale);
        context
    }

    pub(super) fn normalize_routine_decision_trace_candidate_score_for_gui(
        mut score: RoutineDecisionTraceCandidateScore,
    ) -> Option<RoutineDecisionTraceCandidateScore> {
        score.routine_id = score.routine_id.trim().to_string();
        if score.routine_id.is_empty() {
            return None;
        }
        score.routine_title = score
            .routine_title
            .take()
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty());
        score.routine_family = score.routine_family.trim().to_string();
        score.estimated_time_hours = score
            .estimated_time_hours
            .filter(|value| value.is_finite() && *value >= 0.0);
        score.estimated_cost = score
            .estimated_cost
            .filter(|value| value.is_finite() && *value >= 0.0);
        score.local_fit_score = score
            .local_fit_score
            .filter(|value| value.is_finite() && *value >= 0.0 && *value <= 1.0);
        score.composite_meta_score = score.composite_meta_score.filter(|value| value.is_finite());
        score.routine_family_alignment_bonus = score
            .routine_family_alignment_bonus
            .filter(|value| value.is_finite());
        let mut sources = vec![];
        for source in std::mem::take(&mut score.routine_family_alignment_sources) {
            Self::push_unique_trace_token(&mut sources, &source);
        }
        score.routine_family_alignment_sources = sources;
        Some(score)
    }

    pub(super) fn normalize_macro_template_suggestion_for_gui(
        mut suggestion: MacroTemplateSuggestion,
    ) -> Option<MacroTemplateSuggestion> {
        suggestion.macro_kind = suggestion.macro_kind.trim().to_string();
        if suggestion.macro_kind.is_empty() {
            return None;
        }
        suggestion.template_name = suggestion.template_name.trim().to_string();
        if suggestion.template_name.is_empty() {
            return None;
        }
        suggestion.description = suggestion
            .description
            .take()
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty());
        suggestion.details_url = suggestion
            .details_url
            .take()
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty());
        if !suggestion.score.is_finite() || suggestion.score < 0.0 {
            suggestion.score = 0.0;
        }
        let normalize_vec = |values: &mut Vec<String>| {
            let mut normalized = vec![];
            for value in std::mem::take(values) {
                Self::push_unique_trace_token(&mut normalized, &value);
            }
            *values = normalized;
        };
        normalize_vec(&mut suggestion.matched_routine_families);
        normalize_vec(&mut suggestion.matched_terms);
        normalize_vec(&mut suggestion.rationale);
        Some(suggestion)
    }

    pub(super) fn normalize_routine_decision_trace_for_gui(
        mut trace: RoutineDecisionTrace,
    ) -> Option<RoutineDecisionTrace> {
        let schema = trace.schema.trim();
        if schema.is_empty() {
            trace.schema = ROUTINE_DECISION_TRACE_SCHEMA.to_string();
        } else if !schema.eq_ignore_ascii_case(ROUTINE_DECISION_TRACE_SCHEMA) {
            return None;
        } else {
            trace.schema = ROUTINE_DECISION_TRACE_SCHEMA.to_string();
        }
        trace.trace_id = trace.trace_id.trim().to_string();
        if trace.trace_id.is_empty() {
            return None;
        }
        trace.source = trace.source.trim().to_string();
        if trace.source.is_empty() {
            trace.source = "gui_routine_assistant".to_string();
        }
        trace.status = trace.status.trim().to_string();
        if trace.status.is_empty() {
            trace.status = "draft".to_string();
        }
        trace.goal_text = trace.goal_text.trim().to_string();
        trace.query_text = trace.query_text.trim().to_string();
        if trace.created_at_unix_ms == 0 {
            trace.created_at_unix_ms = trace.updated_at_unix_ms;
        }
        if trace.updated_at_unix_ms == 0 {
            trace.updated_at_unix_ms = trace.created_at_unix_ms;
        }

        let normalize_opt = |value: &mut Option<String>| {
            *value = value
                .take()
                .map(|v| v.trim().to_string())
                .filter(|v| !v.is_empty());
        };
        normalize_opt(&mut trace.selected_routine_id);
        normalize_opt(&mut trace.selected_routine_title);
        normalize_opt(&mut trace.selected_routine_family);
        normalize_opt(&mut trace.macro_instance_id);
        normalize_opt(&mut trace.execution_error);
        trace.routine_preference_context = trace
            .routine_preference_context
            .take()
            .map(Self::normalize_routine_preference_context_for_gui);

        let mut normalized_candidates = vec![];
        for token in std::mem::take(&mut trace.candidate_routine_ids) {
            Self::push_unique_trace_token(&mut normalized_candidates, &token);
        }
        trace.candidate_routine_ids = normalized_candidates;

        let mut candidate_scores: Vec<RoutineDecisionTraceCandidateScore> = vec![];
        let mut seen_candidate_ids: HashSet<String> = HashSet::new();
        for row in std::mem::take(&mut trace.candidate_planning_scores) {
            let Some(row) = Self::normalize_routine_decision_trace_candidate_score_for_gui(row)
            else {
                continue;
            };
            if !seen_candidate_ids.insert(row.routine_id.to_ascii_lowercase()) {
                continue;
            }
            candidate_scores.push(row);
        }
        candidate_scores.sort_by(|left, right| {
            right
                .passes_guardrails
                .cmp(&left.passes_guardrails)
                .then_with(|| {
                    right
                        .composite_meta_score
                        .unwrap_or(f64::NEG_INFINITY)
                        .total_cmp(&left.composite_meta_score.unwrap_or(f64::NEG_INFINITY))
                })
                .then_with(|| left.routine_family.cmp(&right.routine_family))
                .then_with(|| left.routine_id.cmp(&right.routine_id))
        });
        trace.candidate_planning_scores = candidate_scores;

        let mut normalized_alternatives = vec![];
        for token in std::mem::take(&mut trace.alternatives_presented) {
            Self::push_unique_trace_token(&mut normalized_alternatives, &token);
        }
        trace.alternatives_presented = normalized_alternatives;

        let mut macro_suggestions: Vec<MacroTemplateSuggestion> = vec![];
        let mut seen_macro_keys: HashSet<String> = HashSet::new();
        for row in std::mem::take(&mut trace.macro_suggestions) {
            let Some(row) = Self::normalize_macro_template_suggestion_for_gui(row) else {
                continue;
            };
            let key = format!(
                "{}\u{1f}{}",
                row.macro_kind.to_ascii_lowercase(),
                row.template_name.to_ascii_lowercase()
            );
            if !seen_macro_keys.insert(key) {
                continue;
            }
            macro_suggestions.push(row);
        }
        macro_suggestions.sort_by(|left, right| {
            right
                .score
                .total_cmp(&left.score)
                .then_with(|| left.macro_kind.cmp(&right.macro_kind))
                .then_with(|| left.template_name.cmp(&right.template_name))
        });
        trace.macro_suggestions = macro_suggestions;

        let mut normalized_questions: Vec<RoutineDecisionTraceDisambiguationQuestion> = vec![];
        let mut used_question_ids: HashMap<String, usize> = HashMap::new();
        for mut row in std::mem::take(&mut trace.disambiguation_questions_presented) {
            row.question_id = row.question_id.trim().to_string();
            row.question_text = row.question_text.trim().to_string();
            if row.question_text.is_empty() {
                continue;
            }
            let base_id = if row.question_id.is_empty() {
                Self::routine_assistant_disambiguation_question_id_from_text(&row.question_text)
            } else {
                row.question_id.clone()
            };
            let count = used_question_ids.entry(base_id.clone()).or_insert(0);
            *count += 1;
            row.question_id = if *count == 1 {
                base_id
            } else {
                format!("{}_{}", base_id, *count)
            };
            if normalized_questions.iter().any(|existing| {
                existing.question_id.eq_ignore_ascii_case(&row.question_id)
                    || existing
                        .question_text
                        .eq_ignore_ascii_case(&row.question_text)
            }) {
                continue;
            }
            normalized_questions.push(row);
        }
        trace.disambiguation_questions_presented = normalized_questions;

        let mut normalized_answers_by_question: BTreeMap<String, String> = BTreeMap::new();
        for mut row in std::mem::take(&mut trace.disambiguation_answers) {
            row.question_id = row.question_id.trim().to_string();
            row.answer_text = row.answer_text.trim().to_string();
            if row.question_id.is_empty() || row.answer_text.is_empty() {
                continue;
            }
            normalized_answers_by_question.insert(row.question_id, row.answer_text);
        }
        trace.disambiguation_answers = normalized_answers_by_question
            .into_iter()
            .map(
                |(question_id, answer_text)| RoutineDecisionTraceDisambiguationAnswer {
                    question_id,
                    answer_text,
                },
            )
            .collect();

        let mut normalized_op_ids = vec![];
        for token in std::mem::take(&mut trace.emitted_operation_ids) {
            Self::push_unique_trace_token(&mut normalized_op_ids, &token);
        }
        trace.emitted_operation_ids = normalized_op_ids;

        let mut normalized_bindings: BTreeMap<String, String> = BTreeMap::new();
        for (key, value) in std::mem::take(&mut trace.bindings_snapshot) {
            let key = key.trim().to_string();
            let value = value.trim().to_string();
            if key.is_empty() || value.is_empty() {
                continue;
            }
            normalized_bindings.insert(key, value);
        }
        trace.bindings_snapshot = normalized_bindings;

        let mut preflight_history: Vec<RoutineDecisionTracePreflightSnapshot> = vec![];
        for mut snapshot in std::mem::take(&mut trace.preflight_history) {
            Self::normalize_routine_assistant_preflight_snapshot(&mut snapshot);
            preflight_history.push(snapshot);
        }
        trace.preflight_history = preflight_history;

        if let Some(snapshot) = trace.preflight_snapshot.as_mut() {
            Self::normalize_routine_assistant_preflight_snapshot(snapshot);
        }
        if trace.preflight_history.is_empty() {
            if let Some(snapshot) = trace.preflight_snapshot.clone() {
                trace.preflight_history.push(snapshot);
            }
        }
        trace.preflight_snapshot = trace.preflight_history.last().cloned();

        let mut comparisons: Vec<RoutineDecisionTraceComparison> = vec![];
        for mut row in std::mem::take(&mut trace.comparisons) {
            row.left_routine_id = row.left_routine_id.trim().to_string();
            row.right_routine_id = row.right_routine_id.trim().to_string();
            if row.left_routine_id.is_empty() || row.right_routine_id.is_empty() {
                continue;
            }
            if comparisons.iter().any(|existing| {
                existing.left_routine_id == row.left_routine_id
                    && existing.right_routine_id == row.right_routine_id
            }) {
                continue;
            }
            comparisons.push(row);
        }
        trace.comparisons = comparisons;

        let mut export_events: Vec<RoutineDecisionTraceExportEvent> = vec![];
        for mut event in std::mem::take(&mut trace.export_events) {
            event.run_bundle_path = event.run_bundle_path.trim().to_string();
            if event.run_bundle_path.is_empty() {
                continue;
            }
            if export_events.iter().any(|existing| {
                existing.run_bundle_path == event.run_bundle_path
                    && existing.exported_at_unix_ms == event.exported_at_unix_ms
            }) {
                continue;
            }
            export_events.push(event);
        }
        export_events.sort_by(|left, right| {
            left.exported_at_unix_ms
                .cmp(&right.exported_at_unix_ms)
                .then_with(|| left.run_bundle_path.cmp(&right.run_bundle_path))
        });
        trace.export_events = export_events;
        Some(trace)
    }

    pub(super) fn normalize_routine_decision_trace_store_for_gui(
        store: RoutineDecisionTraceStore,
    ) -> RoutineDecisionTraceStore {
        let mut by_trace_id: HashMap<String, RoutineDecisionTrace> = HashMap::new();
        for trace in store.traces {
            let Some(normalized) = Self::normalize_routine_decision_trace_for_gui(trace) else {
                continue;
            };
            let should_replace = by_trace_id
                .get(&normalized.trace_id)
                .map(|existing| {
                    (
                        normalized.updated_at_unix_ms,
                        normalized.created_at_unix_ms,
                        normalized.trace_id.as_str(),
                    ) > (
                        existing.updated_at_unix_ms,
                        existing.created_at_unix_ms,
                        existing.trace_id.as_str(),
                    )
                })
                .unwrap_or(true);
            if should_replace {
                by_trace_id.insert(normalized.trace_id.clone(), normalized);
            }
        }
        let mut traces = by_trace_id.into_values().collect::<Vec<_>>();
        traces.sort_by(|left, right| {
            left.created_at_unix_ms
                .cmp(&right.created_at_unix_ms)
                .then_with(|| left.trace_id.cmp(&right.trace_id))
        });
        RoutineDecisionTraceStore {
            schema: ROUTINE_DECISION_TRACE_STORE_SCHEMA.to_string(),
            traces,
        }
    }

    pub(super) fn load_routine_decision_trace_store_from_state(&self) -> RoutineDecisionTraceStore {
        let raw = self
            .engine
            .read()
            .unwrap()
            .state()
            .metadata
            .get(ROUTINE_DECISION_TRACES_METADATA_KEY)
            .cloned();
        let Some(raw) = raw else {
            return RoutineDecisionTraceStore {
                schema: ROUTINE_DECISION_TRACE_STORE_SCHEMA.to_string(),
                traces: vec![],
            };
        };
        if let Ok(mut store) = serde_json::from_value::<RoutineDecisionTraceStore>(raw.clone()) {
            if store.schema.trim().is_empty() {
                store.schema = ROUTINE_DECISION_TRACE_STORE_SCHEMA.to_string();
            }
            if !store
                .schema
                .trim()
                .eq_ignore_ascii_case(ROUTINE_DECISION_TRACE_STORE_SCHEMA)
            {
                return RoutineDecisionTraceStore {
                    schema: ROUTINE_DECISION_TRACE_STORE_SCHEMA.to_string(),
                    traces: vec![],
                };
            }
            return Self::normalize_routine_decision_trace_store_for_gui(store);
        }
        let traces = serde_json::from_value::<Vec<RoutineDecisionTrace>>(raw).unwrap_or_default();
        Self::normalize_routine_decision_trace_store_for_gui(RoutineDecisionTraceStore {
            schema: ROUTINE_DECISION_TRACE_STORE_SCHEMA.to_string(),
            traces,
        })
    }

    pub(super) fn persist_routine_decision_trace_store_to_state(
        &mut self,
        store: RoutineDecisionTraceStore,
    ) {
        let normalized = Self::normalize_routine_decision_trace_store_for_gui(store);
        let Ok(value) = serde_json::to_value(&normalized) else {
            return;
        };
        let mut engine = self.engine.write().unwrap();
        let state = engine.state_mut();
        if state.metadata.get(ROUTINE_DECISION_TRACES_METADATA_KEY) == Some(&value) {
            return;
        }
        state
            .metadata
            .insert(ROUTINE_DECISION_TRACES_METADATA_KEY.to_string(), value);
    }

    pub(super) fn persist_routine_assistant_decision_trace(&mut self) {
        let Some(active_trace) = self.routine_assistant_decision_trace.clone() else {
            return;
        };
        let Some(active_trace) = Self::normalize_routine_decision_trace_for_gui(active_trace)
        else {
            return;
        };
        self.routine_assistant_decision_trace = Some(active_trace.clone());
        let mut store = self.load_routine_decision_trace_store_from_state();
        let mut replaced = false;
        for trace in &mut store.traces {
            if trace.trace_id == active_trace.trace_id {
                *trace = active_trace.clone();
                replaced = true;
                break;
            }
        }
        if !replaced {
            store.traces.push(active_trace);
        }
        self.persist_routine_decision_trace_store_to_state(store);
    }

    pub(super) fn ensure_routine_assistant_decision_trace_started(&mut self) {
        if self.routine_assistant_decision_trace.is_some() {
            return;
        }
        let now = Self::now_unix_ms();
        let selected_routine = self.routine_assistant_selected_routine();
        let (routine_preference_context, candidate_planning_scores, macro_suggestions) =
            self.routine_assistant_planning_trace_artifacts(selected_routine.as_ref());
        self.routine_assistant_preference_context = routine_preference_context.clone();
        self.routine_assistant_macro_suggestions = macro_suggestions.clone();
        let trace = RoutineDecisionTrace {
            schema: ROUTINE_DECISION_TRACE_SCHEMA.to_string(),
            trace_id: self.next_routine_assistant_trace_id(),
            source: "gui_routine_assistant".to_string(),
            status: "draft".to_string(),
            created_at_unix_ms: now,
            updated_at_unix_ms: now,
            goal_text: self.routine_assistant_goal.trim().to_string(),
            query_text: self.routine_assistant_query.trim().to_string(),
            candidate_routine_ids: self.routine_assistant_candidate_ids_snapshot(),
            routine_preference_context,
            candidate_planning_scores,
            macro_suggestions,
            ..RoutineDecisionTrace::default()
        };
        self.routine_assistant_decision_trace = Some(trace);
        self.persist_routine_assistant_decision_trace();
    }

    pub(super) fn update_routine_assistant_decision_trace<F>(&mut self, updater: F)
    where
        F: FnOnce(&mut RoutineDecisionTrace),
    {
        self.ensure_routine_assistant_decision_trace_started();
        let goal_text = self.routine_assistant_goal.trim().to_string();
        let query_text = self.routine_assistant_query.trim().to_string();
        let candidate_routine_ids = self.routine_assistant_candidate_ids_snapshot();
        let selected_routine = self.routine_assistant_selected_routine();
        let (routine_preference_context, candidate_planning_scores, macro_suggestions) =
            self.routine_assistant_planning_trace_artifacts(selected_routine.as_ref());
        self.routine_assistant_preference_context = routine_preference_context.clone();
        self.routine_assistant_macro_suggestions = macro_suggestions.clone();
        let now = Self::now_unix_ms();
        if let Some(trace) = self.routine_assistant_decision_trace.as_mut() {
            updater(trace);
            trace.goal_text = goal_text;
            trace.query_text = query_text;
            trace.candidate_routine_ids = candidate_routine_ids;
            trace.routine_preference_context = routine_preference_context;
            trace.candidate_planning_scores = candidate_planning_scores;
            trace.macro_suggestions = macro_suggestions;
            trace.updated_at_unix_ms = now;
        }
        self.persist_routine_assistant_decision_trace();
    }

    pub(super) fn maybe_mark_routine_assistant_trace_aborted(&mut self) {
        let should_mark = self
            .routine_assistant_decision_trace
            .as_ref()
            .map(|trace| {
                !matches!(
                    trace.status.as_str(),
                    "executed" | "execution_failed" | "aborted" | "exported"
                )
            })
            .unwrap_or(false);
        if !should_mark {
            return;
        }
        self.update_routine_assistant_decision_trace(|trace| {
            trace.status = "aborted".to_string();
        });
    }

    pub(super) fn routine_assistant_capture_selected_routine(
        trace: &mut RoutineDecisionTrace,
        routine: Option<&CloningRoutineCatalogRow>,
    ) {
        if let Some(routine) = routine {
            trace.selected_routine_id = Some(routine.routine_id.trim().to_string());
            trace.selected_routine_title = Some(routine.title.trim().to_string());
            trace.selected_routine_family = Some(routine.family.trim().to_string());
        }
    }

    pub(super) fn routine_assistant_preflight_snapshot_from_output(
        output: &serde_json::Value,
    ) -> Option<RoutineDecisionTracePreflightSnapshot> {
        let can_execute = output
            .get("can_execute")
            .and_then(|value| value.as_bool())?;
        let preflight = output.get("preflight")?;
        let warnings = preflight
            .get("warnings")
            .and_then(|value| value.as_array())
            .map(|rows| {
                rows.iter()
                    .filter_map(|row| row.as_str())
                    .map(str::trim)
                    .filter(|row| !row.is_empty())
                    .map(|row| row.to_string())
                    .collect::<Vec<_>>()
            })
            .unwrap_or_default();
        let errors = preflight
            .get("errors")
            .and_then(|value| value.as_array())
            .map(|rows| {
                rows.iter()
                    .filter_map(|row| row.as_str())
                    .map(str::trim)
                    .filter(|row| !row.is_empty())
                    .map(|row| row.to_string())
                    .collect::<Vec<_>>()
            })
            .unwrap_or_default();
        let contract_source = preflight
            .get("contract_source")
            .and_then(|value| value.as_str())
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_string());
        Some(RoutineDecisionTracePreflightSnapshot {
            can_execute,
            warnings,
            errors,
            contract_source,
        })
    }

    pub(super) fn collect_op_ids_from_json(value: &serde_json::Value, op_ids: &mut Vec<String>) {
        match value {
            serde_json::Value::Object(map) => {
                if let Some(op_id) = map.get("op_id").and_then(|value| value.as_str()) {
                    Self::push_unique_trace_token(op_ids, op_id);
                }
                for nested in map.values() {
                    Self::collect_op_ids_from_json(nested, op_ids);
                }
            }
            serde_json::Value::Array(rows) => {
                for row in rows {
                    Self::collect_op_ids_from_json(row, op_ids);
                }
            }
            _ => {}
        }
    }

    pub(super) fn routine_assistant_emitted_op_ids_from_execute_output(
        output: &serde_json::Value,
    ) -> Vec<String> {
        let mut op_ids: Vec<String> = vec![];
        if let Some(run) = output.get("run") {
            Self::collect_op_ids_from_json(run, &mut op_ids);
        } else {
            Self::collect_op_ids_from_json(output, &mut op_ids);
        }
        op_ids
    }

    pub(super) fn list_cloning_routines(
        &mut self,
        family: Option<&str>,
        status: Option<&str>,
        query: Option<&str>,
    ) -> std::result::Result<Vec<CloningRoutineCatalogRow>, String> {
        let command = ShellCommand::RoutinesList {
            catalog_path: Some(DEFAULT_CLONING_ROUTINE_CATALOG_PATH.to_string()),
            family: family
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(|value| value.to_string()),
            status: status
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(|value| value.to_string()),
            tag: None,
            seq_id: self.routine_assistant_construct_reasoning_seq_id(),
            query: query
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(|value| value.to_string()),
        };
        let (output, _) = self.execute_shared_shell_command_json(&command)?;
        let routines_json = output
            .get("routines")
            .cloned()
            .unwrap_or_else(|| serde_json::json!([]));
        serde_json::from_value::<Vec<CloningRoutineCatalogRow>>(routines_json)
            .map_err(|e| format!("Could not parse routine catalog output: {e}"))
    }

    pub(super) fn refresh_routine_assistant_candidates(&mut self) {
        let query = self.routine_assistant_query.trim().to_string();
        let fallback = self.routine_assistant_goal.trim().to_string();
        let effective_query = if query.is_empty() { fallback } else { query };
        match self.list_cloning_routines(
            None,
            None,
            if effective_query.is_empty() {
                None
            } else {
                Some(effective_query.as_str())
            },
        ) {
            Ok(rows) => {
                self.routine_assistant_candidates = rows;
                if !self.routine_assistant_selected_routine_id.trim().is_empty()
                    && !self.routine_assistant_candidates.iter().any(|row| {
                        row.routine_id
                            .eq_ignore_ascii_case(self.routine_assistant_selected_routine_id.trim())
                    })
                {
                    self.routine_assistant_selected_routine_id.clear();
                    self.routine_assistant_compare_routine_id.clear();
                    self.routine_assistant_bindings.clear();
                    self.routine_assistant_disambiguation_answers.clear();
                    self.routine_assistant_explain_output = None;
                    self.routine_assistant_compare_output = None;
                    self.routine_assistant_preflight_output = None;
                    self.routine_assistant_execute_output = None;
                }
                self.routine_assistant_status = format!(
                    "Routine Assistant: loaded {} candidate routine(s)",
                    self.routine_assistant_candidates.len()
                );
                let selected = self.routine_assistant_selected_routine();
                let bindings_snapshot = self.routine_assistant_bindings_snapshot();
                self.update_routine_assistant_decision_trace(|trace| {
                    trace.status = "draft".to_string();
                    trace.bindings_snapshot = bindings_snapshot;
                    Self::routine_assistant_capture_selected_routine(trace, selected.as_ref());
                });
            }
            Err(err) => {
                self.routine_assistant_candidates.clear();
                self.routine_assistant_status =
                    format!("Routine Assistant: could not list routines: {err}");
                self.update_routine_assistant_decision_trace(|trace| {
                    trace.status = "draft".to_string();
                });
            }
        }
    }

    pub(super) fn routine_assistant_selected_routine(&self) -> Option<CloningRoutineCatalogRow> {
        let selected_id = self.routine_assistant_selected_routine_id.trim();
        if selected_id.is_empty() {
            return None;
        }
        self.routine_assistant_candidates
            .iter()
            .find(|row| row.routine_id.eq_ignore_ascii_case(selected_id))
            .cloned()
    }

    pub(super) fn routine_assistant_input_port_ids(
        routine: &CloningRoutineCatalogRow,
    ) -> Vec<String> {
        routine
            .input_ports
            .iter()
            .filter_map(|port| {
                port.get("port_id")
                    .and_then(|value| value.as_str())
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .map(|value| value.to_string())
            })
            .collect::<Vec<_>>()
    }

    pub(super) fn routine_assistant_sequence_port_ids(
        routine: &CloningRoutineCatalogRow,
    ) -> Vec<String> {
        routine
            .input_ports
            .iter()
            .filter_map(|port| {
                let kind = port
                    .get("kind")
                    .and_then(|value| value.as_str())
                    .map(str::trim)
                    .unwrap_or("");
                if !kind.eq_ignore_ascii_case("sequence") {
                    return None;
                }
                port.get("port_id")
                    .and_then(|value| value.as_str())
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .map(|value| value.to_string())
            })
            .collect::<Vec<_>>()
    }

    pub(super) fn routine_assistant_bound_sequence_topologies_for_routine(
        &self,
        routine: &CloningRoutineCatalogRow,
    ) -> Vec<RoutineAssistantBoundSequenceTopology> {
        let sequence_ports = Self::routine_assistant_sequence_port_ids(routine)
            .into_iter()
            .collect::<HashSet<_>>();
        if sequence_ports.is_empty() {
            return vec![];
        }
        let Ok(engine) = self.engine.read() else {
            return vec![];
        };
        self.routine_assistant_bindings
            .iter()
            .filter_map(|(port_id, value)| {
                if !sequence_ports.contains(port_id) {
                    return None;
                }
                let seq_id = value.trim();
                if seq_id.is_empty() {
                    return None;
                }
                let dna = engine.state().sequences.get(seq_id)?;
                Some(RoutineAssistantBoundSequenceTopology {
                    port_id: port_id.clone(),
                    seq_id: seq_id.to_string(),
                    circular: dna.is_circular(),
                    length_bp: dna.len(),
                })
            })
            .collect::<Vec<_>>()
    }

    pub(super) fn routine_assistant_sequence_topology_for_seq_id(
        &self,
        seq_id: &str,
    ) -> Option<(bool, usize)> {
        let compact = seq_id.trim();
        if compact.is_empty() {
            return None;
        }
        let engine = self.engine.read().ok()?;
        let dna = engine.state().sequences.get(compact)?;
        Some((dna.is_circular(), dna.len()))
    }

    pub(super) fn routine_assistant_is_gibson_family(routine: &CloningRoutineCatalogRow) -> bool {
        if routine.family.trim().eq_ignore_ascii_case("gibson") {
            return true;
        }
        let routine_id = routine.routine_id.to_ascii_lowercase();
        let template_name = routine.template_name.to_ascii_lowercase();
        routine_id.contains("gibson") || template_name.contains("gibson")
    }

    pub(super) fn routine_assistant_gibson_circular_blocking_error(
        binding: &RoutineAssistantBoundSequenceTopology,
    ) -> String {
        format!(
            "Gibson requires linear fragments: binding '{}' on port '{}' is circular ({} bp).",
            binding.seq_id, binding.port_id, binding.length_bp
        )
    }

    pub(super) fn routine_assistant_gibson_circular_blocking_preflight_output(
        &self,
        routine: &CloningRoutineCatalogRow,
        binding: &RoutineAssistantBoundSequenceTopology,
    ) -> serde_json::Value {
        let error = Self::routine_assistant_gibson_circular_blocking_error(binding);
        serde_json::json!({
            "schema": "gentle.macro_template_preflight.v1",
            "can_execute": false,
            "routine_id": routine.routine_id,
            "template_name": routine.template_name,
            "preflight": {
                "contract_source": "routine_assistant.gibson_linearization_guard.v1",
                "errors": [error],
                "warnings": [
                    "Use 'Linearize Vector...' to create a linear branch before preflight/execute."
                ]
            }
        })
    }

    pub(super) fn routine_assistant_gibson_circular_binding_for_routine(
        &self,
        routine: &CloningRoutineCatalogRow,
    ) -> Option<RoutineAssistantBoundSequenceTopology> {
        if !Self::routine_assistant_is_gibson_family(routine) {
            return None;
        }
        let circular_inputs = self
            .routine_assistant_bound_sequence_topologies_for_routine(routine)
            .into_iter()
            .filter(|row| row.circular)
            .collect::<Vec<_>>();
        if circular_inputs.is_empty() {
            return None;
        }
        for preferred_port in ["vector_seq_id", "backbone_seq_id", "right_seq_id"] {
            if let Some(row) = circular_inputs
                .iter()
                .find(|row| row.port_id.eq_ignore_ascii_case(preferred_port))
            {
                return Some(row.clone());
            }
        }
        if let Some(row) = circular_inputs.iter().find(|row| {
            let lower = row.port_id.to_ascii_lowercase();
            lower.contains("vector") || lower.contains("backbone")
        }) {
            return Some(row.clone());
        }
        circular_inputs.into_iter().next()
    }

    pub(super) fn render_routine_assistant_planning_context_strip(&self, ui: &mut Ui) {
        let Some(context) = self.routine_assistant_preference_context.as_ref() else {
            return;
        };
        if context.helper_profile_id.is_none()
            && context.construct_reasoning_seq_id.is_none()
            && context.effective_preferred_routine_families.is_empty()
            && context.variant_effect_tags.is_empty()
            && context.variant_suggested_assay_ids.is_empty()
            && context.rationale.is_empty()
        {
            return;
        }
        ui.group(|ui| {
            ui.strong("Planning Context");
            if let Some(seq_id) = context.construct_reasoning_seq_id.as_deref() {
                ui.small(format!("construct reasoning: {seq_id}"));
            }
            if let Some(helper_profile_id) = context.helper_profile_id.as_deref() {
                ui.small(format!(
                    "helper profile: {} [{}]",
                    helper_profile_id, context.helper_resolution_status
                ));
            }
            if !context.effective_preferred_routine_families.is_empty() {
                ui.small(format!(
                    "preferred routine families: {}",
                    context.effective_preferred_routine_families.join(", ")
                ));
            }
            if !context.variant_effect_tags.is_empty() {
                ui.small(format!(
                    "variant effect tags: {}",
                    context.variant_effect_tags.join(", ")
                ));
            }
            if !context.variant_suggested_assay_ids.is_empty() {
                ui.small(format!(
                    "suggested variant assays: {}",
                    context.variant_suggested_assay_ids.join(", ")
                ));
            }
            if let Some(line) = context.rationale.first() {
                ui.small(line);
            }
        });
    }

    pub(super) fn render_routine_assistant_macro_suggestions(&self, ui: &mut Ui) {
        if self.routine_assistant_macro_suggestions.is_empty() {
            return;
        }
        ui.group(|ui| {
            ui.strong("Suggested Macros");
            for suggestion in &self.routine_assistant_macro_suggestions {
                ui.horizontal_wrapped(|ui| {
                    ui.label(format!(
                        "{}: {} (score {:.2})",
                        suggestion.macro_kind, suggestion.template_name, suggestion.score
                    ));
                    if let Some(details_url) = suggestion.details_url.as_deref() {
                        ui.hyperlink_to("docs", details_url);
                    }
                });
                if let Some(description) = suggestion.description.as_deref() {
                    ui.small(description);
                }
                if !suggestion.matched_routine_families.is_empty() {
                    ui.small(format!(
                        "matched families: {}",
                        suggestion.matched_routine_families.join(", ")
                    ));
                }
                if let Some(line) = suggestion.rationale.first() {
                    ui.small(line);
                }
                ui.add_space(4.0);
            }
        });
    }

    pub(super) fn render_routine_assistant_gibson_linearization_notice(
        &mut self,
        ui: &mut Ui,
        routine: &CloningRoutineCatalogRow,
    ) {
        let Some(binding) = self.routine_assistant_gibson_circular_binding_for_routine(routine)
        else {
            return;
        };
        let error = Self::routine_assistant_gibson_circular_blocking_error(&binding);
        ui.group(|ui| {
            ui.colored_label(egui::Color32::from_rgb(190, 70, 70), error);
            ui.small(
                "One-click fix: create a branched copy, force linear topology, and re-bind this input.",
            );
            if ui
                .button("Linearize Vector...")
                .on_hover_text(
                    "Create a branched linear copy of the bound circular sequence and rebind this Gibson input to that copy.",
                )
                .clicked()
            {
                let port_id = binding.port_id.clone();
                let seq_id = binding.seq_id.clone();
                match self.routine_assistant_linearize_binding_sequence(&port_id, &seq_id) {
                    Ok(new_id) => {
                        self.routine_assistant_status = format!(
                            "Routine Assistant: linearized '{}' as '{}' and rebound '{}'",
                            seq_id, new_id, port_id
                        );
                    }
                    Err(err) => {
                        self.routine_assistant_status =
                            format!("Routine Assistant linearization failed: {err}");
                    }
                }
            }
        });
    }

    pub(super) fn routine_assistant_linearize_binding_sequence(
        &mut self,
        port_id: &str,
        seq_id: &str,
    ) -> std::result::Result<String, String> {
        let compact_port = port_id.trim();
        if compact_port.is_empty() {
            return Err("Linearize Vector requires a non-empty binding port".to_string());
        }
        let compact_seq = seq_id.trim();
        if compact_seq.is_empty() {
            return Err("Linearize Vector requires a non-empty sequence ID".to_string());
        }
        let (exists, is_circular) = {
            let engine = self
                .engine
                .read()
                .map_err(|_| "Engine lock poisoned while checking sequence topology".to_string())?;
            match engine.state().sequences.get(compact_seq) {
                Some(dna) => (true, dna.is_circular()),
                None => (false, false),
            }
        };
        if !exists {
            return Err(format!(
                "Linearize Vector could not find sequence '{}'",
                compact_seq
            ));
        }
        if !is_circular {
            return Err(format!(
                "Sequence '{}' is already linear; no linearization needed",
                compact_seq
            ));
        }

        let suggested_id = format!("{}_linear", compact_seq);
        let branch_result = {
            let mut engine = self
                .engine
                .write()
                .map_err(|_| "Engine lock poisoned while branching sequence".to_string())?;
            engine
                .apply(Operation::Branch {
                    input: compact_seq.to_string(),
                    output_id: Some(suggested_id.clone()),
                })
                .map_err(|e| format!("Linearize Vector branch failed: {}", e.message))?
        };
        self.lineage_cache_valid = false;
        if branch_result.created_seq_ids.is_empty() {
            return Err("Linearize Vector branch operation did not change state".to_string());
        }
        let created_id = branch_result
            .created_seq_ids
            .first()
            .cloned()
            .unwrap_or(suggested_id);

        {
            let mut engine = self
                .engine
                .write()
                .map_err(|_| "Engine lock poisoned while updating topology".to_string())?;
            engine
                .apply(Operation::SetTopology {
                    seq_id: created_id.clone(),
                    circular: false,
                })
                .map_err(|e| {
                    format!(
                        "Linearize Vector could not set linear topology for '{}': {}",
                        created_id, e.message
                    )
                })?;
        }
        self.lineage_cache_valid = false;

        self.routine_assistant_bindings
            .insert(compact_port.to_string(), created_id.clone());
        self.routine_assistant_preflight_output = None;
        self.routine_assistant_execute_output = None;
        Ok(created_id)
    }

    pub(super) fn sync_routine_assistant_bindings_for_selected(&mut self) {
        let Some(routine) = self.routine_assistant_selected_routine() else {
            self.routine_assistant_bindings.clear();
            self.routine_assistant_disambiguation_answers.clear();
            return;
        };
        let allowed = Self::routine_assistant_input_port_ids(&routine)
            .into_iter()
            .collect::<HashSet<_>>();
        self.routine_assistant_bindings
            .retain(|key, _| allowed.contains(key));
        for key in allowed {
            self.routine_assistant_bindings.entry(key).or_default();
        }
    }

    pub(super) fn routine_assistant_bindings_compact(&self) -> HashMap<String, String> {
        self.routine_assistant_bindings
            .iter()
            .filter_map(|(key, value)| {
                let compact = value.trim();
                if compact.is_empty() {
                    None
                } else {
                    Some((key.clone(), compact.to_string()))
                }
            })
            .collect::<HashMap<_, _>>()
    }

    pub(super) fn load_routine_assistant_explain(&mut self) {
        let selected_id = self
            .routine_assistant_selected_routine_id
            .trim()
            .to_string();
        if selected_id.is_empty() {
            self.routine_assistant_status =
                "Routine Assistant: select a primary routine first".to_string();
            return;
        }
        let command = ShellCommand::RoutinesExplain {
            catalog_path: Some(DEFAULT_CLONING_ROUTINE_CATALOG_PATH.to_string()),
            routine_id: selected_id.clone(),
            seq_id: self.routine_assistant_construct_reasoning_seq_id(),
        };
        match self.execute_shared_shell_command_json(&command) {
            Ok((output, _)) => {
                self.routine_assistant_explain_output = Some(output.clone());
                if self.routine_assistant_compare_routine_id.trim().is_empty() {
                    if let Some(alt_id) = output
                        .get("alternatives")
                        .and_then(|value| value.as_array())
                        .and_then(|rows| rows.first())
                        .and_then(|row| row.get("routine_id"))
                        .and_then(|value| value.as_str())
                    {
                        self.routine_assistant_compare_routine_id = alt_id.trim().to_string();
                    }
                }
                self.routine_assistant_status =
                    format!("Routine Assistant: loaded explanation for '{selected_id}'");
                let selected = self.routine_assistant_selected_routine();
                let mut alternatives: Vec<String> = vec![];
                if let Some(rows) = output
                    .get("alternatives")
                    .and_then(|value| value.as_array())
                {
                    for row in rows {
                        if let Some(routine_id) =
                            row.get("routine_id").and_then(|value| value.as_str())
                        {
                            Self::push_unique_trace_token(&mut alternatives, routine_id);
                        }
                    }
                }
                let bindings_snapshot = self.routine_assistant_bindings_snapshot();
                let disambiguation_questions =
                    Self::routine_assistant_disambiguation_questions_from_output(&output);
                self.sync_routine_assistant_disambiguation_answers_for_questions(
                    &disambiguation_questions,
                    &[],
                );
                let disambiguation_answers = self
                    .routine_assistant_disambiguation_answers_snapshot(&disambiguation_questions);
                self.update_routine_assistant_decision_trace(|trace| {
                    trace.status = "draft".to_string();
                    trace.alternatives_presented = alternatives;
                    trace.disambiguation_questions_presented = disambiguation_questions.clone();
                    trace.disambiguation_answers = disambiguation_answers;
                    trace.bindings_snapshot = bindings_snapshot;
                    Self::routine_assistant_capture_selected_routine(trace, selected.as_ref());
                });
            }
            Err(err) => {
                self.routine_assistant_status = format!("Routine Assistant explain failed: {err}");
                self.update_routine_assistant_decision_trace(|trace| {
                    trace.status = "draft".to_string();
                });
            }
        }
    }

    pub(super) fn load_routine_assistant_compare(&mut self) {
        let left = self
            .routine_assistant_selected_routine_id
            .trim()
            .to_string();
        let right = self.routine_assistant_compare_routine_id.trim().to_string();
        if left.is_empty() || right.is_empty() {
            self.routine_assistant_status =
                "Routine Assistant: select both primary and comparison routines".to_string();
            return;
        }
        let command = ShellCommand::RoutinesCompare {
            catalog_path: Some(DEFAULT_CLONING_ROUTINE_CATALOG_PATH.to_string()),
            left_routine_id: left.clone(),
            right_routine_id: right.clone(),
            seq_id: self.routine_assistant_construct_reasoning_seq_id(),
        };
        match self.execute_shared_shell_command_json(&command) {
            Ok((output, _)) => {
                self.routine_assistant_compare_output = Some(output);
                self.routine_assistant_status =
                    format!("Routine Assistant: compared '{left}' vs '{right}'");
                let selected = self.routine_assistant_selected_routine();
                let disambiguation_questions =
                    self.routine_assistant_effective_disambiguation_questions();
                self.sync_routine_assistant_disambiguation_answers_for_questions(
                    &disambiguation_questions,
                    &[],
                );
                let disambiguation_answers = self
                    .routine_assistant_disambiguation_answers_snapshot(&disambiguation_questions);
                self.update_routine_assistant_decision_trace(|trace| {
                    trace.status = "draft".to_string();
                    Self::routine_assistant_capture_selected_routine(trace, selected.as_ref());
                    Self::merge_routine_assistant_disambiguation_questions(
                        &mut trace.disambiguation_questions_presented,
                        disambiguation_questions.clone(),
                    );
                    trace.disambiguation_answers = disambiguation_answers;
                    if !trace.comparisons.iter().any(|row| {
                        row.left_routine_id.eq_ignore_ascii_case(&left)
                            && row.right_routine_id.eq_ignore_ascii_case(&right)
                    }) {
                        trace.comparisons.push(RoutineDecisionTraceComparison {
                            left_routine_id: left.clone(),
                            right_routine_id: right.clone(),
                        });
                    }
                });
            }
            Err(err) => {
                self.routine_assistant_status = format!("Routine Assistant compare failed: {err}");
                self.update_routine_assistant_decision_trace(|trace| {
                    trace.status = "draft".to_string();
                });
            }
        }
    }

    pub(super) fn ensure_routine_assistant_template_imported(
        &mut self,
        routine: &CloningRoutineCatalogRow,
    ) -> std::result::Result<(), String> {
        let Some(path) = routine
            .template_path
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
        else {
            return Err(format!(
                "Routine '{}' has no template_path configured",
                routine.routine_id
            ));
        };
        let command = ShellCommand::MacrosTemplateImport {
            path: path.to_string(),
        };
        self.execute_shared_shell_command_json(&command).map(|_| ())
    }

    pub(super) fn run_routine_assistant_preflight(&mut self) {
        let Some(routine) = self.routine_assistant_selected_routine() else {
            self.routine_assistant_status =
                "Routine Assistant: select a routine before preflight".to_string();
            return;
        };
        if let Some(binding) = self.routine_assistant_gibson_circular_binding_for_routine(&routine)
        {
            self.routine_assistant_preflight_output = Some(
                self.routine_assistant_gibson_circular_blocking_preflight_output(
                    &routine, &binding,
                ),
            );
            self.routine_assistant_execute_output = None;
            self.routine_assistant_stage = RoutineAssistantStage::Preflight;
            self.routine_assistant_status = format!(
                "Routine Assistant preflight blocked: '{}' on '{}' is circular",
                binding.seq_id, binding.port_id
            );
            let preflight_snapshot = self
                .routine_assistant_preflight_output
                .as_ref()
                .and_then(Self::routine_assistant_preflight_snapshot_from_output);
            let bindings_snapshot = self.routine_assistant_bindings_snapshot();
            self.update_routine_assistant_decision_trace(|trace| {
                trace.status = "preflight_failed".to_string();
                trace.bindings_snapshot = bindings_snapshot;
                Self::routine_assistant_commit_preflight_snapshot(trace, preflight_snapshot);
                Self::routine_assistant_capture_selected_routine(trace, Some(&routine));
            });
            return;
        }
        if let Err(err) = self.ensure_routine_assistant_template_imported(&routine) {
            self.routine_assistant_status = format!("Routine Assistant preflight failed: {err}");
            let bindings_snapshot = self.routine_assistant_bindings_snapshot();
            self.update_routine_assistant_decision_trace(|trace| {
                trace.status = "preflight_failed".to_string();
                trace.bindings_snapshot = bindings_snapshot;
                Self::routine_assistant_commit_preflight_snapshot(trace, None);
                trace.execution_error = Some(err.clone());
                Self::routine_assistant_capture_selected_routine(trace, Some(&routine));
            });
            return;
        }
        let command = ShellCommand::MacrosTemplateRun {
            name: routine.template_name.clone(),
            bindings: self.routine_assistant_bindings_compact(),
            transactional: false,
            validate_only: true,
        };
        match self.execute_shared_shell_command_json(&command) {
            Ok((output, _)) => {
                self.routine_assistant_preflight_output = Some(output.clone());
                self.routine_assistant_execute_output = None;
                self.routine_assistant_stage = RoutineAssistantStage::Preflight;
                let can_execute = output
                    .get("can_execute")
                    .and_then(|value| value.as_bool())
                    .unwrap_or(false);
                self.routine_assistant_status = if can_execute {
                    "Routine Assistant: preflight passed".to_string()
                } else {
                    "Routine Assistant: preflight reported blocking errors".to_string()
                };
                let preflight_snapshot =
                    Self::routine_assistant_preflight_snapshot_from_output(&output);
                let bindings_snapshot = self.routine_assistant_bindings_snapshot();
                self.update_routine_assistant_decision_trace(|trace| {
                    trace.status = if can_execute {
                        "ready".to_string()
                    } else {
                        "preflight_failed".to_string()
                    };
                    trace.bindings_snapshot = bindings_snapshot;
                    Self::routine_assistant_commit_preflight_snapshot(trace, preflight_snapshot);
                    trace.execution_error = None;
                    Self::routine_assistant_capture_selected_routine(trace, Some(&routine));
                });
            }
            Err(err) => {
                self.routine_assistant_status =
                    format!("Routine Assistant preflight failed: {err}");
                let bindings_snapshot = self.routine_assistant_bindings_snapshot();
                self.update_routine_assistant_decision_trace(|trace| {
                    trace.status = "preflight_failed".to_string();
                    trace.bindings_snapshot = bindings_snapshot;
                    Self::routine_assistant_commit_preflight_snapshot(trace, None);
                    trace.execution_error = Some(err.clone());
                    Self::routine_assistant_capture_selected_routine(trace, Some(&routine));
                });
            }
        }
    }

    pub(super) fn run_routine_assistant_execute(&mut self) {
        let Some(routine) = self.routine_assistant_selected_routine() else {
            self.routine_assistant_status =
                "Routine Assistant: select a routine before execution".to_string();
            return;
        };
        if let Some(binding) = self.routine_assistant_gibson_circular_binding_for_routine(&routine)
        {
            self.routine_assistant_preflight_output = Some(
                self.routine_assistant_gibson_circular_blocking_preflight_output(
                    &routine, &binding,
                ),
            );
            self.routine_assistant_execute_output = None;
            self.routine_assistant_stage = RoutineAssistantStage::Preflight;
            self.routine_assistant_status = format!(
                "Routine Assistant execution blocked: '{}' on '{}' is circular",
                binding.seq_id, binding.port_id
            );
            let preflight_snapshot = self
                .routine_assistant_preflight_output
                .as_ref()
                .and_then(Self::routine_assistant_preflight_snapshot_from_output);
            let bindings_snapshot = self.routine_assistant_bindings_snapshot();
            self.update_routine_assistant_decision_trace(|trace| {
                trace.status = "preflight_failed".to_string();
                trace.bindings_snapshot = bindings_snapshot;
                Self::routine_assistant_commit_preflight_snapshot(trace, preflight_snapshot);
                trace.execution_attempted = false;
                trace.execution_success = Some(false);
                trace.transactional = Some(true);
                trace.execution_error = Some("execution blocked by preflight guard".to_string());
                Self::routine_assistant_capture_selected_routine(trace, Some(&routine));
            });
            return;
        }
        if let Err(err) = self.ensure_routine_assistant_template_imported(&routine) {
            self.routine_assistant_status = format!("Routine Assistant execution failed: {err}");
            let bindings_snapshot = self.routine_assistant_bindings_snapshot();
            self.update_routine_assistant_decision_trace(|trace| {
                trace.status = "execution_failed".to_string();
                trace.bindings_snapshot = bindings_snapshot;
                trace.execution_attempted = false;
                trace.execution_success = Some(false);
                trace.transactional = Some(true);
                trace.execution_error = Some(err.clone());
                Self::routine_assistant_capture_selected_routine(trace, Some(&routine));
            });
            return;
        }
        let command = ShellCommand::MacrosTemplateRun {
            name: routine.template_name.clone(),
            bindings: self.routine_assistant_bindings_compact(),
            transactional: true,
            validate_only: false,
        };
        match self.execute_shared_shell_command_json(&command) {
            Ok((output, _)) => {
                let emitted_operation_ids =
                    Self::routine_assistant_emitted_op_ids_from_execute_output(&output);
                let macro_instance_id = output
                    .get("macro_instance_id")
                    .and_then(|value| value.as_str())
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .map(|value| value.to_string());
                let preflight_snapshot =
                    Self::routine_assistant_preflight_snapshot_from_output(&output);
                self.routine_assistant_execute_output = Some(output);
                self.routine_assistant_stage = RoutineAssistantStage::ExecuteAndExport;
                self.routine_assistant_status =
                    "Routine Assistant: transactional run completed".to_string();
                let bindings_snapshot = self.routine_assistant_bindings_snapshot();
                self.update_routine_assistant_decision_trace(|trace| {
                    trace.status = "executed".to_string();
                    trace.bindings_snapshot = bindings_snapshot;
                    Self::routine_assistant_commit_preflight_snapshot(trace, preflight_snapshot);
                    trace.execution_attempted = true;
                    trace.execution_success = Some(true);
                    trace.transactional = Some(true);
                    trace.macro_instance_id = macro_instance_id;
                    trace.emitted_operation_ids = emitted_operation_ids;
                    trace.execution_error = None;
                    Self::routine_assistant_capture_selected_routine(trace, Some(&routine));
                });
            }
            Err(err) => {
                self.routine_assistant_status =
                    format!("Routine Assistant execution failed: {err}");
                let bindings_snapshot = self.routine_assistant_bindings_snapshot();
                self.update_routine_assistant_decision_trace(|trace| {
                    trace.status = "execution_failed".to_string();
                    trace.bindings_snapshot = bindings_snapshot;
                    trace.execution_attempted = true;
                    trace.execution_success = Some(false);
                    trace.transactional = Some(true);
                    trace.execution_error = Some(err.clone());
                    Self::routine_assistant_capture_selected_routine(trace, Some(&routine));
                });
            }
        }
    }

    pub(super) fn export_routine_assistant_run_bundle(&mut self) {
        let Some(path) = rfd::FileDialog::new()
            .set_file_name("run_bundle.routine_assistant.json")
            .add_filter("JSON", &["json"])
            .save_file()
        else {
            self.routine_assistant_status =
                "Routine Assistant: run-bundle export canceled".to_string();
            return;
        };
        let path_text = path.display().to_string();
        let command = ShellCommand::ExportRunBundle {
            output: path_text.clone(),
            run_id: None,
        };
        match self.execute_shared_shell_command_json(&command) {
            Ok(_) => {
                self.routine_assistant_status =
                    format!("Routine Assistant: exported run bundle to '{path_text}'");
                self.update_routine_assistant_decision_trace(|trace| {
                    trace.status = "exported".to_string();
                    trace.export_events.push(RoutineDecisionTraceExportEvent {
                        run_bundle_path: path_text.clone(),
                        exported_at_unix_ms: Self::now_unix_ms(),
                    });
                });
            }
            Err(err) => {
                self.routine_assistant_status = format!(
                    "Routine Assistant: could not export run bundle '{}': {}",
                    path_text, err
                );
                self.update_routine_assistant_decision_trace(|trace| {
                    trace.status = "execution_failed".to_string();
                    trace.execution_error = Some(format!("run-bundle export failed: {}", err));
                });
            }
        }
    }

    pub(super) fn routine_assistant_can_execute(&self) -> bool {
        self.routine_assistant_preflight_output
            .as_ref()
            .and_then(|value| value.get("can_execute"))
            .and_then(|value| value.as_bool())
            .unwrap_or(false)
    }
}
