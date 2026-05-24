//! External-service GUI workspace backed by shared shell contracts.
//!
//! The GUI intentionally does not know vendor business logic. It renders the
//! active provider catalog, lets users edit one provider-neutral request JSON,
//! and calls the same `services ...` shell commands used by CLI/agent/ClawBio.

use super::*;
use serde_json::json;

#[derive(Clone, Debug)]
struct ExternalServiceCapabilityChoice {
    service_kind: String,
    display_name: String,
}

#[derive(Clone, Debug)]
struct ExternalServiceProviderChoice {
    provider: String,
    display_name: String,
    support_status: String,
    capabilities: Vec<ExternalServiceCapabilityChoice>,
}

impl GENtleApp {
    pub(super) fn open_external_services_dialog(&mut self) {
        let was_open = self.external_services_ui.show_panel;
        self.external_services_ui.show_panel = true;
        if self.external_services_ui.provider_catalog_output.is_none() {
            self.refresh_external_services_provider_catalog();
        }
        if self.external_services_ui.request_json.trim().is_empty() {
            self.reset_external_services_request_from_selection();
        }
        if self.external_services_ui.quote_output_dir.trim().is_empty() {
            self.external_services_ui.quote_output_dir =
                self.default_external_service_quote_output_dir();
        }
        self.mark_window_open_or_focus(Self::external_services_viewport_id(), was_open);
    }

    fn external_service_provider_choices(&self) -> Vec<ExternalServiceProviderChoice> {
        let Some(catalog) = self.external_services_ui.provider_catalog_output.as_ref() else {
            return vec![];
        };
        catalog
            .get("providers")
            .and_then(serde_json::Value::as_array)
            .into_iter()
            .flatten()
            .filter_map(|provider| {
                let provider_id = provider.get("provider")?.as_str()?.to_string();
                let display_name = provider
                    .get("display_name")
                    .and_then(serde_json::Value::as_str)
                    .unwrap_or(provider_id.as_str())
                    .to_string();
                let support_status = provider
                    .get("support_status")
                    .and_then(serde_json::Value::as_str)
                    .unwrap_or_default()
                    .to_string();
                let capabilities = provider
                    .get("capabilities")
                    .and_then(serde_json::Value::as_array)
                    .into_iter()
                    .flatten()
                    .filter_map(|capability| {
                        let service_kind = capability.get("service_kind")?.as_str()?.to_string();
                        let display_name = capability
                            .get("display_name")
                            .and_then(serde_json::Value::as_str)
                            .unwrap_or(service_kind.as_str())
                            .to_string();
                        Some(ExternalServiceCapabilityChoice {
                            service_kind,
                            display_name,
                        })
                    })
                    .collect::<Vec<_>>();
                Some(ExternalServiceProviderChoice {
                    provider: provider_id,
                    display_name,
                    support_status,
                    capabilities,
                })
            })
            .collect()
    }

    fn selected_external_service_provider_choice(&self) -> Option<ExternalServiceProviderChoice> {
        let selected = self.external_services_ui.selected_provider.trim();
        self.external_service_provider_choices()
            .into_iter()
            .find(|provider| provider.provider == selected)
            .or_else(|| self.external_service_provider_choices().into_iter().next())
    }

    fn selected_external_service_capability_choice(
        &self,
    ) -> Option<ExternalServiceCapabilityChoice> {
        let selected = self.external_services_ui.selected_service_kind.trim();
        self.selected_external_service_provider_choice()
            .and_then(|provider| {
                provider
                    .capabilities
                    .iter()
                    .find(|capability| capability.service_kind == selected)
                    .cloned()
                    .or_else(|| provider.capabilities.first().cloned())
            })
    }

    fn safe_external_service_output_segment(value: &str) -> String {
        let mut out = String::new();
        let mut previous_was_sep = false;
        for ch in value.chars() {
            let mapped = if ch.is_ascii_alphanumeric() {
                ch.to_ascii_lowercase()
            } else {
                '_'
            };
            if mapped == '_' {
                if previous_was_sep {
                    continue;
                }
                previous_was_sep = true;
            } else {
                previous_was_sep = false;
            }
            out.push(mapped);
        }
        let trimmed = out.trim_matches('_');
        if trimmed.is_empty() {
            "handoff".to_string()
        } else {
            trimmed.to_string()
        }
    }

    fn default_external_service_quote_output_dir(&self) -> String {
        let provider = Self::safe_external_service_output_segment(
            self.external_services_ui.selected_provider.as_str(),
        );
        let service_kind = Self::safe_external_service_output_segment(
            self.external_services_ui.selected_service_kind.as_str(),
        );
        format!("artifacts/external_services/{provider}_{service_kind}_handoff")
    }

    fn ensure_external_services_selection_from_catalog(&mut self) {
        let providers = self.external_service_provider_choices();
        if providers.is_empty() {
            return;
        }
        let provider_exists = providers
            .iter()
            .any(|provider| provider.provider == self.external_services_ui.selected_provider);
        if !provider_exists {
            self.external_services_ui.selected_provider = providers[0].provider.clone();
        }
        let capabilities = providers
            .iter()
            .find(|provider| provider.provider == self.external_services_ui.selected_provider)
            .map(|provider| provider.capabilities.clone())
            .unwrap_or_default();
        if capabilities.is_empty() {
            self.external_services_ui.selected_service_kind.clear();
            return;
        }
        let capability_exists = capabilities.iter().any(|capability| {
            capability.service_kind == self.external_services_ui.selected_service_kind
        });
        if !capability_exists {
            self.external_services_ui.selected_service_kind = capabilities[0].service_kind.clone();
        }
    }

    fn refresh_external_services_provider_catalog(&mut self) {
        let command = ShellCommand::ServicesProvidersList;
        match self.execute_shared_shell_command_json(&command) {
            Ok((output, _)) => {
                let provider_count = output
                    .get("providers")
                    .and_then(serde_json::Value::as_array)
                    .map(Vec::len)
                    .unwrap_or(0);
                self.external_services_ui.provider_catalog_output = Some(output);
                self.ensure_external_services_selection_from_catalog();
                self.external_services_ui.status =
                    format!("Loaded {provider_count} external-service provider(s)");
            }
            Err(err) => {
                self.external_services_ui.status =
                    format!("Could not load external-service provider catalog: {err}");
            }
        }
    }

    fn refresh_external_services_provider_doctor(&mut self) {
        let command = ShellCommand::ServicesProvidersDoctor {
            catalog_path: None,
            output: None,
        };
        match self.execute_shared_shell_command_json(&command) {
            Ok((output, _)) => {
                let errors = output
                    .get("error_count")
                    .and_then(serde_json::Value::as_u64)
                    .unwrap_or(0);
                let warnings = output
                    .get("warning_count")
                    .and_then(serde_json::Value::as_u64)
                    .unwrap_or(0);
                self.external_services_ui.provider_doctor_output = Some(output);
                self.external_services_ui.status = format!(
                    "Provider config doctor complete: {errors} error(s), {warnings} warning(s)"
                );
            }
            Err(err) => {
                self.external_services_ui.status =
                    format!("External-service provider doctor failed: {err}");
            }
        }
    }

    fn default_external_service_request_json(&self) -> String {
        let provider = self
            .selected_external_service_provider_choice()
            .map(|row| row.provider)
            .unwrap_or_else(|| self.external_services_ui.selected_provider.clone());
        let service_kind = self
            .selected_external_service_capability_choice()
            .map(|row| row.service_kind)
            .unwrap_or_else(|| self.external_services_ui.selected_service_kind.clone());
        let lower_service = service_kind.to_ascii_lowercase();
        let source_target = if lower_service.contains("protein") {
            json!({
                "kind": "inline_protein",
                "name": "protein_candidate",
                "protein_sequence": "MSEQUENCE"
            })
        } else if lower_service.contains("fragment") {
            json!({
                "kind": "inline_dna",
                "name": "fragment_1",
                "sequence": "ATGGCTGCTGCTTAA"
            })
        } else {
            json!({
                "kind": "inline_dna",
                "name": "oligo_1",
                "sequence": "ACGTACGTACGTACGTACGT"
            })
        };
        let request = json!({
            "schema": "gentle.external_service_request.v1",
            "provider": provider,
            "service_kind": service_kind,
            "source_target": source_target,
            "delivery_options": {
                "delivery_form": "dry",
                "purification": "confirm_before_submission"
            },
            "return_spec": {
                "requested_payloads": ["quote_metadata", "handoff_bundle"],
                "redact_commercial_fields": true,
                "prefer_artifact_bundle": true
            }
        });
        serde_json::to_string_pretty(&request).unwrap_or_else(|_| "{}".to_string())
    }

    pub(super) fn reset_external_services_request_from_selection(&mut self) {
        self.ensure_external_services_selection_from_catalog();
        self.external_services_ui.request_json = self.default_external_service_request_json();
        self.external_services_ui.preflight_output = None;
        self.external_services_ui.quote_output = None;
        self.external_services_ui.selected_quote_payload = 0;
        self.external_services_ui.quote_output_dir =
            self.default_external_service_quote_output_dir();
        self.external_services_ui.status = format!(
            "Prepared editable request template for {} / {}",
            self.external_services_ui.selected_provider,
            self.external_services_ui.selected_service_kind
        );
    }

    pub(super) fn run_external_services_preflight(&mut self) {
        let command = ShellCommand::ServicesProjectPreflight {
            request_json: self.external_services_ui.request_json.clone(),
        };
        match self.execute_shared_shell_command_json(&command) {
            Ok((output, _)) => {
                let eligible = output
                    .get("eligible")
                    .and_then(serde_json::Value::as_bool)
                    .unwrap_or(false);
                let blocking_count = output
                    .get("blocking_issues")
                    .and_then(serde_json::Value::as_array)
                    .map(Vec::len)
                    .unwrap_or(0);
                self.external_services_ui.preflight_output = Some(output);
                self.external_services_ui.status = format!(
                    "Preflight complete: eligible={eligible}, blocking issues={blocking_count}"
                );
            }
            Err(err) => {
                self.external_services_ui.status =
                    format!("External-service preflight failed: {err}");
            }
        }
    }

    fn run_external_services_quote(&mut self) {
        let command = ShellCommand::ServicesProjectQuote {
            request_json: self.external_services_ui.request_json.clone(),
            output_dir: None,
        };
        match self.execute_shared_shell_command_json(&command) {
            Ok((output, _)) => {
                let quote_status = output
                    .get("quote_status")
                    .and_then(serde_json::Value::as_str)
                    .unwrap_or("unknown")
                    .to_string();
                let payload_count = output
                    .pointer("/service_ready_bundle/inline_payloads")
                    .and_then(serde_json::Value::as_array)
                    .map(Vec::len)
                    .unwrap_or(0);
                self.external_services_ui.quote_output = Some(output);
                self.external_services_ui.selected_quote_payload = 0;
                self.external_services_ui.status = format!(
                    "Quote handoff complete: status={quote_status}, inline payloads={payload_count}"
                );
            }
            Err(err) => {
                self.external_services_ui.status =
                    format!("External-service quote handoff failed: {err}");
            }
        }
    }

    pub(super) fn export_external_services_quote_bundle(&mut self) {
        if self.external_services_ui.quote_output_dir.trim().is_empty() {
            self.external_services_ui.quote_output_dir =
                self.default_external_service_quote_output_dir();
        }
        let output_dir = self
            .external_services_ui
            .quote_output_dir
            .trim()
            .to_string();
        let command = ShellCommand::ServicesProjectQuote {
            request_json: self.external_services_ui.request_json.clone(),
            output_dir: Some(output_dir.clone()),
        };
        match self.execute_shared_shell_command_json(&command) {
            Ok((output, _)) => {
                let quote_status = output
                    .get("quote_status")
                    .and_then(serde_json::Value::as_str)
                    .unwrap_or("unknown")
                    .to_string();
                let file_count = output
                    .pointer("/service_ready_bundle/local_files")
                    .and_then(serde_json::Value::as_array)
                    .map(Vec::len)
                    .unwrap_or(0);
                self.external_services_ui.quote_output = Some(output);
                self.external_services_ui.selected_quote_payload = 0;
                self.external_services_ui.status = format!(
                    "Exported quote handoff bundle to {output_dir}: status={quote_status}, files={file_count}"
                );
            }
            Err(err) => {
                self.external_services_ui.status =
                    format!("External-service quote bundle export failed: {err}");
            }
        }
    }

    fn render_external_services_provider_picker(&mut self, ui: &mut Ui) {
        let providers = self.external_service_provider_choices();
        if providers.is_empty() {
            ui.small("No provider catalog loaded yet.");
            return;
        }
        let mut provider_changed = false;
        ui.horizontal_wrapped(|ui| {
            ui.label("Provider");
            egui::ComboBox::from_id_salt("external_services_provider_picker")
                .selected_text(
                    providers
                        .iter()
                        .find(|row| row.provider == self.external_services_ui.selected_provider)
                        .map(|row| row.display_name.as_str())
                        .unwrap_or(self.external_services_ui.selected_provider.as_str()),
                )
                .show_ui(ui, |ui| {
                    for provider in &providers {
                        let label = format!("{} ({})", provider.display_name, provider.provider);
                        if ui
                            .selectable_value(
                                &mut self.external_services_ui.selected_provider,
                                provider.provider.clone(),
                                label,
                            )
                            .changed()
                        {
                            provider_changed = true;
                        }
                    }
                });

            if let Some(provider) = providers
                .iter()
                .find(|row| row.provider == self.external_services_ui.selected_provider)
            {
                ui.small(format!("status: {}", provider.support_status));
            }
        });
        if provider_changed {
            self.reset_external_services_request_from_selection();
        }

        let capabilities = providers
            .iter()
            .find(|row| row.provider == self.external_services_ui.selected_provider)
            .map(|row| row.capabilities.clone())
            .unwrap_or_default();
        if capabilities.is_empty() {
            ui.small("Selected provider has no configured capabilities.");
            return;
        }
        let mut service_kind_changed = false;
        ui.horizontal_wrapped(|ui| {
            ui.label("Service kind");
            egui::ComboBox::from_id_salt("external_services_service_kind_picker")
                .selected_text(
                    capabilities
                        .iter()
                        .find(|row| {
                            row.service_kind == self.external_services_ui.selected_service_kind
                        })
                        .map(|row| row.display_name.as_str())
                        .unwrap_or(self.external_services_ui.selected_service_kind.as_str()),
                )
                .show_ui(ui, |ui| {
                    for capability in &capabilities {
                        let label =
                            format!("{} ({})", capability.display_name, capability.service_kind);
                        if ui
                            .selectable_value(
                                &mut self.external_services_ui.selected_service_kind,
                                capability.service_kind.clone(),
                                label,
                            )
                            .changed()
                        {
                            service_kind_changed = true;
                        }
                    }
                });
        });
        if service_kind_changed {
            self.reset_external_services_request_from_selection();
        }
    }

    fn render_json_preview(ui: &mut Ui, title: &str, value: &Option<serde_json::Value>) {
        let Some(value) = value else {
            ui.small(format!("{title}: not run yet"));
            return;
        };
        egui::CollapsingHeader::new(title)
            .default_open(false)
            .show(ui, |ui| {
                let text = serde_json::to_string_pretty(value)
                    .unwrap_or_else(|err| format!("Could not format JSON: {err}"));
                egui::ScrollArea::vertical()
                    .id_salt(format!("external_services_json_preview_{title}"))
                    .max_height(220.0)
                    .show(ui, |ui| {
                        ui.monospace(text);
                    });
            });
    }

    fn render_quote_payloads(&mut self, ui: &mut Ui) {
        let Some(quote) = self.external_services_ui.quote_output.as_ref() else {
            ui.small("Quote handoff: not run yet");
            return;
        };
        let payloads = quote
            .pointer("/service_ready_bundle/inline_payloads")
            .and_then(serde_json::Value::as_array)
            .cloned()
            .unwrap_or_default();
        if payloads.is_empty() {
            ui.small("Quote handoff produced no inline payloads.");
            return;
        }
        if self.external_services_ui.selected_quote_payload >= payloads.len() {
            self.external_services_ui.selected_quote_payload = payloads.len() - 1;
        }
        ui.horizontal_wrapped(|ui| {
            ui.label("Payload");
            egui::ComboBox::from_id_salt("external_services_quote_payload_picker")
                .selected_text(
                    payloads
                        .get(self.external_services_ui.selected_quote_payload)
                        .and_then(|payload| payload.get("payload_kind"))
                        .and_then(serde_json::Value::as_str)
                        .unwrap_or("payload"),
                )
                .show_ui(ui, |ui| {
                    for (idx, payload) in payloads.iter().enumerate() {
                        let kind = payload
                            .get("payload_kind")
                            .and_then(serde_json::Value::as_str)
                            .unwrap_or("payload");
                        ui.selectable_value(
                            &mut self.external_services_ui.selected_quote_payload,
                            idx,
                            kind,
                        );
                    }
                });
        });
        let payload = &payloads[self.external_services_ui.selected_quote_payload];
        if let Some(description) = payload
            .get("description")
            .and_then(serde_json::Value::as_str)
            .filter(|value| !value.trim().is_empty())
        {
            ui.small(description);
        }
        let text = payload
            .get("text")
            .and_then(serde_json::Value::as_str)
            .unwrap_or("");
        egui::ScrollArea::vertical()
            .id_salt("external_services_quote_payload_text")
            .max_height(260.0)
            .show(ui, |ui| {
                ui.monospace(text);
            });
        let local_files = quote
            .pointer("/service_ready_bundle/local_files")
            .and_then(serde_json::Value::as_array)
            .cloned()
            .unwrap_or_default();
        if !local_files.is_empty() {
            ui.separator();
            ui.strong("Bundle files");
            for file in &local_files {
                let kind = file
                    .get("artifact_kind")
                    .and_then(serde_json::Value::as_str)
                    .unwrap_or("artifact");
                let path = file
                    .get("path")
                    .and_then(serde_json::Value::as_str)
                    .unwrap_or("-");
                ui.monospace(format!("{kind}: {path}"));
            }
        }
    }

    fn render_external_services_contents(&mut self, ui: &mut Ui) {
        let close_requested = self.render_specialist_window_nav_with_close(
            ui,
            Some(("Close", "Close the External Services window (Cmd/Ctrl+W)")),
        );
        if close_requested {
            self.external_services_ui.show_panel = false;
        }

        ui.heading("External Services");
        ui.label(
            "Provider-neutral quote/handoff preparation. This window calls the same shared-shell services routes used by CLI, agents, MCP, and ClawBio.",
        );
        ui.small(
            "No WOP scraping, API submission, credential lookup, PO, shipping, or billing persistence is performed here.",
        );
        ui.separator();

        ui.horizontal_wrapped(|ui| {
            if ui
                .button("Refresh Providers")
                .on_hover_text("Run `services providers list`")
                .clicked()
            {
                self.refresh_external_services_provider_catalog();
            }
            if ui
                .button("Provider Config Doctor")
                .on_hover_text("Run `services providers doctor` against the overlay discovery chain")
                .clicked()
            {
                self.refresh_external_services_provider_doctor();
            }
            if ui
                .button("Use Selected Template")
                .on_hover_text("Replace the editable request JSON with a provider-neutral template for the selected provider/service kind")
                .clicked()
            {
                self.reset_external_services_request_from_selection();
            }
        });
        if !self.external_services_ui.status.trim().is_empty() {
            ui.monospace(self.external_services_ui.status.clone());
        }

        ui.separator();
        self.render_external_services_provider_picker(ui);

        ui.separator();
        ui.strong("Editable service request JSON");
        ui.small("Keep this provider-neutral. Provider-specific mapping is resolved by the shared services layer.");
        egui::ScrollArea::vertical()
            .id_salt("external_services_request_json_scroll")
            .max_height(260.0)
            .show(ui, |ui| {
                ui.add(
                    egui::TextEdit::multiline(&mut self.external_services_ui.request_json)
                        .desired_width(f32::INFINITY)
                        .desired_rows(12)
                        .code_editor(),
                );
            });
        ui.horizontal_wrapped(|ui| {
            if ui
                .button("Preflight")
                .on_hover_text("Run `services project-preflight` for the editable request")
                .clicked()
            {
                self.run_external_services_preflight();
            }
            if ui
                .button("Prepare Quote Handoff")
                .on_hover_text(
                    "Run `services project-quote` for the editable request; no order is submitted",
                )
                .clicked()
            {
                self.run_external_services_quote();
            }
        });
        ui.horizontal_wrapped(|ui| {
            ui.label("Output dir");
            ui.add(
                egui::TextEdit::singleline(&mut self.external_services_ui.quote_output_dir)
                    .desired_width(420.0),
            )
            .on_hover_text(
                "Directory for generated quote handoff files; defaults under ignored artifacts/",
            );
            if ui
                .button("Export Handoff Bundle")
                .on_hover_text("Run `services project-quote --output-dir DIR` and write bundle files locally; no order is submitted")
                .clicked()
            {
                self.export_external_services_quote_bundle();
            }
        });

        ui.separator();
        ui.strong("Handoff payload preview");
        self.render_quote_payloads(ui);

        ui.separator();
        Self::render_json_preview(
            ui,
            "Provider Config Doctor JSON",
            &self.external_services_ui.provider_doctor_output,
        );
        Self::render_json_preview(
            ui,
            "Preflight JSON",
            &self.external_services_ui.preflight_output,
        );
        Self::render_json_preview(ui, "Quote JSON", &self.external_services_ui.quote_output);
        Self::render_json_preview(
            ui,
            "Provider Catalog JSON",
            &self.external_services_ui.provider_catalog_output,
        );
    }

    pub(super) fn render_external_services_dialog(&mut self, ctx: &egui::Context) {
        if !self.external_services_ui.show_panel {
            return;
        }
        let mut open = self.external_services_ui.show_panel;
        let viewport_id = Self::external_services_viewport_id();
        let spec = self.hosted_window_spec_for_viewport(
            "External Services",
            egui::Id::new(("external_services_hosted_window", viewport_id)),
            viewport_id,
            Vec2::new(920.0, 720.0),
            Vec2::new(620.0, 420.0),
        );
        crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
            self.render_external_services_contents(ui);
        });
        self.clear_viewport_foreground_request_after_render(viewport_id);
        self.finalize_viewport_open_probe(viewport_id, "External Services");
        self.external_services_ui.show_panel = open && self.external_services_ui.show_panel;
    }
}
