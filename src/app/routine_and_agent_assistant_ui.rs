//! Routine Assistant and Agent Assistant GUI helpers.
//!
//! This module is a move-only extraction from `app.rs`: it keeps the
//! intertwined Routine Assistant and Agent Assistant dialog/rendering helpers
//! close to `GENtleApp` while reducing the top-level app monolith.

use super::*;
use crate::agent_bridge::AgentSuggestedCommand;

impl GENtleApp {
    const AGENT_MODEL_SELECTION_REQUIRED_MESSAGE: &'static str =
        "Connection established, please select the model to use in the drop-down box.";
    const AGENT_SCREENSHOT_CAPTURE_TIMEOUT: Duration = Duration::from_secs(15);
    const AGENT_DISCARDED_SCREENSHOT_CAPTURE_LIMIT: usize = 64;

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

    pub(super) fn agent_initial_actions_visible(&self) -> bool {
        self.agent_conversation.turns.is_empty() && self.agent_last_invocation.is_none()
    }

    fn encode_agent_help_png(image: &egui::ColorImage) -> Result<Vec<u8>, String> {
        let [width, height] = image.size;
        if width == 0 || height == 0 || image.pixels.len() != width.saturating_mul(height) {
            return Err("Captured GENtle view has invalid pixel dimensions".to_string());
        }
        let rgba = image
            .pixels
            .iter()
            .flat_map(|pixel| pixel.to_array())
            .collect::<Vec<_>>();
        let buffer = image::RgbaImage::from_raw(width as u32, height as u32, rgba)
            .ok_or_else(|| "Could not construct captured GENtle image".to_string())?;
        let mut encoded = Cursor::new(Vec::new());
        image::DynamicImage::ImageRgba8(buffer)
            .write_to(&mut encoded, image::ImageFormat::Png)
            .map_err(|error| format!("Could not encode captured GENtle view as PNG: {error}"))?;
        Ok(encoded.into_inner())
    }

    fn agent_help_sha256(bytes: &[u8]) -> String {
        ring::digest::digest(&ring::digest::SHA256, bytes)
            .as_ref()
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect()
    }

    fn prepare_agent_help_attachment(
        capture: crate::agent_help::AgentHelpCapturedImage,
    ) -> Result<AgentPendingImageAttachment, String> {
        let png_bytes = Self::encode_agent_help_png(&capture.image)?;
        if png_bytes.len() > 20 * 1024 * 1024 {
            return Err(format!(
                "Captured image is too large for Agent Assistant ({} MiB; limit 20 MiB)",
                png_bytes.len() / (1024 * 1024)
            ));
        }
        let mut temp_file = tempfile::Builder::new()
            .prefix("gentle-agent-help-")
            .suffix(".png")
            .tempfile()
            .map_err(|error| format!("Could not create temporary screenshot file: {error}"))?;
        temp_file
            .write_all(&png_bytes)
            .map_err(|error| format!("Could not write temporary screenshot file: {error}"))?;
        temp_file
            .flush()
            .map_err(|error| format!("Could not flush temporary screenshot file: {error}"))?;
        let path = temp_file
            .path()
            .canonicalize()
            .map_err(|error| format!("Could not resolve temporary screenshot path: {error}"))?;
        let [pixel_width, pixel_height] = capture.image.size;
        let file_name = path
            .file_name()
            .and_then(|value| value.to_str())
            .unwrap_or("gentle-agent-help.png")
            .to_string();
        let request = AgentRequestAttachment {
            schema: AGENT_ATTACHMENT_SCHEMA.to_string(),
            id: format!("agent_help_{}", capture.request_id),
            kind: "image".to_string(),
            file_name,
            mime_type: "image/png".to_string(),
            path: path.to_string_lossy().to_string(),
            byte_len: png_bytes.len() as u64,
            sha256: Self::agent_help_sha256(&png_bytes),
            source_window_title: Some(capture.window_title),
            capture_backend: Some(capture.backend),
            pixel_width: Some(pixel_width),
            pixel_height: Some(pixel_height),
        };
        Ok(AgentPendingImageAttachment {
            request,
            png_bytes: Arc::<[u8]>::from(png_bytes),
            temp_file: Arc::new(temp_file),
        })
    }

    fn agent_help_prompt(window_title: &str) -> String {
        crate::i18n::tr("agent.help_prompt").replace("{window}", window_title)
    }

    fn agent_screenshot_project_generation(&self) -> AgentScreenshotProjectGeneration {
        AgentScreenshotProjectGeneration {
            engine_identity: Arc::as_ptr(&self.engine) as usize,
            structural_revision: self
                .engine
                .read()
                .map(|engine| engine.structural_revision())
                .unwrap_or(u64::MAX),
        }
    }

    fn agent_screenshot_targets(&self) -> Vec<OpenWindowEntry> {
        let entries = self
            .collect_open_window_entries()
            .into_iter()
            .filter(|entry| entry.viewport_id != Self::agent_assistant_viewport_id())
            .collect::<Vec<_>>();
        let mut key_counts = HashMap::<u64, usize>::new();
        for entry in &entries {
            *key_counts.entry(entry.native_menu_key).or_default() += 1;
        }
        entries
            .into_iter()
            .filter(|entry| key_counts.get(&entry.native_menu_key) == Some(&1))
            .collect()
    }

    fn resolve_agent_screenshot_target(
        &self,
        native_menu_key: u64,
    ) -> Result<OpenWindowEntry, String> {
        let mut matches = self
            .agent_screenshot_targets()
            .into_iter()
            .filter(|entry| entry.native_menu_key == native_menu_key);
        let Some(entry) = matches.next() else {
            return Err(
                "The selected GENtle window is closed, unavailable, or ambiguous. Ask the agent again if a new screenshot is still needed."
                    .to_string(),
            );
        };
        if matches.next().is_some() {
            return Err(
                "The selected GENtle window is ambiguous, so no screenshot was captured."
                    .to_string(),
            );
        }
        Ok(entry)
    }

    fn agent_screenshot_turn_is_present(
        &self,
        system_id: &str,
        response_completed_at_unix_ms: u128,
        request_id: &str,
    ) -> bool {
        self.agent_conversation.turns.iter().any(|turn| {
            turn.system_id == system_id
                && turn.completed_at_unix_ms == response_completed_at_unix_ms
                && turn
                    .response
                    .screenshot_request
                    .as_ref()
                    .is_some_and(|request| request.id == request_id)
        })
    }

    fn remember_discarded_agent_screenshot_capture(&mut self, request_id: u64) {
        if !self
            .agent_discarded_screenshot_capture_ids
            .contains(&request_id)
        {
            self.agent_discarded_screenshot_capture_ids
                .push_back(request_id);
        }
        while self.agent_discarded_screenshot_capture_ids.len()
            > Self::AGENT_DISCARDED_SCREENSHOT_CAPTURE_LIMIT
        {
            self.agent_discarded_screenshot_capture_ids.pop_front();
        }
    }

    fn consume_discarded_agent_screenshot_capture(&mut self, request_id: u64) -> bool {
        let Some(position) = self
            .agent_discarded_screenshot_capture_ids
            .iter()
            .position(|candidate| *candidate == request_id)
        else {
            return false;
        };
        self.agent_discarded_screenshot_capture_ids.remove(position);
        true
    }

    fn invalidate_agent_screenshot_state(&mut self, reason: Option<&str>) {
        let had_state = self.agent_screenshot_consent.take().is_some();
        if let Some(capture) = self.agent_screenshot_capture.take() {
            self.remember_discarded_agent_screenshot_capture(capture.capture_request_id);
            if let Some(reason) = reason {
                self.agent_help_capture_failure = Some(AgentHelpCaptureFailure {
                    request_id: capture.capture_request_id,
                    window_title: capture.source_window_title,
                    kind: crate::agent_help::AgentHelpCaptureFailureKind::CaptureFailed,
                    message: reason.to_string(),
                });
            }
        }
        if (had_state || self.agent_help_capture_failure.is_some())
            && let Some(reason) = reason
        {
            self.agent_status = reason.to_string();
        }
    }

    pub(super) fn activate_agent_screenshot_consent(
        &mut self,
        request: AgentScreenshotRequest,
        system_id: String,
        system_label: String,
        response_completed_at_unix_ms: u128,
    ) {
        self.invalidate_agent_screenshot_state(None);
        let targets = self.agent_screenshot_targets();
        let selected_window_key = self
            .active_window_menu_key
            .filter(|key| targets.iter().any(|entry| entry.native_menu_key == *key))
            .or_else(|| {
                targets
                    .iter()
                    .find(|entry| entry.viewport_id == ViewportId::ROOT)
                    .map(|entry| entry.native_menu_key)
            })
            .or_else(|| targets.first().map(|entry| entry.native_menu_key));
        self.agent_screenshot_consent = Some(AgentScreenshotConsent {
            request,
            system_id,
            system_label,
            response_completed_at_unix_ms,
            project_generation: self.agent_screenshot_project_generation(),
            selected_window_key,
        });
    }

    fn agent_screenshot_consent_binding_error(
        &self,
        consent: &AgentScreenshotConsent,
    ) -> Option<String> {
        if self.agent_task.is_some() {
            return Some(
                "A new agent request is running, so the earlier screenshot consent request expired."
                    .to_string(),
            );
        }
        if self.agent_system_id.trim() != consent.system_id {
            return Some(
                "The selected agent system changed, so the earlier screenshot consent request expired."
                    .to_string(),
            );
        }
        if self.agent_screenshot_project_generation() != consent.project_generation {
            return Some(
                "The project changed, so the earlier screenshot consent request expired."
                    .to_string(),
            );
        }
        if !self.agent_screenshot_turn_is_present(
            &consent.system_id,
            consent.response_completed_at_unix_ms,
            &consent.request.id,
        ) {
            return Some(
                "The originating agent response is no longer active, so its screenshot request expired."
                    .to_string(),
            );
        }
        if let Some(key) = consent.selected_window_key
            && let Err(error) = self.resolve_agent_screenshot_target(key)
        {
            return Some(error);
        }
        None
    }

    fn agent_screenshot_capture_binding_error(
        &self,
        capture: &AgentScreenshotCapture,
    ) -> Option<String> {
        if self.agent_task.is_some() {
            return Some(
                "An agent request started before capture completed, so the screenshot was discarded."
                    .to_string(),
            );
        }
        if self.agent_system_id.trim() != capture.system_id {
            return Some(
                "The selected agent system changed before capture completed, so the screenshot was discarded."
                    .to_string(),
            );
        }
        if self.agent_screenshot_project_generation() != capture.project_generation {
            return Some(
                "The project changed before capture completed, so the screenshot was discarded."
                    .to_string(),
            );
        }
        if !self.agent_screenshot_turn_is_present(
            &capture.system_id,
            capture.response_completed_at_unix_ms,
            &capture.agent_request_id,
        ) {
            return Some(
                "The originating agent response changed before capture completed, so the screenshot was discarded."
                    .to_string(),
            );
        }
        match self.resolve_agent_screenshot_target(capture.target_window_key) {
            Ok(entry) if entry.viewport_id == capture.target_viewport_id => {}
            Ok(_) => {
                return Some(
                    "The selected GENtle window changed identity before capture completed, so the screenshot was discarded."
                        .to_string(),
                );
            }
            Err(error) => return Some(error),
        }
        if capture.started.elapsed() > Self::AGENT_SCREENSHOT_CAPTURE_TIMEOUT {
            return Some(
                "The approved screenshot capture timed out. No image was attached or sent."
                    .to_string(),
            );
        }
        None
    }

    pub(super) fn validate_agent_screenshot_state(&mut self) {
        let invalid_reason = self
            .agent_screenshot_consent
            .as_ref()
            .and_then(|consent| self.agent_screenshot_consent_binding_error(consent))
            .or_else(|| {
                self.agent_screenshot_capture
                    .as_ref()
                    .and_then(|capture| self.agent_screenshot_capture_binding_error(capture))
            });
        if let Some(reason) = invalid_reason {
            self.invalidate_agent_screenshot_state(Some(&reason));
        }
    }

    pub(super) fn decline_agent_screenshot_request(&mut self) {
        if self.agent_screenshot_consent.take().is_some() {
            self.agent_status =
                "Screenshot request declined. No image was captured or attached.".to_string();
        }
    }

    pub(super) fn approve_agent_screenshot_request(&mut self, ctx: &egui::Context) {
        let Some(consent) = self.agent_screenshot_consent.take() else {
            self.agent_status =
                "This screenshot request is no longer active. No image was captured.".to_string();
            return;
        };
        if let Some(error) = self.agent_screenshot_consent_binding_error(&consent) {
            self.agent_status = error;
            return;
        }
        let Some(system) = self
            .agent_systems
            .iter()
            .find(|system| system.id == consent.system_id)
        else {
            self.agent_status =
                "The originating agent system is unavailable. No image was captured.".to_string();
            return;
        };
        if !system.supports_image_attachments {
            self.agent_status = format!(
                "Agent system '{}' cannot receive image attachments. No image was captured.",
                consent.system_label
            );
            return;
        }
        let Some(target_key) = consent.selected_window_key else {
            self.agent_status =
                "Select one registered GENtle content window before allowing a screenshot."
                    .to_string();
            self.agent_screenshot_consent = Some(consent);
            return;
        };
        let target = match self.resolve_agent_screenshot_target(target_key) {
            Ok(target) => target,
            Err(error) => {
                self.agent_status = error;
                return;
            }
        };
        let capture_viewport_id = if ctx.embed_viewports()
            && self
                .embedded_window_layer_id_for_viewport(target.viewport_id)
                .is_some()
        {
            ViewportId::ROOT
        } else {
            target.viewport_id
        };
        let capture_request_id = crate::agent_help::request_egui_viewport_capture_for(
            ctx,
            capture_viewport_id,
            target.title.clone(),
        );
        self.agent_screenshot_capture = Some(AgentScreenshotCapture {
            capture_request_id,
            agent_request_id: consent.request.id.clone(),
            system_id: consent.system_id,
            system_label: consent.system_label,
            response_completed_at_unix_ms: consent.response_completed_at_unix_ms,
            project_generation: consent.project_generation,
            target_window_key: target.native_menu_key,
            target_viewport_id: target.viewport_id,
            capture_viewport_id,
            source_window_title: target.title.clone(),
            started: Instant::now(),
        });
        self.agent_status = format!(
            "Capturing one user-approved screenshot from '{}'. It will be previewed locally before it can be sent.",
            target.title
        );
        ctx.request_repaint_after(Duration::from_millis(100));
    }

    fn agent_screenshot_followup_prompt(request_id: &str, window_title: &str) -> String {
        let bounded_title = window_title.chars().take(200).collect::<String>();
        format!(
            "This is the user-approved screenshot for agent screenshot request '{request_id}' from the registered GENtle window '{bounded_title}'. Inspect only the attached image, distinguish visible evidence from inference, and answer the original visual question."
        )
    }

    fn render_agent_screenshot_consent_card(&mut self, ui: &mut egui::Ui) {
        self.validate_agent_screenshot_state();
        if let Some(capture) = self.agent_screenshot_capture.clone() {
            ui.group(|ui| {
                ui.strong(self.tr("agent.screenshot_request.title"));
                ui.add(egui::Spinner::new());
                ui.label(
                    self.tr("agent.screenshot_request.capturing")
                        .replace("{window}", &capture.source_window_title),
                );
                ui.small(self.tr("agent.screenshot_request.preview_notice"));
            });
            return;
        }

        let Some(consent) = self.agent_screenshot_consent.clone() else {
            return;
        };
        let targets = self.agent_screenshot_targets();
        if self
            .agent_screenshot_consent
            .as_ref()
            .and_then(|state| state.selected_window_key)
            .is_some_and(|selected| {
                !targets
                    .iter()
                    .any(|entry| entry.native_menu_key == selected)
            })
        {
            self.invalidate_agent_screenshot_state(Some(
                "The selected GENtle window closed, so the screenshot request expired.",
            ));
            return;
        }

        let supports_images = self
            .agent_systems
            .iter()
            .find(|system| system.id == consent.system_id)
            .is_some_and(|system| system.supports_image_attachments);
        let selected_window_key = self
            .agent_screenshot_consent
            .as_ref()
            .and_then(|state| state.selected_window_key);
        let selected_text = selected_window_key
            .and_then(|selected| {
                targets
                    .iter()
                    .find(|entry| entry.native_menu_key == selected)
            })
            .map(|entry| entry.title.clone())
            .unwrap_or_else(|| self.tr("agent.screenshot_request.choose_window"));
        let mut selected_after = selected_window_key;
        let mut allow_clicked = false;
        let mut decline_clicked = false;
        let can_allow = supports_images
            && self.agent_task.is_none()
            && selected_window_key.is_some()
            && self
                .agent_screenshot_consent_binding_error(&consent)
                .is_none();

        ui.group(|ui| {
            ui.strong(self.tr("agent.screenshot_request.title"));
            ui.horizontal_wrapped(|ui| {
                ui.small(self.tr("agent.screenshot_request.system"));
                ui.monospace(format!("{} ({})", consent.system_label, consent.system_id));
            });
            ui.label(
                self.tr("agent.screenshot_request.reason")
                    .replace("{reason}", &consent.request.reason),
            );
            ui.horizontal_wrapped(|ui| {
                ui.label(self.tr("agent.screenshot_request.window"));
                egui::ComboBox::from_id_salt((
                    "agent_screenshot_target",
                    consent.response_completed_at_unix_ms,
                    consent.request.id.as_str(),
                ))
                .selected_text(selected_text)
                .show_ui(ui, |ui| {
                    for target in &targets {
                        ui.selectable_value(
                            &mut selected_after,
                            Some(target.native_menu_key),
                            format!("{} - {}", target.title, target.detail),
                        );
                    }
                });
            });
            if targets.is_empty() {
                ui.colored_label(
                    egui::Color32::from_rgb(180, 70, 45),
                    self.tr("agent.screenshot_request.no_windows"),
                );
            }
            if !supports_images {
                ui.colored_label(
                    egui::Color32::from_rgb(180, 70, 45),
                    self.tr("agent.screenshot_request.unsupported")
                        .replace("{system}", &consent.system_label),
                );
            }
            ui.small(self.tr("agent.screenshot_request.preview_notice"));
            ui.horizontal(|ui| {
                let allow_response = ui.add_enabled(
                    can_allow,
                    egui::Button::new(self.tr("agent.screenshot_request.allow_once")),
                );
                let decline_response = ui.button(self.tr("agent.screenshot_request.decline"));
                #[cfg(feature = "gui-test-support")]
                {
                    crate::gui_test_support::register_response(
                        &allow_response,
                        "agent.screenshot.allow_once",
                        "window.agent_assistant",
                        None,
                        crate::gui_test_support::GuiTestWidgetKind::Button,
                        false,
                    );
                    crate::gui_test_support::register_response(
                        &decline_response,
                        "agent.screenshot.decline",
                        "window.agent_assistant",
                        None,
                        crate::gui_test_support::GuiTestWidgetKind::Button,
                        false,
                    );
                }
                allow_clicked = allow_response.clicked();
                decline_clicked = decline_response.clicked();
            });
        });

        if let Some(state) = self.agent_screenshot_consent.as_mut() {
            state.selected_window_key = selected_after;
        }
        if decline_clicked {
            self.decline_agent_screenshot_request();
        } else if allow_clicked {
            self.approve_agent_screenshot_request(ui.ctx());
        }
    }

    fn render_agent_help_attachment_panel(&mut self, ui: &mut egui::Ui) {
        if let Some(attachment) = self.agent_pending_image_attachment.clone() {
            let supports_images = self
                .selected_agent_system()
                .map(|system| system.supports_image_attachments)
                .unwrap_or(false);
            ui.group(|ui| {
                ui.horizontal_wrapped(|ui| {
                    ui.strong("Attached screenshot");
                    ui.small(
                        attachment
                            .request
                            .source_window_title
                            .as_deref()
                            .unwrap_or("GENtle window"),
                    );
                    if let (Some(width), Some(height)) = (
                        attachment.request.pixel_width,
                        attachment.request.pixel_height,
                    ) {
                        ui.small(format!("{width} x {height} px"));
                    }
                    ui.small(format!(
                        "{} KiB | {}",
                        attachment.request.byte_len.div_ceil(1024),
                        attachment
                            .request
                            .capture_backend
                            .as_deref()
                            .unwrap_or("capture")
                    ));
                });
                ui.add(
                    egui::Image::from_bytes(
                        format!(
                            "bytes://gentle-agent-help-{}.png",
                            attachment.request.sha256
                        ),
                        attachment.png_bytes.clone(),
                    )
                    .max_width(ui.available_width())
                    .max_height(260.0)
                    .shrink_to_fit(),
                );
                ui.small(
                    "The screenshot remains local until you click Ask agent. Only this selected image is attached.",
                );
                if !supports_images {
                    ui.colored_label(
                        egui::Color32::from_rgb(180, 70, 45),
                        "The selected agent system does not support image attachments. Choose an image-capable system or remove the screenshot.",
                    );
                }
                let remove_response = ui
                    .add_enabled(
                        self.agent_task.is_none(),
                        egui::Button::new("Remove screenshot"),
                    );
                #[cfg(feature = "gui-test-support")]
                crate::gui_test_support::register_response(
                    &remove_response,
                    "agent.screenshot.remove",
                    "window.agent_assistant",
                    None,
                    crate::gui_test_support::GuiTestWidgetKind::Button,
                    false,
                );
                if remove_response.clicked() {
                    self.agent_pending_image_attachment = None;
                    self.agent_status = "Pending screenshot removed".to_string();
                }
            });
        }

        if let Some(failure) = self.agent_help_capture_failure.clone() {
            ui.group(|ui| {
                ui.colored_label(
                    egui::Color32::from_rgb(180, 70, 45),
                    format!("Screenshot unavailable: {}", failure.message),
                );
                ui.small(format!("Requested for '{}'.", failure.window_title));
                if matches!(
                    failure.kind,
                    crate::agent_help::AgentHelpCaptureFailureKind::PermissionRequired
                        | crate::agent_help::AgentHelpCaptureFailureKind::RestartRequired
                ) && ui.button("Open Screen Recording settings").clicked()
                {
                    self.agent_status = match crate::agent_help::open_macos_screen_recording_settings()
                    {
                        Ok(()) => {
                            "Opened macOS Screen Recording settings. Restart GENtle after granting access."
                                .to_string()
                        }
                        Err(error) => error,
                    };
                }
                ui.small(
                    "The normal Agent help click captures only GENtle's drawn viewport and does not require Screen Recording permission; native full-window capture is optional.",
                );
            });
        }
    }

    pub(super) fn poll_agent_help_capture_events(&mut self, ctx: &egui::Context) {
        self.validate_agent_screenshot_state();
        if let Some(capture) = &self.agent_screenshot_capture {
            crate::agent_help::collect_egui_capture_events_for(ctx, capture.capture_viewport_id);
            ctx.request_repaint_after(Duration::from_millis(100));
        }
        self.validate_agent_screenshot_state();
        for event in take_capture_events() {
            match event {
                AgentHelpCaptureEvent::Captured(capture) => {
                    let request_id = capture.request_id;
                    let window_title = capture.window_title.clone();
                    if self.consume_discarded_agent_screenshot_capture(request_id) {
                        continue;
                    }
                    let agent_requested = self
                        .agent_screenshot_capture
                        .as_ref()
                        .is_some_and(|pending| pending.capture_request_id == request_id);
                    if agent_requested {
                        let pending = self
                            .agent_screenshot_capture
                            .take()
                            .expect("matching screenshot capture exists");
                        if let Some(message) = self.agent_screenshot_capture_binding_error(&pending)
                        {
                            self.agent_pending_image_attachment = None;
                            self.agent_help_capture_failure = Some(AgentHelpCaptureFailure {
                                request_id,
                                window_title,
                                kind: crate::agent_help::AgentHelpCaptureFailureKind::CaptureFailed,
                                message: message.clone(),
                            });
                            self.agent_status = message;
                            self.open_agent_assistant_dialog();
                            continue;
                        }
                        match Self::prepare_agent_help_attachment(capture) {
                            Ok(attachment) => {
                                self.agent_pending_image_attachment = Some(attachment);
                                self.agent_help_capture_failure = None;
                                self.agent_prompt = Self::agent_screenshot_followup_prompt(
                                    &pending.agent_request_id,
                                    &pending.source_window_title,
                                );
                                self.agent_status = format!(
                                    "One screenshot from '{}' requested by '{}' is attached locally. Review the preview and prompt, then click Ask Agent to send it.",
                                    pending.source_window_title, pending.system_label
                                );
                            }
                            Err(message) => {
                                self.agent_pending_image_attachment = None;
                                self.agent_help_capture_failure = Some(AgentHelpCaptureFailure {
                                    request_id,
                                    window_title: pending.source_window_title,
                                    kind: crate::agent_help::AgentHelpCaptureFailureKind::CaptureFailed,
                                    message: message.clone(),
                                });
                                self.agent_status = message;
                            }
                        }
                        self.open_agent_assistant_dialog();
                        continue;
                    }
                    self.invalidate_agent_screenshot_state(Some(
                        "A user-invoked Agent Help capture superseded the pending agent screenshot request.",
                    ));
                    match Self::prepare_agent_help_attachment(capture) {
                        Ok(attachment) => {
                            self.agent_pending_image_attachment = Some(attachment);
                            self.agent_help_capture_failure = None;
                            self.agent_prompt = Self::agent_help_prompt(&window_title);
                            self.agent_status = format!(
                                "Screenshot from '{window_title}' is attached locally. Review it and your prompt before asking the agent."
                            );
                        }
                        Err(message) => {
                            self.agent_pending_image_attachment = None;
                            self.agent_help_capture_failure = Some(AgentHelpCaptureFailure {
                                request_id,
                                window_title: window_title.clone(),
                                kind: crate::agent_help::AgentHelpCaptureFailureKind::CaptureFailed,
                                message: message.clone(),
                            });
                            self.agent_status = message;
                        }
                    }
                    self.open_agent_assistant_dialog();
                }
                AgentHelpCaptureEvent::Failed(failure) => {
                    if self.consume_discarded_agent_screenshot_capture(failure.request_id) {
                        continue;
                    }
                    if self
                        .agent_screenshot_capture
                        .as_ref()
                        .is_some_and(|pending| pending.capture_request_id == failure.request_id)
                    {
                        self.agent_screenshot_capture = None;
                    } else {
                        self.invalidate_agent_screenshot_state(Some(
                            "A user-invoked Agent Help capture superseded the pending agent screenshot request.",
                        ));
                    }
                    self.agent_pending_image_attachment = None;
                    self.agent_status = failure.message.clone();
                    self.agent_help_capture_failure = Some(failure);
                    self.open_agent_assistant_dialog();
                }
            }
        }
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

    pub(super) fn agent_status_mentions_openai_quota(status: &str) -> bool {
        let lower = status.to_ascii_lowercase();
        lower.contains("openai")
            && (lower.contains("insufficient quota")
                || lower.contains("insufficient_quota")
                || lower.contains("quota or billing"))
            && (status.contains(OPENAI_USAGE_URL) || status.contains(OPENAI_BILLING_URL))
    }

    pub(super) fn render_openai_quota_links(&self, ui: &mut egui::Ui) {
        ui.horizontal_wrapped(|ui| {
            ui.small("OpenAI quota links:");
            ui.hyperlink_to("Usage", OPENAI_USAGE_URL);
            ui.hyperlink_to("Billing", OPENAI_BILLING_URL);
        });
    }

    pub(super) fn render_agent_status_message(
        &self,
        ui: &mut egui::Ui,
        status: &str,
        monospace: bool,
    ) {
        if monospace {
            ui.monospace(status.to_string());
        } else {
            ui.small(status.to_string());
        }
        if Self::agent_status_mentions_openai_quota(status) {
            self.render_openai_quota_links(ui);
        }
    }

    pub(super) fn agent_response_clipboard_payload(invocation: &AgentInvocationOutcome) -> String {
        let raw = invocation.raw_stdout.trim();
        if !raw.is_empty() {
            return raw.to_string();
        }
        serde_json::to_string_pretty(&invocation.response)
            .unwrap_or_else(|_| invocation.response.assistant_message.clone())
    }

    fn render_agent_web_research(&self, ui: &mut egui::Ui, response: &AgentResponse) {
        let Some(research) = response.web_research.as_ref() else {
            return;
        };
        if research.searches.is_empty() && research.pages.is_empty() && research.warnings.is_empty()
        {
            return;
        }
        egui::CollapsingHeader::new(format!(
            "{} ({} / {})",
            self.tr("agent.web_sources"),
            research.searches.len(),
            research.pages.len()
        ))
        .default_open(false)
        .show(ui, |ui| {
            for search in &research.searches {
                ui.small(format!(
                    "{}: {}",
                    self.tr("agent.web_sources.search"),
                    search.query
                ));
            }
            let mut displayed_urls = BTreeSet::new();
            for page in &research.pages {
                let url = page.final_url.trim();
                if url.is_empty() || !displayed_urls.insert(url.to_string()) {
                    continue;
                }
                ui.hyperlink_to(
                    page.title
                        .as_deref()
                        .map(str::trim)
                        .filter(|title| !title.is_empty())
                        .unwrap_or(url),
                    url,
                );
            }
            if research.pages.is_empty() {
                for result in research
                    .searches
                    .iter()
                    .flat_map(|search| search.results.iter())
                {
                    let url = result.url.trim();
                    if url.is_empty() || !displayed_urls.insert(url.to_string()) {
                        continue;
                    }
                    ui.hyperlink_to(result.title.trim(), url);
                }
            }
            for warning in &research.warnings {
                ui.colored_label(
                    egui::Color32::from_rgb(180, 120, 50),
                    format!("{}: {warning}", self.tr("agent.web_sources.warning")),
                );
            }
        });
    }

    pub(super) fn agent_preflight_summary_status(preflight: &AgentSystemPreflight) -> &'static str {
        if !preflight.available {
            return "unavailable";
        }
        match preflight
            .live_probe
            .as_ref()
            .map(|probe| probe.status_class)
        {
            None | Some(AgentLiveProbeStatusClass::Ok) => "ok",
            Some(AgentLiveProbeStatusClass::ProviderError)
            | Some(AgentLiveProbeStatusClass::UnsupportedTransport) => "live_warning",
            Some(_) => "live_failed",
        }
    }

    pub(super) fn agent_preflight_overall_label(
        preflight: &AgentSystemPreflight,
    ) -> (&'static str, egui::Color32) {
        match Self::agent_preflight_summary_status(preflight) {
            "ok" => ("Ready", egui::Color32::from_rgb(60, 140, 80)),
            "live_warning" => ("Live test warning", egui::Color32::from_rgb(180, 120, 50)),
            "live_failed" => ("Live test failed", egui::Color32::from_rgb(190, 70, 70)),
            _ => ("Unavailable", egui::Color32::from_rgb(190, 70, 70)),
        }
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

    pub(super) fn refresh_agent_token_file_credentials(&mut self) {
        let home = env::var_os("HOME").map(PathBuf::from);
        self.agent_token_file_credentials = load_agent_token_file_credentials(home.as_deref());
        self.agent_token_file_credentials_loaded = true;
    }

    pub(super) fn select_agent_system_and_reset_setup(&mut self, system_id: &str) -> bool {
        if self.agent_system_id == system_id {
            return false;
        }
        self.agent_system_id = system_id.to_string();
        self.clear_agent_preflight_output();
        self.clear_agent_model_discovery_snapshot();
        true
    }

    fn select_agent_system_and_persist_setup(&mut self, system_id: &str) {
        if !self.select_agent_system_and_reset_setup(system_id) {
            return;
        }
        if let Err(err) = self.persist_agent_system_selection_to_disk() {
            self.agent_status = format!("Could not persist selected agent system: {err}");
        }
    }

    pub(super) fn selected_agent_discovered_model(&self) -> Option<String> {
        let mut normalized_models = Vec::new();
        for model in &self.agent_discovered_models {
            if let Some(model) = normalize_agent_model_name(model)
                && !normalized_models.iter().any(|item| item == &model)
            {
                normalized_models.push(model);
            }
        }
        if normalized_models.is_empty() {
            return None;
        }
        if let Some(picked) = normalize_agent_model_name(&self.agent_discovered_model_pick)
            && normalized_models.iter().any(|item| item == &picked)
        {
            return Some(picked);
        }
        if normalized_models.len() == 1 {
            return normalized_models.into_iter().next();
        }
        None
    }

    pub(super) fn should_render_agent_model_selector(&self, system: &AgentSystemSpec) -> bool {
        agent_system_supports_model_discovery(system)
            && (matches!(system.transport, AgentSystemTransport::ExternalJsonStdio)
                || !self.agent_discovered_models.is_empty())
    }

    pub(super) fn agent_model_selection_prompt(
        &self,
        system: &AgentSystemSpec,
    ) -> Option<&'static str> {
        if !matches!(
            system.transport,
            AgentSystemTransport::NativeOpenai
                | AgentSystemTransport::NativeAnthropic
                | AgentSystemTransport::NativeMistral
                | AgentSystemTransport::NativeOpenaiCompat
        ) {
            return None;
        }
        if normalize_agent_model_name(self.agent_model_override.trim()).is_some()
            || system
                .model
                .as_deref()
                .and_then(normalize_agent_model_name)
                .is_some()
            || self.selected_agent_discovered_model().is_some()
        {
            return None;
        }
        let mut discovered_count = 0usize;
        let mut normalized_models = Vec::new();
        for model in &self.agent_discovered_models {
            if let Some(model) = normalize_agent_model_name(model)
                && !normalized_models.iter().any(|item| item == &model)
            {
                normalized_models.push(model);
                discovered_count += 1;
            }
        }
        (discovered_count > 1).then_some(Self::AGENT_MODEL_SELECTION_REQUIRED_MESSAGE)
    }

    pub(super) fn selected_agent_session_env_overrides(
        &self,
        system: &AgentSystemSpec,
    ) -> Result<HashMap<String, String>, String> {
        let mut overrides = HashMap::new();
        let session_api_key = self.agent_openai_api_key.trim();
        if let Some((key_env, _token_path)) = agent_api_key_source(system.transport) {
            if !session_api_key.is_empty() {
                overrides.insert(key_env.to_string(), session_api_key.to_string());
            } else if std::env::var(key_env)
                .ok()
                .map(|value| value.trim().is_empty())
                .unwrap_or(true)
                && let Some(value) = self
                    .agent_token_file_credentials
                    .get(key_env)
                    .and_then(AgentTokenFileCredential::value)
            {
                overrides.insert(key_env.to_string(), value.to_string());
            }
        } else if matches!(
            system.transport,
            AgentSystemTransport::NativeOpenaiCompat | AgentSystemTransport::ExternalJsonStdio
        ) && !is_codex_local_agent_system(system)
            && !is_pi_local_agent_system(system)
            && !session_api_key.is_empty()
        {
            overrides.insert(OPENAI_API_KEY_ENV.to_string(), session_api_key.to_string());
        }
        if system.supports_web_research && self.agent_allow_web_research {
            overrides.insert(AGENT_ALLOW_WEB_RESEARCH_ENV.to_string(), "1".to_string());
        }
        let override_base_url = self.agent_base_url_override.trim();
        if !override_base_url.is_empty()
            && matches!(
                system.transport,
                AgentSystemTransport::NativeOpenai
                    | AgentSystemTransport::NativeAnthropic
                    | AgentSystemTransport::NativeMistral
                    | AgentSystemTransport::NativeOpenaiCompat
            )
        {
            overrides.insert(
                AGENT_BASE_URL_ENV.to_string(),
                override_base_url.to_string(),
            );
        }
        let selected_discovered_model = self.selected_agent_discovered_model();
        let override_model = normalize_agent_model_name(self.agent_model_override.trim())
            .or(selected_discovered_model);
        if let Some(override_model) = override_model
            && agent_system_supports_model_selection(system)
        {
            overrides.insert(AGENT_MODEL_ENV.to_string(), override_model);
        }
        if let Some(timeout_override) = self.parse_agent_timeout_seconds()?
            && matches!(
                system.transport,
                AgentSystemTransport::ExternalJsonStdio
                    | AgentSystemTransport::NativeOpenai
                    | AgentSystemTransport::NativeAnthropic
                    | AgentSystemTransport::NativeMistral
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
                    | AgentSystemTransport::NativeMistral
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
                    | AgentSystemTransport::NativeMistral
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
                    | AgentSystemTransport::NativeMistral
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
                    | AgentSystemTransport::NativeMistral
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
                | AgentSystemTransport::NativeMistral
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
            AgentSystemTransport::NativeMistral => GUI_MISTRAL_DEFAULT_BASE_URL.to_string(),
            AgentSystemTransport::NativeOpenaiCompat => {
                GUI_OPENAI_COMPAT_DEFAULT_BASE_URL.to_string()
            }
            _ => return None,
        })
    }

    pub(super) fn selected_agent_base_url_placeholder(&self) -> String {
        let Some(system) = self.selected_agent_system() else {
            return GUI_OPENAI_COMPAT_DEFAULT_BASE_URL.to_string();
        };
        if !matches!(
            system.transport,
            AgentSystemTransport::NativeOpenai
                | AgentSystemTransport::NativeAnthropic
                | AgentSystemTransport::NativeMistral
                | AgentSystemTransport::NativeOpenaiCompat
        ) {
            return GUI_OPENAI_COMPAT_DEFAULT_BASE_URL.to_string();
        }
        if let Some(catalog_base_url) = system
            .base_url
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            return catalog_base_url.to_string();
        }
        match system.transport {
            AgentSystemTransport::NativeOpenai => GUI_OPENAI_DEFAULT_BASE_URL.to_string(),
            AgentSystemTransport::NativeAnthropic => GUI_ANTHROPIC_DEFAULT_BASE_URL.to_string(),
            AgentSystemTransport::NativeMistral => GUI_MISTRAL_DEFAULT_BASE_URL.to_string(),
            AgentSystemTransport::NativeOpenaiCompat => {
                GUI_OPENAI_COMPAT_DEFAULT_BASE_URL.to_string()
            }
            _ => GUI_OPENAI_COMPAT_DEFAULT_BASE_URL.to_string(),
        }
    }

    pub(super) fn selected_agent_model_discovery_source_key(
        &self,
        system: &AgentSystemSpec,
    ) -> Option<String> {
        if matches!(system.transport, AgentSystemTransport::ExternalJsonStdio)
            && agent_system_supports_model_discovery(system)
        {
            let source = if is_pi_local_agent_system(system) {
                "local-pi-model-list"
            } else {
                "local-codex-model-cache"
            };
            return Some(format!(
                "{}|{}|{}",
                system.id,
                system.transport.as_str(),
                source
            ));
        }
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
    ) -> String {
        if system.transport == AgentSystemTransport::NativeOpenaiCompat {
            return if self.agent_openai_api_key.trim().is_empty() {
                "optional-no-key".to_string()
            } else {
                "session-key".to_string()
            };
        }
        let Some((env_key, token_path)) = agent_api_key_source(system.transport) else {
            return "not-applicable".to_string();
        };
        if !self.agent_openai_api_key.trim().is_empty() {
            return "session-key".to_string();
        }
        if std::env::var(env_key)
            .ok()
            .map(|value| !value.trim().is_empty())
            .unwrap_or(false)
        {
            return match system.transport {
                AgentSystemTransport::NativeAnthropic => "env-anthropic-api-key".to_string(),
                AgentSystemTransport::NativeMistral => "env-mistral-api-key".to_string(),
                _ => "env-openai-api-key".to_string(),
            };
        }
        if self
            .agent_token_file_credentials
            .get(env_key)
            .and_then(AgentTokenFileCredential::value)
            .is_some()
        {
            return format!("file:{token_path}");
        }
        "no-key".to_string()
    }

    fn agent_credential_text(&self, key: &str, replacements: &[(&str, &str)]) -> String {
        let mut text = self.tr(key);
        for (name, value) in replacements {
            text = text.replace(&format!("{{{name}}}"), value);
        }
        text
    }

    fn selected_agent_credential_messages(&self, system: &AgentSystemSpec) -> Vec<(String, bool)> {
        if is_codex_local_agent_system(system) {
            return vec![(self.tr("agent.credential.codex_login"), false)];
        }
        if is_pi_local_agent_system(system) {
            return vec![(self.tr("agent.credential.pi_login"), false)];
        }
        if system.transport == AgentSystemTransport::NativeOpenaiCompat {
            return vec![(self.tr("agent.credential.optional"), false)];
        }
        let Some((env_key, token_path)) = agent_api_key_source(system.transport) else {
            return Vec::new();
        };
        if !self.agent_openai_api_key.trim().is_empty() {
            return vec![(self.tr("agent.credential.session"), false)];
        }
        if std::env::var(env_key)
            .ok()
            .map(|value| !value.trim().is_empty())
            .unwrap_or(false)
        {
            return vec![(
                self.agent_credential_text("agent.credential.environment", &[("env", env_key)]),
                false,
            )];
        }
        let Some(file) = self.agent_token_file_credentials.get(env_key) else {
            return vec![(
                self.agent_credential_text(
                    "agent.credential.missing",
                    &[("env", env_key), ("path", token_path)],
                ),
                true,
            )];
        };
        let mut messages = Vec::new();
        match file.status {
            AgentTokenFileStatus::Loaded => {
                messages.push((
                    self.agent_credential_text(
                        "agent.credential.file",
                        &[("path", file.display_path), ("env", env_key)],
                    ),
                    false,
                ));
                if system.transport == AgentSystemTransport::NativeAnthropic
                    && let Some(warning) = file.value().and_then(anthropic_api_key_kind_warning)
                {
                    messages.push((warning.to_string(), true));
                }
            }
            AgentTokenFileStatus::Missing => messages.push((
                self.agent_credential_text(
                    "agent.credential.missing",
                    &[("env", env_key), ("path", file.display_path)],
                ),
                true,
            )),
            AgentTokenFileStatus::Empty
            | AgentTokenFileStatus::Invalid
            | AgentTokenFileStatus::Unreadable => messages.push((
                self.agent_credential_text(
                    "agent.credential.file_problem",
                    &[("path", file.display_path), ("reason", &file.detail)],
                ),
                true,
            )),
        }
        if file.broadly_readable {
            messages.push((
                self.agent_credential_text(
                    "agent.credential.permissions",
                    &[("path", file.display_path)],
                ),
                true,
            ));
        }
        messages
    }

    pub(super) fn agent_model_discovery_failure_hint(error: &str) -> Option<&'static str> {
        let lower = error.to_ascii_lowercase();
        let auth_failed = lower.contains("401")
            || lower.contains("403")
            || lower.contains("unauthorized")
            || lower.contains("invalid_api_key")
            || lower.contains("incorrect api key")
            || lower.contains("authentication_error");
        if auth_failed && (lower.contains("mistral") || lower.contains("la plateforme")) {
            return Some(MISTRAL_API_KEY_AUTH_HINT);
        }
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
            return Some(
                "Authentication failed. Use an OpenAI Platform API key for OPENAI_API_KEY; ChatGPT/Codex subscription tokens are not OpenAI API keys.",
            );
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
        is_pi_local_agent_system(system)
            || matches!(
                system.transport,
                AgentSystemTransport::NativeOpenai
                    | AgentSystemTransport::NativeAnthropic
                    | AgentSystemTransport::NativeMistral
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

    pub(super) fn render_copyable_command_line(ui: &mut Ui, label: &str, command: &str) -> bool {
        let mut copied = false;
        ui.horizontal_wrapped(|ui| {
            if !label.trim().is_empty() {
                ui.small(label);
            }
            ui.monospace(command);
            if ui
                .small_button("⧉")
                .on_hover_text("Copy command to clipboard")
                .clicked()
            {
                ui.ctx().copy_text(command.to_string());
                copied = true;
            }
        });
        copied
    }

    pub(super) fn agent_response_sanity_warnings(
        invocation: &AgentInvocationOutcome,
    ) -> Vec<String> {
        let prompt = invocation
            .request
            .get("prompt")
            .and_then(|value| value.as_str())
            .unwrap_or_default();
        Self::agent_response_sanity_warnings_for_prompt(prompt, &invocation.response)
    }

    pub(super) fn agent_suggestion_run_blocker(
        command: &str,
        execution: AgentExecutionIntent,
    ) -> Option<String> {
        if execution == AgentExecutionIntent::Chat {
            return Some("This suggestion is explanatory (execution=chat).".to_string());
        }
        let command = command.trim();
        if command.is_empty() {
            return Some("The suggested command is empty.".to_string());
        }
        parse_shell_line(command).err().map(|err| {
            format!(
                "Invalid GENtle command: {} Use /help to inspect supported commands.",
                Self::compact_agent_validation_message(&err.to_string())
            )
        })
    }

    pub(super) fn agent_suggestion_precondition_blocker(
        &self,
        expr: Option<&serde_json::Value>,
    ) -> Option<String> {
        let expr = expr?;
        let expression = match serde_json::from_value::<crate::engine::FactExpression>(expr.clone())
        {
            Ok(expression) => expression,
            Err(err) => {
                return Some(format!(
                    "Cannot verify preconditions: invalid fact expression ({err})."
                ));
            }
        };
        let evaluation = match self.engine.read() {
            Ok(engine) => engine.evaluate_fact_expression(&expression, &[]),
            Err(_) => {
                return Some(
                    "Cannot verify preconditions while the project state is unavailable."
                        .to_string(),
                );
            }
        };
        if evaluation.truth == crate::engine::FactTruth::Satisfied {
            None
        } else {
            Some(format!(
                "Waiting for preconditions: {}.",
                crate::agent_bridge::agent_fact_readiness_label(&evaluation)
            ))
        }
    }

    pub(super) fn agent_suggestion_live_blocker(
        &self,
        suggestion: &AgentSuggestedCommand,
    ) -> Option<String> {
        Self::agent_suggestion_run_blocker(&suggestion.command, suggestion.execution).or_else(
            || self.agent_suggestion_precondition_blocker(suggestion.precondition_expr.as_ref()),
        )
    }

    pub(super) fn compact_agent_validation_message(message: &str) -> String {
        let cutoff = [
            " Supported GENtle-local alternatives:",
            " Supported commands:",
            " Details:",
        ]
        .into_iter()
        .filter_map(|marker| message.find(marker))
        .min()
        .unwrap_or(message.len());
        message[..cutoff].trim().to_string()
    }

    pub(super) fn agent_response_sanity_warnings_for_prompt(
        prompt: &str,
        response: &AgentResponse,
    ) -> Vec<String> {
        let mut warnings = Vec::new();
        let chat_only_count = response
            .suggested_commands
            .iter()
            .filter(|suggestion| suggestion.execution == AgentExecutionIntent::Chat)
            .count();
        if chat_only_count == response.suggested_commands.len()
            && !response.suggested_commands.is_empty()
        {
            warnings.push(
                "All suggestions are marked execution=chat, so GENtle will not run them. Runnable suggestions should use execution=ask."
                    .to_string(),
            );
        } else if chat_only_count > 0 {
            warnings.push(format!(
                "{chat_only_count} suggestion(s) are marked execution=chat, so GENtle will not run those rows."
            ));
        }

        for (idx, suggestion) in response.suggested_commands.iter().enumerate() {
            let index_1based = idx + 1;
            let command = suggestion.command.trim();
            if command.is_empty() {
                warnings.push(format!("Suggestion #{index_1based} has an empty command."));
                continue;
            }
            if Self::agent_command_has_placeholder(command) {
                warnings.push(format!(
                    "Suggestion #{index_1based} contains placeholder/help syntax rather than an executable command."
                ));
            }
            let parsed_command = parse_shell_line(command);
            if suggestion.execution != AgentExecutionIntent::Chat
                && let Err(err) = &parsed_command
            {
                warnings.push(format!(
                    "Suggestion #{index_1based} is not parseable by GENtle: {err}"
                ));
            }
            if suggestion.execution == AgentExecutionIntent::Auto
                && parsed_command.as_ref().is_ok_and(|parsed| {
                    matches!(
                        parsed,
                        ShellCommand::HistoryUndo | ShellCommand::HistoryRedo
                    )
                })
            {
                warnings.push(format!(
                    "Suggestion #{index_1based} cannot auto-run: {AGENT_HISTORY_CONFIRMATION_REQUIRED}. Click Run to confirm it explicitly."
                ));
            }
            if command.eq_ignore_ascii_case("/list") {
                let described_text = format!(
                    "{}\n{}\n{}",
                    response.assistant_message,
                    suggestion.title.clone().unwrap_or_default(),
                    suggestion.rationale.clone().unwrap_or_default()
                )
                .to_ascii_lowercase();
                if (described_text.contains("file")
                    || described_text.contains("folder")
                    || described_text.contains("director"))
                    && (described_text.contains("list") || described_text.contains("/list"))
                {
                    warnings.push(
                        "Suggestion #{index_1based} describes /list like a filesystem command; in GENtle, /list reports project state and loaded sequences."
                            .to_string(),
                    );
                }
            }
        }

        let lower_prompt = prompt.to_ascii_lowercase();
        let lower_reply = format!(
            "{}\n{}",
            response.assistant_message,
            response.questions.join("\n")
        )
        .to_ascii_lowercase();
        if !prompt.trim().is_empty()
            && (lower_reply.contains("what would you like to do")
                || lower_reply.contains("what would you like me to do"))
        {
            warnings.push(
                "The reply asks what to do even though the prompt already specified a task; the model likely ignored task context."
                    .to_string(),
            );
        }
        if Self::prompt_looks_like_retrieval_task(&lower_prompt)
            && !response
                .suggested_commands
                .iter()
                .any(|suggestion| Self::command_looks_like_retrieval(&suggestion.command))
        {
            warnings.push(
                "The prompt looks like a public-database or gene-retrieval task, but no suggestion uses a retrieval command such as /fetch, ensembl-gene, or genomes genes/extract-gene."
                    .to_string(),
            );
        }
        warnings
    }

    fn agent_command_has_placeholder(command: &str) -> bool {
        if command.contains('[')
            || command.contains(']')
            || command.contains('<')
            || command.contains('>')
            || command.contains("...")
        {
            return true;
        }
        command.split_whitespace().any(|token| {
            let bare = token.trim_matches(|ch: char| {
                matches!(ch, '\'' | '"' | ',' | ':' | ';' | '(' | ')' | '`')
            });
            matches!(
                bare,
                "ACCESSION"
                    | "CHR"
                    | "DNA"
                    | "END"
                    | "ENTRY_ID"
                    | "GENOME_ID"
                    | "ID"
                    | "MODEL"
                    | "PATH"
                    | "QUERY"
                    | "SEQ_ID"
                    | "SPECIES"
                    | "START"
                    | "SYSTEM_ID"
                    | "TEXT"
            )
        })
    }

    fn prompt_looks_like_retrieval_task(lower_prompt: &str) -> bool {
        (lower_prompt.contains("retrieve")
            || lower_prompt.contains("fetch")
            || lower_prompt.contains("database")
            || lower_prompt.contains("public database")
            || lower_prompt.contains("ensembl")
            || lower_prompt.contains("genbank")
            || lower_prompt.contains("ncbi")
            || lower_prompt.contains("uniprot"))
            && (lower_prompt.contains("gene")
                || lower_prompt.contains("isoform")
                || lower_prompt.contains("sequence")
                || lower_prompt.contains("protein")
                || lower_prompt.contains("fus"))
    }

    fn command_looks_like_retrieval(command: &str) -> bool {
        let lower = command.trim().to_ascii_lowercase();
        lower.starts_with("/fetch ")
            || lower.starts_with("ensembl-gene ")
            || lower.starts_with("ensembl-region ")
            || lower.starts_with("genbank ")
            || lower.starts_with("ncbi ")
            || lower.starts_with("uniprot ")
            || lower.starts_with("genomes genes ")
            || lower.starts_with("genomes extract-gene ")
            || lower.starts_with("helpers genes ")
            || lower.starts_with("helpers extract-gene ")
            || lower.contains("dbsnp")
    }

    pub(super) fn agent_prompt_direct_shell_command(prompt: &str) -> Option<&str> {
        let trimmed = prompt.trim();
        if trimmed.contains('\n') {
            return None;
        }
        if agent_path_is_supported_local_document(std::path::Path::new(trimmed)) {
            return None;
        }
        if trimmed.starts_with('/') || matches!(trimmed, "capabilities" | "help" | "state-summary")
        {
            Some(trimmed)
        } else {
            None
        }
    }

    fn agent_prompt_bare_absolute_path_hint(prompt: &str) -> Option<String> {
        let trimmed = prompt.trim();
        if trimmed.is_empty() || trimmed.contains('\n') || trimmed.split_whitespace().count() != 1 {
            return None;
        }
        if !std::path::Path::new(trimmed).is_absolute() {
            return None;
        }
        let without_root = trimmed.trim_start_matches('/');
        if !without_root.contains('/') && !without_root.contains('.') {
            return None;
        }
        Some(format!(
            "Prompt command looks like a bare absolute file path. GENtle's Agent Assistant does not use Ollama-style `/path/to/file` attachments; use `/open file {trimmed}` for sequence files, or import/attach the file through the GENtle GUI."
        ))
    }

    pub(super) fn agent_preflight_next_actions(preflight: &AgentSystemPreflight) -> Vec<String> {
        let key_hint = match preflight.transport.as_str() {
            transport if transport == AgentSystemTransport::NativeAnthropic.as_str() => {
                format!("Paste an Anthropic API key or set {ANTHROPIC_API_KEY_ENV}.")
            }
            transport if transport == AgentSystemTransport::NativeMistral.as_str() => {
                format!("Paste a Mistral API key or set {MISTRAL_API_KEY_ENV}.")
            }
            _ => format!(
                "Paste a session key or set {OPENAI_API_KEY_ENV}; ChatGPT/Codex subscriptions are not OpenAI API keys."
            ),
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
                    } else if preflight.transport == AgentSystemTransport::NativeMistral.as_str() {
                        vec![MISTRAL_API_KEY_AUTH_HINT.to_string()]
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
                AgentLiveProbeStatusClass::EndpointUnreachable => {
                    if preflight.system_id == "pi_local_stdio" {
                        vec![
                            "Install Pi, add it to PATH, or set PI_BIN to the Pi executable, then run Test Setup again."
                                .to_string(),
                        ]
                    } else {
                        vec![
                            "Start the local server or correct Base URL override, then run Test Setup again."
                                .to_string(),
                        ]
                    }
                }
                AgentLiveProbeStatusClass::UnsupportedTransport => vec![
                    "Use this setup check as config-only validation for this non-HTTP transport."
                        .to_string(),
                ],
                AgentLiveProbeStatusClass::ProviderError => {
                    if preflight.system_id == "pi_local_stdio" {
                        vec![
                            "Update Pi or check PI_BIN; GENtle's Pi bridge requires the documented no-tools/no-session flags."
                                .to_string(),
                        ]
                    } else {
                        vec![
                            "Inspect the provider response; model discovery must return JSON with model ids."
                                .to_string(),
                        ]
                    }
                }
            };
        }

        let mut actions = Vec::new();
        if preflight.warnings.iter().any(|warning| {
            warning.contains(OPENAI_API_KEY_ENV)
                || warning.contains(ANTHROPIC_API_KEY_ENV)
                || warning.contains(MISTRAL_API_KEY_ENV)
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
        if !agent_system_supports_model_discovery(system) {
            return;
        }
        let discovery_source =
            if matches!(system.transport, AgentSystemTransport::ExternalJsonStdio) {
                if is_pi_local_agent_system(system) {
                    "installed Pi CLI".to_string()
                } else {
                    "local Codex model metadata cache".to_string()
                }
            } else {
                let Some(base_url) = self.selected_agent_runtime_base_url(system) else {
                    return;
                };
                base_url
            };
        let Some(source_key) = self.selected_agent_model_discovery_source_key(system) else {
            return;
        };
        if !force {
            if let Some(task) = &self.agent_model_discovery_task
                && task.source_key == source_key
            {
                return;
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
        self.agent_model_discovery_status =
            format!("Discovering models from {discovery_source} (auth={key_label}) ...");
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
        let runtime_frame = Self::push_runtime_external_tool_frame(
            "agent model discovery",
            format!("{} from {}", system_id, discovery_source),
        );
        self.agent_model_discovery_task = Some(AgentModelDiscoveryTask {
            started: Instant::now(),
            source_key: source_key.clone(),
            runtime_frame,
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
                let status = Self::agent_preflight_summary_status(&preflight);
                let live_status = preflight
                    .live_probe
                    .as_ref()
                    .map(|probe| format!(", live={}", probe.status_class.as_str()))
                    .unwrap_or_default();
                self.agent_status = format!(
                    "Agent setup preflight: {} (status={}, transport={}{})",
                    selected_system.id, status, preflight.transport, live_status
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
        if self.agent_pending_image_attachment.is_some()
            && !selected_system.supports_image_attachments
        {
            self.agent_status = format!(
                "Selected agent system '{}' does not support image attachments. Choose an image-capable system or remove the screenshot.",
                selected_system.label
            );
            return;
        }
        let (available, reason) = self.selected_agent_system_availability(&selected_system);
        if !available {
            if let Some(prompt) = self.agent_model_selection_prompt(&selected_system) {
                self.agent_status = prompt.to_string();
            } else {
                self.agent_status = format!(
                    "Selected agent system is unavailable: {}",
                    reason.unwrap_or_else(|| "unknown reason".to_string())
                );
            }
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
                self.agent_status = self
                    .agent_model_selection_prompt(&selected_system)
                    .unwrap_or(
                        "Model is unspecified. Discover models and select one, or set Model override.",
                    )
                    .to_string();
                return;
            }
        }

        let include_state_summary = self.agent_include_state_summary;
        let conversation = self.agent_conversation.clone();
        let attachments = self
            .agent_pending_image_attachment
            .iter()
            .map(|attachment| attachment.request.clone())
            .collect::<Vec<_>>();
        let attachment_summaries = attachments
            .iter()
            .map(AgentAttachmentSummary::from)
            .collect::<Vec<_>>();
        let attachment_files = self
            .agent_pending_image_attachment
            .iter()
            .map(|attachment| attachment.temp_file.clone())
            .collect::<Vec<_>>();
        let worker_attachment_files = attachment_files.clone();
        let recent_project_paths = self.recent_project_paths.clone();
        let current_project_path = self.current_project_path.clone();
        let engine = self.engine.clone();
        let catalog_path = self.agent_catalog_path.trim().to_string();
        let job_id = self.alloc_background_job_id();
        let (tx, rx) = mpsc::channel::<AgentAskTaskMessage>();
        self.invalidate_agent_screenshot_state(Some(
            "A new agent request started, so the earlier screenshot request expired.",
        ));
        self.agent_last_command_output = None;
        self.agent_status = if let Some(timeout) = timeout_seconds {
            format!(
                "Starting agent '{}' in background (timeout={}s, retries={})",
                system_id,
                timeout,
                max_retries.unwrap_or(2)
            )
        } else {
            format!("Starting agent '{}' in background", system_id)
        };
        self.push_job_event(
            BackgroundJobKind::AgentAssist,
            BackgroundJobEventPhase::Started,
            Some(job_id),
            format!("Agent request started for system '{}'", system_id),
        );
        let runtime_frame = Self::push_runtime_background_job_frame(
            BackgroundJobKind::AgentAssist,
            job_id,
            format!("agent system '{system_id}'"),
        );
        self.agent_task = Some(AgentAskTask {
            job_id,
            prompt: prompt.clone(),
            attachment_summaries,
            _attachment_files: attachment_files,
            started: Instant::now(),
            runtime_frame,
            receiver: rx,
        });
        std::thread::spawn(move || {
            let _attachment_files = worker_attachment_files;
            let _ = tx.send(AgentAskTaskMessage::Status {
                job_id,
                message: "Building recent-project, tutorial, and Configuration context".to_string(),
            });
            let gui_context = GENtleApp::build_agent_gui_context_from(
                &recent_project_paths,
                current_project_path.as_deref(),
            );
            let request_context = if include_state_summary {
                let _ = tx.send(AgentAskTaskMessage::Status {
                    job_id,
                    message: "Building project summary and fact context for agent request"
                        .to_string(),
                });
                engine
                    .read()
                    .map(|guard| guard.clone_without_history())
                    .ok()
                    .map(|snapshot| {
                        let state_summary = snapshot.summarize_state();
                        let introspection =
                            build_agent_introspection_context(&snapshot.project_fact_graph());
                        (state_summary, introspection)
                    })
            } else {
                None
            };
            let _ = tx.send(AgentAskTaskMessage::Status {
                job_id,
                message: format!("Contacting agent system '{}'", system_id),
            });
            let result = invoke_agent_support_with_gui_context_and_attachments(
                Some(catalog_path.as_str()),
                &system_id,
                &prompt,
                request_context.as_ref().map(|(summary, _)| summary),
                request_context
                    .as_ref()
                    .map(|(_, introspection)| introspection),
                Some(&conversation),
                Some(&gui_context),
                &attachments,
                if env_overrides.is_empty() {
                    None
                } else {
                    Some(&env_overrides)
                },
            );
            let _ = tx.send(AgentAskTaskMessage::Done { job_id, result });
        });
    }

    pub(super) fn load_agent_conversation_from_state(&mut self) {
        let stored = self
            .engine
            .read()
            .ok()
            .and_then(|engine| {
                engine
                    .state()
                    .metadata
                    .get(AGENT_CONVERSATION_METADATA_KEY)
                    .cloned()
            })
            .and_then(|value| serde_json::from_value::<AgentConversation>(value).ok())
            .filter(|conversation| {
                conversation.schema == crate::agent_bridge::AGENT_CONVERSATION_SCHEMA
            })
            .map(AgentConversation::normalize)
            .unwrap_or_default();
        self.agent_conversation = stored;
    }

    fn persist_agent_conversation_to_state(&self) {
        let value = if self.agent_conversation.turns.is_empty() {
            None
        } else {
            serde_json::to_value(&self.agent_conversation).ok()
        };
        self.persist_project_metadata_values(&[(AGENT_CONVERSATION_METADATA_KEY, value)]);
    }

    pub(super) fn clear_agent_conversation(&mut self) {
        self.invalidate_agent_screenshot_state(None);
        self.agent_conversation = AgentConversation::default();
        self.agent_last_invocation = None;
        self.agent_pending_image_attachment = None;
        self.agent_help_capture_failure = None;
        self.persist_agent_conversation_to_state();
        self.agent_status = "Conversation and pending attachment cleared".to_string();
    }

    pub(super) fn execute_agent_suggested_command(
        &mut self,
        index_1based: usize,
        command_text: &str,
        trigger: &str,
    ) {
        let source_label = format!("Suggestion #{index_1based}");
        self.execute_agent_shell_command_from_ui(
            index_1based,
            &source_label,
            command_text,
            trigger,
        );
    }

    pub(super) fn execute_agent_suggestion(
        &mut self,
        index_1based: usize,
        suggestion: &AgentSuggestedCommand,
        trigger: &str,
    ) {
        if let Some(reason) = self.agent_suggestion_live_blocker(suggestion) {
            let source_label = format!("Suggestion #{index_1based}");
            self.agent_status = format!("{source_label} not run: {reason}");
            self.agent_execution_log.push(AgentCommandExecutionRecord {
                index_1based,
                command: suggestion.command.trim().to_string(),
                trigger: trigger.to_string(),
                ok: false,
                state_changed: false,
                summary: reason,
                executed_at_unix_ms: Self::now_unix_ms(),
            });
            if self.agent_execution_log.len() > 100 {
                let drain = self.agent_execution_log.len() - 100;
                self.agent_execution_log.drain(0..drain);
            }
            return;
        }
        self.execute_agent_suggested_command(index_1based, &suggestion.command, trigger);
    }

    pub(super) fn execute_agent_prompt_command(&mut self, command_text: &str) {
        self.execute_agent_shell_command_from_ui(0, "Prompt command", command_text, "prompt");
    }

    fn execute_agent_shell_command_from_ui(
        &mut self,
        index_1based: usize,
        source_label: &str,
        command_text: &str,
        trigger: &str,
    ) {
        self.agent_last_command_output = None;
        let trimmed = command_text.trim();
        if trimmed.is_empty() {
            self.agent_status = format!("{source_label} is empty");
            return;
        }
        if trigger == "prompt"
            && let Some(hint) = Self::agent_prompt_bare_absolute_path_hint(trimmed)
        {
            self.agent_status = hint.clone();
            self.agent_execution_log.push(AgentCommandExecutionRecord {
                index_1based,
                command: trimmed.to_string(),
                trigger: trigger.to_string(),
                ok: false,
                state_changed: false,
                summary: hint,
                executed_at_unix_ms: Self::now_unix_ms(),
            });
            if self.agent_execution_log.len() > 100 {
                let drain = self.agent_execution_log.len() - 100;
                self.agent_execution_log.drain(0..drain);
            }
            return;
        }
        let command = match parse_shell_line(trimmed) {
            Ok(command) => command,
            Err(err) => {
                self.agent_status = format!("{source_label} parse error: {err}");
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
        if trigger == "auto"
            && matches!(
                command,
                ShellCommand::HistoryUndo | ShellCommand::HistoryRedo
            )
        {
            let summary = AGENT_HISTORY_CONFIRMATION_REQUIRED.to_string();
            self.agent_status = format!("{source_label} rejected: {summary}");
            self.agent_execution_log.push(AgentCommandExecutionRecord {
                index_1based,
                command: trimmed.to_string(),
                trigger: trigger.to_string(),
                ok: false,
                state_changed: false,
                summary,
                executed_at_unix_ms: Self::now_unix_ms(),
            });
            return;
        }
        if matches!(
            command,
            ShellCommand::AgentsAsk { .. }
                | ShellCommand::AgentsPlan { .. }
                | ShellCommand::AgentsExecutePlan { .. }
        ) {
            self.agent_status = format!(
                "{source_label} rejected: agent-to-agent 'agents ...' commands are blocked"
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
        if matches!(
            command,
            ShellCommand::HistoryUndo | ShellCommand::HistoryRedo
        ) {
            let state_changed = match command {
                ShellCommand::HistoryUndo => self.undo_last_operation(),
                ShellCommand::HistoryRedo => self.redo_last_operation(),
                _ => unreachable!("history transition match is exhaustive"),
            };
            let summary = self.app_status.clone();
            self.agent_status = format!("{source_label}: {summary}");
            self.agent_execution_log.push(AgentCommandExecutionRecord {
                index_1based,
                command: trimmed.to_string(),
                trigger: trigger.to_string(),
                ok: state_changed,
                state_changed,
                summary,
                executed_at_unix_ms: Self::now_unix_ms(),
            });
            if self.agent_execution_log.len() > 100 {
                let drain = self.agent_execution_log.len() - 100;
                self.agent_execution_log.drain(0..drain);
            }
            return;
        }
        let suppress_auto_open = Self::agent_command_suppresses_auto_open(trimmed, &command);
        if let Some(summary) = self.try_apply_shell_ui_intent(&command) {
            self.agent_status = format!("{source_label}: {summary}");
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
                let mut opened_seq_ids = if matches!(command, ShellCommand::LoadFile { .. }) {
                    Self::agent_sequence_ids_from_shell_output(&run.output)
                } else {
                    Vec::new()
                };
                let mut extra_summary: Option<String> = match &command {
                    ShellCommand::Help { topic, .. } => {
                        self.open_help_doc(HelpDoc::Shell);
                        if !topic.is_empty() {
                            self.help_search_query = topic.join(" ");
                            self.help_search_selected = 0;
                            self.refresh_help_search_matches();
                        }
                        Some("opened Help > Shell Commands".to_string())
                    }
                    ShellCommand::StateSummary => {
                        Some("showing the current project summary below".to_string())
                    }
                    _ => None,
                };
                let mut effective_state_changed = run.state_changed;
                if let ShellCommand::EnsemblGeneFetch { entry_id, .. } = &command {
                    match self.import_agent_ensembl_gene_fetch_result(entry_id.as_deref(), &run) {
                        Ok(imported_seq_ids) if !imported_seq_ids.is_empty() => {
                            effective_state_changed = true;
                            opened_seq_ids.extend(imported_seq_ids.iter().cloned());
                            extra_summary = Some(format!(
                                "imported Ensembl gene sequence {}",
                                imported_seq_ids.join(", ")
                            ));
                        }
                        Ok(_) => {
                            extra_summary = Some(
                                "stored Ensembl gene metadata; no sequence was imported"
                                    .to_string(),
                            );
                        }
                        Err(err) => {
                            extra_summary = Some(format!(
                                "stored Ensembl gene metadata; sequence import failed: {err}"
                            ));
                        }
                    }
                }
                let mut opened_seq_ids_unique = Vec::new();
                for seq_id in opened_seq_ids {
                    if !opened_seq_ids_unique.contains(&seq_id) {
                        opened_seq_ids_unique.push(seq_id);
                    }
                }
                if !suppress_auto_open {
                    for seq_id in &opened_seq_ids_unique {
                        self.open_sequence_window(seq_id);
                    }
                }
                let summary = if let Some(extra) = extra_summary {
                    if opened_seq_ids_unique.is_empty() {
                        format!("success - {extra}")
                    } else if suppress_auto_open {
                        format!(
                            "success - {extra}; did not open {}",
                            opened_seq_ids_unique.join(", ")
                        )
                    } else {
                        format!(
                            "success - {extra}; opened {}",
                            opened_seq_ids_unique.join(", ")
                        )
                    }
                } else if effective_state_changed {
                    if opened_seq_ids_unique.is_empty() {
                        "success - executed (state changed)".to_string()
                    } else if suppress_auto_open {
                        format!(
                            "success - executed (state changed; did not open {})",
                            opened_seq_ids_unique.join(", ")
                        )
                    } else {
                        format!(
                            "success - executed (state changed; opened {})",
                            opened_seq_ids_unique.join(", ")
                        )
                    }
                } else {
                    "success - executed".to_string()
                };
                self.agent_status = format!("{source_label}: {summary}");
                self.agent_execution_log.push(AgentCommandExecutionRecord {
                    index_1based,
                    command: trimmed.to_string(),
                    trigger: trigger.to_string(),
                    ok: true,
                    state_changed: effective_state_changed,
                    summary,
                    executed_at_unix_ms: Self::now_unix_ms(),
                });
                self.agent_last_command_output = Some(AgentCommandOutput {
                    command: trimmed.to_string(),
                    output: run.output,
                    state_changed: effective_state_changed,
                });
            }
            Err(err) => {
                self.agent_status = format!("{source_label} failed: {err}");
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

    pub(super) fn agent_command_suppresses_auto_open(
        command_text: &str,
        command: &ShellCommand,
    ) -> bool {
        if !matches!(command, ShellCommand::EnsemblGeneFetch { .. }) {
            return false;
        }
        split_shell_words(command_text)
            .map(|tokens| tokens.iter().any(|token| token == "--no-open"))
            .unwrap_or(false)
    }

    fn agent_sequence_ids_from_shell_output(output: &serde_json::Value) -> Vec<String> {
        let Some(result) = output.get("result") else {
            return Vec::new();
        };
        ["created_seq_ids", "changed_seq_ids"]
            .into_iter()
            .filter_map(|key| result.get(key).and_then(|ids| ids.as_array()))
            .flat_map(|ids| ids.iter().filter_map(|value| value.as_str()))
            .map(str::to_string)
            .collect()
    }

    fn agent_ensembl_gene_entry_id_from_fetch_run(
        explicit_entry_id: Option<&str>,
        run: &ShellRunResult,
    ) -> Option<String> {
        explicit_entry_id
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(str::to_string)
            .or_else(|| {
                run.output
                    .get("result")
                    .and_then(|result| result.get("messages"))
                    .and_then(|messages| messages.as_array())
                    .and_then(|messages| {
                        messages.iter().find_map(|message| {
                            let message = message.as_str()?;
                            message
                                .strip_prefix("Fetched Ensembl gene '")
                                .and_then(|rest| {
                                    rest.split_once('\'').map(|(entry_id, _)| entry_id)
                                })
                                .map(str::to_string)
                        })
                    })
            })
    }

    pub(super) fn import_agent_ensembl_gene_fetch_result(
        &mut self,
        explicit_entry_id: Option<&str>,
        run: &ShellRunResult,
    ) -> Result<Vec<String>, String> {
        let Some(entry_id) =
            Self::agent_ensembl_gene_entry_id_from_fetch_run(explicit_entry_id, run)
        else {
            return Ok(Vec::new());
        };
        let import_result = {
            let mut guard = self.engine.write().unwrap();
            guard
                .apply(Operation::ImportEnsemblGeneSequence {
                    entry_id: entry_id.clone(),
                    output_id: Some(entry_id.clone()),
                })
                .map_err(|err| err.to_string())?
        };
        self.lineage_cache_valid = false;
        let output = serde_json::json!({ "result": import_result });
        Ok(Self::agent_sequence_ids_from_shell_output(&output))
    }

    pub(super) fn try_apply_shell_ui_intent(&mut self, command: &ShellCommand) -> Option<String> {
        if let ShellCommand::UiRecentProject { item_id } = command {
            return Some(self.apply_recent_project_intent(item_id));
        }
        if let ShellCommand::UiTutorialProject { chapter_id } = command {
            return Some(self.apply_tutorial_project_intent(chapter_id));
        }
        if let ShellCommand::UiConfiguration { action, section } = command {
            return Some(self.apply_configuration_intent(*action, *section));
        }
        if let ShellCommand::UiSequenceWindow { action, seq_id } = command {
            return Some(self.apply_sequence_window_intent(*action, seq_id));
        }
        if let ShellCommand::UiSequenceSelection {
            seq_id,
            start_0based,
            end_0based_exclusive,
        } = command
        {
            return Some(self.apply_sequence_selection_intent(
                seq_id,
                *start_0based,
                *end_0based_exclusive,
            ));
        }
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
        if matches!(action, UiIntentAction::Close) {
            return Some(self.apply_close_ui_intent_target(*target));
        }
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
            UiIntentTarget::OpenSequence => self.prompt_open_sequence(),
            UiIntentTarget::RecentProject => {
                self.app_status =
                    "Recent-project UI intent requires an item id from the current GUI host context"
                        .to_string();
            }
            UiIntentTarget::TutorialProject => {
                self.app_status =
                    "Tutorial-project UI intent requires a chapter id from the tutorial catalog"
                        .to_string();
            }
            UiIntentTarget::Configuration => self.open_configuration_dialog(),
            UiIntentTarget::PreparedReferences => self.open_reference_genome_inspector_dialog(),
            UiIntentTarget::PrepareReferenceGenome => self.open_reference_genome_prepare_dialog(),
            UiIntentTarget::RetrieveGenomeSequence => self.open_reference_genome_retrieve_dialog(),
            UiIntentTarget::BlastGenomeSequence => self.open_reference_genome_blast_dialog(),
            UiIntentTarget::ImportGenomeTrack => self.open_genome_bed_track_dialog(),
            UiIntentTarget::FeatureLocationEditor => self.open_feature_location_editor(),
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

    pub(super) fn recent_project_agent_item_id(path: &str) -> String {
        let normalized = Self::normalize_project_path(path);
        let digest = Self::digest_hex(&normalized);
        format!("recent-{}", &digest[..16.min(digest.len())])
    }

    pub(super) fn build_agent_gui_context_from(
        recent_project_paths: &[String],
        current_project_path: Option<&str>,
    ) -> AgentGuiContext {
        let current_project_path = current_project_path.map(Self::normalize_project_path);
        let recent_projects = recent_project_paths
            .iter()
            .take(crate::agent_bridge::AGENT_GUI_RECENT_PROJECT_LIMIT)
            .enumerate()
            .map(|(index, path)| {
                let parsed = Path::new(path);
                let file_name = parsed
                    .file_name()
                    .map(|value| value.to_string_lossy().to_string())
                    .unwrap_or_else(|| "saved project".to_string());
                let parent_label = parsed
                    .parent()
                    .and_then(Path::file_name)
                    .map(|value| value.to_string_lossy().to_string())
                    .unwrap_or_default();
                let display_label = if parent_label.is_empty() {
                    file_name.clone()
                } else {
                    format!("{file_name} ({parent_label})")
                };
                let metadata = fs::metadata(parsed).ok();
                let modified_at_unix_ms = metadata
                    .as_ref()
                    .and_then(|metadata| metadata.modified().ok())
                    .and_then(|modified| modified.duration_since(UNIX_EPOCH).ok())
                    .map(|duration| duration.as_millis().min(u64::MAX as u128) as u64);
                let item_id = Self::recent_project_agent_item_id(path);
                AgentGuiRecentProject {
                    item_id: item_id.clone(),
                    display_label,
                    file_name,
                    parent_label,
                    list_position: index + 1,
                    exists: parsed.is_file(),
                    byte_count: metadata.as_ref().map(std::fs::Metadata::len),
                    modified_at_unix_ms,
                    current_project: current_project_path
                        .as_deref()
                        .is_some_and(|current| current == Self::normalize_project_path(path)),
                    open_command: format!("ui open recent-project {item_id}"),
                }
            })
            .collect::<Vec<_>>();

        let mut context = AgentGuiContext {
            host_available: true,
            recent_project_count: recent_projects.len(),
            recent_projects,
            ..AgentGuiContext::default()
        };
        match Self::load_tutorial_project_entries() {
            Ok(entries) => {
                context.tutorial_project_count = entries.len();
                context.tutorial_projects = entries
                    .into_iter()
                    .take(AGENT_GUI_TUTORIAL_PROJECT_LIMIT)
                    .map(|entry| {
                        let display_label = Self::tutorial_display_label(
                            entry.decimal_id.as_deref(),
                            Some(entry.chapter_order),
                            &entry.chapter_title,
                        );
                        let chapter_id = entry.chapter_id;
                        AgentGuiTutorialProject {
                            open_command: format!("ui open tutorial-project {chapter_id}"),
                            chapter_id,
                            decimal_id: entry.decimal_id,
                            display_label,
                            title: entry.chapter_title,
                            summary: entry.chapter_summary,
                            group: entry.group_label,
                            tier: entry.tier.as_str().to_string(),
                            example_id: entry.example.id,
                            online: entry.example.test_mode == ExampleTestMode::Online,
                            review_status: entry.review_status,
                            review_stale: entry.review_stale,
                        }
                    })
                    .collect();
                context.included_tutorial_project_count = context.tutorial_projects.len();
                context.omitted_tutorial_project_count = context
                    .tutorial_project_count
                    .saturating_sub(context.included_tutorial_project_count);
                context.tutorial_projects_truncated = context.omitted_tutorial_project_count > 0;
            }
            Err(err) => context
                .warnings
                .push(format!("Could not load the GUI tutorial catalog: {err}")),
        }
        context.configuration_sections = UiConfigurationSection::all()
            .iter()
            .copied()
            .map(|section| AgentGuiConfigurationSection {
                section_id: section.as_str().to_string(),
                title: section.title().to_string(),
                detail: section.detail().to_string(),
                open_command: format!("ui open configuration {}", section.as_str()),
            })
            .collect();
        context
    }

    fn apply_recent_project_intent(&mut self, item_id: &str) -> String {
        let item_id = item_id.trim();
        let Some(path) = self
            .recent_project_paths
            .iter()
            .find(|path| Self::recent_project_agent_item_id(path) == item_id)
            .cloned()
        else {
            return format!(
                "ui intent open 'recent-project' found no current recent-project item '{item_id}'; ask the agent to list the current GUI context again"
            );
        };
        if !Path::new(&path).is_file() {
            return format!(
                "ui intent open 'recent-project' found item '{item_id}', but its saved project file is missing"
            );
        }
        let label = Self::recent_project_menu_label(&path);
        self.request_project_action(ProjectAction::OpenPath(path));
        format!("ui intent open 'recent-project' ({label})")
    }

    fn apply_tutorial_project_intent(&mut self, chapter_id: &str) -> String {
        let chapter_id = chapter_id.trim();
        let entries = match Self::load_tutorial_project_entries() {
            Ok(entries) => entries,
            Err(err) => {
                return format!(
                    "ui intent open 'tutorial-project' could not load the tutorial catalog: {err}"
                );
            }
        };
        let Some(entry) = entries.iter().find(|entry| entry.chapter_id == chapter_id) else {
            return format!(
                "ui intent open 'tutorial-project' found no current chapter '{chapter_id}'"
            );
        };
        let title = entry.chapter_title.clone();
        self.request_project_action(ProjectAction::OpenTutorialChapter(chapter_id.to_string()));
        format!("ui intent open 'tutorial-project' '{title}' ({chapter_id})")
    }

    fn apply_configuration_intent(
        &mut self,
        action: UiIntentAction,
        section: UiConfigurationSection,
    ) -> String {
        if matches!(action, UiIntentAction::Close) {
            let was_open = self.show_configuration_dialog;
            self.show_configuration_dialog = false;
            return if was_open {
                "ui intent close 'configuration'".to_string()
            } else {
                "ui intent close 'configuration' requested; target was already closed".to_string()
            };
        }
        let tab = match section {
            UiConfigurationSection::ExternalApplications => ConfigurationTab::ExternalApplications,
            UiConfigurationSection::AgentSystems => ConfigurationTab::AgentSystems,
            UiConfigurationSection::Microarrays => ConfigurationTab::Microarrays,
            UiConfigurationSection::Graphics => ConfigurationTab::Graphics,
            UiConfigurationSection::Language => ConfigurationTab::Language,
        };
        if matches!(section, UiConfigurationSection::AgentSystems) {
            self.refresh_agent_token_file_credentials();
        }
        self.open_configuration_dialog_for_tab(tab);
        format!(
            "ui intent {} 'configuration' section '{}'",
            action.as_str(),
            section.as_str()
        )
    }

    fn apply_close_ui_intent_target(&mut self, target: UiIntentTarget) -> String {
        let was_open = match target {
            UiIntentTarget::OpenSequence => {
                return "ui intent close 'open-sequence' is not applicable; use ui close sequence-window SEQ_ID for DNA viewers".to_string();
            }
            UiIntentTarget::RecentProject | UiIntentTarget::TutorialProject => {
                return format!(
                    "ui intent close '{}' is not applicable; opening a project is an action, not a persistent dialog",
                    target.as_str()
                );
            }
            UiIntentTarget::Configuration => {
                let was_open = self.show_configuration_dialog;
                self.show_configuration_dialog = false;
                was_open
            }
            UiIntentTarget::PreparedReferences => {
                let was_open = self.show_reference_genome_inspector_dialog;
                self.show_reference_genome_inspector_dialog = false;
                was_open
            }
            UiIntentTarget::PrepareReferenceGenome | UiIntentTarget::PrepareHelperGenome => {
                let was_open = self.show_reference_genome_prepare_dialog;
                self.show_reference_genome_prepare_dialog = false;
                was_open
            }
            UiIntentTarget::RetrieveGenomeSequence | UiIntentTarget::RetrieveHelperSequence => {
                let was_open = self.show_reference_genome_retrieve_dialog;
                self.show_reference_genome_retrieve_dialog = false;
                was_open
            }
            UiIntentTarget::BlastGenomeSequence | UiIntentTarget::BlastHelperSequence => {
                let was_open = self.show_reference_genome_blast_dialog;
                self.show_reference_genome_blast_dialog = false;
                was_open
            }
            UiIntentTarget::ImportGenomeTrack => {
                let was_open = self.show_genome_bed_track_dialog;
                self.show_genome_bed_track_dialog = false;
                was_open
            }
            UiIntentTarget::FeatureLocationEditor => self.close_feature_location_editor(),
            UiIntentTarget::PcrDesign => {
                let was_open = self.show_pcr_design_dialog;
                self.show_pcr_design_dialog = false;
                was_open
            }
            UiIntentTarget::SequencingConfirmation => {
                let was_open = self.show_sequencing_confirmation_dialog;
                self.show_sequencing_confirmation_dialog = false;
                was_open
            }
            UiIntentTarget::AgentAssistant => {
                let was_open = self.show_agent_assistant_dialog;
                self.show_agent_assistant_dialog = false;
                was_open
            }
        };
        if was_open {
            format!("ui intent close '{}'", target.as_str())
        } else {
            format!(
                "ui intent close '{}' requested; target was already closed",
                target.as_str()
            )
        }
    }

    fn apply_sequence_window_intent(&mut self, action: UiIntentAction, seq_id: &str) -> String {
        match action {
            UiIntentAction::Open | UiIntentAction::Focus => {
                self.apply_open_or_focus_sequence_window_intent(action, seq_id)
            }
            UiIntentAction::Close => self.apply_close_sequence_window_intent(seq_id),
        }
    }

    fn apply_open_or_focus_sequence_window_intent(
        &mut self,
        action: UiIntentAction,
        seq_id: &str,
    ) -> String {
        let seq_id = seq_id.trim();
        if seq_id.is_empty() {
            return format!(
                "ui intent {} 'sequence-window' requires seq_id",
                action.as_str()
            );
        }
        let sequence_loaded = self
            .engine
            .read()
            .unwrap()
            .state()
            .sequences
            .contains_key(seq_id);
        if !sequence_loaded {
            return format!(
                "ui intent {} 'sequence-window' found no loaded sequence {seq_id}",
                action.as_str()
            );
        }
        let existing_viewport = self.find_open_sequence_viewport_id(seq_id);
        let existing_is_opening = existing_viewport
            .and_then(|viewport_id| self.windows.get(&viewport_id))
            .and_then(|window| window.read().ok())
            .is_some_and(|window| window.is_sequence_opening());
        let pending_is_opening = self
            .new_windows
            .iter()
            .any(|window| window.sequence_id().as_deref() == Some(seq_id));
        self.open_sequence_window(seq_id);
        if existing_is_opening || pending_is_opening {
            format!(
                "DNA Sequence Viewer for '{seq_id}' is still opening (sequence record kept loaded)"
            )
        } else if existing_viewport.is_some() {
            format!(
                "Focused the existing DNA Sequence Viewer for '{seq_id}' (sequence record kept loaded)"
            )
        } else {
            format!(
                "Opening DNA Sequence Viewer for '{seq_id}' (queued; sequence record kept loaded)"
            )
        }
    }

    fn apply_close_sequence_window_intent(&mut self, seq_id: &str) -> String {
        let seq_id = seq_id.trim();
        let mut close_requested = false;
        if let Some(viewport_id) = self.find_open_sequence_viewport_id(seq_id) {
            if let Ok(mut to_close) = self.windows_to_close.write() {
                if !to_close.contains(&viewport_id) {
                    to_close.push(viewport_id);
                }
                close_requested = true;
            }
            self.pending_focus_viewports.retain(|id| *id != viewport_id);
            if close_requested {
                self.process_window_close_queue();
            }
        }
        let pending_before = self.new_windows.len();
        self.new_windows
            .retain(|window| window.sequence_id().as_deref() != Some(seq_id));
        let removed_pending = self.new_windows.len() != pending_before;
        let sequence_loaded = self
            .engine
            .read()
            .unwrap()
            .state()
            .sequences
            .contains_key(seq_id);
        if close_requested || removed_pending {
            format!(
                "ui intent close 'sequence-window' (seq_id={seq_id}; sequence record kept loaded)"
            )
        } else if sequence_loaded {
            format!(
                "ui intent close 'sequence-window' found no open window for seq_id={seq_id}; sequence record kept loaded"
            )
        } else {
            format!("ui intent close 'sequence-window' found no open or loaded sequence {seq_id}")
        }
    }

    fn current_sequence_window_selection_range(&self, seq_id: &str) -> Option<(usize, usize)> {
        if let Some(viewport_id) = self.find_open_sequence_viewport_id(seq_id)
            && let Some(window) = self.windows.get(&viewport_id)
            && let Ok(window) = window.read()
        {
            return window.selection_range_0based();
        }
        self.new_windows
            .iter()
            .find(|window| window.sequence_id().as_deref() == Some(seq_id))
            .and_then(|window| window.selection_range_0based())
    }

    fn set_sequence_window_selection_range(
        &mut self,
        seq_id: &str,
        start_0based: usize,
        end_0based_exclusive: usize,
    ) -> Result<(usize, usize), String> {
        if let Some(viewport_id) = self.find_open_sequence_viewport_id(seq_id) {
            let Some(window) = self.windows.get(&viewport_id).cloned() else {
                return Err(format!("Open sequence window for {seq_id} disappeared"));
            };
            window
                .write()
                .map_err(|_| "Sequence window lock poisoned while setting selection".to_string())?
                .set_selection_range_0based(start_0based, end_0based_exclusive)?;
            self.queue_focus_viewport(viewport_id);
            return Ok((start_0based, end_0based_exclusive));
        }
        if let Some(window) = self.find_pending_sequence_window_mut(seq_id) {
            window.set_selection_range_0based(start_0based, end_0based_exclusive)?;
            return Ok((start_0based, end_0based_exclusive));
        }
        let dna = self
            .engine
            .read()
            .map_err(|_| "Engine lock poisoned while opening sequence window".to_string())?
            .state()
            .sequences
            .get(seq_id)
            .cloned()
            .ok_or_else(|| format!("No loaded sequence {seq_id}"))?;
        let mut window = Window::new_dna(dna, seq_id.to_string(), self.engine.clone());
        window.set_selection_range_0based(start_0based, end_0based_exclusive)?;
        self.new_windows.push(window);
        Ok((start_0based, end_0based_exclusive))
    }

    fn apply_sequence_selection_intent(
        &mut self,
        seq_id: &str,
        start_0based: Option<usize>,
        end_0based_exclusive: Option<usize>,
    ) -> String {
        let seq_id = seq_id.trim();
        if seq_id.is_empty() {
            return "ui selection 'sequence-window' requires seq_id".to_string();
        }
        match (start_0based, end_0based_exclusive) {
            (Some(start), Some(end)) => {
                match self.set_sequence_window_selection_range(seq_id, start, end) {
                    Ok((start, end)) => format!(
                        "ui selection 'sequence-window' set seq_id={seq_id} range={start}..{end} (0-based, end-exclusive; sequence record kept loaded)"
                    ),
                    Err(err) => {
                        format!("ui selection 'sequence-window' failed for seq_id={seq_id}: {err}")
                    }
                }
            }
            (None, None) => {
                if let Some((start, end)) = self.current_sequence_window_selection_range(seq_id) {
                    format!(
                        "ui selection 'sequence-window' seq_id={seq_id} range={start}..{end} (0-based, end-exclusive)"
                    )
                } else {
                    let sequence_loaded = self
                        .engine
                        .read()
                        .map(|engine| engine.state().sequences.contains_key(seq_id))
                        .unwrap_or(false);
                    if sequence_loaded {
                        format!(
                            "ui selection 'sequence-window' found no current selection for loaded seq_id={seq_id}"
                        )
                    } else {
                        format!("ui selection 'sequence-window' found no loaded sequence {seq_id}")
                    }
                }
            }
            _ => {
                "ui selection 'sequence-window' requires both start_0based and end_0based_exclusive"
                    .to_string()
            }
        }
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
                self.execute_agent_suggestion(idx + 1, suggestion, "auto");
            }
        }
    }

    pub(super) fn poll_agent_assistant_task(&mut self, ctx: &egui::Context) {
        if self.agent_task.is_none() {
            return;
        }
        ctx.request_repaint_after(Duration::from_millis(100));
        let mut done: Option<(u64, Result<AgentInvocationOutcome, String>)> = None;
        let mut latest_status: Option<(String, f64)> = None;
        if let Some(task) = &self.agent_task {
            loop {
                match task.receiver.try_recv() {
                    Ok(AgentAskTaskMessage::Status { job_id, message }) => {
                        if job_id == task.job_id {
                            task.runtime_frame.update_phase("agent_status");
                            task.runtime_frame.update_detail(message.clone());
                            latest_status = Some((message, task.started.elapsed().as_secs_f64()));
                        }
                    }
                    Ok(AgentAskTaskMessage::Done { job_id, result }) => {
                        if job_id == task.job_id {
                            done = Some((job_id, result));
                        }
                        break;
                    }
                    Err(mpsc::TryRecvError::Empty) => break,
                    Err(mpsc::TryRecvError::Disconnected) => {
                        done = Some((task.job_id, Err("Agent worker disconnected".to_string())));
                        break;
                    }
                }
            }
        }
        if let Some((message, elapsed)) = latest_status
            && done.is_none()
        {
            self.agent_status = format!("{message} ({elapsed:.1}s)");
        }
        if let Some((job_id, outcome)) = done {
            let elapsed = self
                .agent_task
                .as_ref()
                .map(|task| task.started.elapsed().as_secs_f64())
                .unwrap_or(0.0);
            let completed_prompt = self
                .agent_task
                .as_ref()
                .map(|task| task.prompt.clone())
                .unwrap_or_default();
            let completed_attachments = self
                .agent_task
                .as_ref()
                .map(|task| task.attachment_summaries.clone())
                .unwrap_or_default();
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
                    let completed_at_unix_ms = Self::now_unix_ms();
                    self.agent_conversation.push_turn(AgentConversationTurn {
                        user_message: completed_prompt,
                        response: response.clone(),
                        attachments: completed_attachments,
                        system_id: invocation.system_id.clone(),
                        system_label: invocation.system_label.clone(),
                        completed_at_unix_ms,
                    });
                    self.persist_agent_conversation_to_state();
                    if let Some(request) = response.screenshot_request.clone() {
                        self.activate_agent_screenshot_consent(
                            request,
                            invocation.system_id.clone(),
                            invocation.system_label.clone(),
                            completed_at_unix_ms,
                        );
                    }
                    self.agent_last_invocation = Some(invocation);
                    self.agent_pending_image_attachment = None;
                    self.agent_help_capture_failure = None;
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
                Err(mpsc::TryRecvError::Empty) => {
                    task.runtime_frame.update_phase("waiting_for_models");
                }
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
                        if let Some(model) = self.selected_agent_discovered_model() {
                            self.agent_discovered_model_pick = model;
                        } else {
                            self.agent_discovered_model_pick.clear();
                        }
                        self.agent_model_discovery_status =
                            if self.agent_discovered_model_pick.trim().is_empty() {
                                Self::AGENT_MODEL_SELECTION_REQUIRED_MESSAGE.to_string()
                            } else {
                                format!(
                                    "Discovered {} model(s) in {:.1}s",
                                    self.agent_discovered_models.len(),
                                    elapsed
                                )
                            };
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
                    && let Some(estimate) = planning.get("estimate")
                {
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
                            && !errors.is_empty()
                        {
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
                        if let Some(warnings) =
                            preflight.get("warnings").and_then(|value| value.as_array())
                            && !warnings.is_empty()
                        {
                            ui.label("warnings");
                            for warning in warnings {
                                if let Some(text) = warning.as_str() {
                                    ui.small(format!("- {text}"));
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
            self.clear_viewport_foreground_request_after_render(viewport_id);
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
                crate::egui_compat::show_hosted_window(&mut *ctx, &spec, &mut open, |ui| {
                    close_requested = self.render_routine_assistant_contents(ui);
                });
                if close_requested {
                    open = false;
                }
            } else {
                let mut close_requested = false;
                crate::egui_compat::show_central_panel(
                    &mut *ctx,
                    egui::CentralPanel::default(),
                    |ui| {
                        close_requested = self.render_routine_assistant_contents(ui);
                    },
                );
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

    pub(super) fn render_agent_configuration_tab(&mut self, ui: &mut Ui) {
        self.refresh_agent_system_catalog();
        if !self.agent_token_file_credentials_loaded {
            self.refresh_agent_token_file_credentials();
        }
        ui.heading(self.tr("configuration.agent.heading"));
        ui.label(self.tr("configuration.agent.description"));
        ui.small(self.tr("configuration.agent.session_note"));
        ui.add_space(8.0);
        let mut preflight_inputs_changed = false;
        let mut requested_agent_system_id: Option<String> = None;
        let selected_system_text = self
            .agent_systems
            .iter()
            .find(|system| system.id == self.agent_system_id)
            .map(|system| format!("{} ({})", system.label, system.id))
            .unwrap_or_else(|| self.tr("agent.choose_system"));
        ui.horizontal(|ui| {
            ui.label(self.tr("agent.system"));
            egui::ComboBox::from_id_salt("agent_system_combo")
                .selected_text(selected_system_text)
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
            self.select_agent_system_and_persist_setup(&system_id);
        }
        if !self.agent_systems.is_empty() {
            egui::CollapsingHeader::new(self.tr("agent.quick_start.title"))
                .default_open(false)
                .show(ui, |ui| {
                    ui.small(self.tr("agent.quick_start.description"));
                    ui.horizontal_wrapped(|ui| {
                        if let Some(openai_system_id) =
                            preferred_openai_agent_system_id(&self.agent_systems)
                            && ui
                                .button(self.tr("agent.quick_start.openai"))
                                .on_hover_text(
                                    "Select the native OpenAI agent profile and use OPENAI_API_KEY for requests",
                                )
                                .clicked()
                            {
                                self.select_agent_system_and_persist_setup(&openai_system_id);
                                self.agent_base_url_override.clear();
                                self.agent_model_override.clear();
                                self.agent_discovered_model_pick.clear();
                                self.agent_status = self.tr("agent.status.selected_openai");
                            }
                        if let Some(anthropic_system_id) =
                            preferred_anthropic_agent_system_id(&self.agent_systems)
                            && ui
                                .button(self.tr("agent.quick_start.claude"))
                                .on_hover_text(
                                    "Select the native Anthropic Claude profile and use ANTHROPIC_API_KEY for requests",
                                )
                                .clicked()
                            {
                                self.select_agent_system_and_persist_setup(&anthropic_system_id);
                                self.agent_base_url_override.clear();
                                self.agent_model_override.clear();
                                self.agent_discovered_model_pick.clear();
                                self.agent_status = self.tr("agent.status.selected_claude");
                            }
                        if let Some(mistral_system_id) =
                            preferred_mistral_agent_system_id(&self.agent_systems)
                            && ui
                                .button(self.tr("agent.quick_start.mistral"))
                                .on_hover_text(
                                    "Select the native Mistral profile and use MISTRAL_API_KEY for requests",
                                )
                                .clicked()
                            {
                                self.select_agent_system_and_persist_setup(&mistral_system_id);
                                self.agent_base_url_override.clear();
                                self.agent_model_override.clear();
                                self.agent_discovered_model_pick.clear();
                                self.agent_status = self.tr("agent.status.selected_mistral");
                            }
                        if let Some(local_system_id) =
                            preferred_local_agent_system_id(&self.agent_systems)
                            && ui
                                .button(self.tr("agent.quick_start.local"))
                                .on_hover_text(
                                    "Select a local OpenAI-compatible endpoint such as Msty MLX, Msty gateway, Ollama, or Jan",
                                )
                                .clicked()
                            {
                                self.select_agent_system_and_persist_setup(&local_system_id);
                                self.agent_base_url_override.clear();
                                self.agent_model_override.clear();
                                self.agent_discovered_model_pick.clear();
                                self.agent_status = self.tr("agent.status.selected_local");
                            }
                        if self.agent_systems.iter().any(|system| system.id == "builtin_echo")
                            && ui
                                .button(self.tr("agent.quick_start.demo"))
                                .on_hover_text(
                                    "Select the offline demo assistant that never contacts a remote service",
                                )
                                .clicked()
                        {
                            self.select_agent_system_and_persist_setup("builtin_echo");
                            self.agent_status = self.tr("agent.status.selected_demo");
                        }
                    });
                    ui.small(self.tr("agent.quick_start.cloud_note"));
                });
            egui::CollapsingHeader::new(self.tr("agent.catalog"))
                .default_open(false)
                .show(ui, |ui| {
                    ui.horizontal(|ui| {
                        ui.text_edit_singleline(&mut self.agent_catalog_path);
                        if ui
                            .button(self.tr("button.browse"))
                            .on_hover_text("Browse filesystem and fill this path")
                            .clicked()
                            && let Some(path) = rfd::FileDialog::new()
                                .add_filter("JSON", &["json"])
                                .pick_file()
                        {
                            self.agent_catalog_path = path.display().to_string();
                            self.agent_catalog_loaded_path.clear();
                            self.refresh_agent_system_catalog();
                        }
                    });
                    if !self.agent_catalog_error.is_empty() {
                        ui.colored_label(
                            egui::Color32::from_rgb(190, 70, 70),
                            format!("Catalog error: {}", self.agent_catalog_error),
                        );
                    }
                });
            ui.group(|ui| {
                ui.strong(self.tr("agent.external_mcp.title"));
                ui.small(self.tr("agent.external_mcp.route"));
                ui.small(self.tr("agent.external_mcp.subscriptions"));
                let state_path = self.external_agent_mcp_state_path();
                ui.small(format!(
                    "{}: {}{}",
                    self.tr("agent.external_mcp.state_path"),
                    state_path,
                    if self.current_project_path.is_some() {
                        format!(" ({})", self.tr("agent.external_mcp.active_project"))
                    } else {
                        format!(" ({})", self.tr("agent.external_mcp.default_state_path"))
                    }
                ));
                let command = self.external_agent_mcp_command_snippet();
                if Self::render_copyable_command_line(ui, "", &command) {
                    self.agent_status = "Copied external MCP command".to_string();
                }
            });
        }
        if let Some(system) = self.selected_agent_system() {
            let (available, reason) = self.selected_agent_system_availability(&system);
            if let Some(description) = system.description.as_deref() {
                let trimmed = description.trim();
                if !trimmed.is_empty() {
                    ui.small(trimmed);
                }
            }
            ui.small(format!("transport: {}", system.transport.as_str()));
            if !system.command.is_empty() {
                let command = system.command.join(" ");
                if Self::render_copyable_command_line(ui, "command:", &command) {
                    self.agent_status = "Copied selected agent command".to_string();
                }
            }
            if agent_system_supports_model_selection(&system) {
                if let Some(source_key) = self.selected_agent_model_discovery_source_key(&system) {
                    if self.agent_model_discovery_source_key != source_key {
                        self.clear_agent_model_discovery_snapshot();
                        self.agent_model_discovery_source_key = source_key;
                        self.clear_agent_preflight_output();
                    }
                    if normalize_agent_model_name(self.agent_model_override.trim()).is_none()
                        && self.agent_model_discovery_task.is_none()
                        && self.agent_discovered_models.is_empty()
                    {
                        self.start_agent_model_discovery_task(&system, false);
                    }
                }
                if !matches!(system.transport, AgentSystemTransport::ExternalJsonStdio) {
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
                        AgentSystemTransport::NativeMistral => {
                            GUI_MISTRAL_DEFAULT_MODEL.to_string()
                        }
                        AgentSystemTransport::ExternalJsonStdio => {
                            if is_pi_local_agent_system(&system) {
                                "Pi default".to_string()
                            } else {
                                "Codex default".to_string()
                            }
                        }
                        _ => OPENAI_COMPAT_UNSPECIFIED_MODEL.to_string(),
                    });
                let model_override = normalize_agent_model_name(self.agent_model_override.trim());
                if let Some(model_override) = model_override {
                    ui.small(format!("model: {model_override} (session override)"));
                } else if let Some(discovered_model) = self.selected_agent_discovered_model() {
                    ui.small(format!(
                        "model: {} (selected discovered model)",
                        discovered_model
                    ));
                } else {
                    ui.small(format!("model: {catalog_model}"));
                }
            } else {
                self.clear_agent_model_discovery_snapshot();
            }
            if !available {
                if let Some(prompt) = self.agent_model_selection_prompt(&system) {
                    ui.small(prompt);
                } else {
                    ui.colored_label(
                        egui::Color32::from_rgb(190, 70, 70),
                        format!(
                            "Unavailable: {}",
                            reason.unwrap_or_else(|| "unknown reason".to_string())
                        ),
                    );
                }
            }
            if system.id == "builtin_echo" {
                ui.small(
                    "Built-in Echo is only a demo bridge. Use one of the quick-start buttons above for a real model-backed assistant.",
                );
            }
        } else if self.agent_systems.is_empty() {
            ui.small(self.tr("agent.no_systems_loaded"));
        }
        if let Some(selected_system) = self.selected_agent_system() {
            let credential_messages = self.selected_agent_credential_messages(&selected_system);
            let key_field_relevant = selected_system.id != "builtin_echo"
                && !is_codex_local_agent_system(&selected_system)
                && !is_pi_local_agent_system(&selected_system);
            if key_field_relevant {
                let (key_label, key_hint) = match selected_system.transport {
                    AgentSystemTransport::NativeAnthropic => {
                        (self.tr("agent.key.anthropic"), "sk-ant-...")
                    }
                    AgentSystemTransport::NativeMistral => {
                        (self.tr("agent.key.mistral"), "api key")
                    }
                    AgentSystemTransport::NativeOpenaiCompat => {
                        (self.tr("agent.key.openai_compatible"), "optional")
                    }
                    _ => (self.tr("agent.key.openai"), "sk-..."),
                };
                ui.horizontal_wrapped(|ui| {
                    ui.label(key_label);
                    let response = ui.add(
                        egui::TextEdit::singleline(&mut self.agent_openai_api_key)
                            .password(true)
                            .hint_text(key_hint),
                    );
                    preflight_inputs_changed |= response.changed();
                    if ui
                        .button(self.tr("agent.clear_key"))
                        .on_hover_text("Clear session-only API key override")
                        .clicked()
                    {
                        self.agent_openai_api_key.clear();
                        preflight_inputs_changed = true;
                    }
                    if agent_api_key_source(selected_system.transport).is_some()
                        && ui
                            .button(self.tr("agent.credential.reload"))
                            .on_hover_text("Read provider token files again")
                            .clicked()
                    {
                        self.refresh_agent_token_file_credentials();
                        preflight_inputs_changed = true;
                    }
                });
            }
            for (message, warning) in credential_messages {
                if warning {
                    ui.colored_label(egui::Color32::from_rgb(180, 100, 40), message);
                } else {
                    ui.small(message);
                }
            }
        }
        let base_url_placeholder = self.selected_agent_base_url_placeholder();
        ui.horizontal(|ui| {
            ui.label(self.tr("agent.base_url_override"));
            let response = ui.add(
                egui::TextEdit::singleline(&mut self.agent_base_url_override)
                    .hint_text(base_url_placeholder),
            );
            preflight_inputs_changed |= response.changed();
            if ui
                .button(self.tr("agent.clear_url"))
                .on_hover_text("Clear session-only base URL override")
                .clicked()
            {
                self.agent_base_url_override.clear();
                preflight_inputs_changed = true;
            }
        });
        ui.horizontal(|ui| {
            ui.label(self.tr("agent.model_override"));
            let response = ui.add(
                egui::TextEdit::singleline(&mut self.agent_model_override).hint_text("unspecified"),
            );
            preflight_inputs_changed |= response.changed();
            if ui
                .button(self.tr("agent.clear_model"))
                .on_hover_text("Clear session-only model override")
                .clicked()
            {
                self.agent_model_override.clear();
                self.agent_discovered_model_pick.clear();
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
                .button(self.tr("agent.clear_timeout"))
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
                .button(self.tr("agent.clear_http_timeouts"))
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
                .button(self.tr("agent.clear_limits"))
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
                    .button(self.tr("agent.test_setup"))
                    .on_hover_text(
                        "Validate the current system, key, endpoint, model, and runtime settings without sending a prompt",
                    )
                    .clicked()
                {
                    self.run_agent_preflight_probe();
                }
                if self.agent_preflight_output.is_some()
                    && ui
                        .button(self.tr("agent.clear_test"))
                        .on_hover_text("Clear the latest setup-preflight snapshot")
                        .clicked()
                {
                    self.clear_agent_preflight_output();
                }
                if agent_system_supports_model_discovery(&system) {
                    let discovery_hover = if matches!(
                        system.transport,
                        AgentSystemTransport::ExternalJsonStdio
                    ) {
                        if is_pi_local_agent_system(&system) {
                            "Read selectable provider/model ids from the installed Pi CLI"
                        } else {
                            "Read selectable model ids from the local Codex model metadata cache"
                        }
                    } else {
                        "Query local/server model list from current base URL"
                    };
                    if ui
                        .button(self.tr("agent.discover_models"))
                        .on_hover_text(discovery_hover)
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
            if !agent_system_supports_model_discovery(&system) {
                self.clear_agent_model_discovery_snapshot();
            } else {
                if self.should_render_agent_model_selector(&system) {
                    let codex_local = is_codex_local_agent_system(&system);
                    let pi_local = is_pi_local_agent_system(&system);
                    let local_default = if pi_local {
                        "Pi default"
                    } else {
                        "Codex default"
                    };
                    let previous_pick = self.agent_discovered_model_pick.clone();
                    ui.horizontal(|ui| {
                        ui.label(self.tr("agent.discovered_model"));
                        egui::ComboBox::from_id_salt("agent_discovered_model_combo")
                            .selected_text(if self.agent_discovered_model_pick.trim().is_empty() {
                                if codex_local || pi_local {
                                    local_default.to_string()
                                } else {
                                    self.tr("agent.choose_model")
                                }
                            } else {
                                self.agent_discovered_model_pick.clone()
                            })
                            .show_ui(ui, |ui| {
                                if codex_local || pi_local {
                                    ui.selectable_value(
                                        &mut self.agent_discovered_model_pick,
                                        String::new(),
                                        local_default,
                                    );
                                    ui.separator();
                                }
                                for model in &self.agent_discovered_models {
                                    ui.selectable_value(
                                        &mut self.agent_discovered_model_pick,
                                        model.clone(),
                                        model,
                                    );
                                }
                            });
                    });
                    if previous_pick != self.agent_discovered_model_pick {
                        self.agent_model_override.clear();
                        preflight_inputs_changed = true;
                    }
                    if codex_local {
                        ui.small(
                            "Models are read from the logged-in Codex CLI metadata cache. Codex default leaves model choice to the CLI.",
                        );
                    } else if pi_local {
                        ui.small(
                            "Models are read from 'pi --list-models'. Pi default leaves model choice to Pi; authenticate providers in Pi with /login.",
                        );
                    } else {
                        ui.small(self.tr("agent.discovered_model_note"));
                    }
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
                ui.strong(self.tr("agent.setup_preflight"));
                let (overall_label, color) = Self::agent_preflight_overall_label(preflight);
                ui.colored_label(
                    color,
                    format!(
                        "{} (static={}, transport={})",
                        overall_label,
                        if preflight.available { "available" } else { "unavailable" },
                        preflight.transport,
                    ),
                );
                if let Some(reason) = preflight.availability_reason.as_deref()
                    && !reason.trim().is_empty() {
                        ui.small(format!("detail: {}", reason.trim()));
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
                    ui.strong(self.tr("agent.live_probe"));
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
                    if Self::agent_status_mentions_openai_quota(&live.message)
                        || (live.status_class == AgentLiveProbeStatusClass::QuotaOrBilling
                            && preflight.transport == AgentSystemTransport::NativeOpenai.as_str())
                    {
                        self.render_openai_quota_links(ui);
                    }
                    match live.probe_kind {
                        AgentLiveProbeKind::CommandShape => ui.small(format!(
                            "probe=command_shape | executable_reachable={} | command_shape_ok={}",
                            live.reachable,
                            live.status_class == AgentLiveProbeStatusClass::Ok
                        )),
                        AgentLiveProbeKind::ModelDiscovery => ui.small(format!(
                            "probe=model_discovery | reachable={} | auth_ok={} | model_list_ok={} | selected_model_seen={}",
                            live.reachable,
                            live.auth_ok,
                            live.model_list_ok,
                            live.selected_model_seen
                        )),
                    };
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
                    ui.strong(self.tr("agent.next_action"));
                    for action in next_actions {
                        ui.small(action);
                    }
                }
            });
        }
        ui.small(
            "Session only: if set, this key overrides the selected provider API key env var for agent requests started from this GUI window.",
        );
        ui.small(
            "Session only: Base URL override applies to native_openai/native_anthropic/native_mistral/native_openai_compat. For local roots (e.g. http://localhost:11964), GENtle tries /chat/completions and /v1/chat/completions on that same base URL; Msty MLX model servers commonly use http://localhost:11973/v1.",
        );
        ui.small(
            "Session only: Model override applies to native providers and Codex Local, and maps to GENTLE_AGENT_MODEL. The Codex bridge forwards it as codex --model. Value 'unspecified' means no override.",
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
    }

    pub(super) fn agent_sequence_object_commands(seq_id: &str) -> Vec<AgentObjectCommand> {
        let seq_id = Self::shell_quote_command_arg(seq_id);
        vec![
            AgentObjectCommand {
                label: "Open sequence",
                detail: "Open this project sequence in its DNA/RNA window.",
                command: format!("/open sequence-window {seq_id}"),
            },
            AgentObjectCommand {
                label: "List annotations",
                detail: "Return a read-only structured list of annotations on this sequence.",
                command: format!("features query {seq_id} --limit 100"),
            },
            AgentObjectCommand {
                label: "Find restriction sites",
                detail: "Run the shared read-only restriction-site scan on this sequence.",
                command: format!("features restriction-scan {seq_id}"),
            },
        ]
    }

    pub(super) fn agent_container_object_command(
        container_id: &str,
        declared_contents_exclusive: bool,
    ) -> AgentObjectCommand {
        let container_id = Self::shell_quote_command_arg(container_id);
        let (label, detail, exclusive) = if declared_contents_exclusive {
            (
                "Mark contents as a measured subset",
                "Match the subset interpretation offered by the main project overview.",
                false,
            )
        } else {
            (
                "Mark contents as exhaustive",
                "Match the exhaustive-content interpretation offered by the main project overview.",
                true,
            )
        };
        AgentObjectCommand {
            label,
            detail,
            command: format!("containers set-exclusive {container_id} {exclusive}"),
        }
    }

    fn render_agent_object_command_menu(
        ui: &mut egui::Ui,
        object_label: &str,
        commands: &[AgentObjectCommand],
        pending_command: &mut Option<String>,
    ) {
        ui.set_min_width(440.0);
        ui.strong(format!("Actions for {object_label}"));
        ui.small("These are the same parser-validated commands used by GUI Shell and CLI Shell.");
        ui.separator();
        for action in commands {
            ui.label(egui::RichText::new(action.label).strong());
            ui.small(action.detail);
            ui.monospace(&action.command);
            ui.horizontal(|ui| {
                if ui.button("Run").clicked() {
                    *pending_command = Some(action.command.clone());
                    ui.close();
                }
                if ui.button("Copy command").clicked() {
                    ui.ctx().copy_text(action.command.clone());
                }
            });
            ui.separator();
        }
        ui.small("Use /help for the complete command catalog.");
    }

    fn render_agent_project_state_summary(
        &mut self,
        ui: &mut egui::Ui,
        summary: &EngineStateSummary,
    ) -> Option<String> {
        let mut pending_command = None;
        ui.label(format!(
            "Current project: {} sequence(s), {} container(s), {} arrangement(s).",
            summary.sequence_count, summary.container_count, summary.arrangement_count
        ));
        ui.horizontal_wrapped(|ui| {
            if ui
                .button("Show complete project overview")
                .on_hover_text(
                    "Focus the main window's table/graph with sequences, analyses, containers, arrangements, and their context menus",
                )
                .clicked()
            {
                self.focus_project_overview_target(ProjectOverviewTarget::Lineage);
                self.queue_focus_viewport(egui::ViewportId::ROOT);
                self.agent_status = "Focused the complete project overview in the main window"
                    .to_string();
            }
            ui.small("Right-click a sequence or container below for its contextual commands.");
        });
        if summary.sequence_count == 0
            && summary.container_count == 0
            && summary.arrangement_count == 0
        {
            ui.small("No sequences, containers, or arrangements are currently loaded.");
            return pending_command;
        }

        if !summary.sequences.is_empty() {
            egui::CollapsingHeader::new(format!("Sequences ({})", summary.sequence_count))
                .default_open(true)
                .show(ui, |ui| {
                    for sequence in summary.sequences.iter().take(50) {
                        let name = sequence
                            .name
                            .as_deref()
                            .map(str::trim)
                            .filter(|name| !name.is_empty() && *name != sequence.id.as_str())
                            .map(|name| format!(" | {name}"))
                            .unwrap_or_default();
                        let response = ui.add(
                            egui::Label::new(format!(
                                "{}{} | {} bp | {}",
                                sequence.id,
                                name,
                                sequence.length,
                                if sequence.circular {
                                    "circular"
                                } else {
                                    "linear"
                                }
                            ))
                            .wrap(),
                        );
                        response
                            .on_hover_text("Right-click for commands that operate on this sequence")
                            .context_menu(|ui| {
                                Self::render_agent_object_command_menu(
                                    ui,
                                    &sequence.id,
                                    &Self::agent_sequence_object_commands(&sequence.id),
                                    &mut pending_command,
                                );
                            });
                    }
                    if summary.sequences.len() > 50 {
                        ui.small(format!(
                            "... and {} more sequence(s); use Copy JSON for the complete list.",
                            summary.sequences.len() - 50
                        ));
                    }
                });
        }
        if !summary.containers.is_empty() {
            egui::CollapsingHeader::new(format!("Containers ({})", summary.container_count))
                .default_open(false)
                .show(ui, |ui| {
                    for container in summary.containers.iter().take(50) {
                        let response = ui.add(
                            egui::Label::new(format!(
                                "{} | {} | {} member(s)",
                                container.id, container.kind, container.member_count
                            ))
                            .wrap(),
                        );
                        response
                            .on_hover_text(
                                "Right-click for commands that operate on this container",
                            )
                            .context_menu(|ui| {
                                Self::render_agent_object_command_menu(
                                    ui,
                                    &container.id,
                                    &[Self::agent_container_object_command(
                                        &container.id,
                                        container.declared_contents_exclusive,
                                    )],
                                    &mut pending_command,
                                );
                            });
                    }
                    if summary.containers.len() > 50 {
                        ui.small(format!(
                            "... and {} more container(s); use Copy JSON for the complete list.",
                            summary.containers.len() - 50
                        ));
                    }
                });
        }
        if !summary.arrangements.is_empty() {
            egui::CollapsingHeader::new(format!("Arrangements ({})", summary.arrangement_count))
                .default_open(false)
                .show(ui, |ui| {
                    for arrangement in summary.arrangements.iter().take(50) {
                        ui.add(
                            egui::Label::new(format!(
                                "{} | {} | {} lane(s)",
                                arrangement.id, arrangement.mode, arrangement.lane_count
                            ))
                            .wrap(),
                        )
                        .on_hover_text(
                            "Use Show complete project overview for arrangement actions",
                        );
                    }
                    if summary.arrangements.len() > 50 {
                        ui.small(format!(
                            "... and {} more arrangement(s); use Copy JSON for the complete list.",
                            summary.arrangements.len() - 50
                        ));
                    }
                });
        }
        pending_command
    }

    fn render_agent_command_output(&mut self, ui: &mut egui::Ui, result: &AgentCommandOutput) {
        let pretty_output = serde_json::to_string_pretty(&result.output)
            .unwrap_or_else(|_| result.output.to_string());
        ui.separator();
        let mut pending_command = None;
        ui.group(|ui| {
            ui.horizontal_wrapped(|ui| {
                ui.strong("Command result");
                ui.monospace(result.command.trim());
                ui.small(if result.state_changed {
                    "project changed"
                } else {
                    "read only"
                });
                if ui
                    .button("Copy JSON")
                    .on_hover_text("Copy the complete structured result from this local command")
                    .clicked()
                {
                    ui.ctx().copy_text(pretty_output.clone());
                    self.agent_status =
                        format!("Copied structured result for {}", result.command.trim());
                }
            });

            if let Ok(summary) = serde_json::from_value::<EngineStateSummary>(result.output.clone())
            {
                pending_command = self.render_agent_project_state_summary(ui, &summary);
            } else if result.output.get("help").is_some()
                || result.output.get("help_markdown").is_some()
            {
                ui.label("The shared-shell command reference is open in Help > Shell Commands.");
                ui.small("Use the Help window's Find field to locate a command or topic.");
            } else {
                ui.small("Structured command output:");
                let mut visible_output = pretty_output;
                egui::ScrollArea::vertical()
                    .id_salt("agent_local_command_output_scroll")
                    .max_height(220.0)
                    .auto_shrink([false, true])
                    .show(ui, |ui| {
                        ui.add(
                            egui::TextEdit::multiline(&mut visible_output)
                                .code_editor()
                                .interactive(false)
                                .desired_rows(8)
                                .desired_width(ui.available_width()),
                        );
                    });
            }
        });
        if let Some(command) = pending_command {
            self.execute_agent_shell_command_from_ui(
                0,
                "Project item action",
                &command,
                "list context menu",
            );
        }
    }

    pub(super) fn render_agent_assistant_contents(&mut self, ui: &mut Ui) -> bool {
        #[cfg(feature = "gui-test-support")]
        crate::gui_test_support::register_rect(
            ui.ctx().clone(),
            "window.agent_assistant",
            "window.agent_assistant",
            None,
            crate::gui_test_support::GuiTestWidgetKind::Window,
            ui.max_rect(),
            true,
            true,
            true,
            Some(if self.agent_task.is_some() {
                "running"
            } else {
                "ready"
            }),
        );
        self.refresh_agent_system_catalog();
        let mut close_requested = false;
        let close_hover = Self::specialist_window_close_hover_text("Agent Assistant");
        let close_label = self.tr("button.close");
        if self.render_specialist_window_nav_with_close(
            ui,
            Some((close_label.as_str(), close_hover.as_str())),
        ) {
            close_requested = true;
        }
        ui.label(self.tr("agent.description"));
        let selected_system = self.selected_agent_system();
        let selected_available = selected_system
            .as_ref()
            .map(|system| self.selected_agent_system_availability(system).0)
            .unwrap_or(false);
        ui.group(|ui| {
            ui.horizontal_wrapped(|ui| {
                if let Some(system) = selected_system.as_ref() {
                    ui.strong(format!("{} ({})", system.label, system.id));
                    let model = normalize_agent_model_name(self.agent_model_override.trim())
                        .or_else(|| self.selected_agent_discovered_model())
                        .or_else(|| system.model.as_deref().and_then(normalize_agent_model_name));
                    if let Some(model) = model {
                        ui.small(format!("model: {model}"));
                    }
                } else {
                    ui.colored_label(
                        egui::Color32::from_rgb(190, 70, 70),
                        self.tr("agent.choose_system"),
                    );
                }
                if ui
                    .button(self.tr("agent.configure"))
                    .on_hover_text(self.tr("agent.configure.tooltip"))
                    .clicked()
                {
                    self.open_configuration_agent_systems_dialog();
                }
            });
            ui.small(self.tr("agent.configuration_separated_note"));
        });
        let include_state_summary_label = self.tr("agent.include_state_summary");
        let auto_run_suggestions_label = self.tr("agent.auto_run_suggestions");
        let allow_web_research_label = self.tr("agent.allow_web_research");
        let allow_web_research_tooltip = self.tr("agent.allow_web_research.tooltip");
        ui.horizontal_wrapped(|ui| {
            ui.checkbox(
                &mut self.agent_include_state_summary,
                include_state_summary_label,
            )
            .on_hover_text(self.tr("agent.include_state_summary.tooltip"));
            ui.checkbox(&mut self.agent_allow_auto_exec, auto_run_suggestions_label)
                .on_hover_text(self.tr("agent.auto_run_suggestions.tooltip"));
            if selected_system
                .as_ref()
                .is_some_and(|system| system.supports_web_research)
            {
                ui.checkbox(&mut self.agent_allow_web_research, allow_web_research_label)
                    .on_hover_text(allow_web_research_tooltip);
            }
        });
        self.render_agent_help_attachment_panel(ui);
        self.render_agent_screenshot_consent_card(ui);
        if self.agent_initial_actions_visible() {
            let mut open_document = false;
            let mut open_configuration = false;
            ui.group(|ui| {
                ui.strong(self.tr("agent.initial_actions.title"));
                ui.small(self.tr("agent.initial_actions.description"));
                ui.horizontal_wrapped(|ui| {
                    if ui
                        .button(self.tr("agent.initial_actions.open_previous"))
                        .on_hover_text(self.tr("agent.initial_actions.open_previous.tooltip"))
                        .clicked()
                    {
                        open_document = true;
                    }
                    if ui
                        .button(self.tr("agent.initial_actions.configure"))
                        .on_hover_text(self.tr("agent.initial_actions.configure.tooltip"))
                        .clicked()
                    {
                        open_configuration = true;
                    }
                });
            });
            if open_document {
                self.prompt_open_sequence();
            }
            if open_configuration {
                self.open_configuration_dialog();
            }
        }
        if !self.agent_conversation.turns.is_empty() {
            let turns = self.agent_conversation.turns.clone();
            let mut copied_response = false;
            egui::CollapsingHeader::new(format!(
                "{} ({})",
                self.tr("agent.conversation"),
                turns.len()
            ))
            .default_open(true)
            .show(ui, |ui| {
                ui.small(self.tr("agent.conversation.project_storage_note"));
                egui::ScrollArea::vertical()
                    .id_salt("agent_conversation_scroll")
                    .max_height(320.0)
                    .auto_shrink([false, true])
                    .show(ui, |ui| {
                        for (index, turn) in turns.iter().enumerate() {
                            if index > 0 {
                                ui.separator();
                            }
                            ui.strong(self.tr("agent.conversation.you"));
                            ui.add(egui::Label::new(turn.user_message.trim()).wrap());
                            for attachment in &turn.attachments {
                                let dimensions =
                                    match (attachment.pixel_width, attachment.pixel_height) {
                                        (Some(width), Some(height)) => {
                                            format!(" | {width} x {height} px")
                                        }
                                        _ => String::new(),
                                    };
                                ui.small(format!(
                                    "Attached image: {}{} | source {}",
                                    attachment.file_name,
                                    dimensions,
                                    attachment
                                        .source_window_title
                                        .as_deref()
                                        .unwrap_or("GENtle window")
                                ));
                            }
                            ui.horizontal_wrapped(|ui| {
                                ui.strong(if turn.system_label.trim().is_empty() {
                                    turn.system_id.as_str()
                                } else {
                                    turn.system_label.as_str()
                                });
                                if ui
                                    .small_button(self.tr("agent.conversation.copy"))
                                    .on_hover_text(
                                        "Copy this stored agent response as structured JSON",
                                    )
                                    .clicked()
                                    && let Ok(payload) =
                                        serde_json::to_string_pretty(&turn.response)
                                {
                                    ui.ctx().copy_text(payload);
                                    copied_response = true;
                                }
                            });
                            ui.add(egui::Label::new(turn.response.assistant_message.trim()).wrap());
                            self.render_agent_web_research(ui, &turn.response);
                            for question in &turn.response.questions {
                                ui.add(
                                    egui::Label::new(
                                        egui::RichText::new(format!(
                                            "Question: {}",
                                            question.trim()
                                        ))
                                        .small(),
                                    )
                                    .wrap(),
                                );
                            }
                            if let Some(request) = &turn.response.screenshot_request {
                                ui.small(format!(
                                    "Screenshot request {}: {} (one-shot approval is not stored)",
                                    request.id, request.reason
                                ));
                            }
                        }
                    });
            });
            if copied_response {
                self.agent_status = "Copied stored agent response JSON".to_string();
            }
        }
        if !agent_prompt_template_options()
            .iter()
            .any(|(id, _)| *id == self.agent_prompt_template_id)
        {
            self.agent_prompt_template_id = AGENT_PROMPT_TEMPLATE_DEFAULT_ID.to_string();
        }
        ui.horizontal(|ui| {
            ui.label(self.tr("agent.prompt_template"))
                .on_hover_text(self.tr("agent.prompt_template.tooltip"));
            let template_response = egui::ComboBox::from_id_salt("agent_prompt_template_combo")
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
            template_response
                .response
                .on_hover_text(self.tr("agent.prompt_template.tooltip"));
            if ui
                .button(self.tr("agent.insert"))
                .on_hover_text(
                    "Replace current prompt with selected template and apply its request defaults",
                )
                .clicked()
            {
                self.agent_prompt =
                    agent_prompt_template_text(&self.agent_prompt_template_id).to_string();
                self.agent_include_state_summary =
                    agent_prompt_template_includes_state_summary_by_default(
                        &self.agent_prompt_template_id,
                    );
            }
            if ui
                .button(self.tr("agent.append"))
                .on_hover_text(
                    "Append selected template below current prompt and apply its request defaults",
                )
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
                self.agent_include_state_summary =
                    agent_prompt_template_includes_state_summary_by_default(
                        &self.agent_prompt_template_id,
                    );
            }
        });
        ui.horizontal_wrapped(|ui| {
            ui.label(self.tr("agent.prompt"));
            if ui
                .button(self.tr("agent.attach_document"))
                .on_hover_text(self.tr("agent.attach_document.tooltip"))
                .clicked()
                && let Some(path) = rfd::FileDialog::new()
                    .add_filter(
                        "Text documents",
                        &[
                            "md", "markdown", "txt", "rst", "log", "json", "toml", "yaml", "yml",
                            "csv", "tsv",
                        ],
                    )
                    .pick_file()
            {
                let path = fs::canonicalize(&path).unwrap_or(path);
                if !self.agent_prompt.trim().is_empty() {
                    self.agent_prompt.push('\n');
                }
                self.agent_prompt.push_str(&format!(
                    "Use local text document `{}` as reference context.",
                    path.display()
                ));
            }
        });
        let local_document_paths = agent_explicit_local_document_paths(&self.agent_prompt);
        if !local_document_paths.is_empty() {
            ui.colored_label(
                egui::Color32::from_rgb(165, 105, 35),
                self.tr("agent.local_documents.notice"),
            );
            ui.horizontal_wrapped(|ui| {
                ui.small(self.tr("agent.local_documents.detected"));
                for path in &local_document_paths {
                    ui.monospace(path.display().to_string());
                }
            });
        }
        let prompt_edit_id = ui.make_persistent_id("agent_assistant_prompt_edit");
        let prompt_submit_shortcut = ui.memory(|memory| memory.has_focus(prompt_edit_id))
            && ui.input_mut(|input| {
                input.consume_shortcut(&egui::KeyboardShortcut::new(
                    egui::Modifiers::COMMAND,
                    egui::Key::Enter,
                )) || input.consume_shortcut(&egui::KeyboardShortcut::new(
                    egui::Modifiers::CTRL,
                    egui::Key::Enter,
                ))
            });
        ui.add(
            egui::TextEdit::multiline(&mut self.agent_prompt)
                .id(prompt_edit_id)
                .desired_rows(6)
                .desired_width(f32::INFINITY),
        );
        let mut running = self.agent_task.is_some();
        if running && Self::consume_command_or_ctrl_shortcut(ui.ctx(), egui::Key::Period) {
            self.request_agent_task_cancel("agent assistant shortcut");
            running = self.agent_task.is_some();
        }
        let direct_prompt_command =
            Self::agent_prompt_direct_shell_command(&self.agent_prompt).map(str::to_string);
        let selected_supports_pending_attachment = self.agent_pending_image_attachment.is_none()
            || selected_system
                .as_ref()
                .map(|system| system.supports_image_attachments)
                .unwrap_or(false);
        let can_submit_prompt = !running
            && self.agent_screenshot_capture.is_none()
            && (selected_available || direct_prompt_command.is_some())
            && (direct_prompt_command.is_some() || selected_supports_pending_attachment);
        if prompt_submit_shortcut && can_submit_prompt {
            if let Some(command) = direct_prompt_command.as_deref() {
                self.execute_agent_prompt_command(command);
            } else {
                self.start_agent_assistant_request();
            }
        }
        ui.horizontal(|ui| {
            let ask_button_text = if direct_prompt_command.is_some() {
                "Run command".to_string()
            } else {
                self.tr("agent.ask_agent")
            };
            let ask_hover_text = if direct_prompt_command.is_some() {
                "Run this GENtle Agent Assistant slash command locally (Command/Ctrl+Return)"
            } else {
                "Send prompt to selected agent system (Command/Ctrl+Return)"
            };
            let ask_response = ui
                .add_enabled(can_submit_prompt, egui::Button::new(ask_button_text))
                .on_hover_text(ask_hover_text);
            #[cfg(feature = "gui-test-support")]
            crate::gui_test_support::register_response(
                &ask_response,
                "agent.ask",
                "window.agent_assistant",
                None,
                crate::gui_test_support::GuiTestWidgetKind::Button,
                false,
            );
            if ask_response.clicked() {
                if let Some(command) = direct_prompt_command.as_deref() {
                    self.execute_agent_prompt_command(command);
                } else {
                    self.start_agent_assistant_request();
                }
            }
            if ui
                .add_enabled(
                    !running,
                    egui::Button::new(self.tr("agent.clear_conversation")),
                )
                .on_hover_text(
                    "Clear the project-stored conversation, latest response, and agent status",
                )
                .clicked()
            {
                self.clear_agent_conversation();
            }
            if ui
                .button(self.tr("agent.clear_execution_log"))
                .on_hover_text("Clear local execution history for agent suggestions")
                .clicked()
            {
                self.agent_execution_log.clear();
                self.agent_last_command_output = None;
            }
        });
        let mut stop_agent_request = false;
        if let Some(task) = &self.agent_task {
            ui.horizontal(|ui| {
                ui.add(egui::Spinner::new());
                ui.label(format!(
                    "Agent request running ({:.1}s)",
                    task.started.elapsed().as_secs_f32()
                ));
                if ui
                    .button("Stop")
                    .on_hover_text(
                        "Stop waiting for the current agent request (Command/Ctrl+Period)",
                    )
                    .clicked()
                {
                    stop_agent_request = true;
                }
            });
        }
        if stop_agent_request {
            self.request_agent_task_cancel("agent assistant");
        }
        if self.sequence_window_open_in_progress() {
            ui.horizontal_wrapped(|ui| {
                ui.label(
                    egui::RichText::new(self.tr("sequence.window.opening"))
                        .strong()
                        .color(egui::Color32::from_rgb(40, 120, 170)),
                );
                ui.small(self.tr("agent.sequence_window.opening.detail"));
            });
        }
        if !self.agent_status.is_empty() {
            ui.separator();
            self.render_agent_status_message(ui, &self.agent_status, true);
        }

        if let Some(output) = self.agent_last_command_output.clone() {
            self.render_agent_command_output(ui, &output);
        }

        if let Some(invocation) = self.agent_last_invocation.clone() {
            ui.separator();
            ui.horizontal_wrapped(|ui| {
                ui.label(format!(
                    "Latest response from {} ({})",
                    invocation.system_label, invocation.system_id
                ));
                if ui
                    .button("Copy Response JSON")
                    .on_hover_text("Copy the latest agent response JSON to the clipboard")
                    .clicked()
                {
                    let payload = Self::agent_response_clipboard_payload(&invocation);
                    ui.ctx().copy_text(payload);
                    self.agent_status = "Copied latest agent response JSON".to_string();
                }
            });
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
            let sanity_warnings = Self::agent_response_sanity_warnings(&invocation);
            if !sanity_warnings.is_empty() {
                ui.group(|ui| {
                    ui.strong("Response sanity checks");
                    for warning in &sanity_warnings {
                        let warning = Self::compact_agent_validation_message(warning);
                        ui.add(
                            egui::Label::new(
                                egui::RichText::new(&warning)
                                    .color(egui::Color32::from_rgb(180, 110, 25)),
                            )
                            .wrap(),
                        );
                    }
                    ui.add(
                        egui::Label::new(
                            egui::RichText::new(
                                "These checks are deterministic GENtle validation hints, not a second AI judgment.",
                            )
                            .small(),
                        )
                        .wrap(),
                    );
                });
            }
            if !invocation.response.assistant_message.trim().is_empty() {
                ui.group(|ui| {
                    ui.strong("Agent message");
                    ui.add(egui::Label::new(invocation.response.assistant_message.trim()).wrap());
                });
            }
            self.render_agent_web_research(ui, &invocation.response);
            if !invocation.response.questions.is_empty() {
                ui.group(|ui| {
                    ui.strong("Agent questions");
                    for question in &invocation.response.questions {
                        ui.add(egui::Label::new(format!("- {}", question)).wrap());
                    }
                });
            }
            if invocation.response.suggested_commands.is_empty() {
                ui.small("No executable suggestions in this reply.");
            } else {
                ui.separator();
                ui.strong("Suggested commands");
                let mut run_request: Option<(usize, AgentSuggestedCommand)> = None;
                for (idx, suggestion) in invocation.response.suggested_commands.iter().enumerate() {
                    let index_1based = idx + 1;
                    let command_blocker = Self::agent_suggestion_run_blocker(
                        &suggestion.command,
                        suggestion.execution,
                    );
                    let precondition_blocker = self.agent_suggestion_precondition_blocker(
                        suggestion.precondition_expr.as_ref(),
                    );
                    let run_blocker = command_blocker
                        .clone()
                        .or_else(|| precondition_blocker.clone());
                    ui.group(|ui| {
                        if precondition_blocker.is_some() {
                            ui.disable();
                        }
                        ui.horizontal_wrapped(|ui| {
                            ui.strong(format!("#{index_1based}"));
                            let run_response = ui
                                .add_enabled(run_blocker.is_none(), egui::Button::new("Run"))
                                .on_hover_text(run_blocker.as_deref().unwrap_or(
                                    "Execute this suggestion through GENtle's shared shell parser/executor",
                                ));
                            if run_response.clicked() {
                                run_request = Some((index_1based, suggestion.clone()));
                            }
                            ui.strong(
                                suggestion
                                    .title
                                    .as_deref()
                                    .unwrap_or("Suggested GENtle command"),
                            );
                            ui.small(format!("mode: {}", suggestion.execution.as_str()));
                        });
                        if let Some(reason) = &run_blocker {
                            let reason_color = if precondition_blocker.is_some()
                                || suggestion.execution == AgentExecutionIntent::Chat
                            {
                                ui.visuals().weak_text_color()
                            } else {
                                egui::Color32::from_rgb(190, 70, 70)
                            };
                            ui.add(
                                egui::Label::new(
                                    egui::RichText::new(reason).color(reason_color),
                                )
                                .wrap(),
                            );
                        }
                        let command_text = egui::RichText::new(suggestion.command.trim())
                            .monospace()
                            .color(if command_blocker.is_some()
                                && suggestion.execution != AgentExecutionIntent::Chat
                            {
                                egui::Color32::from_rgb(190, 70, 70)
                            } else {
                                ui.visuals().text_color()
                            });
                        ui.add(egui::Label::new(command_text).wrap());

                        let mut details = Vec::new();
                        if !suggestion.preconditions.is_empty() {
                            details.push(format!(
                                "Preconditions: {}",
                                suggestion.preconditions.join("; ")
                            ));
                        }
                        if let Some(expr) = &suggestion.precondition_expr {
                            if let Some(readiness) = self.agent_suggestion_fact_readiness(expr) {
                                details.push(format!("Readiness: {readiness}"));
                            }
                            if let Ok(expr_json) = serde_json::to_string(expr) {
                                details.push(format!("Precondition logic: {expr_json}"));
                            }
                        }
                        if !suggestion.expected_outcomes.is_empty() {
                            details.push(format!(
                                "Expected outcomes: {}",
                                suggestion.expected_outcomes.join("; ")
                            ));
                        }
                        if !suggestion.expected_effects.is_empty()
                            && let Ok(effects_json) =
                                serde_json::to_string(&suggestion.expected_effects)
                        {
                            details.push(format!("Expected effects: {effects_json}"));
                        }
                        if let Some(rationale) = suggestion.rationale.as_deref()
                            && !rationale.trim().is_empty()
                        {
                            details.push(format!("Rationale: {}", rationale.trim()));
                        }
                        for detail in details {
                            ui.add(
                                egui::Label::new(egui::RichText::new(detail).small()).wrap(),
                            );
                        }
                    });
                }
                if let Some((index_1based, suggestion)) = run_request {
                    self.execute_agent_suggestion(index_1based, &suggestion, "manual");
                }
            }
            if !invocation.raw_stderr.trim().is_empty() {
                ui.separator();
                ui.strong("Agent stderr");
                let mut stderr = invocation.raw_stderr.clone();
                egui::ScrollArea::horizontal()
                    .id_salt("agent_stderr_horizontal_scroll")
                    .max_height(120.0)
                    .show(ui, |ui| {
                        ui.add(
                            egui::TextEdit::multiline(&mut stderr)
                                .desired_rows(4)
                                .desired_width(1200.0),
                        );
                    });
            }
        }

        if !self.agent_execution_log.is_empty() {
            ui.separator();
            ui.strong("Execution log");
            egui::ScrollArea::vertical()
                .id_salt("agent_execution_log_scroll")
                .max_height(180.0)
                .auto_shrink([false, true])
                .show(ui, |ui| {
                    scroll_input_policy::apply_scrollarea_keyboard_navigation(
                        ui,
                        scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                    );
                    for entry in self.agent_execution_log.iter().rev() {
                        let source = if entry.index_1based == 0 {
                            "prompt".to_string()
                        } else {
                            format!("#{}", entry.index_1based)
                        };
                        ui.add(
                            egui::Label::new(format!(
                                "{} [{}] {} | {} | changed={} | t={}",
                                source,
                                entry.trigger,
                                if entry.ok { "ok" } else { "error" },
                                entry.command,
                                entry.state_changed,
                                entry.executed_at_unix_ms
                            ))
                            .wrap(),
                        );
                        ui.add(
                            egui::Label::new(egui::RichText::new(&entry.summary).small()).wrap(),
                        );
                    }
                });
        }
        close_requested
    }

    pub(super) fn render_agent_assistant_contents_scrollable(
        &mut self,
        ui: &mut Ui,
        id_salt: &'static str,
    ) -> egui::containers::scroll_area::ScrollAreaOutput<bool> {
        window_backdrop::paint_window_backdrop(
            ui,
            WindowBackdropKind::AgentAssistant,
            &self.window_backdrops,
        );
        with_window_content_inset(ui, |ui| {
            egui::ScrollArea::vertical()
                .id_salt(id_salt)
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    scroll_input_policy::apply_scrollarea_keyboard_navigation(
                        ui,
                        scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                    );
                    self.render_agent_assistant_contents(ui)
                })
        })
    }

    pub(super) fn render_agent_assistant_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_agent_assistant_dialog {
            return;
        }
        let mut open = self.show_agent_assistant_dialog;
        let viewport_id = Self::agent_assistant_viewport_id();
        let title = self.tr("agent.title");
        let spec = self
            .hosted_window_spec_for_viewport(
                title.clone(),
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
                close_requested = self
                    .render_agent_assistant_contents_scrollable(ui, "agent_assistant_main_scroll")
                    .inner;
            });
            self.clear_viewport_foreground_request_after_render(viewport_id);
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
            title,
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
                crate::egui_compat::show_hosted_window(
                    &mut *ctx,
                    &viewport_spec,
                    &mut open,
                    |ui| {
                        close_requested = self
                            .render_agent_assistant_contents_scrollable(
                                ui,
                                "agent_assistant_main_scroll",
                            )
                            .inner;
                    },
                );
                if close_requested {
                    open = false;
                }
            } else {
                let mut close_requested = false;
                crate::egui_compat::show_central_panel(
                    &mut *ctx,
                    egui::CentralPanel::default(),
                    |ui| {
                        close_requested = self
                            .render_agent_assistant_contents_scrollable(
                                ui,
                                "agent_assistant_main_scroll",
                            )
                            .inner;
                    },
                );

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
        if trace.preflight_history.is_empty()
            && let Some(snapshot) = trace.preflight_snapshot.clone()
        {
            trace.preflight_history.push(snapshot);
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
        self.persist_project_metadata_values(&[(
            ROUTINE_DECISION_TRACES_METADATA_KEY,
            Some(value),
        )]);
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

    pub(super) fn agent_suggestion_fact_readiness(
        &self,
        expr: &serde_json::Value,
    ) -> Option<String> {
        let expression =
            serde_json::from_value::<crate::engine::FactExpression>(expr.clone()).ok()?;
        let engine = self.engine.read().ok()?;
        let evaluation = engine.evaluate_fact_expression(&expression, &[]);
        Some(crate::agent_bridge::agent_fact_readiness_label(&evaluation))
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
                if self.routine_assistant_compare_routine_id.trim().is_empty()
                    && let Some(alt_id) = output
                        .get("alternatives")
                        .and_then(|value| value.as_array())
                        .and_then(|rows| rows.first())
                        .and_then(|row| row.get("routine_id"))
                        .and_then(|value| value.as_str())
                {
                    self.routine_assistant_compare_routine_id = alt_id.trim().to_string();
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
