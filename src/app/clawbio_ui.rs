use super::clawbio_bridge::{
    CLAWBIO_POLL_INTERVAL, CancelToken, ClawBioResultView, ClawBioRunSpec, ClawBioTransport,
    DEFAULT_CLAWBIO_SKILL_ALIAS, SubprocessTransport, clawbio_run_root, default_clawbio_bin,
};
use super::*;
use serde_json::Value;
use std::{path::PathBuf, process::Command, sync::mpsc, thread};

pub(super) struct ClawBioPanelState {
    pub(super) show_panel: bool,
    skill_alias: String,
    clawbio_bin: String,
    task: Option<ClawBioTask>,
    status: String,
    last_request_json: String,
    last_result: Option<ClawBioResultView>,
    last_error: Option<String>,
}

impl Default for ClawBioPanelState {
    fn default() -> Self {
        Self {
            show_panel: false,
            skill_alias: DEFAULT_CLAWBIO_SKILL_ALIAS.to_string(),
            clawbio_bin: default_clawbio_bin(),
            task: None,
            status: String::new(),
            last_request_json: String::new(),
            last_result: None,
            last_error: None,
        }
    }
}

struct ClawBioTask {
    started: Instant,
    cancel: CancelToken,
    receiver: mpsc::Receiver<Result<ClawBioResultView, String>>,
}

pub(super) fn run_clawbio_on_worker(
    transport: Box<dyn ClawBioTransport + Send>,
    spec: ClawBioRunSpec,
    cancel: CancelToken,
    tx: mpsc::Sender<Result<ClawBioResultView, String>>,
) {
    thread::spawn(move || {
        let result = transport.dispatch(spec, cancel);
        let _ = tx.send(result);
    });
}

impl GENtleApp {
    pub(super) fn open_clawbio_dialog(&mut self) {
        let was_open = self.clawbio_panel.show_panel;
        self.clawbio_panel.show_panel = true;
        self.refresh_clawbio_request_preview();
        self.mark_window_open_or_focus(Self::clawbio_viewport_id(), was_open);
    }

    pub(super) fn poll_clawbio_task(&mut self, ctx: &egui::Context) {
        let mut finished = None;
        if let Some(task) = &self.clawbio_panel.task {
            match task.receiver.try_recv() {
                Ok(result) => finished = Some(result),
                Err(mpsc::TryRecvError::Empty) => {
                    ctx.request_repaint_after(CLAWBIO_POLL_INTERVAL);
                }
                Err(mpsc::TryRecvError::Disconnected) => {
                    finished = Some(Err("ClawBio worker disconnected.".to_string()));
                }
            }
        }
        if let Some(result) = finished {
            self.clawbio_panel.task = None;
            match result {
                Ok(result) => {
                    self.clawbio_panel.status = "ClawBio result loaded.".to_string();
                    self.clawbio_panel.last_error = None;
                    self.clawbio_panel.last_result = Some(result);
                }
                Err(err) => {
                    self.clawbio_panel.status = "ClawBio run failed.".to_string();
                    self.clawbio_panel.last_error = Some(err);
                }
            }
        }
    }

    pub(super) fn render_clawbio_dialog(&mut self, ctx: &egui::Context) {
        if !self.clawbio_panel.show_panel {
            return;
        }
        let mut open = self.clawbio_panel.show_panel;
        let viewport_id = Self::clawbio_viewport_id();
        let spec = self.hosted_window_spec_for_viewport(
            "ClawBio",
            Self::hosted_clawbio_window_id(),
            viewport_id,
            Vec2::new(920.0, 700.0),
            Vec2::new(620.0, 420.0),
        );
        crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
            self.render_clawbio_contents(ui);
        });
        self.clear_viewport_foreground_request_after_render(viewport_id);
        self.finalize_viewport_open_probe(viewport_id, "ClawBio");
        self.clawbio_panel.show_panel = open && self.clawbio_panel.show_panel;
    }

    fn render_clawbio_contents(&mut self, ui: &mut Ui) {
        if self.render_specialist_window_nav_with_close(ui, Some(("Close", "Close ClawBio"))) {
            self.clawbio_panel.show_panel = false;
        }
        ui.heading("ClawBio");
        ui.horizontal(|ui| {
            ui.label("Skill alias");
            ui.text_edit_singleline(&mut self.clawbio_panel.skill_alias);
        });
        ui.horizontal(|ui| {
            ui.label("clawbio.py");
            ui.text_edit_singleline(&mut self.clawbio_panel.clawbio_bin);
        });
        ui.horizontal(|ui| {
            let running = self.clawbio_panel.task.is_some();
            if ui
                .add_enabled(!running, egui::Button::new("Run"))
                .on_hover_text("Run ClawBio with the current sequence/project context")
                .clicked()
            {
                self.start_clawbio_request(self.current_clawbio_request());
            }
            if ui
                .add_enabled(!running, egui::Button::new("Refresh Context"))
                .clicked()
            {
                self.refresh_clawbio_request_preview();
            }
            if let Some(task) = &self.clawbio_panel.task {
                ui.spinner();
                ui.label(format!("Running for {}s", task.started.elapsed().as_secs()));
                if ui.button("Cancel").clicked() {
                    task.cancel.cancel();
                    self.clawbio_panel.status = "Cancelling ClawBio run...".to_string();
                }
            }
        });
        if !self.clawbio_panel.status.is_empty() {
            ui.small(self.clawbio_panel.status.as_str());
        }
        if let Some(err) = &self.clawbio_panel.last_error {
            ui.colored_label(egui::Color32::RED, err);
        }
        egui::CollapsingHeader::new("Request context")
            .default_open(false)
            .show(ui, |ui| {
                ui.add(
                    egui::TextEdit::multiline(&mut self.clawbio_panel.last_request_json)
                        .desired_rows(8)
                        .code_editor(),
                );
            });
        if let Some(result) = self.clawbio_panel.last_result.clone() {
            self.render_clawbio_result(ui, &result);
        }
    }

    fn render_clawbio_result(&mut self, ui: &mut Ui, result: &ClawBioResultView) {
        ui.separator();
        ui.label(clawbio_lifecycle_label(result));
        for line in &result.chat_summary_lines {
            ui.label(line);
        }
        ui.horizontal(|ui| {
            if ui.button("Reveal output...").clicked() {
                reveal_path(&result.output_dir);
            }
            ui.monospace(result.output_dir.as_str());
        });
        if let Some(report) = &result.report_md_path {
            if ui.button("Open report.md").clicked() {
                reveal_path(report);
            }
        }
        if ui.button("Open result.json").clicked() {
            reveal_path(&result.result_json_path);
        }
        for artifact in &result.preferred_artifacts {
            let path = clawbio_artifact_path(artifact, &result.output_dir);
            ui.horizontal(|ui| {
                if let Some(path) = path.as_ref() {
                    if ui.button(clawbio_artifact_label(artifact)).clicked() {
                        reveal_path(path);
                    }
                    ui.monospace(path);
                } else {
                    ui.label(clawbio_artifact_label(artifact));
                }
            });
        }
        for (idx, action) in result.suggested_actions.iter().enumerate() {
            let request = action.get("request").cloned();
            let enabled = clawbio_action_has_request(action) && self.clawbio_panel.task.is_none();
            ui.horizontal(|ui| {
                if ui
                    .add_enabled(
                        enabled,
                        egui::Button::new(clawbio_action_label(action, idx + 1)),
                    )
                    .clicked()
                {
                    if let Some(request) = request.clone() {
                        self.start_clawbio_request(request);
                    }
                }
                if !clawbio_action_has_request(action) {
                    ui.small("No nested request supplied by ClawBio.");
                }
            });
        }
        egui::CollapsingHeader::new("Raw result.json")
            .default_open(false)
            .show(ui, |ui| {
                let mut raw = serde_json::to_string_pretty(&result.raw_result)
                    .unwrap_or_else(|_| "{}".into());
                ui.add(
                    egui::TextEdit::multiline(&mut raw)
                        .desired_rows(10)
                        .code_editor(),
                );
            });
    }

    fn current_clawbio_request(&self) -> Value {
        self.clawbio_context_request(Some(self.effective_clawbio_skill_alias().as_str()))
    }

    fn refresh_clawbio_request_preview(&mut self) {
        self.clawbio_panel.last_request_json =
            serde_json::to_string_pretty(&self.current_clawbio_request()).unwrap_or_default();
    }

    fn start_clawbio_request(&mut self, request: Value) {
        if self.clawbio_panel.task.is_some() {
            return;
        }
        let spec = ClawBioRunSpec {
            clawbio_bin: self.effective_clawbio_bin(),
            skill_alias: self.effective_clawbio_skill_alias(),
            request: request.clone(),
            run_root: clawbio_run_root(),
        };
        self.clawbio_panel.last_request_json =
            serde_json::to_string_pretty(&request).unwrap_or_default();
        let (tx, rx) = mpsc::channel();
        let cancel = CancelToken::default();
        run_clawbio_on_worker(Box::new(SubprocessTransport), spec, cancel.clone(), tx);
        self.clawbio_panel.task = Some(ClawBioTask {
            started: Instant::now(),
            cancel,
            receiver: rx,
        });
        self.clawbio_panel.status = "Running ClawBio...".to_string();
        self.clawbio_panel.last_error = None;
    }

    fn effective_clawbio_skill_alias(&self) -> String {
        let trimmed = self.clawbio_panel.skill_alias.trim();
        if trimmed.is_empty() {
            DEFAULT_CLAWBIO_SKILL_ALIAS.to_string()
        } else {
            trimmed.to_string()
        }
    }

    fn effective_clawbio_bin(&self) -> String {
        let trimmed = self.clawbio_panel.clawbio_bin.trim();
        if trimmed.is_empty() {
            default_clawbio_bin()
        } else {
            trimmed.to_string()
        }
    }
}

pub(super) fn clawbio_lifecycle_label(result: &ClawBioResultView) -> String {
    let state = result.workflow_state.as_ref();
    let label = state
        .and_then(|v| v.get("state_label"))
        .and_then(Value::as_str)
        .unwrap_or("ClawBio result");
    let lifecycle = state
        .and_then(|v| v.get("lifecycle"))
        .and_then(Value::as_str)
        .or_else(|| result.raw_result.get("status").and_then(Value::as_str))
        .unwrap_or("unknown");
    format!("{label}: {lifecycle}")
}

pub(super) fn clawbio_artifact_label(artifact: &Value) -> String {
    artifact
        .get("label")
        .or_else(|| artifact.get("artifact_id"))
        .or_else(|| artifact.get("artifact_kind"))
        .and_then(Value::as_str)
        .unwrap_or("Artifact")
        .to_string()
}

pub(super) fn clawbio_artifact_path(artifact: &Value, output_dir: &str) -> Option<String> {
    artifact
        .get("path")
        .or_else(|| artifact.get("bundle_path"))
        .or_else(|| artifact.get("declared_path"))
        .and_then(Value::as_str)
        .map(|path| {
            let path = PathBuf::from(path);
            if path.is_absolute() {
                path.display().to_string()
            } else {
                PathBuf::from(output_dir).join(path).display().to_string()
            }
        })
}

pub(super) fn clawbio_action_label(action: &Value, index_1based: usize) -> String {
    let label = action
        .get("label")
        .and_then(Value::as_str)
        .unwrap_or("Action");
    match action.get("estimate").and_then(Value::as_str) {
        Some(estimate) if !estimate.trim().is_empty() => {
            format!("{index_1based}. {label} ({estimate})")
        }
        _ => format!("{index_1based}. {label}"),
    }
}

pub(super) fn clawbio_action_has_request(action: &Value) -> bool {
    matches!(action.get("request"), Some(Value::Object(_)))
}

fn reveal_path(path: &str) {
    let mut command = if cfg!(target_os = "macos") {
        let mut command = Command::new("open");
        command.arg(path);
        command
    } else if cfg!(target_os = "windows") {
        let mut command = Command::new("explorer");
        command.arg(path);
        command
    } else {
        let mut command = Command::new("xdg-open");
        command.arg(path);
        command
    };
    let _ = command.spawn();
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_view() -> ClawBioResultView {
        ClawBioResultView {
            output_dir: "/tmp/clawbio".to_string(),
            result_json_path: "/tmp/clawbio/result.json".to_string(),
            report_md_path: None,
            workflow_state: Some(serde_json::json!({
                "state_label": "Ready",
                "lifecycle": "ready"
            })),
            chat_summary_lines: vec![],
            preferred_artifacts: vec![],
            suggested_actions: vec![],
            raw_result: serde_json::json!({"status": "ok"}),
        }
    }

    #[test]
    fn clawbio_panel_defaults_to_gentle_cloning_alias() {
        let state = ClawBioPanelState::default();
        assert_eq!(state.skill_alias, DEFAULT_CLAWBIO_SKILL_ALIAS);
        assert!(!state.clawbio_bin.trim().is_empty());
    }

    #[test]
    fn clawbio_render_helpers_format_result_artifacts_and_actions() {
        let view = sample_view();
        assert_eq!(clawbio_lifecycle_label(&view), "Ready: ready");
        let artifact = serde_json::json!({"artifact_id": "fig", "path": "fig.svg"});
        assert_eq!(clawbio_artifact_label(&artifact), "fig");
        assert_eq!(
            clawbio_artifact_path(&artifact, "/tmp/run").unwrap(),
            "/tmp/run/fig.svg"
        );
        let action = serde_json::json!({
            "label": "Continue",
            "estimate": "~5s",
            "request": {"mode": "skill-info"}
        });
        assert_eq!(clawbio_action_label(&action, 2), "2. Continue (~5s)");
        assert!(clawbio_action_has_request(&action));
        assert!(!clawbio_action_has_request(
            &serde_json::json!({"label": "No request"})
        ));
    }
}
