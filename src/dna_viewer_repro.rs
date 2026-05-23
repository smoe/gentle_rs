//! Diagnostic DNA-viewer repro harness used to bisect embedded-window layout.
//!
//! This module intentionally keeps the biology/UI implementation in the real
//! sequence-window types and exposes only a small mode switcher for local GUI
//! diagnosis. It is not part of the stable adapter contract.

use crate::{
    dna_sequence::{DNAsequence, load_from_file},
    egui_compat::{self, HostedWindowSpec},
    engine::GentleEngine,
    main_area_dna::MainAreaDna,
    window_dna::WindowDna,
};
use eframe::{Frame, egui};
use std::{
    path::Path,
    str::FromStr,
    sync::{Arc, RwLock},
};

const DEFAULT_HOSTED_SIZE: egui::Vec2 = egui::vec2(1200.0, 860.0);
const MIN_HOSTED_SIZE: egui::Vec2 = egui::vec2(560.0, 360.0);
const DEFAULT_SYNTHETIC_LEN_BP: usize = 120_000;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum DnaViewerReproMode {
    HostedShell,
    HostedWindow,
    HostedMain,
    HostedTop,
    HostedMiddle,
    HostedSequence,
    HostedTopMiddle,
    HostedMiddleSequence,
    DirectWindow,
    DirectMain,
    DirectMiddle,
}

impl DnaViewerReproMode {
    pub const ALL: [Self; 11] = [
        Self::HostedShell,
        Self::HostedWindow,
        Self::HostedMain,
        Self::HostedTop,
        Self::HostedMiddle,
        Self::HostedSequence,
        Self::HostedTopMiddle,
        Self::HostedMiddleSequence,
        Self::DirectWindow,
        Self::DirectMain,
        Self::DirectMiddle,
    ];

    pub fn name(self) -> &'static str {
        match self {
            Self::HostedShell => "hosted-shell",
            Self::HostedWindow => "hosted-window",
            Self::HostedMain => "hosted-main",
            Self::HostedTop => "hosted-top",
            Self::HostedMiddle => "hosted-middle",
            Self::HostedSequence => "hosted-sequence",
            Self::HostedTopMiddle => "hosted-top-middle",
            Self::HostedMiddleSequence => "hosted-middle-sequence",
            Self::DirectWindow => "direct-window",
            Self::DirectMain => "direct-main",
            Self::DirectMiddle => "direct-middle",
        }
    }

    pub fn description(self) -> &'static str {
        match self {
            Self::HostedShell => "hosted egui::Window shell with diagnostics only",
            Self::HostedWindow => "real WindowDna::update_embedded path inside hosted shell",
            Self::HostedMain => {
                "MainAreaDna::render_inside_without_auxiliary_windows inside hosted shell"
            }
            Self::HostedTop => "top controls only inside hosted shell",
            Self::HostedMiddle => "map/sidebar middle body only inside hosted shell",
            Self::HostedSequence => "sequence strip only inside hosted shell",
            Self::HostedTopMiddle => "top controls plus middle body inside hosted shell",
            Self::HostedMiddleSequence => "middle body plus sequence strip inside hosted shell",
            Self::DirectWindow => {
                "real WindowDna::update_embedded path directly in root CentralPanel"
            }
            Self::DirectMain => {
                "MainAreaDna::render_inside_without_auxiliary_windows directly in root CentralPanel"
            }
            Self::DirectMiddle => "map/sidebar middle body directly in root CentralPanel",
        }
    }

    fn uses_hosted_window(self) -> bool {
        matches!(
            self,
            Self::HostedShell
                | Self::HostedWindow
                | Self::HostedMain
                | Self::HostedTop
                | Self::HostedMiddle
                | Self::HostedSequence
                | Self::HostedTopMiddle
                | Self::HostedMiddleSequence
        )
    }
}

impl FromStr for DnaViewerReproMode {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        let normalized = value.trim().to_ascii_lowercase().replace('_', "-");
        match normalized.as_str() {
            "hosted-shell" | "shell" => Ok(Self::HostedShell),
            "hosted-window" | "window" | "full" | "hosted-full" => Ok(Self::HostedWindow),
            "hosted-main" | "main" => Ok(Self::HostedMain),
            "hosted-top" | "top" => Ok(Self::HostedTop),
            "hosted-middle" | "middle" => Ok(Self::HostedMiddle),
            "hosted-sequence" | "sequence" | "bottom" => Ok(Self::HostedSequence),
            "hosted-top-middle" | "top-middle" => Ok(Self::HostedTopMiddle),
            "hosted-middle-sequence" | "middle-sequence" => Ok(Self::HostedMiddleSequence),
            "direct-window" => Ok(Self::DirectWindow),
            "direct-main" => Ok(Self::DirectMain),
            "direct-middle" => Ok(Self::DirectMiddle),
            _ => Err(format!(
                "unknown DNA viewer repro mode '{value}'. Valid modes: {}",
                mode_names().join(", ")
            )),
        }
    }
}

pub fn mode_names() -> Vec<&'static str> {
    DnaViewerReproMode::ALL
        .iter()
        .map(|mode| mode.name())
        .collect()
}

pub fn mode_help() -> String {
    DnaViewerReproMode::ALL
        .iter()
        .map(|mode| format!("  {:24} {}", mode.name(), mode.description()))
        .collect::<Vec<_>>()
        .join("\n")
}

pub fn parse_mode_list(values: &[String]) -> Result<Vec<DnaViewerReproMode>, String> {
    let mut modes = Vec::new();
    for value in values {
        for part in value.split(',') {
            let trimmed = part.trim();
            if trimmed.is_empty() {
                continue;
            }
            let mode = DnaViewerReproMode::from_str(trimmed)?;
            if !modes.contains(&mode) {
                modes.push(mode);
            }
        }
    }
    if modes.is_empty() {
        modes.push(DnaViewerReproMode::HostedWindow);
    }
    Ok(modes)
}

pub fn load_repro_sequence(
    path: Option<&Path>,
    synthetic_len_bp: Option<usize>,
) -> Result<(DNAsequence, String), String> {
    if let Some(path) = path {
        let path_text = path.to_string_lossy();
        let dna = load_from_file(&path_text)
            .map_err(|err| format!("could not load '{}': {err}", path.display()))?;
        let label = dna
            .name()
            .as_deref()
            .filter(|name| !name.trim().is_empty())
            .map(ToOwned::to_owned)
            .unwrap_or_else(|| path.display().to_string());
        return Ok((dna, label));
    }

    let len = synthetic_len_bp.unwrap_or(DEFAULT_SYNTHETIC_LEN_BP).max(1);
    let pattern = b"ACGTGCAATTCG";
    let mut sequence = String::with_capacity(len);
    for idx in 0..len {
        sequence.push(pattern[idx % pattern.len()] as char);
    }
    let mut dna = DNAsequence::from_sequence(&sequence)
        .map_err(|err| format!("could not create synthetic sequence: {err}"))?;
    dna.set_name(format!("synthetic_{len}bp"));
    Ok((dna, format!("synthetic {len} bp")))
}

pub struct DnaViewerReproApp {
    modes: Vec<DnaViewerReproMode>,
    active_mode_index: usize,
    sequence_label: String,
    sequence_len_bp: usize,
    hosted_open: bool,
    window_shell: WindowDna,
    main_area: MainAreaDna,
}

impl DnaViewerReproApp {
    pub fn new(dna: DNAsequence, sequence_label: String, modes: Vec<DnaViewerReproMode>) -> Self {
        let sequence_len_bp = dna.len();
        let engine = Arc::new(RwLock::new(GentleEngine::new()));
        let mut window_shell = WindowDna::new(dna.clone(), sequence_label.clone(), engine.clone());
        window_shell.set_window_scope_id("dna_viewer_repro_window".to_string());
        let mut main_area = MainAreaDna::new(dna, Some(sequence_label.clone()), Some(engine));
        main_area.set_window_scope_id("dna_viewer_repro_main".to_string());
        Self {
            modes: if modes.is_empty() {
                vec![DnaViewerReproMode::HostedWindow]
            } else {
                modes
            },
            active_mode_index: 0,
            sequence_label,
            sequence_len_bp,
            hosted_open: true,
            window_shell,
            main_area,
        }
    }

    fn active_mode(&self) -> DnaViewerReproMode {
        self.modes
            .get(self.active_mode_index)
            .copied()
            .unwrap_or(DnaViewerReproMode::HostedWindow)
    }

    fn render_controls(&mut self, ui: &mut egui::Ui) {
        ui.horizontal_wrapped(|ui| {
            ui.heading("DNA viewer repro");
            ui.separator();
            ui.label(format!(
                "{} ({} bp)",
                self.sequence_label, self.sequence_len_bp
            ));
            ui.separator();
            if ui.button("reopen hosted window").clicked() {
                self.hosted_open = true;
            }
        });
        ui.horizontal_wrapped(|ui| {
            ui.label("mode:");
            for index in 0..self.modes.len() {
                let mode = self.modes[index];
                if ui
                    .selectable_label(index == self.active_mode_index, mode.name())
                    .on_hover_text(mode.description())
                    .clicked()
                {
                    self.active_mode_index = index;
                    if mode.uses_hosted_window() {
                        self.hosted_open = true;
                    }
                }
            }
        });
    }

    fn render_root_surface(&self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();
        let content_rect = ctx.content_rect();
        let safe = egui_compat::hosted_window_safe_rect(&ctx);
        let max_inner =
            egui_compat::hosted_window_max_inner_size(safe, MIN_HOSTED_SIZE, egui::Vec2::ZERO);
        let hosted_area_rect = if self.active_mode().uses_hosted_window() {
            let stable_id = hosted_stable_id(self.active_mode());
            ctx.memory(|mem| mem.area_rect(stable_id))
        } else {
            None
        };
        ui.heading("Root surface");
        ui.label("Resize this app window and the hosted DNA window. The numbers below are deliberately live so a visual desync has coordinates attached to it.");
        ui.separator();
        egui::Grid::new("dna_viewer_repro_root_geometry")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                ui.label("content rect");
                ui.monospace(format_rect(content_rect));
                ui.end_row();
                ui.label("host safe rect");
                ui.monospace(format_rect(safe));
                ui.end_row();
                ui.label("host max inner");
                ui.monospace(format!("{:.1} x {:.1}", max_inner.x, max_inner.y));
                ui.end_row();
                ui.label("active mode");
                ui.monospace(self.active_mode().name());
                ui.end_row();
            });
        if self.active_mode().uses_hosted_window() {
            ui.label(format!(
                "hosted area rect: {}",
                hosted_area_rect
                    .map(format_rect)
                    .unwrap_or_else(|| "n/a".to_string())
            ));
        }
        ui.separator();
    }

    fn render_hosted_window(&mut self, ctx: &egui::Context, mode: DnaViewerReproMode) {
        let title = format!("DNA viewer repro: {}", mode.name());
        let spec = HostedWindowSpec::new(
            title,
            hosted_stable_id(mode),
            DEFAULT_HOSTED_SIZE,
            MIN_HOSTED_SIZE,
        )
        .initial_pos(Some(egui::pos2(42.0, 86.0)))
        .drag_margin(egui::Vec2::ZERO)
        .resizable(true);
        let mut hosted_open = self.hosted_open;
        egui_compat::show_hosted_window(ctx, &spec, &mut hosted_open, |ui| {
            self.render_mode_body(ui, mode);
        });
        self.hosted_open = hosted_open;
    }

    fn render_mode_body(&mut self, ui: &mut egui::Ui, mode: DnaViewerReproMode) {
        self.render_body_geometry_header(ui, mode);
        match mode {
            DnaViewerReproMode::HostedShell => {
                ui.separator();
                ui.label("Hosted shell only. If this cannot resize, the problem is outside the DNA viewer body.");
                fill_diagnostic_canvas(ui);
            }
            DnaViewerReproMode::HostedWindow | DnaViewerReproMode::DirectWindow => {
                self.window_shell.update_embedded(ui);
            }
            DnaViewerReproMode::HostedMain | DnaViewerReproMode::DirectMain => {
                self.main_area.render_inside_without_auxiliary_windows(ui);
            }
            DnaViewerReproMode::HostedTop => {
                self.main_area.render_top_panel(ui);
            }
            DnaViewerReproMode::HostedMiddle | DnaViewerReproMode::DirectMiddle => {
                let ctx = ui.ctx().clone();
                self.main_area.render_middle(&ctx, ui);
            }
            DnaViewerReproMode::HostedSequence => {
                self.main_area.render_sequence(ui);
            }
            DnaViewerReproMode::HostedTopMiddle => {
                egui_compat::show_top_panel_inside(
                    ui,
                    egui::Panel::top("dna_repro_top_middle_top"),
                    |ui| self.main_area.render_top_panel(ui),
                );
                egui_compat::show_central_panel_inside(ui, egui::CentralPanel::default(), |ui| {
                    let ctx = ui.ctx().clone();
                    self.main_area.render_middle(&ctx, ui);
                });
            }
            DnaViewerReproMode::HostedMiddleSequence => {
                egui_compat::show_bottom_panel_inside(
                    ui,
                    egui::Panel::bottom("dna_repro_middle_sequence_bottom")
                        .default_size(180.0)
                        .resizable(true),
                    |ui| self.main_area.render_sequence(ui),
                );
                egui_compat::show_central_panel_inside(ui, egui::CentralPanel::default(), |ui| {
                    let ctx = ui.ctx().clone();
                    self.main_area.render_middle(&ctx, ui);
                });
            }
        }
    }

    fn render_body_geometry_header(&self, ui: &mut egui::Ui, mode: DnaViewerReproMode) {
        let available = ui.available_rect_before_wrap();
        ui.horizontal_wrapped(|ui| {
            ui.label("body mode");
            ui.monospace(mode.name());
            ui.separator();
            ui.label("available");
            ui.monospace(format_rect(available));
            ui.separator();
            ui.label("clip");
            ui.monospace(format_rect(ui.clip_rect()));
        });
    }
}

impl eframe::App for DnaViewerReproApp {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut Frame) {
        self.render_controls(ui);
        let mode = self.active_mode();
        if mode.uses_hosted_window() {
            self.render_root_surface(ui);
            if self.hosted_open {
                let ctx = ui.ctx().clone();
                self.render_hosted_window(&ctx, mode);
            }
        } else {
            ui.separator();
            self.render_mode_body(ui, mode);
        }
    }
}

fn hosted_stable_id(mode: DnaViewerReproMode) -> egui::Id {
    egui::Id::new(("dna_viewer_repro_hosted", mode.name()))
}

fn format_rect(rect: egui::Rect) -> String {
    format!(
        "min=({:.1},{:.1}) max=({:.1},{:.1}) size={:.1}x{:.1}",
        rect.min.x,
        rect.min.y,
        rect.max.x,
        rect.max.y,
        rect.width(),
        rect.height()
    )
}

fn fill_diagnostic_canvas(ui: &mut egui::Ui) {
    let size = egui::vec2(
        ui.available_width().max(240.0),
        ui.available_height().max(180.0),
    );
    let (response, painter) = ui.allocate_painter(size, egui::Sense::hover());
    let rect = response.rect;
    painter.rect_filled(rect, 8.0, egui::Color32::from_rgb(248, 250, 252));
    painter.rect_stroke(
        rect,
        8.0,
        egui::Stroke::new(1.0, egui::Color32::from_rgb(148, 163, 184)),
        egui::StrokeKind::Middle,
    );
    painter.text(
        rect.center(),
        egui::Align2::CENTER_CENTER,
        "hosted shell diagnostic canvas",
        egui::FontId::proportional(18.0),
        egui::Color32::from_rgb(71, 85, 105),
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_mode_list_accepts_repeated_and_comma_separated_modes() {
        let modes = parse_mode_list(&[
            "hosted-shell,hosted-window".to_string(),
            "direct-middle".to_string(),
            "hosted-shell".to_string(),
        ])
        .expect("modes");
        assert_eq!(
            modes,
            vec![
                DnaViewerReproMode::HostedShell,
                DnaViewerReproMode::HostedWindow,
                DnaViewerReproMode::DirectMiddle,
            ]
        );
    }

    #[test]
    fn load_repro_sequence_can_create_synthetic_sequence() {
        let (dna, label) = load_repro_sequence(None, Some(37)).expect("synthetic DNA");
        assert_eq!(dna.len(), 37);
        assert_eq!(label, "synthetic 37 bp");
    }
}
