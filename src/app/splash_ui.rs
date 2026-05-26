//! Startup splash screen rendering for the root GUI.

use super::*;
use crate::icons::SPLASH_SCREEN;

const SPLASH_AUTO_HIDE_AFTER: Duration = Duration::from_millis(2_200);
const SPLASH_REPAINT_INTERVAL: Duration = Duration::from_millis(80);

impl GENtleApp {
    pub(super) fn splash_should_render_at(&self, now: Instant) -> bool {
        !self.splash_dismissed
            && now.saturating_duration_since(self.splash_started_at) < SPLASH_AUTO_HIDE_AFTER
    }

    pub(super) fn dismiss_splash_screen(&mut self) {
        self.splash_dismissed = true;
    }

    pub(super) fn render_splash_screen(&mut self, ctx: &egui::Context) {
        let now = Instant::now();
        if !self.splash_should_render_at(now) {
            self.dismiss_splash_screen();
            return;
        }
        if ctx.input(|input| {
            input.pointer.any_click()
                || input.key_pressed(Key::Escape)
                || input.key_pressed(Key::Enter)
                || input.key_pressed(Key::Space)
        }) {
            self.dismiss_splash_screen();
            ctx.request_repaint();
            return;
        }

        ctx.request_repaint_after(SPLASH_REPAINT_INTERVAL);
        crate::egui_compat::show_central_panel(
            ctx,
            egui::CentralPanel::default().frame(egui::Frame::NONE),
            |ui| {
                let rect = ui.max_rect();
                let palette = SciencePalette::default();
                let fill = egui::Color32::from_rgb(239, 237, 230);
                let accent = palette.sequence_node;
                let reverse = palette.reverse_strand;
                ui.painter().rect_filled(rect, 0.0, fill);

                let dna_y = rect.center().y + 120.0;
                let left = rect.left() + rect.width() * 0.18;
                let right = rect.right() - rect.width() * 0.18;
                let width = (right - left).max(160.0);
                let step = (width / 9.0).max(24.0);
                for idx in 0..=9 {
                    let x = left + step * idx as f32;
                    let phase = idx as f32 * 0.65;
                    let top = dna_y + phase.sin() * 18.0;
                    let bottom = dna_y - phase.sin() * 18.0;
                    ui.painter().line_segment(
                        [egui::pos2(x, top), egui::pos2(x, bottom)],
                        egui::Stroke::new(1.0, egui::Color32::from_gray(165)),
                    );
                }
                let mut points_a = Vec::new();
                let mut points_b = Vec::new();
                let segments = 80;
                for idx in 0..=segments {
                    let t = idx as f32 / segments as f32;
                    let x = left + width * t;
                    let wave = (t * std::f32::consts::TAU * 2.2).sin() * 18.0;
                    points_a.push(egui::pos2(x, dna_y + wave));
                    points_b.push(egui::pos2(x, dna_y - wave));
                }
                ui.painter()
                    .add(egui::Shape::line(points_a, egui::Stroke::new(2.0, accent)));
                ui.painter()
                    .add(egui::Shape::line(points_b, egui::Stroke::new(2.0, reverse)));

                ui.vertical_centered(|ui| {
                    let top_space = (ui.available_height() - 330.0).max(28.0) * 0.45;
                    ui.add_space(top_space);
                    ui.add(
                        SPLASH_SCREEN
                            .clone()
                            .fit_to_exact_size(Vec2::new(112.0, 112.0)),
                    );
                    ui.add_space(14.0);
                    ui.heading("GENtle");
                    ui.add_space(4.0);
                    ui.label(format!(
                        "Molecular biology workbench · v{}",
                        about::GENTLE_DISPLAY_VERSION
                    ));
                    ui.add_space(18.0);
                    ui.spinner();
                    ui.add_space(8.0);
                    ui.small("Initialising workspace");
                    ui.add_space(6.0);
                    ui.small("Click or press Esc to continue");
                });
            },
        );
    }
}
