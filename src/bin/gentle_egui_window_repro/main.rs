//! Minimal native-viewport repro harness for macOS/egui window lifecycle bugs.
//!
//! This binary intentionally avoids GENtle biology/engine state and focuses on
//! one thing only:
//! - a root window with a continuously animated surface,
//! - a button that opens another native window,
//! - each child window exposes the same functionality again.
//!
//! The goal is to reproduce stale/respawn/twin-window behavior under resize or
//! maximize without the full DNA-sequence window stack getting in the way.

use eframe::egui::{self, Color32, Pos2, Sense, Stroke, Vec2};
use eframe::{Frame, NativeOptions, Renderer};
use std::collections::BTreeMap;
use std::sync::{Arc, Mutex, RwLock};
use std::time::Duration;

const DEFAULT_WINDOW_SIZE: [f32; 2] = [920.0, 640.0];
const MIN_WINDOW_SIZE: [f32; 2] = [480.0, 320.0];
const CASCADE_START: f32 = 64.0;
const CASCADE_STEP: f32 = 36.0;
const REPAINT_DELAY_MS: u64 = 16;

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct WindowActions {
    open_child: bool,
    close_self: bool,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum ChildViewportMode {
    Immediate,
    Deferred,
    Embedded,
}

impl ChildViewportMode {
    fn as_str(self) -> &'static str {
        match self {
            Self::Immediate => "immediate native viewport",
            Self::Deferred => "deferred native viewport",
            Self::Embedded => "embedded egui::Window",
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum ReproAction {
    OpenChild,
    CloseViewport(egui::ViewportId),
}

#[derive(Debug)]
struct ReproWindowState {
    serial: usize,
    opened_children: usize,
}

impl ReproWindowState {
    fn new(serial: usize) -> Self {
        Self {
            serial,
            opened_children: 0,
        }
    }

    fn title(&self) -> String {
        format!("GENtle egui repro [{}]", self.serial)
    }

    fn render_contents(
        &mut self,
        ui: &mut egui::Ui,
        allow_close: bool,
        ctx: &egui::Context,
        animate: bool,
        allow_nested_children: bool,
    ) -> WindowActions {
        if animate {
            ctx.request_repaint_after(Duration::from_millis(REPAINT_DELAY_MS));
        }
        let mut actions = WindowActions::default();
        let window_size = ctx.content_rect().size();
        let t = ctx.input(|i| i.time) as f32;

        ui.heading(self.title());
        ui.horizontal_wrapped(|ui| {
            if allow_nested_children {
                if ui
                    .button("Open another window")
                    .on_hover_text("Spawn another native viewport with the same controls")
                    .clicked()
                {
                    self.opened_children = self.opened_children.saturating_add(1);
                    actions.open_child = true;
                }
            }
            if allow_close
                && ui
                    .button("Close")
                    .on_hover_text("Close this native child window")
                    .clicked()
            {
                actions.close_self = true;
            }
            ui.separator();
            ui.small(format!(
                "viewport {:.0} x {:.0} px | spawned children {}{}",
                window_size.x,
                window_size.y,
                self.opened_children,
                if animate {
                    format!(" | t={:.2}s", t)
                } else {
                    String::new()
                }
            ));
        });
        ui.small(
            "Use this app for resize/maximize testing. If a child window respawns, duplicates, or leaves a stale twin, the issue is likely in native viewport lifecycle handling rather than GENtle biology logic.",
        );
        ui.separator();

        let canvas_height = (ui.available_height() - 8.0).max(220.0);
        let desired = Vec2::new(ui.available_width().max(320.0), canvas_height);
        let (response, painter) = ui.allocate_painter(desired, Sense::hover());
        let rect = response.rect;
        let bg = Color32::from_rgb(250, 246, 236);
        painter.rect_filled(rect, 12.0, bg);
        painter.rect_stroke(
            rect,
            12.0,
            Stroke::new(1.0, Color32::from_rgb(186, 180, 164)),
            egui::StrokeKind::Middle,
        );

        if animate {
            let stripe_count = 10usize;
            for idx in 0..stripe_count {
                let frac = idx as f32 / stripe_count as f32;
                let y = egui::lerp(rect.top()..=rect.bottom(), frac);
                let phase = t * 0.9 + self.serial as f32 * 0.11 + idx as f32 * 0.23;
                let bend = phase.sin() * rect.width() * 0.07;
                let start = Pos2::new(rect.left() + 18.0, y);
                let mid = Pos2::new(rect.center().x + bend, y + phase.cos() * 8.0);
                let end = Pos2::new(rect.right() - 18.0, y);
                painter.line_segment(
                    [start, mid],
                    Stroke::new(1.6, Color32::from_rgb(120, 127, 136)),
                );
                painter.line_segment(
                    [mid, end],
                    Stroke::new(1.6, Color32::from_rgb(120, 127, 136)),
                );
            }

            let x = egui::lerp(
                rect.left() + 24.0..=rect.right() - 24.0,
                (t * 0.8 + self.serial as f32 * 0.19).sin() * 0.5 + 0.5,
            );
            let y = egui::lerp(
                rect.top() + 24.0..=rect.bottom() - 24.0,
                (t * 1.15 + self.serial as f32 * 0.13).cos() * 0.5 + 0.5,
            );
            let color = match self.serial % 4 {
                0 => Color32::from_rgb(30, 64, 175),
                1 => Color32::from_rgb(180, 83, 9),
                2 => Color32::from_rgb(21, 128, 61),
                _ => Color32::from_rgb(190, 24, 93),
            };
            painter.circle_filled(Pos2::new(x, y), 16.0, color);
        } else {
            painter.text(
                rect.center(),
                egui::Align2::CENTER_CENTER,
                "static surface",
                egui::FontId::proportional(24.0),
                Color32::from_rgb(71, 85, 105),
            );
        }
        painter.text(
            rect.left_top() + Vec2::new(18.0, 16.0),
            egui::Align2::LEFT_TOP,
            format!("window {}", self.serial),
            egui::FontId::proportional(18.0),
            Color32::from_rgb(55, 65, 81),
        );

        actions
    }
}

struct ReproApp {
    child_mode: ChildViewportMode,
    animate_children: bool,
    cascade_new_windows: bool,
    nested_children: bool,
    next_serial: usize,
    root_window: ReproWindowState,
    child_windows: BTreeMap<egui::ViewportId, Arc<RwLock<ReproWindowState>>>,
    pending_window_positions: BTreeMap<egui::ViewportId, Pos2>,
    pending_actions: Arc<Mutex<Vec<ReproAction>>>,
}

impl ReproApp {
    fn new(_cc: &eframe::CreationContext<'_>) -> Self {
        Self {
            child_mode: ChildViewportMode::Immediate,
            animate_children: true,
            cascade_new_windows: true,
            nested_children: true,
            next_serial: 1,
            root_window: ReproWindowState::new(0),
            child_windows: BTreeMap::new(),
            pending_window_positions: BTreeMap::new(),
            pending_actions: Arc::new(Mutex::new(Vec::new())),
        }
    }

    fn child_viewport_id(serial: usize) -> egui::ViewportId {
        egui::ViewportId::from_hash_of(("gentle_egui_window_repro_child", serial))
    }

    fn cascade_position(index: usize) -> Pos2 {
        Pos2::new(
            CASCADE_START + index as f32 * CASCADE_STEP,
            CASCADE_START + index as f32 * CASCADE_STEP,
        )
    }

    fn spawn_child_window(&mut self) {
        let serial = self.next_serial;
        self.next_serial += 1;
        let viewport_id = Self::child_viewport_id(serial);
        self.child_windows.insert(
            viewport_id,
            Arc::new(RwLock::new(ReproWindowState::new(serial))),
        );
        if self.cascade_new_windows {
            self.pending_window_positions.insert(
                viewport_id,
                Self::cascade_position(serial.saturating_sub(1)),
            );
        }
        eprintln!("gentle_egui_window_repro: opened child window {serial} ({viewport_id:?})");
    }

    fn reset_children(&mut self) {
        self.child_windows.clear();
        self.pending_window_positions.clear();
        if let Ok(mut actions) = self.pending_actions.lock() {
            actions.clear();
        }
    }

    fn push_action(&self, action: ReproAction) {
        if let Ok(mut actions) = self.pending_actions.lock() {
            actions.push(action);
        }
    }

    fn render_child_window_immediate(
        &mut self,
        ctx: &egui::Context,
        viewport_id: egui::ViewportId,
    ) {
        let Some(window_state) = self.child_windows.get(&viewport_id).cloned() else {
            return;
        };
        let Some(title) = window_state.read().ok().map(|w| w.title()) else {
            return;
        };
        let builder = egui::ViewportBuilder::default()
            .with_title(title)
            .with_inner_size(DEFAULT_WINDOW_SIZE)
            .with_min_inner_size(MIN_WINDOW_SIZE);
        let pending_actions = self.pending_actions.clone();
        let animate = self.animate_children;
        let allow_nested_children = self.nested_children;

        ctx.show_viewport_immediate(viewport_id, builder, |ui, class| {
            if !matches!(
                class,
                egui::ViewportClass::Immediate | egui::ViewportClass::EmbeddedWindow
            ) {
                return;
            }

            let window_ctx = ui.ctx().clone();
            let native_close_requested = window_ctx.input(|i| i.viewport().close_requested());
            let mut actions = WindowActions::default();
            if let Ok(mut window) = window_state.write() {
                actions =
                    window.render_contents(ui, true, &window_ctx, animate, allow_nested_children);
            }

            if let Ok(mut queued) = pending_actions.lock() {
                if actions.open_child {
                    queued.push(ReproAction::OpenChild);
                }
                if actions.close_self || native_close_requested {
                    queued.push(ReproAction::CloseViewport(viewport_id));
                }
            }
        });

        if let Some(position) = self.pending_window_positions.remove(&viewport_id) {
            ctx.send_viewport_cmd_to(viewport_id, egui::ViewportCommand::OuterPosition(position));
        }
    }

    fn render_child_window_deferred(&mut self, ctx: &egui::Context, viewport_id: egui::ViewportId) {
        let Some(window_state) = self.child_windows.get(&viewport_id).cloned() else {
            return;
        };
        let Some(title) = window_state.read().ok().map(|w| w.title()) else {
            return;
        };
        let builder = egui::ViewportBuilder::default()
            .with_title(title)
            .with_inner_size(DEFAULT_WINDOW_SIZE)
            .with_min_inner_size(MIN_WINDOW_SIZE);
        let pending_actions = self.pending_actions.clone();
        let animate = self.animate_children;
        let allow_nested_children = self.nested_children;

        ctx.show_viewport_deferred(viewport_id, builder, move |ui, class| {
            if !matches!(
                class,
                egui::ViewportClass::Deferred | egui::ViewportClass::EmbeddedWindow
            ) {
                return;
            }
            let window_ctx = ui.ctx().clone();
            let native_close_requested = window_ctx.input(|i| i.viewport().close_requested());
            let mut actions = WindowActions::default();
            if let Ok(mut window) = window_state.write() {
                actions =
                    window.render_contents(ui, true, &window_ctx, animate, allow_nested_children);
            }
            if let Ok(mut queued) = pending_actions.lock() {
                if actions.open_child {
                    queued.push(ReproAction::OpenChild);
                }
                if actions.close_self || native_close_requested {
                    queued.push(ReproAction::CloseViewport(viewport_id));
                }
            }
        });

        if let Some(position) = self.pending_window_positions.remove(&viewport_id) {
            ctx.send_viewport_cmd_to(viewport_id, egui::ViewportCommand::OuterPosition(position));
        }
    }

    fn render_child_window_embedded(
        &mut self,
        ctx: &egui::Context,
        viewport_id: egui::ViewportId,
        state: Arc<RwLock<ReproWindowState>>,
    ) {
        let Some(title) = state.read().ok().map(|w| w.title()) else {
            return;
        };
        let mut open = true;
        egui::Window::new(title)
            .id(egui::Id::new((
                "gentle_egui_window_repro_embedded",
                viewport_id,
            )))
            .open(&mut open)
            .default_size(Vec2::from(DEFAULT_WINDOW_SIZE))
            .min_width(MIN_WINDOW_SIZE[0])
            .min_height(MIN_WINDOW_SIZE[1])
            .resizable(true)
            .show(ctx, |ui| {
                let window_ctx = ui.ctx().clone();
                if let Ok(mut window) = state.write() {
                    let actions = window.render_contents(
                        ui,
                        true,
                        &window_ctx,
                        self.animate_children,
                        self.nested_children,
                    );
                    if actions.open_child {
                        self.push_action(ReproAction::OpenChild);
                    }
                    if actions.close_self {
                        self.push_action(ReproAction::CloseViewport(viewport_id));
                    }
                }
            });
        if !open {
            self.push_action(ReproAction::CloseViewport(viewport_id));
        }
    }

    fn apply_pending_actions(&mut self) {
        let actions = if let Ok(mut queued) = self.pending_actions.lock() {
            std::mem::take(&mut *queued)
        } else {
            Vec::new()
        };
        for action in actions {
            match action {
                ReproAction::OpenChild => self.spawn_child_window(),
                ReproAction::CloseViewport(viewport_id) => {
                    self.child_windows.remove(&viewport_id);
                    self.pending_window_positions.remove(&viewport_id);
                }
            }
        }
    }
}

impl eframe::App for ReproApp {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut Frame) {
        let ctx = ui.ctx().clone();
        ui.horizontal_wrapped(|ui| {
            ui.label("child window mode");
            let old_mode = self.child_mode;
            egui::ComboBox::from_id_salt("gentle_egui_window_repro_mode")
                .selected_text(self.child_mode.as_str())
                .show_ui(ui, |ui| {
                    ui.selectable_value(
                        &mut self.child_mode,
                        ChildViewportMode::Immediate,
                        ChildViewportMode::Immediate.as_str(),
                    );
                    ui.selectable_value(
                        &mut self.child_mode,
                        ChildViewportMode::Deferred,
                        ChildViewportMode::Deferred.as_str(),
                    );
                    ui.selectable_value(
                        &mut self.child_mode,
                        ChildViewportMode::Embedded,
                        ChildViewportMode::Embedded.as_str(),
                    );
                });
            if self.child_mode != old_mode {
                self.reset_children();
            }
            ui.small("Changing mode resets existing child windows so each test starts cleanly.");
        });
        ui.horizontal_wrapped(|ui| {
            let mut reset_needed = false;
            reset_needed |= ui
                .checkbox(&mut self.animate_children, "animate child surfaces")
                .changed();
            reset_needed |= ui
                .checkbox(&mut self.cascade_new_windows, "cascade initial position")
                .changed();
            reset_needed |= ui
                .checkbox(&mut self.nested_children, "child can open child")
                .changed();
            if reset_needed {
                self.reset_children();
            }
            ui.small("Changing reduction knobs resets existing child windows.");
        });
        ui.separator();
        let root_actions = self
            .root_window
            .render_contents(ui, false, &ctx, true, true);
        if root_actions.open_child {
            self.spawn_child_window();
        }

        let child_entries = self
            .child_windows
            .iter()
            .map(|(id, state)| (*id, state.clone()))
            .collect::<Vec<_>>();
        for (viewport_id, state) in child_entries {
            match self.child_mode {
                ChildViewportMode::Immediate => {
                    self.render_child_window_immediate(&ctx, viewport_id)
                }
                ChildViewportMode::Deferred => self.render_child_window_deferred(&ctx, viewport_id),
                ChildViewportMode::Embedded => {
                    self.render_child_window_embedded(&ctx, viewport_id, state)
                }
            }
        }
        self.apply_pending_actions();
    }
}

fn main() -> eframe::Result {
    let options = NativeOptions {
        persist_window: false,
        renderer: Renderer::Wgpu,
        viewport: egui::ViewportBuilder::default()
            .with_title("GENtle egui window repro")
            .with_inner_size(DEFAULT_WINDOW_SIZE)
            .with_min_inner_size(MIN_WINDOW_SIZE),
        ..Default::default()
    };
    eframe::run_native(
        "GENtle egui window repro",
        options,
        Box::new(|cc| Ok(Box::new(ReproApp::new(cc)))),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cascade_position_uses_fixed_step() {
        assert_eq!(ReproApp::cascade_position(0), Pos2::new(64.0, 64.0));
        assert_eq!(ReproApp::cascade_position(3), Pos2::new(172.0, 172.0));
    }

    #[test]
    fn child_viewport_id_is_stable_for_serial() {
        assert_eq!(
            ReproApp::child_viewport_id(7),
            ReproApp::child_viewport_id(7)
        );
        assert_ne!(
            ReproApp::child_viewport_id(7),
            ReproApp::child_viewport_id(8)
        );
    }
}
