//! Operation-history window and undo/redo UI.
//!
//! This module is the first low-risk GUI extraction from the monolithic
//! `app.rs`: it keeps the history panel close to `GENtleApp` while consuming the
//! engine-owned history summary shared with shell/CLI adapters.

use super::*;

impl GENtleApp {
    pub(super) fn render_history_panel(&mut self, ctx: &egui::Context) {
        if !self.history_ui.show_panel {
            return;
        }
        let mut open = self.history_ui.show_panel;
        let (history_summary, history_rows, lock_failed) = match self.engine.read() {
            Ok(engine) => {
                let rows = engine
                    .operation_log()
                    .iter()
                    .rev()
                    .take(120)
                    .map(|record| {
                        (
                            record.result.op_id.clone(),
                            record.run_id.clone(),
                            Self::summarize_operation(&record.op),
                        )
                    })
                    .collect::<Vec<_>>();
                (engine.history_summary(), rows, false)
            }
            Err(_) => (
                crate::engine::EngineHistorySummary::default(),
                Vec::new(),
                true,
            ),
        };
        if lock_failed {
            self.app_status = "Operation History unavailable: engine lock poisoned".to_string();
        }

        let builder = egui::ViewportBuilder::default()
            .with_title("Operation History")
            .with_inner_size([820.0, 520.0])
            .with_min_inner_size([560.0, 320.0]);
        let viewport_id = Self::history_viewport_id();
        let render_history_in_foreground = self.viewport_foreground_requested(viewport_id);
        if ctx.embed_viewports() {
            let mut render_contents = |ui: &mut Ui| {
                ui.label("Operation-level history with multi-level undo/redo");
                ui.horizontal(|ui| {
                    if ui
                        .add_enabled(
                            history_summary.undo_count > 0,
                            egui::Button::new(match &history_summary.next_undo {
                                Some(next) => format!("Undo {}", next.operation),
                                None => "Undo".to_string(),
                            }),
                        )
                        .on_hover_text("Undo the most recent operation-level state transition")
                        .clicked()
                    {
                        self.undo_last_operation();
                    }
                    if ui
                        .add_enabled(
                            history_summary.redo_count > 0,
                            egui::Button::new(match &history_summary.next_redo {
                                Some(next) => format!("Redo {}", next.operation),
                                None => "Redo".to_string(),
                            }),
                        )
                        .on_hover_text("Redo the most recently undone operation-level transition")
                        .clicked()
                    {
                        self.redo_last_operation();
                    }
                    ui.small(format!(
                        "undo available: {} | redo available: {} | limit: {}",
                        history_summary.undo_count,
                        history_summary.redo_count,
                        history_summary.history_limit
                    ));
                });
                ui.small(format!(
                    "journal rows: {} | history is session-local and not saved in project files",
                    history_summary.operation_log_count
                ));
                ui.separator();
                if history_rows.is_empty() {
                    ui.small("No operations recorded yet.");
                } else {
                    egui::ScrollArea::vertical()
                        .id_salt(OPERATION_HISTORY_SCROLL_ID)
                        .show(ui, |ui| {
                            scroll_input_policy::apply_scrollarea_keyboard_navigation(
                                ui,
                                scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                            );
                            for (op_id, run_id, summary) in &history_rows {
                                ui.monospace(format!("[{op_id}] run={run_id}"));
                                ui.small(summary);
                                ui.separator();
                            }
                        });
                }
            };
            let spec = crate::egui_compat::HostedWindowSpec::new(
                "Operation History",
                egui::Id::new("Operation History"),
                Vec2::new(820.0, 520.0),
                Vec2::new(560.0, 320.0),
            )
            .foreground(render_history_in_foreground);
            crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| render_contents(ui));
            self.history_ui.show_panel = open;
            return;
        }
        ctx.show_viewport_immediate(viewport_id, builder, |ctx, class| {
            self.note_viewport_focus_if_active(ctx, viewport_id);
            let mut render_contents = |ui: &mut Ui| {
                ui.label("Operation-level history with multi-level undo/redo");
                ui.horizontal(|ui| {
                    if ui
                        .add_enabled(
                            history_summary.undo_count > 0,
                            egui::Button::new(match &history_summary.next_undo {
                                Some(next) => format!("Undo {}", next.operation),
                                None => "Undo".to_string(),
                            }),
                        )
                        .on_hover_text("Undo the most recent operation-level state transition")
                        .clicked()
                    {
                        self.undo_last_operation();
                    }
                    if ui
                        .add_enabled(
                            history_summary.redo_count > 0,
                            egui::Button::new(match &history_summary.next_redo {
                                Some(next) => format!("Redo {}", next.operation),
                                None => "Redo".to_string(),
                            }),
                        )
                        .on_hover_text("Redo the most recently undone operation-level transition")
                        .clicked()
                    {
                        self.redo_last_operation();
                    }
                    ui.small(format!(
                        "undo available: {} | redo available: {} | limit: {}",
                        history_summary.undo_count,
                        history_summary.redo_count,
                        history_summary.history_limit
                    ));
                });
                ui.small(format!(
                    "journal rows: {} | history is session-local and not saved in project files",
                    history_summary.operation_log_count
                ));
                ui.separator();
                if history_rows.is_empty() {
                    ui.small("No operations recorded yet.");
                } else {
                    egui::ScrollArea::vertical()
                        .id_salt(OPERATION_HISTORY_SCROLL_ID)
                        .show(ui, |ui| {
                            scroll_input_policy::apply_scrollarea_keyboard_navigation(
                                ui,
                                scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                            );
                            for (op_id, run_id, summary) in &history_rows {
                                ui.monospace(format!("[{op_id}] run={run_id}"));
                                ui.small(summary);
                                ui.separator();
                            }
                        });
                }
            };

            if class == egui::ViewportClass::EmbeddedWindow {
                let spec = crate::egui_compat::HostedWindowSpec::new(
                    "Operation History",
                    egui::Id::new("Operation History"),
                    Vec2::new(820.0, 520.0),
                    Vec2::new(560.0, 320.0),
                )
                .foreground(render_history_in_foreground);
                crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
                    render_contents(ui)
                });
            } else {
                crate::egui_compat::show_central_panel(ctx, egui::CentralPanel::default(), |ui| {
                    render_contents(ui);
                });
                if Self::viewport_close_requested_or_shortcut(ctx) {
                    open = false;
                }
            }
        });
        self.history_ui.show_panel = open;
    }
}
