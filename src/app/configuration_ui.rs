//! Configuration dialog shell and tab dispatcher.
//!
//! Detailed controls remain in `app.rs` for now; this module extracts the
//! window/content scaffold so the monolithic app update surface keeps shrinking
//! without changing the existing Apply/Cancel model.

use super::*;

impl GENtleApp {
    pub(super) fn render_configuration_contents(&mut self, ui: &mut Ui) {
        let has_unapplied_changes = self.configuration_has_unapplied_changes();
        window_backdrop::paint_window_backdrop(
            ui,
            WindowBackdropKind::Configuration,
            &self.window_backdrops,
        );
        with_window_content_inset(ui, |ui| {
            self.render_specialist_window_nav(ui);
            ui.horizontal_wrapped(|ui| {
                if ui
                    .selectable_label(
                        self.configuration_tab == ConfigurationTab::ExternalApplications,
                        "External Applications",
                    )
                    .clicked()
                {
                    self.configuration_tab = ConfigurationTab::ExternalApplications;
                }
                if ui
                    .selectable_label(
                        self.configuration_tab == ConfigurationTab::Graphics,
                        "Graphics",
                    )
                    .clicked()
                {
                    self.configuration_tab = ConfigurationTab::Graphics;
                }
                ui.separator();
                if has_unapplied_changes {
                    ui.colored_label(egui::Color32::from_rgb(185, 95, 25), "Unapplied changes");
                }
                if ui
                    .button("Close")
                    .on_hover_text(if has_unapplied_changes {
                        "Close configuration dialog (unapplied changes will be discarded)"
                    } else {
                        "Close configuration dialog"
                    })
                    .clicked()
                {
                    self.show_configuration_dialog = false;
                }
            });
            ui.separator();
            ui.with_layout(egui::Layout::bottom_up(egui::Align::Min), |ui| {
                ui.horizontal(|ui| {
                    if has_unapplied_changes {
                        ui.colored_label(egui::Color32::from_rgb(185, 95, 25), "Unapplied changes");
                    }
                    if ui
                        .button("Cancel")
                        .on_hover_text(if has_unapplied_changes {
                            "Discard unapplied configuration changes and close"
                        } else {
                            "Close configuration dialog"
                        })
                        .clicked()
                    {
                        self.sync_configuration_from_runtime();
                        self.configuration_status = if has_unapplied_changes {
                            "Discarded unapplied configuration changes".to_string()
                        } else {
                            "Closed configuration dialog".to_string()
                        };
                        self.show_configuration_dialog = false;
                    }
                    if ui
                        .add_enabled(has_unapplied_changes, egui::Button::new("Apply"))
                        .on_hover_text("Apply all unapplied configuration changes")
                        .clicked()
                    {
                        self.apply_pending_configuration_changes();
                    }
                });
                if !self.configuration_status.trim().is_empty() {
                    ui.separator();
                    ui.monospace(self.configuration_status.clone());
                }
                ui.separator();
                egui::ScrollArea::vertical()
                    .id_salt("configuration_dialog_scroll")
                    .auto_shrink([false, false])
                    .show(ui, |ui| {
                        scroll_input_policy::apply_scrollarea_keyboard_navigation(
                            ui,
                            scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                        );
                        match self.configuration_tab {
                            ConfigurationTab::ExternalApplications => {
                                self.render_configuration_external_tab(ui);
                            }
                            ConfigurationTab::Graphics => {
                                self.render_configuration_graphics_tab(ui);
                            }
                        }
                    });
            });
        });
    }

    pub(super) fn render_configuration_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_configuration_dialog {
            return;
        }
        let viewport_id = Self::configuration_viewport_id();
        let spec = self.hosted_window_spec_for_viewport(
            "Configuration",
            Self::hosted_configuration_window_id(),
            viewport_id,
            Vec2::new(720.0, 540.0),
            Vec2::new(460.0, 320.0),
        );
        if ctx.embed_viewports() {
            let mut open = self.show_configuration_dialog;
            let render_started = Instant::now();
            crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
                self.render_configuration_contents(ui)
            });
            self.clear_viewport_foreground_request_after_render(viewport_id);
            self.note_slow_open_phase(
                viewport_id,
                "Configuration first-frame render",
                render_started.elapsed().as_millis(),
            );
            self.show_configuration_dialog =
                Self::reconcile_embedded_window_open_state(self.show_configuration_dialog, open);
            self.finalize_viewport_open_probe(viewport_id, "Configuration");
            return;
        }
        let builder = egui::ViewportBuilder::default()
            .with_title("Configuration")
            .with_inner_size([720.0, 540.0])
            .with_min_inner_size([460.0, 320.0]);
        ctx.show_viewport_immediate(viewport_id, builder, |ctx, class| {
            self.note_viewport_focus_if_active(ctx, viewport_id);
            if class == egui::ViewportClass::EmbeddedWindow {
                let mut open = self.show_configuration_dialog;
                let render_started = Instant::now();
                let min_size = Vec2::new(460.0, 320.0);
                let spec = crate::egui_compat::HostedWindowSpec::new(
                    "Configuration",
                    Self::hosted_configuration_window_id(),
                    Vec2::new(720.0, 540.0),
                    min_size,
                )
                .foreground(
                    self.pending_focus_viewports.contains(&viewport_id)
                        || self
                            .pending_viewport_focus_timestamps
                            .contains_key(&viewport_id),
                );
                crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
                    self.render_configuration_contents(ui)
                });
                self.note_slow_open_phase(
                    viewport_id,
                    "Configuration first-frame render",
                    render_started.elapsed().as_millis(),
                );
                self.show_configuration_dialog = Self::reconcile_embedded_window_open_state(
                    self.show_configuration_dialog,
                    open,
                );
                return;
            }

            let render_started = Instant::now();
            crate::egui_compat::show_central_panel(
                ctx,
                egui::CentralPanel::default().frame(egui::Frame::NONE),
                |ui| {
                    self.render_configuration_contents(ui);
                },
            );
            self.note_slow_open_phase(
                viewport_id,
                "Configuration first-frame render",
                render_started.elapsed().as_millis(),
            );

            if Self::viewport_close_requested_or_shortcut(ctx) {
                self.show_configuration_dialog = false;
            }
        });
        self.finalize_viewport_open_probe(viewport_id, "Configuration");
    }
}
