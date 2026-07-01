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
            let external_tab_label = self.tr("configuration.tab.external_applications");
            let graphics_tab_label = self.tr("configuration.tab.graphics");
            let language_tab_label = self.tr("configuration.tab.language");
            let unapplied_changes_label = self.tr("configuration.status.unapplied_changes");
            let close_label = self.tr("button.close");
            ui.horizontal_wrapped(|ui| {
                if ui
                    .selectable_label(
                        self.configuration_tab == ConfigurationTab::ExternalApplications,
                        external_tab_label,
                    )
                    .clicked()
                {
                    self.configuration_tab = ConfigurationTab::ExternalApplications;
                }
                if ui
                    .selectable_label(
                        self.configuration_tab == ConfigurationTab::Graphics,
                        graphics_tab_label,
                    )
                    .clicked()
                {
                    self.configuration_tab = ConfigurationTab::Graphics;
                }
                if ui
                    .selectable_label(
                        self.configuration_tab == ConfigurationTab::Language,
                        language_tab_label,
                    )
                    .clicked()
                {
                    self.configuration_tab = ConfigurationTab::Language;
                }
                ui.separator();
                if has_unapplied_changes {
                    ui.colored_label(
                        egui::Color32::from_rgb(185, 95, 25),
                        unapplied_changes_label.clone(),
                    );
                }
                if ui
                    .button(close_label)
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
                    let cancel_label = self.tr("button.cancel");
                    let apply_label = self.tr("button.apply");
                    if has_unapplied_changes {
                        ui.colored_label(
                            egui::Color32::from_rgb(185, 95, 25),
                            self.tr("configuration.status.unapplied_changes"),
                        );
                    }
                    if ui
                        .button(cancel_label)
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
                        .add_enabled(has_unapplied_changes, egui::Button::new(apply_label))
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
                            ConfigurationTab::Language => {
                                self.render_configuration_language_tab(ui);
                            }
                        }
                    });
            });
        });
    }

    fn render_configuration_language_tab(&mut self, ui: &mut Ui) {
        ui.heading(self.tr("configuration.language.heading"));
        ui.label(self.tr("configuration.language.description"));
        ui.add_space(8.0);

        ui.horizontal(|ui| {
            ui.label(self.tr("configuration.language.selector"));
            let before = self.configuration_ui_language;
            egui::ComboBox::from_id_salt("configuration_language_selector")
                .selected_text(self.configuration_ui_language.label())
                .show_ui(ui, |ui| {
                    for language in UiLanguage::ALL {
                        ui.selectable_value(
                            &mut self.configuration_ui_language,
                            language,
                            language.label(),
                        );
                    }
                });
            if self.configuration_ui_language != before {
                self.configuration_language_dirty =
                    self.configuration_ui_language != self.ui_language;
            }
        });

        ui.small(format!(
            "{}: {} | {}: {}",
            self.tr("configuration.language.active"),
            self.i18n.language().label(),
            self.tr("configuration.language.selected"),
            self.configuration_ui_language.label()
        ));
        ui.small(self.tr("configuration.language.apply_note"));
        ui.small(self.tr("configuration.language.experimental_note"));
    }

    pub(super) fn render_configuration_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_configuration_dialog {
            return;
        }
        let viewport_id = Self::configuration_viewport_id();
        let title = self.tr("dialog.configuration.title");
        let spec = self.hosted_window_spec_for_viewport(
            title.clone(),
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
            .with_title(title.clone())
            .with_inner_size([720.0, 540.0])
            .with_min_inner_size([460.0, 320.0]);
        ctx.show_viewport_immediate(viewport_id, builder, |ctx, class| {
            self.note_viewport_focus_if_active(ctx, viewport_id);
            if class == egui::ViewportClass::EmbeddedWindow {
                let mut open = self.show_configuration_dialog;
                let render_started = Instant::now();
                let min_size = Vec2::new(460.0, 320.0);
                let spec = crate::egui_compat::HostedWindowSpec::new(
                    title.clone(),
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
                crate::egui_compat::show_hosted_window(&mut *ctx, &spec, &mut open, |ui| {
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
                &mut *ctx,
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
