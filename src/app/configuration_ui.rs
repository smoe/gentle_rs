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
            let agent_tab_label = self.tr("configuration.tab.agent_systems");
            let microarray_tab_label = self.tr("configuration.tab.microarrays");
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
                        self.configuration_tab == ConfigurationTab::AgentSystems,
                        agent_tab_label,
                    )
                    .clicked()
                {
                    self.configuration_tab = ConfigurationTab::AgentSystems;
                }
                if ui
                    .selectable_label(
                        self.configuration_tab == ConfigurationTab::Microarrays,
                        microarray_tab_label,
                    )
                    .clicked()
                {
                    self.configuration_tab = ConfigurationTab::Microarrays;
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
                    let tab_is_read_only_or_immediate = matches!(
                        self.configuration_tab,
                        ConfigurationTab::AgentSystems | ConfigurationTab::Microarrays
                    ) && !has_unapplied_changes;
                    let cancel_label = if tab_is_read_only_or_immediate {
                        self.tr("button.close")
                    } else {
                        self.tr("button.cancel")
                    };
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
                    if !tab_is_read_only_or_immediate
                        && ui
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
                        ui.with_layout(egui::Layout::top_down(egui::Align::Min), |ui| {
                            match self.configuration_tab {
                                ConfigurationTab::ExternalApplications => {
                                    self.render_configuration_external_tab(ui);
                                }
                                ConfigurationTab::AgentSystems => {
                                    self.render_agent_configuration_tab(ui);
                                }
                                ConfigurationTab::Microarrays => {
                                    self.render_configuration_microarrays_tab(ui);
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

    fn render_microarray_resource_row(ui: &mut Ui, label: &str, path: &Path, exists: bool) {
        ui.vertical(|ui| {
            let (status, color) = if exists {
                ("present", egui::Color32::from_rgb(20, 140, 45))
            } else {
                ("missing", egui::Color32::from_rgb(180, 83, 9))
            };
            ui.horizontal_wrapped(|ui| {
                ui.colored_label(color, status);
                ui.label(label);
                if ui.small_button("Copy path").clicked() {
                    ui.ctx().copy_text(path.display().to_string());
                }
            });
            ui.add(
                egui::Label::new(egui::RichText::new(path.display().to_string()).monospace())
                    .wrap(),
            );
        });
    }

    fn render_configuration_microarrays_tab(&mut self, ui: &mut Ui) {
        use crate::microarray_setup::{
            CLARIOM_D_HUMAN_SUPPORT_README, CLARIOM_D_PROBESET_URL, CLARIOM_D_TRANSCRIPT_URL,
            E_MTAB_14704_LANDING_URL, MicroarraySetupDiscovery,
        };

        let resources = MicroarraySetupDiscovery::discover_checkout();
        ui.heading(self.tr("configuration.microarray.heading"));
        ui.label(self.tr("configuration.microarray.description"));
        ui.small(
            "This page configures and explains reusable platform resources. Import, projection, and interpretation remain sequence-specific actions in Array setup.",
        );
        ui.add_space(8.0);

        egui::Frame::group(ui.style()).show(ui, |ui| {
            ui.strong("Choose the route that matches the files you actually have");
            ui.label(
                "1. Prepared GENtle output: select a directory containing region_intensity_chrom_order.csv, normalized_feature_matrix_manifest.json, and provenance.json, then inspect and project it.",
            );
            ui.label(
                "2. Completed APT/R tables: provide a normalized probeset summary matrix plus a NetAffx-style coordinate annotation table. Sample metadata and a PM-probe matrix are optional.",
            );
            ui.label(
                "3. Raw CEL files: obtain CEL files from your own experiment or a public repository, add sample metadata and platform annotation, then create and review an arrays probe-regions preflight before allowing R/APT execution.",
            );
        });

        ui.add_space(8.0);
        ui.heading("Expected sources and file roles");
        ui.label(
            "Raw experiment input: .CEL files plus IDF/SDRF or another sample table. The bundled TP73 example is E-MTAB-14704 from EMBL-EBI BioStudies/ArrayExpress.",
        );
        ui.hyperlink_to("Open E-MTAB-14704 at EMBL-EBI", E_MTAB_14704_LANDING_URL);
        ui.label(
            "Platform coordinates: the Clariom D Human na36 hg38 probeset and transcript annotation ZIPs are manual Thermo Fisher/NetAffx downloads. They are login-walled and GENtle never downloads them silently.",
        );
        ui.horizontal_wrapped(|ui| {
            ui.hyperlink_to("Probeset annotation ZIP", CLARIOM_D_PROBESET_URL);
            ui.hyperlink_to("Transcript annotation ZIP", CLARIOM_D_TRANSCRIPT_URL);
        });
        let support_readme = resources.checkout_root.join(CLARIOM_D_HUMAN_SUPPORT_README);
        Self::render_microarray_resource_row(
            ui,
            "local staging instructions",
            &support_readme,
            support_readme.is_file(),
        );

        ui.add_space(8.0);
        ui.heading("Detected in this checkout");
        ui.small(format!(
            "Selected bounded checkout root: {}. GENtle compares only the build checkout and current working checkout; it does not search unrelated folders.",
            resources.checkout_root.display()
        ));
        Self::render_microarray_resource_row(
            ui,
            resources.vendor_probeset_zip.role,
            &resources.vendor_probeset_zip.path,
            resources.vendor_probeset_zip.exists,
        );
        Self::render_microarray_resource_row(
            ui,
            resources.vendor_transcript_zip.role,
            &resources.vendor_transcript_zip.path,
            resources.vendor_transcript_zip.exists,
        );
        Self::render_microarray_resource_row(
            ui,
            resources.publication_idf.role,
            &resources.publication_idf.path,
            resources.publication_idf.exists,
        );
        Self::render_microarray_resource_row(
            ui,
            resources.publication_sdrf.role,
            &resources.publication_sdrf.path,
            resources.publication_sdrf.exists,
        );
        ui.label(format!(
            "E-MTAB-14704 raw CEL files present: {} (the published TP73 example expects 9)",
            resources.publication_cel_paths.len()
        ));
        Self::render_microarray_resource_row(
            ui,
            resources.local_summary_tsv.role,
            &resources.local_summary_tsv.path,
            resources.local_summary_tsv.exists,
        );
        Self::render_microarray_resource_row(
            ui,
            resources.local_annotation_tsv.role,
            &resources.local_annotation_tsv.path,
            resources.local_annotation_tsv.exists,
        );
        Self::render_microarray_resource_row(
            ui,
            resources.local_metadata_tsv.role,
            &resources.local_metadata_tsv.path,
            resources.local_metadata_tsv.exists,
        );
        if resources.has_local_analysis_tables() {
            ui.small(
                egui::RichText::new(
                    "These are full-platform analysis source tables, not a coordinate-complete direct-import subset. Use the raw-data probe-region preflight to select TP73 before helper output is projected.",
                )
                .color(egui::Color32::from_rgb(180, 83, 9)),
            );
        }
        ui.vertical(|ui| {
            let ready = resources.has_synthetic_output();
            let color = if ready {
                egui::Color32::from_rgb(20, 140, 45)
            } else {
                egui::Color32::from_rgb(180, 83, 9)
            };
            ui.horizontal_wrapped(|ui| {
                ui.colored_label(color, if ready { "present" } else { "incomplete" });
                ui.label(format!(
                    "synthetic TP73 demonstration output: {}/{} required files",
                    resources.synthetic_output_files_present,
                    resources.synthetic_output_required_file_count()
                ));
            });
            ui.add(
                egui::Label::new(
                    egui::RichText::new(resources.synthetic_output_dir.display().to_string())
                        .monospace(),
                )
                .wrap(),
            );
        });
        ui.small(
            egui::RichText::new(
                "The synthetic output is for learning and deterministic tests only; it is not biological evidence.",
            )
            .color(egui::Color32::from_rgb(180, 83, 9)),
        );
        ui.small(
            "Array setup in a DNA window reads these same checks, suggests a TP73 preflight for local CEL files, and offers only validated prepared output for direct selection; it never replaces paths you selected yourself.",
        );
    }

    pub(super) fn render_configuration_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_configuration_dialog {
            return;
        }
        let viewport_id = Self::configuration_viewport_id();
        let title = self.tr("dialog.configuration.title");
        let (preferred_size, min_size) = if matches!(
            self.configuration_tab,
            ConfigurationTab::AgentSystems | ConfigurationTab::Microarrays
        ) {
            (Vec2::new(980.0, 720.0), Vec2::new(620.0, 420.0))
        } else {
            (Vec2::new(720.0, 540.0), Vec2::new(460.0, 320.0))
        };
        let spec = self.hosted_window_spec_for_viewport(
            title.clone(),
            Self::hosted_configuration_window_id(),
            viewport_id,
            preferred_size,
            min_size,
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
            .with_inner_size(preferred_size)
            .with_min_inner_size(min_size);
        ctx.show_viewport_immediate(viewport_id, builder, |ctx, class| {
            self.note_viewport_focus_if_active(ctx, viewport_id);
            if class == egui::ViewportClass::EmbeddedWindow {
                let mut open = self.show_configuration_dialog;
                let render_started = Instant::now();
                let spec = crate::egui_compat::HostedWindowSpec::new(
                    title.clone(),
                    Self::hosted_configuration_window_id(),
                    preferred_size,
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
