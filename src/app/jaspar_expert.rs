//! App-level JASPAR expert window.
//!
//! This module keeps the JASPAR registry browser and expert detail surface out
//! of the already-large `app.rs` file while still reusing the shared
//! engine-owned expert payload.
//!
//! Look here for:
//! - the app-level `JASPAR Expert` window
//! - motif filtering/selection helpers
//! - simple in-window sequence-logo and histogram painters

use super::*;
use crate::engine::{JasparCatalogRow, JasparScoreDistributionPanel};

impl GENtleApp {
    fn jaspar_catalog_rows_local_fallback(&self) -> Vec<JasparCatalogRow> {
        tf_motifs::list_motif_summaries()
            .into_iter()
            .map(|row| JasparCatalogRow {
                motif_id: row.id,
                motif_name: row.name,
                consensus_iupac: row.consensus_iupac,
                motif_length_bp: row.motif_length_bp,
                remote_summary: None,
            })
            .collect()
    }

    fn filtered_jaspar_catalog_rows(&self) -> Vec<JasparCatalogRow> {
        let filter = self.jaspar_expert_filter.trim().to_ascii_uppercase();
        self.jaspar_catalog_report
            .as_ref()
            .map(|report| report.rows.clone())
            .unwrap_or_else(|| self.jaspar_catalog_rows_local_fallback())
            .into_iter()
            .filter(|row| {
                if filter.is_empty() {
                    return true;
                }
                row.motif_id.to_ascii_uppercase().contains(&filter)
                    || row
                        .motif_name
                        .as_deref()
                        .unwrap_or("")
                        .to_ascii_uppercase()
                        .contains(&filter)
                    || row.consensus_iupac.to_ascii_uppercase().contains(&filter)
            })
            .collect()
    }

    fn merge_jaspar_catalog_report(&mut self, report: JasparCatalogReport) {
        if !report.include_remote_metadata {
            self.jaspar_catalog_report = Some(report);
            return;
        }
        let registry_entry_count = report.registry_entry_count;
        let existing_rows = self
            .jaspar_catalog_report
            .as_ref()
            .map(|report| report.rows.clone())
            .unwrap_or_else(|| self.jaspar_catalog_rows_local_fallback());
        let mut by_id = existing_rows
            .into_iter()
            .map(|row| (row.motif_id.clone(), row))
            .collect::<BTreeMap<_, _>>();
        for row in report.rows.clone() {
            by_id.insert(row.motif_id.clone(), row);
        }
        let rows = by_id.into_values().collect::<Vec<_>>();
        self.jaspar_catalog_report = Some(JasparCatalogReport {
            rows: rows.clone(),
            returned_entry_count: rows.len(),
            registry_entry_count: registry_entry_count.max(rows.len()),
            ..report
        });
    }

    fn refresh_jaspar_catalog(&mut self, include_remote_metadata: bool) {
        let result = self.engine.write().expect("Engine lock poisoned").apply(
            Operation::ListJasparCatalog {
                filter: if include_remote_metadata {
                    Some(self.jaspar_expert_filter.trim().to_string())
                } else {
                    None
                },
                limit: None,
                include_remote_metadata,
                path: None,
            },
        );
        match result {
            Ok(result) => {
                if let Some(report) = result.jaspar_catalog_report {
                    let fetched_count = report.returned_entry_count;
                    let warning_count = report.warnings.len();
                    self.merge_jaspar_catalog_report(report);
                    self.jaspar_expert_status = if include_remote_metadata {
                        if warning_count == 0 {
                            format!(
                                "Fetched remote JASPAR metadata for {} visible catalog entr{}.",
                                fetched_count,
                                if fetched_count == 1 { "y" } else { "ies" }
                            )
                        } else {
                            format!(
                                "Fetched remote JASPAR metadata for {} visible catalog entr{} with {} warning(s).",
                                fetched_count,
                                if fetched_count == 1 { "y" } else { "ies" },
                                warning_count
                            )
                        }
                    } else {
                        format!(
                            "Loaded local JASPAR catalog with {} entries.",
                            fetched_count
                        )
                    };
                }
            }
            Err(err) => {
                self.jaspar_expert_status = err.to_string();
            }
        }
    }

    fn refresh_jaspar_expert_view(&mut self) {
        let motif = self.jaspar_expert_selected_motif_id.trim().to_string();
        if motif.is_empty() {
            self.jaspar_expert_status = "Choose one JASPAR entry first.".to_string();
            self.jaspar_expert_view = None;
            return;
        }
        let random_sequence_length_bp = match self.jaspar_expert_random_length_bp.trim().parse() {
            Ok(value) => value,
            Err(e) => {
                self.jaspar_expert_status = format!(
                    "Invalid random background length '{}': {e}",
                    self.jaspar_expert_random_length_bp.trim()
                );
                return;
            }
        };
        let random_seed = match self.jaspar_expert_random_seed.trim().parse() {
            Ok(value) => value,
            Err(e) => {
                self.jaspar_expert_status = format!(
                    "Invalid random seed '{}': {e}",
                    self.jaspar_expert_random_seed.trim()
                );
                return;
            }
        };

        let result = self.engine.write().expect("Engine lock poisoned").apply(
            Operation::InspectJasparEntry {
                motif: motif.clone(),
                random_sequence_length_bp,
                random_seed,
                include_remote_metadata: self.jaspar_expert_fetch_remote_metadata,
                path: None,
            },
        );
        match result {
            Ok(result) => {
                self.jaspar_expert_view = result.jaspar_entry_expert_view;
                if let Some(view) = self.jaspar_expert_view.as_ref() {
                    if let Some(remote_metadata) = view.remote_metadata.clone() {
                        let mut row = self
                            .jaspar_catalog_report
                            .as_ref()
                            .and_then(|report| {
                                report
                                    .rows
                                    .iter()
                                    .find(|row| row.motif_id == view.motif_id)
                                    .cloned()
                            })
                            .unwrap_or(JasparCatalogRow {
                                motif_id: view.motif_id.clone(),
                                motif_name: view.motif_name.clone(),
                                consensus_iupac: view.consensus_iupac.clone(),
                                motif_length_bp: view.motif_length_bp,
                                remote_summary: None,
                            });
                        row.remote_summary = Some(GentleEngine::jaspar_remote_metadata_summary(
                            &remote_metadata,
                        ));
                        let base =
                            self.jaspar_catalog_report
                                .as_ref()
                                .cloned()
                                .unwrap_or_else(|| JasparCatalogReport {
                                    schema: "gentle.jaspar_catalog.v1".to_string(),
                                    generated_at_unix_ms: 0,
                                    op_id: None,
                                    run_id: None,
                                    filter: None,
                                    limit: None,
                                    include_remote_metadata: true,
                                    registry_entry_count: tf_motifs::list_motif_summaries().len(),
                                    returned_entry_count: 0,
                                    rows: self.jaspar_catalog_rows_local_fallback(),
                                    warnings: vec![],
                                });
                        self.merge_jaspar_catalog_report(JasparCatalogReport {
                            rows: vec![row],
                            returned_entry_count: 1,
                            ..base
                        });
                    }
                }
                let mut status_parts = vec![];
                if !result.messages.is_empty() {
                    status_parts.push(result.messages.join(" "));
                }
                if !result.warnings.is_empty() {
                    status_parts.push(format!("warnings: {}", result.warnings.join(" | ")));
                }
                self.jaspar_expert_status = if status_parts.is_empty() {
                    format!("Loaded JASPAR expert view for '{}'.", motif)
                } else {
                    status_parts.join("  ")
                };
            }
            Err(err) => {
                self.jaspar_expert_view = None;
                self.jaspar_expert_status = err.to_string();
            }
        }
    }

    fn jaspar_logo_base_color(base: char) -> egui::Color32 {
        match base {
            'A' => egui::Color32::from_rgb(34, 139, 34),
            'C' => egui::Color32::from_rgb(30, 144, 255),
            'G' => egui::Color32::from_rgb(255, 140, 0),
            'T' => egui::Color32::from_rgb(220, 20, 60),
            _ => egui::Color32::GRAY,
        }
    }

    fn render_jaspar_logo(ui: &mut Ui, view: &JasparEntryExpertView) {
        let column_width = 28.0_f32;
        let width = (view.columns.len() as f32 * column_width + 60.0).max(280.0);
        let desired = egui::Vec2::new(width, 220.0);
        let (rect, _) = ui.allocate_exact_size(desired, egui::Sense::hover());
        let painter = ui.painter_at(rect);
        let margin = 32.0_f32;
        let baseline = rect.bottom() - margin;
        let left = rect.left() + margin;
        let usable_height = (rect.height() - margin * 1.8).max(1.0);
        let bits_to_px = usable_height / 2.0;

        painter.rect_stroke(
            rect,
            4.0,
            egui::Stroke::new(1.0, egui::Color32::from_gray(90)),
            egui::StrokeKind::Inside,
        );
        painter.line_segment(
            [
                egui::pos2(left, baseline),
                egui::pos2(rect.right() - 8.0, baseline),
            ],
            egui::Stroke::new(1.0, egui::Color32::LIGHT_GRAY),
        );
        painter.line_segment(
            [
                egui::pos2(left, baseline),
                egui::pos2(left, rect.top() + 8.0),
            ],
            egui::Stroke::new(1.0, egui::Color32::LIGHT_GRAY),
        );
        painter.text(
            egui::pos2(rect.left() + 10.0, baseline),
            egui::Align2::LEFT_BOTTOM,
            "0",
            egui::FontId::proportional(11.0),
            egui::Color32::LIGHT_GRAY,
        );
        painter.text(
            egui::pos2(rect.left() + 10.0, baseline - usable_height),
            egui::Align2::LEFT_CENTER,
            "2 bits",
            egui::FontId::proportional(11.0),
            egui::Color32::LIGHT_GRAY,
        );

        for column in &view.columns {
            let x = left + (column.position_1based as f32 - 0.5) * column_width;
            let mut rows = vec![
                ('A', column.a_logo_bits),
                ('C', column.c_logo_bits),
                ('G', column.g_logo_bits),
                ('T', column.t_logo_bits),
            ];
            rows.sort_by(|left, right| left.1.total_cmp(&right.1));
            let mut used_height = 0.0_f32;
            for (base, bits) in rows {
                if bits <= 0.0001 {
                    continue;
                }
                let height_px = (bits as f32 * bits_to_px).max(10.0);
                painter.text(
                    egui::pos2(x, baseline - used_height),
                    egui::Align2::CENTER_BOTTOM,
                    base.to_string(),
                    egui::FontId::monospace(height_px),
                    Self::jaspar_logo_base_color(base),
                );
                used_height += height_px * 0.86;
            }
            painter.text(
                egui::pos2(x, baseline + 10.0),
                egui::Align2::CENTER_TOP,
                column.position_1based.to_string(),
                egui::FontId::proportional(10.0),
                egui::Color32::LIGHT_GRAY,
            );
        }
    }

    fn render_jaspar_distribution_histogram(ui: &mut Ui, panel: &JasparScoreDistributionPanel) {
        let desired = egui::Vec2::new(ui.available_width().max(360.0), 170.0);
        let (rect, _) = ui.allocate_exact_size(desired, egui::Sense::hover());
        let painter = ui.painter_at(rect);
        let max_count = panel
            .histogram_bins
            .iter()
            .map(|row| row.count)
            .max()
            .unwrap_or(1)
            .max(1) as f32;
        let margin = 28.0_f32;
        let plot_rect = egui::Rect::from_min_max(
            egui::pos2(rect.left() + margin, rect.top() + 8.0),
            egui::pos2(rect.right() - 10.0, rect.bottom() - margin),
        );
        painter.rect_stroke(
            rect,
            4.0,
            egui::Stroke::new(1.0, egui::Color32::from_gray(90)),
            egui::StrokeKind::Inside,
        );
        painter.line_segment(
            [
                egui::pos2(plot_rect.left(), plot_rect.bottom()),
                egui::pos2(plot_rect.right(), plot_rect.bottom()),
            ],
            egui::Stroke::new(1.0, egui::Color32::LIGHT_GRAY),
        );
        painter.line_segment(
            [
                egui::pos2(plot_rect.left(), plot_rect.top()),
                egui::pos2(plot_rect.left(), plot_rect.bottom()),
            ],
            egui::Stroke::new(1.0, egui::Color32::LIGHT_GRAY),
        );
        if !panel.histogram_bins.is_empty() {
            let bar_width = plot_rect.width() / panel.histogram_bins.len() as f32;
            for (idx, bin) in panel.histogram_bins.iter().enumerate() {
                let height = if max_count <= 0.0 {
                    0.0
                } else {
                    plot_rect.height() * (bin.count as f32 / max_count)
                };
                let bar = egui::Rect::from_min_max(
                    egui::pos2(
                        plot_rect.left() + idx as f32 * bar_width,
                        plot_rect.bottom() - height,
                    ),
                    egui::pos2(
                        plot_rect.left() + (idx as f32 + 1.0) * bar_width - 1.0,
                        plot_rect.bottom(),
                    ),
                );
                painter.rect_filled(bar, 0.0, egui::Color32::from_rgb(88, 129, 255));
            }
            if let Some(first) = panel.histogram_bins.first() {
                painter.text(
                    egui::pos2(plot_rect.left(), rect.bottom() - 8.0),
                    egui::Align2::LEFT_BOTTOM,
                    format!("{:.2}", first.start_score),
                    egui::FontId::proportional(10.0),
                    egui::Color32::LIGHT_GRAY,
                );
            }
            if let Some(last) = panel.histogram_bins.last() {
                painter.text(
                    egui::pos2(plot_rect.right(), rect.bottom() - 8.0),
                    egui::Align2::RIGHT_BOTTOM,
                    format!("{:.2}", last.end_score),
                    egui::FontId::proportional(10.0),
                    egui::Color32::LIGHT_GRAY,
                );
            }
        }
        painter.text(
            egui::pos2(rect.left() + 6.0, plot_rect.top()),
            egui::Align2::LEFT_TOP,
            format!("max {}", max_count as usize),
            egui::FontId::proportional(10.0),
            egui::Color32::LIGHT_GRAY,
        );
    }

    pub(super) fn render_jaspar_expert_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_jaspar_expert_dialog {
            return;
        }

        if self.jaspar_catalog_report.is_none() {
            self.refresh_jaspar_catalog(false);
        }
        let entries = self.filtered_jaspar_catalog_rows();
        if self.jaspar_expert_selected_motif_id.trim().is_empty() {
            self.jaspar_expert_selected_motif_id = entries
                .first()
                .map(|row| row.motif_id.clone())
                .unwrap_or_default();
        }
        let selected_still_visible = entries
            .iter()
            .any(|row| row.motif_id == self.jaspar_expert_selected_motif_id);
        if !selected_still_visible {
            self.jaspar_expert_selected_motif_id = entries
                .first()
                .map(|row| row.motif_id.clone())
                .unwrap_or_default();
        }
        if self.jaspar_expert_view.is_none() && !self.jaspar_expert_selected_motif_id.is_empty() {
            self.refresh_jaspar_expert_view();
        }

        let mut open = self.show_jaspar_expert_dialog;
        egui::Window::new("JASPAR Expert")
            .open(&mut open)
            .resizable(true)
            .default_size(egui::Vec2::new(1280.0, 900.0))
            .min_size(egui::Vec2::new(920.0, 640.0))
            .show(ctx, |ui| {
                ui.label("Inspect local JASPAR entries through GENtle’s own matrix/scoring path, with optional remote species metadata from the JASPAR REST API.");
                ui.horizontal(|ui| {
                    ui.label("Filter");
                    ui.text_edit_singleline(&mut self.jaspar_expert_filter);
                    ui.separator();
                    ui.label("Random background bp");
                    ui.add(
                        egui::TextEdit::singleline(&mut self.jaspar_expert_random_length_bp)
                            .desired_width(90.0),
                    );
                    ui.label("Seed");
                    ui.add(
                        egui::TextEdit::singleline(&mut self.jaspar_expert_random_seed)
                            .desired_width(120.0),
                    );
                    ui.checkbox(
                        &mut self.jaspar_expert_fetch_remote_metadata,
                        "Fetch JASPAR species metadata",
                    );
                    if ui.button("Reload catalog").clicked() {
                        self.refresh_jaspar_catalog(false);
                    }
                    if ui
                        .button("Fetch visible species")
                        .on_hover_text(
                            "Request remote JASPAR metadata for the currently visible filtered subset and merge species summaries into the catalog table.",
                        )
                        .clicked()
                    {
                        self.refresh_jaspar_catalog(true);
                    }
                    if ui.button("Inspect selected").clicked() {
                        self.refresh_jaspar_expert_view();
                    }
                });
                if !self.jaspar_expert_status.trim().is_empty() {
                    ui.small(self.jaspar_expert_status.clone());
                }
                ui.separator();
                ui.columns(2, |columns| {
                    columns[0].vertical(|ui| {
                        let total_entries = self
                            .jaspar_catalog_report
                            .as_ref()
                            .map(|report| report.registry_entry_count)
                            .unwrap_or(entries.len());
                        ui.heading(format!(
                            "Available entries (showing {} / {})",
                            entries.len(),
                            total_entries
                        ));
                        egui::ScrollArea::both()
                            .auto_shrink([false, false])
                            .max_height(720.0)
                            .show(ui, |ui| {
                                egui::Grid::new("jaspar_catalog_grid")
                                    .striped(true)
                                    .min_col_width(52.0)
                                    .show(ui, |ui| {
                                        ui.strong("ID");
                                        ui.strong("Name");
                                        ui.strong("Len");
                                        ui.strong("Consensus");
                                        ui.strong("Species");
                                        ui.strong("Collection");
                                        ui.end_row();
                                        for row in &entries {
                                            let selected = self.jaspar_expert_selected_motif_id
                                                == row.motif_id;
                                            if ui
                                                .selectable_label(selected, row.motif_id.clone())
                                                .clicked()
                                            {
                                                self.jaspar_expert_selected_motif_id =
                                                    row.motif_id.clone();
                                                self.refresh_jaspar_expert_view();
                                            }
                                            ui.label(
                                                row.motif_name
                                                    .as_deref()
                                                    .filter(|value| !value.trim().is_empty())
                                                    .unwrap_or(""),
                                            );
                                            ui.label(row.motif_length_bp.to_string());
                                            ui.monospace(row.consensus_iupac.clone());
                                            ui.label(
                                                row.remote_summary
                                                    .as_ref()
                                                    .map(|remote| {
                                                        if remote.species_count == 0 {
                                                            "0".to_string()
                                                        } else {
                                                            remote.species_count.to_string()
                                                        }
                                                    })
                                                    .unwrap_or_else(|| "—".to_string()),
                                            );
                                            ui.label(
                                                row.remote_summary
                                                    .as_ref()
                                                    .and_then(|remote| remote.collection.as_deref())
                                                    .unwrap_or(""),
                                            );
                                            ui.end_row();
                                        }
                                    });
                                if let Some(report) = self.jaspar_catalog_report.as_ref() {
                                    if !report.warnings.is_empty() {
                                        ui.separator();
                                        ui.small(format!(
                                            "Catalog warnings: {}",
                                            report.warnings.join(" | ")
                                        ));
                                    }
                                }
                            });
                    });
                    columns[1].vertical(|ui| {
                        let Some(view) = self.jaspar_expert_view.as_ref() else {
                            ui.colored_label(
                                egui::Color32::from_rgb(190, 70, 70),
                                "No JASPAR expert view is loaded yet.",
                            );
                            return;
                        };
                        ui.heading(format!(
                            "{}{}",
                            view.motif_id,
                            view.motif_name
                                .as_deref()
                                .map(|name| format!(" — {name}"))
                                .unwrap_or_default()
                        ));
                        ui.horizontal_wrapped(|ui| {
                            ui.label(format!("consensus: {}", view.consensus_iupac));
                            ui.separator();
                            ui.label(format!("length: {} bp", view.motif_length_bp));
                            ui.separator();
                            ui.label(format!("registry entries: {}", view.registry_entry_count));
                        });
                        if let Some(remote) = view.remote_metadata.as_ref() {
                            ui.separator();
                            ui.heading("Remote metadata");
                            ui.horizontal_wrapped(|ui| {
                                if let Some(collection) = remote.collection.as_deref() {
                                    ui.label(format!("collection: {collection}"));
                                }
                                if let Some(tax_group) = remote.tax_group.as_deref() {
                                    ui.label(format!("tax group: {tax_group}"));
                                }
                                if let Some(tf_class) = remote.tf_class.as_deref() {
                                    ui.label(format!("class: {tf_class}"));
                                }
                                if let Some(tf_family) = remote.tf_family.as_deref() {
                                    ui.label(format!("family: {tf_family}"));
                                }
                                if let Some(data_type) = remote.data_type.as_deref() {
                                    ui.label(format!("data type: {data_type}"));
                                }
                            });
                            if !remote.species_assignments.is_empty() {
                                ui.label("Species assignments:");
                                for species in &remote.species_assignments {
                                    let mut line = species.scientific_name.clone();
                                    if let Some(tax_id) = species.tax_id.as_deref() {
                                        line.push_str(&format!(" (tax_id={tax_id})"));
                                    }
                                    if let Some(common_name) = species.common_name.as_deref() {
                                        line.push_str(&format!(" — {common_name}"));
                                    }
                                    ui.small(line);
                                }
                            }
                            ui.small(format!("source: {}", remote.source_url));
                        } else if view.include_remote_metadata {
                            ui.small("Remote metadata was requested but no species/class metadata could be attached.");
                        }
                        if !view.warnings.is_empty() {
                            ui.separator();
                            ui.colored_label(
                                egui::Color32::from_rgb(200, 140, 60),
                                format!("Warnings: {}", view.warnings.join(" | ")),
                            );
                        }
                        ui.separator();
                        ui.heading("Sequence logo");
                        Self::render_jaspar_logo(ui, view);
                        ui.separator();
                        ui.heading("Count matrix");
                        egui::ScrollArea::both()
                            .auto_shrink([false, false])
                            .max_height(220.0)
                            .show(ui, |ui| {
                                egui::Grid::new("jaspar_matrix_grid")
                                    .striped(true)
                                    .show(ui, |ui| {
                                        ui.strong("Pos");
                                        ui.strong("A");
                                        ui.strong("C");
                                        ui.strong("G");
                                        ui.strong("T");
                                        ui.strong("Total");
                                        ui.strong("IC bits");
                                        ui.strong("Dominant");
                                        ui.end_row();
                                        for column in &view.columns {
                                            ui.label(column.position_1based.to_string());
                                            ui.label(format!("{:.1}", column.a_count));
                                            ui.label(format!("{:.1}", column.c_count));
                                            ui.label(format!("{:.1}", column.g_count));
                                            ui.label(format!("{:.1}", column.t_count));
                                            ui.label(format!("{:.1}", column.total_count));
                                            ui.label(format!("{:.2}", column.information_content_bits));
                                            ui.label(column.dominant_base.clone());
                                            ui.end_row();
                                        }
                                    });
                            });
                        ui.separator();
                        ui.heading("Score distributions");
                        for panel in &view.score_panels {
                            egui::CollapsingHeader::new(panel.label.clone())
                                .default_open(true)
                                .show(ui, |ui| {
                                    ui.horizontal_wrapped(|ui| {
                                        ui.label(format!(
                                            "max: {} ({:.2}, q={:.3})",
                                            panel.maximizing_sequence,
                                            panel.maximizing_score,
                                            panel.maximizing_quantile
                                        ));
                                        ui.separator();
                                        ui.label(format!(
                                            "min: {} ({:.2}, q={:.3})",
                                            panel.minimizing_sequence,
                                            panel.minimizing_score,
                                            panel.minimizing_quantile
                                        ));
                                    });
                                    ui.small(format!(
                                        "mean={:.2} stddev={:.2} p05={:.2} p50={:.2} p95={:.2} positive_fraction={:.3}",
                                        panel.distribution.mean_score,
                                        panel.distribution.stddev_score,
                                        panel.distribution.p05_score,
                                        panel.distribution.p50_score,
                                        panel.distribution.p95_score,
                                        panel.distribution.positive_fraction
                                    ));
                                    Self::render_jaspar_distribution_histogram(ui, panel);
                                });
                        }
                    });
                });
            });
        self.show_jaspar_expert_dialog = open;
    }
}

#[cfg(test)]
mod tests {
    use super::GENtleApp;

    #[test]
    fn open_jaspar_expert_dialog_seeds_default_selection() {
        let mut app = GENtleApp::default();
        app.jaspar_expert_selected_motif_id.clear();
        app.open_jaspar_expert_dialog();
        assert!(app.show_jaspar_expert_dialog);
        assert!(!app.jaspar_expert_selected_motif_id.trim().is_empty());
    }

    #[test]
    fn refresh_jaspar_expert_view_loads_selected_entry() {
        let mut app = GENtleApp::default();
        app.jaspar_expert_selected_motif_id = "SP1".to_string();
        app.jaspar_expert_random_length_bp = "512".to_string();
        app.jaspar_expert_random_seed = "7".to_string();
        app.jaspar_expert_fetch_remote_metadata = false;
        app.refresh_jaspar_expert_view();
        assert_eq!(
            app.jaspar_expert_view
                .as_ref()
                .and_then(|view| view.motif_name.as_deref()),
            Some("SP1")
        );
    }

    #[test]
    fn refresh_jaspar_catalog_loads_local_rows() {
        let mut app = GENtleApp::default();
        app.refresh_jaspar_catalog(false);
        let report = app
            .jaspar_catalog_report
            .as_ref()
            .expect("catalog report should load");
        assert!(report.returned_entry_count > 0);
        assert!(
            report
                .rows
                .iter()
                .any(|row| row.motif_name.as_deref() == Some("SP1"))
        );
    }
}
