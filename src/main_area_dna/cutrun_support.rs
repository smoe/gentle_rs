//! CUT&RUN regulatory-support inspection helpers for `MainAreaDna`.
//!
//! This module keeps the DNA-window CUT&RUN UI as a thin presentation layer over
//! the shared `InspectCutRunRegulatorySupport` engine report. It deliberately
//! does not add GUI-only scoring or biological interpretation logic.

use super::*;

impl MainAreaDna {
    pub(super) fn cutrun_evidence_source_kind_label(
        kind: CutRunRegulatoryEvidenceSourceKind,
    ) -> &'static str {
        match kind {
            CutRunRegulatoryEvidenceSourceKind::Dataset => "dataset",
            CutRunRegulatoryEvidenceSourceKind::ReadReport => "read report",
        }
    }

    pub(super) fn cutrun_support_strength_label(strength: CutRunSupportStrength) -> &'static str {
        match strength {
            CutRunSupportStrength::Weak => "weak",
            CutRunSupportStrength::Moderate => "moderate",
            CutRunSupportStrength::Strong => "strong",
        }
    }

    pub(super) fn cutrun_tfbs_confirmation_label(
        status: CutRunRegulatoryTfbsConfirmationStatus,
    ) -> &'static str {
        match status {
            CutRunRegulatoryTfbsConfirmationStatus::Confirmed => "confirmed",
            CutRunRegulatoryTfbsConfirmationStatus::Unconfirmed => "unconfirmed",
        }
    }

    pub(super) fn cutrun_motif_context_scope_label(scope: CutRunMotifContextScope) -> &'static str {
        match scope {
            CutRunMotifContextScope::InsideWindow => "inside",
            CutRunMotifContextScope::NeighborWindow => "nearby",
        }
    }

    pub(super) fn cutrun_occupancy_interpretation_label(
        interpretation: CutRunMotifAbsentOccupancyInterpretation,
    ) -> &'static str {
        match interpretation {
            CutRunMotifAbsentOccupancyInterpretation::ContextSupportedByOtherMotifs => {
                "context motifs"
            }
            CutRunMotifAbsentOccupancyInterpretation::MotifPoorSupported => "motif poor",
        }
    }

    pub(super) fn collect_cutrun_regulatory_support_request(
        &self,
    ) -> Result<
        (
            String,
            Vec<String>,
            Vec<String>,
            Option<usize>,
            Option<usize>,
            usize,
            Vec<String>,
        ),
        String,
    > {
        let seq_id = self.seq_id.clone().unwrap_or_default();
        if seq_id.trim().is_empty() {
            return Err("No active genome-anchored sequence".to_string());
        }
        let dataset_ids = Self::parse_ids(&self.cutrun_regulatory_dataset_ids);
        let read_report_ids = Self::parse_ids(&self.cutrun_regulatory_read_report_ids);
        if dataset_ids.is_empty() && read_report_ids.is_empty() {
            return Err(
                "Provide at least one CUT&RUN dataset id or saved read-report id".to_string(),
            );
        }
        let promoter_start = Self::parse_optional_usize_text(
            &self.cutrun_regulatory_promoter_start_0based,
            "CUT&RUN promoter start",
        )?;
        let promoter_end = Self::parse_optional_usize_text(
            &self.cutrun_regulatory_promoter_end_0based_exclusive,
            "CUT&RUN promoter end",
        )?;
        let neighbor_window_bp = Self::parse_positive_usize_text(
            &self.cutrun_regulatory_neighbor_window_bp,
            "CUT&RUN neighbor window bp",
        )?;
        let species_filters = Self::parse_ids(&self.cutrun_regulatory_species_filters);
        Ok((
            seq_id,
            dataset_ids,
            read_report_ids,
            promoter_start,
            promoter_end,
            neighbor_window_bp,
            species_filters,
        ))
    }

    pub(super) fn summarize_cutrun_regulatory_support_report_status(
        report: &CutRunRegulatorySupportReport,
    ) -> String {
        format!(
            "CUT&RUN regulatory support for '{}' merged {} source(s), {} support window(s), {} confirmed TFBS, {} unconfirmed TFBS, {} motif-absent supported window(s)",
            report.seq_id,
            report.evidence_sources.len(),
            report.support_windows.len(),
            report.confirmed_tfbs_rows.len(),
            report.unconfirmed_tfbs_rows.len(),
            report.motif_absent_supported_windows.len()
        )
    }

    pub(super) fn inspect_cutrun_regulatory_support_for_active_sequence(&mut self) {
        let (
            seq_id,
            dataset_ids,
            read_report_ids,
            promoter_search_start_0based,
            promoter_search_end_0based_exclusive,
            neighbor_window_bp,
            species_filters,
        ) = match self.collect_cutrun_regulatory_support_request() {
            Ok(request) => request,
            Err(message) => {
                self.op_status = message;
                return;
            }
        };
        let Some(result) = self.apply_operation_with_feedback_and_result(
            Operation::InspectCutRunRegulatorySupport {
                seq_id,
                dataset_ids,
                read_report_ids,
                promoter_search_start_0based,
                promoter_search_end_0based_exclusive,
                neighbor_window_bp,
                species_filters,
                path: None,
            },
        ) else {
            return;
        };
        if let Some(report) = result.cutrun_regulatory_support {
            self.op_status = Self::summarize_cutrun_regulatory_support_report_status(&report);
            self.cached_cutrun_regulatory_support = Some(report);
            self.op_error_popup = None;
        }
    }

    pub(super) fn seed_latest_cutrun_read_report_for_active_sequence(&mut self) {
        let seq_id = self.seq_id.clone().unwrap_or_default();
        if seq_id.trim().is_empty() {
            self.op_status = "No active sequence for CUT&RUN read-report lookup".to_string();
            return;
        }
        let Some(engine) = self.engine.as_ref() else {
            self.op_status = "No engine attached".to_string();
            return;
        };
        let reports = match engine.read() {
            Ok(guard) => guard.list_cutrun_read_reports(Some(seq_id.as_str())),
            Err(_) => {
                self.op_status =
                    "Engine lock poisoned while listing CUT&RUN read reports".to_string();
                return;
            }
        };
        if let Some(report) = reports.first() {
            self.cutrun_regulatory_read_report_ids = report.report_id.clone();
            self.cached_cutrun_regulatory_support = None;
            self.save_engine_ops_state();
            self.op_status = format!(
                "Using latest CUT&RUN read report '{}' for '{}'",
                report.report_id, seq_id
            );
        } else {
            self.op_status = format!("No saved CUT&RUN read reports found for '{}'", seq_id);
        }
    }

    pub(super) fn render_cutrun_regulatory_support_summary_panel(&mut self, ui: &mut egui::Ui) {
        let navigate_to = {
            let Some(report) = self.cached_cutrun_regulatory_support.as_ref() else {
                ui.small(
                    egui::RichText::new(
                        "No CUT&RUN regulatory-support report cached yet. Provide dataset and/or read-report ids, then run the shared engine inspection.",
                    )
                    .color(egui::Color32::from_rgb(100, 116, 139)),
                );
                return;
            };
            let mut navigate_to: Option<(usize, usize, String)> = None;
            ui.group(|ui| {
                ui.label(egui::RichText::new("CUT&RUN regulatory support").strong());
                ui.small(
                    egui::RichText::new(format!(
                        "{} source(s) | {} support window(s) | {} confirmed TFBS | {} unconfirmed TFBS | {} motif-absent strong window(s)",
                        report.evidence_sources.len(),
                        report.support_windows.len(),
                        report.confirmed_tfbs_rows.len(),
                        report.unconfirmed_tfbs_rows.len(),
                        report.motif_absent_supported_windows.len()
                    ))
                    .color(egui::Color32::from_rgb(71, 85, 105)),
                );
                ui.small(
                    egui::RichText::new(format!(
                        "seq {} | promoter span {}..{} | neighbor window {} bp | motif-context quantile {:.2}",
                        report.seq_id,
                        report.promoter_search_start_0based,
                        report.promoter_search_end_0based_exclusive,
                        report.neighbor_window_bp,
                        report.motif_context_min_llr_quantile
                    ))
                    .color(egui::Color32::from_rgb(100, 116, 139)),
                );
                if !report.species_filters.is_empty() {
                    ui.small(
                        egui::RichText::new(format!(
                            "species filters: {}",
                            report.species_filters.join(", ")
                        ))
                        .color(egui::Color32::from_rgb(100, 116, 139)),
                    );
                }
                if !report.warnings.is_empty() {
                    ui.collapsing("Warnings", |ui| {
                        for warning in &report.warnings {
                            ui.small(
                                egui::RichText::new(warning)
                                    .color(egui::Color32::from_rgb(180, 83, 9)),
                            );
                        }
                    });
                }
                if !report.evidence_sources.is_empty() {
                    ui.collapsing("Evidence sources", |ui| {
                        egui::Grid::new(("cutrun_evidence_sources", report.seq_id.as_str()))
                            .num_columns(4)
                            .striped(true)
                            .show(ui, |ui| {
                                ui.small(egui::RichText::new("kind").strong());
                                ui.small(egui::RichText::new("id").strong());
                                ui.small(egui::RichText::new("factor").strong());
                                ui.small(egui::RichText::new("species").strong());
                                ui.end_row();
                                for source in &report.evidence_sources {
                                    ui.small(Self::cutrun_evidence_source_kind_label(
                                        source.source_kind,
                                    ));
                                    ui.small(&source.source_id);
                                    ui.small(source.target_factor.as_deref().unwrap_or("-"));
                                    ui.small(source.species.as_deref().unwrap_or("-"));
                                    ui.end_row();
                                }
                            });
                    });
                }
                ui.collapsing("Support windows", |ui| {
                    if report.support_windows.is_empty() {
                        ui.small("No support windows were derived from the selected evidence.");
                        return;
                    }
                    let shown = report.support_windows.iter().take(20).collect::<Vec<_>>();
                    egui::Grid::new(("cutrun_support_windows", report.seq_id.as_str()))
                        .num_columns(6)
                        .striped(true)
                        .show(ui, |ui| {
                            ui.small(egui::RichText::new("window").strong());
                            ui.small(egui::RichText::new("local").strong());
                            ui.small(egui::RichText::new("strength").strong());
                            ui.small(egui::RichText::new("peaks").strong());
                            ui.small(egui::RichText::new("frags").strong());
                            ui.small(egui::RichText::new("cuts").strong());
                            ui.end_row();
                            for window in &shown {
                                let response = ui
                                    .link(&window.window_id)
                                    .on_hover_text("Select and center this support window");
                                if response.clicked() {
                                    navigate_to = Some((
                                        window.local_start_0based,
                                        window.local_end_0based_exclusive,
                                        format!("CUT&RUN support {}", window.window_id),
                                    ));
                                }
                                ui.small(format!(
                                    "{}..{}",
                                    window.local_start_0based, window.local_end_0based_exclusive
                                ));
                                ui.small(Self::cutrun_support_strength_label(
                                    window.support_strength,
                                ));
                                ui.small(window.overlapping_peak_count.to_string());
                                ui.small(window.supporting_fragment_count.to_string());
                                ui.small(window.cut_site_count.to_string());
                                ui.end_row();
                            }
                        });
                    if report.support_windows.len() > shown.len() {
                        ui.small(
                            egui::RichText::new(format!(
                                "showing first {} of {} support windows",
                                shown.len(),
                                report.support_windows.len()
                            ))
                            .color(egui::Color32::from_rgb(100, 116, 139)),
                        );
                    }
                });
                ui.collapsing("TFBS confirmation", |ui| {
                    let combined = report
                        .confirmed_tfbs_rows
                        .iter()
                        .chain(report.unconfirmed_tfbs_rows.iter())
                        .take(24)
                        .collect::<Vec<_>>();
                    if combined.is_empty() {
                        ui.small("No TFBS feature rows were present in the inspected span.");
                        return;
                    }
                    egui::Grid::new(("cutrun_tfbs_rows", report.seq_id.as_str()))
                        .num_columns(7)
                        .striped(true)
                        .show(ui, |ui| {
                            ui.small(egui::RichText::new("TF").strong());
                            ui.small(egui::RichText::new("local").strong());
                            ui.small(egui::RichText::new("status").strong());
                            ui.small(egui::RichText::new("support").strong());
                            ui.small(egui::RichText::new("peaks").strong());
                            ui.small(egui::RichText::new("frags").strong());
                            ui.small(egui::RichText::new("cuts").strong());
                            ui.end_row();
                            for row in &combined {
                                let label = row
                                    .motif_label
                                    .as_deref()
                                    .or(row.motif_id.as_deref())
                                    .unwrap_or(&row.feature_label);
                                let response = ui
                                    .link(label)
                                    .on_hover_text("Select and center this TFBS row");
                                if response.clicked() {
                                    navigate_to = Some((
                                        row.local_start_0based,
                                        row.local_end_0based_exclusive,
                                        format!("CUT&RUN TFBS {label}"),
                                    ));
                                }
                                ui.small(format!(
                                    "{}..{}{}",
                                    row.local_start_0based,
                                    row.local_end_0based_exclusive,
                                    row.strand
                                ));
                                ui.small(Self::cutrun_tfbs_confirmation_label(
                                    row.confirmation_status,
                                ));
                                ui.small(
                                    row.strongest_support_strength
                                        .map(Self::cutrun_support_strength_label)
                                        .unwrap_or("-"),
                                );
                                ui.small(row.overlapping_peak_count.to_string());
                                ui.small(row.supporting_fragment_count.to_string());
                                ui.small(row.cut_site_count.to_string());
                                ui.end_row();
                            }
                        });
                    let total = report
                        .confirmed_tfbs_rows
                        .len()
                        .saturating_add(report.unconfirmed_tfbs_rows.len());
                    if total > combined.len() {
                        ui.small(
                            egui::RichText::new(format!(
                                "showing first {} of {} TFBS rows",
                                combined.len(),
                                total
                            ))
                            .color(egui::Color32::from_rgb(100, 116, 139)),
                        );
                    }
                });
                ui.collapsing("Motif-absent supported windows", |ui| {
                    if report.motif_absent_supported_windows.is_empty() {
                        ui.small("No strong support windows lacked the selected target motif.");
                    } else {
                        egui::Grid::new(("cutrun_motif_absent", report.seq_id.as_str()))
                            .num_columns(5)
                            .striped(true)
                            .show(ui, |ui| {
                                ui.small(egui::RichText::new("window").strong());
                                ui.small(egui::RichText::new("local").strong());
                                ui.small(egui::RichText::new("interpretation").strong());
                                ui.small(egui::RichText::new("inside motifs").strong());
                                ui.small(egui::RichText::new("nearby motifs").strong());
                                ui.end_row();
                                for window in report.motif_absent_supported_windows.iter().take(16)
                                {
                                    let response = ui
                                        .link(&window.window_id)
                                        .on_hover_text("Select and center this motif-absent window");
                                    if response.clicked() {
                                        navigate_to = Some((
                                            window.local_start_0based,
                                            window.local_end_0based_exclusive,
                                            format!("CUT&RUN motif-absent {}", window.window_id),
                                        ));
                                    }
                                    ui.small(format!(
                                        "{}..{}",
                                        window.local_start_0based,
                                        window.local_end_0based_exclusive
                                    ));
                                    ui.small(Self::cutrun_occupancy_interpretation_label(
                                        window.occupancy_interpretation,
                                    ));
                                    ui.small(Self::cutrun_context_hit_preview(
                                        &window.motifs_inside_window,
                                    ));
                                    ui.small(Self::cutrun_context_hit_preview(
                                        &window.motifs_in_neighbor_window,
                                    ));
                                    ui.end_row();
                                }
                            });
                    }
                    let summaries = report
                        .common_motifs_inside_supported_windows
                        .iter()
                        .chain(report.common_motifs_near_supported_windows.iter())
                        .take(12)
                        .collect::<Vec<_>>();
                    if !summaries.is_empty() {
                        ui.separator();
                        ui.small(egui::RichText::new("Recurring motif context").strong());
                        for row in summaries {
                            let label = row.motif_label.as_deref().unwrap_or(&row.motif_id);
                            ui.small(format!(
                                "{} {}: {}/{} windows, mean {:.2}, max {:.2}",
                                Self::cutrun_motif_context_scope_label(row.context_scope),
                                label,
                                row.window_count,
                                report.motif_absent_supported_windows.len().max(1),
                                row.mean_best_score,
                                row.max_best_score
                            ));
                        }
                    }
                });
            });
            navigate_to
        };
        if let Some((start, end_exclusive, context_label)) = navigate_to {
            let _ = self.inspect_sequence_span_0based(start, end_exclusive, &context_label);
        }
    }

    pub(super) fn cutrun_context_hit_preview(
        hits: &[crate::engine::CutRunMotifContextHit],
    ) -> String {
        let preview = hits
            .iter()
            .take(3)
            .map(|hit| {
                format!(
                    "{}({})",
                    hit.motif_label.as_deref().unwrap_or(&hit.motif_id),
                    hit.hit_count
                )
            })
            .collect::<Vec<_>>();
        if preview.is_empty() {
            "-".to_string()
        } else {
            preview.join(", ")
        }
    }

    pub(super) fn default_cached_cutrun_regulatory_support_json_file_name(
        report: &CutRunRegulatorySupportReport,
    ) -> String {
        let seq_id = if report.seq_id.trim().is_empty() {
            "cutrun_regulatory_support"
        } else {
            report.seq_id.trim()
        };
        format!(
            "{}_cutrun_regulatory_support_{}-{}.json",
            Self::sanitize_export_name_component(seq_id, "cutrun_regulatory_support"),
            report.promoter_search_start_0based,
            report.promoter_search_end_0based_exclusive
        )
    }

    pub(super) fn export_cached_cutrun_regulatory_support_json(&mut self) {
        let Some(report) = self.cached_cutrun_regulatory_support.clone() else {
            self.op_status =
                "No cached CUT&RUN regulatory-support report available for JSON export".to_string();
            return;
        };
        let default_name = Self::default_cached_cutrun_regulatory_support_json_file_name(&report);
        let Some(path) = rfd::FileDialog::new()
            .set_file_name(&default_name)
            .save_file()
        else {
            self.op_status = "CUT&RUN regulatory-support export canceled".to_string();
            return;
        };
        let dataset_ids = report
            .evidence_sources
            .iter()
            .filter_map(|source| source.dataset_id.clone())
            .collect::<Vec<_>>();
        let read_report_ids = report
            .evidence_sources
            .iter()
            .filter_map(|source| source.report_id.clone())
            .collect::<Vec<_>>();
        let result = self.apply_operation_with_feedback_and_result(
            Operation::InspectCutRunRegulatorySupport {
                seq_id: report.seq_id.clone(),
                dataset_ids,
                read_report_ids,
                promoter_search_start_0based: Some(report.promoter_search_start_0based),
                promoter_search_end_0based_exclusive: Some(
                    report.promoter_search_end_0based_exclusive,
                ),
                neighbor_window_bp: report.neighbor_window_bp,
                species_filters: report.species_filters.clone(),
                path: Some(path.display().to_string()),
            },
        );
        if let Some(updated_report) = result.and_then(|row| row.cutrun_regulatory_support) {
            self.op_status =
                Self::summarize_cutrun_regulatory_support_report_status(&updated_report);
            self.cached_cutrun_regulatory_support = Some(updated_report);
        }
    }
}
