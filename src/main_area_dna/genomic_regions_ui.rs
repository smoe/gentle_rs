//! Portable genomic-region manager for the DNA sequence window.
//!
//! The GUI only stages user intent and presents engine-owned records. Genome
//! projection, coordinate conversion, provenance, digesting, and persistence
//! remain in `GentleEngine` so the same behavior is available headlessly.

use super::*;

#[derive(Debug, Clone)]
enum GenomicRegionManagerAction {
    Refresh,
    SavePendingSelection,
    CopyHuman(gentle_protocol::GenomicRegionOfInterest),
    CopyBed(gentle_protocol::GenomicRegionOfInterest),
    CopyJson(gentle_protocol::GenomicRegionOfInterest),
    ImportJson,
    ImportBed,
    ExportJson(String),
    ExportBed(String),
}

impl MainAreaDna {
    pub(super) fn open_genomic_region_manager(&mut self, selection: Option<(usize, usize)>) {
        self.show_genomic_region_manager = true;
        if let Some((start, end)) = selection {
            self.genomic_region_pending_selection = Some((start.min(end), start.max(end)));
            if self.genomic_region_new_label.trim().is_empty() {
                self.genomic_region_new_label = format!(
                    "Selected region {}-{}",
                    start.min(end).saturating_add(1),
                    start.max(end)
                );
            }
        }
        self.refresh_genomic_region_store_cache();
    }

    fn refresh_genomic_region_store_cache(&mut self) {
        let Some(engine) = self.engine.as_ref() else {
            self.genomic_region_status = "No engine attached".to_string();
            self.genomic_region_store_cache = None;
            return;
        };
        match engine
            .try_read()
            .map_err(|_| "Engine is busy; try refreshing again".to_string())
            .and_then(|guard| {
                guard
                    .genomic_region_store_snapshot()
                    .map_err(|error| error.message)
            }) {
            Ok(store) => self.genomic_region_store_cache = Some(store),
            Err(error) => {
                self.genomic_region_status = error;
                self.genomic_region_store_cache = None;
            }
        }
    }

    fn apply_genomic_region_operation(&mut self, operation: Operation) -> bool {
        let Some(result) = self.apply_operation_with_feedback_and_result(operation) else {
            return false;
        };
        let Some(report) = result.genomic_region_operation else {
            self.genomic_region_status =
                "GENtle returned no genomic-region operation report".to_string();
            return false;
        };
        self.genomic_region_status = if let Some(region) = report.region.as_ref() {
            format!(
                "{}: {} ({})",
                report.action, region.region_id, region.content_sha256
            )
        } else if !report.written_artifacts.is_empty() {
            format!(
                "{}: wrote {} artifact(s)",
                report.action,
                report.written_artifacts.len()
            )
        } else {
            report.action.clone()
        };
        self.refresh_genomic_region_store_cache();
        true
    }

    fn save_pending_genomic_region_selection(&mut self) {
        let Some((start, end)) = self.genomic_region_pending_selection else {
            self.genomic_region_status = "No non-empty sequence selection is staged".to_string();
            return;
        };
        if start >= end {
            self.genomic_region_status = "The staged sequence selection is empty".to_string();
            return;
        }
        let Some(seq_id) = self.seq_id.clone().filter(|value| !value.trim().is_empty()) else {
            self.genomic_region_status = "No active sequence is available".to_string();
            return;
        };
        let label = self.genomic_region_new_label.trim();
        let request = gentle_protocol::GenomicRegionCaptureRequest {
            set_id: self.genomic_region_default_set_id.trim().to_string(),
            set_label: (!self.genomic_region_default_set_label.trim().is_empty())
                .then(|| self.genomic_region_default_set_label.trim().to_string()),
            label: (!label.is_empty()).then(|| label.to_string()),
            purpose: self.genomic_region_new_purpose,
            source: gentle_protocol::GenomicRegionCaptureSource::SequenceSelection {
                seq_id,
                local_start_0based: start as u64,
                local_end_0based_exclusive: end as u64,
                strand: gentle_protocol::GenomicRegionStrand::Unstranded,
                reference_override: None,
            },
            collision_policy: gentle_protocol::GenomicRegionCollisionPolicy::Reject,
            ..Default::default()
        };
        if self.apply_genomic_region_operation(Operation::CaptureGenomicRegion { request }) {
            self.genomic_region_pending_selection = None;
            self.genomic_region_new_label.clear();
        }
    }

    pub(super) fn capture_cutrun_support_window_as_region(
        &mut self,
        report: CutRunRegulatorySupportReport,
        window_id: String,
    ) {
        let Some(anchor) = self.active_sequence_anchor_summary() else {
            self.genomic_region_status =
                "An exact genome anchor is required to save a CUT&RUN support window".to_string();
            self.show_genomic_region_manager = true;
            return;
        };
        let reference = gentle_protocol::GenomicRegionReference {
            assembly_name: anchor.genome_id,
            contig_name: anchor.chromosome,
            ..Default::default()
        };
        let request = gentle_protocol::GenomicRegionCaptureRequest {
            set_id: self.genomic_region_default_set_id.trim().to_string(),
            set_label: Some(self.genomic_region_default_set_label.trim().to_string()),
            label: Some(format!("CUT&RUN {window_id}")),
            purpose: gentle_protocol::GenomicRegionPurpose::OccupancyRegion,
            source: gentle_protocol::GenomicRegionCaptureSource::CutrunSupportWindow {
                report: Box::new(report),
                window_id,
                reference,
            },
            collision_policy: gentle_protocol::GenomicRegionCollisionPolicy::Reject,
            ..Default::default()
        };
        self.show_genomic_region_manager = true;
        let _ = self.apply_genomic_region_operation(Operation::CaptureGenomicRegion { request });
    }

    pub(super) fn capture_gene_locus_ensembl_region(
        &mut self,
        report: &GeneLocusEvidenceDisplayReport,
        row: gentle_protocol::GeneLocusEnsemblRegulationFeatureRow,
    ) {
        let Some(ensembl) = report.ensembl_regulation.as_ref() else {
            self.genomic_region_status =
                "This locus report has no verified Ensembl Regulation source binding".to_string();
            self.show_genomic_region_manager = true;
            return;
        };
        let Some(source_binding) = ensembl.source_binding.clone() else {
            self.genomic_region_status =
                "This locus report has no verified Ensembl Regulation source binding".to_string();
            self.show_genomic_region_manager = true;
            return;
        };
        let Some(anchor) = self.active_sequence_anchor_summary() else {
            self.genomic_region_status =
                "An exact genome anchor is required to save this Ensembl feature".to_string();
            self.show_genomic_region_manager = true;
            return;
        };
        let reference = gentle_protocol::GenomicRegionReference {
            species_scientific_name: Some(source_binding.source.species_scientific_name.clone()),
            taxon_id: Some(source_binding.source.taxon_id),
            assembly_name: source_binding.source.assembly_name.clone(),
            assembly_accession: Some(source_binding.source.assembly_accession.clone())
                .filter(|value| !value.trim().is_empty()),
            contig_name: anchor.chromosome,
            ..Default::default()
        };
        let request = gentle_protocol::GenomicRegionCaptureRequest {
            set_id: self.genomic_region_default_set_id.trim().to_string(),
            set_label: Some(self.genomic_region_default_set_label.trim().to_string()),
            label: Some(format!("{} {}", row.feature_type, row.feature_id)),
            purpose: gentle_protocol::GenomicRegionPurpose::CandidateCisRegulatoryRegion,
            source:
                gentle_protocol::GenomicRegionCaptureSource::GeneLocusEnsemblRegulatoryFeature {
                    row,
                    source_binding,
                    reference,
                    seq_id: Some(report.seq_id.clone()),
                    interval_kind: gentle_protocol::GenomicRegionEnsemblIntervalKind::Core,
                    evidence_statement: ensembl.evidence_statement.clone(),
                    non_claims: ensembl.non_claims.clone(),
                },
            collision_policy: gentle_protocol::GenomicRegionCollisionPolicy::Reject,
            ..Default::default()
        };
        self.show_genomic_region_manager = true;
        let _ = self.apply_genomic_region_operation(Operation::CaptureGenomicRegion { request });
    }

    fn set_region_overlay_enabled(&mut self, set_id: &str, enabled: bool) {
        let mut ids = Self::parse_ids(&self.splicing_locus_region_set_ids)
            .into_iter()
            .collect::<BTreeSet<_>>();
        if enabled {
            ids.insert(set_id.to_string());
        } else {
            ids.remove(set_id);
        }
        self.splicing_locus_region_set_ids = ids.into_iter().collect::<Vec<_>>().join(", ");
    }

    fn export_genomic_region_set_json(&mut self, set_id: &str) {
        let Some(path) = rfd::FileDialog::new()
            .set_file_name(format!("{set_id}.genomic-regions.json"))
            .save_file()
        else {
            self.genomic_region_status = "JSON export canceled".to_string();
            return;
        };
        let request = gentle_protocol::GenomicRegionExportRequest {
            set_id: set_id.to_string(),
            json_path: Some(path.display().to_string()),
            ..Default::default()
        };
        let _ = self.apply_genomic_region_operation(Operation::ExportGenomicRegionSet { request });
    }

    fn import_genomic_region_set_json(&mut self) {
        let Some(path) = rfd::FileDialog::new()
            .add_filter("GENtle genomic region set", &["json"])
            .pick_file()
        else {
            self.genomic_region_status = "JSON import canceled".to_string();
            return;
        };
        let request = gentle_protocol::GenomicRegionImportRequest {
            path: path.display().to_string(),
            format: gentle_protocol::GenomicRegionImportFormat::Json,
            collision_policy: gentle_protocol::GenomicRegionCollisionPolicy::Reject,
            ..Default::default()
        };
        let _ = self.apply_genomic_region_operation(Operation::ImportGenomicRegionSet { request });
    }

    fn import_genomic_region_set_bed(&mut self) {
        let Some(bed_path) = rfd::FileDialog::new()
            .add_filter("BED6 genomic regions", &["bed"])
            .pick_file()
        else {
            self.genomic_region_status = "BED import canceled".to_string();
            return;
        };
        let mut manifest_dialog = rfd::FileDialog::new()
            .add_filter("GENtle BED manifest", &["json"])
            .set_file_name(format!(
                "{}.manifest.json",
                bed_path
                    .file_name()
                    .and_then(|value| value.to_str())
                    .unwrap_or("regions.bed")
            ));
        if let Some(parent) = bed_path.parent() {
            manifest_dialog = manifest_dialog.set_directory(parent);
        }
        let Some(manifest_path) = manifest_dialog.pick_file() else {
            self.genomic_region_status = "BED manifest selection canceled".to_string();
            return;
        };
        let request = gentle_protocol::GenomicRegionImportRequest {
            path: bed_path.display().to_string(),
            format: gentle_protocol::GenomicRegionImportFormat::Bed,
            manifest_path: Some(manifest_path.display().to_string()),
            collision_policy: gentle_protocol::GenomicRegionCollisionPolicy::Reject,
            ..Default::default()
        };
        let _ = self.apply_genomic_region_operation(Operation::ImportGenomicRegionSet { request });
    }

    fn export_genomic_region_set_bed(&mut self, set_id: &str) {
        let Some(path) = rfd::FileDialog::new()
            .set_file_name(format!("{set_id}.bed"))
            .save_file()
        else {
            self.genomic_region_status = "BED export canceled".to_string();
            return;
        };
        let manifest_path = format!("{}.manifest.json", path.display());
        let request = gentle_protocol::GenomicRegionExportRequest {
            set_id: set_id.to_string(),
            bed_path: Some(path.display().to_string()),
            manifest_path: Some(manifest_path),
            ..Default::default()
        };
        let _ = self.apply_genomic_region_operation(Operation::ExportGenomicRegionSet { request });
    }

    pub(super) fn render_genomic_region_manager(&mut self, ctx: &egui::Context) {
        if !self.show_genomic_region_manager {
            return;
        }
        let store = self.genomic_region_store_cache.clone().unwrap_or_default();
        let selected_overlay_ids = Self::parse_ids(&self.splicing_locus_region_set_ids)
            .into_iter()
            .collect::<BTreeSet<_>>();
        #[cfg(feature = "gui-test-support")]
        let subject_scope = crate::gui_test_support::pseudonymous_subject_scope(&[self
            .seq_id
            .as_deref()
            .unwrap_or("unnamed")]);
        let mut action: Option<GenomicRegionManagerAction> = None;
        let mut overlay_change: Option<(String, bool)> = None;
        let mut open = self.show_genomic_region_manager;
        let spec = crate::egui_compat::HostedWindowSpec::new(
            "Saved genomic regions",
            egui::Id::new(("genomic_region_manager", self.panel_scope_key())),
            Vec2::new(980.0, 700.0),
            Vec2::new(720.0, 480.0),
        );
        crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
            #[cfg(feature = "gui-test-support")]
            crate::gui_test_support::register_rect(
                ui.ctx().clone(),
                crate::tutorial_gui_semantics::WINDOW_GENOMIC_REGIONS,
                crate::tutorial_gui_semantics::WINDOW_GENOMIC_REGIONS,
                Some(&subject_scope),
                crate::gui_test_support::GuiTestWidgetKind::Window,
                ui.max_rect(),
                true,
                true,
                true,
                Some("ready"),
            );
            egui::ScrollArea::both()
                .id_salt(("genomic_region_manager_scroll", self.panel_scope_key()))
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    scroll_input_policy::apply_scrollarea_keyboard_navigation(
                        ui,
                        scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                    );
                    ui.set_min_width(930.0);
                    ui.heading("Portable genomic regions of interest");
                    ui.small(
                        "Canonical coordinates are assembly-bound, 0-based and half-open. Human copy text is explicitly 1-based and inclusive; BED remains 0-based and half-open.",
                    );
                    ui.small(
                        "A purpose or evidence source records why a span was retained. It does not establish causal regulation, biochemical affinity, or a target-gene assignment.",
                    );
                    ui.separator();
                    ui.horizontal_wrapped(|ui| {
                        ui.label("Set id");
                        ui.text_edit_singleline(&mut self.genomic_region_default_set_id);
                        ui.label("Set label");
                        ui.text_edit_singleline(&mut self.genomic_region_default_set_label);
                        let refresh = ui.button("Refresh");
                        #[cfg(feature = "gui-test-support")]
                        crate::gui_test_support::register_response(
                            &refresh,
                            crate::tutorial_gui_semantics::GENOMIC_REGION_REFRESH,
                            crate::tutorial_gui_semantics::WINDOW_GENOMIC_REGIONS,
                            Some(&subject_scope),
                            crate::gui_test_support::GuiTestWidgetKind::Button,
                            false,
                        );
                        if refresh.clicked() {
                            action = Some(GenomicRegionManagerAction::Refresh);
                        }
                        let import_json = ui.button("Import JSON...");
                        #[cfg(feature = "gui-test-support")]
                        crate::gui_test_support::register_response(
                            &import_json,
                            crate::tutorial_gui_semantics::GENOMIC_REGION_IMPORT_JSON,
                            crate::tutorial_gui_semantics::WINDOW_GENOMIC_REGIONS,
                            Some(&subject_scope),
                            crate::gui_test_support::GuiTestWidgetKind::Button,
                            false,
                        );
                        if import_json.clicked() {
                            action = Some(GenomicRegionManagerAction::ImportJson);
                        }
                        let import_bed = ui.button("Import BED + manifest...");
                        #[cfg(feature = "gui-test-support")]
                        crate::gui_test_support::register_response(
                            &import_bed,
                            crate::tutorial_gui_semantics::GENOMIC_REGION_IMPORT_BED,
                            crate::tutorial_gui_semantics::WINDOW_GENOMIC_REGIONS,
                            Some(&subject_scope),
                            crate::gui_test_support::GuiTestWidgetKind::Button,
                            false,
                        );
                        if import_bed.clicked() {
                            action = Some(GenomicRegionManagerAction::ImportBed);
                        }
                    });
                    if let Some((start, end)) = self.genomic_region_pending_selection {
                        ui.group(|ui| {
                            ui.label(egui::RichText::new("Staged sequence selection").strong());
                            ui.label(format!(
                                "local {}..{} (0-based, end-exclusive; {} bp)",
                                start,
                                end,
                                end.saturating_sub(start)
                            ));
                            ui.horizontal_wrapped(|ui| {
                                ui.label("Label");
                                ui.text_edit_singleline(&mut self.genomic_region_new_label);
                                egui::ComboBox::from_id_salt("genomic_region_new_purpose")
                                    .selected_text(self.genomic_region_new_purpose.as_str())
                                    .show_ui(ui, |ui| {
                                        for purpose in [
                                            gentle_protocol::GenomicRegionPurpose::CandidateCisRegulatoryRegion,
                                            gentle_protocol::GenomicRegionPurpose::OccupancyRegion,
                                            gentle_protocol::GenomicRegionPurpose::PromoterRegion,
                                            gentle_protocol::GenomicRegionPurpose::ReporterCandidate,
                                            gentle_protocol::GenomicRegionPurpose::Other,
                                        ] {
                                            ui.selectable_value(
                                                &mut self.genomic_region_new_purpose,
                                                purpose,
                                                purpose.as_str(),
                                            );
                                        }
                                    });
                                let save = ui.button("Save/share region");
                                #[cfg(feature = "gui-test-support")]
                                crate::gui_test_support::register_response(
                                    &save,
                                    crate::tutorial_gui_semantics::GENOMIC_REGION_SAVE_PENDING,
                                    crate::tutorial_gui_semantics::WINDOW_GENOMIC_REGIONS,
                                    Some(&subject_scope),
                                    crate::gui_test_support::GuiTestWidgetKind::Button,
                                    false,
                                );
                                if save.clicked() {
                                    action = Some(GenomicRegionManagerAction::SavePendingSelection);
                                }
                            });
                        });
                    }
                    ui.separator();
                    if store.sets.is_empty() {
                        ui.label("No saved genomic regions in this project.");
                    }
                    for set in &store.sets {
                        ui.horizontal_wrapped(|ui| {
                            ui.heading(set.label.as_deref().unwrap_or(&set.set_id));
                            ui.monospace(format!("{} region(s)", set.regions.len()));
                            let mut enabled = selected_overlay_ids.contains(&set.set_id);
                            if ui
                                .checkbox(&mut enabled, "Show in locus figure")
                                .on_hover_text(
                                    "Include this saved set when the gene-locus report is next composed",
                                )
                                .changed()
                            {
                                overlay_change = Some((set.set_id.clone(), enabled));
                            }
                            let export_json = ui.button("Export JSON...");
                            #[cfg(feature = "gui-test-support")]
                            crate::gui_test_support::register_response(
                                &export_json,
                                crate::tutorial_gui_semantics::GENOMIC_REGION_EXPORT_JSON,
                                crate::tutorial_gui_semantics::WINDOW_GENOMIC_REGIONS,
                                Some(&crate::gui_test_support::pseudonymous_subject_scope(&[
                                    self.seq_id.as_deref().unwrap_or("unnamed"),
                                    &set.set_id,
                                ])),
                                crate::gui_test_support::GuiTestWidgetKind::Button,
                                false,
                            );
                            if export_json.clicked() {
                                action = Some(GenomicRegionManagerAction::ExportJson(
                                    set.set_id.clone(),
                                ));
                            }
                            let export_bed = ui.button("Export BED + manifest...");
                            #[cfg(feature = "gui-test-support")]
                            crate::gui_test_support::register_response(
                                &export_bed,
                                crate::tutorial_gui_semantics::GENOMIC_REGION_EXPORT_BED,
                                crate::tutorial_gui_semantics::WINDOW_GENOMIC_REGIONS,
                                Some(&crate::gui_test_support::pseudonymous_subject_scope(&[
                                    self.seq_id.as_deref().unwrap_or("unnamed"),
                                    &set.set_id,
                                ])),
                                crate::gui_test_support::GuiTestWidgetKind::Button,
                                false,
                            );
                            if export_bed.clicked() {
                                action = Some(GenomicRegionManagerAction::ExportBed(
                                    set.set_id.clone(),
                                ));
                            }
                        });
                        ui.small(format!("set digest {}", set.content_sha256));
                        egui::Grid::new(("genomic_region_table", set.set_id.as_str()))
                            .num_columns(8)
                            .striped(true)
                            .spacing(Vec2::new(12.0, 4.0))
                            .show(ui, |ui| {
                                for heading in [
                                    "ID", "label", "assembly", "coordinate", "strand", "method",
                                    "evidence", "copy",
                                ] {
                                    ui.small(egui::RichText::new(heading).strong());
                                }
                                ui.end_row();
                                for region in &set.regions {
                                    ui.monospace(&region.region_id);
                                    ui.label(region.label.as_deref().unwrap_or("-"));
                                    ui.label(&region.interval.reference.assembly_name);
                                    ui.monospace(format!(
                                        "{}:{}-{} (1-based)",
                                        region.interval.reference.contig_name,
                                        region.interval.start_0based.saturating_add(1),
                                        region.interval.end_0based_exclusive
                                    ));
                                    ui.label(region.interval.strand.human_value());
                                    ui.label(region.selection_method.as_str());
                                    let availability = region
                                        .evidence
                                        .iter()
                                        .map(|item| item.availability)
                                        .find(|status| {
                                            *status
                                                != gentle_protocol::GenomicRegionEvidenceAvailability::Available
                                        })
                                        .unwrap_or(
                                            gentle_protocol::GenomicRegionEvidenceAvailability::Available,
                                        );
                                    ui.label(format!(
                                        "{} ({})",
                                        availability.as_str(),
                                        region.evidence.len()
                                    ));
                                    ui.horizontal(|ui| {
                                        let human = ui.small_button("1-based").on_hover_text(
                                            "Copy human coordinates (1-based inclusive)",
                                        );
                                        #[cfg(feature = "gui-test-support")]
                                        crate::gui_test_support::register_response(
                                            &human,
                                            crate::tutorial_gui_semantics::GENOMIC_REGION_COPY_HUMAN,
                                            crate::tutorial_gui_semantics::WINDOW_GENOMIC_REGIONS,
                                            Some(&crate::gui_test_support::pseudonymous_subject_scope(&[
                                                self.seq_id.as_deref().unwrap_or("unnamed"),
                                                &set.set_id,
                                                &region.region_id,
                                            ])),
                                            crate::gui_test_support::GuiTestWidgetKind::Button,
                                            false,
                                        );
                                        if human.clicked() {
                                            action = Some(GenomicRegionManagerAction::CopyHuman(
                                                region.clone(),
                                            ));
                                        }
                                        let bed = ui.small_button("BED").on_hover_text(
                                            "Copy BED row (0-based half-open)",
                                        );
                                        #[cfg(feature = "gui-test-support")]
                                        crate::gui_test_support::register_response(
                                            &bed,
                                            crate::tutorial_gui_semantics::GENOMIC_REGION_COPY_BED,
                                            crate::tutorial_gui_semantics::WINDOW_GENOMIC_REGIONS,
                                            Some(&crate::gui_test_support::pseudonymous_subject_scope(&[
                                                self.seq_id.as_deref().unwrap_or("unnamed"),
                                                &set.set_id,
                                                &region.region_id,
                                            ])),
                                            crate::gui_test_support::GuiTestWidgetKind::Button,
                                            false,
                                        );
                                        if bed.clicked() {
                                            action = Some(GenomicRegionManagerAction::CopyBed(
                                                region.clone(),
                                            ));
                                        }
                                        let json = ui.small_button("JSON").on_hover_text(
                                            "Copy canonical ROI JSON (lossless)",
                                        );
                                        #[cfg(feature = "gui-test-support")]
                                        crate::gui_test_support::register_response(
                                            &json,
                                            crate::tutorial_gui_semantics::GENOMIC_REGION_COPY_JSON,
                                            crate::tutorial_gui_semantics::WINDOW_GENOMIC_REGIONS,
                                            Some(&crate::gui_test_support::pseudonymous_subject_scope(&[
                                                self.seq_id.as_deref().unwrap_or("unnamed"),
                                                &set.set_id,
                                                &region.region_id,
                                            ])),
                                            crate::gui_test_support::GuiTestWidgetKind::Button,
                                            false,
                                        );
                                        if json.clicked() {
                                            action = Some(GenomicRegionManagerAction::CopyJson(
                                                region.clone(),
                                            ));
                                        }
                                    });
                                    ui.end_row();
                                }
                            });
                        ui.separator();
                    }
                    if !self.genomic_region_status.trim().is_empty() {
                        ui.small(&self.genomic_region_status);
                    }
                });
        });
        self.show_genomic_region_manager = open;
        if let Some((set_id, enabled)) = overlay_change {
            self.set_region_overlay_enabled(&set_id, enabled);
        }
        match action {
            Some(GenomicRegionManagerAction::Refresh) => self.refresh_genomic_region_store_cache(),
            Some(GenomicRegionManagerAction::SavePendingSelection) => {
                self.save_pending_genomic_region_selection()
            }
            Some(GenomicRegionManagerAction::CopyHuman(region)) => {
                ctx.copy_text(GentleEngine::genomic_region_human_copy(&region));
                self.genomic_region_status =
                    "Copied human coordinates (1-based inclusive)".to_string();
            }
            Some(GenomicRegionManagerAction::CopyBed(region)) => {
                ctx.copy_text(GentleEngine::genomic_region_bed_row(&region));
                self.genomic_region_status = "Copied BED row (0-based half-open)".to_string();
            }
            Some(GenomicRegionManagerAction::CopyJson(region)) => {
                match GentleEngine::genomic_region_canonical_json(&region) {
                    Ok(json) => {
                        ctx.copy_text(json);
                        self.genomic_region_status = "Copied canonical ROI JSON".to_string();
                    }
                    Err(error) => self.genomic_region_status = error.message,
                }
            }
            Some(GenomicRegionManagerAction::ImportJson) => self.import_genomic_region_set_json(),
            Some(GenomicRegionManagerAction::ImportBed) => self.import_genomic_region_set_bed(),
            Some(GenomicRegionManagerAction::ExportJson(set_id)) => {
                self.export_genomic_region_set_json(&set_id)
            }
            Some(GenomicRegionManagerAction::ExportBed(set_id)) => {
                self.export_genomic_region_set_bed(&set_id)
            }
            None => {}
        }
    }
}
