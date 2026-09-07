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
    OpenConservation {
        set_id: String,
        region: gentle_protocol::GenomicRegionOfInterest,
    },
}

/// A colour edit that is still being dragged.
///
/// Committing on every reported change would apply one engine operation, and
/// therefore one full project checkpoint, per rendered frame.
#[derive(Debug, Clone)]
pub(super) struct StagedGenomicRegionColor {
    pub(super) set_id: String,
    pub(super) region_id: String,
    pub(super) color_hex: String,
}

impl MainAreaDna {
    pub(super) fn open_genomic_region_conservation_by_id(&mut self, set_id: &str, region_id: &str) {
        if self.genomic_region_store_cache.is_none() {
            self.refresh_genomic_region_store_cache();
        }
        let region = self
            .genomic_region_store_cache
            .as_ref()
            .and_then(|store| store.sets.iter().find(|set| set.set_id == set_id))
            .and_then(|set| {
                set.regions
                    .iter()
                    .find(|region| region.region_id == region_id)
            })
            .cloned();
        let Some(region) = region else {
            self.genomic_region_status = format!(
                "Could not find saved region '{set_id}/{region_id}' for conservation analysis"
            );
            self.show_genomic_region_manager = true;
            return;
        };
        self.open_genomic_region_conservation(set_id, &region);
    }

    fn open_genomic_region_conservation(
        &mut self,
        set_id: &str,
        region: &gentle_protocol::GenomicRegionOfInterest,
    ) {
        let changed = self.genomic_region_conservation_set_id != set_id
            || self.genomic_region_conservation_region_id != region.region_id
            || self.genomic_region_conservation_expected_sha256.as_deref()
                != Some(region.content_sha256.as_str());
        self.genomic_region_conservation_set_id = set_id.to_string();
        self.genomic_region_conservation_region_id = region.region_id.clone();
        self.genomic_region_conservation_expected_sha256 = Some(region.content_sha256.clone());
        self.genomic_region_conservation_query_genome_id = region
            .local_projection
            .as_ref()
            .map(|projection| projection.source_genome_id.clone())
            .filter(|value| !value.trim().is_empty())
            .unwrap_or_else(|| region.interval.reference.assembly_name.clone());
        if changed {
            self.genomic_region_conservation_report = None;
            self.genomic_region_conservation_selected_block_id = None;
            self.genomic_region_conservation_evidence_region_ids.clear();
            self.genomic_region_conservation_module_assessment = None;
            self.genomic_region_conservation_progress = None;
        }
        self.genomic_region_conservation_status =
            "Ready to inspect validated local genomic indexes".to_string();
        self.show_genomic_region_conservation = true;
    }

    pub(super) fn open_genomic_region_manager(&mut self, selection: Option<(usize, usize)>) {
        self.show_genomic_region_manager = true;
        if let Some((start, end)) = selection {
            self.genomic_region_pending_locus_report = None;
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

    pub(super) fn save_pending_genomic_region_selection(&mut self) {
        if let Some(report) = &self.genomic_region_pending_locus_report
            && let Err(error) = self.verify_splicing_locus_binding(report)
        {
            self.genomic_region_status = error;
            return;
        }
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
            display_color_hex: Some(self.genomic_region_new_color_hex.clone()),
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
            self.genomic_region_pending_locus_report = None;
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
            display_color_hex: Some("#DC2626".to_string()),
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
            display_color_hex: Some("#7C3AED".to_string()),
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

    /// Apply a staged colour once the pointer is released.
    ///
    /// Presentation updates are ordinary engine operations, so one per dragged
    /// frame would clone the whole project per frame and evict real undo
    /// history. Keyboard-driven edits commit immediately, since no pointer is
    /// held during them.
    fn commit_staged_genomic_region_color(&mut self, ctx: &egui::Context) {
        if ctx.input(|input| input.pointer.any_down()) {
            return;
        }
        let Some(staged) = self.genomic_region_pending_color.take() else {
            return;
        };
        let Some(region) = self
            .genomic_region_store_cache
            .as_ref()
            .and_then(|store| store.sets.iter().find(|set| set.set_id == staged.set_id))
            .and_then(|set| {
                set.regions
                    .iter()
                    .find(|region| region.region_id == staged.region_id)
            })
            .cloned()
        else {
            self.genomic_region_status =
                "The recoloured genomic region is no longer in the loaded store".to_string();
            return;
        };
        if region.display_color_hex.as_deref() == Some(staged.color_hex.as_str()) {
            return;
        }
        let request = gentle_protocol::GenomicRegionUpdateRequest {
            set_id: staged.set_id,
            region_id: region.region_id,
            label: region.label,
            description: region.description,
            display_color_hex: Some(staged.color_hex),
            notes: region.notes,
        };
        let _ = self
            .apply_genomic_region_operation(Operation::UpdateGenomicRegionPresentation { request });
    }

    pub(super) fn render_genomic_region_manager(&mut self, ctx: &egui::Context) {
        if !self.show_genomic_region_manager {
            self.commit_staged_genomic_region_color(ctx);
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
        let mut staged_color: Option<StagedGenomicRegionColor> = None;
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
                                ui.label("Colour");
                                ui.add(
                                    egui::TextEdit::singleline(
                                        &mut self.genomic_region_new_color_hex,
                                    )
                                    .desired_width(76.0),
                                );
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
                            .num_columns(10)
                            .striped(true)
                            .spacing(Vec2::new(12.0, 4.0))
                            .show(ui, |ui| {
                                for heading in [
                                    "ID", "label", "assembly", "coordinate", "strand", "method",
                                    "evidence", "colour", "copy", "analyse",
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
                                    // The picker reports a change every frame
                                    // it is dragged, and each engine operation
                                    // captures a full project checkpoint, so
                                    // stage the value here and commit it once
                                    // the pointer is released.
                                    let staged_hex = self
                                        .genomic_region_pending_color
                                        .as_ref()
                                        .filter(|staged| {
                                            staged.set_id == set.set_id
                                                && staged.region_id == region.region_id
                                        })
                                        .map(|staged| staged.color_hex.clone());
                                    let mut color = MainAreaDna::locus_inspector_color(
                                        staged_hex
                                            .as_deref()
                                            .or(region.display_color_hex.as_deref()),
                                        egui::Color32::from_rgb(194, 65, 12),
                                    );
                                    if egui::color_picker::color_edit_button_srgba(
                                        ui,
                                        &mut color,
                                        egui::color_picker::Alpha::Opaque,
                                    )
                                    .changed()
                                    {
                                        staged_color = Some(StagedGenomicRegionColor {
                                            set_id: set.set_id.clone(),
                                            region_id: region.region_id.clone(),
                                            color_hex: format!(
                                                "#{:02X}{:02X}{:02X}",
                                                color.r(),
                                                color.g(),
                                                color.b()
                                            ),
                                        });
                                    }
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
                                    if ui
                                        .small_button("Conservation...")
                                        .on_hover_text(
                                            "Compare this saved region with validated local genomic BLAST indexes",
                                        )
                                        .clicked()
                                    {
                                        action = Some(
                                            GenomicRegionManagerAction::OpenConservation {
                                                set_id: set.set_id.clone(),
                                                region: region.clone(),
                                            },
                                        );
                                    }
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
        if let Some(staged) = staged_color {
            self.genomic_region_pending_color = Some(staged);
        }
        self.commit_staged_genomic_region_color(ctx);
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
            Some(GenomicRegionManagerAction::OpenConservation { set_id, region }) => {
                self.open_genomic_region_conservation(&set_id, &region)
            }
            None => {}
        }
    }

    fn start_genomic_region_homology_screen(&mut self) {
        if self.genomic_region_conservation_task.is_some() {
            self.genomic_region_conservation_status =
                "A conservation screen is already running".to_string();
            return;
        }
        let Some(engine) = self.engine.clone() else {
            self.genomic_region_conservation_status = "No engine is attached".to_string();
            return;
        };
        let request = gentle_protocol::GenomicRegionHomologyScreenRequest {
            set_id: self.genomic_region_conservation_set_id.clone(),
            region_id: self.genomic_region_conservation_region_id.clone(),
            expected_region_content_sha256: self
                .genomic_region_conservation_expected_sha256
                .clone(),
            query_genome_id: (!self
                .genomic_region_conservation_query_genome_id
                .trim()
                .is_empty())
            .then(|| {
                self.genomic_region_conservation_query_genome_id
                    .trim()
                    .to_string()
            }),
            ..Default::default()
        };
        let (sender, receiver) = mpsc::channel::<GenomicRegionHomologyTaskMessage>();
        self.genomic_region_conservation_task = Some(GenomicRegionHomologyTask {
            started: Instant::now(),
            receiver: Arc::new(Mutex::new(receiver)),
        });
        self.genomic_region_conservation_progress = None;
        self.genomic_region_conservation_status =
            "Preparing an immutable project snapshot...".to_string();
        std::thread::spawn(move || {
            let progress_sender = sender.clone();
            let result = crate::background_engine::execute_read_only_operation_on_engine_snapshot_with_progress(
                &engine,
                Operation::ScreenGenomicRegionHomology {
                    request,
                    path: None,
                },
                move |progress| match progress {
                    OperationProgress::GenomicRegionHomology(progress) => progress_sender
                        .send(GenomicRegionHomologyTaskMessage::Progress(progress))
                        .is_ok(),
                    _ => true,
                },
            )
            .and_then(|mut result| {
                result
                    .genomic_region_homology
                    .take()
                    .map(|report| *report)
                    .ok_or_else(|| {
                        EngineError::new(
                            ErrorCode::Internal,
                            "Conservation operation returned no typed homology report",
                        )
                    })
            });
            let _ = sender.send(GenomicRegionHomologyTaskMessage::Done(result));
        });
    }

    pub(super) fn poll_genomic_region_homology_task(&mut self, ctx: &egui::Context) {
        let Some(task) = self.genomic_region_conservation_task.as_ref() else {
            return;
        };
        let started = task.started;
        let receiver = Arc::clone(&task.receiver);
        let mut done = None;
        let mut disconnected = false;
        match receiver.lock() {
            Ok(receiver) => loop {
                match receiver.try_recv() {
                    Ok(GenomicRegionHomologyTaskMessage::Progress(progress)) => {
                        self.genomic_region_conservation_status = if progress.target_count > 0 {
                            format!(
                                "{} {}/{}: {} ({:.1}s)",
                                progress.phase,
                                progress.target_ordinal,
                                progress.target_count,
                                progress.detail,
                                started.elapsed().as_secs_f32()
                            )
                        } else {
                            format!(
                                "{}: {} ({:.1}s)",
                                progress.phase,
                                progress.detail,
                                started.elapsed().as_secs_f32()
                            )
                        };
                        self.genomic_region_conservation_progress = Some(progress);
                    }
                    Ok(GenomicRegionHomologyTaskMessage::Done(result)) => {
                        done = Some(result);
                        break;
                    }
                    Err(TryRecvError::Empty) => break,
                    Err(TryRecvError::Disconnected) => {
                        disconnected = true;
                        break;
                    }
                }
            },
            Err(_) => {
                done = Some(Err(EngineError::new(
                    ErrorCode::Internal,
                    "Conservation worker result channel is unavailable",
                )));
            }
        }
        if disconnected && done.is_none() {
            done = Some(Err(EngineError::new(
                ErrorCode::Internal,
                "Conservation worker disconnected before returning a report",
            )));
        }
        let Some(result) = done else {
            ctx.request_repaint_after(Duration::from_millis(150));
            return;
        };
        self.genomic_region_conservation_task = None;
        match result {
            Ok(report) => {
                self.genomic_region_conservation_status = format!(
                    "Completed in {:.1}s: {} target(s), {} locus/loci, {} conserved block(s)",
                    started.elapsed().as_secs_f32(),
                    report.targets.len(),
                    report.loci.len(),
                    report.conserved_blocks.len()
                );
                self.genomic_region_conservation_selected_block_id = report
                    .conserved_blocks
                    .first()
                    .map(|block| block.block_id.clone());
                self.genomic_region_conservation_report = Some(Arc::new(report));
                self.genomic_region_conservation_module_assessment = None;
            }
            Err(error) => {
                self.genomic_region_conservation_status = format!(
                    "Conservation screen failed after {:.1}s: {}",
                    started.elapsed().as_secs_f32(),
                    error.message
                );
            }
        }
        ctx.request_repaint();
    }

    fn open_genomic_region_homology_json(&mut self) {
        let Some(path) = rfd::FileDialog::new()
            .add_filter("GENtle homology report", &["json"])
            .pick_file()
        else {
            self.genomic_region_conservation_status = "Open report canceled".to_string();
            return;
        };
        let loaded = std::fs::read(&path)
            .map_err(|error| format!("Could not read {}: {error}", path.display()))
            .and_then(|payload| {
                serde_json::from_slice::<gentle_protocol::GenomicRegionHomologyScreenReport>(
                    &payload,
                )
                .map_err(|error| format!("Could not parse {}: {error}", path.display()))
            })
            .and_then(|report| {
                crate::engine::validate_genomic_region_homology_report(&report)
                    .map_err(|error| error.message)?;
                if report.query.set_id != self.genomic_region_conservation_set_id
                    || report.query.region.region_id != self.genomic_region_conservation_region_id
                    || self
                        .genomic_region_conservation_expected_sha256
                        .as_deref()
                        .is_some_and(|expected| {
                            expected != report.query.region.content_sha256.as_str()
                        })
                {
                    return Err(
                        "The report is bound to a different or changed saved genomic region"
                            .to_string(),
                    );
                }
                Ok(report)
            });
        match loaded {
            Ok(report) => {
                self.genomic_region_conservation_selected_block_id = report
                    .conserved_blocks
                    .first()
                    .map(|block| block.block_id.clone());
                self.genomic_region_conservation_module_assessment = None;
                self.genomic_region_conservation_status = format!(
                    "Opened {}: {} target(s), {} locus/loci, {} conserved block(s)",
                    path.display(),
                    report.targets.len(),
                    report.loci.len(),
                    report.conserved_blocks.len()
                );
                self.genomic_region_conservation_report = Some(Arc::new(report));
            }
            Err(error) => {
                self.genomic_region_conservation_status = error;
            }
        }
    }

    fn export_genomic_region_homology_json(&mut self) {
        let Some(report) = self.genomic_region_conservation_report.as_ref() else {
            self.genomic_region_conservation_status =
                "Run the conservation screen before exporting".to_string();
            return;
        };
        let Some(path) = rfd::FileDialog::new()
            .set_file_name(format!("{}.homology.json", report.query.region.region_id))
            .save_file()
        else {
            self.genomic_region_conservation_status = "JSON export canceled".to_string();
            return;
        };
        match serde_json::to_vec_pretty(report.as_ref())
            .map_err(|error| error.to_string())
            .and_then(|payload| std::fs::write(&path, payload).map_err(|error| error.to_string()))
        {
            Ok(()) => {
                self.genomic_region_conservation_status = format!("Wrote {}", path.display());
            }
            Err(error) => {
                self.genomic_region_conservation_status =
                    format!("Could not export homology JSON: {error}");
            }
        }
    }

    fn export_genomic_region_homology_svg(&mut self) {
        let Some(report) = self.genomic_region_conservation_report.as_ref().cloned() else {
            self.genomic_region_conservation_status =
                "Run the conservation screen before exporting".to_string();
            return;
        };
        let Some(path) = rfd::FileDialog::new()
            .set_file_name(format!("{}.homology.svg", report.query.region.region_id))
            .save_file()
        else {
            self.genomic_region_conservation_status = "SVG export canceled".to_string();
            return;
        };
        if self
            .apply_operation_with_feedback_and_result(Operation::RenderGenomicRegionHomologySvg {
                report: Box::new(report.as_ref().clone()),
                path: path.display().to_string(),
            })
            .is_some()
        {
            self.genomic_region_conservation_status = format!("Wrote {}", path.display());
        }
    }

    fn save_selected_conserved_block(&mut self) {
        let Some(report) = self.genomic_region_conservation_report.as_ref().cloned() else {
            return;
        };
        let Some(block_id) = self
            .genomic_region_conservation_selected_block_id
            .as_deref()
        else {
            self.genomic_region_conservation_status = "Select a conserved block first".to_string();
            return;
        };
        let Some(block) = report
            .conserved_blocks
            .iter()
            .find(|block| block.block_id == block_id)
        else {
            self.genomic_region_conservation_status =
                "The selected conserved block is no longer present".to_string();
            return;
        };
        let mut local_projection = report.query.region.local_projection.clone();
        if let Some(projection) = local_projection.as_mut() {
            let block_start = block.query_start_0based as u64;
            let block_end = block.query_end_0based_exclusive as u64;
            if projection.local_strand == gentle_protocol::GenomicRegionStrand::Minus {
                let parent_end = projection.local_end_0based_exclusive;
                projection.local_start_0based = parent_end.saturating_sub(block_end);
                projection.local_end_0based_exclusive = parent_end.saturating_sub(block_start);
            } else {
                let parent_start = projection.local_start_0based;
                projection.local_start_0based = parent_start.saturating_add(block_start);
                projection.local_end_0based_exclusive = parent_start.saturating_add(block_end);
            }
        }
        let request = gentle_protocol::GenomicRegionCreateRequest {
            set_id: report.query.set_id.clone(),
            label: Some(format!(
                "Conserved {} block {}-{}",
                block.support_class.as_str(),
                block.query_start_0based.saturating_add(1),
                block.query_end_0based_exclusive
            )),
            description: Some(
                "Query-referenced exact-support block from a content-bound local homology screen"
                    .to_string(),
            ),
            interval: gentle_protocol::GenomicRegionInterval {
                reference: report.query.region.interval.reference.clone(),
                start_0based: block.genomic_start_0based,
                end_0based_exclusive: block.genomic_end_0based_exclusive,
                strand: block.genomic_strand,
                coordinate_convention:
                    gentle_protocol::GenomicRegionCoordinateConvention::ZeroBasedHalfOpen,
            },
            local_projection,
            purpose: gentle_protocol::GenomicRegionPurpose::ReporterCandidate,
            display_color_hex: Some("#15803D".to_string()),
            selection_method: gentle_protocol::GenomicRegionSelectionMethod::HomologyConservedBlock,
            evidence: vec![gentle_protocol::GenomicRegionEvidenceReference {
                evidence_id: format!("homology:{}", block.block_id),
                source_kind: "genomic_region_homology".to_string(),
                source_id: report.schema.clone(),
                source_sha256: Some(report.content_sha256.clone()),
                feature_or_window_id: Some(block.block_id.clone()),
                availability: gentle_protocol::GenomicRegionEvidenceAvailability::Available,
                evidence_statement: format!(
                    "Exact query support from {} genome(s) in the {} evidence class",
                    block.supporting_genome_ids.len(),
                    block.support_class.as_str()
                ),
                non_claims: report.non_claims.clone(),
                ..Default::default()
            }],
            notes: vec![
                "Saving this block preserves a reporter-fragment hypothesis; it does not assert autonomous regulatory function."
                    .to_string(),
            ],
            collision_policy: gentle_protocol::GenomicRegionCollisionPolicy::Reject,
            ..Default::default()
        };
        if self.apply_genomic_region_operation(Operation::CreateGenomicRegion { request }) {
            self.genomic_region_conservation_status =
                "Saved the selected block as a new portable genomic region".to_string();
        }
    }

    fn assess_selected_conservation_evidence(&mut self) {
        let Some(report) = self.genomic_region_conservation_report.as_ref().cloned() else {
            self.genomic_region_conservation_status =
                "Run the conservation screen before assessing modules".to_string();
            return;
        };
        if self
            .genomic_region_conservation_evidence_region_ids
            .is_empty()
        {
            self.genomic_region_conservation_status =
                "Select at least one saved evidence region inside the query".to_string();
            return;
        }
        let Some(set) = self.genomic_region_store_cache.as_ref().and_then(|store| {
            store
                .sets
                .iter()
                .find(|set| set.set_id == report.query.set_id)
        }) else {
            self.genomic_region_conservation_status =
                "The saved-region set is unavailable; refresh the region manager".to_string();
            return;
        };
        let query_interval = &report.query.region.interval;
        let mut spans = vec![];
        for region_id in &self.genomic_region_conservation_evidence_region_ids {
            let Some(region) = set
                .regions
                .iter()
                .find(|region| &region.region_id == region_id)
            else {
                self.genomic_region_conservation_status = format!(
                    "Selected evidence region '{region_id}' is no longer present; refresh and select again"
                );
                return;
            };
            if region.interval.reference != query_interval.reference
                || region.interval.start_0based < query_interval.start_0based
                || region.interval.end_0based_exclusive > query_interval.end_0based_exclusive
            {
                self.genomic_region_conservation_status = format!(
                    "Selected evidence region '{}' is not fully contained in the query",
                    region.region_id
                );
                return;
            }
            let (query_start, query_end) =
                if query_interval.strand == gentle_protocol::GenomicRegionStrand::Minus {
                    (
                        query_interval
                            .end_0based_exclusive
                            .saturating_sub(region.interval.end_0based_exclusive),
                        query_interval
                            .end_0based_exclusive
                            .saturating_sub(region.interval.start_0based),
                    )
                } else {
                    (
                        region
                            .interval
                            .start_0based
                            .saturating_sub(query_interval.start_0based),
                        region
                            .interval
                            .end_0based_exclusive
                            .saturating_sub(query_interval.start_0based),
                    )
                };
            let Ok(query_start_0based) = usize::try_from(query_start) else {
                self.genomic_region_conservation_status =
                    "Evidence coordinates exceed this platform's range".to_string();
                return;
            };
            let Ok(query_end_0based_exclusive) = usize::try_from(query_end) else {
                self.genomic_region_conservation_status =
                    "Evidence coordinates exceed this platform's range".to_string();
                return;
            };
            spans.push(gentle_protocol::PromoterModuleEvidenceSpan {
                evidence_id: format!("saved_region:{}:{}", set.set_id, region.region_id),
                evidence_kind: region.purpose.as_str().to_string(),
                query_start_0based,
                query_end_0based_exclusive,
                required: true,
                source_id: set.set_id.clone(),
                source_sha256: Some(region.content_sha256.clone()),
                evidence_statement: region
                    .description
                    .clone()
                    .or_else(|| region.label.clone())
                    .unwrap_or_else(|| "Caller-selected saved genomic evidence span".to_string()),
            });
        }
        let request = gentle_protocol::PromoterModuleAssessmentRequest {
            homology_report: Box::new(report.as_ref().clone()),
            selected_evidence_spans: spans,
            ..Default::default()
        };
        let Some(mut result) = self.apply_operation_with_feedback_and_result(
            Operation::AssessPromoterConservedModules {
                request,
                path: None,
            },
        ) else {
            self.genomic_region_conservation_status = self.op_status.clone();
            return;
        };
        let Some(assessment) = result.promoter_module_assessment.take() else {
            self.genomic_region_conservation_status =
                "Module assessment returned no typed report".to_string();
            return;
        };
        self.genomic_region_conservation_status = format!(
            "Module assessment: {} ({} selected block(s))",
            assessment.hypothesis.as_str(),
            assessment.selected_block_ids.len()
        );
        self.genomic_region_conservation_module_assessment = Some(Arc::new(*assessment));
    }

    pub(super) fn render_genomic_region_conservation_workspace(&mut self, ctx: &egui::Context) {
        if !self.show_genomic_region_conservation {
            return;
        }
        let report = self.genomic_region_conservation_report.clone();
        let module_assessment = self.genomic_region_conservation_module_assessment.clone();
        let evidence_regions = report
            .as_ref()
            .and_then(|report| {
                self.genomic_region_store_cache.as_ref().and_then(|store| {
                    store
                        .sets
                        .iter()
                        .find(|set| set.set_id == report.query.set_id)
                        .map(|set| {
                            set.regions
                                .iter()
                                .filter(|region| {
                                    region.region_id != report.query.region.region_id
                                        && region.interval.reference
                                            == report.query.region.interval.reference
                                        && region.interval.start_0based
                                            >= report.query.region.interval.start_0based
                                        && region.interval.end_0based_exclusive
                                            <= report.query.region.interval.end_0based_exclusive
                                })
                                .cloned()
                                .collect::<Vec<_>>()
                        })
                })
            })
            .unwrap_or_default();
        let running = self.genomic_region_conservation_task.is_some();
        let mut open = self.show_genomic_region_conservation;
        let mut run = false;
        let mut open_report = false;
        let mut export_json = false;
        let mut export_svg = false;
        let mut save_block = false;
        let mut assess_modules = false;
        let _subject_scope = crate::tutorial_gui_semantics::pseudonymous_subject_scope(&[
            self.seq_id.as_deref().unwrap_or("unnamed"),
            self.genomic_region_conservation_set_id.as_str(),
            self.genomic_region_conservation_region_id.as_str(),
        ]);
        let spec = crate::egui_compat::HostedWindowSpec::new(
            "Conservation",
            egui::Id::new(("genomic_region_conservation", self.panel_scope_key())),
            Vec2::new(1120.0, 780.0),
            Vec2::new(760.0, 520.0),
        );
        crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
            #[cfg(feature = "gui-test-support")]
            crate::gui_test_support::register_rect(
                ui.ctx().clone(),
                crate::tutorial_gui_semantics::WINDOW_REGION_CONSERVATION,
                crate::tutorial_gui_semantics::WINDOW_REGION_CONSERVATION,
                Some(&_subject_scope),
                crate::gui_test_support::GuiTestWidgetKind::Window,
                ui.max_rect(),
                true,
                true,
                true,
                Some(if running { "running" } else { "ready" }),
            );
            ui.horizontal_wrapped(|ui| {
                ui.heading("Genomic region conservation");
                ui.monospace(format!(
                    "{} / {}",
                    self.genomic_region_conservation_set_id,
                    self.genomic_region_conservation_region_id
                ));
            });
            ui.small(
                "Only validated local genomic BLAST indexes are searched. Similarity is shown separately for expected orthologs, unassigned cross-species loci, and same-genome alternatives.",
            );
            ui.horizontal_wrapped(|ui| {
                ui.label("Query genome");
                ui.add(
                    egui::TextEdit::singleline(
                        &mut self.genomic_region_conservation_query_genome_id,
                    )
                    .desired_width(180.0),
                );
                let run_response = ui.add_enabled(!running, egui::Button::new("Run local screen"));
                #[cfg(feature = "gui-test-support")]
                crate::gui_test_support::register_response(
                    &run_response,
                    crate::tutorial_gui_semantics::REGION_CONSERVATION_RUN,
                    crate::tutorial_gui_semantics::WINDOW_REGION_CONSERVATION,
                    Some(&_subject_scope),
                    crate::gui_test_support::GuiTestWidgetKind::Button,
                    false,
                );
                if run_response.clicked() {
                    run = true;
                }
                let open_report_response = ui.button("Open report...");
                #[cfg(feature = "gui-test-support")]
                crate::gui_test_support::register_response(
                    &open_report_response,
                    crate::tutorial_gui_semantics::REGION_CONSERVATION_OPEN_REPORT,
                    crate::tutorial_gui_semantics::WINDOW_REGION_CONSERVATION,
                    Some(&_subject_scope),
                    crate::gui_test_support::GuiTestWidgetKind::Button,
                    false,
                );
                if open_report_response.clicked() {
                    open_report = true;
                }
                let export_json_response =
                    ui.add_enabled(report.is_some(), egui::Button::new("Export JSON..."));
                #[cfg(feature = "gui-test-support")]
                crate::gui_test_support::register_response(
                    &export_json_response,
                    crate::tutorial_gui_semantics::REGION_CONSERVATION_EXPORT_JSON,
                    crate::tutorial_gui_semantics::WINDOW_REGION_CONSERVATION,
                    Some(&_subject_scope),
                    crate::gui_test_support::GuiTestWidgetKind::Button,
                    false,
                );
                if export_json_response.clicked() {
                    export_json = true;
                }
                let export_svg_response =
                    ui.add_enabled(report.is_some(), egui::Button::new("Export SVG..."));
                #[cfg(feature = "gui-test-support")]
                crate::gui_test_support::register_response(
                    &export_svg_response,
                    crate::tutorial_gui_semantics::REGION_CONSERVATION_EXPORT_SVG,
                    crate::tutorial_gui_semantics::WINDOW_REGION_CONSERVATION,
                    Some(&_subject_scope),
                    crate::gui_test_support::GuiTestWidgetKind::Button,
                    false,
                );
                if export_svg_response.clicked() {
                    export_svg = true;
                }
            });
            if running {
                ui.spinner();
            }
            if !self.genomic_region_conservation_status.trim().is_empty() {
                let _status = ui.small(&self.genomic_region_conservation_status);
                #[cfg(feature = "gui-test-support")]
                crate::gui_test_support::register_response(
                    &_status,
                    crate::tutorial_gui_semantics::REGION_CONSERVATION_STATUS,
                    crate::tutorial_gui_semantics::WINDOW_REGION_CONSERVATION,
                    Some(&_subject_scope),
                    crate::gui_test_support::GuiTestWidgetKind::Status,
                    false,
                );
            }
            ui.separator();
            let Some(report) = report.as_ref() else {
                ui.label(
                    "Run the screen to inspect local similarity. No index is downloaded or prepared automatically.",
                );
                return;
            };
            egui::ScrollArea::both()
                .id_salt(("genomic_region_conservation_scroll", self.panel_scope_key()))
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.set_min_width(1040.0);
                    ui.horizontal_wrapped(|ui| {
                        ui.strong(format!("{} bp query", report.query.sequence.len()));
                        ui.label(format!("{} target(s)", report.targets.len()));
                        ui.label(format!("{} retained locus/loci", report.loci.len()));
                        ui.label(format!("{} exact-support block(s)", report.conserved_blocks.len()));
                        ui.monospace(&report.content_sha256);
                    });
                    for warning in &report.warnings {
                        ui.colored_label(egui::Color32::from_rgb(180, 83, 9), warning);
                    }
                    egui::CollapsingHeader::new("Target readiness")
                        .default_open(!report.targets.is_empty())
                        .show(ui, |ui| {
                            egui::Grid::new("region_homology_target_grid")
                                .num_columns(5)
                                .striped(true)
                                .show(ui, |ui| {
                                    for heading in ["genome", "role", "status", "HSPs", "loci"] {
                                        ui.strong(heading);
                                    }
                                    ui.end_row();
                                    for target in &report.targets {
                                        ui.monospace(&target.target.genome_id);
                                        ui.label(target.target.role.as_str());
                                        ui.label(target.status.as_str());
                                        ui.monospace(format!(
                                            "{}/{}",
                                            target.accepted_hsp_count, target.raw_hsp_count
                                        ));
                                        ui.monospace(target.retained_locus_count.to_string());
                                        ui.end_row();
                                    }
                                });
                        });
                    ui.separator();
                    ui.strong("Conserved blocks");
                    if report.conserved_blocks.is_empty() {
                        ui.label("No exact-support block met the configured minimum length.");
                    }
                    for block in &report.conserved_blocks {
                        let selected = self
                            .genomic_region_conservation_selected_block_id
                            .as_deref()
                            == Some(block.block_id.as_str());
                        if ui
                            .selectable_label(
                                selected,
                                format!(
                                    "{}  {}-{}  {} bp  {}/{} available genome(s)",
                                    block.support_class.as_str(),
                                    block.query_start_0based.saturating_add(1),
                                    block.query_end_0based_exclusive,
                                    block
                                        .query_end_0based_exclusive
                                        .saturating_sub(block.query_start_0based),
                                    block.supporting_genome_ids.len(),
                                    block.available_genome_ids.len()
                                ),
                            )
                            .clicked()
                        {
                            self.genomic_region_conservation_selected_block_id =
                                Some(block.block_id.clone());
                        }
                    }
                    if let Some(block_id) = self
                        .genomic_region_conservation_selected_block_id
                        .as_deref()
                        && let Some(block) = report
                            .conserved_blocks
                            .iter()
                            .find(|block| block.block_id == block_id)
                    {
                        ui.horizontal_wrapped(|ui| {
                            ui.monospace(format!(
                                "{}:{}-{}",
                                report.query.region.interval.reference.contig_name,
                                block.genomic_start_0based.saturating_add(1),
                                block.genomic_end_0based_exclusive
                            ));
                            ui.label(format!("support {:.1}%", block.support_fraction * 100.0));
                            let save_response = ui.button("Save as portable region");
                            #[cfg(feature = "gui-test-support")]
                            crate::gui_test_support::register_response(
                                &save_response,
                                crate::tutorial_gui_semantics::REGION_CONSERVATION_SAVE_BLOCK,
                                crate::tutorial_gui_semantics::WINDOW_REGION_CONSERVATION,
                                Some(&_subject_scope),
                                crate::gui_test_support::GuiTestWidgetKind::Button,
                                false,
                            );
                            if save_response.clicked() {
                                save_block = true;
                            }
                        });
                    }
                    ui.separator();
                    ui.strong("Reporter-module hypothesis");
                    ui.small(
                        "Select saved evidence spans that the candidate fragment must retain. GENtle evaluates their geometry against explicitly supported ortholog blocks and same-genome similarity; it does not infer function from conservation.",
                    );
                    if evidence_regions.is_empty() {
                        ui.label(
                            "No other saved region is fully contained in this query. Save CUT&RUN, motif, Ensembl/SCREEN, or reporter-candidate spans first.",
                        );
                    } else {
                        for region in &evidence_regions {
                            let mut selected = self
                                .genomic_region_conservation_evidence_region_ids
                                .contains(&region.region_id);
                            let label = region.label.as_deref().unwrap_or(&region.region_id);
                            if ui
                                .checkbox(
                                    &mut selected,
                                    format!("{} ({})", label, region.purpose.as_str()),
                                )
                                .changed()
                            {
                                if selected {
                                    self.genomic_region_conservation_evidence_region_ids
                                        .insert(region.region_id.clone());
                                } else {
                                    self.genomic_region_conservation_evidence_region_ids
                                        .remove(&region.region_id);
                                }
                                self.genomic_region_conservation_module_assessment = None;
                            }
                        }
                    }
                    let assess_response = ui.add_enabled(
                        !self
                            .genomic_region_conservation_evidence_region_ids
                            .is_empty(),
                        egui::Button::new("Assess selected evidence"),
                    );
                    #[cfg(feature = "gui-test-support")]
                    crate::gui_test_support::register_response(
                        &assess_response,
                        crate::tutorial_gui_semantics::REGION_CONSERVATION_ASSESS_MODULES,
                        crate::tutorial_gui_semantics::WINDOW_REGION_CONSERVATION,
                        Some(&_subject_scope),
                        crate::gui_test_support::GuiTestWidgetKind::Button,
                        false,
                    );
                    if assess_response.clicked() {
                        assess_modules = true;
                    }
                    if let Some(assessment) = module_assessment.as_ref() {
                        ui.horizontal_wrapped(|ui| {
                            ui.label("Result");
                            ui.strong(assessment.hypothesis.as_str());
                            ui.monospace(&assessment.content_sha256);
                        });
                        for rule in &assessment.decision_trace {
                            ui.label(format!(
                                "{} {}: {}",
                                if rule.satisfied { "PASS" } else { "NO" },
                                rule.rule_id,
                                rule.detail
                            ));
                        }
                        for non_claim in &assessment.non_claims {
                            ui.small(egui::RichText::new(non_claim).italics());
                        }
                    }
                    ui.separator();
                    ui.strong("Query-referenced alignment");
                    ui.small(
                        ". exact match; letters are substitutions; - is a target deletion; blanks have no accepted HSP. Target insertions never add columns and remain in JSON provenance.",
                    );
                    let query_len = report.query.sequence.len();
                    for start in (0..query_len).step_by(100) {
                        let end = (start + 100).min(query_len);
                        ui.monospace(format!("query {}..{}", start + 1, end));
                        ui.horizontal(|ui| {
                            ui.monospace(format!("{:<28}", "QUERY"));
                            ui.monospace(&report.query.sequence[start..end]);
                        });
                        for row in report
                            .alignment_rows
                            .iter()
                            .filter(|row| {
                                row.locus_class
                                    != gentle_protocol::GenomicRegionHomologyLocusClass::Query
                            })
                            .take(50)
                        {
                            ui.horizontal(|ui| {
                                ui.monospace(format!(
                                    "{:<28}",
                                    format!("{}:{}", row.target_genome_id, row.subject_id)
                                ));
                                ui.monospace(&row.query_projection[start..end]);
                            });
                        }
                        ui.add_space(5.0);
                    }
                    if report.alignment_rows.len() > 51 {
                        ui.small(format!(
                            "{} additional retained rows are available in JSON/SVG.",
                            report.alignment_rows.len() - 51
                        ));
                    }
                    for non_claim in &report.non_claims {
                        ui.small(egui::RichText::new(non_claim).italics());
                    }
                });
        });
        self.show_genomic_region_conservation = open;
        if run {
            self.start_genomic_region_homology_screen();
        }
        if open_report {
            self.open_genomic_region_homology_json();
        }
        if export_json {
            self.export_genomic_region_homology_json();
        }
        if export_svg {
            self.export_genomic_region_homology_svg();
        }
        if save_block {
            self.save_selected_conserved_block();
        }
        if assess_modules {
            self.assess_selected_conservation_evidence();
        }
    }
}
