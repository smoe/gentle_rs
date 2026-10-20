//! Thin measured-gel forms and frozen-result presentation. Numerical work and
//! file validation stay in the shared engine; proposals require explicit review.

use super::*;

pub(super) struct DetectionIdentity {
    pub(super) draft: GelImageAnalysisRequest,
    pub(super) request: GelDetectRequest,
    pub(super) revision: u64,
}

pub(super) struct PendingSuggestions {
    pub(super) identity: DetectionIdentity,
    pub(super) suggestions: GelBandSuggestions,
}

pub(super) struct ComparisonInspection {
    revision: u64,
    mode: GelComparisonDisplay,
    report: Arc<GelComparisonReport>,
    #[cfg(test)]
    pub(super) svg: Arc<[u8]>,
    pub(super) pixels: Option<egui::ColorImage>,
    pub(super) texture: Option<egui::TextureHandle>,
    image_uri: String,
}

impl ComparisonInspection {
    /// Worker only: replay, project, rasterize and decode before GUI adoption.
    pub(super) fn new(
        revision: u64,
        report: Arc<GelComparisonReport>,
        mode: GelComparisonDisplay,
    ) -> Result<Self, String> {
        crate::gel_workbench::verify_comparison(&report).map_err(|e| e.message)?;
        let svg = gentle_render::gel_workbench::render_comparison_profile_svg(&report, mode)?;
        let pixels = Some(plot_pixels(&svg)?);
        let image_uri = format!(
            "gel-comparison-{}",
            crate::digest_utils::sha256_hex_bytes(svg.as_bytes())
        );
        Ok(Self {
            revision,
            mode,
            report,
            #[cfg(test)]
            svg: svg.into_bytes().into(),
            pixels,
            texture: None,
            image_uri,
        })
    }

    pub(super) fn install_texture(&mut self, ctx: &egui::Context) {
        if let Some(pixels) = self.pixels.take() {
            self.texture =
                Some(ctx.load_texture(&self.image_uri, pixels, egui::TextureOptions::LINEAR));
        }
    }
}

fn plot_pixels(svg: &str) -> Result<egui::ColorImage, String> {
    let rendered = crate::svg_png::render_svg_to_png_bytes(svg, Default::default())?;
    let rgba = image::load_from_memory_with_format(&rendered.bytes, image::ImageFormat::Png)
        .map_err(|e| e.to_string())?
        .into_rgba8();
    Ok(egui::ColorImage::from_rgba_unmultiplied(
        [rgba.width() as usize, rgba.height() as usize],
        rgba.as_raw(),
    ))
}

pub(super) struct AlgebraInspection {
    revision: u64,
    reference_profile_id: Option<String>,
    pub(super) presentation: GelAlgebraPresentation,
    #[cfg(test)]
    pub(super) svg: Arc<[u8]>,
    pub(super) pixels: Option<egui::ColorImage>,
    pub(super) texture: Option<egui::TextureHandle>,
    image_uri: String,
}

impl AlgebraInspection {
    pub(super) fn new(
        revision: u64,
        presentation: GelAlgebraPresentation,
        reference_profile_id: Option<String>,
    ) -> Result<Self, String> {
        let svg = gentle_render::gel_workbench::render_algebra_svg(&presentation);
        let pixels = if presentation.verification == GelAlgebraVerification::Verified {
            Some(plot_pixels(&svg)?)
        } else {
            None
        };
        let image_uri = format!(
            "gel-algebra-{}",
            crate::digest_utils::sha256_hex_bytes(svg.as_bytes())
        );
        Ok(Self {
            revision,
            reference_profile_id,
            presentation,
            #[cfg(test)]
            svg: svg.into_bytes().into(),
            pixels,
            texture: None,
            image_uri,
        })
    }

    pub(super) fn install_texture(&mut self, ctx: &egui::Context) {
        if let Some(pixels) = self.pixels.take() {
            self.texture =
                Some(ctx.load_texture(&self.image_uri, pixels, egui::TextureOptions::LINEAR));
        }
    }
}

pub(super) struct WorkbenchUi {
    pdf_object_selected: bool,
    pdf_object: u32,
    native_thumbnail: bool,
    pub(super) detection_lanes: BTreeSet<String>,
    pub(super) dark_bands: bool,
    pub(super) prominence: f64,
    distance: u32,
    smoothing: u32,
    pub(super) pending: Option<PendingSuggestions>,
    pub(super) selected_suggestions: BTreeSet<String>,
    pub(super) analysis_id: String,
    pub(super) lane_id: String,
    profile_id: String,
    pub(super) run_id: String,
    pub(super) extraction_dark: bool,
    pub(super) background: f64,
    pub(super) left: String,
    pub(super) right: String,
    pub(super) comparison_id: String,
    pub(super) range: [f64; 2],
    threshold_enabled: bool,
    threshold: f64,
    pub(super) coefficients: [f64; 2],
    pub(super) algebra_id: String,
    pub(super) saved_comparison: String,
    pub(super) comparison_display: GelComparisonDisplay,
    pub(super) comparison_inspection: Option<ComparisonInspection>,
    pub(super) selected_algebra: String,
    pub(super) algebra_reference: String,
    pub(super) algebra_inspection: Option<AlgebraInspection>,
    pub(super) export_analyses: Vec<String>,
    pub(super) export_comparisons: Vec<String>,
    pub(super) export_algebra: Vec<String>,
    export_title: String,
    export_path: String,
    advanced: String,
}

impl Default for WorkbenchUi {
    fn default() -> Self {
        Self {
            pdf_object_selected: false,
            pdf_object: 1,
            native_thumbnail: false,
            detection_lanes: BTreeSet::new(),
            dark_bands: false,
            prominence: 0.05,
            distance: 5,
            smoothing: 1,
            pending: None,
            selected_suggestions: BTreeSet::new(),
            analysis_id: String::new(),
            lane_id: String::new(),
            profile_id: "profile_1".into(),
            run_id: String::new(),
            extraction_dark: false,
            background: 0.0,
            left: String::new(),
            right: String::new(),
            comparison_id: "comparison_1".into(),
            range: [50.0, 5000.0],
            threshold_enabled: false,
            threshold: 200.0,
            coefficients: [1.0, -1.0],
            algebra_id: "residual_1".into(),
            saved_comparison: String::new(),
            comparison_display: GelComparisonDisplay::SignalDensityPerSize,
            comparison_inspection: None,
            selected_algebra: String::new(),
            algebra_reference: String::new(),
            algebra_inspection: None,
            export_analyses: Vec::new(),
            export_comparisons: Vec::new(),
            export_algebra: Vec::new(),
            export_title: "Measured gel comparison".into(),
            export_path: String::new(),
            advanced: String::new(),
        }
    }
}

#[derive(Default)]
pub(super) struct Inventory {
    pub(super) analyses: BTreeMap<String, Arc<GelImageAnalysisReport>>,
    pub(super) store: GelWorkbenchStore,
}

impl Inventory {
    fn read(shared: &Arc<RwLock<GentleEngine>>) -> Self {
        let engine = shared.read().unwrap();
        let gels = &engine.state().gel_images;
        Self {
            analyses: gels.analyses.clone(),
            store: gels.workbench.clone(),
        }
    }
}

fn valid_range(range: [f64; 2]) -> bool {
    range.iter().all(|x| x.is_finite() && *x > 0.0) && range[0] < range[1]
}

fn compatible(left: &GelProfile, right: &GelProfile) -> bool {
    left.size_kind == right.size_kind
        && left.evidence_kind == right.evidence_kind
        && left.signal_units == right.signal_units
        && left.measure == right.measure
}

impl WorkbenchUi {
    pub(super) fn comparison_matches(
        &self,
        inspection: &ComparisonInspection,
        revision: u64,
    ) -> bool {
        inspection.revision == revision
            && inspection.report.request.comparison_id == self.saved_comparison
            && inspection.mode == self.comparison_display
    }

    pub(super) fn invalidate_comparison(&mut self, revision: u64) {
        if self
            .comparison_inspection
            .as_ref()
            .is_some_and(|p| !self.comparison_matches(p, revision))
        {
            self.comparison_inspection = None;
        }
    }

    fn comparison_display_controls(&mut self, ui: &mut Ui) {
        ui.horizontal_wrapped(|ui| {
            ui.label("Comparison view / PDF (session only)");
            ui.radio_value(
                &mut self.comparison_display,
                GelComparisonDisplay::SignalDensityPerSize,
                "Signal per size width",
            );
            ui.radio_value(
                &mut self.comparison_display,
                GelComparisonDisplay::BinIntegral,
                "Legacy bin integrals",
            );
        });
        if self
            .comparison_inspection
            .as_ref()
            .is_some_and(|p| p.mode != self.comparison_display)
        {
            self.comparison_inspection = None;
        }
    }

    pub(super) fn algebra_matches(&self, inspection: &AlgebraInspection, revision: u64) -> bool {
        inspection.revision == revision
            && inspection.presentation.profile_id == self.selected_algebra
            && inspection.reference_profile_id.as_deref() == self.algebra_reference_id()
    }

    pub(super) fn algebra_reference_id(&self) -> Option<&str> {
        (!self.algebra_reference.is_empty()).then_some(self.algebra_reference.as_str())
    }

    pub(super) fn invalidate_algebra(&mut self, revision: u64) {
        if self
            .algebra_inspection
            .as_ref()
            .is_some_and(|p| !self.algebra_matches(p, revision))
        {
            self.algebra_inspection = None;
        }
    }

    pub(super) fn clear_suggestions(&mut self) {
        self.pending = None;
        self.selected_suggestions.clear();
    }

    fn pair<'a>(&self, inventory: &'a Inventory) -> Option<(&'a GelProfile, &'a GelProfile)> {
        Some((
            inventory.store.profiles.get(&self.left)?.as_ref(),
            inventory.store.profiles.get(&self.right)?.as_ref(),
        ))
    }

    pub(super) fn extraction_request(&self, inventory: &Inventory) -> Option<GelWorkbenchRequest> {
        let report = inventory.analyses.get(&self.analysis_id)?;
        if !report
            .request
            .lanes
            .iter()
            .any(|lane| lane.id == self.lane_id && lane.id != report.request.ladder.lane_id)
            || self.profile_id.trim().is_empty()
            || self.run_id.trim().is_empty()
            || inventory.store.profiles.contains_key(&self.profile_id)
            || !self.background.is_finite()
            || !(0.0..=1.0).contains(&self.background)
        {
            return None;
        }
        Some(GelWorkbenchRequest::ExtractProfile(
            GelExtractProfileRequest {
                profile_id: self.profile_id.clone(),
                analysis_id: self.analysis_id.clone(),
                lane_id: self.lane_id.clone(),
                run_id: self.run_id.clone(),
                dark_bands: self.extraction_dark,
                background: self.background,
                masked_ranges: vec![],
            },
        ))
    }

    pub(super) fn compare_request(&self, inventory: &Inventory) -> Option<GelWorkbenchRequest> {
        let (left, right) = self.pair(inventory)?;
        if !compatible(left, right)
            || left.evidence_kind == GelEvidenceKind::SignedResidual
            || !valid_range(self.range)
            || self.comparison_id.trim().is_empty()
            || inventory
                .store
                .comparisons
                .contains_key(&self.comparison_id)
            || (self.threshold_enabled
                && (!self.threshold.is_finite()
                    || self.threshold < self.range[0]
                    || self.threshold > self.range[1]))
        {
            return None;
        }
        Some(GelWorkbenchRequest::Compare(GelCompareRequest {
            comparison_id: self.comparison_id.clone(),
            left_profile: self.left.clone(),
            right_profile: self.right.clone(),
            range: self.range,
            short_threshold: self.threshold_enabled.then_some(self.threshold),
        }))
    }

    pub(super) fn algebra_request(&self, inventory: &Inventory) -> Option<GelWorkbenchRequest> {
        let (left, right) = self.pair(inventory)?;
        if !compatible(left, right)
            || !valid_range(self.range)
            || self.algebra_id.trim().is_empty()
            || inventory.store.profiles.contains_key(&self.algebra_id)
            || !self
                .coefficients
                .iter()
                .all(|c| c.is_finite() && c.abs() <= 1e6)
        {
            return None;
        }
        Some(GelWorkbenchRequest::Algebra(GelAlgebraRequest {
            output_id: self.algebra_id.clone(),
            range: self.range,
            terms: vec![
                GelAlgebraTerm {
                    profile_id: self.left.clone(),
                    coefficient: self.coefficients[0],
                },
                GelAlgebraTerm {
                    profile_id: self.right.clone(),
                    coefficient: self.coefficients[1],
                },
            ],
        }))
    }

    pub(super) fn export_request(&self, inventory: &Inventory) -> Option<GelWorkbenchRequest> {
        if self.export_title.trim().is_empty()
            || (self.export_analyses.is_empty()
                && self.export_comparisons.is_empty()
                && self.export_algebra.is_empty())
            || self.export_analyses.len() > 8
            || self.export_comparisons.len() > 8
            || self.export_algebra.len() > 8
            || (self.algebra_reference_id().is_some()
                && (self.export_algebra.len() != 1
                    || self.export_algebra[0] != self.selected_algebra
                    || !inventory
                        .store
                        .profiles
                        .get(&self.algebra_reference)
                        .is_some_and(|p| p.evidence_kind != GelEvidenceKind::SignedResidual)))
            || !self.export_algebra.iter().all(|id| {
                inventory
                    .store
                    .profiles
                    .get(id)
                    .is_some_and(|p| p.evidence_kind == GelEvidenceKind::SignedResidual)
            })
            || !self
                .export_analyses
                .iter()
                .all(|id| inventory.analyses.contains_key(id))
            || !self
                .export_comparisons
                .iter()
                .all(|id| inventory.store.comparisons.contains_key(id))
        {
            return None;
        }
        Some(GelWorkbenchRequest::Export(GelWorkbenchExportRequest {
            path: self.export_path.clone(),
            analysis_ids: self.export_analyses.clone(),
            comparison_ids: self.export_comparisons.clone(),
            algebra_ids: self.export_algebra.clone(),
            comparison_display: (!self.export_comparisons.is_empty())
                .then_some(self.comparison_display),
            algebra_reference_profile_id: self.algebra_reference_id().map(str::to_owned),
            title: self.export_title.clone(),
        }))
    }
}

/// Dialogs and extension validation run on the same worker as engine I/O.
pub(super) fn prepare_file_request(request: &mut GelWorkbenchRequest) -> Result<(), String> {
    match request {
        GelWorkbenchRequest::ImportSource(source) => {
            if source.path.trim().is_empty() {
                source.path = rfd::FileDialog::new()
                    .set_title("Import PDF raster / native gel thumbnail")
                    .add_filter("Gel source", &["pdf", "D5000", "d5000"])
                    .pick_file()
                    .ok_or("Source selection canceled")?
                    .display()
                    .to_string();
            }
            let extension = std::path::Path::new(&source.path)
                .extension()
                .and_then(|x| x.to_str())
                .unwrap_or("")
                .to_ascii_lowercase();
            if extension == "pdf" && source.pdf_object.is_none_or(|object| object == 0) {
                return Err("Select an explicit PDF raster object number before importing a PDF. No image is chosen by page order.".into());
            }
            if extension == "d5000" && !source.allow_native_thumbnail {
                return Err("D5000 requires explicit display-only thumbnail consent. Encrypted numerical data and complete raw import are unsupported.".into());
            }
        }
        GelWorkbenchRequest::ImportProfileFile(source) => {
            if source.path.trim().is_empty() {
                source.path = rfd::FileDialog::new()
                    .set_title("Import canonical integrated-bin CSV")
                    .add_filter("Integrated-bin CSV", &["csv"])
                    .pick_file()
                    .ok_or("Numeric source selection canceled")?
                    .display()
                    .to_string();
            }
        }
        GelWorkbenchRequest::Export(export) => {
            if export.path.trim().is_empty() {
                export.path = rfd::FileDialog::new()
                    .set_title("Export saved gel report (new file only)")
                    .set_file_name("measured-gel-report.pdf")
                    .add_filter("PDF report", &["pdf"])
                    .save_file()
                    .ok_or("Report export selection canceled")?
                    .display()
                    .to_string();
            }
            let extension = std::path::Path::new(&export.path)
                .extension()
                .and_then(|x| x.to_str())
                .unwrap_or("");
            if !extension.eq_ignore_ascii_case("pdf") {
                return Err("Choose a new .pdf path for the workbench report. Existing files are never overwritten.".into());
            }
        }
        _ => {}
    }
    Ok(())
}

impl GelImageEditor {
    pub(super) fn detection_request(&self) -> Option<GelDetectRequest> {
        let draft = self.draft.as_ref()?;
        let lanes: Vec<_> = draft
            .lanes
            .iter()
            .filter(|lane| {
                lane.id != draft.ladder.lane_id && self.workbench.detection_lanes.contains(&lane.id)
            })
            .cloned()
            .collect();
        if lanes.is_empty() {
            return None;
        }
        Some(GelDetectRequest {
            image_id: draft.image_id.clone(),
            lanes,
            migration: draft.migration,
            dark_bands: self.workbench.dark_bands,
            minimum_prominence: self.workbench.prominence,
            minimum_distance_px: self.workbench.distance,
            smoothing_radius_px: self.workbench.smoothing,
        })
    }

    pub(super) fn suggestions_match(&self, pending: &PendingSuggestions, revision: u64) -> bool {
        pending.identity.revision == revision
            && self.draft.as_ref() == Some(&pending.identity.draft)
            && self.detection_request().as_ref() == Some(&pending.identity.request)
            && pending.suggestions.request == pending.identity.request
            && pending.suggestions.image_sha256 == pending.identity.draft.image_sha256
            && self.image.as_ref().is_some_and(|image| {
                image.descriptor.image_id == pending.identity.draft.image_id
                    && image.descriptor.sha256 == pending.suggestions.image_sha256
            })
    }

    pub(super) fn accept_suggestions(&mut self, shared: &Arc<RwLock<GentleEngine>>) {
        self.bind_project(shared);
        self.sync_draft(shared);
        let mut engine = shared.write().unwrap();
        let Some(pending) = &self.workbench.pending else {
            return;
        };
        if !self.suggestions_match(pending, engine.mutation_revision()) {
            self.workbench.clear_suggestions();
            self.status =
                "Draft or request changed; proposals discarded without applying marks.".into();
            return;
        }
        let selected: Vec<_> = pending
            .suggestions
            .bands
            .iter()
            .filter(|b| self.workbench.selected_suggestions.contains(&b.id))
            .cloned()
            .collect();
        if selected.is_empty() {
            return;
        }
        let mut draft = self.draft.clone().unwrap();
        if draft.sample_bands.len().saturating_add(selected.len()) > 2048 {
            self.status = "Accept fewer proposals: at most 2048 sample marks are supported.".into();
            return;
        }
        for mut band in selected {
            if band.lane_id == draft.ladder.lane_id
                || !pending
                    .identity
                    .request
                    .lanes
                    .iter()
                    .any(|lane| lane.id == band.lane_id && contains(lane, band.center))
            {
                self.status =
                    "Proposal does not belong to a selected sample rectangle; nothing applied."
                        .into();
                return;
            }
            if draft
                .sample_bands
                .iter()
                .any(|b| b.lane_id == band.lane_id && b.center == band.center)
            {
                continue;
            }
            band.id = fresh_id("sample", draft.sample_bands.iter().map(|b| b.id.clone()));
            band.label = format!("Reviewed peak {}", band.id);
            draft.sample_bands.push(band);
        }
        match engine.apply(Operation::SaveGelImageDraft {
            request: Arc::new(draft.clone()),
        }) {
            Ok(_) => {
                self.draft = Some(draft);
                self.stored_draft = engine
                    .state()
                    .gel_images
                    .drafts
                    .get(&self.image_id)
                    .cloned();
                self.workbench.clear_suggestions();
                self.status = "Reviewed sample marks saved as a draft. Ladder unchanged; Analyze to validate sizing.".into();
            }
            Err(error) => self.status = format!("No proposal marks saved: {}", error.message),
        }
    }

    pub(super) fn native_source_controls(
        &mut self,
        ui: &mut Ui,
        shared: &Arc<RwLock<GentleEngine>>,
    ) {
        ui.collapsing("PDF / D5000 source import", |ui| {
            ui.label("Uses the image path above, or opens a file dialog. PDF import extracts an embedded raster, not a rendered page.");
            ui.horizontal_wrapped(|ui| {
                ui.checkbox(&mut self.workbench.pdf_object_selected, "Select PDF raster object number");
                ui.add_enabled(self.workbench.pdf_object_selected, egui::DragValue::new(&mut self.workbench.pdf_object).range(1..=u32::MAX));
            });
            ui.checkbox(&mut self.workbench.native_thumbnail, "Allow native D5000 thumbnail (display-only)");
            ui.small("Encrypted numerical traces are unsupported. A thumbnail is not complete raw-data import, detector fluorescence, concentration or yield.");
            if ui.button("Import PDF / native source...").clicked() {
                let image_id = fresh_id("gel", shared.read().unwrap().state().gel_images.images.keys().cloned());
                self.start(shared, ui.ctx(), GelJob::Workbench(GelWorkbenchRequest::ImportSource(GelSourceImportRequest {
                    image_id, path: self.import_path.clone(),
                    pdf_object: self.workbench.pdf_object_selected.then_some(self.workbench.pdf_object),
                    allow_native_thumbnail: self.workbench.native_thumbnail,
                })));
            }
        });
        let source = shared
            .read()
            .unwrap()
            .state()
            .gel_images
            .workbench
            .sources
            .get(&self.image_id)
            .cloned();
        if let Some(source) = source {
            ui.collapsing("Imported source provenance", |ui| {
                ui.label(format!("{} | {}", source.source_name, source.extraction));
                ui.small(format!("Source SHA-256: {}", source.source_sha256));
                ui.small(format!(
                    "Extracted image SHA-256: {}",
                    source.extracted_image_sha256
                ));
                for warning in &source.warnings {
                    ui.small(warning);
                }
            });
        }
    }

    pub(super) fn workbench_controls(&mut self, ui: &mut Ui, shared: &Arc<RwLock<GentleEngine>>) {
        let inventory = Inventory::read(shared);
        ui.separator();
        ui.heading("Measured profiles and saved comparisons");
        ui.small("Explicit engine requests only. Saved results are inspected below, never re-extracted or rescored by repainting this window.");
        self.detection_controls(ui, shared);
        self.extraction_controls(ui, shared, &inventory);
        self.comparison_controls(ui, shared, &inventory);
        self.saved_comparison_controls(ui, shared, &inventory);
        self.algebra_inspection_controls(ui, shared, &inventory);
        self.workbench_export_controls(ui, shared, &inventory);
        self.advanced_workbench_controls(ui, shared);
    }

    fn detection_controls(&mut self, ui: &mut Ui, shared: &Arc<RwLock<GentleEngine>>) {
        ui.collapsing("Propose sample bands in manual lane rectangles", |ui| {
            ui.label("Select the sample rectangles to scan. The confirmed ladder lane is excluded; ladder sizes are never autoassigned.");
            if let Some(draft) = &self.draft {
                for lane in draft.lanes.iter().filter(|lane| lane.id != draft.ladder.lane_id) {
                    let mut selected = self.workbench.detection_lanes.contains(&lane.id);
                    if ui.checkbox(&mut selected, format!("{}: {}", lane.id, lane.label)).changed() {
                        if selected { self.workbench.detection_lanes.insert(lane.id.clone()); }
                        else { self.workbench.detection_lanes.remove(&lane.id); }
                    }
                }
            }
            polarity_controls(ui, &mut self.workbench.dark_bands);
            ui.horizontal_wrapped(|ui| {
                ui.add(egui::DragValue::new(&mut self.workbench.prominence).speed(0.01).range(0.001..=1.0).prefix("Prominence "));
                ui.add(egui::DragValue::new(&mut self.workbench.distance).range(1..=1000).prefix("Min distance px "));
                ui.add(egui::DragValue::new(&mut self.workbench.smoothing).range(0..=32).prefix("Smoothing radius px "));
            });
            let request = self.detection_request();
            if ui.add_enabled(request.is_some(), egui::Button::new("Detect sample band proposals")).clicked()
                && let (Some(request), Some(draft)) = (request, self.draft.clone()) {
                self.workbench.clear_suggestions();
                self.start(shared, ui.ctx(), GelJob::Detect { request, draft: Box::new(draft) });
            }
            let revision = shared.read().unwrap().mutation_revision();
            if self.workbench.pending.as_ref().is_some_and(|pending| !self.suggestions_match(pending, revision)) {
                self.workbench.clear_suggestions();
                self.status = "Draft or detector settings changed; old proposals discarded.".into();
            }
            if let Some(pending) = &self.workbench.pending {
                ui.label(format!("{} unconfirmed proposals; select those reviewed on the image", pending.suggestions.bands.len()));
                egui::ScrollArea::vertical().id_salt("gel_peak_proposals").max_height(180.0).show(ui, |ui| {
                    for band in &pending.suggestions.bands {
                        let mut selected = self.workbench.selected_suggestions.contains(&band.id);
                        if ui.checkbox(&mut selected, format!("{} / {}: ({:.1}, {:.1})", band.lane_id, band.id, band.center.x, band.center.y)).changed() {
                            if selected { self.workbench.selected_suggestions.insert(band.id.clone()); }
                            else { self.workbench.selected_suggestions.remove(&band.id); }
                        }
                    }
                });
                for warning in &pending.suggestions.warnings { ui.small(warning); }
                if ui.add_enabled(!self.workbench.selected_suggestions.is_empty(), egui::Button::new("Accept selected proposals as sample marks")).clicked() {
                    self.accept_suggestions(shared);
                }
                if ui.button("Discard proposals").clicked() { self.workbench.clear_suggestions(); }
            }
        });
    }

    fn extraction_controls(
        &mut self,
        ui: &mut Ui,
        shared: &Arc<RwLock<GentleEngine>>,
        inventory: &Inventory,
    ) {
        ui.collapsing("Extract a profile from saved calibration", |ui| {
            ui.label("Uses the chosen immutable report's geometry and calibration, not the current editable draft.");
            let previous = self.workbench.analysis_id.clone();
            select_id(ui, "gel_extract_analysis", &mut self.workbench.analysis_id, "Choose saved calibration", inventory.analyses.iter().map(|(id, r)| (id.clone(), format!("{id} | image {}", r.image.image_id))));
            if previous != self.workbench.analysis_id { self.workbench.lane_id.clear(); }
            let report = inventory.analyses.get(&self.workbench.analysis_id);
            let lanes = report.into_iter().flat_map(|r| r.request.lanes.iter().filter(|lane| lane.id != r.request.ladder.lane_id))
                .map(|lane| (lane.id.clone(), format!("{}: {}", lane.id, lane.label))).collect::<Vec<_>>();
            select_id(ui, "gel_extract_lane", &mut self.workbench.lane_id, "Choose sample lane", lanes);
            ui.horizontal_wrapped(|ui| {
                ui.label("New profile ID"); ui.text_edit_singleline(&mut self.workbench.profile_id);
                ui.label("Run label / ID"); ui.text_edit_singleline(&mut self.workbench.run_id);
            });
            polarity_controls(ui, &mut self.workbench.extraction_dark);
            ui.add(egui::DragValue::new(&mut self.workbench.background).range(0.0..=1.0).speed(0.01).prefix("Background (0..1, after polarity) "));
            ui.small("Background is user-specified in polarity-adjusted display intensity. Extraction subtracts it and floors individual pixel signals at zero. No background is fitted. This is image-weighted shape, not detector fluorescence or abundance.");
            ui.small("This form applies no size masks. To exclude reviewed marker intervals, use ExtractProfile masked_ranges in Advanced JSON; masked pixel bins remain missing, not zero.");
            let request = self.workbench.extraction_request(inventory);
            if ui.add_enabled(request.is_some(), egui::Button::new("Extract and save profile")).clicked()
                && let Some(request) = request {
                self.start(shared, ui.ctx(), GelJob::Workbench(request));
            }
        });
    }

    fn comparison_controls(
        &mut self,
        ui: &mut Ui,
        shared: &Arc<RwLock<GentleEngine>>,
        inventory: &Inventory,
    ) {
        ui.group(|ui| {
            ui.label("Compare any two saved profiles, including different images / runs");
            select_id(ui, "gel_profile_left", &mut self.workbench.left, "Choose left profile", profile_choices(inventory));
            select_id(ui, "gel_profile_right", &mut self.workbench.right, "Choose right profile", profile_choices(inventory));
            if inventory.store.profiles.is_empty() { ui.small("No saved profiles. Extract one above or import canonical numeric bins in Advanced."); }
            if let Some((left, right)) = self.workbench.pair(inventory) {
                ui.small(format!("Size: {} / {}; signal: {} / {}", left.size_kind.unit(), right.size_kind.unit(), left.signal_units, right.signal_units));
                if !compatible(left, right) { ui.colored_label(REFERENCE_COLOR, "Incompatible evidence kind, size kind, signal units or integration measure. No conversion is inferred."); }
            }
            ui.horizontal_wrapped(|ui| {
                ui.label("Size range (chosen profiles' bp or kDa)");
                ui.add(egui::DragValue::new(&mut self.workbench.range[0]).speed(1.0).prefix("From "));
                ui.add(egui::DragValue::new(&mut self.workbench.range[1]).speed(1.0).prefix("To "));
                ui.checkbox(&mut self.workbench.threshold_enabled, "Short-size threshold");
                ui.add_enabled(self.workbench.threshold_enabled, egui::DragValue::new(&mut self.workbench.threshold));
            });
            ui.horizontal_wrapped(|ui| {
                ui.label("New comparison ID"); ui.text_edit_singleline(&mut self.workbench.comparison_id);
                let request = self.workbench.compare_request(inventory);
                if ui.add_enabled(request.is_some(), egui::Button::new("Compare and save result")).clicked() && let Some(request) = request {
                    self.start(shared, ui.ctx(), GelJob::Workbench(request));
                }
            });
            ui.small("Requires compatible non-residual profiles, a new ID, and a positive increasing range. Coverage/weight validity is checked by the engine. p-value is not computed.");
            ui.collapsing("Signed profile algebra (explicit coefficients)", |ui| {
                ui.label("Uses the same left/right selections and size range. Values remain signed, unnormalized and unclipped.");
                ui.horizontal_wrapped(|ui| {
                    ui.add(egui::DragValue::new(&mut self.workbench.coefficients[0]).speed(0.1).prefix("Left coefficient "));
                    ui.add(egui::DragValue::new(&mut self.workbench.coefficients[1]).speed(0.1).prefix("Right coefficient "));
                    ui.label("New output profile ID"); ui.text_edit_singleline(&mut self.workbench.algebra_id);
                });
                ui.label(format!("Output = {} x left + {} x right", self.workbench.coefficients[0], self.workbench.coefficients[1]));
                let request = self.workbench.algebra_request(inventory);
                if ui.add_enabled(request.is_some(), egui::Button::new("Save signed linear combination")).clicked() && let Some(request) = request {
                    self.start(shared, ui.ctx(), GelJob::Workbench(request));
                }
                ui.small("Negative residuals and missing coverage are retained. A signed profile is not a probability distribution or gene-specific enrichment claim.");
            });
        });
    }

    fn saved_comparison_controls(
        &mut self,
        ui: &mut Ui,
        shared: &Arc<RwLock<GentleEngine>>,
        inventory: &Inventory,
    ) {
        select_id(
            ui,
            "gel_saved_comparison",
            &mut self.workbench.saved_comparison,
            "Inspect saved comparison",
            inventory
                .store
                .comparisons
                .keys()
                .map(|id| (id.clone(), id.clone())),
        );
        self.workbench.comparison_display_controls(ui);
        self.workbench
            .invalidate_comparison(shared.read().unwrap().mutation_revision());
        if ui
            .add_enabled(
                self.job.is_none()
                    && inventory
                        .store
                        .comparisons
                        .contains_key(&self.workbench.saved_comparison),
                egui::Button::new("Verify and display saved comparison"),
            )
            .clicked()
        {
            self.workbench.comparison_inspection = None;
            self.start(
                shared,
                ui.ctx(),
                GelJob::ComparisonPlot {
                    comparison_id: self.workbench.saved_comparison.clone(),
                    mode: self.workbench.comparison_display,
                },
            );
        }
        if let Some(cached) = &self.workbench.comparison_inspection {
            ui.small("Replay verified: internal consistency, not scientific acceptance. No new extraction or statistics.");
            if let Some(texture) = &cached.texture {
                ui.add(
                    egui::Image::new(texture)
                        .max_width(ui.available_width())
                        .shrink_to_fit(),
                );
            }
            comparison_summary(ui, &cached.report);
        } else if !self.workbench.saved_comparison.is_empty() {
            ui.small("Comparison not prepared for the selected mode. Verify to show curves, summaries and QQ; any earlier-mode result is discarded.");
        }
        ui.collapsing("Saved profile provenance / signed outputs", |ui| {
            for profile in inventory.store.profiles.values() {
                ui.push_id(&profile.id, |ui| {
                    ui.collapsing(format!("{}: {}", profile.id, profile.label), |ui| {
                        profile_provenance(ui, profile);
                        ui.label(format!("{} integrated bins ({} missing)", profile.weights.len(), profile.weights.iter().filter(|w| w.is_none()).count()));
                        if profile.evidence_kind == GelEvidenceKind::SignedResidual {
                            let minimum = profile.weights.iter().flatten().copied().reduce(f64::min);
                            let maximum = profile.weights.iter().flatten().copied().reduce(f64::max);
                            ui.label(format!("Signed bin weights: {minimum:?} to {maximum:?}; not clipped or normalized"));
                        }
                        if let Some(source) = inventory.store.numeric_sources.get(&profile.id) {
                            ui.small(format!("Retained numeric source: {} | SHA-256 {}", source.source_name, source.sha256));
                        }
                    });
                });
            }
        });
    }

    fn algebra_inspection_controls(
        &mut self,
        ui: &mut Ui,
        shared: &Arc<RwLock<GentleEngine>>,
        inventory: &Inventory,
    ) {
        ui.separator();
        ui.label("Inspect saved signed algebra (no CDF / QQ / significance)");
        select_id(
            ui,
            "gel_algebra_inspect",
            &mut self.workbench.selected_algebra,
            "Select a derived result",
            inventory
                .store
                .profiles
                .values()
                .filter(|p| p.evidence_kind == GelEvidenceKind::SignedResidual)
                .map(|p| (p.id.clone(), p.id.clone())),
        );
        ui.label(
            "Optional custom-only reference (caller-selected; separate Y axis, no arithmetic)",
        );
        select_id(
            ui,
            "gel_algebra_reference",
            &mut self.workbench.algebra_reference,
            "None selected",
            std::iter::once((String::new(), "None (no reference panel)".into())).chain(
                inventory
                    .store
                    .profiles
                    .values()
                    .filter(|p| p.evidence_kind != GelEvidenceKind::SignedResidual)
                    .map(|p| {
                        (
                            p.id.clone(),
                            format!(
                                "{} | {} | run {} / lane {}",
                                p.id, p.label, p.run_id, p.lane_id
                            ),
                        )
                    }),
            ),
        );
        ui.small("Unscaled native-bin signal per size width. Same size axis, independently labelled Y axis; not comparable amplitudes or verified custom priming. No coefficient is selected or fitted.");
        self.workbench
            .invalidate_algebra(shared.read().unwrap().mutation_revision());
        let available = inventory
            .store
            .profiles
            .get(&self.workbench.selected_algebra)
            .is_some_and(|p| p.evidence_kind == GelEvidenceKind::SignedResidual);
        if ui
            .add_enabled(
                available,
                egui::Button::new("Verify and display saved algebra"),
            )
            .clicked()
        {
            self.workbench.algebra_inspection = None;
            self.start(
                shared,
                ui.ctx(),
                GelJob::Workbench(GelWorkbenchRequest::InspectAlgebra {
                    profile_id: self.workbench.selected_algebra.clone(),
                    reference_profile_id: self.workbench.algebra_reference_id().map(str::to_owned),
                }),
            );
        }
        if let Some(cached) = &self.workbench.algebra_inspection {
            let report = &cached.presentation;
            ui.label(format!(
                "Verification: {:?} (internal consistency, not scientific acceptance)",
                report.verification
            ));
            if let Some(reason) = &report.reason {
                ui.label(reason);
            }
            if report.verification == GelAlgebraVerification::Verified {
                if let Some(texture) = &cached.texture {
                    ui.add(
                        egui::Image::new(texture)
                            .max_width(ui.available_width())
                            .shrink_to_fit(),
                    );
                }
                ui.collapsing(
                    "Complete recorded inputs, coefficients and warnings",
                    |ui| {
                        if let Some(request) = &report.request {
                            for term in &request.terms {
                                ui.label(format!("{} x {}", term.coefficient, term.profile_id));
                            }
                        }
                        ui.small(format!(
                            "Algebra input digest: {}",
                            report.input_digest.as_deref().unwrap_or("unavailable")
                        ));
                        for input in &report.inputs {
                            ui.push_id(&input.id, |ui| {
                                ui.collapsing(&input.id, |ui| profile_provenance(ui, input));
                            });
                        }
                        if let Some(reference) = &report.reference {
                            ui.label("Custom-only reference (caller-selected, not an operand)");
                            profile_provenance(ui, &reference.profile);
                        }
                        for warning in &report.warnings {
                            ui.small(warning);
                        }
                    },
                );
            } else {
                ui.small("No verified plot or report export. Raw legacy fields remain available under saved profile provenance.");
            }
        }
    }

    fn workbench_export_controls(
        &mut self,
        ui: &mut Ui,
        shared: &Arc<RwLock<GentleEngine>>,
        inventory: &Inventory,
    ) {
        ui.collapsing("Export selected saved evidence as PDF", |ui| {
            ui.label("Select saved analyses in explicit before / after order. Select saved comparisons separately. No adjacent-lane or filename pairing is inferred.");
            ordered_selection(ui, "gel_pdf_analyses", &mut self.workbench.export_analyses,
                inventory.analyses.keys().cloned().collect(), true);
            ui.label("Saved comparison IDs (export order)");
            ordered_selection(ui, "gel_pdf_comparisons", &mut self.workbench.export_comparisons,
                inventory.store.comparisons.keys().cloned().collect(), false);
            self.workbench.comparison_display_controls(ui);
            ui.label("Saved algebra IDs (append after existing pages; verification required)");
            ordered_selection(ui, "gel_pdf_algebra", &mut self.workbench.export_algebra,
                inventory.store.profiles.values().filter(|p| p.evidence_kind == GelEvidenceKind::SignedResidual).map(|p| p.id.clone()).collect(), false);
            if self.workbench.algebra_reference_id().is_some() {
                ui.small("The selected custom-only reference is exported only with exactly the inspected algebra ID selected. Clear the reference to export other/multiple algebra results. No reference pairing is inferred.");
            }
            ui.horizontal_wrapped(|ui| {
                ui.label("Report title"); ui.text_edit_singleline(&mut self.workbench.export_title);
                ui.label("New PDF path (optional)"); ui.text_edit_singleline(&mut self.workbench.export_path);
            });
            let request = self.workbench.export_request(inventory);
            if ui.add_enabled(request.is_some(), egui::Button::new("Export stored report PDF...")).clicked() && let Some(request) = request {
                self.start(shared, ui.ctx(), GelJob::Workbench(request));
            }
            ui.small("Up to 8 of each record kind. Algebra export refuses mismatching/unverifiable inputs and excessive provenance pages. A file dialog runs in the background if no path is given. Existing files are refused. Export is external: no project edit or undo; closing does not cancel it.");
        });
    }

    fn advanced_workbench_controls(&mut self, ui: &mut Ui, shared: &Arc<RwLock<GentleEngine>>) {
        ui.collapsing("Advanced: typed workbench JSON / numeric import", |ui| {
            ui.label("One GelWorkbenchRequest JSON object. Numeric CSV: action=import_profile_file; columns left_edge,right_edge,weight; empty weights mean missing. Weights must already be integrated bins, never guessed point-sample integration.");
            ui.horizontal(|ui| {
                if ui.button("Canonical CSV request template").clicked() {
                    self.workbench.advanced = serde_json::to_string_pretty(&GelWorkbenchRequest::ImportProfileFile(GelProfileFileRequest {
                        path: String::new(), profile_id: "numeric_profile_1".into(), label: "Numeric integrated bins".into(),
                        run_id: String::new(), lane_id: String::new(), size_kind: GelSizeKind::LinearDnaBp,
                        signal_units: String::new(), measure: String::new(),
                    })).unwrap_or_default();
                }
            });
            ui.small("Supply the run, lane, units and integration measure explicitly. Original CSV bytes are retained in the project; caller-declared calibration is not independently verified. Imports are not complete encrypted D5000 decoding.");
            ui.add(egui::TextEdit::multiline(&mut self.workbench.advanced).desired_rows(8).desired_width(f32::INFINITY));
            if ui.add_enabled(!self.workbench.advanced.trim().is_empty(), egui::Button::new("Run explicit JSON request")).clicked() {
                match serde_json::from_str::<GelWorkbenchRequest>(&self.workbench.advanced) {
                    Ok(GelWorkbenchRequest::DetectBands(request)) => {
                        if self.detection_request().as_ref() == Some(&request) && let Some(draft) = self.draft.clone() {
                            self.workbench.clear_suggestions();
                            self.start(shared, ui.ctx(), GelJob::Detect { request, draft: Box::new(draft) });
                        } else {
                            self.status = "Detection JSON must match the selected manual sample rectangles and current detector controls; ladder proposals are not allowed.".into();
                        }
                    }
                    Ok(request) => {
                        if let GelWorkbenchRequest::InspectAlgebra { profile_id, reference_profile_id } = &request {
                            self.workbench.selected_algebra = profile_id.clone();
                            self.workbench.algebra_reference = reference_profile_id.clone().unwrap_or_default();
                            self.workbench.algebra_inspection = None;
                        }
                        self.start(shared, ui.ctx(), GelJob::Workbench(request));
                    },
                    Err(error) => self.status = format!("Request not submitted: {error}"),
                }
            }
        });
    }
}

fn polarity_controls(ui: &mut Ui, dark: &mut bool) {
    ui.horizontal_wrapped(|ui| {
        ui.label("Signal polarity");
        ui.radio_value(dark, false, "Bright bands on dark background");
        ui.radio_value(dark, true, "Dark bands on bright background");
    });
}

fn select_id(
    ui: &mut Ui,
    salt: &str,
    selected: &mut String,
    placeholder: &str,
    choices: impl IntoIterator<Item = (String, String)>,
) {
    let choices: Vec<_> = choices.into_iter().collect();
    let text = choices
        .iter()
        .find(|(id, _)| id.as_str() == selected.as_str())
        .map(|(_, label)| label.as_str())
        .unwrap_or(placeholder);
    ui.add_enabled_ui(!choices.is_empty(), |ui| {
        egui::ComboBox::from_id_salt(salt)
            .selected_text(text)
            .show_ui(ui, |ui| {
                for (id, label) in &choices {
                    ui.selectable_value(selected, id.clone(), label);
                }
            });
    });
}

fn profile_choices(inventory: &Inventory) -> Vec<(String, String)> {
    inventory
        .store
        .profiles
        .values()
        .map(|p| {
            (
                p.id.clone(),
                if p.evidence_kind == GelEvidenceKind::SignedResidual {
                    format!("{} | derived algebra (verify before reporting)", p.id)
                } else {
                    format!(
                        "{} | {} | run {} / lane {} | {:?}",
                        p.id, p.label, p.run_id, p.lane_id, p.evidence_kind
                    )
                },
            )
        })
        .collect()
}

fn ordered_selection(
    ui: &mut Ui,
    salt: &str,
    selected: &mut Vec<String>,
    choices: Vec<String>,
    analyses: bool,
) {
    ui.push_id(salt, |ui| {
        if choices.is_empty() {
            ui.small("No saved records available.");
        }
        for id in choices {
            let mut checked = selected.contains(&id);
            if ui.checkbox(&mut checked, &id).changed() {
                if checked {
                    selected.push(id);
                } else {
                    selected.retain(|value| value != &id);
                }
            }
        }
        let mut movement = None;
        let mut removal = None;
        for (i, id) in selected.iter().enumerate() {
            ui.horizontal(|ui| {
                let role = if analyses {
                    match i {
                        0 => "Before",
                        1 => "After",
                        _ => "Additional",
                    }
                } else {
                    "Comparison"
                };
                ui.label(format!("{}. {role}: {id}", i + 1));
                if ui.add_enabled(i > 0, egui::Button::new("Up")).clicked() {
                    movement = Some((i, i - 1));
                }
                if ui
                    .add_enabled(i + 1 < selected.len(), egui::Button::new("Down"))
                    .clicked()
                {
                    movement = Some((i, i + 1));
                }
                if ui.button("Remove").clicked() {
                    removal = Some(i);
                }
            });
        }
        if let Some((a, b)) = movement {
            selected.swap(a, b);
        }
        if let Some(i) = removal {
            selected.remove(i);
        }
    });
}

fn profile_provenance(ui: &mut Ui, profile: &GelProfile) {
    if profile.evidence_kind == GelEvidenceKind::SignedResidual {
        ui.label(format!(
            "{} | derived algebra | size {}",
            profile.label,
            profile.size_kind.unit()
        ));
        ui.small("Legacy inherited run/lane fields are not the derived result's identity; inspect its recorded operands.");
    } else {
        ui.label(format!(
            "{} | run {} / lane {} | {:?} | size {}",
            profile.label,
            profile.run_id,
            profile.lane_id,
            profile.evidence_kind,
            profile.size_kind.unit()
        ));
    }
    ui.small(format!(
        "Signal units: {} | integration measure: {}",
        profile.signal_units, profile.measure
    ));
    ui.small(format!(
        "{}: {}",
        if profile.evidence_kind == GelEvidenceKind::SignedResidual {
            "Algebra input digest (requires replay verification)"
        } else {
            "Source SHA-256"
        },
        profile.source_sha256
    ));
    for step in &profile.processing {
        ui.small(step);
    }
    for warning in &profile.warnings {
        ui.small(warning);
    }
}

fn comparison_summary(ui: &mut Ui, report: &GelComparisonReport) {
    ui.label(format!(
        "Saved comparison: {} | {}",
        report.request.comparison_id, report.schema
    ));
    ui.small(format!(
        "Frozen inputs: {} / {}; requested range {} to {} {}",
        report.left.id,
        report.right.id,
        report.request.range[0],
        report.request.range[1],
        report.left.size_kind.unit()
    ));
    ui.small("The form above does not change this saved result. Summaries are conditional on shared valid coverage.");
    ui.small(format!(
        "Frozen evidence: {:?} | signal units: {} | measure: {}",
        report.left.evidence_kind, report.left.signal_units, report.left.measure
    ));
    egui::Grid::new("gel_quantitative_summary")
        .striped(true)
        .show(ui, |ui| {
            ui.label("Saved quantity");
            ui.label("Left");
            ui.label("Right");
            ui.end_row();
            for (name, a, b) in [
                (
                    "Total integrated signal",
                    report.left_summary.total_weight,
                    report.right_summary.total_weight,
                ),
                (
                    "Q10 size",
                    report.left_summary.q10,
                    report.right_summary.q10,
                ),
                (
                    "Median size",
                    report.left_summary.median,
                    report.right_summary.median,
                ),
                (
                    "Q90 size",
                    report.left_summary.q90,
                    report.right_summary.q90,
                ),
            ] {
                ui.label(name);
                ui.label(format!("{a:.4}"));
                ui.label(format!("{b:.4}"));
                ui.end_row();
            }
            if let Some(threshold) = report.request.short_threshold {
                ui.label(format!(
                    "Fraction below {threshold} {}",
                    report.left.size_kind.unit()
                ));
                ui.label(format_optional(
                    report.left_summary.fraction_below_threshold,
                ));
                ui.label(format_optional(
                    report.right_summary.fraction_below_threshold,
                ));
                ui.end_row();
            }
        });
    ui.label(format!(
        "Descriptive KS distance: {:.5} | Wasserstein distance: {:.5} {}",
        report.ks_distance,
        report.wasserstein_distance,
        report.left.size_kind.unit()
    ));
    ui.label("p-value: not computed. Correlated signal bins / pixels are not independent fragment samples or biological replicates.");
    qq_plot(ui, report);
    ui.collapsing("Frozen profile provenance and excluded intervals", |ui| {
        ui.label("Left");
        profile_provenance(ui, &report.left);
        ui.label("Right");
        profile_provenance(ui, &report.right);
        ui.label(format!(
            "{} excluded size intervals (missing in either profile, not replaced with zero)",
            report.excluded_intervals.len()
        ));
        egui::ScrollArea::vertical()
            .id_salt("gel_comparison_exclusions")
            .max_height(100.0)
            .show(ui, |ui| {
                for [a, b] in &report.excluded_intervals {
                    ui.small(format!("{a:.4} to {b:.4}"));
                }
            });
    });
    for warning in &report.warnings {
        ui.small(warning);
    }
}

fn format_optional(value: Option<f64>) -> String {
    value
        .map(|x| format!("{x:.4}"))
        .unwrap_or_else(|| "Not available".into())
}

/// A shared equal-scale size axis, x=left and y=right, with y increasing upward.
/// Out-of-range/nonfinite coordinates are omitted, never clipped into agreement.
pub(super) fn qq_position(
    rect: egui::Rect,
    range: [f64; 2],
    left: f64,
    right: f64,
) -> Option<Pos2> {
    if !valid_range(range)
        || !rect.is_finite()
        || !rect.is_positive()
        || !left.is_finite()
        || !right.is_finite()
        || left < range[0]
        || left > range[1]
        || right < range[0]
        || right > range[1]
    {
        return None;
    }
    let extent = range[1] - range[0];
    Some(egui::pos2(
        rect.left() + ((left - range[0]) / extent) as f32 * rect.width(),
        rect.bottom() - ((right - range[0]) / extent) as f32 * rect.height(),
    ))
}

fn qq_plot(ui: &mut Ui, report: &GelComparisonReport) {
    ui.label("QQ plot: equal probability quantiles of saved weighted profiles");
    let side = (ui.available_width() - 100.0).clamp(80.0, 360.0);
    let (outer, _) =
        ui.allocate_exact_size(egui::vec2(side + 100.0, side + 76.0), egui::Sense::hover());
    let rect =
        egui::Rect::from_min_size(outer.min + egui::vec2(76.0, 24.0), egui::vec2(side, side));
    let painter = ui.painter_at(outer);
    let ink = ui.visuals().text_color();
    let stroke = egui::Stroke::new(1.0, ink);
    painter.line_segment([rect.left_top(), rect.left_bottom()], stroke);
    painter.line_segment([rect.left_bottom(), rect.right_bottom()], stroke);
    painter.line_segment(
        [rect.left_bottom(), rect.right_top()],
        egui::Stroke::new(1.0, Color32::GRAY),
    );
    let range = report.request.range;
    for i in 0..=4 {
        let f = i as f32 / 4.0;
        let value = range[0] + f64::from(f) * (range[1] - range[0]);
        let x = rect.left() + f * side;
        let y = rect.bottom() - f * side;
        painter.line_segment(
            [
                egui::pos2(x, rect.bottom()),
                egui::pos2(x, rect.bottom() + 4.0),
            ],
            stroke,
        );
        painter.line_segment(
            [egui::pos2(rect.left() - 4.0, y), egui::pos2(rect.left(), y)],
            stroke,
        );
        painter.text(
            egui::pos2(x, rect.bottom() + 6.0),
            egui::Align2::CENTER_TOP,
            format!("{value:.1}"),
            egui::FontId::proportional(11.0),
            ink,
        );
        painter.text(
            egui::pos2(rect.left() - 6.0, y),
            egui::Align2::RIGHT_CENTER,
            format!("{value:.1}"),
            egui::FontId::proportional(11.0),
            ink,
        );
    }
    for point in &report.qq {
        if let Some(position) = qq_position(rect, range, point.left, point.right) {
            painter.circle_filled(position, 2.5, SAMPLE_COLOR);
        }
    }
    painter.text(
        rect.left_top() - egui::vec2(0.0, 8.0),
        egui::Align2::LEFT_BOTTOM,
        format!("Right size ({})", report.right.size_kind.unit()),
        egui::FontId::proportional(12.0),
        ink,
    );
    painter.text(
        egui::pos2(rect.center().x, rect.bottom() + 34.0),
        egui::Align2::CENTER_TOP,
        format!(
            "Left size ({}) | gray line: identity",
            report.left.size_kind.unit()
        ),
        egui::FontId::proportional(12.0),
        ink,
    );
}
