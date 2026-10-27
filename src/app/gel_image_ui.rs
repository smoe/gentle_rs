//! Measured-gel workspace. Screen geometry is display-only; the shared
//! engine owns persistent drafts, calibration, immutable reports and exports.

use super::*;
use base64::{Engine as _, engine::general_purpose::STANDARD};
use egui::Color32;
use gentle_protocol::gel_image::*;
use gentle_protocol::gel_workbench::*;
use std::sync::{Weak, mpsc};

#[path = "gel_image_ui/workbench.rs"]
mod workbench;
use workbench::{DetectionIdentity, PendingSuggestions, WorkbenchUi};

const REFERENCE_COLOR: Color32 = Color32::from_rgb(210, 125, 20);
const SAMPLE_COLOR: Color32 = Color32::from_rgb(0, 130, 145);

#[derive(Clone, Copy, Default, PartialEq, Eq)]
enum MarkTool {
    #[default]
    Lane,
    Reference,
    Sample,
    Move,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum BandSelection {
    Reference(usize),
    Sample(usize),
}

/// Pixel centers, not thumbnail pixels or image-edge coordinates.
#[derive(Clone, Copy)]
struct ImageTransform {
    rect: egui::Rect,
    width: u32,
    height: u32,
}

impl ImageTransform {
    fn to_screen(self, point: GelImagePoint) -> Pos2 {
        egui::pos2(
            self.rect.left() + ((point.x + 0.5) / self.width as f64) as f32 * self.rect.width(),
            self.rect.top() + ((point.y + 0.5) / self.height as f64) as f32 * self.rect.height(),
        )
    }

    fn to_original(self, point: Pos2) -> Option<GelImagePoint> {
        if !self.rect.contains(point) || self.rect.width() <= 0.0 || self.rect.height() <= 0.0 {
            return None;
        }
        Some(GelImagePoint {
            x: (((point.x - self.rect.left()) / self.rect.width()) as f64 * self.width as f64
                - 0.5)
                .clamp(0.0, (self.width - 1) as f64),
            y: (((point.y - self.rect.top()) / self.rect.height()) as f64 * self.height as f64
                - 0.5)
                .clamp(0.0, (self.height - 1) as f64),
        })
    }
}

fn contains(lane: &GelImageLane, point: GelImagePoint) -> bool {
    point.x >= lane.min.x && point.x <= lane.max.x && point.y >= lane.min.y && point.y <= lane.max.y
}

fn fresh_id(prefix: &str, occupied: impl Iterator<Item = String>) -> String {
    let occupied: BTreeSet<_> = occupied.collect();
    (1..)
        .map(|n| format!("{prefix}_{n}"))
        .find(|id| !occupied.contains(id))
        .expect("unbounded identifier space")
}

fn empty_draft(image: &GelImageDescriptor) -> GelImageAnalysisRequest {
    GelImageAnalysisRequest {
        report_id: String::new(),
        image_id: image.image_id.clone(),
        image_sha256: image.sha256.clone(),
        size_kind: GelSizeKind::LinearDnaBp,
        migration: GelMigrationDirection::Down,
        lanes: vec![],
        ladder: GelImageLadder {
            lane_id: String::new(),
            label: "Custom ladder".into(),
            source: String::new(),
            gel_system: None,
            prestained: false,
            bands: vec![],
        },
        sample_bands: vec![],
    }
}

struct GelJobResult {
    detached: Option<crate::engine::DetachedEngineExecution>,
    preview: Option<(Arc<GelImageRecord>, egui::ColorImage)>,
    report: Option<Arc<GelImageAnalysisReport>>,
    source: Option<Arc<GelSourceSummary>>,
    suggestions: Option<PendingSuggestions>,
    workbench: Option<Box<GelWorkbenchResult>>,
    algebra_inspection: Option<workbench::AlgebraInspection>,
    comparison_inspection: Option<workbench::ComparisonInspection>,
    message: String,
}

enum GelJob {
    Import {
        path: String,
        image_id: String,
        tiff: bool,
    },
    Preview(Arc<GelImageRecord>),
    Analyze(GelImageAnalysisRequest),
    Export(GelImageExportRequest),
    Workbench(GelWorkbenchRequest),
    ComparisonPlot {
        comparison_id: String,
        mode: GelComparisonDisplay,
    },
    Detect {
        request: GelDetectRequest,
        draft: Box<GelImageAnalysisRequest>,
    },
}

impl GelJob {
    fn is_external_export(&self) -> bool {
        matches!(
            self,
            Self::Export(_) | Self::Workbench(GelWorkbenchRequest::Export(_))
        )
    }
}

/// Decode only in workers, never during an idle frame. Saved previews are not
/// trusted evidence: rebuilding from the original also verifies its digest.
fn prepare_preview(
    record: Arc<GelImageRecord>,
) -> Result<(Arc<GelImageRecord>, egui::ColorImage), String> {
    let png = crate::gel_image::gel_image_export_preview(&record).map_err(|e| e.message)?;
    let bytes = STANDARD.decode(png).map_err(|e| e.to_string())?;
    let rgba = image::load_from_memory_with_format(&bytes, image::ImageFormat::Png)
        .map_err(|e| e.to_string())?
        .into_rgba8();
    let pixels = egui::ColorImage::from_rgba_unmultiplied(
        [rgba.width() as usize, rgba.height() as usize],
        rgba.as_raw(),
    );
    Ok((record, pixels))
}

fn execute_job(shared: &Arc<RwLock<GentleEngine>>, job: GelJob) -> Result<GelJobResult, String> {
    if let GelJob::Preview(record) = job {
        return Ok(GelJobResult {
            detached: None,
            preview: Some(prepare_preview(record)?),
            report: None,
            source: None,
            suggestions: None,
            workbench: None,
            algebra_inspection: None,
            comparison_inspection: None,
            message: "Image ready. Marks use original pixel coordinates.".into(),
        });
    }
    let mut detached = shared
        .read()
        .map_err(|_| "Engine lock unavailable")?
        .fork_detached_execution();
    let revision = detached.engine().mutation_revision();
    if let GelJob::ComparisonPlot {
        comparison_id,
        mode,
    } = job
    {
        let report = detached
            .engine()
            .state()
            .gel_images
            .workbench
            .comparisons
            .get(&comparison_id)
            .cloned()
            .ok_or("Saved comparison missing")?;
        let inspection = workbench::ComparisonInspection::new(revision, report, mode)?;
        return Ok(GelJobResult {
            detached: None,
            preview: None,
            report: None,
            source: None,
            suggestions: None,
            workbench: None,
            algebra_inspection: None,
            comparison_inspection: Some(inspection),
            message: "Saved comparison verified and displayed; no project change.".into(),
        });
    }
    let mut detection = None;
    let (op, import_id, mutates) = match job {
        GelJob::Import {
            path,
            image_id,
            tiff,
        } => {
            let path = if path.trim().is_empty() {
                rfd::FileDialog::new()
                    .set_title("Import a measured gel image")
                    .add_filter("Gel image", &["png", "jpg", "jpeg", "tif", "tiff"])
                    .pick_file()
                    .ok_or("Image selection canceled")?
                    .display()
                    .to_string()
            } else {
                path
            };
            (
                Operation::ImportGelImage {
                    request: GelImageImportRequest {
                        image_id: image_id.clone(),
                        path,
                        tiff_page: tiff.then_some(0),
                    },
                },
                Some(image_id),
                true,
            )
        }
        GelJob::Analyze(request) => (
            Operation::AnalyzeGelImage {
                request: Box::new(request),
            },
            None,
            true,
        ),
        GelJob::Export(mut request) => {
            if request.path.trim().is_empty() {
                let extension = match request.format {
                    GelImageExportFormat::Svg => "svg",
                    GelImageExportFormat::Tsv => "tsv",
                    GelImageExportFormat::Json => "json",
                };
                request.path = rfd::FileDialog::new()
                    .set_title("Export measured gel (choose a new file)")
                    .set_file_name(format!("{}.{}", request.report_id, extension))
                    .add_filter("Report", &[extension])
                    .save_file()
                    .ok_or("Export selection canceled")?
                    .display()
                    .to_string();
            }
            (Operation::ExportGelImageAnalysis { request }, None, false)
        }
        GelJob::Workbench(mut request) => {
            workbench::prepare_file_request(&mut request)?;
            let import_id = match &request {
                GelWorkbenchRequest::ImportSource(source) => Some(source.image_id.clone()),
                _ => None,
            };
            let mutates = request.mutates();
            (
                Operation::GelWorkbench {
                    request: Box::new(request),
                },
                import_id,
                mutates,
            )
        }
        GelJob::Detect { request, draft } => {
            detection = Some(DetectionIdentity {
                draft: *draft,
                request: request.clone(),
                revision,
            });
            (
                Operation::GelWorkbench {
                    request: Box::new(GelWorkbenchRequest::DetectBands(request)),
                },
                None,
                false,
            )
        }
        GelJob::Preview(_) | GelJob::ComparisonPlot { .. } => unreachable!(),
    };
    let reference_profile_id = match &op {
        Operation::GelWorkbench { request } => match request.as_ref() {
            GelWorkbenchRequest::InspectAlgebra {
                reference_profile_id,
                ..
            } => reference_profile_id.clone(),
            _ => None,
        },
        _ => None,
    };
    let mut result = detached.engine_mut().apply(op).map_err(|e| e.message)?;
    let preview = import_id
        .map(|id| prepare_preview(detached.engine().state().gel_images.images[&id].clone()))
        .transpose()?;
    let source = result
        .gel_workbench
        .as_ref()
        .and_then(|r| r.source.clone())
        .map(Arc::new);
    let suggestions = detection
        .zip(
            result
                .gel_workbench
                .as_mut()
                .and_then(|r| r.suggestions.take()),
        )
        .map(|(identity, suggestions)| PendingSuggestions {
            identity,
            suggestions,
        });
    let message = result
        .messages
        .iter()
        .chain(&result.warnings)
        .cloned()
        .collect::<Vec<_>>()
        .join("\n");
    let algebra_inspection = result
        .gel_workbench
        .as_mut()
        .and_then(|r| r.algebra.take())
        .map(|presentation| {
            workbench::AlgebraInspection::new(revision, presentation, reference_profile_id)
        })
        .transpose()?;
    Ok(GelJobResult {
        detached: mutates.then_some(detached),
        preview,
        report: result.gel_image_analysis.map(|r| Arc::from(*r)),
        source,
        suggestions,
        workbench: result.gel_workbench,
        algebra_inspection,
        comparison_inspection: None,
        message,
    })
}

pub(super) struct GelImageEditor {
    pub(super) open: bool,
    bound_engine: Option<Weak<RwLock<GentleEngine>>>,
    bound_instance_id: Option<u64>,
    job: Option<mpsc::Receiver<Result<GelJobResult, String>>>,
    export_running: bool,
    comparison_plot_running: bool,
    canceled: bool,
    image_id: String,
    image: Option<Arc<GelImageRecord>>,
    texture: Option<egui::TextureHandle>,
    draft: Option<GelImageAnalysisRequest>,
    stored_draft: Option<Arc<GelImageAnalysisRequest>>,
    report: Option<Arc<GelImageAnalysisReport>>,
    import_path: String,
    tiff: bool,
    export_path: String,
    export_format: GelImageExportFormat,
    tool: MarkTool,
    selected: Option<BandSelection>,
    selected_lane: String,
    reference_size: f64,
    zoom: f32,
    drag_start: Option<Pos2>,
    status: String,
    help: bool,
    workbench: WorkbenchUi,
}

impl Default for GelImageEditor {
    fn default() -> Self {
        Self {
            open: false,
            bound_engine: None,
            bound_instance_id: None,
            job: None,
            export_running: false,
            comparison_plot_running: false,
            canceled: false,
            image_id: String::new(),
            image: None,
            texture: None,
            draft: None,
            stored_draft: None,
            report: None,
            import_path: String::new(),
            tiff: false,
            export_path: String::new(),
            export_format: GelImageExportFormat::Svg,
            tool: MarkTool::Lane,
            selected: None,
            selected_lane: String::new(),
            reference_size: 0.0,
            zoom: 1.0,
            drag_start: None,
            status: "Import a gel, draw the ladder lane, then mark its known bands.".into(),
            help: true,
            workbench: WorkbenchUi::default(),
        }
    }
}

impl GelImageEditor {
    fn start(&mut self, shared: &Arc<RwLock<GentleEngine>>, ctx: &egui::Context, job: GelJob) {
        if self.job.is_some() {
            return;
        }
        self.export_running = job.is_external_export();
        self.comparison_plot_running = matches!(job, GelJob::ComparisonPlot { .. });
        self.canceled = false;
        self.status = match &job {
            GelJob::Import { .. } => "Loading and preserving original image...",
            GelJob::Preview(_) => "Verifying original and preparing preview...",
            GelJob::Analyze(_) => "Validating assignments and calibrating band sizes...",
            GelJob::Export(_) => "Validating and writing export (no overwrite)...",
            GelJob::Workbench(GelWorkbenchRequest::Export(_)) => {
                "Writing saved workbench report (external file; no overwrite)..."
            }
            GelJob::Workbench(_) => "Running measured-gel workbench request...",
            GelJob::ComparisonPlot { .. } => {
                "Verifying saved comparison and preparing selected display..."
            }
            GelJob::Detect { .. } => {
                "Finding unconfirmed peaks in selected sample lane rectangles..."
            }
        }
        .into();
        let (tx, rx) = mpsc::channel();
        self.job = Some(rx);
        let shared = shared.clone();
        let ctx = ctx.clone();
        let worker = std::thread::Builder::new()
            .name("gel-image-workbench".into())
            .stack_size(16 * 1024 * 1024)
            .spawn(move || {
                let result = execute_job(&shared, job);
                let _ = tx.send(result);
                ctx.request_repaint();
                ctx.request_repaint_of(ViewportId::ROOT);
            });
        if let Err(error) = worker {
            self.job = None;
            self.export_running = false;
            self.comparison_plot_running = false;
            self.status = format!("Could not start gel worker: {error}");
        }
    }

    fn bind_project(&mut self, shared: &Arc<RwLock<GentleEngine>>) {
        let instance_id = shared.read().unwrap().instance_id();
        if self
            .bound_engine
            .as_ref()
            .is_some_and(|old| !old.ptr_eq(&Arc::downgrade(shared)))
            || self.bound_instance_id.is_some_and(|old| old != instance_id)
        {
            // Dropping the receiver prevents a late result committing into a new project.
            if let Some(pending) = self.job.take() {
                std::thread::spawn(move || drop(pending));
            }
            *self = Self::default();
        }
        self.bound_engine = Some(Arc::downgrade(shared));
        self.bound_instance_id = Some(instance_id);
    }

    fn poll(&mut self, shared: &Arc<RwLock<GentleEngine>>, ctx: &egui::Context) {
        self.bind_project(shared);
        self.sync_draft(shared);
        let revision = shared.read().unwrap().mutation_revision();
        self.workbench.invalidate_algebra(revision);
        self.workbench.invalidate_comparison(revision);
        if !self.open {
            self.workbench.algebra_inspection = None;
            self.workbench.comparison_inspection = None;
            // Closing ends this inspection, even if reopened before it finishes.
            if self.comparison_plot_running {
                self.canceled = true;
            }
        }
        if self
            .workbench
            .pending
            .as_ref()
            .is_some_and(|pending| !self.suggestions_match(pending, revision))
        {
            self.workbench.clear_suggestions();
            self.status =
                "Project, draft or detector settings changed; old proposals discarded.".into();
        }
        let outcome = match self.job.as_ref().map(|rx| rx.try_recv()) {
            Some(Ok(result)) => result,
            Some(Err(mpsc::TryRecvError::Disconnected)) => {
                Err("Gel worker stopped without a result".into())
            }
            _ => return,
        };
        self.job = None;
        self.export_running = false;
        self.comparison_plot_running = false;
        if self.canceled {
            self.canceled = false;
            self.status = "Pending result discarded. No project change was applied.".into();
            std::thread::spawn(move || drop(outcome));
            return;
        }
        match outcome {
            Err(error) => self.status = error,
            Ok(mut result) => {
                if let Some(mut inspection) = result.comparison_inspection.take() {
                    if self.open
                        && self.workbench.comparison_matches(
                            &inspection,
                            shared.read().unwrap().mutation_revision(),
                        )
                    {
                        inspection.install_texture(ctx);
                        self.workbench.comparison_inspection = Some(inspection);
                        self.status = result.message;
                    } else {
                        self.status = "Project, comparison or display mode changed; stale inspection discarded.".into();
                    }
                    return;
                }
                if let Some(mut inspection) = result.algebra_inspection.take() {
                    let engine = shared.read().unwrap();
                    if self.open
                        && self
                            .workbench
                            .algebra_matches(&inspection, engine.mutation_revision())
                    {
                        inspection.install_texture(ctx);
                        self.workbench.algebra_inspection = Some(inspection);
                        self.status = result.message;
                    } else {
                        self.status = "Project or selected algebra result changed; stale inspection discarded.".into();
                    }
                    return;
                }
                if let Some(mut detached) = result.detached.take() {
                    let commit = shared
                        .write()
                        .unwrap()
                        .commit_detached_execution(&mut detached);
                    match commit {
                        Ok(old) => {
                            self.workbench.algebra_inspection = None;
                            self.workbench.comparison_inspection = None;
                            std::thread::spawn(move || drop(old));
                        }
                        Err(error) => {
                            self.status =
                                format!("Project changed; result not applied. {}", error.message);
                            std::thread::spawn(move || drop(detached));
                            return;
                        }
                    }
                }
                if result.preview.as_ref().is_some_and(|(record, _)| {
                    !shared
                        .read()
                        .unwrap()
                        .state()
                        .gel_images
                        .images
                        .get(&record.descriptor.image_id)
                        .is_some_and(|current| Arc::ptr_eq(current, record))
                }) {
                    self.status = "Image changed or was removed; stale preview discarded.".into();
                    std::thread::spawn(move || drop(result));
                    return;
                }
                if let Some((record, pixels)) = result.preview {
                    self.image_id = record.descriptor.image_id.clone();
                    self.texture = Some(ctx.load_texture(
                        format!("gel:{}", record.descriptor.sha256),
                        pixels,
                        egui::TextureOptions::LINEAR,
                    ));
                    self.image = Some(record.clone());
                    self.draft = Some(empty_draft(&record.descriptor));
                    self.stored_draft = None;
                    self.selected = None;
                    self.selected_lane.clear();
                    self.report = None;
                    self.zoom = 1.0;
                    self.workbench.clear_suggestions();
                    self.workbench.detection_lanes.clear();
                }
                if let Some(report) = result.report {
                    self.report = Some(report);
                }
                self.status = result.message;
                if let Some(source) = result.source {
                    self.status
                        .push_str(&format!("\nSource extraction: {}", source.extraction));
                }
                if let Some(output) = result.workbench {
                    if let Some(profile) = output.profile {
                        if profile.evidence_kind == GelEvidenceKind::SignedResidual {
                            self.workbench.selected_algebra = profile.id.clone();
                            self.status.push_str(&format!("\nSaved derived algebra profile: {}. Inspect to verify; inherited run/lane fields are not its identity.", profile.id));
                        } else {
                            self.status.push_str(&format!(
                                "\nSaved profile: {} (run {}, lane {}, {:?})",
                                profile.id, profile.run_id, profile.lane_id, profile.evidence_kind
                            ));
                        }
                    }
                    if let Some(comparison) = output.comparison {
                        self.workbench.saved_comparison = comparison.request.comparison_id;
                    }
                    if let Some(hash) = output.exported_sha256 {
                        self.status =
                            format!("External report written (no project commit). SHA-256: {hash}");
                    }
                }
                if let Some(pending) = result.suggestions {
                    let revision = shared.read().unwrap().mutation_revision();
                    if self.suggestions_match(&pending, revision) {
                        self.workbench.selected_suggestions.clear();
                        self.workbench.pending = Some(pending);
                        self.status = "Unconfirmed proposals ready. Select and explicitly accept sample marks; no ladder assignments were made.".into();
                    } else {
                        self.workbench.clear_suggestions();
                        self.status = "Draft, image, project or detection request changed; stale proposals discarded.".into();
                    }
                }
            }
        }
    }

    fn save_draft(&mut self, shared: &Arc<RwLock<GentleEngine>>) {
        let Some(request) = &self.draft else {
            return;
        };
        let mut engine = shared.write().unwrap();
        match engine.apply(Operation::SaveGelImageDraft {
            request: Arc::new(request.clone()),
        }) {
            Ok(_) => {
                self.stored_draft = engine
                    .state()
                    .gel_images
                    .drafts
                    .get(&self.image_id)
                    .cloned();
                self.status = "Assignments updated in project. Save Project to retain them on disk; Analyze to validate sizing.".into();
            }
            Err(error) => self.status = format!("Draft not saved: {}", error.message),
        }
    }

    fn sync_draft(&mut self, shared: &Arc<RwLock<GentleEngine>>) {
        let engine = shared.read().unwrap();
        let store = &engine.state().gel_images;
        if let Some(image) = &self.image {
            if !store
                .images
                .get(&self.image_id)
                .is_some_and(|current| Arc::ptr_eq(current, image))
            {
                self.image = None;
                self.texture = None;
                self.draft = None;
                self.stored_draft = None;
                self.report = None;
                self.workbench.clear_suggestions();
                self.workbench.detection_lanes.clear();
                self.status =
                    "The selected image changed or was removed. Select an image to reload it."
                        .into();
                return;
            }
            let current = store.drafts.get(&self.image_id);
            let unchanged = match (current, &self.stored_draft) {
                (Some(a), Some(b)) => Arc::ptr_eq(a, b),
                (None, None) => true,
                _ => false,
            };
            if !unchanged {
                self.draft = Some(
                    current
                        .map(|r| (**r).clone())
                        .unwrap_or_else(|| empty_draft(&image.descriptor)),
                );
                self.stored_draft = current.cloned();
                self.selected = None;
                self.drag_start = None;
                self.workbench.clear_suggestions();
            }
        }
        if self
            .report
            .as_ref()
            .is_some_and(|r| !store.analyses.contains_key(&r.request.report_id))
        {
            self.report = None;
        }
    }

    fn contents(&mut self, ui: &mut Ui, shared: &Arc<RwLock<GentleEngine>>) {
        help_topics::button(ui, HelpTopic::GelImage);
        self.sync_draft(shared);
        #[cfg(feature = "gui-test-support")]
        crate::gui_test_support::register_rect(
            ui.ctx().clone(),
            "window.gel_image",
            "window.gel_image",
            None,
            crate::gui_test_support::GuiTestWidgetKind::Window,
            ui.max_rect(),
            true,
            true,
            true,
            Some(if self.job.is_some() {
                "running"
            } else {
                "ready"
            }),
        );
        ui.horizontal(|ui| {
            if ui.button("Help").clicked() {
                self.help = !self.help;
            }
            ui.label("Measured Gel Workbench");
        });
        if self.help {
            ui.group(|ui| {
                ui.label("1. Import PNG/JPEG or explicitly select TIFF page 0. Originals are never altered.");
                ui.label("2. Draw lane rectangles. The first lane is the ladder. Select a different ladder only before marking references.");
                ui.label("3. Enter the ladder source and each known size, then click its band center. Mark sample bands in sample lanes.");
                ui.label("4. Analyze to validate and save a sizing report. Review out-of-range results and warnings, then export.");
                ui.label("5. Select sample rectangles for peak proposals; review and accept explicitly. Extract a profile from a saved calibration, compare saved profiles across runs, then export selected saved evidence.");
                ui.small("Save Project retains images, drafts, profiles and comparisons. Zoom/scroll do not change coordinates. DNA: linear fragments only; protein: apparent SDS mass. Display intensity is not concentration, molecule count or yield. No smile correction or inferred ladder identities.");
            });
        }
        ui.horizontal(|ui| {
            ui.label(&self.status);
            if self.job.is_some() {
                ui.spinner();
                if self.export_running {
                    ui.small("External export continues even if this window closes or the project changes; it cannot be canceled or undone.");
                }
                if !self.export_running && !self.canceled && ui.button("Cancel pending result").clicked() {
                    self.canceled = true;
                    self.status = "Discarding result; waiting for background work to finish. No project change will be applied.".into();
                }
            }
        });
        ui.add_enabled_ui(self.job.is_none(), |ui| {
            ui.horizontal(|ui| {
                let (undo,redo) = { let engine = shared.read().unwrap(); (engine.undo_available(),engine.redo_available()) };
                let action = if ui.add_enabled(undo > 0,egui::Button::new("Undo project edit")).clicked() { Some(false) }
                    else if ui.add_enabled(redo > 0,egui::Button::new("Redo project edit")).clicked() { Some(true) } else { None };
                if let Some(redo) = action {
                    let result = { let mut engine = shared.write().unwrap(); if redo { engine.redo_last_operation() } else { engine.undo_last_operation() } };
                    self.status = match result { Ok(_) => "Project history updated.".into(), Err(e) => e.message };
                    self.sync_draft(shared);
                }
            });
            self.source_controls(ui, shared);
            self.native_source_controls(ui, shared);
            if self.texture.is_some() {
            let before = self.draft.clone();
            ui.separator();
            self.assignment_controls(ui);
            ui.horizontal(|ui| {
                for (tool, label) in [(MarkTool::Lane, "Draw lane"), (MarkTool::Reference, "Ladder band"), (MarkTool::Sample, "Sample band"), (MarkTool::Move, "Select / move band")] {
                    if semantic_button(ui, label, match tool { MarkTool::Lane => "gel.tool.lane", MarkTool::Reference => "gel.tool.reference", MarkTool::Sample => "gel.tool.sample", MarkTool::Move => "gel.tool.move" }, self.tool == tool).clicked() { self.tool = tool; self.drag_start = None; }
                }
                ui.add(egui::Slider::new(&mut self.zoom, 0.5..=8.0).text("Zoom"));
                if ui.button("Fit").clicked() { self.zoom = 1.0; }
            });
            ui.small("Orange: confirmed ladder sizes. Teal: sample marks. Drag to draw a lane or move a selected band; scroll to pan the zoomed image.");
            self.canvas(ui);
            self.band_editor(ui);
            if before != self.draft { self.save_draft(shared); }
            self.analysis_controls(ui, shared);
            }
            self.workbench_controls(ui, shared);
        });
    }

    fn source_controls(&mut self, ui: &mut Ui, shared: &Arc<RwLock<GentleEngine>>) {
        ui.horizontal(|ui| {
            ui.label("Image path (optional)"); ui.text_edit_singleline(&mut self.import_path);
            ui.checkbox(&mut self.tiff, "TIFF page 0 only").on_hover_text("Explicitly confirm first-page-only import for TIFF; leave unchecked for PNG/JPEG.");
            if semantic_button(ui, "Import image...", "gel.import", false).clicked() {
                let id = fresh_id("gel", shared.read().unwrap().state().gel_images.images.keys().cloned());
                self.start(shared, ui.ctx(), GelJob::Import { path: self.import_path.clone(), image_id: id, tiff: self.tiff });
            }
        });
        let (images, reports) = {
            let engine = shared.read().unwrap();
            (
                engine
                    .state()
                    .gel_images
                    .images
                    .values()
                    .cloned()
                    .collect::<Vec<_>>(),
                engine
                    .state()
                    .gel_images
                    .analyses
                    .values()
                    .filter(|r| r.image.image_id == self.image_id)
                    .cloned()
                    .collect::<Vec<_>>(),
            )
        };
        ui.horizontal(|ui| {
            egui::ComboBox::from_id_salt("gel_images")
                .selected_text(if self.image_id.is_empty() {
                    "Select imported image"
                } else {
                    &self.image_id
                })
                .show_ui(ui, |ui| {
                    for record in images {
                        if ui
                            .selectable_label(
                                self.image_id == record.descriptor.image_id,
                                format!(
                                    "{}: {}",
                                    record.descriptor.image_id, record.descriptor.source_name
                                ),
                            )
                            .clicked()
                        {
                            self.start(shared, ui.ctx(), GelJob::Preview(record));
                        }
                    }
                });
            egui::ComboBox::from_id_salt("gel_reports")
                .selected_text("Reopen a saved report")
                .show_ui(ui, |ui| {
                    for report in reports {
                        if ui.button(&report.request.report_id).clicked() {
                            self.draft = Some(report.request.clone());
                            self.report = Some(report);
                            self.selected = None;
                            self.save_draft(shared);
                            ui.close();
                        }
                    }
                });
            if let Some(image) = &self.image {
                ui.small(format!(
                    "{} x {} original pixels | {}",
                    image.descriptor.width, image.descriptor.height, image.descriptor.format
                ));
            }
        });
    }

    fn assignment_controls(&mut self, ui: &mut Ui) {
        let Some(draft) = &mut self.draft else {
            return;
        };
        ui.horizontal(|ui| {
            egui::ComboBox::from_id_salt("gel_units")
                .selected_text(draft.size_kind.unit())
                .show_ui(ui, |ui| {
                    ui.selectable_value(
                        &mut draft.size_kind,
                        GelSizeKind::LinearDnaBp,
                        "Linear DNA (bp)",
                    );
                    ui.selectable_value(
                        &mut draft.size_kind,
                        GelSizeKind::SdsProteinKda,
                        "SDS protein (kDa)",
                    );
                });
            egui::ComboBox::from_id_salt("gel_migration")
                .selected_text(format!("Migration: {:?}", draft.migration))
                .show_ui(ui, |ui| {
                    for direction in [
                        GelMigrationDirection::Down,
                        GelMigrationDirection::Up,
                        GelMigrationDirection::Right,
                        GelMigrationDirection::Left,
                    ] {
                        ui.selectable_value(
                            &mut draft.migration,
                            direction,
                            format!("{direction:?}"),
                        );
                    }
                });
            ui.label("Ladder name");
            ui.text_edit_singleline(&mut draft.ladder.label);
        });
        ui.horizontal(|ui| {
            ui.label("Ladder source / catalog / lot");
            ui.text_edit_singleline(&mut draft.ladder.source);
            ui.checkbox(&mut draft.ladder.prestained, "Prestained protein ladder");
        });
        if draft.size_kind == GelSizeKind::SdsProteinKda {
            ui.horizontal(|ui| {
                ui.label("Gel system (required if prestained)");
                let text = draft.ladder.gel_system.get_or_insert_with(String::new);
                ui.text_edit_singleline(text);
            });
        }
        ui.horizontal(|ui| {
            egui::ComboBox::from_id_salt("gel_lane")
                .selected_text(if self.selected_lane.is_empty() {
                    "Select lane"
                } else {
                    &self.selected_lane
                })
                .show_ui(ui, |ui| {
                    for lane in &draft.lanes {
                        ui.selectable_value(&mut self.selected_lane, lane.id.clone(), &lane.label);
                    }
                });
            if let Some(lane) = draft.lanes.iter_mut().find(|l| l.id == self.selected_lane) {
                ui.label("Lane label");
                ui.text_edit_singleline(&mut lane.label);
                if ui
                    .add_enabled(
                        draft.ladder.bands.is_empty()
                            && !draft.sample_bands.iter().any(|b| b.lane_id == lane.id),
                        egui::Button::new("Use as ladder"),
                    )
                    .clicked()
                {
                    draft.ladder.lane_id = lane.id.clone();
                }
                if ui.button("Delete lane and its marks").clicked() {
                    draft
                        .sample_bands
                        .retain(|b| b.lane_id != self.selected_lane);
                    if draft.ladder.lane_id == self.selected_lane {
                        draft.ladder.bands.clear();
                        draft.ladder.lane_id.clear();
                    }
                    draft.lanes.retain(|l| l.id != self.selected_lane);
                    self.selected_lane.clear();
                    self.selected = None;
                }
            }
            ui.label(format!("Ladder: {}", draft.ladder.lane_id));
        });
        if self.tool == MarkTool::Reference {
            ui.horizontal(|ui| {
                ui.label(format!("Next known size ({})", draft.size_kind.unit()));
                ui.add(egui::DragValue::new(&mut self.reference_size).range(0.0..=1e9).speed(1.0));
                ui.small("Enter a positive value from your ladder documentation before clicking. Reset after each mark.");
            });
        }
    }

    fn canvas(&mut self, ui: &mut Ui) -> Option<ImageTransform> {
        let (Some(image), Some(texture), Some(draft)) =
            (&self.image, &self.texture, self.draft.clone())
        else {
            return None;
        };
        let width = image.descriptor.width;
        let height = image.descriptor.height;
        let texture_id = texture.id();
        let fit = (ui.available_width() / width as f32).min(440.0 / height as f32);
        let size = egui::vec2(width as f32, height as f32) * fit * self.zoom;
        let output = egui::ScrollArea::both()
            .id_salt("gel_canvas_scroll")
            .max_height(480.0)
            .auto_shrink([false, false])
            .show(ui, |ui| {
                let (rect, response) = ui.allocate_exact_size(size, egui::Sense::click_and_drag());
                #[cfg(feature = "gui-test-support")]
                crate::gui_test_support::register_response(
                    &response,
                    "gel.canvas",
                    "window.gel_image",
                    None,
                    crate::gui_test_support::GuiTestWidgetKind::Row,
                    false,
                );
                let transform = ImageTransform {
                    rect,
                    width,
                    height,
                };
                let painter = ui.painter_at(rect);
                painter.image(
                    texture_id,
                    rect,
                    egui::Rect::from_min_max(egui::pos2(0.0, 0.0), egui::pos2(1.0, 1.0)),
                    Color32::WHITE,
                );
                for lane in &draft.lanes {
                    let color = if lane.id == draft.ladder.lane_id {
                        REFERENCE_COLOR
                    } else {
                        SAMPLE_COLOR
                    };
                    let lane_rect = egui::Rect::from_two_pos(
                        transform.to_screen(lane.min),
                        transform.to_screen(lane.max),
                    );
                    painter.rect_stroke(
                        lane_rect,
                        0.0,
                        egui::Stroke::new(1.5, color),
                        egui::StrokeKind::Inside,
                    );
                    painter.text(
                        lane_rect.left_top(),
                        egui::Align2::LEFT_TOP,
                        &lane.label,
                        egui::FontId::proportional(12.0),
                        color,
                    );
                }
                let mut centers = vec![];
                for (i, band) in draft.ladder.bands.iter().enumerate() {
                    centers.push((
                        BandSelection::Reference(i),
                        band.center,
                        format!("{} {}", band.size, draft.size_kind.unit()),
                        REFERENCE_COLOR,
                    ));
                }
                for (i, band) in draft.sample_bands.iter().enumerate() {
                    centers.push((
                        BandSelection::Sample(i),
                        band.center,
                        band.label.clone(),
                        SAMPLE_COLOR,
                    ));
                }
                for (selected, center, label, color) in &centers {
                    let p = transform.to_screen(*center);
                    painter.circle_stroke(
                        p,
                        if self.selected == Some(*selected) {
                            7.0
                        } else {
                            4.0
                        },
                        egui::Stroke::new(2.0, *color),
                    );
                    painter.text(
                        p + egui::vec2(8.0, 0.0),
                        egui::Align2::LEFT_CENTER,
                        label,
                        egui::FontId::proportional(12.0),
                        *color,
                    );
                }
                if let Some(pending) = &self.workbench.pending
                    && self.suggestions_match(pending, pending.identity.revision)
                {
                    for band in &pending.suggestions.bands {
                        let p = transform.to_screen(band.center);
                        let selected = self.workbench.selected_suggestions.contains(&band.id);
                        painter.circle_stroke(
                            p,
                            if selected { 9.0 } else { 6.0 },
                            egui::Stroke::new(2.0, Color32::YELLOW),
                        );
                        painter.text(
                            p + egui::vec2(9.0, -9.0),
                            egui::Align2::LEFT_BOTTOM,
                            format!("? {}", band.id),
                            egui::FontId::proportional(11.0),
                            Color32::YELLOW,
                        );
                    }
                }
                if let Some(p) = response.hover_pos().and_then(|p| transform.to_original(p)) {
                    response
                        .clone()
                        .on_hover_text(format!("Original pixel: x={:.1}, y={:.1}", p.x, p.y));
                }
                if response.drag_started() {
                    self.drag_start = ui.input(|i| i.pointer.press_origin());
                    if self.tool == MarkTool::Move {
                        self.selected = self
                            .drag_start
                            .and_then(|p| nearest_band(&centers, transform, p));
                    }
                }
                if let (Some(start), Some(end)) = (self.drag_start, response.interact_pointer_pos())
                {
                    if self.tool == MarkTool::Lane {
                        painter.rect_stroke(
                            egui::Rect::from_two_pos(start, end),
                            0.0,
                            egui::Stroke::new(2.0, Color32::WHITE),
                            egui::StrokeKind::Inside,
                        );
                    } else if self.tool == MarkTool::Move && self.selected.is_some() {
                        painter.circle_stroke(end, 7.0, egui::Stroke::new(2.0, Color32::WHITE));
                    }
                }
                if response.drag_stopped() {
                    let points = self
                        .drag_start
                        .take()
                        .zip(response.interact_pointer_pos())
                        .and_then(|(a, b)| transform.to_original(a).zip(transform.to_original(b)));
                    if let Some((a, b)) = points {
                        if self.tool == MarkTool::Lane {
                            self.add_lane(a, b);
                        } else if self.tool == MarkTool::Move {
                            self.move_selected(b);
                        }
                    }
                } else if response.clicked()
                    && let Some(position) = response.interact_pointer_pos()
                {
                    if self.tool == MarkTool::Move {
                        self.selected = nearest_band(&centers, transform, position);
                    } else if let Some(point) = transform.to_original(position) {
                        self.add_band(point);
                    }
                }
                transform
            });
        Some(output.inner)
    }

    fn add_lane(&mut self, a: GelImagePoint, b: GelImagePoint) {
        let Some(draft) = &mut self.draft else {
            return;
        };
        if (a.x - b.x).abs() < 1.0 || (a.y - b.y).abs() < 1.0 || draft.lanes.len() >= 256 {
            return;
        }
        let id = fresh_id("lane", draft.lanes.iter().map(|l| l.id.clone()));
        draft.lanes.push(GelImageLane {
            id: id.clone(),
            label: id.clone(),
            min: GelImagePoint {
                x: a.x.min(b.x),
                y: a.y.min(b.y),
            },
            max: GelImagePoint {
                x: a.x.max(b.x),
                y: a.y.max(b.y),
            },
        });
        if draft.ladder.lane_id.is_empty() {
            draft.ladder.lane_id = id.clone();
        }
        self.selected_lane = id;
    }

    fn add_band(&mut self, center: GelImagePoint) {
        let Some(draft) = &mut self.draft else {
            return;
        };
        match self.tool {
            MarkTool::Reference => {
                if !self.reference_size.is_finite() || self.reference_size <= 0.0 {
                    self.status =
                        "Enter the known positive ladder size before marking its band.".into();
                    return;
                }
                if !draft
                    .lanes
                    .iter()
                    .any(|l| l.id == draft.ladder.lane_id && contains(l, center))
                {
                    self.status = "Mark reference bands inside the ladder lane.".into();
                    return;
                }
                if draft.ladder.bands.len() >= 2048 {
                    return;
                }
                let id = fresh_id("reference", draft.ladder.bands.iter().map(|b| b.id.clone()));
                draft.ladder.bands.push(GelLadderBand {
                    id,
                    center,
                    size: self.reference_size,
                });
                self.reference_size = 0.0;
                self.selected = Some(BandSelection::Reference(draft.ladder.bands.len() - 1));
            }
            MarkTool::Sample => {
                let matching: Vec<_> = draft
                    .lanes
                    .iter()
                    .filter(|l| l.id != draft.ladder.lane_id && contains(l, center))
                    .collect();
                let lane = matching
                    .iter()
                    .find(|l| l.id == self.selected_lane)
                    .copied()
                    .or_else(|| (matching.len() == 1).then(|| matching[0]));
                let Some(lane) = lane else {
                    self.status = "Click inside a sample lane; select the intended lane if rectangles overlap.".into();
                    return;
                };
                if draft.sample_bands.len() >= 2048 {
                    return;
                }
                let id = fresh_id("sample", draft.sample_bands.iter().map(|b| b.id.clone()));
                draft.sample_bands.push(GelSampleBand {
                    id: id.clone(),
                    lane_id: lane.id.clone(),
                    label: id,
                    center,
                    position_half_width_px: None,
                });
                self.selected = Some(BandSelection::Sample(draft.sample_bands.len() - 1));
            }
            _ => {}
        }
    }

    fn move_selected(&mut self, center: GelImagePoint) {
        let Some(draft) = &mut self.draft else {
            return;
        };
        match self.selected {
            Some(BandSelection::Reference(i)) => {
                if draft
                    .lanes
                    .iter()
                    .any(|l| l.id == draft.ladder.lane_id && contains(l, center))
                    && let Some(band) = draft.ladder.bands.get_mut(i)
                {
                    band.center = center;
                }
            }
            Some(BandSelection::Sample(i)) => {
                if let Some(band) = draft.sample_bands.get_mut(i)
                    && draft
                        .lanes
                        .iter()
                        .any(|l| l.id == band.lane_id && contains(l, center))
                {
                    band.center = center;
                }
            }
            None => {}
        }
    }

    fn band_editor(&mut self, ui: &mut Ui) {
        let Some(draft) = &mut self.draft else {
            return;
        };
        ui.collapsing(
            format!(
                "Assignments: {} ladder bands, {} sample bands",
                draft.ladder.bands.len(),
                draft.sample_bands.len()
            ),
            |ui| {
                egui::ScrollArea::vertical()
                    .id_salt("gel_assignments")
                    .max_height(170.0)
                    .show(ui, |ui| {
                        for (i, b) in draft.ladder.bands.iter().enumerate() {
                            if ui
                                .selectable_label(
                                    self.selected == Some(BandSelection::Reference(i)),
                                    format!(
                                        "{}: {} {} at ({:.1}, {:.1})",
                                        b.id,
                                        b.size,
                                        draft.size_kind.unit(),
                                        b.center.x,
                                        b.center.y
                                    ),
                                )
                                .clicked()
                            {
                                self.selected = Some(BandSelection::Reference(i));
                            }
                        }
                        for (i, b) in draft.sample_bands.iter().enumerate() {
                            if ui
                                .selectable_label(
                                    self.selected == Some(BandSelection::Sample(i)),
                                    format!(
                                        "{} / {} at ({:.1}, {:.1})",
                                        b.lane_id, b.label, b.center.x, b.center.y
                                    ),
                                )
                                .clicked()
                            {
                                self.selected = Some(BandSelection::Sample(i));
                            }
                        }
                    });
            },
        );
        ui.horizontal(|ui| {
            match self.selected {
                Some(BandSelection::Reference(i)) if i < draft.ladder.bands.len() => {
                    let band = &mut draft.ladder.bands[i]; ui.label(&band.id);
                    ui.add(egui::DragValue::new(&mut band.size).range(0.0..=1e9).prefix("Size ").suffix(format!(" {}", draft.size_kind.unit())));
                    point_controls(ui, &mut band.center);
                    if ui.button("Delete selected band").clicked() { draft.ladder.bands.remove(i); self.selected = None; }
                }
                Some(BandSelection::Sample(i)) if i < draft.sample_bands.len() => {
                    let band = &mut draft.sample_bands[i]; ui.text_edit_singleline(&mut band.label);
                    point_controls(ui, &mut band.center);
                    let mut has_width = band.position_half_width_px.is_some();
                    if ui.checkbox(&mut has_width, "Localization +/- px").changed() { band.position_half_width_px = has_width.then_some(1.0); }
                    if let Some(width) = &mut band.position_half_width_px { ui.add(egui::DragValue::new(width).range(0.0..=1e6)); }
                    if ui.button("Delete selected band").clicked() { draft.sample_bands.remove(i); self.selected = None; }
                }
                _ => { ui.small("Select a band on the image or in Assignments to correct its position/size or remove it."); }
            }
        });
    }

    fn analysis_controls(&mut self, ui: &mut Ui, shared: &Arc<RwLock<GentleEngine>>) {
        ui.separator();
        if semantic_button(ui, "Analyze and save report", "gel.analyze", false).clicked()
            && let Some(mut request) = self.draft.clone()
        {
            request.report_id = fresh_id(
                "gel_analysis",
                shared
                    .read()
                    .unwrap()
                    .state()
                    .gel_images
                    .analyses
                    .keys()
                    .cloned(),
            );
            self.draft = Some(request.clone());
            self.save_draft(shared);
            self.start(shared, ui.ctx(), GelJob::Analyze(request));
        }
        let Some(report) = self.report.clone() else {
            return;
        };
        let current = self.draft.as_ref() == Some(&report.request);
        ui.label(format!("Saved report: {}", report.request.report_id));
        if !current {
            ui.colored_label(REFERENCE_COLOR, "Assignments changed: the report below belongs to earlier marks. Analyze again before exporting this view.");
        }
        calibration_plot(ui, &report);
        egui::Grid::new("gel_results").striped(true).show(ui, |ui| {
            ui.label("Band");
            ui.label(format!("Size ({})", report.request.size_kind.unit()));
            ui.label("Localization bounds / status");
            ui.end_row();
            for estimate in &report.estimates {
                ui.label(&estimate.label);
                ui.label(
                    estimate
                        .estimated_size
                        .map(|s| format!("{s:.2}"))
                        .unwrap_or_else(|| "Not sized".into()),
                );
                ui.label(match estimate.status {
                    GelBandSizingStatus::OutsideCalibratedRange => {
                        "Outside ladder range; no extrapolation".into()
                    }
                    GelBandSizingStatus::Interpolated => estimate
                        .localization_size_bounds
                        .map(|[low, high]| format!("{low:.2} to {high:.2} (localization only)"))
                        .unwrap_or_else(|| {
                            "Within ladder range; localization bounds not specified/available"
                                .into()
                        }),
                });
                ui.end_row();
            }
        });
        for warning in &report.warnings {
            ui.small(warning);
        }
        for estimate in &report.estimates {
            for warning in &estimate.warnings {
                ui.small(format!("{}: {warning}", estimate.label));
            }
        }
        ui.add_enabled_ui(current && self.job.is_none(), |ui| {
            ui.horizontal(|ui| {
                egui::ComboBox::from_id_salt("gel_export_format")
                    .selected_text(format!("{:?}", self.export_format))
                    .show_ui(ui, |ui| {
                        for format in [
                            GelImageExportFormat::Svg,
                            GelImageExportFormat::Tsv,
                            GelImageExportFormat::Json,
                        ] {
                            ui.selectable_value(
                                &mut self.export_format,
                                format,
                                format!("{format:?}"),
                            );
                        }
                    });
                ui.label("New output path (optional)");
                ui.text_edit_singleline(&mut self.export_path);
                if semantic_button(ui, "Export...", "gel.export", false).clicked() {
                    self.start(
                        shared,
                        ui.ctx(),
                        GelJob::Export(GelImageExportRequest {
                            report_id: report.request.report_id.clone(),
                            path: self.export_path.clone(),
                            format: self.export_format,
                        }),
                    );
                }
            });
        });
    }
}

fn point_controls(ui: &mut Ui, point: &mut GelImagePoint) {
    ui.add(egui::DragValue::new(&mut point.x).speed(0.2).prefix("x "));
    ui.add(egui::DragValue::new(&mut point.y).speed(0.2).prefix("y "));
}

fn nearest_band(
    centers: &[(BandSelection, GelImagePoint, String, Color32)],
    transform: ImageTransform,
    position: Pos2,
) -> Option<BandSelection> {
    centers
        .iter()
        .map(|(id, p, _, _)| (*id, transform.to_screen(*p).distance(position)))
        .filter(|(_, d)| *d <= 12.0)
        .min_by(|a, b| a.1.total_cmp(&b.1))
        .map(|(id, _)| id)
}

fn semantic_button(ui: &mut Ui, label: &str, id: &'static str, selected: bool) -> egui::Response {
    let response = ui.add(egui::Button::new(label).selected(selected));
    #[cfg(feature = "gui-test-support")]
    crate::gui_test_support::register_response(
        &response,
        id,
        "window.gel_image",
        None,
        crate::gui_test_support::GuiTestWidgetKind::Button,
        false,
    );
    let _ = id;
    response
}

/// Visualize the already validated piecewise-log calibration; do not refit it.
fn calibration_plot(ui: &mut Ui, report: &GelImageAnalysisReport) {
    ui.collapsing("Calibration: migration distance vs log10(size)", |ui| {
        let bands = &report.calibration_bands;
        if bands.len() < 2 { return; }
        let (rect, _) = ui.allocate_exact_size(egui::vec2(ui.available_width().min(680.0), 180.0), egui::Sense::hover());
        let plot = rect.shrink2(egui::vec2(65.0, 25.0));
        let migration = report.request.migration;
        let start = migration.coordinate(bands[0].center); let end = migration.coordinate(bands[bands.len()-1].center);
        let high = bands[0].size.log10(); let low = bands[bands.len()-1].size.log10();
        if end <= start || high <= low { return; }
        let points: Vec<_> = bands.iter().map(|b| egui::pos2(
            plot.left() + ((migration.coordinate(b.center)-start)/(end-start)) as f32 * plot.width(),
            plot.bottom() - ((b.size.log10()-low)/(high-low)) as f32 * plot.height(),
        )).collect();
        let painter = ui.painter();
        painter.line_segment([plot.left_top(),plot.left_bottom()], egui::Stroke::new(1.0,Color32::GRAY));
        painter.line_segment([plot.left_bottom(),plot.right_bottom()], egui::Stroke::new(1.0,Color32::GRAY));
        painter.add(egui::Shape::line(points.clone(), egui::Stroke::new(1.5,REFERENCE_COLOR)));
        for point in points { painter.circle_filled(point, 3.0,REFERENCE_COLOR); }
        for (value,y) in [(high,plot.top()),((high+low)/2.0,plot.center().y),(low,plot.bottom())] {
            painter.text(egui::pos2(plot.left()-5.0,y), egui::Align2::RIGHT_CENTER, format!("{value:.2}"), egui::FontId::proportional(11.0),ui.visuals().text_color());
        }
        painter.text(egui::pos2(plot.left(),plot.bottom()+6.0), egui::Align2::LEFT_TOP, "0 px", egui::FontId::proportional(11.0),ui.visuals().text_color());
        painter.text(egui::pos2(plot.right(),plot.bottom()+6.0), egui::Align2::RIGHT_TOP, format!("{:.1} px from first reference",end-start), egui::FontId::proportional(11.0),ui.visuals().text_color());
        ui.small(format!("Y: log10({}); X: pixels along migration. Straight segments interpolate between confirmed references only.", report.request.size_kind.unit()));
    });
}

impl GENtleApp {
    pub(super) fn gel_image_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("GENtle Measured Gel Image")
    }
    pub(super) fn gel_image_window_id() -> egui::Id {
        egui::Id::new("hosted_gel_image_editor")
    }

    pub(super) fn open_gel_image_editor(&mut self) {
        self.gel_image_editor.bind_project(&self.engine);
        let was_open = self.gel_image_editor.open;
        self.gel_image_editor.open = true;
        self.mark_window_open_or_focus(Self::gel_image_viewport_id(), was_open);
    }

    pub(super) fn render_gel_image_editor(&mut self, ctx: &egui::Context) {
        self.gel_image_editor.poll(&self.engine, ctx);
        if !self.gel_image_editor.open {
            return;
        }
        let viewport_id = Self::gel_image_viewport_id();
        let spec = self.hosted_window_spec_for_viewport(
            "Gel Image Analysis",
            Self::gel_image_window_id(),
            viewport_id,
            egui::vec2(1120.0, 900.0),
            egui::vec2(780.0, 600.0),
        );
        let mut open = true;
        let engine = self.engine.clone();
        if ctx.embed_viewports() {
            crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
                egui::ScrollArea::vertical()
                    .id_salt("gel_editor_scroll")
                    .show(ui, |ui| self.gel_image_editor.contents(ui, &engine));
            });
            self.clear_viewport_foreground_request_after_render(viewport_id);
        } else {
            let builder = crate::egui_compat::viewport_builder_for_hosted_window(&spec);
            ctx.show_viewport_immediate(viewport_id, builder, |ctx, class| {
                self.note_viewport_focus_if_active(ctx, viewport_id);
                if class == egui::ViewportClass::EmbeddedWindow {
                    crate::egui_compat::show_hosted_window(&mut *ctx, &spec, &mut open, |ui| {
                        egui::ScrollArea::vertical()
                            .id_salt("gel_editor_scroll")
                            .show(ui, |ui| self.gel_image_editor.contents(ui, &engine));
                    });
                } else {
                    crate::egui_compat::show_central_panel(
                        &mut *ctx,
                        egui::CentralPanel::default(),
                        |ui| {
                            egui::ScrollArea::vertical()
                                .id_salt("gel_editor_scroll")
                                .show(ui, |ui| self.gel_image_editor.contents(ui, &engine));
                        },
                    );
                }
                if Self::viewport_close_requested_or_shortcut(ctx) {
                    open = false;
                }
            });
        }
        self.gel_image_editor.open = open;
    }
}

#[cfg(test)]
#[path = "gel_image_ui_tests.rs"]
mod tests;
