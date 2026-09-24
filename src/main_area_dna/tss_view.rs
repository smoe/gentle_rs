//! Native TSS lanes over the shared annotated-window presentation model.

use super::*;
use crate::tss_sequence_view::{TssLaneKind, TssSequenceView, TssViewFeature, TssViewLane};
use std::sync::mpsc::{Receiver, TryRecvError};

mod local_scoring;
use local_scoring::LocalScoreState;

type LoadedView = Result<Arc<TssSequenceView>, String>;
pub(super) type SvgExportReceiver = Arc<Mutex<Receiver<Result<String, String>>>>;

#[derive(Clone, Debug)]
struct ProfileLoad {
    source: Arc<TssSequenceView>,
    receiver: Arc<std::sync::Mutex<Receiver<LoadedView>>>,
}

#[derive(Clone, Debug)]
pub(super) struct TssUiState {
    key: Option<(u64, usize)>,
    detected: bool,
    pending: Option<Arc<std::sync::Mutex<Receiver<LoadedView>>>>,
    document: Option<LoadedView>,
    annotations: Option<Arc<TssSequenceView>>,
    profile_load: Option<ProfileLoad>,
    requested_profile_path: Option<std::path::PathBuf>,
    profile_error: Option<String>,
    local_scores: LocalScoreState,
    structures: bool,
    signals: bool,
    motifs: bool,
    traces: bool,
    local_traces: bool,
    imported: bool,
    other: bool,
    filter: String,
    selected: Option<TssViewFeature>,
}

impl Default for TssUiState {
    fn default() -> Self {
        Self {
            key: None,
            detected: false,
            pending: None,
            document: None,
            annotations: None,
            profile_load: None,
            requested_profile_path: None,
            profile_error: None,
            local_scores: LocalScoreState::default(),
            structures: true,
            signals: true,
            motifs: true,
            traces: true,
            local_traces: true,
            imported: true,
            other: true,
            filter: String::new(),
            selected: None,
        }
    }
}

impl TssUiState {
    fn visible_lanes(&self, view: &TssSequenceView) -> Vec<usize> {
        let filter = self.filter.to_lowercase();
        view.lanes
            .iter()
            .enumerate()
            .filter(|(_, lane)| {
                let enabled = match lane.kind {
                    TssLaneKind::Structure => self.structures,
                    TssLaneKind::Signal => self.signals,
                    TssLaneKind::Motif => self.motifs,
                    TssLaneKind::ScoreTrace => self.traces,
                    TssLaneKind::LocalScoreTrace => self.local_traces,
                    TssLaneKind::ImportedMotif => self.imported,
                    TssLaneKind::Other => self.other,
                };
                enabled
                    && (filter.is_empty()
                        || lane.label.to_lowercase().contains(&filter)
                        || lane.id.to_lowercase().contains(&filter))
            })
            .map(|(i, _)| i)
            .collect()
    }
}

fn color(kind: TssLaneKind) -> egui::Color32 {
    match kind {
        TssLaneKind::Structure => egui::Color32::from_rgb(55, 130, 170),
        TssLaneKind::Signal => egui::Color32::from_rgb(165, 65, 100),
        TssLaneKind::Motif => egui::Color32::from_rgb(0, 140, 115),
        TssLaneKind::ScoreTrace => egui::Color32::from_rgb(30, 105, 185),
        TssLaneKind::LocalScoreTrace => egui::Color32::from_rgb(0, 130, 100),
        TssLaneKind::ImportedMotif => egui::Color32::from_rgb(180, 100, 25),
        TssLaneKind::Other => egui::Color32::from_rgb(180, 125, 30),
    }
}

impl MainAreaDna {
    fn tss_svg_snapshot(
        &mut self,
        profile: ViewSvgExportProfile,
    ) -> Result<
        (
            Arc<TssSequenceView>,
            crate::tss_sequence_view::TssViewSvgOptions,
        ),
        String,
    > {
        self.refresh_tss_recognition();
        if self.tss_ui.profile_load.is_some()
            || self.tss_ui.pending.is_some()
            || self.tss_ui.requested_profile_path.is_some()
        {
            return Err("Wait for the TSS document/report to finish loading before export".into());
        }
        if self.tss_ui.local_scores.running() {
            return Err("Wait for local TSS scoring, or cancel it, before export".into());
        }
        let view = self
            .tss_ui
            .document
            .as_ref()
            .ok_or("TSS document is not ready")?
            .as_ref()
            .map_err(Clone::clone)?
            .clone();
        let (start, span, length) = self.current_linear_viewport();
        let layout =
            Self::view_svg_export_layout(profile, self.last_linear_map_width_px, true, false);
        let (start, span) =
            Self::expand_linear_export_window(start, span, length, layout.viewport_span_multiplier);
        let options = crate::tss_sequence_view::TssViewSvgOptions {
            start_0based: start,
            end_0based_exclusive: start.saturating_add(span).min(length),
            lane_indices: self.tss_ui.visible_lanes(&view),
            width_px: layout.canvas_width_px as u32,
            print_size_mm: layout.print_size_mm,
        };
        Ok((view, options))
    }

    pub(super) fn export_tss_view_svg(
        &mut self,
        profile: ViewSvgExportProfile,
        ctx: &egui::Context,
    ) {
        if self.tss_svg_export.is_some() {
            self.op_status = "A TSS SVG export is already running".into();
            return;
        }
        let (view, options) = match self.tss_svg_snapshot(profile) {
            Ok(snapshot) => snapshot,
            Err(e) => {
                self.op_status = e;
                return;
            }
        };
        let Some(path) = rfd::FileDialog::new()
            .set_file_name(format!("tss.{}.svg", profile.file_stem_suffix()))
            .add_filter("SVG", &["svg"])
            .save_file()
        else {
            self.op_status = "Export canceled".into();
            return;
        };
        let (send, receiver) = std::sync::mpsc::channel();
        self.tss_svg_export = Some(Arc::new(Mutex::new(receiver)));
        self.op_status = "Exporting the selected TSS view snapshot...".into();
        let ctx = ctx.clone();
        std::thread::spawn(move || {
            let result =
                crate::tss_sequence_view::write_tss_view_svg(&view, &options, &path).map(|hash| {
                    format!(
                        "Exported TSS SVG snapshot: {} (SHA-256 {hash})",
                        path.display()
                    )
                });
            let _ = send.send(result);
            ctx.request_repaint();
        });
    }

    pub(super) fn poll_tss_svg_export(&mut self, ctx: &egui::Context) {
        let Some(receiver) = &self.tss_svg_export else {
            return;
        };
        let result = receiver
            .lock()
            .map_err(|_| TryRecvError::Disconnected)
            .and_then(|r| r.try_recv());
        match result {
            Ok(result) => {
                self.op_status = result.unwrap_or_else(|e| format!("TSS SVG export failed: {e}"));
                self.tss_svg_export = None;
            }
            Err(TryRecvError::Empty) => {
                ctx.request_repaint_after(std::time::Duration::from_millis(100))
            }
            Err(_) => {
                self.op_status = "TSS SVG export worker stopped".into();
                self.tss_svg_export = None;
            }
        }
    }

    pub(super) fn tss_display_title(&self) -> Option<&str> {
        self.tss_ui
            .document
            .as_ref()?
            .as_ref()
            .ok()
            .map(|view| view.title.as_str())
    }

    fn refresh_tss_recognition(&mut self) {
        let Ok(dna) = self.dna.read() else { return };
        let key = (dna.feature_generation(), dna.len());
        if self.tss_ui.key != Some(key) {
            self.tss_ui = TssUiState {
                key: Some(key),
                detected: TssSequenceView::recognizes(&dna),
                ..Default::default()
            };
        }
    }

    pub(super) fn tss_view_available(&mut self) -> bool {
        self.refresh_tss_recognition();
        self.tss_ui.detected
    }

    /// Shared UI intents and the toolbar enter the same read-only DNA display mode.
    pub(crate) fn set_tss_view(&mut self, enabled: bool) -> Result<(), String> {
        if !enabled && self.primary_map_mode != PrimaryMapMode::Tss {
            return Err("TSS view is not active; the current map was left unchanged".into());
        }
        if enabled && (!self.tss_view_available() || self.is_circular()) {
            return Err("Active DNA is not a linear GENtle annotated TSS window. Open its EMBL/GenBank export first.".into());
        }
        self.primary_map_mode = if enabled {
            PrimaryMapMode::Tss
        } else {
            PrimaryMapMode::Standard
        };
        self.save_engine_ops_state();
        Ok(())
    }

    /// Queue the same profile attachment used by the native file-picker path.
    ///
    /// Validation remains inside `TssSequenceView::load_profile`; this adapter
    /// neither scores motifs nor queries external evidence.
    pub(crate) fn queue_tss_profile(&mut self, path: std::path::PathBuf) -> Result<(), String> {
        if path.as_os_str().is_empty() {
            return Err("TSS profile report path must not be empty".into());
        }
        self.set_tss_view(true)?;
        self.tss_ui.requested_profile_path = Some(path);
        Ok(())
    }

    fn poll_tss_view(&mut self, ctx: &egui::Context) {
        self.refresh_tss_recognition();
        if self.tss_ui.document.is_none() && self.tss_ui.pending.is_none() {
            // No engine lock, scoring, file I/O, or feature decoding in the paint callback.
            let dna = self.dna.clone();
            let ctx = ctx.clone();
            let (send, recv) = std::sync::mpsc::channel();
            self.tss_ui.pending = Some(Arc::new(std::sync::Mutex::new(recv)));
            std::thread::spawn(move || {
                let result = dna
                    .read()
                    .map_err(|_| "Could not read TSS sequence".to_string())
                    .and_then(|dna| TssSequenceView::from_dna(&dna))
                    .map(Arc::new);
                let _ = send.send(result);
                ctx.request_repaint();
            });
        }
        if let Some(receiver) = self.tss_ui.pending.clone() {
            let result = receiver
                .lock()
                .map_err(|_| TryRecvError::Disconnected)
                .and_then(|r| r.try_recv());
            match result {
                Ok(document) => {
                    self.tss_ui.annotations = document.as_ref().ok().cloned();
                    self.tss_ui.document = Some(document);
                    self.tss_ui.pending = None;
                }
                Err(TryRecvError::Disconnected) => {
                    self.tss_ui.document = Some(Err("TSS inspection worker stopped".into()));
                    self.tss_ui.pending = None;
                }
                Err(TryRecvError::Empty) => {
                    ctx.request_repaint_after(std::time::Duration::from_millis(100));
                }
            }
        }
    }

    fn load_tss_profile(&mut self, path: std::path::PathBuf, ctx: &egui::Context) {
        let Some(Ok(source)) = self.tss_ui.document.clone() else {
            return;
        };
        // Attaching/detaching a report must preserve independently computed lanes.
        let annotations = source.clone();
        let (send, receiver) = std::sync::mpsc::channel();
        self.tss_ui.profile_error = None;
        self.tss_ui.profile_load = Some(ProfileLoad {
            source,
            receiver: Arc::new(std::sync::Mutex::new(receiver)),
        });
        let ctx = ctx.clone();
        std::thread::spawn(move || {
            let _ = send.send(annotations.load_profile(&path).map(Arc::new));
            ctx.request_repaint();
        });
    }

    fn poll_tss_profile(&mut self, ctx: &egui::Context) {
        let Some(load) = self.tss_ui.profile_load.clone() else {
            return;
        };
        let result = load
            .receiver
            .lock()
            .map_err(|_| TryRecvError::Disconnected)
            .and_then(|r| r.try_recv());
        if matches!(result, Err(TryRecvError::Empty)) {
            ctx.request_repaint_after(std::time::Duration::from_millis(100));
            return;
        }
        self.tss_ui.profile_load = None;
        if !self
            .tss_ui
            .document
            .as_ref()
            .and_then(|d| d.as_ref().ok())
            .is_some_and(|current| Arc::ptr_eq(current, &load.source))
        {
            return; // A replaced/edited document cannot receive an old attachment.
        }
        match result {
            Ok(Ok(view)) => {
                self.tss_ui.document = Some(Ok(view));
                self.tss_ui.selected = None;
            }
            Ok(Err(error)) => self.tss_ui.profile_error = Some(error),
            Err(_) => self.tss_ui.profile_error = Some("TSS report worker stopped".into()),
        }
    }

    pub(super) fn render_primary_tss_map_ui(&mut self, ui: &mut egui::Ui) {
        // The primary-map host may use a horizontal layout for its other map
        // implementations.  The TSS evidence view is a document-like stack;
        // establish that layout explicitly so controls and lanes cannot be
        // pushed beyond the right edge of the native viewport.
        ui.vertical(|ui| self.render_primary_tss_map_contents(ui));
    }

    fn render_primary_tss_map_contents(&mut self, ui: &mut egui::Ui) {
        self.poll_tss_view(ui.ctx());
        self.poll_tss_profile(ui.ctx());
        self.poll_tss_local_scores(ui.ctx());
        if self.tss_ui.profile_load.is_none()
            && !self.tss_ui.local_scores.running()
            && self
                .tss_ui
                .document
                .as_ref()
                .is_some_and(|document| document.is_ok())
            && let Some(path) = self.tss_ui.requested_profile_path.take()
        {
            self.load_tss_profile(path, ui.ctx());
        }
        let Some(document) = self.tss_ui.document.clone() else {
            ui.spinner();
            ui.label("Checking TSS sequence binding and grouping evidence...");
            return;
        };
        let view = match document {
            Ok(view) => view,
            Err(error) => {
                ui.colored_label(ui.visuals().error_fg_color, error);
                ui.label("The standard DNA map remains available; no evidence has been changed.");
                return;
            }
        };
        ui.heading(&view.title);
        ui.small("Annotation-derived TSS | transcript-oriented genomic DNA, not a spliced transcript or reporter construct");
        ui.horizontal_wrapped(|ui| {
            if ui.add_enabled(self.tss_ui.profile_load.is_none() && !self.tss_ui.local_scores.running(), egui::Button::new("Load TSS profile report...")).on_hover_text("Select report.json from the same TSS SVG/PDF bundle. Validates the report, exact reference, TSS geometry and sequence hash; does not rescore or query DuckDB.").clicked()
                && let Some(path) = rfd::FileDialog::new().add_filter("TSS profile JSON", &["json"]).pick_file()
            {
                self.load_tss_profile(path, ui.ctx());
            }
            if self.tss_ui.profile_load.is_some() {
                ui.spinner();
                ui.label("Validating report and preparing evidence...");
                if ui.button("Cancel attachment").on_hover_text("Ignore this worker's result; its bounded file read may still finish in the background.").clicked() {
                    self.tss_ui.profile_load = None;
                }
            } else if view.profile.is_some() && ui.add_enabled(!self.tss_ui.local_scores.running(), egui::Button::new("Detach report")).clicked() {
                let mut detached = (*view).clone();
                detached.clear_profile();
                self.tss_ui.document = Some(Ok(Arc::new(detached)));
                self.tss_ui.selected = None;
                self.tss_ui.profile_error = None;
            }
        });
        if let Some(error) = &self.tss_ui.profile_error {
            ui.colored_label(ui.visuals().error_fg_color, error);
        }
        if let Some(profile) = &view.profile {
            ui.small(format!(
                "Report panel: {} | producer {} | no new scoring or database query",
                profile.panel_id, profile.producer_revision
            ));
        } else {
            ui.small("Annotated file only: load report.json for complete TFBS curves and attached DuckDB hits. These cannot be reconstructed from stored peaks.");
        }
        self.render_tss_local_scoring(ui, &view);
        ui.horizontal_wrapped(|ui| {
            ui.checkbox(&mut self.tss_ui.structures, "Exons / CDS");
            ui.checkbox(&mut self.tss_ui.signals, "CUT&RUN / chromatin");
            ui.checkbox(&mut self.tss_ui.motifs, "Stored motif peaks");
            ui.checkbox(&mut self.tss_ui.traces, "Report TFBS curves");
            ui.checkbox(&mut self.tss_ui.local_traces, "Locally computed curves");
            ui.checkbox(&mut self.tss_ui.imported, "DuckDB peaks");
            ui.checkbox(&mut self.tss_ui.other, "Other annotations");
            ui.label("Filter lanes");
            ui.add(egui::TextEdit::singleline(&mut self.tss_ui.filter).desired_width(150.0));
        });
        ui.small("Rose: supplied signal; green: stored peaks (negative values hidden); blue curves: report scores (+ solid / - dashed); amber triangles: imported raw scores (+ up / - down). Local and imported scores use separate scales. Click evidence to select its DNA span.");
        ui.collapsing("Provenance, limitations and missing data", |ui| {
            for warning in &view.warnings {
                ui.label(warning);
            }
            ui.label(&view.provenance);
            if let Some(local) = &view.local_scoring {
                ui.label(format!("Locally computed: {} | report SHA-256 {} | source {} | cache {} | producer {}", local.request.score_kind.as_str(), local.report_sha256, local.source_binding_sha256, local.cache_key_sha256, local.producer_revision));
                ui.label("Model scores are not experimental binding or luciferase activity. Missing/ambiguous windows are gaps, not zero scores. No network or DuckDB query.");
            }
            if let Some(profile) = &view.profile {
                ui.label(format!("Report file SHA-256: {}", profile.file_sha256));
                ui.label("Validated report contents and sequence binding, not a full receipt audit or independent reference authentication. Original annotated features remain unchanged; reporter constructs are not supplied by this attachment.");
                for warning in &profile.warnings {
                    ui.label(warning);
                }
            }
        });

        let lanes = self.tss_ui.visible_lanes(&view);
        ui.small(format!("{} / {} lanes visible. Hover for source intervals; displayed interval ends are not individual read ends.", lanes.len(), view.lanes.len()));
        let (start, span, _) = self.current_linear_viewport();
        let start = start.min(view.geometry.length().unwrap() - 1);
        let end = (start + span)
            .min(view.geometry.length().unwrap())
            .max(start + 1);
        let width = ui.available_width().max(240.0);
        let left = (width * 0.23).clamp(90.0, 190.0);
        let right = (width * 0.25).clamp(90.0, 240.0);
        let (axis, _) = ui.allocate_exact_size(egui::vec2(width, 60.0), egui::Sense::hover());
        let plot_left = axis.left() + left;
        let plot_right = (axis.right() - right).max(plot_left + 50.0);
        let x = |pos: usize| {
            plot_left
                + (pos as f64 - start as f64) as f32 / (end - start) as f32
                    * (plot_right - plot_left)
        };
        let painter = ui.painter();
        let text_color = ui.visuals().text_color();
        painter.text(
            axis.left_top(),
            egui::Align2::LEFT_TOP,
            "TSS-relative bp\nGenomic base\nLocal base",
            egui::FontId::monospace(10.0),
            text_color,
        );
        let mut ticks = vec![];
        if (start..end).contains(&view.geometry.upstream_bp) {
            ticks.push(view.geometry.upstream_bp);
        }
        for pos in [start, end - 1] {
            if ticks
                .iter()
                .all(|existing| (x(*existing) - x(pos)).abs() >= 70.0)
            {
                ticks.push(pos);
            }
        }
        ticks.sort_unstable();
        ticks.dedup();
        for pos in ticks {
            let align = if pos == start {
                egui::Align2::LEFT_TOP
            } else if pos == end - 1 {
                egui::Align2::RIGHT_TOP
            } else {
                egui::Align2::CENTER_TOP
            };
            painter.text(
                egui::pos2(x(pos), axis.top()),
                align,
                format!(
                    "{:+}\n{}\n{}",
                    pos as i64 - view.geometry.upstream_bp as i64,
                    view.geometry.genomic_at(pos).unwrap(),
                    pos + 1
                ),
                egui::FontId::monospace(10.0),
                text_color,
            );
        }
        let mut selected = None;
        egui::ScrollArea::vertical()
            .id_salt(("tss_lanes", self.panel_scope_key()))
            .max_height((ui.available_height() - 115.0).max(120.0))
            .show_rows(ui, 148.0, lanes.len(), |ui, rows| {
                for row in rows {
                    let index = lanes[row];
                    if let Some(feature) = paint_lane(
                        ui,
                        &view,
                        &view.lanes[index],
                        width,
                        left,
                        right,
                        start,
                        end,
                    ) {
                        selected = Some(feature);
                    }
                }
            });
        if let Some(f) = selected {
            if let Err(error) = self.set_selection_range_0based(f.start, f.end) {
                self.op_status = error;
            }
            self.tss_ui.selected = Some(f);
        }
        if let Some(f) = self.tss_ui.selected.clone() {
            ui.separator();
            ui.label(format!(
                "{} | {} to {}{}",
                f.label,
                view.coordinate_label(f.start),
                view.coordinate_label(f.end - 1),
                if f.clipped { " | clipped" } else { "" }
            ));
            ui.horizontal(|ui| {
                if ui.button("Inspect in standard DNA map").clicked() {
                    self.primary_map_mode = PrimaryMapMode::Standard;
                    self.show_sequence = true;
                    self.map_dna.select_feature(f.feature_id);
                    self.save_engine_ops_state();
                }
                ui.label("Full annotation").on_hover_text(&f.details);
            });
        }
    }
}

fn paint_lane(
    ui: &mut egui::Ui,
    view: &TssSequenceView,
    lane: &TssViewLane,
    width: f32,
    left: f32,
    right: f32,
    start: usize,
    end: usize,
) -> Option<TssViewFeature> {
    let (rect, response) = ui.allocate_exact_size(egui::vec2(width, 144.0), egui::Sense::click());
    let plot = egui::Rect::from_min_max(
        egui::pos2(rect.left() + left, rect.top() + 9.0),
        egui::pos2(
            (rect.right() - right).max(rect.left() + left + 50.0),
            rect.bottom() - 12.0,
        ),
    );
    let painter = ui.painter().with_clip_rect(rect.intersect(ui.clip_rect()));
    let foreground = ui.visuals().text_color();
    let font = egui::FontId::proportional(11.0);
    let label = painter.layout(
        lane.label.clone(),
        font.clone(),
        foreground,
        (left - 45.0).max(50.0),
    );
    painter.galley(rect.left_top(), label, foreground);
    painter.rect_filled(plot, 2.0, ui.visuals().faint_bg_color);
    if lane.trace.is_some() {
        return paint_trace(ui, view, lane, plot, rect, response, start, end);
    }
    let imported = lane.kind == TssLaneKind::ImportedMotif;
    let middle = lane.kind == TssLaneKind::Motif
        || imported
        || lane
            .features
            .iter()
            .any(|f| f.score.is_some_and(|s| s < 0.0));
    let baseline = if middle {
        plot.center().y
    } else {
        plot.bottom()
    };
    let amplitude = if middle {
        plot.height() / 2.0
    } else {
        plot.height()
    };
    let x = |p: usize| {
        plot.left() + (p as f64 - start as f64) as f32 / (end - start) as f32 * plot.width()
    };
    painter.line_segment(
        [
            egui::pos2(plot.left(), baseline),
            egui::pos2(plot.right(), baseline),
        ],
        egui::Stroke::new(0.5, foreground),
    );
    if lane.kind == TssLaneKind::Signal || lane.kind == TssLaneKind::Motif || imported {
        painter.text(
            plot.left_top() - egui::vec2(3.0, 0.0),
            egui::Align2::RIGHT_TOP,
            format!("{:.2}", lane.scale_max),
            egui::FontId::monospace(9.0),
            foreground,
        );
        painter.text(
            egui::pos2(plot.left() - 3.0, baseline),
            egui::Align2::RIGHT_BOTTOM,
            format!("{:.2}", if imported { lane.scale_min } else { 0.0 }),
            egui::FontId::monospace(9.0),
            foreground,
        );
        if middle {
            let bottom_max = if lane.kind == TssLaneKind::Motif || imported {
                lane.scale_max
            } else {
                -lane.scale_max
            };
            painter.text(
                plot.left_bottom() - egui::vec2(3.0, 0.0),
                egui::Align2::RIGHT_BOTTOM,
                format!("{bottom_max:.2}"),
                egui::FontId::monospace(9.0),
                foreground,
            );
        }
    }
    let tss = view.geometry.upstream_bp;
    if (start..end).contains(&tss) {
        painter.line_segment(
            [
                egui::pos2(x(tss), plot.top()),
                egui::pos2(x(tss), plot.bottom()),
            ],
            egui::Stroke::new(1.0, foreground),
        );
    }
    let p = ui.painter().with_clip_rect(plot.intersect(ui.clip_rect()));
    let pointer = response.hover_pos();
    let mut hit = None;
    let mut visible = 0;
    let mut shown_bounds: Option<(usize, usize)> = None;
    for (i, f) in lane.features.iter().enumerate() {
        if f.end <= start
            || f.start >= end
            || (lane.kind == TssLaneKind::Motif && !f.score.is_some_and(|v| v >= 0.0))
        {
            continue;
        }
        visible += 1;
        // Bound paint work without fabricating or interpolating the remaining evidence.
        if visible > 10_000 {
            continue;
        }
        let bounds = shown_bounds.get_or_insert((f.start.max(start), f.end.min(end) - 1));
        bounds.0 = bounds.0.min(f.start.max(start));
        bounds.1 = bounds.1.max(f.end.min(end) - 1);
        let mut feature_rect;
        if lane.kind == TssLaneKind::Motif || imported {
            let height = if imported {
                ((f.score.unwrap_or(lane.scale_min) - lane.scale_min)
                    / (lane.scale_max - lane.scale_min)) as f32
                    * amplitude
            } else {
                (f.score.unwrap_or(0.0) / lane.scale_max) as f32 * amplitude
            }
            .max(2.0);
            let y = baseline + if f.reverse { height } else { -height };
            let pos = egui::pos2(
                if imported {
                    (x(f.start) + x(f.end)) / 2.0
                } else {
                    x(f.start)
                },
                y,
            );
            let offset = if f.reverse { -4.0 } else { 4.0 };
            p.line_segment(
                [egui::pos2(pos.x, baseline), pos],
                egui::Stroke::new(0.6, color(lane.kind).gamma_multiply(0.5)),
            );
            let (a, b) = if imported {
                (x(f.start), x(f.end))
            } else {
                (pos.x - 3.0, pos.x + 3.0)
            };
            p.add(egui::Shape::convex_polygon(
                vec![
                    pos,
                    egui::pos2(a, if imported { baseline } else { pos.y + offset }),
                    egui::pos2(
                        b.max(a + 2.0),
                        if imported { baseline } else { pos.y + offset },
                    ),
                ],
                color(lane.kind),
                egui::Stroke::NONE,
            ));
            feature_rect = if imported {
                egui::Rect::from_min_max(
                    egui::pos2(a, y.min(baseline)),
                    egui::pos2(b.max(a + 2.0), y.max(baseline)),
                )
            } else {
                egui::Rect::from_center_size(pos, egui::vec2(8.0, 10.0))
            };
        } else if lane.kind == TssLaneKind::Signal {
            let y = baseline - (f.score.unwrap_or(0.0) / lane.scale_max) as f32 * amplitude;
            feature_rect = egui::Rect::from_min_max(
                egui::pos2(x(f.start), y.min(baseline)),
                egui::pos2(x(f.end), y.max(baseline) + 1.0),
            );
            if f.score.is_some() {
                p.rect_filled(feature_rect, 0.0, color(lane.kind).gamma_multiply(0.7));
            } else {
                p.line_segment(
                    [feature_rect.left_bottom(), feature_rect.right_bottom()],
                    egui::Stroke::new(2.0, color(TssLaneKind::Other)),
                );
            }
        } else {
            let thick = f.label.contains("CDS segment");
            feature_rect = egui::Rect::from_min_max(
                egui::pos2(x(f.start), plot.center().y - if thick { 10.0 } else { 5.0 }),
                egui::pos2(
                    x(f.end).max(x(f.start) + 2.0),
                    plot.center().y + if thick { 10.0 } else { 5.0 },
                ),
            );
            p.rect_filled(
                feature_rect,
                1.0,
                color(lane.kind).gamma_multiply(if thick { 1.0 } else { 0.6 }),
            );
        }
        feature_rect = feature_rect.intersect(plot);
        if pointer.is_some_and(|pos| feature_rect.expand(2.0).contains(pos)) {
            hit = Some(i);
        }
    }
    let count = if visible > 10_000 {
        format!("10,000/{visible} drawn; zoom in")
    } else {
        format!("{visible}/{} intervals", lane.features.len())
    };
    let summary = if visible == 0 {
        format!(
            "{}\nNo displayed intervals here\nNot a measured zero",
            lane.state
        )
    } else {
        let (first, last) = shown_bounds.unwrap();
        let units = match lane.units.as_str() {
            "llr_background_tail_log10" => "-log10 background tail P",
            "llr_log2" => "LLR (bits)",
            "source signal (not read endpoints)" => "source signal",
            "exon / coding segment / translation marker" => "exon / CDS context",
            _ => &lane.units,
        };
        format!(
            "{count}\n{units}\ngenomic {}..{}\nlocal {}..{}\n{}",
            view.geometry.genomic_at(first).unwrap(),
            view.geometry.genomic_at(last).unwrap(),
            first + 1,
            last + 1,
            lane.state
        )
    };
    let text = painter.layout(
        summary,
        egui::FontId::proportional(10.0),
        foreground,
        right - 8.0,
    );
    painter.galley(egui::pos2(plot.right() + 5.0, rect.top()), text, foreground);
    let clicked = response.clicked();
    if let Some(index) = hit {
        let f = &lane.features[index];
        response.on_hover_text(format!(
            "{}\n{}\n{}\n{}\n{}",
            f.label,
            view.coordinate_label(f.start),
            view.coordinate_label(f.end - 1),
            f.details,
            lane.details
        ));
        if clicked {
            return Some(f.clone());
        }
    } else {
        response.on_hover_text(format!("{}\n{}\n{}", lane.id, lane.units, lane.details));
    }
    None
}

/// Curves use window starts on the existing local axis, never strand-flipped x values.
fn paint_trace(
    ui: &mut egui::Ui,
    view: &TssSequenceView,
    lane: &TssViewLane,
    plot: egui::Rect,
    row: egui::Rect,
    response: egui::Response,
    start: usize,
    end: usize,
) -> Option<TssViewFeature> {
    let trace = lane.trace.as_ref()?;
    let painter = ui.painter().with_clip_rect(row.intersect(ui.clip_rect()));
    let p = painter.with_clip_rect(plot.intersect(ui.clip_rect()));
    let foreground = ui.visuals().text_color();
    let x = |pos: usize| {
        plot.left() + (pos as f64 - start as f64) as f32 / (end - start) as f32 * plot.width()
    };
    let y = |value: f64| {
        plot.bottom()
            - ((value - lane.scale_min) / (lane.scale_max - lane.scale_min)) as f32 * plot.height()
    };
    for value in [
        lane.scale_min,
        (lane.scale_min + lane.scale_max) / 2.0,
        lane.scale_max,
    ] {
        painter.text(
            egui::pos2(plot.left() - 3.0, y(value)),
            egui::Align2::RIGHT_CENTER,
            format!("{value:.2}"),
            egui::FontId::monospace(9.0),
            foreground,
        );
        p.line_segment(
            [
                egui::pos2(plot.left(), y(value)),
                egui::pos2(plot.right(), y(value)),
            ],
            egui::Stroke::new(0.4, foreground.gamma_multiply(0.3)),
        );
    }
    let visible_end = end.min(trace.forward.len());
    let bounded = visible_end.saturating_sub(start) <= 10_000;
    if visible_end < end {
        p.rect_filled(
            egui::Rect::from_min_max(
                egui::pos2(x(visible_end.max(start)), plot.top()),
                plot.right_bottom(),
            ),
            0.0,
            egui::Color32::GRAY.gamma_multiply(0.22),
        );
    }
    let mut valid = 0;
    if bounded {
        for (reverse, scores) in [(false, &trace.forward), (true, &trace.reverse)] {
            let stroke = egui::Stroke::new(
                1.2,
                if reverse {
                    egui::Color32::from_rgb(170, 70, 105)
                } else {
                    color(lane.kind)
                },
            );
            let mut path = Vec::new();
            let flush = |path: &mut Vec<egui::Pos2>| {
                if path.len() == 1 {
                    p.circle_filled(path[0], 1.5, stroke.color);
                } else if reverse {
                    p.extend(egui::Shape::dashed_line(path, stroke, 5.0, 3.0));
                } else if !path.is_empty() {
                    p.add(egui::Shape::line(path.clone(), stroke));
                }
                path.clear();
            };
            for pos in start..visible_end {
                match scores[pos] {
                    Some(raw) => {
                        valid += 1;
                        let display = if trace.clip_negative {
                            raw.max(0.0)
                        } else {
                            raw
                        };
                        path.push(egui::pos2(x(pos), y(display)));
                    }
                    None => {
                        flush(&mut path);
                        let (top, bottom) = if reverse {
                            (plot.center().y, plot.bottom())
                        } else {
                            (plot.top(), plot.center().y)
                        };
                        p.rect_filled(
                            egui::Rect::from_min_max(
                                egui::pos2(x(pos), top),
                                egui::pos2(x(pos + 1), bottom),
                            ),
                            0.0,
                            egui::Color32::from_rgb(215, 163, 55).gamma_multiply(0.3),
                        );
                    }
                }
            }
            flush(&mut path);
        }
    }
    if (start..end).contains(&view.geometry.upstream_bp) {
        p.line_segment(
            [
                egui::pos2(x(view.geometry.upstream_bp), plot.top()),
                egui::pos2(x(view.geometry.upstream_bp), plot.bottom()),
            ],
            egui::Stroke::new(1.0, foreground),
        );
    }
    let status = if !bounded {
        "More than 10,000 window starts visible; zoom in to draw curves".to_string()
    } else if valid == 0 {
        "No evaluable windows here; NOT zero signal".to_string()
    } else if trace.range_is_fallback {
        format!("{valid} valid strand-windows; displayed values are zero")
    } else {
        format!("{valid} valid strand-windows here")
    };
    let summary = format!(
        "{}\n{status}{}\n+ solid / - dashed\nAmber: unavailable; grey: no full motif window\n{}",
        lane.units,
        if trace.range_is_fallback {
            "; axis 0..1 is a display fallback"
        } else {
            ""
        },
        lane.state
    );
    let galley = painter.layout(
        summary,
        egui::FontId::proportional(10.0),
        foreground,
        (row.right() - plot.right() - 8.0).max(50.0),
    );
    painter.galley(
        egui::pos2(plot.right() + 5.0, row.top()),
        galley,
        foreground,
    );
    if let Some(pointer) = response.hover_pos().filter(|p| plot.contains(*p)) {
        let pos = (start
            + (((pointer.x - plot.left()) / plot.width()) * (end - start) as f32).floor() as usize)
            .min(end - 1);
        let selection = trace_selection(view, lane, pos);
        let clicked = response.clicked();
        let details = selection
            .as_ref()
            .map(|s| s.details.clone())
            .unwrap_or_else(|| {
                format!(
                    "{}\nUnavailable window or no complete {}-bp motif footprint; NOT a zero score",
                    view.coordinate_label(pos),
                    trace.motif_length_bp
                )
            });
        response.on_hover_text(format!("{details}\n{}", lane.details));
        if clicked {
            return selection;
        }
    } else {
        response.on_hover_text(format!("{}\n{}\n{}", lane.units, lane.state, lane.details));
    }
    None
}

fn trace_selection(
    view: &TssSequenceView,
    lane: &TssViewLane,
    pos: usize,
) -> Option<TssViewFeature> {
    let trace = lane.trace.as_ref()?;
    let forward = trace.forward.get(pos).copied().flatten();
    let reverse = trace.reverse.get(pos).copied().flatten();
    if forward.is_none() && reverse.is_none() {
        return None;
    }
    let value = |score: Option<f64>| {
        score
            .map(|s| s.to_string())
            .unwrap_or_else(|| "unavailable".into())
    };
    Some(TssViewFeature {
        feature_id: None,
        start: pos,
        end: pos + trace.motif_length_bp,
        reverse: false,
        clipped: false,
        label: format!("{} | motif-window footprint", lane.label),
        details: format!(
            "{}\n{}\n{} + {} / - {} [{}]\nGenomic strand of local +: {}; of local -: {}. {}",
            view.coordinate_label(pos),
            view.coordinate_label(pos + trace.motif_length_bp - 1),
            if lane.kind == TssLaneKind::LocalScoreTrace {
                "Computed"
            } else {
                "Raw"
            },
            value(forward),
            value(reverse),
            lane.units,
            view.geometry.strand.as_str(),
            view.geometry.strand.opposite().as_str(),
            lane.state
        ),
        score: None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tss_svg_snapshot_matches_lane_filters_and_refuses_stale_or_loading_views() {
        let (dna, report) = crate::tss_sequence_view::profile_fixture(true);
        let mut area = MainAreaDna::new(dna.clone(), None, None);
        area.tss_view_available();
        let view = Arc::new(
            TssSequenceView::from_dna(&dna)
                .unwrap()
                .with_profile(&report)
                .unwrap(),
        );
        area.tss_ui.document = Some(Ok(view.clone()));
        area.tss_ui.traces = false;
        area.tss_ui.filter = "MA0001.1".into();
        let (_, options) = area.tss_svg_snapshot(ViewSvgExportProfile::Screen).unwrap();
        assert_eq!(options.lane_indices, area.tss_ui.visible_lanes(&view));
        assert!(!options.lane_indices.is_empty());
        assert!(
            options
                .lane_indices
                .iter()
                .all(|&i| view.lanes[i].kind == TssLaneKind::ImportedMotif)
        );
        assert!(
            crate::tss_sequence_view::render_tss_view_svg(&view, &options)
                .unwrap()
                .contains("data-lane-id=")
        );
        let (_, receiver) = std::sync::mpsc::channel();
        area.tss_ui.profile_load = Some(ProfileLoad {
            source: view,
            receiver: Arc::new(Mutex::new(receiver)),
        });
        assert!(
            area.tss_svg_snapshot(ViewSvgExportProfile::Screen)
                .unwrap_err()
                .contains("finish loading")
        );
        area.replace_loaded_sequence(crate::tss_sequence_view::tests::fixture(false));
        assert!(area.tss_svg_snapshot(ViewSvgExportProfile::Screen).is_err());
        assert!(area.tfbs_task.is_none());
    }

    #[test]
    fn tss_profile_native_frame_preserves_raw_scores_and_window_start_selection() {
        for minus in [false, true] {
            let (dna, report) = crate::tss_sequence_view::profile_fixture(minus);
            let base = Arc::new(TssSequenceView::from_dna(&dna).unwrap());
            let document = Arc::new(base.with_profile(&report).unwrap());
            let mut area = MainAreaDna::new(dna, None, None);
            area.tss_view_available();
            area.tss_ui.annotations = Some(base);
            area.tss_ui.document = Some(Ok(document.clone()));
            let curve = document
                .lanes
                .iter()
                .find(|l| l.kind == TssLaneKind::ScoreTrace)
                .unwrap();
            assert!(
                trace_selection(&document, curve, 1).is_none(),
                "unavailable is not zero"
            );
            assert!(
                trace_selection(&document, curve, 4).is_none(),
                "no full motif at the end"
            );
            let selection = trace_selection(&document, curve, 0).unwrap();
            assert_eq!(
                (selection.start, selection.end, selection.feature_id),
                (0, 3, None)
            );
            assert!(selection.details.contains("Raw + -2 / - 4 [llr_bits]"));
            let zero_curve = document
                .lanes
                .iter()
                .filter(|l| l.kind == TssLaneKind::ScoreTrace)
                .nth(2)
                .unwrap();
            assert!(
                trace_selection(&document, zero_curve, 4).is_some(),
                "real zero remains selectable"
            );
            let ctx = egui::Context::default();
            for _ in 0..2 {
                ctx.begin_pass(egui::RawInput {
                    screen_rect: Some(egui::Rect::from_min_size(
                        egui::Pos2::ZERO,
                        egui::vec2(1100.0, 1400.0),
                    )),
                    ..Default::default()
                });
                crate::egui_compat::show_central_panel_for_test_context(
                    &ctx,
                    egui::CentralPanel::default(),
                    |ui| area.render_primary_tss_map_ui(ui),
                );
                assert!(!crate::egui_compat::end_test_pass(&ctx).shapes.is_empty());
                assert!(Arc::ptr_eq(
                    area.tss_ui.document.as_ref().unwrap().as_ref().unwrap(),
                    &document
                ));
                assert!(area.tss_ui.profile_load.is_none());
                assert!(
                    area.tfbs_task.is_none(),
                    "painting must not start scoring or DuckDB"
                );
            }
        }
    }

    #[test]
    fn tss_profile_background_loading_does_not_use_panel_cache() {
        let (dna, report) = crate::tss_sequence_view::profile_fixture(false);
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("report.json");
        std::fs::write(&path, serde_json::to_vec(&report).unwrap()).unwrap();
        let mut area = MainAreaDna::new(dna.clone(), None, None);
        area.tss_view_available();
        area.tss_ui.document = Some(Ok(Arc::new(TssSequenceView::from_dna(&dna).unwrap())));
        let ctx = egui::Context::default();
        area.load_tss_profile(path, &ctx);
        let deadline = Instant::now() + std::time::Duration::from_secs(5);
        while area.tss_ui.profile_load.is_some() && Instant::now() < deadline {
            area.poll_tss_profile(&ctx);
            std::thread::sleep(std::time::Duration::from_millis(1));
        }
        assert!(area.tss_ui.profile_load.is_none());
        assert!(
            area.tss_ui.profile_error.is_none(),
            "{:?}",
            area.tss_ui.profile_error
        );
        assert!(
            area.tss_ui
                .document
                .as_ref()
                .unwrap()
                .as_ref()
                .unwrap()
                .profile
                .is_some()
        );
        assert!(area.cached_tfbs_score_tracks.is_none());
        assert!(area.cached_genomic_motif_evidence.is_none());
        assert_eq!(
            area.dna.read().unwrap().get_forward_string(),
            dna.get_forward_string()
        );
    }

    #[test]
    fn shared_profile_intent_uses_the_native_attachment_path() {
        let (dna, report) = crate::tss_sequence_view::profile_fixture(false);
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("report.json");
        std::fs::write(&path, serde_json::to_vec(&report).unwrap()).unwrap();
        let mut area = MainAreaDna::new(dna.clone(), None, None);
        area.tss_view_available();
        area.tss_ui.document = Some(Ok(Arc::new(TssSequenceView::from_dna(&dna).unwrap())));
        area.queue_tss_profile(path).unwrap();
        assert_eq!(area.primary_map_mode, PrimaryMapMode::Tss);

        let ctx = egui::Context::default();
        let deadline = Instant::now() + std::time::Duration::from_secs(5);
        while area
            .tss_ui
            .document
            .as_ref()
            .and_then(|document| document.as_ref().ok())
            .is_none_or(|document| document.profile.is_none())
            && Instant::now() < deadline
        {
            ctx.begin_pass(egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(
                    egui::Pos2::ZERO,
                    egui::vec2(1000.0, 900.0),
                )),
                ..Default::default()
            });
            crate::egui_compat::show_central_panel_for_test_context(
                &ctx,
                egui::CentralPanel::default(),
                |ui| area.render_primary_tss_map_ui(ui),
            );
            let _ = crate::egui_compat::end_test_pass(&ctx);
            std::thread::sleep(std::time::Duration::from_millis(1));
        }
        assert!(area.tss_ui.requested_profile_path.is_none());
        assert!(area.tss_ui.profile_load.is_none());
        assert!(area.tss_ui.profile_error.is_none());
        assert!(
            area.tss_ui
                .document
                .as_ref()
                .unwrap()
                .as_ref()
                .unwrap()
                .profile
                .is_some()
        );
        assert!(area.cached_tfbs_score_tracks.is_none());
        assert!(area.cached_genomic_motif_evidence.is_none());
    }

    #[test]
    fn tss_profile_worker_failure_or_stale_result_keeps_current_document() {
        let (dna, report) = crate::tss_sequence_view::profile_fixture(false);
        let base = Arc::new(TssSequenceView::from_dna(&dna).unwrap());
        let prepared = Arc::new(base.with_profile(&report).unwrap());
        let mut area = MainAreaDna::new(dna, None, None);
        let ctx = egui::Context::default();
        for fail in [false, true] {
            let (send, receiver) = std::sync::mpsc::channel();
            area.tss_ui.document = Some(Ok(prepared.clone()));
            area.tss_ui.profile_load = Some(ProfileLoad {
                source: if fail { prepared.clone() } else { base.clone() },
                receiver: Arc::new(Mutex::new(receiver)),
            });
            send.send(if fail {
                Err("bad report".into())
            } else {
                Ok(base.clone())
            })
            .unwrap();
            area.poll_tss_profile(&ctx);
            assert!(Arc::ptr_eq(
                area.tss_ui.document.as_ref().unwrap().as_ref().unwrap(),
                &prepared
            ));
        }
        assert_eq!(area.tss_ui.profile_error.as_deref(), Some("bad report"));
        area.replace_loaded_sequence(crate::tss_sequence_view::tests::fixture(true));
        assert!(area.tss_ui.document.is_none());
        assert!(area.tss_ui.profile_load.is_none());
        assert!(area.tss_ui.annotations.is_none());
    }

    #[test]
    fn tss_view_switches_without_mutating_sequence_and_invalidates_on_replacement() {
        let dna = crate::tss_sequence_view::tests::fixture(false);
        let mut area = MainAreaDna::new(dna.clone(), None, None);
        area.set_tss_view(true).unwrap();
        assert_eq!(area.primary_map_mode, PrimaryMapMode::Tss);
        assert_eq!(
            area.dna.read().unwrap().get_forward_string(),
            dna.get_forward_string()
        );
        area.set_tss_view(false).unwrap();
        assert_eq!(area.primary_map_mode, PrimaryMapMode::Standard);
        area.primary_map_mode = PrimaryMapMode::Splicing;
        assert!(area.set_tss_view(false).is_err());
        assert_eq!(area.primary_map_mode, PrimaryMapMode::Splicing);
        area.replace_loaded_sequence(DNAsequence::from_sequence("ACGT").unwrap());
        assert!(area.set_tss_view(true).is_err());
    }

    #[test]
    fn tss_view_native_frame_keeps_cached_document_and_renders_both_strands() {
        for minus in [false, true] {
            let dna = crate::tss_sequence_view::tests::fixture(minus);
            let mut area = MainAreaDna::new(dna.clone(), None, None);
            area.tss_view_available();
            let document = Arc::new(TssSequenceView::from_dna(&dna).unwrap());
            area.tss_ui.document = Some(Ok(document.clone()));
            let ctx = egui::Context::default();
            for _ in 0..2 {
                ctx.begin_pass(egui::RawInput {
                    screen_rect: Some(egui::Rect::from_min_size(
                        egui::Pos2::ZERO,
                        egui::vec2(900.0, 700.0),
                    )),
                    ..Default::default()
                });
                crate::egui_compat::show_central_panel_for_test_context(
                    &ctx,
                    egui::CentralPanel::default(),
                    |ui| area.render_primary_tss_map_ui(ui),
                );
                let output = crate::egui_compat::end_test_pass(&ctx);
                assert!(!output.shapes.is_empty());
                assert!(Arc::ptr_eq(
                    area.tss_ui.document.as_ref().unwrap().as_ref().unwrap(),
                    &document
                ));
                assert!(area.tss_ui.pending.is_none());
            }
        }
    }

    #[test]
    fn tss_view_establishes_a_vertical_document_layout_inside_horizontal_host() {
        let dna = crate::tss_sequence_view::tests::fixture(false);
        let mut area = MainAreaDna::new(dna.clone(), None, None);
        area.tss_view_available();
        area.tss_ui.document = Some(Ok(Arc::new(TssSequenceView::from_dna(&dna).unwrap())));
        let ctx = egui::Context::default();
        let mut used = egui::Rect::NOTHING;
        ctx.begin_pass(egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(
                egui::Pos2::ZERO,
                egui::vec2(900.0, 700.0),
            )),
            ..Default::default()
        });
        crate::egui_compat::show_central_panel_for_test_context(
            &ctx,
            egui::CentralPanel::default(),
            |ui| {
                ui.horizontal(|ui| {
                    area.render_primary_tss_map_ui(ui);
                    used = ui.min_rect();
                });
            },
        );
        let _ = crate::egui_compat::end_test_pass(&ctx);
        assert!(used.height() > 250.0, "TSS document collapsed to {used:?}");
        assert!(
            used.width() <= 900.0,
            "TSS document escaped viewport: {used:?}"
        );
    }
}
