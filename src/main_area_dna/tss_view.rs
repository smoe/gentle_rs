//! Native TSS lanes over the shared annotated-window presentation model.

use super::*;
use crate::tss_sequence_view::{TssLaneKind, TssSequenceView, TssViewLane};
use std::sync::mpsc::{Receiver, TryRecvError};

type LoadedView = Result<Arc<TssSequenceView>, String>;

#[derive(Clone, Debug)]
pub(super) struct TssUiState {
    key: Option<(u64, usize)>,
    detected: bool,
    pending: Option<Arc<std::sync::Mutex<Receiver<LoadedView>>>>,
    document: Option<LoadedView>,
    structures: bool,
    signals: bool,
    motifs: bool,
    other: bool,
    filter: String,
    selected: Option<(usize, usize)>,
}

impl Default for TssUiState {
    fn default() -> Self {
        Self {
            key: None,
            detected: false,
            pending: None,
            document: None,
            structures: true,
            signals: true,
            motifs: true,
            other: true,
            filter: String::new(),
            selected: None,
        }
    }
}

fn color(kind: TssLaneKind) -> egui::Color32 {
    match kind {
        TssLaneKind::Structure => egui::Color32::from_rgb(55, 130, 170),
        TssLaneKind::Signal => egui::Color32::from_rgb(165, 65, 100),
        TssLaneKind::Motif => egui::Color32::from_rgb(0, 140, 115),
        TssLaneKind::Other => egui::Color32::from_rgb(180, 125, 30),
    }
}

impl MainAreaDna {
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

    pub(super) fn render_primary_tss_map_ui(&mut self, ui: &mut egui::Ui) {
        self.poll_tss_view(ui.ctx());
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
            ui.checkbox(&mut self.tss_ui.structures, "Exons / CDS");
            ui.checkbox(&mut self.tss_ui.signals, "CUT&RUN / chromatin");
            ui.checkbox(&mut self.tss_ui.motifs, "Stored motif peaks");
            ui.checkbox(&mut self.tss_ui.other, "Other annotations");
            ui.label("Filter lanes");
            ui.add(egui::TextEdit::singleline(&mut self.tss_ui.filter).desired_width(150.0));
        });
        ui.small("Blue: structure; rose: supplied signal; green triangles: predicted motif (+ above / - below). No new scoring. Motif values below 0 are hidden, not deleted. Click a feature to select its exact DNA span.");
        ui.collapsing("Provenance, limitations and missing data", |ui| {
            for warning in &view.warnings {
                ui.label(warning);
            }
            ui.label(&view.provenance);
        });

        let filter = self.tss_ui.filter.to_lowercase();
        let lanes: Vec<usize> = view
            .lanes
            .iter()
            .enumerate()
            .filter(|(_, lane)| {
                let enabled = match lane.kind {
                    TssLaneKind::Structure => self.tss_ui.structures,
                    TssLaneKind::Signal => self.tss_ui.signals,
                    TssLaneKind::Motif => self.tss_ui.motifs,
                    TssLaneKind::Other => self.tss_ui.other,
                };
                enabled
                    && (filter.is_empty()
                        || lane.label.to_lowercase().contains(&filter)
                        || lane.id.to_lowercase().contains(&filter))
            })
            .map(|(i, _)| i)
            .collect();
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
            .show_rows(ui, 90.0, lanes.len(), |ui, rows| {
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
                        selected = Some((index, feature));
                    }
                }
            });
        if let Some((lane, feature)) = selected {
            self.tss_ui.selected = Some((lane, feature));
            let f = &view.lanes[lane].features[feature];
            if let Err(error) = self.set_selection_range_0based(f.start, f.end) {
                self.op_status = error;
            }
        }
        if let Some((lane, feature)) = self.tss_ui.selected {
            let f = &view.lanes[lane].features[feature];
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
                    self.map_dna.select_feature(Some(f.feature_id));
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
) -> Option<usize> {
    let (rect, response) = ui.allocate_exact_size(egui::vec2(width, 86.0), egui::Sense::click());
    let plot = egui::Rect::from_min_max(
        egui::pos2(rect.left() + left, rect.top() + 9.0),
        egui::pos2(
            (rect.right() - right).max(rect.left() + left + 50.0),
            rect.bottom() - 12.0,
        ),
    );
    let painter = ui.painter();
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
    let middle = lane.kind == TssLaneKind::Motif
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
    if lane.kind == TssLaneKind::Signal || lane.kind == TssLaneKind::Motif {
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
            "0",
            egui::FontId::monospace(9.0),
            foreground,
        );
        if middle {
            let bottom_max = if lane.kind == TssLaneKind::Motif {
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
        if lane.kind == TssLaneKind::Motif {
            let height = (f.score.unwrap_or(0.0) / lane.scale_max) as f32 * amplitude;
            let y = baseline + if f.reverse { height } else { -height };
            let pos = egui::pos2(x(f.start), y);
            let offset = if f.reverse { -4.0 } else { 4.0 };
            p.line_segment(
                [egui::pos2(pos.x, baseline), pos],
                egui::Stroke::new(0.6, color(lane.kind).gamma_multiply(0.5)),
            );
            p.add(egui::Shape::convex_polygon(
                vec![
                    pos,
                    pos + egui::vec2(-3.0, offset),
                    pos + egui::vec2(3.0, offset),
                ],
                color(lane.kind),
                egui::Stroke::NONE,
            ));
            feature_rect = egui::Rect::from_center_size(pos, egui::vec2(8.0, 10.0));
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
            "{count}\n{units}\ngenomic {}..{}\nlocal {}..{}",
            view.geometry.genomic_at(first).unwrap(),
            view.geometry.genomic_at(last).unwrap(),
            first + 1,
            last + 1
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
            return Some(index);
        }
    } else {
        response.on_hover_text(format!("{}\n{}\n{}", lane.id, lane.units, lane.details));
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

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
}
