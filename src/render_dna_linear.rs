use crate::{
    dna_display::{DnaDisplay, Selection, TfbsDisplayCriteria},
    dna_sequence::DNAsequence,
    open_reading_frame::OpenReadingFrame,
    render_dna::RenderDna,
    render_dna::RestrictionEnzymePosition,
};
use eframe::egui::{
    self, Align2, Color32, FontFamily, FontId, PointerState, Pos2, Rect, Stroke, Vec2,
};
use gb_io::seq::{Feature, Location};
use std::collections::{BTreeSet, HashMap};
use std::sync::{Arc, RwLock};

const BASELINE_STROKE: f32 = 2.0;
const FEATURE_HEIGHT: f32 = 12.0;
const FEATURE_GAP: f32 = 18.0;
const BASELINE_MARGIN: f32 = 30.0;
const LABEL_ROW_HEIGHT: f32 = 12.0;
const LABEL_CHAR_WIDTH: f32 = 6.5;
const ORF_HEIGHT: f32 = 6.0;
const GC_STRIP_HEIGHT: f32 = 6.0;
const METHYLATION_TICK: f32 = 10.0;
const RE_LABEL_BASE_OFFSET: f32 = 76.0;
const RE_SITE_MAX_BP_PER_PX: f32 = 200.0;
const RE_LABEL_MAX_BP_PER_PX: f32 = 30.0;
const FEATURE_LABEL_MAX_BP_PER_PX: f32 = 120.0;
const METHYLATION_MAX_BP_PER_PX: f32 = 40.0;

#[derive(Debug, Clone, Copy)]
struct LinearViewport {
    start: usize,
    end: usize,
    span: usize,
}

#[derive(Debug, Clone, Copy)]
struct LinearDetailLevel {
    show_feature_labels: bool,
    show_restriction_sites: bool,
    show_restriction_labels: bool,
    show_methylation_sites: bool,
}

#[derive(Debug, Clone)]
struct FeaturePosition {
    feature_number: usize,
    from: usize,
    to: usize,
    label: String,
    color: Color32,
    is_pointy: bool,
    is_reverse: bool,
    rect: Rect,
    label_pos: Pos2,
}

impl FeaturePosition {
    fn contains(&self, pos: Pos2) -> bool {
        self.rect.contains(pos)
    }
}

#[derive(Debug, Clone)]
pub struct RenderDnaLinear {
    dna: Arc<RwLock<DNAsequence>>,
    display: Arc<RwLock<DnaDisplay>>,
    sequence_length: usize,
    area: Rect,
    features: Vec<FeaturePosition>,
    restriction_enzyme_sites: Vec<RestrictionEnzymePosition>,
    selected_feature_number: Option<usize>,
    hovered_feature_number: Option<usize>,
    hover_enzyme: Option<RestrictionEnzymePosition>,
}

impl RenderDnaLinear {
    pub fn new(dna: Arc<RwLock<DNAsequence>>, display: Arc<RwLock<DnaDisplay>>) -> Self {
        let sequence_length = dna.read().map(|d| d.len()).unwrap_or(0);
        Self {
            dna,
            display,
            sequence_length,
            area: Rect::NOTHING,
            features: vec![],
            restriction_enzyme_sites: vec![],
            selected_feature_number: None,
            hovered_feature_number: None,
            hover_enzyme: None,
        }
    }

    pub fn area(&self) -> &Rect {
        &self.area
    }

    fn is_rect_usable(rect: Rect) -> bool {
        rect.min.x.is_finite()
            && rect.min.y.is_finite()
            && rect.max.x.is_finite()
            && rect.max.y.is_finite()
            && rect.width() > 0.0
            && rect.height() > 0.0
    }

    fn layout_needs_recomputing(&self) -> bool {
        self.display
            .read()
            .map(|d| d.update_layout().update_map_dna())
            .unwrap_or(false)
    }

    fn layout_was_updated(&self) {
        if let Ok(mut display) = self.display.write() {
            display.update_layout_mut().map_dna_updated();
        }
    }

    fn baseline_y(&self) -> f32 {
        self.area.center().y
    }

    fn viewport(&self) -> LinearViewport {
        if self.sequence_length == 0 {
            return LinearViewport {
                start: 0,
                end: 0,
                span: 0,
            };
        }
        let (start_bp, span_bp) = self
            .display
            .read()
            .map(|display| {
                (
                    display.linear_view_start_bp(),
                    display.linear_view_span_bp(),
                )
            })
            .unwrap_or((0, self.sequence_length));
        let span = if span_bp == 0 || span_bp > self.sequence_length {
            self.sequence_length
        } else {
            span_bp
        };
        let max_start = self.sequence_length.saturating_sub(span);
        let start = start_bp.min(max_start);
        let end = start.saturating_add(span).min(self.sequence_length);
        LinearViewport { start, end, span }
    }

    fn bp_to_x(&self, bp: usize, viewport: LinearViewport) -> f32 {
        if viewport.span == 0 {
            return self.area.left();
        }
        let frac = (bp.saturating_sub(viewport.start)) as f32 / viewport.span as f32;
        self.area.left() + frac * self.area.width()
    }

    fn x_to_bp(&self, x: f32, viewport: LinearViewport) -> usize {
        if viewport.span == 0 {
            return 0;
        }
        let frac = ((x - self.area.left()) / self.area.width().max(1.0)).clamp(0.0, 1.0);
        let bp = viewport.start + (frac * viewport.span as f32).floor() as usize;
        bp.min(viewport.end.saturating_sub(1))
    }

    fn normalize_pos(&self, pos: isize) -> usize {
        if self.sequence_length == 0 {
            return 0;
        }
        let len = self.sequence_length as isize;
        (((pos % len) + len) % len) as usize
    }

    fn range_overlap(
        a_start: usize,
        a_end_exclusive: usize,
        b_start: usize,
        b_end_exclusive: usize,
    ) -> Option<(usize, usize)> {
        let start = a_start.max(b_start);
        let end = a_end_exclusive.min(b_end_exclusive);
        (start < end).then_some((start, end))
    }

    fn bp_per_px(&self, viewport: LinearViewport) -> f32 {
        if viewport.span == 0 {
            return 0.0;
        }
        viewport.span as f32 / self.area.width().max(1.0)
    }

    fn detail_level(&self, viewport: LinearViewport) -> LinearDetailLevel {
        let bp_per_px = self.bp_per_px(viewport);
        LinearDetailLevel {
            show_feature_labels: bp_per_px <= FEATURE_LABEL_MAX_BP_PER_PX,
            show_restriction_sites: bp_per_px <= RE_SITE_MAX_BP_PER_PX,
            show_restriction_labels: bp_per_px <= RE_LABEL_MAX_BP_PER_PX,
            show_methylation_sites: bp_per_px <= METHYLATION_MAX_BP_PER_PX,
        }
    }

    fn draw_feature(
        feature: &Feature,
        show_cds_features: bool,
        show_gene_features: bool,
        show_mrna_features: bool,
        show_tfbs: bool,
        tfbs_display_criteria: TfbsDisplayCriteria,
        hidden_feature_kinds: &BTreeSet<String>,
    ) -> bool {
        if RenderDna::is_source_feature(feature) {
            return false;
        }
        let feature_kind = feature.kind.to_string().to_ascii_uppercase();
        if hidden_feature_kinds.contains(&feature_kind) {
            return false;
        }
        if !RenderDna::feature_passes_kind_filter(
            feature,
            show_cds_features,
            show_gene_features,
            show_mrna_features,
        ) {
            return false;
        }
        if RenderDna::is_tfbs_feature(feature) {
            if !show_tfbs {
                return false;
            }
            return RenderDna::tfbs_feature_passes_display_filter(feature, tfbs_display_criteria);
        }
        true
    }

    fn estimate_label_width(label: &str) -> f32 {
        let chars = label.chars().count().max(1) as f32;
        chars * LABEL_CHAR_WIDTH
    }

    fn allocate_lane(lane_ends: &mut Vec<f32>, start: f32, end: f32, padding: f32) -> usize {
        for (i, lane_end) in lane_ends.iter_mut().enumerate() {
            if start >= *lane_end + padding {
                *lane_end = end;
                return i;
            }
        }
        lane_ends.push(end);
        lane_ends.len() - 1
    }

    fn tick_step(span: usize) -> usize {
        if span <= 1 {
            return 1;
        }
        let raw = (span as f64 / 8.0).max(1.0);
        let base = 10_f64.powf(raw.log10().floor());
        for scale in [1.0, 2.0, 5.0, 10.0] {
            let step = (base * scale).round() as usize;
            if (step as f64) >= raw {
                return step.max(1);
            }
        }
        span.max(1)
    }

    fn collect_location_strands(location: &Location, reverse: bool, strands: &mut Vec<bool>) {
        match location {
            Location::Range(_, _) | Location::Between(_, _) => strands.push(reverse),
            Location::Complement(inner) => Self::collect_location_strands(inner, !reverse, strands),
            Location::Join(parts)
            | Location::Order(parts)
            | Location::Bond(parts)
            | Location::OneOf(parts) => {
                for part in parts {
                    Self::collect_location_strands(part, reverse, strands);
                }
            }
            Location::External(_, maybe_loc) => {
                if let Some(loc) = maybe_loc {
                    Self::collect_location_strands(loc, reverse, strands);
                }
            }
            Location::Gap(_) => {}
        }
    }

    fn feature_is_reverse(feature: &Feature) -> bool {
        let mut strands = Vec::new();
        Self::collect_location_strands(&feature.location, false, &mut strands);
        if strands.is_empty() {
            false
        } else {
            strands.iter().filter(|is_reverse| **is_reverse).count() > strands.len() / 2
        }
    }

    fn layout_features(&mut self, viewport: LinearViewport) {
        self.features.clear();
        if self.sequence_length == 0 {
            return;
        }

        #[derive(Clone)]
        struct Seed {
            feature_number: usize,
            from: usize,
            to: usize,
            x1: f32,
            x2: f32,
            label: String,
            color: Color32,
            is_pointy: bool,
            is_reverse: bool,
        }
        impl Seed {
            fn span(&self) -> usize {
                self.to.saturating_sub(self.from).saturating_add(1)
            }
        }

        #[derive(Clone)]
        struct PositionedSeed {
            seed: Seed,
            feature_lane: usize,
        }

        let mut seeds: Vec<Seed> = Vec::new();
        let (
            show_cds_features,
            show_gene_features,
            show_mrna_features,
            show_tfbs,
            tfbs_display_criteria,
            hidden_feature_kinds,
        ) = self
            .display
            .read()
            .map(|display| {
                (
                    display.show_cds_features(),
                    display.show_gene_features(),
                    display.show_mrna_features(),
                    display.show_tfbs(),
                    display.tfbs_display_criteria(),
                    display.hidden_feature_kinds().clone(),
                )
            })
            .unwrap_or((
                true,
                true,
                true,
                false,
                TfbsDisplayCriteria::default(),
                BTreeSet::new(),
            ));
        let features = self
            .dna
            .read()
            .map(|dna| dna.features().to_owned())
            .unwrap_or_default();
        for (feature_number, feature) in features.iter().enumerate() {
            if !Self::draw_feature(
                feature,
                show_cds_features,
                show_gene_features,
                show_mrna_features,
                show_tfbs,
                tfbs_display_criteria,
                &hidden_feature_kinds,
            ) {
                continue;
            }

            let (raw_from, raw_to) = match feature.location.find_bounds() {
                Ok(bounds) => bounds,
                Err(_) => continue,
            };
            if raw_from < 0 || raw_to < 0 {
                continue;
            }

            let mut from = raw_from as usize;
            let mut to = raw_to as usize;
            if to < from {
                std::mem::swap(&mut to, &mut from);
            }
            if from >= self.sequence_length {
                continue;
            }
            if to >= self.sequence_length {
                to = self.sequence_length.saturating_sub(1);
            }

            let feature_end_exclusive = to.saturating_add(1).min(self.sequence_length);
            let Some((visible_start, visible_end)) =
                Self::range_overlap(from, feature_end_exclusive, viewport.start, viewport.end)
            else {
                continue;
            };
            let x1 = self.bp_to_x(visible_start, viewport).max(self.area.left());
            let x2 = self
                .bp_to_x(visible_end, viewport)
                .max(x1 + 1.0)
                .min(self.area.right());
            let label = RenderDna::feature_name(feature);

            seeds.push(Seed {
                feature_number,
                from,
                to,
                x1,
                x2,
                label,
                color: RenderDna::feature_color(feature),
                is_pointy: RenderDna::is_feature_pointy(feature),
                is_reverse: Self::feature_is_reverse(feature),
            });
        }

        let mut feature_lanes_top: Vec<f32> = vec![];
        let mut feature_lanes_bottom: Vec<f32> = vec![];
        let mut lane_seed: Vec<PositionedSeed> = Vec::with_capacity(seeds.len());
        let mut lane_order = seeds.clone();
        lane_order.sort_by(|a, b| {
            b.span()
                .cmp(&a.span())
                .then_with(|| a.x1.total_cmp(&b.x1))
                .then_with(|| a.feature_number.cmp(&b.feature_number))
        });

        for seed in lane_order {
            let feature_lane = if seed.is_reverse {
                Self::allocate_lane(&mut feature_lanes_bottom, seed.x1, seed.x2, 3.0)
            } else {
                Self::allocate_lane(&mut feature_lanes_top, seed.x1, seed.x2, 3.0)
            };
            lane_seed.push(PositionedSeed { seed, feature_lane });
        }

        lane_seed.sort_by(|a, b| {
            a.seed
                .x1
                .total_cmp(&b.seed.x1)
                .then_with(|| a.seed.feature_number.cmp(&b.seed.feature_number))
        });

        for item in lane_seed {
            let seed = item.seed;
            let feature_lane = item.feature_lane;
            let center_y = if seed.is_reverse {
                self.baseline_y() + BASELINE_MARGIN + FEATURE_GAP * feature_lane as f32
            } else {
                self.baseline_y() - BASELINE_MARGIN - FEATURE_GAP * feature_lane as f32
            };
            let rect = Rect::from_min_max(
                Pos2::new(seed.x1, center_y - FEATURE_HEIGHT / 2.0),
                Pos2::new(seed.x2, center_y + FEATURE_HEIGHT / 2.0),
            );

            let label_pos = if seed.is_reverse {
                Pos2::new(seed.x1, rect.bottom() + 2.0)
            } else {
                Pos2::new(seed.x1, rect.top() - 2.0)
            };

            self.features.push(FeaturePosition {
                feature_number: seed.feature_number,
                from: seed.from,
                to: seed.to,
                label: seed.label,
                color: seed.color,
                is_pointy: seed.is_pointy,
                is_reverse: seed.is_reverse,
                rect,
                label_pos,
            });
        }
    }

    fn get_re_site_for_positon(&self, pos: Pos2) -> Option<RestrictionEnzymePosition> {
        self.restriction_enzyme_sites
            .iter()
            .find(|site| site.area.contains(pos))
            .cloned()
    }

    fn get_clicked_feature(&self, pos: Pos2) -> Option<&FeaturePosition> {
        self.features.iter().find(|feature| feature.contains(pos))
    }

    pub fn selected_feature_number(&self) -> Option<usize> {
        self.selected_feature_number
    }

    pub fn select_feature(&mut self, feature_number: Option<usize>) {
        self.selected_feature_number = feature_number;
    }

    pub fn coordinate_to_basepair(&self, x: f32) -> Result<i64, String> {
        if !self.area.contains(Pos2::new(x, self.baseline_y())) {
            return Err("Coordinate is outside the DNA visualization area.".to_string());
        }
        Ok(self.x_to_bp(x, self.viewport()) as i64)
    }

    pub fn on_hover(&mut self, pointer_state: PointerState) {
        if let Some(pos) = pointer_state.latest_pos() {
            self.hovered_feature_number = self.get_clicked_feature(pos).map(|f| f.feature_number);
            self.hover_enzyme = self.get_re_site_for_positon(pos);
        }
    }

    pub fn on_click(&mut self, pointer_state: PointerState) {
        if let Some(pos) = pointer_state.latest_pos() {
            self.selected_feature_number = self.get_clicked_feature(pos).map(|f| f.feature_number);
        }
    }

    pub fn on_double_click(&mut self, pointer_state: PointerState) {
        if let Ok(mut display) = self.display.write() {
            display.deselect();
        }

        if let Some(pos) = pointer_state.latest_pos() {
            if let Some(feature) = self.get_clicked_feature(pos) {
                let selection = Selection::new(feature.from, feature.to, self.sequence_length);
                if let Ok(mut display) = self.display.write() {
                    display.select(selection);
                }
            } else if let Some(re_pos) = self.get_re_site_for_positon(pos) {
                let selection = Selection::new(
                    re_pos.key().from() as usize,
                    re_pos.key().to() as usize,
                    self.sequence_length,
                );
                if let Ok(mut display) = self.display.write() {
                    display.select(selection);
                }
            }
        }
    }

    fn draw_backbone(&self, painter: &egui::Painter) {
        let y = self.baseline_y();
        painter.line_segment(
            [
                Pos2::new(self.area.left(), y),
                Pos2::new(self.area.right(), y),
            ],
            Stroke::new(BASELINE_STROKE, Color32::BLACK),
        );
    }

    fn draw_bp_ticks(&self, painter: &egui::Painter, viewport: LinearViewport) {
        if viewport.span == 0 {
            return;
        }
        let y = self.baseline_y();
        let tick = Self::tick_step(viewport.span);
        let font = FontId {
            size: 9.0,
            family: FontFamily::Monospace,
        };

        let mut pos = if viewport.start == 0 {
            0
        } else {
            ((viewport.start / tick) + 1) * tick
        };
        while pos < viewport.end {
            let x = self.bp_to_x(pos, viewport);
            painter.line_segment(
                [Pos2::new(x, y - 4.0), Pos2::new(x, y + 4.0)],
                Stroke::new(1.0, Color32::DARK_GRAY),
            );
            painter.text(
                Pos2::new(x, y + 7.0),
                Align2::CENTER_TOP,
                pos.to_string(),
                font.clone(),
                Color32::DARK_GRAY,
            );
            pos += tick;
        }
    }

    fn draw_name_and_length(&self, painter: &egui::Painter, viewport: LinearViewport) {
        let name = self
            .dna
            .read()
            .ok()
            .and_then(|dna| dna.name().clone())
            .unwrap_or_else(|| "<no name>".to_string());
        painter.text(
            Pos2::new(self.area.left() + 6.0, self.area.top() + 6.0),
            Align2::LEFT_TOP,
            name,
            FontId {
                size: 13.0,
                family: FontFamily::Proportional,
            },
            Color32::BLACK,
        );
        painter.text(
            Pos2::new(self.area.right() - 6.0, self.area.top() + 6.0),
            Align2::RIGHT_TOP,
            if viewport.span > 0 && viewport.span < self.sequence_length {
                format!(
                    "view {}..{} ({} bp) of {} bp",
                    viewport.start.saturating_add(1),
                    viewport.end,
                    viewport.span,
                    self.sequence_length
                )
            } else {
                format!("{} bp", self.sequence_length)
            },
            FontId {
                size: 11.0,
                family: FontFamily::Monospace,
            },
            Color32::DARK_GRAY,
        );
    }

    fn orf_colors() -> HashMap<i32, Color32> {
        let mut colors = HashMap::new();
        colors.insert(-1, Color32::LIGHT_RED);
        colors.insert(-2, Color32::LIGHT_GREEN);
        colors.insert(-3, Color32::LIGHT_BLUE);
        colors.insert(1, Color32::DARK_RED);
        colors.insert(2, Color32::DARK_GREEN);
        colors.insert(3, Color32::DARK_BLUE);
        colors
    }

    fn draw_orf(
        &self,
        painter: &egui::Painter,
        orf: &OpenReadingFrame,
        color: Color32,
        viewport: LinearViewport,
    ) {
        let from = self.normalize_pos(orf.from() as isize);
        let to = self.normalize_pos(orf.to() as isize);
        let start = from.min(to);
        let end_exclusive = from.max(to).saturating_add(1).min(self.sequence_length);
        if start >= end_exclusive {
            return;
        }

        let Some((draw_start, draw_end)) =
            Self::range_overlap(start, end_exclusive, viewport.start, viewport.end)
        else {
            return;
        };

        let x1 = self.bp_to_x(draw_start, viewport);
        let x2 = self.bp_to_x(draw_end, viewport).max(x1 + 1.0);
        let frame_abs = orf.frame().unsigned_abs() as f32;
        let y = if orf.is_reverse() {
            self.baseline_y() + 10.0 + frame_abs * 7.0
        } else {
            self.baseline_y() - 10.0 - frame_abs * 7.0
        };

        let rect = Rect::from_min_max(
            Pos2::new(x1, y - ORF_HEIGHT / 2.0),
            Pos2::new(x2, y + ORF_HEIGHT / 2.0),
        );
        painter.rect_filled(rect, 1.0, color);

        let tip = if orf.is_reverse() {
            Pos2::new(x1 - 5.0, y)
        } else {
            Pos2::new(x2 + 5.0, y)
        };
        let base_x = if orf.is_reverse() { x1 } else { x2 };
        painter.add(egui::Shape::convex_polygon(
            vec![
                Pos2::new(base_x, y - ORF_HEIGHT / 2.0),
                tip,
                Pos2::new(base_x, y + ORF_HEIGHT / 2.0),
            ],
            color,
            Stroke::NONE,
        ));

        let label = format!("ORF {}", orf.frame());
        let label_width = Self::estimate_label_width(&label);
        if (x2 - x1) >= label_width + 6.0 {
            let label_pos = Pos2::new((x1 + x2) / 2.0, y);
            painter.text(
                label_pos,
                Align2::CENTER_CENTER,
                label,
                FontId {
                    size: 8.0,
                    family: FontFamily::Monospace,
                },
                Color32::BLACK,
            );
        }
    }

    fn draw_open_reading_frames(&self, painter: &egui::Painter, viewport: LinearViewport) {
        let show_orfs = self
            .display
            .read()
            .map(|display| display.show_open_reading_frames())
            .unwrap_or(false);
        if !show_orfs {
            return;
        }
        let colors = Self::orf_colors();
        let orfs = self
            .dna
            .read()
            .map(|dna| dna.open_reading_frames().clone())
            .unwrap_or_default();
        for orf in &orfs {
            if let Some(color) = colors.get(&orf.frame()) {
                self.draw_orf(painter, orf, *color, viewport);
            }
        }
    }

    fn draw_gc_contents(&self, painter: &egui::Painter, viewport: LinearViewport) {
        let show_gc = self
            .display
            .read()
            .map(|display| display.show_gc_contents())
            .unwrap_or(false);
        if !show_gc {
            return;
        }
        let y1 = self.baseline_y() - 4.0;
        let y2 = y1 + GC_STRIP_HEIGHT;
        let gc_contents = self
            .dna
            .read()
            .map(|dna| dna.gc_content().clone())
            .unwrap_or_default();
        for region in gc_contents.regions() {
            let region_end_exclusive = region.to().saturating_add(1).min(self.sequence_length);
            let Some((draw_start, draw_end)) = Self::range_overlap(
                region.from(),
                region_end_exclusive,
                viewport.start,
                viewport.end,
            ) else {
                continue;
            };
            let x1 = self.bp_to_x(draw_start, viewport);
            let x2 = self.bp_to_x(draw_end, viewport).max(x1 + 1.0);
            let color = Color32::from_rgb(
                255 - (region.gc() * 255.0) as u8,
                (region.gc() * 255.0) as u8,
                0,
            );
            painter.rect_filled(
                Rect::from_min_max(Pos2::new(x1, y1), Pos2::new(x2, y2)),
                0.0,
                color,
            );
        }
    }

    fn draw_methylation_sites(
        &self,
        painter: &egui::Painter,
        viewport: LinearViewport,
        detail: LinearDetailLevel,
    ) {
        let show_methylation = self
            .display
            .read()
            .map(|display| display.show_methylation_sites())
            .unwrap_or(false);
        if !show_methylation {
            return;
        }
        if !detail.show_methylation_sites {
            return;
        }
        let y = self.baseline_y();
        let sites = self
            .dna
            .read()
            .map(|dna| dna.methylation_sites().clone())
            .unwrap_or_default();
        for site in sites.sites() {
            if *site < viewport.start || *site >= viewport.end {
                continue;
            }
            let x = self.bp_to_x(*site, viewport);
            painter.line_segment(
                [Pos2::new(x, y - METHYLATION_TICK), Pos2::new(x, y - 1.0)],
                Stroke::new(1.0, Color32::DARK_RED),
            );
        }
    }

    fn draw_features(&self, painter: &egui::Painter, detail: LinearDetailLevel) {
        let show_features = self
            .display
            .read()
            .map(|display| display.show_features())
            .unwrap_or(false);
        if !show_features {
            return;
        }

        for feature in &self.features {
            let selected = self.selected_feature_number == Some(feature.feature_number);
            let hovered = self.hovered_feature_number == Some(feature.feature_number);

            painter.rect_filled(feature.rect, 1.5, feature.color);

            if feature.is_pointy {
                let tip_x = if feature.is_reverse {
                    feature.rect.left()
                } else {
                    feature.rect.right()
                };
                let y = feature.rect.center().y;
                let h = feature.rect.height() / 2.0;
                let tip = if feature.is_reverse {
                    Pos2::new(tip_x - 6.0, y)
                } else {
                    Pos2::new(tip_x + 6.0, y)
                };
                painter.add(egui::Shape::convex_polygon(
                    vec![Pos2::new(tip_x, y - h), tip, Pos2::new(tip_x, y + h)],
                    feature.color,
                    Stroke::NONE,
                ));
            }

            if selected || hovered {
                let stroke = if selected {
                    Stroke::new(2.0, Color32::YELLOW)
                } else {
                    Stroke::new(1.5, Color32::WHITE)
                };
                painter.rect_stroke(feature.rect.expand(1.0), 2.0, stroke);
            }

            if !detail.show_feature_labels {
                continue;
            }
            let align = if feature.label_pos.y < self.baseline_y() {
                Align2::LEFT_BOTTOM
            } else {
                Align2::LEFT_TOP
            };
            painter.text(
                feature.label_pos,
                align,
                &feature.label,
                FontId {
                    size: 10.0,
                    family: FontFamily::Monospace,
                },
                Color32::BLACK,
            );
        }
    }

    fn draw_restriction_enzyme_sites(
        &mut self,
        painter: &egui::Painter,
        viewport: LinearViewport,
        detail: LinearDetailLevel,
    ) {
        self.restriction_enzyme_sites.clear();

        if !self
            .display
            .read()
            .map(|display| display.show_restriction_enzyme_sites())
            .unwrap_or(false)
        {
            return;
        }
        if !detail.show_restriction_sites {
            painter.text(
                Pos2::new(self.area.left() + 6.0, self.area.bottom() - 6.0),
                Align2::LEFT_BOTTOM,
                "Restriction sites hidden at this zoom; zoom in to inspect cut sites.",
                FontId {
                    size: 10.0,
                    family: FontFamily::Monospace,
                },
                Color32::DARK_GRAY,
            );
            return;
        }

        let groups = self
            .dna
            .read()
            .map(|dna| dna.restriction_enzyme_groups().clone())
            .unwrap_or_default();

        let mut keys: Vec<_> = groups.keys().cloned().collect();
        keys.sort();
        let mut top_label_lanes: Vec<f32> = vec![];
        let mut bottom_label_lanes: Vec<f32> = vec![];

        for (idx, key) in keys.iter().enumerate() {
            let names = match groups.get(key) {
                Some(names) => names,
                None => continue,
            };
            let pos = self.normalize_pos(key.pos());
            if pos < viewport.start || pos >= viewport.end {
                continue;
            }
            let x = self.bp_to_x(pos, viewport);
            let y = self.baseline_y();
            let color = DnaDisplay::restriction_enzyme_group_color(key.number_of_cuts());

            painter.line_segment(
                [Pos2::new(x, y - 8.0), Pos2::new(x, y + 8.0)],
                Stroke::new(1.0, color),
            );

            if let Some(hovered) = &self.hover_enzyme {
                if hovered.key == *key {
                    painter.rect_filled(hovered.area, 1.0, Color32::LIGHT_YELLOW);
                }
            }
            let tick_rect = Rect::from_center_size(Pos2::new(x, y), Vec2::new(6.0, 18.0));
            let area = if detail.show_restriction_labels {
                let label = names.join(",");
                let label_width = Self::estimate_label_width(&label);
                let label_left = x - label_width / 2.0;
                let label_right = x + label_width / 2.0;
                let place_top = idx % 2 == 0;
                let label_lane = if place_top {
                    Self::allocate_lane(&mut top_label_lanes, label_left, label_right, 6.0)
                } else {
                    Self::allocate_lane(&mut bottom_label_lanes, label_left, label_right, 6.0)
                };
                let label_y = if place_top {
                    y - RE_LABEL_BASE_OFFSET - LABEL_ROW_HEIGHT * label_lane as f32
                } else {
                    y + RE_LABEL_BASE_OFFSET + LABEL_ROW_HEIGHT * label_lane as f32
                };
                let align = if place_top {
                    Align2::CENTER_BOTTOM
                } else {
                    Align2::CENTER_TOP
                };
                let text_rect = painter.text(
                    Pos2::new(x, label_y),
                    align,
                    label,
                    FontId {
                        size: 9.0,
                        family: FontFamily::Monospace,
                    },
                    color,
                );
                text_rect.expand(2.0).union(tick_rect)
            } else {
                tick_rect
            };
            self.restriction_enzyme_sites
                .push(RestrictionEnzymePosition {
                    area,
                    key: key.clone(),
                });
        }
    }

    pub fn render(&mut self, ui: &mut egui::Ui, area: Rect) {
        self.sequence_length = self.dna.read().map(|dna| dna.len()).unwrap_or(0);
        let area_changed = self.area != area;
        self.area = area;
        let viewport = self.viewport();
        let detail = self.detail_level(viewport);

        if (area_changed || self.layout_needs_recomputing()) && Self::is_rect_usable(self.area) {
            self.layout_features(viewport);
            self.layout_was_updated();
        }
        if !Self::is_rect_usable(self.area) {
            return;
        }

        self.draw_name_and_length(ui.painter(), viewport);
        self.draw_gc_contents(ui.painter(), viewport);
        self.draw_methylation_sites(ui.painter(), viewport, detail);
        self.draw_backbone(ui.painter());
        self.draw_bp_ticks(ui.painter(), viewport);
        self.draw_open_reading_frames(ui.painter(), viewport);
        self.draw_features(ui.painter(), detail);
        self.draw_restriction_enzyme_sites(ui.painter(), viewport, detail);
    }
}
