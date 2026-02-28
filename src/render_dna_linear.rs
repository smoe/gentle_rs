//! Linear-map DNA renderer.

use crate::{
    dna_display::{DnaDisplay, Selection, TfbsDisplayCriteria, VcfDisplayCriteria},
    dna_sequence::DNAsequence,
    engine::LinearSequenceLetterLayoutMode,
    feature_location::{collect_location_ranges_usize, feature_is_reverse},
    gc_contents::GcContents,
    iupac_code::IupacCode,
    open_reading_frame::OpenReadingFrame,
    render_dna::RenderDna,
    render_dna::RestrictionEnzymePosition,
};
use eframe::egui::{
    self, Align2, Color32, FontFamily, FontId, PointerState, Pos2, Rect, Stroke, Vec2,
};
use gb_io::seq::Feature;
use std::collections::{BTreeSet, HashMap};
use std::sync::{Arc, RwLock};

const BASELINE_STROKE: f32 = 2.0;
const FEATURE_HEIGHT: f32 = 12.0;
const FEATURE_GAP: f32 = 14.0;
const BASELINE_MARGIN: f32 = 24.0;
const MIN_FEATURE_HEIGHT: f32 = 6.0;
const MIN_FEATURE_GAP: f32 = 4.0;
const MIN_BASELINE_MARGIN: f32 = 8.0;
const BASELINE_TOP_PADDING: f32 = 26.0;
const BASELINE_BOTTOM_PADDING: f32 = 24.0;
const BASELINE_SIDE_MIN_EXTENT: f32 = 20.0;
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
const ORF_MAX_BP_PER_PX: f32 = 18.0;
const SEQUENCE_BASE_TEXT_MIN_FONT_SIZE: f32 = 8.0;
const SEQUENCE_BASE_TEXT_MAX_FONT_SIZE: f32 = 14.0;
const REGULATORY_BASELINE_MARGIN: f32 = 5.0;
const REGULATORY_FEATURE_GAP: f32 = 6.0;
const REGULATORY_FEATURE_HEIGHT: f32 = 6.0;
const REGULATORY_GROUP_GAP: f32 = 5.0;
const MIN_REGULATORY_BASELINE_MARGIN: f32 = 2.0;
const MIN_REGULATORY_FEATURE_GAP: f32 = 3.0;
const MIN_REGULATORY_FEATURE_HEIGHT: f32 = 4.0;
const FEATURE_LANE_PADDING: f32 = 5.0;
const REGULATORY_LANE_PADDING: f32 = 2.0;
const HELICAL_MIN_X_COMPRESSION: f32 = 0.2;
const HELICAL_MAX_Y_PHASE_SCALE: f32 = 0.8;
const HELICAL_MIN_VISUAL_PROGRESS: f32 = 0.35;
const CONDENSED_HELICAL_ROW_SPACING_SCALE: f32 = 1.05;
const CONDENSED_HELICAL_ROW_COUNT: usize = 10;
const CONDENSED_HELICAL_SINGLE_ROW_MIN_STEP: f32 = 12.0;
const CONDENSED_HELICAL_ROW_PAIR_PADDING: f32 = 4.0;
const CONDENSED_HELICAL_STRAND_GAP: f32 = 6.0;
const CONDENSED_FEATURE_CLEARANCE_PADDING: f32 = 8.0;
const FEATURE_INLINE_LABEL_FONT_SIZE: f32 = 10.0;
const FEATURE_EXTERNAL_LABEL_FONT_SIZE: f32 = 11.0;
const FEATURE_EXTERNAL_LABEL_HEIGHT: f32 = 16.0;

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
    show_open_reading_frames: bool,
}

#[derive(Debug, Clone)]
struct FeaturePosition {
    feature_number: usize,
    from: usize,
    to: usize,
    label: String,
    kind_upper: String,
    color: Color32,
    is_regulatory: bool,
    is_pointy: bool,
    is_reverse: bool,
    rect: Rect,
    exon_rects: Vec<Rect>,
    intron_connectors: Vec<[Pos2; 3]>,
}

impl FeaturePosition {
    fn contains(&self, pos: Pos2) -> bool {
        self.exon_rects.iter().any(|rect| rect.contains(pos))
    }
}

#[derive(Debug, Clone, Copy)]
struct SideLaneStyle {
    margin: f32,
    gap: f32,
    height: f32,
}

#[derive(Debug, Clone, Copy)]
struct HelicalLayoutParams {
    enabled: bool,
    mode: LinearSequenceLetterLayoutMode,
    active_max_span: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SequenceBaseRenderMode {
    Off,
    StandardLinear,
    ContinuousHelical,
    Condensed10Row,
}

#[derive(Debug, Clone, Copy)]
struct SequenceBaseRenderStatus {
    active_mode: SequenceBaseRenderMode,
    configured_layout: LinearSequenceLetterLayoutMode,
    standard_max_span: usize,
    helical_enabled: bool,
    active_helical_max_span: usize,
    span: usize,
}

impl SequenceBaseRenderStatus {
    fn bases_visible(&self) -> bool {
        self.active_mode != SequenceBaseRenderMode::Off
    }

    fn active_mode_label(&self) -> &'static str {
        match self.active_mode {
            SequenceBaseRenderMode::Off => "OFF",
            SequenceBaseRenderMode::StandardLinear => "STANDARD",
            SequenceBaseRenderMode::ContinuousHelical => "HELICAL",
            SequenceBaseRenderMode::Condensed10Row => "CONDENSED-10",
        }
    }

    fn configured_layout_label(&self) -> &'static str {
        match self.configured_layout {
            LinearSequenceLetterLayoutMode::ContinuousHelical => "continuous-helical",
            LinearSequenceLetterLayoutMode::Condensed10Row => "condensed-10-row",
        }
    }

    fn note_text(&self) -> Option<String> {
        if self.active_mode == SequenceBaseRenderMode::StandardLinear
            && self.helical_enabled
            && self.active_helical_max_span > 0
        {
            return Some(format!(
                "{} configured; activates when span > {} bp (and <= {} bp)",
                self.configured_layout_label(),
                self.standard_max_span,
                self.active_helical_max_span
            ));
        }
        None
    }

    fn off_reason_text(&self) -> String {
        if self.span == 0 {
            return "inactive: empty viewport".to_string();
        }
        if !self.helical_enabled {
            if self.standard_max_span == 0 {
                return "inactive: standard max is 0 and helical letters are disabled".to_string();
            }
            return format!(
                "inactive: span {} bp > standard max {} bp and helical letters are disabled",
                self.span, self.standard_max_span
            );
        }
        if self.active_helical_max_span == 0 {
            return "inactive: selected layout max span is 0".to_string();
        }
        if self.span > self.active_helical_max_span {
            return format!(
                "inactive: span {} bp > active layout max {} bp",
                self.span, self.active_helical_max_span
            );
        }
        "inactive: suppressed by current thresholds".to_string()
    }
}

#[derive(Debug, Clone)]
pub struct RenderDnaLinear {
    dna: Arc<RwLock<DNAsequence>>,
    display: Arc<RwLock<DnaDisplay>>,
    sequence_length: usize,
    area: Rect,
    baseline_y: f32,
    features: Vec<FeaturePosition>,
    restriction_enzyme_sites: Vec<RestrictionEnzymePosition>,
    selected_feature_number: Option<usize>,
    selected_enzyme: Option<RestrictionEnzymePosition>,
    hovered_feature_number: Option<usize>,
    hover_enzyme: Option<RestrictionEnzymePosition>,
    external_labeled_feature_numbers: BTreeSet<usize>,
}

impl RenderDnaLinear {
    fn helical_layout_params(display: &DnaDisplay) -> HelicalLayoutParams {
        let mode = display.linear_sequence_letter_layout_mode();
        let active_max_span = match mode {
            LinearSequenceLetterLayoutMode::ContinuousHelical => {
                display.linear_sequence_helical_max_view_span_bp()
            }
            LinearSequenceLetterLayoutMode::Condensed10Row => {
                display.linear_sequence_condensed_max_view_span_bp()
            }
        };
        HelicalLayoutParams {
            enabled: display.linear_sequence_helical_letters_enabled(),
            mode,
            active_max_span,
        }
    }

    fn sequence_base_render_status(&self, viewport: LinearViewport) -> SequenceBaseRenderStatus {
        let (standard_max_span, helical) = self
            .display
            .read()
            .map(|display| {
                (
                    display.linear_sequence_base_text_max_view_span_bp(),
                    Self::helical_layout_params(&display),
                )
            })
            .unwrap_or((
                500,
                HelicalLayoutParams {
                    enabled: false,
                    mode: LinearSequenceLetterLayoutMode::ContinuousHelical,
                    active_max_span: 1500,
                },
            ));
        let active_mode = if viewport.span == 0 {
            SequenceBaseRenderMode::Off
        } else if standard_max_span > 0 && viewport.span <= standard_max_span {
            SequenceBaseRenderMode::StandardLinear
        } else if helical.enabled
            && helical.active_max_span > 0
            && viewport.span <= helical.active_max_span
        {
            match helical.mode {
                LinearSequenceLetterLayoutMode::ContinuousHelical => {
                    SequenceBaseRenderMode::ContinuousHelical
                }
                LinearSequenceLetterLayoutMode::Condensed10Row => {
                    SequenceBaseRenderMode::Condensed10Row
                }
            }
        } else {
            SequenceBaseRenderMode::Off
        };
        SequenceBaseRenderStatus {
            active_mode,
            configured_layout: helical.mode,
            standard_max_span,
            helical_enabled: helical.enabled,
            active_helical_max_span: helical.active_max_span,
            span: viewport.span,
        }
    }

    pub fn new(dna: Arc<RwLock<DNAsequence>>, display: Arc<RwLock<DnaDisplay>>) -> Self {
        let sequence_length = dna.read().map(|d| d.len()).unwrap_or(0);
        Self {
            dna,
            display,
            sequence_length,
            area: Rect::NOTHING,
            baseline_y: 0.0,
            features: vec![],
            restriction_enzyme_sites: vec![],
            selected_feature_number: None,
            selected_enzyme: None,
            hovered_feature_number: None,
            hover_enzyme: None,
            external_labeled_feature_numbers: BTreeSet::new(),
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
        if self.baseline_y.is_finite() && self.baseline_y > 0.0 {
            self.baseline_y
        } else {
            self.area.center().y
        }
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
            show_open_reading_frames: bp_per_px <= ORF_MAX_BP_PER_PX,
        }
    }

    fn low_value_feature_min_width_px(bp_per_px: f32) -> f32 {
        if bp_per_px <= 4.0 {
            1.5
        } else if bp_per_px <= 10.0 {
            2.5
        } else if bp_per_px <= 20.0 {
            4.0
        } else {
            7.0
        }
    }

    fn feature_label_color(fill: Color32) -> Color32 {
        let luminance = 0.299 * fill.r() as f32 + 0.587 * fill.g() as f32 + 0.114 * fill.b() as f32;
        if luminance < 145.0 {
            Color32::WHITE
        } else {
            Color32::BLACK
        }
    }

    fn draw_feature(
        feature: &Feature,
        show_cds_features: bool,
        show_gene_features: bool,
        show_mrna_features: bool,
        show_tfbs: bool,
        tfbs_display_criteria: TfbsDisplayCriteria,
        vcf_display_criteria: VcfDisplayCriteria,
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
        if RenderDna::is_vcf_track_feature(feature) {
            return RenderDna::vcf_feature_passes_display_filter(feature, &vcf_display_criteria);
        }
        true
    }

    fn estimate_label_width(label: &str) -> f32 {
        let chars = label.chars().count().max(1) as f32;
        chars * LABEL_CHAR_WIDTH
    }

    fn truncate_label(label: &str, max_chars: usize) -> String {
        if max_chars == 0 {
            return String::new();
        }
        let total = label.chars().count();
        if total <= max_chars {
            return label.to_string();
        }
        let keep = max_chars.saturating_sub(1);
        let mut out = label.chars().take(keep).collect::<String>();
        out.push('…');
        out
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

    fn side_extent(lane_count: usize, style: SideLaneStyle) -> f32 {
        if lane_count == 0 {
            return 0.0;
        }
        style.margin + (lane_count.saturating_sub(1) as f32) * style.gap + style.height * 0.5
    }

    fn fit_side_style(lane_count: usize, target_extent: f32) -> SideLaneStyle {
        Self::fit_lane_style(
            lane_count,
            target_extent,
            SideLaneStyle {
                margin: BASELINE_MARGIN,
                gap: FEATURE_GAP,
                height: FEATURE_HEIGHT,
            },
            MIN_BASELINE_MARGIN,
            MIN_FEATURE_GAP,
            MIN_FEATURE_HEIGHT,
        )
    }

    fn fit_lane_style(
        lane_count: usize,
        target_extent: f32,
        default: SideLaneStyle,
        min_margin: f32,
        min_gap: f32,
        min_height: f32,
    ) -> SideLaneStyle {
        if lane_count == 0 {
            return default;
        }
        let natural = Self::side_extent(lane_count, default);
        if natural <= 0.0 || target_extent <= 0.0 || natural <= target_extent {
            return default;
        }
        let scale = (target_extent / natural).clamp(0.1, 1.0);
        SideLaneStyle {
            margin: (default.margin * scale).max(min_margin),
            gap: (default.gap * scale).max(min_gap),
            height: (default.height * scale).max(min_height),
        }
    }

    fn condensed_row_step(font_size: f32, show_double_strand: bool) -> f32 {
        let single = (font_size * CONDENSED_HELICAL_ROW_SPACING_SCALE)
            .max(CONDENSED_HELICAL_SINGLE_ROW_MIN_STEP);
        if !show_double_strand {
            return single;
        }
        single.max(
            font_size * 2.0 + CONDENSED_HELICAL_STRAND_GAP + CONDENSED_HELICAL_ROW_PAIR_PADDING,
        )
    }

    fn condensed_helical_text_half_extent(&self, viewport: LinearViewport) -> Option<f32> {
        let status = self.sequence_base_render_status(viewport);
        if status.active_mode != SequenceBaseRenderMode::Condensed10Row {
            return None;
        }
        let show_double_strand = self
            .display
            .read()
            .map(|display| display.linear_show_double_strand_bases())
            .unwrap_or(true);
        let px_per_bp = self.area.width().max(1.0) / viewport.span.max(1) as f32;
        let font_size = (px_per_bp * 0.85).clamp(
            SEQUENCE_BASE_TEXT_MIN_FONT_SIZE,
            SEQUENCE_BASE_TEXT_MAX_FONT_SIZE,
        );
        let row_step = Self::condensed_row_step(font_size, show_double_strand);
        let row_offset_max =
            row_step * (CONDENSED_HELICAL_ROW_COUNT.saturating_sub(1) as f32 * 0.5);
        let strand_half_gap = if show_double_strand {
            CONDENSED_HELICAL_STRAND_GAP * 0.5
        } else {
            0.0
        };
        Some(font_size + strand_half_gap + row_offset_max)
    }

    fn layout_features(&mut self, viewport: LinearViewport) {
        self.features.clear();
        self.baseline_y = self.area.center().y;
        if self.sequence_length == 0 {
            return;
        }
        let bp_per_px = self.bp_per_px(viewport);
        let low_value_feature_min_width_px = Self::low_value_feature_min_width_px(bp_per_px);

        #[derive(Clone)]
        struct Seed {
            feature_number: usize,
            from: usize,
            to: usize,
            exon_segments: Vec<(f32, f32)>,
            x1: f32,
            x2: f32,
            label: String,
            kind_upper: String,
            color: Color32,
            is_pointy: bool,
            is_reverse: bool,
            is_regulatory: bool,
        }

        #[derive(Clone)]
        struct PositionedSeed {
            seed: Seed,
            feature_lane: usize,
            lane_side: LaneSide,
        }

        #[derive(Clone, Copy)]
        enum LaneSide {
            Top,
            Bottom,
            RegulatoryTop,
            RegulatoryNearBaseline,
        }

        let mut seeds: Vec<Seed> = Vec::new();
        let (
            show_cds_features,
            show_gene_features,
            show_mrna_features,
            show_tfbs,
            tfbs_display_criteria,
            vcf_display_criteria,
            regulatory_tracks_near_baseline,
            regulatory_feature_max_view_span_bp,
            hidden_feature_kinds,
        ) = self
            .display
            .read()
            .map(|display| {
                (
                    display.show_cds_features_effective(),
                    display.show_gene_features(),
                    display.show_mrna_features(),
                    display.show_tfbs(),
                    display.tfbs_display_criteria(),
                    display.vcf_display_criteria(),
                    display.regulatory_tracks_near_baseline(),
                    display.regulatory_feature_max_view_span_bp(),
                    display.hidden_feature_kinds().clone(),
                )
            })
            .unwrap_or((
                true,
                true,
                true,
                false,
                TfbsDisplayCriteria::default(),
                VcfDisplayCriteria::default(),
                false,
                50_000,
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
                vcf_display_criteria.clone(),
                &hidden_feature_kinds,
            ) {
                continue;
            }
            if !RenderDna::feature_visible_for_view_span(
                feature,
                viewport.span,
                regulatory_feature_max_view_span_bp,
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
            if to > self.sequence_length {
                to = self.sequence_length;
            }
            if to <= from {
                continue;
            }

            let label = RenderDna::feature_name(feature);
            let mut exon_ranges: Vec<(usize, usize)> = vec![];
            collect_location_ranges_usize(&feature.location, &mut exon_ranges);
            if exon_ranges.is_empty() {
                exon_ranges.push((from, to));
            }
            exon_ranges.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

            let mut exon_segments: Vec<(f32, f32)> = Vec::new();
            for (exon_start, exon_end_exclusive) in exon_ranges {
                if exon_start >= self.sequence_length {
                    continue;
                }
                let exon_end_exclusive = exon_end_exclusive.min(self.sequence_length);
                if exon_end_exclusive <= exon_start {
                    continue;
                }
                let Some((visible_start, visible_end)) = Self::range_overlap(
                    exon_start,
                    exon_end_exclusive,
                    viewport.start,
                    viewport.end,
                ) else {
                    continue;
                };
                let seg_x1 = self.bp_to_x(visible_start, viewport).max(self.area.left());
                let seg_x2 = self
                    .bp_to_x(visible_end, viewport)
                    .max(seg_x1 + 1.0)
                    .min(self.area.right());
                exon_segments.push((seg_x1, seg_x2));
            }
            if exon_segments.is_empty() {
                continue;
            }
            exon_segments.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));

            let x1 = exon_segments
                .iter()
                .map(|(sx1, _)| *sx1)
                .fold(f32::INFINITY, f32::min);
            let x2 = exon_segments
                .iter()
                .map(|(_, sx2)| *sx2)
                .fold(f32::NEG_INFINITY, f32::max);
            if !x1.is_finite() || !x2.is_finite() {
                continue;
            }
            let kind = feature.kind.to_string().to_ascii_uppercase();
            let is_high_priority_feature = matches!(kind.as_str(), "CDS" | "GENE" | "MRNA");
            if !is_high_priority_feature && (x2 - x1) < low_value_feature_min_width_px {
                continue;
            }

            seeds.push(Seed {
                feature_number,
                from,
                to,
                exon_segments,
                x1,
                x2,
                label,
                kind_upper: kind.clone(),
                color: RenderDna::feature_color(feature),
                is_pointy: RenderDna::is_feature_pointy(feature),
                is_reverse: feature_is_reverse(feature),
                is_regulatory: RenderDna::is_regulatory_feature(feature),
            });
        }

        let mut feature_lanes_top: Vec<f32> = vec![];
        let mut feature_lanes_bottom: Vec<f32> = vec![];
        let mut feature_lanes_regulatory_top: Vec<f32> = vec![];
        let has_regulatory_features = seeds.iter().any(|seed| seed.is_regulatory);
        let has_restriction_sites = self
            .display
            .read()
            .map(|display| display.show_restriction_enzyme_sites())
            .unwrap_or(false)
            && self
                .dna
                .read()
                .map(|dna| !dna.restriction_enzyme_groups().is_empty())
                .unwrap_or(false);
        let compact_feature_spacing = !has_regulatory_features && !has_restriction_sites;
        let mut lane_seed: Vec<PositionedSeed> = Vec::with_capacity(seeds.len());
        let mut lane_order = seeds.clone();
        lane_order.sort_by(|a, b| {
            a.is_regulatory
                .cmp(&b.is_regulatory)
                .then_with(|| a.x1.total_cmp(&b.x1))
                .then_with(|| (b.x2 - b.x1).total_cmp(&(a.x2 - a.x1)))
                .then_with(|| a.feature_number.cmp(&b.feature_number))
        });

        for seed in lane_order {
            let lane_padding = if seed.is_regulatory {
                REGULATORY_LANE_PADDING
            } else {
                FEATURE_LANE_PADDING
            };
            let (lane_side, feature_lane) = if seed.is_regulatory {
                if regulatory_tracks_near_baseline {
                    (LaneSide::RegulatoryNearBaseline, 0)
                } else {
                    (
                        LaneSide::RegulatoryTop,
                        Self::allocate_lane(
                            &mut feature_lanes_regulatory_top,
                            seed.x1,
                            seed.x2,
                            lane_padding,
                        ),
                    )
                }
            } else if seed.is_reverse {
                (
                    LaneSide::Bottom,
                    Self::allocate_lane(&mut feature_lanes_bottom, seed.x1, seed.x2, lane_padding),
                )
            } else {
                (
                    LaneSide::Top,
                    Self::allocate_lane(&mut feature_lanes_top, seed.x1, seed.x2, lane_padding),
                )
            };
            lane_seed.push(PositionedSeed {
                seed,
                feature_lane,
                lane_side,
            });
        }

        let top_lane_count = feature_lanes_top.len();
        let bottom_lane_count = feature_lanes_bottom.len();
        let regulatory_top_lane_count = feature_lanes_regulatory_top.len();
        let default_style = SideLaneStyle {
            margin: if compact_feature_spacing {
                MIN_BASELINE_MARGIN
            } else {
                BASELINE_MARGIN
            },
            gap: FEATURE_GAP,
            height: FEATURE_HEIGHT,
        };
        let regulatory_top_default_style = SideLaneStyle {
            margin: REGULATORY_BASELINE_MARGIN,
            gap: REGULATORY_FEATURE_GAP,
            height: REGULATORY_FEATURE_HEIGHT,
        };
        let top_regular_natural_extent = Self::side_extent(top_lane_count, default_style);
        let regulatory_top_natural_extent =
            Self::side_extent(regulatory_top_lane_count, regulatory_top_default_style);
        let regulatory_group_gap_natural = if !regulatory_tracks_near_baseline
            && top_lane_count > 0
            && regulatory_top_lane_count > 0
        {
            REGULATORY_GROUP_GAP
        } else {
            0.0
        };
        let top_natural_extent = if regulatory_tracks_near_baseline {
            top_regular_natural_extent.max(BASELINE_SIDE_MIN_EXTENT)
        } else {
            (top_regular_natural_extent
                + regulatory_group_gap_natural
                + regulatory_top_natural_extent)
                .max(BASELINE_SIDE_MIN_EXTENT)
        };
        let bottom_natural_extent =
            Self::side_extent(bottom_lane_count, default_style).max(BASELINE_SIDE_MIN_EXTENT);
        let natural_total = (top_natural_extent + bottom_natural_extent).max(1.0);
        let condensed_text_half_extent = self.condensed_helical_text_half_extent(viewport);
        let condensed_band_total = condensed_text_half_extent.unwrap_or(0.0) * 2.0;
        let available_total = (self.area.height()
            - BASELINE_TOP_PADDING
            - BASELINE_BOTTOM_PADDING
            - condensed_band_total)
            .max(40.0);
        let top_target_extent = available_total * (top_natural_extent / natural_total);
        let bottom_target_extent = available_total * (bottom_natural_extent / natural_total);

        let reserved_half = condensed_text_half_extent.unwrap_or(0.0);
        let baseline = self.area.top() + BASELINE_TOP_PADDING + reserved_half + top_target_extent;
        let baseline_min = self.area.top() + BASELINE_TOP_PADDING + reserved_half;
        let baseline_max = self.area.bottom() - BASELINE_BOTTOM_PADDING - reserved_half;
        self.baseline_y = if baseline_min <= baseline_max {
            baseline.clamp(baseline_min, baseline_max)
        } else {
            self.area.center().y
        };

        let top_scale = if top_natural_extent > 0.0 {
            (top_target_extent / top_natural_extent).clamp(0.1, 1.0)
        } else {
            1.0
        };
        let top_target_regular_extent = if regulatory_tracks_near_baseline {
            top_target_extent
        } else {
            top_regular_natural_extent * top_scale
        };
        let top_target_regulatory_extent = if regulatory_tracks_near_baseline {
            top_target_extent
        } else {
            regulatory_top_natural_extent * top_scale
        };
        let regulatory_group_gap = regulatory_group_gap_natural * top_scale;

        let mut top_style = Self::fit_side_style(top_lane_count, top_target_regular_extent);
        let mut bottom_style = Self::fit_side_style(bottom_lane_count, bottom_target_extent);
        let mut regulatory_top_style = Self::fit_lane_style(
            regulatory_top_lane_count,
            top_target_regulatory_extent,
            regulatory_top_default_style,
            MIN_REGULATORY_BASELINE_MARGIN,
            MIN_REGULATORY_FEATURE_GAP,
            MIN_REGULATORY_FEATURE_HEIGHT,
        );
        let condensed_feature_clearance = self
            .condensed_helical_text_half_extent(viewport)
            .map(|half_extent| half_extent + CONDENSED_FEATURE_CLEARANCE_PADDING);
        if let Some(min_clearance) = condensed_feature_clearance {
            top_style.margin = top_style.margin.max(min_clearance + top_style.height * 0.5);
            bottom_style.margin = bottom_style
                .margin
                .max(min_clearance + bottom_style.height * 0.5);
            regulatory_top_style.margin = regulatory_top_style
                .margin
                .max(min_clearance + regulatory_top_style.height * 0.5);
        }
        let regulatory_near_baseline_margin = condensed_feature_clearance
            .map(|min_clearance| {
                (min_clearance + REGULATORY_FEATURE_HEIGHT * 0.5).max(REGULATORY_BASELINE_MARGIN)
            })
            .unwrap_or(REGULATORY_BASELINE_MARGIN);
        let top_regular_extent = Self::side_extent(top_lane_count, top_style);

        lane_seed.sort_by(|a, b| {
            a.seed
                .x1
                .total_cmp(&b.seed.x1)
                .then_with(|| a.seed.feature_number.cmp(&b.seed.feature_number))
        });

        for item in lane_seed {
            let seed = item.seed;
            let feature_lane = item.feature_lane;
            let side_style = match item.lane_side {
                LaneSide::Top => top_style,
                LaneSide::Bottom => bottom_style,
                LaneSide::RegulatoryTop => regulatory_top_style,
                LaneSide::RegulatoryNearBaseline => SideLaneStyle {
                    margin: regulatory_near_baseline_margin,
                    gap: 0.0,
                    height: REGULATORY_FEATURE_HEIGHT,
                },
            };
            let center_y = match item.lane_side {
                LaneSide::Top => {
                    let mut required_clearance =
                        if regulatory_tracks_near_baseline && has_regulatory_features {
                            let regulatory_outer =
                                REGULATORY_BASELINE_MARGIN + REGULATORY_FEATURE_HEIGHT * 0.5;
                            regulatory_outer + side_style.height * 0.5 + 1.0
                        } else {
                            0.0
                        };
                    if let Some(min_clearance) = condensed_feature_clearance {
                        required_clearance =
                            required_clearance.max(min_clearance + side_style.height * 0.5);
                    }
                    self.baseline_y()
                        - side_style.margin.max(required_clearance)
                        - side_style.gap * feature_lane as f32
                }
                LaneSide::Bottom => {
                    self.baseline_y() + side_style.margin + side_style.gap * feature_lane as f32
                }
                LaneSide::RegulatoryTop => {
                    let regulatory_origin = if regulatory_tracks_near_baseline {
                        self.baseline_y()
                    } else {
                        self.baseline_y() - top_regular_extent - regulatory_group_gap
                    };
                    regulatory_origin - side_style.margin - side_style.gap * feature_lane as f32
                }
                LaneSide::RegulatoryNearBaseline => self.baseline_y() - side_style.margin,
            };
            let lane_is_bottom_side = matches!(item.lane_side, LaneSide::Bottom);
            let suppress_introns = matches!(item.lane_side, LaneSide::RegulatoryNearBaseline);
            let exon_rects: Vec<Rect> = seed
                .exon_segments
                .iter()
                .map(|(x1, x2)| {
                    Rect::from_min_max(
                        Pos2::new(*x1, center_y - side_style.height / 2.0),
                        Pos2::new(*x2, center_y + side_style.height / 2.0),
                    )
                })
                .collect();
            if exon_rects.is_empty() {
                continue;
            }
            let mut rect = exon_rects[0];
            for exon_rect in exon_rects.iter().skip(1) {
                rect = rect.union(*exon_rect);
            }
            let connector_delta = (side_style.height * 0.42).max(3.0);
            let connector_apex_y = if lane_is_bottom_side {
                center_y + connector_delta
            } else {
                center_y - connector_delta
            };
            let mut intron_connectors: Vec<[Pos2; 3]> = Vec::new();
            if !suppress_introns {
                for (left, right) in exon_rects.iter().zip(exon_rects.iter().skip(1)) {
                    if right.left() <= left.right() {
                        continue;
                    }
                    let start = Pos2::new(left.right(), center_y);
                    let end = Pos2::new(right.left(), center_y);
                    let apex_x = (start.x + end.x) * 0.5;
                    let apex = Pos2::new(apex_x, connector_apex_y);
                    intron_connectors.push([start, apex, end]);
                }
            }

            self.features.push(FeaturePosition {
                feature_number: seed.feature_number,
                from: seed.from,
                to: seed.to,
                label: seed.label,
                kind_upper: seed.kind_upper,
                color: seed.color,
                is_regulatory: seed.is_regulatory,
                is_pointy: seed.is_pointy && !seed.is_regulatory,
                is_reverse: seed.is_reverse,
                rect,
                exon_rects,
                intron_connectors,
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

    pub fn hovered_feature_number(&self) -> Option<usize> {
        self.hovered_feature_number
    }

    pub fn selected_restriction_enzyme(&self) -> Option<RestrictionEnzymePosition> {
        self.selected_enzyme.clone()
    }

    pub fn select_restriction_enzyme(&mut self, selected: Option<RestrictionEnzymePosition>) {
        let has_selected = selected.is_some();
        self.selected_enzyme = selected;
        if has_selected {
            self.selected_feature_number = None;
        }
    }

    pub fn select_feature(&mut self, feature_number: Option<usize>) {
        self.selected_feature_number = feature_number;
        if feature_number.is_some() {
            self.selected_enzyme = None;
        }
    }

    pub fn set_external_labeled_feature_numbers(&mut self, feature_numbers: BTreeSet<usize>) {
        self.external_labeled_feature_numbers = feature_numbers;
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
            if let Some(feature) = self.get_clicked_feature(pos) {
                self.selected_feature_number = Some(feature.feature_number);
                self.selected_enzyme = None;
            } else {
                self.selected_feature_number = None;
                self.selected_enzyme = self.get_re_site_for_positon(pos);
            }
        }
    }

    pub fn on_double_click(&mut self, pointer_state: PointerState) {
        if let Ok(mut display) = self.display.write() {
            display.deselect();
        }

        if let Some(pos) = pointer_state.latest_pos() {
            if let Some((feature_number, feature_from, feature_to)) = self
                .get_clicked_feature(pos)
                .map(|feature| (feature.feature_number, feature.from, feature.to))
            {
                self.selected_feature_number = Some(feature_number);
                self.selected_enzyme = None;
                let selection = Selection::new(feature_from, feature_to, self.sequence_length);
                if let Ok(mut display) = self.display.write() {
                    display.select(selection);
                }
            } else if let Some(re_pos) = self.get_re_site_for_positon(pos) {
                self.selected_feature_number = None;
                self.selected_enzyme = Some(re_pos.clone());
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

    fn should_draw_backbone(&self, viewport: LinearViewport) -> bool {
        let status = self.sequence_base_render_status(viewport);
        if status.active_mode == SequenceBaseRenderMode::Condensed10Row {
            return false;
        }
        let hide_when_bases_visible = self
            .display
            .read()
            .map(|display| display.linear_hide_backbone_when_sequence_bases_visible())
            .unwrap_or(false);
        !(hide_when_bases_visible && status.bases_visible())
    }

    fn draw_backbone(&self, painter: &egui::Painter, viewport: LinearViewport) {
        if !self.should_draw_backbone(viewport) {
            return;
        }
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
        let show_sequence_bases = self.should_draw_sequence_bases(viewport);
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
            if !show_sequence_bases {
                painter.text(
                    Pos2::new(x, y + 7.0),
                    Align2::CENTER_TOP,
                    pos.to_string(),
                    font.clone(),
                    Color32::DARK_GRAY,
                );
            }
            pos += tick;
        }
    }

    fn should_draw_sequence_bases(&self, viewport: LinearViewport) -> bool {
        self.sequence_base_render_status(viewport).bases_visible()
    }

    fn sequence_base_color(base: u8) -> Color32 {
        match base.to_ascii_uppercase() {
            b'A' => Color32::from_rgb(20, 120, 40),
            b'C' => Color32::from_rgb(20, 80, 170),
            b'G' => Color32::from_rgb(160, 95, 20),
            b'T' | b'U' => Color32::from_rgb(170, 30, 30),
            _ => Color32::DARK_GRAY,
        }
    }

    fn draw_rotated_base_char(
        painter: &egui::Painter,
        x_center: f32,
        y_top: f32,
        base: char,
        font: &FontId,
        color: Color32,
    ) {
        let galley = painter.layout_no_wrap(base.to_string(), font.clone(), color);
        let size = galley.size();
        // Rotate the full glyph box (not symbol substitution) to preserve font fidelity.
        let pivot = Pos2::new(x_center + size.x * 0.5, y_top + size.y);
        painter.add(
            egui::epaint::TextShape::new(pivot, galley, color).with_angle(std::f32::consts::PI),
        );
    }

    fn draw_sequence_bases(&self, painter: &egui::Painter, viewport: LinearViewport) {
        let render_status = self.sequence_base_render_status(viewport);
        if !render_status.bases_visible() {
            return;
        }
        let Some(seq) = self.dna.read().ok().and_then(|dna| {
            dna.get_inclusive_range_safe(viewport.start..=viewport.end.saturating_sub(1))
        }) else {
            return;
        };
        let selection = self
            .display
            .read()
            .ok()
            .and_then(|display| display.selection());
        let (show_double_strand, reverse_strand_upside_down) = self
            .display
            .read()
            .map(|display| {
                (
                    display.linear_show_double_strand_bases(),
                    display.linear_reverse_strand_use_upside_down_letters(),
                )
            })
            .unwrap_or((true, true));
        let (standard_max_span, helical, helical_phase_offset_bp) = self
            .display
            .read()
            .map(|display| {
                (
                    display.linear_sequence_base_text_max_view_span_bp(),
                    Self::helical_layout_params(&display),
                    display.linear_sequence_helical_phase_offset_bp(),
                )
            })
            .unwrap_or((
                500,
                HelicalLayoutParams {
                    enabled: false,
                    mode: LinearSequenceLetterLayoutMode::ContinuousHelical,
                    active_max_span: 1500,
                },
                0,
            ));
        let use_helical_projection = matches!(
            render_status.active_mode,
            SequenceBaseRenderMode::ContinuousHelical | SequenceBaseRenderMode::Condensed10Row
        );
        let use_condensed_helical =
            render_status.active_mode == SequenceBaseRenderMode::Condensed10Row;
        let helical_t = Self::helical_projection_progress(
            viewport.span,
            standard_max_span,
            helical.active_max_span,
        );
        let helical_visual_t = if helical.enabled
            && helical.active_max_span > 0
            && viewport.span <= helical.active_max_span
        {
            helical_t.max(HELICAL_MIN_VISUAL_PROGRESS)
        } else {
            0.0
        };
        let px_per_bp = self.area.width().max(1.0) / viewport.span.max(1) as f32;
        let font_size = (px_per_bp * 0.85).clamp(
            SEQUENCE_BASE_TEXT_MIN_FONT_SIZE,
            SEQUENCE_BASE_TEXT_MAX_FONT_SIZE,
        );
        let font = FontId {
            size: font_size,
            family: FontFamily::Monospace,
        };
        let baseline = self.baseline_y();
        let strand_half_gap = CONDENSED_HELICAL_STRAND_GAP * 0.5;
        let forward_y = baseline - font_size - strand_half_gap;
        let reverse_y = baseline + strand_half_gap;
        let condensed_row_step = Self::condensed_row_step(font_size, show_double_strand);
        let condensed_row_top_index = CONDENSED_HELICAL_ROW_COUNT.saturating_sub(1);
        let condensed_row_center = condensed_row_top_index as f32 * 0.5;
        let helical_phase_row_divisor = condensed_row_top_index.max(1) as f32;
        let helical_compression =
            Self::sequence_x_compression(render_status.active_mode, helical_t);
        let helical_phase_scale = Self::helical_projection_y_phase_scale(helical_visual_t);
        let x_centerline = self.area.center().x;
        let project_x = |x: f32| -> f32 {
            if use_helical_projection {
                x_centerline + (x - x_centerline) * helical_compression
            } else {
                x
            }
        };
        for (offset, base) in seq.iter().enumerate() {
            let bp = viewport.start + offset;
            let x1 = self.bp_to_x(bp, viewport);
            let x2 = self.bp_to_x(bp.saturating_add(1), viewport);
            let x1_projected = project_x(x1);
            let x2_projected = project_x(x2);
            let x_center = (x1_projected + x2_projected) * 0.5;
            let helical_phase = if use_condensed_helical {
                let row_from_bottom =
                    Self::helical_phase_row_from_bottom(bp, helical_phase_offset_bp);
                let row_from_top = condensed_row_top_index.saturating_sub(row_from_bottom);
                let centered = row_from_top as f32 - condensed_row_center;
                centered * condensed_row_step
            } else if helical_visual_t > 0.0 {
                let row_from_bottom =
                    Self::helical_phase_row_from_bottom(bp, helical_phase_offset_bp);
                let normalized = (row_from_bottom as f32 / helical_phase_row_divisor) - 0.5;
                normalized * font_size * helical_phase_scale
            } else {
                0.0
            };
            let forward_y_bp = forward_y + helical_phase;
            let reverse_y_bp = if use_condensed_helical {
                // Keep both strands phase-aligned in condensed mode to avoid row-crossing overlap.
                reverse_y + helical_phase
            } else {
                reverse_y - helical_phase
            };
            let selection_min_y = forward_y_bp.min(reverse_y_bp) - 1.0;
            let selection_max_y = if show_double_strand {
                forward_y_bp.max(reverse_y_bp) + font_size + 1.0
            } else {
                forward_y_bp + font_size + 1.0
            };
            if let Some(selection) = &selection {
                if selection.contains(bp) {
                    let left = x1_projected.min(x2_projected);
                    let right = x1_projected.max(x2_projected).max(left + 1.0);
                    painter.rect_filled(
                        Rect::from_min_max(
                            Pos2::new(left, selection_min_y),
                            Pos2::new(right, selection_max_y),
                        ),
                        0.0,
                        Color32::from_gray(230),
                    );
                }
            }
            let forward_char = (*base as char).to_ascii_uppercase();
            painter.text(
                Pos2::new(x_center, forward_y_bp),
                Align2::CENTER_TOP,
                forward_char,
                font.clone(),
                Self::sequence_base_color(*base),
            );
            if show_double_strand {
                let complemented_byte = IupacCode::letter_complement(*base).to_ascii_uppercase();
                let complemented = complemented_byte as char;
                let reverse_color = Self::sequence_base_color(complemented_byte);
                if reverse_strand_upside_down {
                    Self::draw_rotated_base_char(
                        painter,
                        x_center,
                        reverse_y_bp,
                        complemented,
                        &font,
                        reverse_color,
                    );
                } else {
                    painter.text(
                        Pos2::new(x_center, reverse_y_bp),
                        Align2::CENTER_TOP,
                        complemented,
                        font.clone(),
                        reverse_color,
                    );
                }
            }
        }
    }

    fn helical_projection_progress(
        span: usize,
        standard_max_span: usize,
        helical_max_span: usize,
    ) -> f32 {
        if helical_max_span == 0 {
            return 0.0;
        }
        let start_span = standard_max_span.max(1).min(helical_max_span);
        if span <= start_span {
            return 0.0;
        }
        let range = helical_max_span.saturating_sub(start_span).max(1) as f32;
        ((span.saturating_sub(start_span)) as f32 / range).clamp(0.0, 1.0)
    }

    fn helical_projection_x_compression(progress: f32) -> f32 {
        let progress = progress.clamp(0.0, 1.0);
        1.0 - (1.0 - HELICAL_MIN_X_COMPRESSION) * progress
    }

    fn sequence_x_compression(render_mode: SequenceBaseRenderMode, helical_t: f32) -> f32 {
        if render_mode == SequenceBaseRenderMode::Condensed10Row {
            return 1.0;
        }
        Self::helical_projection_x_compression(helical_t)
    }

    fn helical_phase_row_from_bottom(bp: usize, offset_bp: usize) -> usize {
        bp.wrapping_add(offset_bp) % CONDENSED_HELICAL_ROW_COUNT
    }

    fn helical_projection_y_phase_scale(progress: f32) -> f32 {
        progress.clamp(0.0, 1.0) * HELICAL_MAX_Y_PHASE_SCALE
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
        let status = self.sequence_base_render_status(viewport);
        let status_color = match status.active_mode {
            SequenceBaseRenderMode::Condensed10Row => Color32::from_rgb(10, 115, 56),
            SequenceBaseRenderMode::ContinuousHelical => Color32::from_rgb(14, 98, 150),
            SequenceBaseRenderMode::StandardLinear => Color32::from_rgb(80, 80, 80),
            SequenceBaseRenderMode::Off => Color32::from_rgb(150, 30, 30),
        };
        painter.text(
            Pos2::new(self.area.right() - 6.0, self.area.top() + 21.0),
            Align2::RIGHT_TOP,
            format!("DNA MODE: {}", status.active_mode_label()),
            FontId {
                size: 10.0,
                family: FontFamily::Monospace,
            },
            status_color,
        );
        painter.text(
            Pos2::new(self.area.right() - 6.0, self.area.top() + 34.0),
            Align2::RIGHT_TOP,
            format!("layout: {}", status.configured_layout_label()),
            FontId {
                size: 9.0,
                family: FontFamily::Monospace,
            },
            Color32::from_gray(95),
        );
        if !status.bases_visible() {
            painter.text(
                Pos2::new(self.area.right() - 6.0, self.area.top() + 47.0),
                Align2::RIGHT_TOP,
                status.off_reason_text(),
                FontId {
                    size: 9.0,
                    family: FontFamily::Monospace,
                },
                Color32::from_gray(110),
            );
        } else if let Some(note) = status.note_text() {
            painter.text(
                Pos2::new(self.area.right() - 6.0, self.area.top() + 47.0),
                Align2::RIGHT_TOP,
                note,
                FontId {
                    size: 9.0,
                    family: FontFamily::Monospace,
                },
                Color32::from_gray(110),
            );
        }
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

    fn draw_open_reading_frames(
        &self,
        painter: &egui::Painter,
        viewport: LinearViewport,
        detail: LinearDetailLevel,
    ) {
        if !detail.show_open_reading_frames {
            return;
        }
        let show_orfs = self
            .display
            .read()
            .map(|display| display.show_open_reading_frames_effective())
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
        let (show_gc, gc_content_bin_size_bp) = self
            .display
            .read()
            .map(|display| (display.show_gc_contents(), display.gc_content_bin_size_bp()))
            .unwrap_or((false, 100));
        if !show_gc {
            return;
        }
        let y1 = self.baseline_y() - 4.0;
        let y2 = y1 + GC_STRIP_HEIGHT;
        let gc_contents = self
            .dna
            .read()
            .map(|dna| {
                GcContents::new_from_sequence_with_bin_size(
                    dna.forward_bytes(),
                    gc_content_bin_size_bp,
                )
            })
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
        let (show_features, external_font_size, external_bg_opacity) = self
            .display
            .read()
            .map(|display| {
                (
                    display.show_features(),
                    display.linear_external_feature_label_font_size(),
                    display.linear_external_feature_label_background_opacity(),
                )
            })
            .unwrap_or((false, FEATURE_EXTERNAL_LABEL_FONT_SIZE, 0.9));
        if !show_features {
            return;
        }
        let external_label_height =
            (external_font_size + 5.0).clamp(FEATURE_EXTERNAL_LABEL_HEIGHT, 28.0);
        let external_label_lane_step = (external_font_size + 2.0).clamp(LABEL_ROW_HEIGHT, 24.0);
        let external_bg_alpha = (external_bg_opacity.clamp(0.0, 1.0) * 255.0).round() as u8;

        let mut top_label_lanes: Vec<f32> = vec![];
        let mut bottom_label_lanes: Vec<f32> = vec![];

        for feature in &self.features {
            let selected = self.selected_feature_number == Some(feature.feature_number);
            let hovered = self.hovered_feature_number == Some(feature.feature_number);

            for exon_rect in &feature.exon_rects {
                painter.rect_filled(*exon_rect, 1.5, feature.color);
            }

            if feature.is_pointy {
                let tip_target = if feature.is_reverse {
                    feature.exon_rects.first().copied().unwrap_or(feature.rect)
                } else {
                    feature.exon_rects.last().copied().unwrap_or(feature.rect)
                };
                let tip_x = if feature.is_reverse {
                    tip_target.left()
                } else {
                    tip_target.right()
                };
                let y = tip_target.center().y;
                let h = tip_target.height() / 2.0;
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
            for connector in &feature.intron_connectors {
                painter.line_segment(
                    [connector[0], connector[1]],
                    Stroke::new(1.0, Color32::DARK_GRAY),
                );
                painter.line_segment(
                    [connector[1], connector[2]],
                    Stroke::new(1.0, Color32::DARK_GRAY),
                );
            }

            if selected || hovered {
                let stroke = if selected {
                    Stroke::new(2.0, Color32::YELLOW)
                } else {
                    Stroke::new(1.5, Color32::WHITE)
                };
                for exon_rect in &feature.exon_rects {
                    painter.rect_stroke(exon_rect.expand(1.0), 2.0, stroke);
                }
            }

            if !detail.show_feature_labels {
                continue;
            }
            if feature.label.trim().is_empty() {
                continue;
            }
            if matches!(feature.kind_upper.as_str(), "MRNA" | "EXON") {
                continue;
            }
            let label = Self::truncate_label(feature.label.trim(), 56);
            let label_rect = feature
                .exon_rects
                .iter()
                .copied()
                .max_by(|a, b| a.width().total_cmp(&b.width()))
                .unwrap_or(feature.rect);
            let label_width = Self::estimate_label_width(&label);
            let inline_possible = !feature.is_regulatory
                && label_rect.height() >= FEATURE_INLINE_LABEL_FONT_SIZE + 2.0
                && label_rect.width() >= label_width + 4.0;
            let force_external = selected
                || hovered
                || self
                    .external_labeled_feature_numbers
                    .contains(&feature.feature_number);
            if inline_possible && !force_external {
                let text_painter = painter.with_clip_rect(label_rect.shrink2(Vec2::new(1.0, 1.0)));
                text_painter.text(
                    label_rect.center(),
                    Align2::CENTER_CENTER,
                    &label,
                    FontId {
                        size: FEATURE_INLINE_LABEL_FONT_SIZE,
                        family: FontFamily::Monospace,
                    },
                    Self::feature_label_color(feature.color),
                );
                continue;
            }
            if !force_external {
                continue;
            }

            let place_top = label_rect.center().y <= self.baseline_y();
            let max_width = (self.area.width() - 8.0).max(44.0);
            let label_width = label_width.max(42.0).min(max_width);
            let left_bound = self.area.left() + 2.0;
            let right_bound = self.area.right() - 2.0;
            let mut label_left = label_rect.center().x - label_width * 0.5;
            if label_left < left_bound {
                label_left = left_bound;
            }
            if label_left + label_width > right_bound {
                label_left = (right_bound - label_width).max(left_bound);
            }
            let label_right = (label_left + label_width).min(right_bound);
            let lane = if place_top {
                Self::allocate_lane(&mut top_label_lanes, label_left, label_right, 6.0)
            } else {
                Self::allocate_lane(&mut bottom_label_lanes, label_left, label_right, 6.0)
            };
            let mut label_center_y = if place_top {
                label_rect.top() - 4.0 - external_label_lane_step * lane as f32
            } else {
                label_rect.bottom() + 4.0 + external_label_lane_step * lane as f32
            };
            let min_y = self.area.top() + external_label_height * 0.5;
            let max_y = self.area.bottom() - external_label_height * 0.5;
            label_center_y = label_center_y.clamp(min_y, max_y);
            let label_bg = Rect::from_center_size(
                Pos2::new((label_left + label_right) * 0.5, label_center_y),
                Vec2::new(
                    (label_right - label_left).max(8.0) + 6.0,
                    external_label_height,
                ),
            );
            let connector_x = label_rect
                .center()
                .x
                .clamp(label_bg.left() + 2.0, label_bg.right() - 2.0);
            let connector_from = Pos2::new(
                connector_x,
                if place_top {
                    label_rect.top()
                } else {
                    label_rect.bottom()
                },
            );
            let connector_to = Pos2::new(
                connector_x,
                if place_top {
                    label_bg.bottom()
                } else {
                    label_bg.top()
                },
            );
            painter.line_segment(
                [connector_from, connector_to],
                Stroke::new(0.8, feature.color.gamma_multiply(0.75)),
            );
            painter.rect_filled(
                label_bg,
                2.0,
                Color32::from_rgba_unmultiplied(252, 252, 252, external_bg_alpha),
            );
            painter.rect_stroke(
                label_bg,
                2.0,
                Stroke::new(1.0, feature.color.gamma_multiply(0.75)),
            );
            painter.text(
                label_bg.center(),
                Align2::CENTER_CENTER,
                label,
                FontId {
                    size: external_font_size,
                    family: FontFamily::Monospace,
                },
                Color32::from_rgb(24, 24, 24),
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
            let mut color = DnaDisplay::restriction_enzyme_group_color(key.number_of_cuts());
            let selected_here = self
                .selected_enzyme
                .as_ref()
                .map(|selected| selected.key == *key)
                .unwrap_or(false);
            if selected_here {
                color = Color32::BLACK;
            }
            let stroke_width = if selected_here { 2.0 } else { 1.0 };

            painter.line_segment(
                [Pos2::new(x, y - 8.0), Pos2::new(x, y + 8.0)],
                Stroke::new(stroke_width, color),
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
        self.draw_backbone(ui.painter(), viewport);
        self.draw_bp_ticks(ui.painter(), viewport);
        self.draw_sequence_bases(ui.painter(), viewport);
        self.draw_open_reading_frames(ui.painter(), viewport, detail);
        self.draw_features(ui.painter(), detail);
        self.draw_restriction_enzyme_sites(ui.painter(), viewport, detail);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gb_io::seq::{FeatureKind, Location};

    fn make_test_feature(location: Location) -> Feature {
        make_test_feature_with_kind("mRNA", location)
    }

    fn make_test_feature_with_kind(kind: &str, location: Location) -> Feature {
        Feature {
            kind: FeatureKind::from(kind),
            location,
            qualifiers: vec![("label".into(), Some("tp73 transcript".to_string()))],
        }
    }

    fn test_renderer_with_feature(feature: Feature, sequence_len: usize) -> RenderDnaLinear {
        let mut dna = DNAsequence::from_sequence(&"A".repeat(sequence_len)).expect("valid DNA");
        dna.features_mut().push(feature);
        let dna = Arc::new(RwLock::new(dna));
        let display = Arc::new(RwLock::new(DnaDisplay::default()));
        let mut renderer = RenderDnaLinear::new(dna, display);
        renderer.area = Rect::from_min_max(Pos2::new(0.0, 0.0), Pos2::new(1200.0, 600.0));
        renderer
    }

    fn test_renderer(sequence_len: usize) -> RenderDnaLinear {
        let dna = Arc::new(RwLock::new(
            DNAsequence::from_sequence(&"A".repeat(sequence_len)).expect("valid DNA"),
        ));
        let display = Arc::new(RwLock::new(DnaDisplay::default()));
        let mut renderer = RenderDnaLinear::new(dna, display);
        renderer.area = Rect::from_min_max(Pos2::new(0.0, 0.0), Pos2::new(1200.0, 600.0));
        renderer
    }

    #[test]
    fn multipart_join_creates_exon_rects_and_intron_connectors() {
        let feature = make_test_feature(Location::Join(vec![
            Location::simple_range(100, 160),
            Location::simple_range(220, 280),
            Location::simple_range(360, 410),
        ]));
        let mut renderer = test_renderer_with_feature(feature, 1000);
        renderer.layout_features(LinearViewport {
            start: 0,
            end: 1000,
            span: 1000,
        });
        assert_eq!(renderer.features.len(), 1);
        let fp = &renderer.features[0];
        assert_eq!(fp.exon_rects.len(), 3);
        assert_eq!(fp.intron_connectors.len(), 2);
        assert!(fp.exon_rects[0].width() > 0.0);
        assert!(fp.exon_rects[1].left() > fp.exon_rects[0].right());
    }

    #[test]
    fn complementary_multipart_feature_is_positioned_on_reverse_lane() {
        let feature = make_test_feature(Location::Complement(Box::new(Location::Join(vec![
            Location::simple_range(120, 180),
            Location::simple_range(240, 310),
        ]))));
        let mut renderer = test_renderer_with_feature(feature, 1000);
        renderer.layout_features(LinearViewport {
            start: 0,
            end: 1000,
            span: 1000,
        });
        assert_eq!(renderer.features.len(), 1);
        let fp = &renderer.features[0];
        assert!(fp.is_reverse);
        assert_eq!(fp.exon_rects.len(), 2);
        assert_eq!(fp.intron_connectors.len(), 1);
        let baseline = renderer.baseline_y();
        assert!(fp.rect.center().y > baseline);
    }

    #[test]
    fn sequence_base_track_is_enabled_at_or_below_500_bp_viewport() {
        let renderer = test_renderer(1000);
        assert!(!renderer.should_draw_sequence_bases(LinearViewport {
            start: 0,
            end: 501,
            span: 501,
        }));
        assert!(renderer.should_draw_sequence_bases(LinearViewport {
            start: 0,
            end: 500,
            span: 500,
        }));
        assert!(renderer.should_draw_sequence_bases(LinearViewport {
            start: 200,
            end: 320,
            span: 120,
        }));
    }

    #[test]
    fn sequence_base_track_respects_configured_max_span() {
        let renderer = test_renderer(1000);
        {
            let mut display = renderer.display.write().expect("display");
            display.set_linear_sequence_base_text_max_view_span_bp(120);
        }
        assert!(!renderer.should_draw_sequence_bases(LinearViewport {
            start: 0,
            end: 121,
            span: 121,
        }));
        assert!(renderer.should_draw_sequence_bases(LinearViewport {
            start: 10,
            end: 130,
            span: 120,
        }));
        {
            let mut display = renderer.display.write().expect("display");
            display.set_linear_sequence_base_text_max_view_span_bp(0);
        }
        assert!(!renderer.should_draw_sequence_bases(LinearViewport {
            start: 0,
            end: 100,
            span: 100,
        }));
    }

    #[test]
    fn sequence_base_track_can_extend_to_helical_span_limit() {
        let renderer = test_renderer(5000);
        {
            let mut display = renderer.display.write().expect("display");
            display.set_linear_sequence_base_text_max_view_span_bp(500);
            display.set_linear_sequence_helical_letters_enabled(true);
            display.set_linear_sequence_helical_max_view_span_bp(2000);
        }
        assert!(renderer.should_draw_sequence_bases(LinearViewport {
            start: 0,
            end: 1500,
            span: 1500,
        }));
        assert!(!renderer.should_draw_sequence_bases(LinearViewport {
            start: 0,
            end: 2100,
            span: 2100,
        }));
    }

    #[test]
    fn sequence_base_track_condensed_mode_uses_condensed_span_limit() {
        let renderer = test_renderer(5000);
        {
            let mut display = renderer.display.write().expect("display");
            display.set_linear_sequence_base_text_max_view_span_bp(500);
            display.set_linear_sequence_helical_letters_enabled(true);
            display.set_linear_sequence_helical_max_view_span_bp(2500);
            display.set_linear_sequence_condensed_max_view_span_bp(1200);
            display.set_linear_sequence_letter_layout_mode(
                LinearSequenceLetterLayoutMode::Condensed10Row,
            );
        }
        assert!(renderer.should_draw_sequence_bases(LinearViewport {
            start: 0,
            end: 1200,
            span: 1200,
        }));
        assert!(!renderer.should_draw_sequence_bases(LinearViewport {
            start: 0,
            end: 1300,
            span: 1300,
        }));
    }

    #[test]
    fn backbone_can_be_hidden_when_sequence_letters_are_visible() {
        let renderer = test_renderer(5000);
        let viewport = LinearViewport {
            start: 0,
            end: 500,
            span: 500,
        };
        assert!(renderer.should_draw_backbone(viewport));
        {
            let mut display = renderer.display.write().expect("display");
            display.set_linear_sequence_base_text_max_view_span_bp(500);
            display.set_linear_hide_backbone_when_sequence_bases_visible(true);
        }
        assert!(!renderer.should_draw_backbone(viewport));
        let wide_viewport = LinearViewport {
            start: 0,
            end: 3000,
            span: 3000,
        };
        assert!(renderer.should_draw_backbone(wide_viewport));
    }

    #[test]
    fn helical_projection_progress_grows_with_view_span() {
        assert_eq!(
            RenderDnaLinear::helical_projection_progress(500, 500, 2000),
            0.0
        );
        assert!((RenderDnaLinear::helical_projection_progress(1250, 500, 2000) - 0.5).abs() < 0.01);
        assert_eq!(
            RenderDnaLinear::helical_projection_progress(2000, 500, 2000),
            1.0
        );
    }

    #[test]
    fn helical_projection_approaches_two_letter_width_per_ten_bp_at_max_span() {
        assert!((RenderDnaLinear::helical_projection_x_compression(0.0) - 1.0).abs() < 0.001);
        assert!((RenderDnaLinear::helical_projection_x_compression(0.5) - 0.6).abs() < 0.001);
        assert!((RenderDnaLinear::helical_projection_x_compression(1.0) - 0.2).abs() < 0.001);
        assert!((RenderDnaLinear::helical_projection_y_phase_scale(0.0) - 0.0).abs() < 0.001);
        assert!((RenderDnaLinear::helical_projection_y_phase_scale(1.0) - 0.8).abs() < 0.001);
    }

    #[test]
    fn condensed_layout_uses_full_width_x_projection() {
        assert!(
            (RenderDnaLinear::sequence_x_compression(SequenceBaseRenderMode::Condensed10Row, 1.0)
                - 1.0)
                .abs()
                < 0.001
        );
        assert!(
            (RenderDnaLinear::sequence_x_compression(
                SequenceBaseRenderMode::ContinuousHelical,
                1.0
            ) - HELICAL_MIN_X_COMPRESSION)
                .abs()
                < 0.001
        );
    }

    #[test]
    fn helical_phase_row_mapping_applies_offset_mod_10() {
        assert_eq!(
            RenderDnaLinear::helical_phase_row_from_bottom(CONDENSED_HELICAL_ROW_COUNT, 0),
            0
        );
        assert_eq!(
            RenderDnaLinear::helical_phase_row_from_bottom(CONDENSED_HELICAL_ROW_COUNT, 3),
            3
        );
        assert_eq!(
            RenderDnaLinear::helical_phase_row_from_bottom(CONDENSED_HELICAL_ROW_COUNT, 9),
            9
        );
        assert_eq!(RenderDnaLinear::helical_phase_row_from_bottom(12, 9), 1);
    }

    #[test]
    fn helical_phase_offset_shifts_seam_without_reordering_rows() {
        let rows_offset_0: Vec<usize> = (0..20)
            .map(|bp| RenderDnaLinear::helical_phase_row_from_bottom(bp, 0))
            .collect();
        let rows_offset_3: Vec<usize> = (0..20)
            .map(|bp| RenderDnaLinear::helical_phase_row_from_bottom(bp, 3))
            .collect();
        let rows_offset_9: Vec<usize> = (0..20)
            .map(|bp| RenderDnaLinear::helical_phase_row_from_bottom(bp, 9))
            .collect();

        assert_eq!(
            rows_offset_0,
            vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        );
        assert_eq!(rows_offset_3[7], 0);
        assert_eq!(rows_offset_9[1], 0);
        for ((base, shifted_3), shifted_9) in rows_offset_0
            .iter()
            .zip(rows_offset_3.iter())
            .zip(rows_offset_9.iter())
        {
            assert_eq!(*shifted_3, (*base + 3) % CONDENSED_HELICAL_ROW_COUNT);
            assert_eq!(*shifted_9, (*base + 9) % CONDENSED_HELICAL_ROW_COUNT);
        }
    }

    #[test]
    fn condensed_helical_rows_use_non_overlapping_step_size() {
        let row_step = RenderDnaLinear::condensed_row_step(SEQUENCE_BASE_TEXT_MAX_FONT_SIZE, true);
        assert!(
            row_step
                >= SEQUENCE_BASE_TEXT_MAX_FONT_SIZE * 2.0
                    + CONDENSED_HELICAL_STRAND_GAP
                    + CONDENSED_HELICAL_ROW_PAIR_PADDING
        );
    }

    #[test]
    fn condensed_single_strand_uses_minimum_row_step() {
        let row_step = RenderDnaLinear::condensed_row_step(SEQUENCE_BASE_TEXT_MIN_FONT_SIZE, false);
        assert!(row_step >= CONDENSED_HELICAL_SINGLE_ROW_MIN_STEP);
    }

    #[test]
    fn condensed_layout_reserves_vertical_band_around_baseline() {
        let mut renderer = test_renderer(5000);
        renderer.area = Rect::from_min_max(Pos2::new(0.0, 0.0), Pos2::new(1200.0, 220.0));
        let viewport = LinearViewport {
            start: 0,
            end: 1500,
            span: 1500,
        };
        {
            let mut display = renderer.display.write().expect("display");
            display.set_linear_sequence_base_text_max_view_span_bp(500);
            display.set_linear_sequence_helical_letters_enabled(true);
            display.set_linear_sequence_helical_max_view_span_bp(2000);
            display.set_linear_sequence_letter_layout_mode(
                LinearSequenceLetterLayoutMode::Condensed10Row,
            );
        }
        renderer.layout_features(viewport);
        let half_extent = renderer
            .condensed_helical_text_half_extent(viewport)
            .expect("condensed extent");
        let baseline = renderer.baseline_y();
        let baseline_min = renderer.area.top() + BASELINE_TOP_PADDING + half_extent;
        let baseline_max = renderer.area.bottom() - BASELINE_BOTTOM_PADDING - half_extent;
        if baseline_min <= baseline_max {
            assert!(baseline >= baseline_min - 0.5);
            assert!(baseline <= baseline_max + 0.5);
        } else {
            assert!((baseline - renderer.area.center().y).abs() <= 0.5);
        }
    }

    #[test]
    fn condensed_layout_keeps_reverse_strand_ordered_below_forward_strand() {
        let font_size = SEQUENCE_BASE_TEXT_MAX_FONT_SIZE;
        let strand_half_gap = CONDENSED_HELICAL_STRAND_GAP * 0.5;
        let forward_y = -font_size - strand_half_gap;
        let reverse_y = strand_half_gap;
        let row_step = RenderDnaLinear::condensed_row_step(font_size, true);
        let row_top_index = CONDENSED_HELICAL_ROW_COUNT.saturating_sub(1);
        let row_center = row_top_index as f32 * 0.5;
        for row_from_top in 0..CONDENSED_HELICAL_ROW_COUNT {
            let centered = row_from_top as f32 - row_center;
            let phase = centered * row_step;
            let forward = forward_y + phase;
            let reverse = reverse_y + phase;
            assert!(reverse > forward);
            assert!((reverse - forward - (font_size + CONDENSED_HELICAL_STRAND_GAP)).abs() < 0.001);
        }
    }

    #[test]
    fn condensed_layout_single_strand_reserves_full_text_band() {
        let feature = make_test_feature(Location::simple_range(120, 260));
        let mut renderer = test_renderer_with_feature(feature, 5000);
        renderer.area = Rect::from_min_max(Pos2::new(0.0, 0.0), Pos2::new(1200.0, 260.0));
        let viewport = LinearViewport {
            start: 0,
            end: 1500,
            span: 1500,
        };
        {
            let mut display = renderer.display.write().expect("display");
            display.set_linear_show_double_strand_bases(false);
            display.set_linear_sequence_base_text_max_view_span_bp(500);
            display.set_linear_sequence_helical_letters_enabled(true);
            display.set_linear_sequence_condensed_max_view_span_bp(2000);
            display.set_linear_sequence_letter_layout_mode(
                LinearSequenceLetterLayoutMode::Condensed10Row,
            );
        }
        renderer.layout_features(viewport);
        let half_extent = renderer
            .condensed_helical_text_half_extent(viewport)
            .expect("condensed extent");
        let minimum_half_extent =
            CONDENSED_HELICAL_SINGLE_ROW_MIN_STEP * CONDENSED_HELICAL_ROW_COUNT as f32 * 0.5;
        assert!(half_extent >= minimum_half_extent - 0.5);
        assert_eq!(renderer.features.len(), 1);
        let clearance = renderer.baseline_y() - renderer.features[0].rect.bottom();
        assert!(clearance >= half_extent + CONDENSED_FEATURE_CLEARANCE_PADDING - 0.5);
    }

    #[test]
    fn condensed_helical_layout_replaces_backbone_line_when_bases_are_visible() {
        let renderer = test_renderer(5000);
        let viewport = LinearViewport {
            start: 0,
            end: 1500,
            span: 1500,
        };
        {
            let mut display = renderer.display.write().expect("display");
            display.set_linear_sequence_base_text_max_view_span_bp(500);
            display.set_linear_sequence_helical_letters_enabled(true);
            display.set_linear_sequence_helical_max_view_span_bp(2000);
            display.set_linear_sequence_letter_layout_mode(
                LinearSequenceLetterLayoutMode::Condensed10Row,
            );
            display.set_linear_hide_backbone_when_sequence_bases_visible(false);
        }
        assert!(renderer.should_draw_sequence_bases(viewport));
        assert!(!renderer.should_draw_backbone(viewport));
    }

    #[test]
    fn condensed_layout_pushes_feature_lanes_outward() {
        let feature = make_test_feature(Location::simple_range(120, 260));
        let mut renderer = test_renderer_with_feature(feature, 5000);
        let viewport = LinearViewport {
            start: 0,
            end: 1500,
            span: 1500,
        };
        {
            let mut display = renderer.display.write().expect("display");
            display.set_linear_sequence_base_text_max_view_span_bp(500);
            display.set_linear_sequence_helical_letters_enabled(true);
            display.set_linear_sequence_helical_max_view_span_bp(2000);
            display.set_linear_sequence_letter_layout_mode(
                LinearSequenceLetterLayoutMode::ContinuousHelical,
            );
        }
        renderer.layout_features(viewport);
        assert_eq!(renderer.features.len(), 1);
        let continuous_baseline = renderer.baseline_y();
        let continuous_clearance = continuous_baseline - renderer.features[0].rect.bottom();

        {
            let mut display = renderer.display.write().expect("display");
            display.set_linear_sequence_letter_layout_mode(
                LinearSequenceLetterLayoutMode::Condensed10Row,
            );
        }
        renderer.layout_features(viewport);
        assert_eq!(renderer.features.len(), 1);
        let condensed_baseline = renderer.baseline_y();
        let condensed_clearance = condensed_baseline - renderer.features[0].rect.bottom();

        assert!(condensed_clearance > continuous_clearance);
        assert!(condensed_clearance >= 20.0);
    }

    #[test]
    fn condensed_layout_pushes_reverse_feature_lanes_outward() {
        let feature = make_test_feature(Location::Complement(Box::new(Location::simple_range(
            120, 260,
        ))));
        let mut renderer = test_renderer_with_feature(feature, 5000);
        let viewport = LinearViewport {
            start: 0,
            end: 1500,
            span: 1500,
        };
        {
            let mut display = renderer.display.write().expect("display");
            display.set_linear_sequence_base_text_max_view_span_bp(500);
            display.set_linear_sequence_helical_letters_enabled(true);
            display.set_linear_sequence_helical_max_view_span_bp(2000);
            display.set_linear_sequence_letter_layout_mode(
                LinearSequenceLetterLayoutMode::ContinuousHelical,
            );
        }
        renderer.layout_features(viewport);
        assert_eq!(renderer.features.len(), 1);
        assert!(renderer.features[0].is_reverse);
        let continuous_baseline = renderer.baseline_y();
        let continuous_clearance = renderer.features[0].rect.top() - continuous_baseline;

        {
            let mut display = renderer.display.write().expect("display");
            display.set_linear_sequence_letter_layout_mode(
                LinearSequenceLetterLayoutMode::Condensed10Row,
            );
        }
        renderer.layout_features(viewport);
        assert_eq!(renderer.features.len(), 1);
        let condensed_baseline = renderer.baseline_y();
        let condensed_clearance = renderer.features[0].rect.top() - condensed_baseline;

        assert!(condensed_clearance > continuous_clearance);
        assert!(condensed_clearance >= 20.0);
    }

    #[test]
    fn condensed_layout_pushes_regulatory_near_baseline_lane_outward() {
        let feature =
            make_test_feature_with_kind("regulatory_region", Location::simple_range(120, 260));
        let mut renderer = test_renderer_with_feature(feature, 5000);
        let viewport = LinearViewport {
            start: 0,
            end: 1500,
            span: 1500,
        };
        {
            let mut display = renderer.display.write().expect("display");
            display.set_regulatory_tracks_near_baseline(true);
            display.set_linear_sequence_base_text_max_view_span_bp(500);
            display.set_linear_sequence_helical_letters_enabled(true);
            display.set_linear_sequence_helical_max_view_span_bp(2000);
            display.set_linear_sequence_letter_layout_mode(
                LinearSequenceLetterLayoutMode::ContinuousHelical,
            );
        }
        renderer.layout_features(viewport);
        assert_eq!(renderer.features.len(), 1);
        assert!(renderer.features[0].is_regulatory);
        let continuous_baseline = renderer.baseline_y();
        let continuous_clearance = continuous_baseline - renderer.features[0].rect.bottom();

        {
            let mut display = renderer.display.write().expect("display");
            display.set_linear_sequence_letter_layout_mode(
                LinearSequenceLetterLayoutMode::Condensed10Row,
            );
        }
        renderer.layout_features(viewport);
        assert_eq!(renderer.features.len(), 1);
        let condensed_baseline = renderer.baseline_y();
        let condensed_clearance = condensed_baseline - renderer.features[0].rect.bottom();

        assert!(condensed_clearance > continuous_clearance);
        assert!(condensed_clearance >= 12.0);
    }
}
