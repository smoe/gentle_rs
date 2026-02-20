use crate::{
    dna_display::{DnaDisplay, Selection, TfbsDisplayCriteria, VcfDisplayCriteria},
    dna_sequence::DNAsequence,
    feature_location::{feature_ranges_sorted_i64, normalize_range, unwrap_ranges_monotonic},
    gc_contents::GcRegion,
    render_dna::RenderDna,
    render_dna::RestrictionEnzymePosition,
    restriction_enzyme::RestrictionEnzymeKey,
};
use eframe::egui::{
    self, Align2, Color32, FontFamily, FontId, PointerState, Pos2, Rect, Shape, Stroke,
};
use gb_io::seq::Feature;
use lazy_static::lazy_static;
use std::{
    collections::{BTreeSet, HashMap},
    sync::{Arc, RwLock},
};

// Defines static stroke styles (BLACK_1, GRAY_1) and a color map (ORF_COLORS) for Open Reading Frames (ORFs).
lazy_static! {
    pub static ref BLACK_1: Stroke = Stroke {
        width: 1.0,
        color: Color32::BLACK,
    };
    pub static ref GRAY_1: Stroke = Stroke {
        width: 1.0,
        color: Color32::GRAY,
    };
    static ref ORF_COLORS: HashMap<i32, Color32> = {
        let mut m = HashMap::new();
        m.insert(-1, Color32::LIGHT_RED);
        m.insert(-2, Color32::LIGHT_GREEN);
        m.insert(-3, Color32::LIGHT_BLUE);
        m.insert(1, Color32::DARK_RED);
        m.insert(2, Color32::DARK_GREEN);
        m.insert(3, Color32::DARK_BLUE);
        m
    };
}

/// Represents the position and attributes of a feature (e.g., gene, CDS) on the circular DNA.
#[derive(Debug, Clone)]
struct FeatureSegmentPosition {
    // FIXME: Clarify if the first position is 0 or 1
    from: i64,
    // FIXME: Clarify if the first position is 0 or 1
    to: i64,
    to_90: i64,
    angle_start: f32,
    angle_stop: f32,
}

/// Represents the position and attributes of a feature (e.g., gene, CDS) on the circular DNA.
#[derive(Debug, Clone)]
struct FeaturePosition {
    feature_number: usize,
    // FIXME: Clarify if the first position is 0 or 1
    from: i64,
    // FIXME: Clarify if the first position is 0 or 1
    to: i64,
    /// Representation as segment - start (degree: 0-360)
    angle_start: f32,
    /// Representation as segment - end (degree: 0-360)
    angle_stop: f32,
    /// Minimal distance to center of circular feature
    inner: f32,
    /// Maximal distance to center of circular feature
    outer: f32,
    to_90: i64,
    /// Directed features shall end with an arrow-like indication of the direction.
    is_pointy: bool,
    color: Color32,
    band: f32,
    label: String,
    segments: Vec<FeatureSegmentPosition>,
    intron_arches: Vec<(i64, i64)>,
    intron_arch_color: Color32,
    intron_arch_width: f32,
    intron_arch_lift_factor: f32,
}

impl FeaturePosition {
    /// Checks if a given angle and distance fall within the bounds of the feature.
    fn contains_angle_distance(&self, angle: f32, distance: f32) -> bool {
        if self.inner > distance || self.outer < distance {
            return false;
        }
        self.segments.iter().any(|segment| {
            if segment.angle_stop < segment.angle_start {
                // Feature extends over zero point
                (angle >= segment.angle_start && angle <= 360.0)
                    || (angle >= 0.0 && angle <= segment.angle_stop)
            } else {
                angle >= segment.angle_start && angle <= segment.angle_stop
            }
        })
    }
}

/// Manages the rendering of circular DNA, including features, enzyme sites, and user interactions.
#[derive(Debug, Clone)]
pub struct RenderDnaCircular {
    dna: Arc<RwLock<DNAsequence>>,
    display: Arc<RwLock<DnaDisplay>>,
    sequence_length: i64,
    area: Rect,
    center: Pos2,
    radius: f32,
    features: Vec<FeaturePosition>,
    restriction_enzyme_sites: Vec<RestrictionEnzymePosition>,
    selected_feature_number: Option<usize>,
    hovered_feature_number: Option<usize>,
    hover_enzyme: Option<RestrictionEnzymePosition>,
}

impl RenderDnaCircular {
    /// Initializes a new RenderDnaCircular instance.
    pub fn new(dna: Arc<RwLock<DNAsequence>>, display: Arc<RwLock<DnaDisplay>>) -> Self {
        Self {
            dna,
            sequence_length: 0,
            display,
            area: Rect::NOTHING,
            center: Pos2::ZERO,
            radius: 0.0,
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

    /// Handles click events to select features.
    pub fn on_click(&mut self, pointer_state: PointerState) {
        if let Some(pos) = pointer_state.latest_pos() {
            self.selected_feature_number = self.get_clicked_feature(pos).map(|f| f.feature_number);
        } else {
            println!("I RenderDnaCircular::on_click: Could not select any feature.");
        }
    }

    /// Handles hover events to highlight features and enzyme sites.
    pub fn on_hover(&mut self, pointer_state: PointerState) {
        if let Some(pos) = pointer_state.latest_pos() {
            self.hovered_feature_number = self.get_clicked_feature(pos).map(|f| f.feature_number);
            self.hover_enzyme = self.get_re_site_for_positon(pos);
        } else {
            println!("I RenderDnaCircular::on_hover: Could not select any feature.");
        }
    }

    /// Handles double-click events to select features or enzyme sites.
    pub fn on_double_click(&mut self, pointer_state: PointerState) {
        self.display.write().unwrap().deselect();
        if let Some(pos) = pointer_state.latest_pos() {
            if let Some(feature) = self.get_clicked_feature(pos) {
                // println!("Double-clicked {:?}", feature);
                let selection = Selection::new(
                    feature.from as usize,
                    feature.to as usize,
                    self.sequence_length as usize,
                );
                self.display.write().unwrap().select(selection);
            } else if let Some(re_pos) = self.get_re_site_for_positon(pos) {
                println!("Double-clicked {re_pos:?}");
                let selection = Selection::new(
                    re_pos.key().from() as usize,
                    re_pos.key().to() as usize,
                    self.sequence_length as usize,
                );
                self.display.write().unwrap().select(selection);
            } else {
                println!("I RenderDnaCircular::on_double_click: Could not select any feature or RestrictionEnzyme.");
            }
        }
    }

    /// Finds the restriction enzyme site at a given position.
    fn get_re_site_for_positon(&self, pos: Pos2) -> Option<RestrictionEnzymePosition> {
        self.restriction_enzyme_sites
            .iter()
            .find(|rep| rep.area.contains(pos))
            .cloned()
    }

    /// Finds the feature at a given position.
    fn get_clicked_feature(&self, pos: Pos2) -> Option<&FeaturePosition> {
        let (angle, distance) = self.get_angle_distance(pos);
        let angle = Self::normalize_angle(angle - 90.0);
        let clicked_features = self
            .features
            .iter()
            .filter(|feature| feature.contains_angle_distance(angle, distance))
            .collect::<Vec<_>>();
        clicked_features.first().map(|f| f.to_owned())
    }

    pub fn set_area(&mut self, area: Rect) {
        self.area = area;
        self.center = self.area.center();
    }

    pub fn selected_feature_number(&self) -> Option<usize> {
        self.selected_feature_number.to_owned()
    }

    pub fn hovered_feature_number(&self) -> Option<usize> {
        self.hovered_feature_number
    }

    pub fn select_feature(&mut self, feature_number: Option<usize>) {
        self.selected_feature_number = feature_number;
    }

    fn layout_needs_recomputing(&self) -> bool {
        self.display
            .read()
            .map(|d| d.update_layout().update_map_dna())
            .unwrap_or(false)
    }

    fn layout_was_updated(&self) {
        self.display
            .write()
            .unwrap()
            .update_layout_mut()
            .map_dna_updated();
    }

    /// Renders the circular DNA visualization
    pub fn render(&mut self, ui: &mut egui::Ui, area: Rect) {
        let area_changed = self.area != area;
        self.area = area;
        self.radius = self.area.width().min(self.area.height()) * 0.35;
        self.center = self.area.center();
        self.sequence_length = self.dna.read().expect("DNA lock poisoned").len() as i64;
        if !Self::is_rect_usable(self.area) {
            return;
        }

        if self.sequence_length <= 0 {
            let painter = ui.painter();
            self.draw_main_label(painter);
            self.draw_bp(painter);
            return;
        }

        if (area_changed || self.layout_needs_recomputing()) && Self::is_rect_usable(self.area) {
            self.layout_features();
            self.layout_was_updated();
        }

        let painter = ui.painter();
        self.draw_selection(painter);
        self.draw_backbone(painter);
        self.draw_gc_contents(painter);
        self.draw_methylation_sites(painter);
        self.draw_main_label(painter);
        self.draw_bp(painter);
        self.draw_open_reading_frames(painter);
        self.draw_restriction_enzyme_sites(painter);
        self.draw_features(painter);
        self.draw_hovered_feature(painter);
    }

    /// Draws a filled circular section
    fn draw_circle_section(&self, start: u64, end: u64, color: Color32, painter: &egui::Painter) {
        let center = self.area.center();
        let radius = self.radius;

        // Generate points to create the filled section
        let num_points = 500; // More points = smoother curve
        let points: Vec<Pos2> = (0..=num_points)
            .map(|i| {
                let pos = start + (end - start) * i / num_points;
                self.pos2xy(pos as i64, radius)
                // let t = i as f32 / num_points as f32;
                // let angle = start_angle + t * (end_angle - start_angle);
                // let x = center.x + radius * angle.cos();
                // let y = center.y + radius * angle.sin();
                // Pos2::new(x, y)
            })
            .collect();

        // Create a vector of points including the center to form a filled shape
        let mut filled_points = vec![center];
        filled_points.extend(points);

        // Draw the filled section
        let fill_color = color.to_owned();
        let stroke = Stroke::new(2.0, color);

        let shape = Shape::convex_polygon(filled_points, fill_color, stroke);

        painter.add(shape);
    }

    /// Draws the selected region on the circular DNA
    fn draw_selection(&self, painter: &egui::Painter) {
        let selection = match self.display.read().unwrap().selection() {
            Some(selection) => selection,
            None => return,
        };
        let parts = selection.parts();
        for part in parts {
            self.draw_circle_section(part.0 as u64, part.1 as u64, Color32::LIGHT_GRAY, painter);
        }
    }

    /// Displays information about the hovered feature
    fn draw_hovered_feature(&self, painter: &egui::Painter) {
        let Some(feature_id) = self.hovered_feature_number else {
            return;
        };
        let Some(fp) = self
            .features
            .iter()
            .find(|feature| feature.feature_number == feature_id)
        else {
            return;
        };
        let feature = self
            .dna
            .read()
            .unwrap()
            .features()
            .get(fp.feature_number)
            .cloned();
        if let Some(feature) = feature {
            if let Ok((from, to)) = feature.location.find_bounds() {
                let text = format!("{}: {}-{}", &fp.label, from, to);
                let font = FontId {
                    size: 12.0,
                    family: FontFamily::Monospace,
                };
                painter.text(
                    self.area.left_bottom(),
                    Align2::LEFT_BOTTOM,
                    text,
                    font.to_owned(),
                    Color32::DARK_GRAY,
                );
            }
        }
    }

    /// Converts a position to an angle and geometric distance (Pythagoras) from the center.
    fn get_angle_distance(&self, pos: Pos2) -> (f32, f32) {
        let diff_x = pos.x - self.center.x;
        let diff_y = pos.y - self.center.y;
        let angle = diff_y.atan2(diff_x) * 180.0 / std::f32::consts::PI + 90.0;
        let angle = Self::normalize_angle(angle);
        let distance = (diff_x.powi(2) + diff_y.powi(2)).sqrt();
        (angle, distance)
    }

    /// Draws an arc with an arrow indicating direction
    fn draw_pointed_arc(
        &self,
        from: i32,
        to: i32,
        radius: f32,
        is_reverse: bool,
        stroke: Stroke,
        painter: &egui::Painter,
    ) {
        // start <= end
        let start = from.min(to);
        let end = from.max(to);
        if is_reverse {
            let step = -10;
            let mut pos = end;
            let mut last_pos = pos;

            // Draw starting point
            let r0 = radius / 75.0;
            let point = self.pos2xy(pos as i64, radius);
            painter.circle_filled(point, r0, stroke.color.to_owned());

            // Draw arc
            while pos > start {
                let point1 = self.pos2xy(last_pos as i64, radius);
                let point2 = self.pos2xy(pos as i64, radius);
                painter.line_segment([point1, point2], stroke.to_owned());
                last_pos = pos;
                pos += step;
            }

            // Draw arrow
            last_pos = start + step * 2;
            let point1 = self.pos2xy(last_pos as i64, radius * 0.98);
            let point2 = self.pos2xy(start as i64, radius);
            painter.line_segment([point1, point2], stroke.to_owned());
            let point1 = self.pos2xy(last_pos as i64, radius * 1.02);
            let point2 = self.pos2xy(start as i64, radius);
            painter.line_segment([point1, point2], stroke.to_owned());
        } else {
            let step = 10;
            let mut pos = start;
            let mut last_pos = pos;

            // Draw starting point
            let r0 = radius / 75.0;
            let point = self.pos2xy(pos as i64, radius);
            painter.circle_filled(point, r0, stroke.color.to_owned());

            // Draw arc
            while pos < end {
                let point1 = self.pos2xy(last_pos as i64, radius);
                let point2 = self.pos2xy(pos as i64, radius);
                painter.line_segment([point1, point2], stroke.to_owned());
                last_pos = pos;
                pos += step;
            }

            // Draw arrow
            last_pos = end - step * 2;
            let point1 = self.pos2xy(last_pos as i64, radius * 0.98);
            let point2 = self.pos2xy(end as i64, radius);
            painter.line_segment([point1, point2], stroke.to_owned());
            let point1 = self.pos2xy(last_pos as i64, radius * 1.02);
            let point2 = self.pos2xy(end as i64, radius);
            painter.line_segment([point1, point2], stroke.to_owned());
        };
    }

    /// Draws Open Reading Frames (ORFs) on the circular DNA
    fn draw_open_reading_frames(&self, painter: &egui::Painter) {
        if !self
            .display
            .read()
            .unwrap()
            .show_open_reading_frames_effective()
        {
            return;
        }
        let orfs = self.dna.read().unwrap().open_reading_frames().to_owned();
        for orf in orfs {
            let color = match ORF_COLORS.get(&orf.frame()) {
                Some(color) => color.to_owned(),
                None => continue,
            };
            let radius = self.radius * 1.1 + self.radius * 0.05 * (orf.frame() as f32);
            let stroke = Stroke::new(1.0, color);
            self.draw_pointed_arc(
                orf.from(),
                orf.to(),
                radius,
                orf.is_reverse(),
                stroke,
                painter,
            );
        }
    }

    /// Draws methylation sites on the circular DNA
    fn draw_methylation_sites(&self, painter: &egui::Painter) {
        if !self.display.read().unwrap().show_methylation_sites() {
            return;
        }
        let radius_lower = self.radius * 0.97;
        let stroke = Stroke::new(1.0, Color32::DARK_RED);
        let methylation_sites = self.dna.read().unwrap().methylation_sites().to_owned();
        for site in methylation_sites.sites() {
            let point1 = self.pos2xy(*site as i64, self.radius);
            let point2 = self.pos2xy(*site as i64, radius_lower);
            painter.line_segment([point1, point2], stroke.to_owned());
        }
    }

    /// Draws GC content regions on the circular DNA.
    fn draw_gc_contents(&self, painter: &egui::Painter) {
        if !self.display.read().unwrap().show_gc_contents() {
            return;
        }
        let radius = self.radius * 2.0 / 3.0;
        let mut last_point = self.pos2xy(0, radius);
        let gc_content = self.dna.read().unwrap().gc_content().to_owned();
        for gc_region in gc_content.regions() {
            last_point = self.draw_gc_arc(gc_region, radius, painter, last_point);
        }
    }

    /// Draws a GC content arc segment
    fn draw_gc_arc(
        &self,
        gc_region: &GcRegion,
        radius: f32,
        painter: &egui::Painter,
        last_point: Pos2,
    ) -> Pos2 {
        let point = self.pos2xy(gc_region.to() as i64, radius);
        let color = Color32::from_rgb(
            255 - (gc_region.gc() * 255.0) as u8,
            (gc_region.gc() * 255.0) as u8,
            0,
        );
        let stroke = Stroke::new(10.0, color);
        painter.line_segment([last_point, point], stroke);
        point
    }

    /// Draws the backbone of the circular DNA with tick marks
    fn draw_backbone(&mut self, painter: &egui::Painter) {
        painter.circle_stroke(self.center.to_owned(), self.radius, BLACK_1.to_owned());
        let mut tick: i64 = 1;
        while tick * 10 < self.sequence_length {
            tick *= 10;
        }

        let font_tick = FontId {
            size: 9.0,
            family: FontFamily::Monospace,
        };

        let mut pos = tick; // Skip 0 point
        while pos < self.sequence_length {
            let p1 = self.pos2xy(pos, self.radius);
            let p2 = self.pos2xy(pos, self.radius * 0.87);
            let p3 = self.pos2xy(pos, self.radius * 0.85);
            painter.line_segment([p1, p2], BLACK_1.to_owned());

            let align = if pos > self.sequence_length / 2 {
                Align2::LEFT_CENTER
            } else {
                Align2::RIGHT_CENTER
            };
            painter.text(
                p3,
                align,
                format!("{pos}"),
                font_tick.to_owned(),
                Color32::BLACK,
            );

            pos += tick;
        }
    }

    /// Lays out features on the circular DNA
    fn layout_features(&mut self) {
        self.features.clear();
        let (
            show_cds_features,
            show_gene_features,
            show_mrna_features,
            show_tfbs,
            tfbs_display_criteria,
            vcf_display_criteria,
            hidden_feature_kinds,
        ) = {
            let display = self.display.read().expect("Display lock poisoned");
            (
                display.show_cds_features(),
                display.show_gene_features(),
                display.show_mrna_features(),
                display.show_tfbs(),
                display.tfbs_display_criteria(),
                display.vcf_display_criteria(),
                display.hidden_feature_kinds().clone(),
            )
        };
        let features = self
            .dna
            .read()
            .expect("DNA lock poisoned")
            .features()
            .to_owned();
        for (feature_number, feature) in features.iter().enumerate() {
            let fp_opt = self.layout_feature_from_location(
                feature,
                show_cds_features,
                show_gene_features,
                show_mrna_features,
                show_tfbs,
                tfbs_display_criteria,
                vcf_display_criteria.clone(),
                &hidden_feature_kinds,
            );
            if let Some(mut fp) = fp_opt {
                fp.feature_number = feature_number;
                self.features.push(fp);
            }
        }
        self.compact_feature_bands();
    }

    fn normalized_feature_range(&self, from: i64, to: i64) -> Option<(i64, i64)> {
        normalize_range(self.sequence_length, from, to)
    }

    fn compact_feature_bands(&mut self) {
        if self.features.is_empty() || self.sequence_length <= 0 {
            return;
        }

        let thickness = self.feature_thickness().max(1.0);
        let lane_gap = (thickness * 0.35).max(1.0);
        let base_offset = (thickness * 0.15).max(0.5);
        let overlap_padding_bp = ((self.sequence_length as f32) * 0.002).ceil() as i64;

        #[derive(Clone)]
        struct Seed {
            feature_index: usize,
            start: i64,
            end: i64,
            span: i64,
        }

        let mut seeds: Vec<Seed> = self
            .features
            .iter()
            .enumerate()
            .filter_map(|(idx, feature)| {
                let (start, end) = self.normalized_feature_range(feature.from, feature.to)?;
                Some(Seed {
                    feature_index: idx,
                    start,
                    end,
                    span: end.saturating_sub(start).max(1),
                })
            })
            .collect();

        seeds.sort_by(|a, b| b.span.cmp(&a.span).then_with(|| a.start.cmp(&b.start)));

        let mut lane_ends: Vec<Vec<(i64, i64)>> = Vec::new();
        let mut lane_for_feature: Vec<usize> = vec![0; self.features.len()];

        'seed_loop: for seed in seeds {
            for (lane_idx, lane_ranges) in lane_ends.iter_mut().enumerate() {
                let collides = lane_ranges.iter().any(|(other_start, other_end)| {
                    !(seed.end + overlap_padding_bp < *other_start
                        || seed.start > *other_end + overlap_padding_bp)
                });
                if !collides {
                    lane_ranges.push((seed.start, seed.end));
                    lane_for_feature[seed.feature_index] = lane_idx;
                    continue 'seed_loop;
                }
            }
            lane_ends.push(vec![(seed.start, seed.end)]);
            lane_for_feature[seed.feature_index] = lane_ends.len() - 1;
        }

        for (feature_idx, feature) in self.features.iter_mut().enumerate() {
            let lane = lane_for_feature.get(feature_idx).copied().unwrap_or(0) as f32;
            let inner = self.radius + base_offset + lane * (thickness + lane_gap);
            feature.inner = inner;
            feature.outer = inner + thickness;
            feature.band = lane + 1.0;
        }
    }

    fn feature_thickness(&self) -> f32 {
        self.radius / 20.0
    }

    fn intron_arch_style(feature: &Feature, base_color: Color32) -> (Color32, f32, f32) {
        if RenderDna::is_mrna_feature(feature) {
            // mRNA gets stronger intron arches because exon/intron structure is central here.
            return (Color32::from_rgb(120, 60, 20), 1.8, 1.05);
        }
        if RenderDna::is_cds_feature(feature) || RenderDna::is_gene_feature(feature) {
            return (base_color.gamma_multiply(0.8), 1.25, 0.75);
        }
        if RenderDna::is_regulatory_feature(feature) {
            return (base_color.gamma_multiply(0.85), 1.1, 0.65);
        }
        (Color32::DARK_GRAY, 1.0, 0.6)
    }

    fn layout_feature_from_location(
        &self,
        feature: &Feature,
        show_cds_features: bool,
        show_gene_features: bool,
        show_mrna_features: bool,
        show_tfbs: bool,
        tfbs_display_criteria: TfbsDisplayCriteria,
        vcf_display_criteria: VcfDisplayCriteria,
        hidden_feature_kinds: &BTreeSet<String>,
    ) -> Option<FeaturePosition> {
        if !Self::draw_feature(
            feature,
            show_cds_features,
            show_gene_features,
            show_mrna_features,
            show_tfbs,
            tfbs_display_criteria,
            vcf_display_criteria,
            hidden_feature_kinds,
        ) {
            return None;
        }
        let seq_len = self.sequence_length;
        if seq_len <= 0 {
            return None;
        }
        let location_ranges = feature_ranges_sorted_i64(feature);
        if location_ranges.is_empty() {
            return None;
        }
        let ranges = unwrap_ranges_monotonic(seq_len, &location_ranges);
        if ranges.is_empty() {
            return None;
        }

        let feature_from = ranges
            .iter()
            .map(|(start, _)| start.rem_euclid(seq_len))
            .min()
            .unwrap_or(0);
        let feature_to = ranges
            .iter()
            .map(|(_, end)| end.rem_euclid(seq_len))
            .max()
            .unwrap_or(feature_from);

        let feature_color = RenderDna::feature_color(feature);
        let (intron_arch_color, intron_arch_width, intron_arch_lift_factor) =
            Self::intron_arch_style(feature, feature_color);
        let mut ret: FeaturePosition = FeaturePosition {
            feature_number: 0,
            from: feature_from,
            to: feature_to,
            angle_start: 0.0,
            angle_stop: 0.0,
            inner: 0.0,
            outer: 0.0,
            to_90: 0,
            is_pointy: RenderDna::is_feature_pointy(feature),
            color: feature_color,
            band: 0.0,
            label: RenderDna::feature_name(feature),
            segments: vec![],
            intron_arches: vec![],
            intron_arch_color,
            intron_arch_width,
            intron_arch_lift_factor,
        };
        // Actual radial packing is computed later for currently visible features.
        ret.inner = self.radius - self.feature_thickness() / 2.0;
        ret.outer = self.radius + self.feature_thickness() / 2.0;
        if ret.inner > ret.outer {
            std::mem::swap(&mut ret.inner, &mut ret.outer);
        }
        ret.segments = ranges
            .iter()
            .enumerate()
            .map(|(idx, (range_start, range_end))| {
                let is_last_segment = idx + 1 == ranges.len();
                let to_90 = if ret.is_pointy && is_last_segment {
                    range_end - (range_end - range_start) / 20
                } else {
                    *range_end
                };
                FeatureSegmentPosition {
                    from: *range_start,
                    to: *range_end,
                    to_90,
                    angle_start: self.angle(*range_start),
                    angle_stop: self.angle(to_90),
                }
            })
            .collect();
        ret.to_90 = ret.segments.last().map(|s| s.to_90).unwrap_or(ret.to);
        ret.angle_start = ret
            .segments
            .first()
            .map(|s| s.angle_start)
            .unwrap_or_default();
        ret.angle_stop = ret
            .segments
            .last()
            .map(|s| s.angle_stop)
            .unwrap_or_default();
        ret.intron_arches = ranges
            .windows(2)
            .filter_map(|pair| {
                let (_left_start, left_end) = pair[0];
                let (right_start, _right_end) = pair[1];
                if right_start > left_end + 1 {
                    Some((left_end, right_start))
                } else {
                    None
                }
            })
            .collect();
        Some(ret)
    }

    fn draw_features(&mut self, painter: &egui::Painter) {
        if !self
            .display
            .read()
            .expect("Display lock poisoned")
            .show_features()
        {
            return;
        }
        for feature in &self.features {
            self.draw_feature_from_range(painter, feature);
        }
    }

    /// Draws a feature from a given range
    fn draw_feature_from_range(&self, painter: &egui::Painter, ret: &FeaturePosition) {
        let selected = self.selected_feature_number == Some(ret.feature_number);
        let hovered = self.hovered_feature_number == Some(ret.feature_number);
        let feature_stroke_width = if selected {
            1.7
        } else if hovered {
            1.3
        } else {
            1.0
        };

        for (segment_idx, segment) in ret.segments.iter().enumerate() {
            let mut feature_points: Vec<Pos2> = vec![];
            feature_points.push(self.pos2xy(segment.from, ret.outer));
            feature_points.push(self.pos2xy(segment.from, ret.inner));

            let points = self.generate_arc(ret.inner, segment.angle_start, segment.angle_stop);
            feature_points.extend(points);

            let is_last_segment = segment_idx + 1 == ret.segments.len();
            if ret.is_pointy && is_last_segment {
                feature_points.push(self.pos2xy(segment.to, (ret.outer + ret.inner) / 2.0));
            }

            feature_points.push(self.pos2xy(segment.to_90, ret.outer));

            let points = self.generate_arc(ret.outer, segment.angle_stop, segment.angle_start);
            feature_points.extend(points);

            let stroke = Stroke {
                width: feature_stroke_width,
                color: ret.color,
            };
            let line = Shape::closed_line(feature_points, stroke);
            painter.add(line);
        }

        let mut connector_radius =
            ret.outer + self.feature_thickness() * ret.intron_arch_lift_factor;
        if selected {
            connector_radius += self.feature_thickness() * 0.15;
        } else if hovered {
            connector_radius += self.feature_thickness() * 0.08;
        }
        let mut connector_color = ret.intron_arch_color;
        let mut connector_width = ret.intron_arch_width;
        if hovered {
            connector_color = connector_color.gamma_multiply(1.15);
            connector_width += 0.2;
        }
        if selected {
            connector_color = Color32::YELLOW;
            connector_width += 0.5;
        }
        for (from, to) in &ret.intron_arches {
            self.draw_intron_arch(
                painter,
                *from,
                *to,
                ret.outer,
                connector_radius,
                connector_color,
                connector_width,
            );
        }

        let font_feature = FontId {
            size: 10.0,
            family: FontFamily::Monospace,
        };

        // Draw feature label
        if !ret.label.trim().is_empty() {
            if let Some(widest_segment) = ret
                .segments
                .iter()
                .max_by_key(|segment| segment.to.saturating_sub(segment.from))
            {
                let label_radius = ret.outer + self.feature_thickness() * 0.35;
                let arc_px = (widest_segment.to.saturating_sub(widest_segment.from) as f32
                    / self.sequence_length.max(1) as f32)
                    * std::f32::consts::TAU
                    * label_radius.max(1.0);
                let label_px = (ret.label.chars().count().max(1) as f32) * 6.0;
                if arc_px >= label_px {
                    let middle = (widest_segment.to + widest_segment.from) / 2;
                    let middle_norm = middle.rem_euclid(self.sequence_length.max(1));
                    let point = self.pos2xy(middle, label_radius);
                    let align = if middle_norm > self.sequence_length / 2 {
                        Align2::RIGHT_CENTER
                    } else {
                        Align2::LEFT_CENTER
                    };
                    painter.text(point, align, ret.label.to_owned(), font_feature, ret.color);
                }
            }
        }
    }

    fn draw_intron_arch(
        &self,
        painter: &egui::Painter,
        from: i64,
        to: i64,
        base_radius: f32,
        arch_radius: f32,
        color: Color32,
        stroke_width: f32,
    ) {
        if self.sequence_length <= 0 || to <= from {
            return;
        }
        let stroke = Stroke::new(stroke_width, color);
        let start_base = self.pos2xy(from, base_radius);
        let start_arch = self.pos2xy(from, arch_radius);
        let end_base = self.pos2xy(to, base_radius);
        let end_arch = self.pos2xy(to, arch_radius);
        painter.line_segment([start_base, start_arch], stroke.to_owned());
        painter.line_segment([end_base, end_arch], stroke.to_owned());

        let arch_points = self.generate_arc_bp(arch_radius, from, to);
        if arch_points.len() >= 2 {
            painter.add(Shape::line(arch_points, stroke));
        }
    }

    /// Generates points for drawing an arc
    fn generate_arc(&self, radius: f32, angle_start: f32, angle_stop: f32) -> Vec<Pos2> {
        if angle_start == 0.0 && angle_stop == 360.0 {
            return vec![];
        }
        let mut points = vec![];
        let n = 100;
        for i in 0..n {
            let angle = angle_start + (angle_stop - angle_start) * (i as f32 / n as f32);
            let angle = angle * std::f32::consts::PI / 180.0;
            let x = self.center.x + radius * angle.cos();
            let y = self.center.y + radius * angle.sin();
            points.push(Pos2::new(x, y));
        }
        points
    }

    fn generate_arc_bp(&self, radius: f32, start_bp: i64, end_bp: i64) -> Vec<Pos2> {
        if self.sequence_length <= 0 || end_bp <= start_bp {
            return vec![];
        }
        let span = end_bp - start_bp;
        let mut n = ((span as f32 / self.sequence_length as f32) * 180.0).ceil() as i64;
        n = n.clamp(8, 180);

        let mut points = Vec::with_capacity((n + 1) as usize);
        for i in 0..=n {
            let pos = start_bp + (span * i) / n;
            points.push(self.pos2xy(pos, radius));
        }
        points
    }

    /// Draws the main label (name) of the DNA sequence
    fn draw_main_label(&self, painter: &egui::Painter) {
        let font_label = FontId {
            size: 20.0,
            family: FontFamily::Proportional,
        };

        let name = match self
            .dna
            .read()
            .expect("DNA lock poisoned")
            .name()
            .as_ref()
            .map(|s| s.to_owned())
        {
            Some(name) => name,
            None => "<no name>".to_string(),
        };

        painter.text(
            self.center.to_owned(),
            Align2::CENTER_BOTTOM,
            name,
            font_label,
            Color32::BLACK,
        );
    }

    /// Draws restriction enzyme sites on the circular DNA.
    fn draw_restriction_enzyme_sites(&mut self, painter: &egui::Painter) {
        self.restriction_enzyme_sites.clear();
        if !self.display.read().unwrap().show_restriction_enzyme_sites() {
            return;
        }
        let font_tick = FontId {
            size: 9.0,
            family: FontFamily::Proportional,
        };

        let mut re_positions: Vec<RestrictionEnzymeKey> = self
            .dna
            .read()
            .unwrap()
            .restriction_enzyme_groups()
            .keys()
            .cloned()
            .collect();
        re_positions.sort();
        let mut last_rect = Rect::NOTHING;
        for restriction_enzyme_key in re_positions {
            let pos = restriction_enzyme_key.pos() as i64;
            let label = self
                .dna
                .read()
                .unwrap()
                .restriction_enzyme_groups()
                .get(&restriction_enzyme_key)
                .unwrap()
                .join(", ");
            let label = if pos < self.sequence_length / 2 {
                format!("{pos} {label}")
            } else {
                format!("{label} {pos}")
            };
            let cuts = restriction_enzyme_key.number_of_cuts();
            let font_color = DnaDisplay::restriction_enzyme_group_color(cuts);

            let p1 = self.pos2xy(pos, self.radius);
            let p2 = self.pos2xy(pos, self.radius * 1.15);
            let mut p3 = self.pos2xy(pos, self.radius * 1.25);
            p3.y = p2.y;
            if pos < self.sequence_length / 2 {
                while p3.y < last_rect.bottom() + 3.0 {
                    p3.y += 1.0;
                }
            } else {
                while p3.y > last_rect.top() - 3.0 {
                    p3.y -= 1.0;
                }
                p3.y = p3.y.max(0.0);
            }
            let mut p4 = self.pos2xy(pos, self.radius * 1.28);
            p4.y = p3.y;
            painter.line_segment([p1, p2], GRAY_1.to_owned());
            painter.line_segment([p2, p3], GRAY_1.to_owned());

            let align = if pos > self.sequence_length / 2 {
                Align2::RIGHT_CENTER
            } else {
                Align2::LEFT_CENTER
            };
            if let Some(he) = &self.hover_enzyme {
                if he.key == restriction_enzyme_key {
                    painter.rect_filled(he.area, 0.0, Color32::LIGHT_YELLOW);
                }
            }
            last_rect = painter.text(p4, align, label, font_tick.to_owned(), font_color);
            self.restriction_enzyme_sites
                .push(RestrictionEnzymePosition {
                    area: last_rect.to_owned(),
                    key: restriction_enzyme_key.to_owned(),
                });
        }
    }

    /// Displays the length of the DNA sequence in base pairs (bp).
    fn draw_bp(&self, painter: &egui::Painter) {
        let font_label = FontId {
            size: 12.0,
            family: FontFamily::Monospace,
        };
        painter.text(
            self.center.to_owned(),
            Align2::CENTER_TOP,
            format!("{} bp", self.sequence_length),
            font_label,
            Color32::BLACK,
        );
    }

    /// Determines if a feature should be drawn, true for all features that are not of kind
    /// "SOURCE"
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

    /// Determines if a feature should be drawn
    fn normalize_angle(angle: f32) -> f32 {
        angle.rem_euclid(360.0)
    }

    /// Converts a position to an angle
    fn angle(&self, pos: i64) -> f32 {
        let denom = self.sequence_length.max(1) as f32;
        Self::normalize_angle(360.0 * (pos as f32) / denom - 90.0)
    }

    /// Converts a position to Cartesian coordinates.
    fn pos2xy(&self, pos: i64, radius: f32) -> Pos2 {
        let angle = self.angle(pos);
        let t = angle * std::f32::consts::PI / 180.0;
        let x = radius * t.cos() + self.center.x;
        let y = radius * t.sin() + self.center.y;
        Pos2 { x, y }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gb_io::seq::{FeatureKind, Location};

    fn make_test_feature(location: Location) -> Feature {
        Feature {
            kind: FeatureKind::from("mRNA"),
            location,
            qualifiers: vec![("label".into(), Some("tp73 isoform".to_string()))],
        }
    }

    fn test_renderer_with_feature(feature: Feature, sequence_len: usize) -> RenderDnaCircular {
        let mut dna = DNAsequence::from_sequence(&"A".repeat(sequence_len)).expect("valid DNA");
        dna.features_mut().push(feature);
        let dna = Arc::new(RwLock::new(dna));
        let display = Arc::new(RwLock::new(DnaDisplay::default()));
        let mut renderer = RenderDnaCircular::new(dna, display);
        renderer.area = Rect::from_min_max(Pos2::new(0.0, 0.0), Pos2::new(1200.0, 900.0));
        renderer.center = renderer.area.center();
        renderer.radius = 220.0;
        renderer.sequence_length = sequence_len as i64;
        renderer
    }

    #[test]
    fn multipart_join_builds_segments_and_arches() {
        let feature = make_test_feature(Location::Join(vec![
            Location::simple_range(100, 150),
            Location::simple_range(210, 250),
            Location::simple_range(300, 340),
        ]));
        let mut renderer = test_renderer_with_feature(feature, 1000);
        renderer.layout_features();
        assert_eq!(renderer.features.len(), 1);
        let fp = &renderer.features[0];
        assert_eq!(fp.segments.len(), 3);
        assert_eq!(fp.intron_arches.len(), 2);
        assert_eq!(fp.intron_arches[0], (150, 210));
    }

    #[test]
    fn layout_features_handles_complement_join_without_dropping_feature() {
        let feature = make_test_feature(Location::Complement(Box::new(Location::Join(vec![
            Location::simple_range(920, 980),
            Location::simple_range(20, 70),
            Location::simple_range(130, 170),
        ]))));
        let mut renderer = test_renderer_with_feature(feature, 1000);
        renderer.layout_features();
        assert_eq!(renderer.features.len(), 1);
        let fp = &renderer.features[0];
        assert_eq!(fp.segments.len(), 3);
        assert_eq!(fp.intron_arches.len(), 2);
        assert!(fp.outer > fp.inner);
    }
}
