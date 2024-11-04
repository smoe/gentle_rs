use crate::{dna_sequence::DNAsequence, render_dna::RenderDnaEnum};
use eframe::egui::{
    self, Align2, Color32, FontFamily, FontId, PointerState, Pos2, Rect, Shape, Stroke,
};
use gb_io::seq::Feature;
use lazy_static::lazy_static;
use std::sync::{Arc, RwLock};

lazy_static! {
    pub static ref BLACK_1: Stroke = Stroke {
        width: 1.0,
        color: Color32::BLACK,
    };
}

#[derive(Debug, Clone)]
struct FeaturePosition {
    feature_number: usize,
    from: i64,
    to: i64,
    angle_start: f32,
    angle_stop: f32,
    inner: f32,
    outer: f32,
    to_90: i64,
    is_pointy: bool,
    color: Color32,
    band: f32,
    label: String,
}

impl FeaturePosition {
    fn contains_point(&self, angle: f32, distance: f32) -> bool {
        self.inner <= distance
            && distance <= self.outer
            && self.angle_start <= angle
            && angle <= self.angle_stop
    }
}

#[derive(Debug)]
pub struct RenderDnaCircular {
    dna: Arc<RwLock<DNAsequence>>,
    sequence_length: i64,
    area: Rect,
    center: Pos2,
    radius: f32,
    features: Vec<FeaturePosition>,
    selected_feature_number: Option<usize>,
}

impl RenderDnaCircular {
    pub fn new(dna: Arc<RwLock<DNAsequence>>) -> Self {
        Self {
            dna,
            sequence_length: 0,
            area: Rect::NOTHING,
            center: Pos2::ZERO,
            radius: 0.0,
            features: vec![],
            selected_feature_number: None,
        }
    }

    pub fn on_click(&mut self, pointer_state: PointerState) {
        if let Some(pos) = pointer_state.latest_pos() {
            let (angle, distance) = self.get_angle_distance(pos);
            let clicked_features = self
                .features
                .iter()
                .filter(|feature| feature.contains_point(angle, distance))
                .collect::<Vec<_>>();
            self.selected_feature_number = clicked_features.first().map(|f| f.feature_number);
        }
    }

    pub fn selected_feature_number(&self) -> Option<usize> {
        self.selected_feature_number.to_owned()
    }

    pub fn select_feature(&mut self, feature_number: Option<usize>) {
        self.selected_feature_number = feature_number;
    }

    pub fn render(&mut self, ui: &mut egui::Ui) {
        self.features.clear();
        self.area = ui.available_rect_before_wrap();
        self.radius = self.area.width().min(self.area.height()) * 0.4;
        self.center = self.area.center();
        self.sequence_length = self.dna.read().expect("DNA lock poisoned").len() as i64;

        let painter = ui.painter();
        self.draw_backbone(painter);
        self.draw_main_label(painter);
        self.draw_bp(painter);
        self.layout_features(); // TODO cache and invalidate on update
        self.draw_features(painter);
    }

    fn get_angle_distance(&self, pos: Pos2) -> (f32, f32) {
        let diff_x = pos.x - self.center.x;
        let diff_y = pos.y - self.center.y;
        let angle = diff_y.atan2(diff_x) * 180.0 / std::f32::consts::PI + 90.0;
        let angle = if angle < 0.0 { angle + 360.0 } else { angle };
        let distance = (diff_x.powi(2) + diff_y.powi(2)).sqrt();
        (angle, distance)
    }

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

    fn layout_features(&mut self) {
        self.features.clear();
        let features = self
            .dna
            .read()
            .expect("DNA lock poisoned")
            .features()
            .to_owned();
        for (feature_number, feature) in features.iter().enumerate() {
            let fp_opt = match feature.location {
                gb_io::seq::Location::Range(from, to) => {
                    self.layout_feature_from_range(feature, from, to)
                }
                gb_io::seq::Location::Between(_, _) => todo!(),
                gb_io::seq::Location::Complement(_) => todo!(),
                gb_io::seq::Location::Join(_) => todo!(),
                gb_io::seq::Location::Order(_) => todo!(),
                gb_io::seq::Location::Bond(_) => todo!(),
                gb_io::seq::Location::OneOf(_) => todo!(),
                gb_io::seq::Location::External(_, _) => todo!(),
                gb_io::seq::Location::Gap(_) => todo!(),
            };
            if let Some(mut fp) = fp_opt {
                fp.feature_number = feature_number;
                self.features.push(fp);
            }
        }
    }

    fn feature_thickness(&self) -> f32 {
        self.radius / 20.0
    }

    fn layout_feature_from_range(
        &self,
        feature: &Feature,
        start: (i64, gb_io::seq::Before),
        end: (i64, gb_io::seq::After),
    ) -> Option<FeaturePosition> {
        if !Self::draw_feature(feature) {
            return None;
        }
        let mut ret: FeaturePosition = FeaturePosition {
            feature_number: 0,
            from: start.0,
            to: end.0,
            angle_start: 0.0,
            angle_stop: 0.0,
            inner: 0.0,
            outer: 0.0,
            to_90: 0,
            is_pointy: Self::is_feature_pointy(feature),
            color: Self::feature_color(feature),
            band: Self::feature_band(feature),
            label: RenderDnaEnum::feature_name(feature),
        };
        if Self::feature_band(feature) == 0.0 {
            ret.inner = self.radius - self.feature_thickness() / 2.0;
            ret.outer = self.radius + self.feature_thickness() / 2.0;
        } else {
            ret.inner = self.radius + Self::feature_band(feature) * self.feature_thickness();
            ret.outer = self.radius + 2.0 * Self::feature_band(feature) * self.feature_thickness();
        }
        ret.to_90 = if ret.is_pointy {
            ret.to - (ret.to - ret.from) / 20
        } else {
            ret.to
        };
        ret.angle_start = self.angle(ret.from);
        ret.angle_stop = self.angle(ret.to_90);
        Some(ret)
    }

    fn draw_features(&mut self, painter: &egui::Painter) {
        for feature in &self.features {
            self.draw_feature_from_range(painter, feature);
        }
    }

    fn draw_feature_from_range(&self, painter: &egui::Painter, ret: &FeaturePosition) {
        let mut feature_points: Vec<Pos2> = vec![];
        feature_points.push(self.pos2xy(ret.from, ret.outer));
        feature_points.push(self.pos2xy(ret.from, ret.inner));

        let points = self.generate_arc(ret.inner, ret.angle_start, ret.angle_stop);
        feature_points.extend(points);

        if ret.is_pointy {
            feature_points.push(self.pos2xy(ret.to, (ret.outer + ret.inner) / 2.0));
        }

        feature_points.push(self.pos2xy(ret.to_90, ret.outer));

        let points = self.generate_arc(ret.outer, ret.angle_stop, ret.angle_start);
        feature_points.extend(points);

        let stroke = Stroke {
            width: 1.0,
            color: ret.color,
        };
        let line = Shape::closed_line(feature_points, stroke);
        painter.add(line);

        let font_feature = FontId {
            size: 10.0,
            family: FontFamily::Monospace,
        };

        // Draw feature label
        let middle = (ret.to + ret.from) / 2;
        let point = self.pos2xy(
            middle,
            ret.outer + ret.band * self.feature_thickness() / 2.0,
        );
        let align = if ret.band < 0.0 {
            // Inside
            if middle < self.sequence_length / 2 {
                Align2::RIGHT_CENTER
            } else {
                Align2::LEFT_CENTER
            }
        } else {
            // Outside
            if middle > self.sequence_length / 2 {
                Align2::RIGHT_CENTER
            } else {
                Align2::LEFT_CENTER
            }
        };
        painter.text(point, align, ret.label.to_owned(), font_feature, ret.color);
    }

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

    fn feature_band(feature: &Feature) -> f32 {
        match feature.kind.to_string().to_ascii_uppercase().as_str() {
            "CDS" => 1.0,
            "GENE" => -1.0,
            _ => 0.0,
        }
    }

    fn draw_feature(feature: &Feature) -> bool {
        let feature_kind = feature.kind.to_string().to_ascii_uppercase();
        feature_kind != "SOURCE"
    }

    fn is_feature_pointy(feature: &Feature) -> bool {
        matches!(
            feature.kind.to_string().to_ascii_uppercase().as_str(),
            "CDS" | "GENE"
        )
    }

    fn feature_color(feature: &Feature) -> Color32 {
        match feature.kind.to_string().to_ascii_uppercase().as_str() {
            "CDS" => Color32::RED,
            "GENE" => Color32::BLUE,
            _ => Color32::GRAY,
        }
    }

    fn angle(&self, pos: i64) -> f32 {
        360.0 * (pos as f32) / (self.sequence_length as f32) - 90.0
    }

    fn pos2xy(&self, pos: i64, radius: f32) -> Pos2 {
        let angle = self.angle(pos);
        let t = angle * std::f32::consts::PI / 180.0;
        let x = radius * t.cos() + self.center.x;
        let y = radius * t.sin() + self.center.y;
        Pos2 { x, y }
    }
}
