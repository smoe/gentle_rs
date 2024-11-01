use crate::dna_sequence::DNAsequence;
use eframe::egui::{
    self, Align2, Color32, FontFamily, FontId, PointerState, Pos2, Rect, Shape, Stroke,
};
use gb_io::seq::Feature;
use lazy_static::lazy_static;
use std::sync::{Arc, Mutex};

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
    dna: Arc<Mutex<DNAsequence>>,
    area: Rect,
    center: Pos2,
    radius: f32,
    features: Vec<FeaturePosition>,
    selected_feature_number: Option<usize>,
}

impl RenderDnaCircular {
    pub fn new(dna: Arc<Mutex<DNAsequence>>) -> Self {
        Self {
            dna,
            area: Rect::NOTHING,
            center: Pos2::ZERO,
            radius: 0.0,
            features: vec![],
            selected_feature_number: None,
        }
    }

    fn get_angle_distance(&self, pos: Pos2) -> (f32, f32) {
        let diff_x = pos.x - self.center.x;
        let diff_y = pos.y - self.center.y;
        let angle = diff_y.atan2(diff_x) * 180.0 / std::f32::consts::PI + 90.0;
        let angle = if angle < 0.0 { angle + 360.0 } else { angle };
        let distance = (diff_x.powi(2) + diff_y.powi(2)).sqrt();
        (angle, distance)
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

    pub fn render(&mut self, ui: &mut egui::Ui) {
        self.features.clear();
        self.area = ui.available_rect_before_wrap();
        self.radius = self.area.width().min(self.area.height()) * 0.4;
        self.center = self.area.center();
        let bp = self.dna.lock().expect("DNA lock poisoned").len() as i64;

        let painter = ui.painter();
        self.draw_backbone(painter, bp);
        self.draw_main_label(painter);
        self.draw_bp(painter, bp);
        self.draw_features(painter, bp);
    }

    fn draw_backbone(&mut self, painter: &egui::Painter, bp: i64) {
        painter.circle_stroke(self.center.to_owned(), self.radius, BLACK_1.to_owned());
        let mut tick: i64 = 1;
        while tick * 10 < bp {
            tick *= 10;
        }

        let font_tick = FontId {
            size: 9.0,
            family: FontFamily::Monospace,
        };

        let mut pos = tick; // Skip 0 point
        while pos < bp {
            let p1 = self.pos2xy(pos, bp, self.radius);
            let p2 = self.pos2xy(pos, bp, self.radius * 0.87);
            let p3 = self.pos2xy(pos, bp, self.radius * 0.85);
            painter.line_segment([p1, p2], BLACK_1.to_owned());

            let align = if pos > bp / 2 {
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

    fn draw_features(&mut self, painter: &egui::Painter, bp: i64) {
        let features = self
            .dna
            .lock()
            .expect("DNA lock poisoned")
            .features()
            .to_owned();
        for (feature_number, feature) in features.iter().enumerate() {
            let fp_opt = match feature.location {
                gb_io::seq::Location::Range(from, to) => {
                    self.draw_feature_from_range(painter, feature, bp, from, to)
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

    fn draw_feature_from_range(
        &self,
        painter: &egui::Painter,
        feature: &Feature,
        bp: i64,
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
        };
        let feature_thickness = self.radius / 20.0;
        if Self::feature_band(feature) == 0.0 {
            ret.inner = self.radius - feature_thickness / 2.0;
            ret.outer = self.radius + feature_thickness / 2.0;
        } else {
            ret.inner = self.radius + Self::feature_band(feature) * feature_thickness;
            ret.outer = self.radius + 2.0 * Self::feature_band(feature) * feature_thickness;
        }

        // let feature_kind = feature.kind.to_string().to_ascii_uppercase();
        // if feature_label != "bla" || feature_kind != "GENE" {
        //     return; // TESTING FIXME
        // }
        // println!("{feature_label}: {}", feature.kind);

        let mut feature_points: Vec<Pos2> = vec![];
        feature_points.push(self.pos2xy(ret.from, bp, ret.outer));
        feature_points.push(self.pos2xy(ret.from, bp, ret.inner));

        let is_pointy = Self::is_feature_pointy(feature);

        let to_90 = if is_pointy {
            ret.to - (ret.to - ret.from) / 20
        } else {
            ret.to
        };
        ret.angle_start = self.angle(ret.from, bp);
        ret.angle_stop = self.angle(to_90, bp);
        let points = self.generate_arc(ret.inner, ret.angle_start, ret.angle_stop);
        feature_points.extend(points);

        if is_pointy {
            feature_points.push(self.pos2xy(ret.to, bp, (ret.outer + ret.inner) / 2.0));
        }

        feature_points.push(self.pos2xy(to_90, bp, ret.outer));

        let points = self.generate_arc(ret.outer, ret.angle_stop, ret.angle_start);
        feature_points.extend(points);

        // let feature_polygon =
        //     Shape::convex_polygon(feature_points, feature_color.to_owned(), BLACK_1.to_owned());
        // painter.add(feature_polygon);

        let stroke = Stroke {
            width: 1.0,
            color: Self::feature_color(feature),
        };
        let line = Shape::closed_line(feature_points, stroke);
        painter.add(line);

        let font_feature = FontId {
            size: 10.0,
            family: FontFamily::Monospace,
        };

        // Draw feature label
        let band = Self::feature_band(feature);
        let feature_label = Self::feature_name(feature);
        let middle = (ret.to + ret.from) / 2;
        let point = self.pos2xy(middle, bp, ret.outer + band * feature_thickness / 2.0);
        let align = if band < 0.0 {
            // Inside
            if middle < bp / 2 {
                Align2::RIGHT_CENTER
            } else {
                Align2::LEFT_CENTER
            }
        } else {
            // Outside
            if middle > bp / 2 {
                Align2::RIGHT_CENTER
            } else {
                Align2::LEFT_CENTER
            }
        };
        painter.text(
            point,
            align,
            feature_label,
            font_feature,
            Self::feature_color(feature),
        );

        Some(ret)
    }

    fn feature_name(feature: &Feature) -> String {
        let mut label_text = String::new();
        for k in ["gene", "product", "standard_name", "protein_id"] {
            label_text = match feature.qualifier_values(k.into()).next() {
                Some(s) => s.to_owned(),
                None => continue,
            };
            break;
        }
        label_text
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
            .lock()
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

    fn draw_bp(&self, painter: &egui::Painter, bp: i64) {
        let font_label = FontId {
            size: 12.0,
            family: FontFamily::Monospace,
        };
        painter.text(
            self.center.to_owned(),
            Align2::CENTER_TOP,
            format!("{bp} bp"),
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

    fn angle(&self, pos: i64, max_pos: i64) -> f32 {
        360.0 * (pos as f32) / (max_pos as f32) - 90.0
    }

    fn pos2xy(&self, pos: i64, max_pos: i64, radius: f32) -> Pos2 {
        let angle = self.angle(pos, max_pos);
        let t = angle * std::f32::consts::PI / 180.0;
        let x = radius * t.cos() + self.center.x;
        let y = radius * t.sin() + self.center.y;
        Pos2 { x, y }
    }
}
