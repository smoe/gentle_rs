use crate::{
    dna_sequence::DNAsequence, render_dna_circular::RenderDnaCircular,
    render_dna_linear::RenderDnaLinear,
};
use eframe::egui::{self, Color32, PointerState, Rect, Response, Sense, Ui, Widget};
use gb_io::seq::Feature;
use std::{
    fmt::Debug,
    sync::{Arc, RwLock},
};

#[derive(Debug, Clone)]
pub enum RenderDna {
    Circular(Arc<RwLock<RenderDnaCircular>>),
    Linear(Arc<RwLock<RenderDnaLinear>>),
}

impl RenderDna {
    pub fn new(dna: Arc<RwLock<DNAsequence>>) -> Self {
        let is_circular = dna.read().unwrap().is_circular();
        match is_circular {
            true => RenderDna::Circular(Arc::new(RwLock::new(RenderDnaCircular::new(dna)))),
            false => RenderDna::Linear(Arc::new(RwLock::new(RenderDnaLinear::new(dna)))),
        }
    }

    fn area(&self) -> Rect {
        match self {
            RenderDna::Circular(renderer) => renderer.read().unwrap().area().to_owned(),
            RenderDna::Linear(renderer) => renderer.read().unwrap().area().to_owned(),
        }
    }

    pub fn is_circular(&self) -> bool {
        match self {
            RenderDna::Circular(_) => true,
            RenderDna::Linear(_) => false,
        }
    }

    pub fn on_click(&self, pointer_state: PointerState) {
        match self {
            RenderDna::Circular(renderer) => renderer.write().unwrap().on_click(pointer_state),
            RenderDna::Linear(renderer) => renderer.write().unwrap().on_click(pointer_state),
        }
    }

    pub fn on_double_click(&self, pointer_state: PointerState) {
        match self {
            RenderDna::Circular(renderer) => {
                renderer.write().unwrap().on_double_click(pointer_state)
            }
            RenderDna::Linear(renderer) => renderer.write().unwrap().on_double_click(pointer_state),
        }
    }

    pub fn get_selected_feature_id(&self) -> Option<usize> {
        match self {
            RenderDna::Circular(renderer) => renderer.read().unwrap().selected_feature_number(),
            RenderDna::Linear(renderer) => renderer.read().unwrap().selected_feature_number(),
        }
    }

    pub fn select_feature(&self, feature_number: Option<usize>) {
        match self {
            RenderDna::Circular(renderer) => {
                renderer.write().unwrap().select_feature(feature_number)
            }
            RenderDna::Linear(renderer) => renderer.write().unwrap().select_feature(feature_number),
        }
    }

    fn render(&self, ui: &mut egui::Ui) {
        match self {
            RenderDna::Circular(renderer) => renderer.write().unwrap().render(ui),
            RenderDna::Linear(renderer) => renderer.write().unwrap().render(ui),
        }
    }

    pub fn is_feature_pointy(feature: &Feature) -> bool {
        matches!(
            feature.kind.to_string().to_ascii_uppercase().as_str(),
            "CDS" | "GENE"
        )
    }

    pub fn feature_color(feature: &Feature) -> Color32 {
        match feature.kind.to_string().to_ascii_uppercase().as_str() {
            "CDS" => Color32::RED,
            "GENE" => Color32::BLUE,
            _ => Color32::GRAY,
        }
    }

    pub fn feature_name(feature: &Feature) -> String {
        let mut label_text = match feature.location.find_bounds() {
            Ok((from, to)) => format!("{from}..{to}"),
            Err(_) => String::new(),
        };
        for k in [
            "name",
            "standard_name",
            "gene",
            "protein_id",
            "product",
            "region_name",
            "bound_moiety",
        ] {
            label_text = match feature.qualifier_values(k.into()).next() {
                Some(s) => s.to_owned(),
                None => continue,
            };
            break;
        }
        label_text
    }
}

impl Widget for RenderDna {
    fn ui(self, ui: &mut Ui) -> Response {
        self.render(ui);
        ui.allocate_response(self.area().size(), Sense::click())
    }
}
