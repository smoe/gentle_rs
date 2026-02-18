use crate::{
    dna_display::{DnaDisplay, TfbsDisplayCriteria},
    dna_sequence::DNAsequence,
    render_dna_circular::RenderDnaCircular,
    render_dna_linear::RenderDnaLinear,
    restriction_enzyme::RestrictionEnzymeKey,
};
use eframe::egui::{self, Color32, PointerState, Rect, Response, Sense, Ui, Widget};
use gb_io::seq::Feature;
use std::{
    fmt::Debug,
    sync::{Arc, RwLock},
};

#[derive(Debug, Clone, PartialEq)]
pub struct RestrictionEnzymePosition {
    pub area: Rect,
    pub key: RestrictionEnzymeKey,
}

impl RestrictionEnzymePosition {
    pub fn key(&self) -> &RestrictionEnzymeKey {
        &self.key
    }
    pub fn area(&self) -> Rect {
        self.area
    }
}

#[derive(Debug, Clone)]
pub enum RenderDna {
    Circular(Arc<RwLock<RenderDnaCircular>>),
    Linear(Arc<RwLock<RenderDnaLinear>>),
}

impl RenderDna {
    pub fn new(dna: Arc<RwLock<DNAsequence>>, display: Arc<RwLock<DnaDisplay>>) -> Self {
        let is_circular = dna.read().unwrap().is_circular();
        match is_circular {
            true => {
                RenderDna::Circular(Arc::new(RwLock::new(RenderDnaCircular::new(dna, display))))
            }
            false => RenderDna::Linear(Arc::new(RwLock::new(RenderDnaLinear::new(dna, display)))),
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

    pub fn on_hover(&self, pointer_state: PointerState) {
        match self {
            RenderDna::Circular(renderer) => renderer.write().unwrap().on_hover(pointer_state),
            RenderDna::Linear(renderer) => renderer.write().unwrap().on_hover(pointer_state),
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
            "TFBS" | "TF_BINDING_SITE" | "PROTEIN_BIND" => Color32::from_rgb(35, 120, 35),
            _ => Color32::GRAY,
        }
    }

    pub fn is_tfbs_feature(feature: &Feature) -> bool {
        matches!(
            feature.kind.to_string().to_ascii_uppercase().as_str(),
            "TFBS" | "TF_BINDING_SITE" | "PROTEIN_BIND"
        )
    }

    fn feature_qualifier_f64(feature: &Feature, key: &str) -> Option<f64> {
        feature
            .qualifier_values(key.into())
            .next()
            .and_then(|v| v.trim().parse::<f64>().ok())
    }

    pub fn tfbs_feature_passes_display_filter(
        feature: &Feature,
        criteria: TfbsDisplayCriteria,
    ) -> bool {
        if !Self::is_tfbs_feature(feature) {
            return true;
        }
        if criteria.use_llr_bits {
            let Some(value) = Self::feature_qualifier_f64(feature, "llr_bits") else {
                return false;
            };
            if value < criteria.min_llr_bits {
                return false;
            }
        }
        if criteria.use_llr_quantile {
            let Some(value) = Self::feature_qualifier_f64(feature, "llr_quantile") else {
                return false;
            };
            if value < criteria.min_llr_quantile {
                return false;
            }
        }
        if criteria.use_true_log_odds_bits {
            let value = Self::feature_qualifier_f64(feature, "true_log_odds_bits")
                .or_else(|| Self::feature_qualifier_f64(feature, "log_odds_ratio_bits"));
            let Some(value) = value else {
                return false;
            };
            if value < criteria.min_true_log_odds_bits {
                return false;
            }
        }
        if criteria.use_true_log_odds_quantile {
            let value = Self::feature_qualifier_f64(feature, "true_log_odds_quantile")
                .or_else(|| Self::feature_qualifier_f64(feature, "log_odds_ratio_quantile"));
            let Some(value) = value else {
                return false;
            };
            if value < criteria.min_true_log_odds_quantile {
                return false;
            }
        }
        true
    }

    pub fn is_source_feature(feature: &Feature) -> bool {
        feature.kind.to_string().to_ascii_uppercase() == "SOURCE"
    }

    pub fn feature_name(feature: &Feature) -> String {
        let mut label_text = match feature.location.find_bounds() {
            Ok((from, to)) => format!("{from}..{to}"),
            Err(_) => String::new(),
        };
        for k in [
            "label",
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
