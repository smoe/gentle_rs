use crate::{
    dna_sequence::DNAsequence, render_dna_circular::RenderDnaCircular,
    render_dna_linear::RenderDnaLinear,
};
use eframe::egui::{self, PointerState};
use gb_io::seq::Feature;
use std::{
    fmt::Debug,
    sync::{Arc, RwLock},
};

#[derive(Debug)]
pub enum RenderDnaEnum {
    Circular(RenderDnaCircular),
    Linear(RenderDnaLinear),
}

impl RenderDnaEnum {
    pub fn new(dna: Arc<RwLock<DNAsequence>>) -> Self {
        let is_circular = dna.read().unwrap().is_circular();
        match is_circular {
            true => RenderDnaEnum::Circular(RenderDnaCircular::new(dna)),
            false => RenderDnaEnum::Linear(RenderDnaLinear::new(dna)),
        }
    }

    pub fn is_circular(&self) -> bool {
        match self {
            RenderDnaEnum::Circular(_) => true,
            RenderDnaEnum::Linear(_) => false,
        }
    }

    pub fn on_click(&mut self, pointer_state: PointerState) {
        match self {
            RenderDnaEnum::Circular(renderer) => renderer.on_click(pointer_state),
            RenderDnaEnum::Linear(renderer) => renderer.on_click(pointer_state),
        }
    }

    pub fn get_selected_feature_id(&self) -> Option<usize> {
        match self {
            RenderDnaEnum::Circular(renderer) => renderer.selected_feature_number(),
            RenderDnaEnum::Linear(renderer) => renderer.selected_feature_number(),
        }
    }

    pub fn select_feature(&mut self, feature_number: Option<usize>) {
        match self {
            RenderDnaEnum::Circular(renderer) => renderer.select_feature(feature_number),
            RenderDnaEnum::Linear(renderer) => renderer.select_feature(feature_number),
        }
    }

    pub fn render(&mut self, ui: &mut egui::Ui) {
        match self {
            RenderDnaEnum::Circular(renderer) => renderer.render(ui),
            RenderDnaEnum::Linear(renderer) => renderer.render(ui),
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
