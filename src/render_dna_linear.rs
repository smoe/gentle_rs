use eframe::egui::{self, PointerState};

use crate::dna_sequence::DNAsequence;
use std::sync::{Arc, Mutex};

#[derive(Debug, Clone)]
struct FeaturePosition {
    _feature_number: usize,
    _from: i64,
    _to: i64,
}

#[derive(Debug)]
pub struct RenderDnaLinear {
    _dna: Arc<Mutex<DNAsequence>>,
    _features: Vec<FeaturePosition>,
    selected_feature_number: Option<usize>,
}

impl RenderDnaLinear {
    pub fn new(dna: Arc<Mutex<DNAsequence>>) -> Self {
        Self {
            _dna: dna,
            _features: vec![],
            selected_feature_number: None,
        }
    }

    pub fn render(&mut self, ui: &mut egui::Ui) {
        ui.heading("Linear DNA");
    }

    pub fn on_click(&mut self, pointer_state: PointerState) {
        if let Some(_pos) = pointer_state.latest_pos() {
            todo!()
        }
    }

    pub fn selected_feature_number(&self) -> Option<usize> {
        self.selected_feature_number.to_owned()
    }

    pub fn select_feature(&mut self, feature_number: Option<usize>) {
        self.selected_feature_number = feature_number;
    }
}
