use eframe::egui::{self, PointerState, Rect};

use crate::{dna_display::DnaDisplay, dna_sequence::DNAsequence};
use std::sync::{Arc, RwLock};

#[derive(Debug, Clone)]
struct FeaturePosition {
    _feature_number: usize,
    _from: i64,
    _to: i64,
}

#[derive(Debug, Clone)]
pub struct RenderDnaLinear {
    area: Rect,
    display: Arc<RwLock<DnaDisplay>>,
    _dna: Arc<RwLock<DNAsequence>>,
    _features: Vec<FeaturePosition>,
    selected_feature_number: Option<usize>,
}

impl RenderDnaLinear {
    pub fn new(dna: Arc<RwLock<DNAsequence>>, display: Arc<RwLock<DnaDisplay>>) -> Self {
        Self {
            area: Rect::NOTHING,
            display,
            _dna: dna,
            _features: vec![],
            selected_feature_number: None,
        }
    }

    pub fn area(&self) -> &Rect {
        &self.area
    }

    fn layout_needs_recomputing(&mut self, ui: &mut egui::Ui) -> bool {
        let mut ret = false;

        // Recompute layout if area has changed
        let new_area = ui.available_rect_before_wrap();
        if self.area != new_area {
            ret = true;
            self.area = new_area;
        }

        // Recompute layout if update flag is set
        ret = ret
            || self
                .display
                .read()
                .unwrap()
                .update_layout()
                .update_map_dna();

        ret
    }

    fn layout_was_updated(&self) {
        self.display
            .write()
            .unwrap()
            .update_layout_mut()
            .map_dna_updated();
    }

    pub fn render(&mut self, ui: &mut egui::Ui) {
        if self.layout_needs_recomputing(ui) {
            // TODO layout recompute
            self.layout_was_updated();
        }

        if self.display.read().unwrap().show_restriction_enzyme_sites() {
            ui.heading("Linear DNA (RE)");
        } else {
            ui.heading("Linear DNA");
        }
    }

    pub fn on_click(&mut self, pointer_state: PointerState) {
        if let Some(_pos) = pointer_state.latest_pos() {
            println!("Clicking on linear DNA sequence not implemented yet");
        }
    }

    pub fn on_double_click(&mut self, pointer_state: PointerState) {
        if let Some(_pos) = pointer_state.latest_pos() {
            println!("Double-clicking on linear DNA sequence not implemented yet");
        }
    }

    pub fn selected_feature_number(&self) -> Option<usize> {
        self.selected_feature_number.to_owned()
    }

    pub fn select_feature(&mut self, feature_number: Option<usize>) {
        self.selected_feature_number = feature_number;
    }
}
