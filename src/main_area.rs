use std::sync::{Arc, Mutex};

use crate::{dna_sequence::DNAsequence, main_area_dna::MainAreaDna};
use eframe::egui;

#[derive(Debug)]
pub enum MainArea {
    Dna(MainAreaDna),
}

impl MainArea {
    pub fn render(&mut self, ctx: &egui::Context, ui: &mut egui::Ui) {
        match self {
            MainArea::Dna(main_area) => {
                main_area.render(ctx, ui);
            }
        }
    }

    pub fn new_dna(dna: Arc<Mutex<DNAsequence>>) -> Self {
        MainArea::Dna(MainAreaDna::new(dna))
    }
}
