use crate::{dna_sequence::DNAsequence, window_dna::WindowDna};
use eframe::egui;

#[derive(Clone, Debug)]
pub enum Window {
    Dna(WindowDna),
}

impl Window {
    pub fn new_dna(dna: DNAsequence) -> Self {
        Self::Dna(WindowDna::new(dna))
    }

    pub fn update(&mut self, ctx: &egui::Context) {
        match self {
            Self::Dna(window) => window.update(ctx),
        }
    }

    pub fn name(&self) -> String {
        match self {
            Self::Dna(window) => window.name(),
        }
    }
}
