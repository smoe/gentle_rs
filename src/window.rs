use crate::{dna_sequence::DNAsequence, engine::GentleEngine, window_dna::WindowDna};
use eframe::egui;
use std::sync::{Arc, RwLock};

#[derive(Clone, Debug)]
pub enum Window {
    Dna(WindowDna),
}

impl Window {
    pub fn new_dna(dna: DNAsequence, seq_id: String, engine: Arc<RwLock<GentleEngine>>) -> Self {
        Self::Dna(WindowDna::new(dna, seq_id, engine))
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
