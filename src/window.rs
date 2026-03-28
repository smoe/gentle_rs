//! Generic GUI window abstraction.

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

    pub fn new_dna_lazy(seq_id: String, engine: Arc<RwLock<GentleEngine>>) -> Self {
        Self::Dna(WindowDna::new_lazy(seq_id, engine))
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

    pub fn set_pool_context(&mut self, pool_seq_ids: Vec<String>) {
        match self {
            Self::Dna(window) => window.set_pool_context(pool_seq_ids),
        }
    }

    pub fn focus_dotplot_analysis(&mut self, dotplot_id: &str) {
        match self {
            Self::Dna(window) => window.focus_dotplot_analysis(dotplot_id),
        }
    }

    pub fn focus_flexibility_track_analysis(&mut self, track_id: &str) {
        match self {
            Self::Dna(window) => window.focus_flexibility_track_analysis(track_id),
        }
    }

    pub fn sequence_id(&self) -> Option<String> {
        match self {
            Self::Dna(window) => window.sequence_id(),
        }
    }

    pub fn selection_range_0based(&self) -> Option<(usize, usize)> {
        match self {
            Self::Dna(window) => window.selection_range_0based(),
        }
    }

    pub fn collect_open_auxiliary_window_entries(&self) -> Vec<(egui::ViewportId, String, String)> {
        match self {
            Self::Dna(window) => window.collect_open_auxiliary_window_entries(),
        }
    }

    pub fn render_pcr_designer_specialist(&mut self, ui: &mut egui::Ui, ctx: &egui::Context) {
        match self {
            Self::Dna(window) => window.render_pcr_designer_specialist(ui, ctx),
        }
    }

    pub fn refresh_from_engine_settings(&mut self) {
        match self {
            Self::Dna(window) => window.refresh_from_engine_settings(),
        }
    }

    pub fn refresh_from_engine_state(&mut self) {
        match self {
            Self::Dna(window) => window.refresh_from_engine_state(),
        }
    }

    pub fn take_close_requested(&mut self) -> bool {
        match self {
            Self::Dna(window) => window.take_close_requested(),
        }
    }
}
