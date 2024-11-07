use crate::{dna_sequence::DNAsequence, main_area_dna::MainAreaDna};
use eframe::egui;

#[derive(Clone, Debug)]
pub struct WindowDna {
    main_area: MainAreaDna,
}

impl WindowDna {
    pub fn new(dna: DNAsequence) -> Self {
        Self {
            main_area: MainAreaDna::new(dna),
        }
    }

    pub fn update(&mut self, ctx: &egui::Context) {
        egui::CentralPanel::default().show(ctx, |ui| self.main_area.render(ctx, ui));
    }

    pub fn name(&self) -> String {
        self.main_area
            .sequence_name()
            .unwrap_or("<Unnamed sequence>".to_string())
    }
}
