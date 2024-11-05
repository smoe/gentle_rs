use crate::{
    dna_sequence::{self, DNAsequence},
    left_panel::LeftPanel,
    main_area::MainArea,
    ENZYMES,
};
use anyhow::{anyhow, Result};
use eframe::egui;
use std::sync::{Arc, RwLock};

#[derive(Default)]
pub struct GENtleApp {
    main_area: Option<MainArea>,
    left_panel: LeftPanel,
    dna_sequence: Vec<Arc<RwLock<DNAsequence>>>,
}

impl GENtleApp {
    pub fn new(cc: &eframe::CreationContext<'_>) -> Self {
        // Customize egui here with cc.egui_ctx.set_fonts and cc.egui_ctx.set_visuals.
        // Restore app state using cc.storage (requires the "persistence" feature).
        // Use the cc.gl (a glow::Context) to create graphics shaders and buffers that you can use
        // for e.g. egui::PaintCallback.

        egui_extras::install_image_loaders(&cc.egui_ctx);
        let mut ret = Self {
            main_area: None,
            left_panel: Default::default(),
            dna_sequence: Self::load_demo_data(), // For testing
        };

        // Select first DNA sequence, if any
        if let Some(dna) = ret.dna_sequence.first() {
            ret.main_area = Some(MainArea::new_dna(dna.clone()));
        }
        ret
    }

    fn load_demo_data() -> Vec<Arc<RwLock<DNAsequence>>> {
        // Load two demo sequences
        let ret = vec![
            Self::load_dna_from_genbank_file("test_files/pGEX-3X.gb").unwrap(),
            Self::load_dna_from_fasta_file("test_files/pGEX_3X.fa").unwrap(),
        ];

        // Set up restriction enzyme sites for first one
        ENZYMES
            .restriction_enzymes()
            .clone_into(ret[0].write().unwrap().re_mut());
        ret[0].write().unwrap().update_re_sites();

        ret
    }

    pub fn dna_sequence(&self) -> &Vec<Arc<RwLock<DNAsequence>>> {
        &self.dna_sequence
    }

    pub fn main_area(&self) -> &Option<MainArea> {
        &self.main_area
    }

    pub fn set_main_area(&mut self, main_area: Option<MainArea>) {
        self.main_area = main_area;
    }

    fn load_dna_from_genbank_file(filename: &str) -> Result<Arc<RwLock<DNAsequence>>> {
        let dna = dna_sequence::DNAsequence::from_genbank_file(filename)?
            .pop()
            .ok_or_else(|| anyhow!("Could not read genbank file {filename}"))?;
        Ok(Arc::new(RwLock::new(dna)))
    }

    fn load_dna_from_fasta_file(filename: &str) -> Result<Arc<RwLock<DNAsequence>>> {
        let dna = dna_sequence::DNAsequence::from_fasta_file(filename)?
            .pop()
            .ok_or_else(|| anyhow!("Could not read fasta file {filename}"))?;
        Ok(Arc::new(RwLock::new(dna)))
    }
}

impl eframe::App for GENtleApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::SidePanel::left("left_panel")
            .resizable(true)
            .default_width(200.0)
            .min_width(200.0)
            .show(ctx, |ui| {
                let mut left_panel = self.left_panel.clone();
                left_panel.render(self, ui);
                self.left_panel = left_panel;
            });
        egui::CentralPanel::default().show(ctx, |ui| match self.main_area {
            Some(ref mut main_area) => {
                main_area.render(ctx, ui);
            }
            None => {
                ui.label("No main area");
            }
        });
    }
}
