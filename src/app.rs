use crate::{
    dna_sequence::{self, DNAsequence},
    left_panel::LeftPanel,
    main_area::MainArea,
    ENZYMES, TRANSLATIONS,
};
use anyhow::{anyhow, Result};
use eframe::egui::{self, menu, Ui};
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
            dna_sequence: vec![],
        };

        // Load two demo sequences
        ret.load_from_file(&cc.egui_ctx, "test_files/pGEX-3X.gb");
        ret.load_from_file(&cc.egui_ctx, "test_files/pGEX_3X.fa");

        // Select first DNA sequence, if any
        if let Some(dna) = ret.dna_sequence.first() {
            ret.main_area = Some(MainArea::new_dna(dna.clone()));
        }
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

    fn add_dna(&mut self, _ctx: &egui::Context, dna_lock: Arc<RwLock<DNAsequence>>) {
        let mut dna = dna_lock.write().unwrap();
        ENZYMES.restriction_enzymes().clone_into(dna.re_mut());
        dna.set_max_re_sites(Some(2));
        dna.update_re_sites();
        drop(dna);

        // Trying to open a window for each DNA sequence
        // let name = match dna.read().unwrap().name() {
        //     Some(name) => name.to_string(),
        //     None => "Unnamed".to_string(),
        // };
        // egui::Window::new(name).show(ctx, |ui| {
        // MainArea::new_dna(dna.clone())
        // ui.label("Hello World!");
        // });

        self.dna_sequence.push(dna_lock.clone());
        self.main_area = Some(MainArea::new_dna(dna_lock));
    }

    fn load_from_file(&mut self, ctx: &egui::Context, path: &str) {
        if let Ok(dna) = Self::load_dna_from_genbank_file(path) {
            self.add_dna(ctx, dna);
        } else if let Ok(dna) = Self::load_dna_from_fasta_file(path) {
            self.add_dna(ctx, dna);
        }
    }

    pub fn render_menu_bar(&mut self, ctx: &egui::Context, ui: &mut Ui) {
        menu::bar(ui, |ui| {
            ui.menu_button(TRANSLATIONS.get("m_file"), |ui| {
                if ui.button(TRANSLATIONS.get("m_open")).clicked() {
                    if let Some(path) = rfd::FileDialog::new().pick_file() {
                        if let Some(path) = Some(path.display().to_string()) {
                            self.load_from_file(ctx, &path);
                        }
                    }
                }
            });
        });
    }
}

impl eframe::App for GENtleApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::TopBottomPanel::top("top").show(ctx, |ui| {
            self.render_menu_bar(ctx, ui);
        });
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
