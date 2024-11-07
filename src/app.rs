use std::{
    collections::HashMap,
    sync::{Arc, RwLock},
};

use crate::{
    dna_sequence::{self, DNAsequence},
    window::Window,
    ENZYMES, TRANSLATIONS,
};
use anyhow::{anyhow, Result};
use eframe::egui::{self, menu, Pos2, Ui, ViewportId};

#[derive(Default)]
pub struct GENtleApp {
    new_windows: Vec<Window>,
    windows: HashMap<ViewportId, Arc<RwLock<Window>>>,
    windows_to_close: Arc<RwLock<Vec<ViewportId>>>,
    viewport_id_counter: usize,
    update_has_run_before: bool,
}

impl GENtleApp {
    pub fn new() -> Self {
        let mut ret = Self::default();

        // Load test sequences
        ret.load_from_file("test_files/pGEX-3X.gb");
        ret.load_from_file("test_files/pGEX_3X.fa");

        ret
    }

    fn load_dna_from_genbank_file(filename: &str) -> Result<DNAsequence> {
        let dna = dna_sequence::DNAsequence::from_genbank_file(filename)?
            .pop()
            .ok_or_else(|| anyhow!("Could not read genbank file {filename}"))?;
        Ok(dna)
    }

    fn load_dna_from_fasta_file(filename: &str) -> Result<DNAsequence> {
        let dna = dna_sequence::DNAsequence::from_fasta_file(filename)?
            .pop()
            .ok_or_else(|| anyhow!("Could not read fasta file {filename}"))?;
        Ok(dna)
    }

    fn new_dna_window(&mut self, mut dna: DNAsequence) {
        ENZYMES.restriction_enzymes().clone_into(dna.re_mut());
        dna.set_max_re_sites(Some(2));
        dna.update_computed_features();

        self.new_windows.push(Window::new_dna(dna));
    }

    fn load_from_file(&mut self, path: &str) {
        if let Ok(dna) = Self::load_dna_from_genbank_file(path) {
            self.new_dna_window(dna);
        } else if let Ok(dna) = Self::load_dna_from_fasta_file(path) {
            self.new_dna_window(dna);
        } else {
            // TODO error, could not load file
        }
    }

    pub fn render_menu_bar(&mut self, ui: &mut Ui) {
        menu::bar(ui, |ui| {
            ui.menu_button(TRANSLATIONS.get("m_file"), |ui| {
                if ui.button(TRANSLATIONS.get("m_open")).clicked() {
                    if let Some(path) = rfd::FileDialog::new().pick_file() {
                        if let Some(path) = Some(path.display().to_string()) {
                            self.load_from_file(&path);
                        }
                    }
                }
            });
        });
    }

    fn show_window(&self, ctx: &egui::Context, id: ViewportId, window: Arc<RwLock<Window>>) {
        let windows_to_close = self.windows_to_close.clone();
        let window_number = self.get_window_number_from_id(id);
        let window_pos = Pos2 {
            x: window_number as f32 * 200.0,
            y: window_number as f32 * 200.0,
        };
        let window_title = window.read().unwrap().name();
        ctx.show_viewport_deferred(
            id,
            egui::ViewportBuilder::default()
                .with_title(window_title)
                // .with_maximized(true),
                .with_position(window_pos),
            move |ctx, class| {
                assert!(
                    class == egui::ViewportClass::Deferred,
                    "This egui backend doesn't support multiple viewports"
                );

                // Draw the window
                window.write().unwrap().update(ctx);

                // "Close window" action
                if ctx.input(|i| i.viewport().close_requested()) {
                    windows_to_close.write().unwrap().push(id);
                }
            },
        );
    }

    fn get_window_number_from_id(&self, id: ViewportId) -> usize {
        let window_number = self
            .windows
            .keys()
            .enumerate()
            .find(|(_num, viewport_id)| **viewport_id == id)
            .unwrap()
            .0;
        window_number
    }
}

impl eframe::App for GENtleApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        if !self.update_has_run_before {
            egui_extras::install_image_loaders(ctx);
            self.update_has_run_before = true;
        }

        // Show menu bar
        egui::TopBottomPanel::top("top").show(ctx, |ui| {
            self.render_menu_bar(ui);
        });

        // Show main window
        egui::CentralPanel::default().show(ctx, |ui| {
            ui.label("Madness??? This! Is! GENtle-Rust!");
        });

        // Open new windows
        let mut new_windows: Vec<Window> = self.new_windows.drain(..).collect();
        for window in new_windows.drain(..) {
            let id = format!("Viewport {}", self.viewport_id_counter);
            let id = ViewportId::from_hash_of(id);
            self.windows.insert(id, Arc::new(RwLock::new(window)));
            self.viewport_id_counter += 1;
        }

        // Close windows
        for id in self.windows_to_close.write().unwrap().drain(..) {
            self.windows.remove(&id);
        }

        // Show windows
        for (id, window) in self.windows.iter() {
            let id = id.to_owned();
            self.show_window(ctx, id, window.clone());
        }
    }
}
