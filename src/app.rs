use std::{
    collections::HashMap,
    sync::{Arc, RwLock},
};

use crate::{
    dna_sequence::{self, DNAsequence},
    engine::{Engine, GentleEngine, Operation},
    window::Window,
    TRANSLATIONS,
};
use anyhow::{anyhow, Result};
use eframe::egui::{self, menu, Pos2, Ui, ViewportId};

pub struct GENtleApp {
    engine: Arc<RwLock<GentleEngine>>,
    new_windows: Vec<Window>,
    windows: HashMap<ViewportId, Arc<RwLock<Window>>>,
    windows_to_close: Arc<RwLock<Vec<ViewportId>>>,
    viewport_id_counter: usize,
    update_has_run_before: bool,
}

impl Default for GENtleApp {
    fn default() -> Self {
        Self {
            engine: Arc::new(RwLock::new(GentleEngine::new())),
            new_windows: vec![],
            windows: HashMap::new(),
            windows_to_close: Arc::new(RwLock::new(vec![])),
            viewport_id_counter: 0,
            update_has_run_before: false,
        }
    }
}

impl GENtleApp {
    pub fn new() -> Self {
        let mut ret = Self::default();

        // Load test sequences
        ret.open_new_window_from_file("test_files/pGEX-3X.gb");   // GenBank entry of complete circular genome
        ret.open_new_window_from_file("test_files/pGEX_3X.fa");   // Sequence only
        ret.open_new_window_from_file("test_files/tp73.ncbi.gb"); // GenBank entry of complex human gene (linear, fragment of genome)

        ret
    }

    fn load_dna_from_genbank_file(filename: &str) -> Result<DNAsequence> {
        let dna = dna_sequence::DNAsequence::from_genbank_file(filename)?
            .pop()
            .ok_or_else(|| anyhow!("Could not read GenBank file {filename}"))?;
        Ok(dna)
    }

    fn load_dna_from_fasta_file(filename: &str) -> Result<DNAsequence> {
        let dna = dna_sequence::DNAsequence::from_fasta_file(filename)?
            .pop()
            .ok_or_else(|| anyhow!("Could not read fasta file {filename}"))?;
        Ok(dna)
    }

    fn new_dna_window(&mut self, seq_id: String, dna: DNAsequence) {
        self.new_windows
            .push(Window::new_dna(dna, seq_id, self.engine.clone()));
    }

    pub fn load_from_file(path: &str) -> Result<DNAsequence> {
        if let Ok(dna) = Self::load_dna_from_genbank_file(path) {
            Ok(dna)
        } else if let Ok(dna) = Self::load_dna_from_fasta_file(path) {
            Ok(dna)
        } else {
            Err(anyhow!("Could not load file '{path}'"))
        }
    }

    fn open_new_window_from_file(&mut self, path: &str) {
        let op = Operation::LoadFile {
            path: path.to_string(),
            as_id: None,
        };
        let load_result = {
            let mut engine = self.engine.write().unwrap();
            engine.apply(op)
        };
        if let Ok(result) = load_result {
            if let Some(seq_id) = result.created_seq_ids.first() {
                let seq_id = seq_id.to_string();
                let dna = self
                    .engine
                    .read()
                    .unwrap()
                    .state()
                    .sequences
                    .get(&seq_id)
                    .cloned();
                if let Some(dna) = dna {
                    self.new_dna_window(seq_id, dna);
                    return;
                }
            }
            // TODO warning: load op ran but no created sequence
        } else {
            // TODO error, could not load file through engine
        }
    }

    pub fn render_menu_bar(&mut self, ui: &mut Ui) {
        menu::bar(ui, |ui| {
            ui.menu_button(TRANSLATIONS.get("m_file"), |ui| {
                if ui.button(TRANSLATIONS.get("m_open")).clicked() {
                    if let Some(path) = rfd::FileDialog::new().pick_file() {
                        if let Some(path) = Some(path.display().to_string()) {
                            self.open_new_window_from_file(&path);
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
