use std::sync::{Arc, Mutex};

use crate::{
    dna_sequence::DNAsequence,
    icons::{ICON_CIRCULAR_LINEAR, ICON_SHOW_MAP, ICON_SHOW_SEQUENCE},
};
use eframe::egui::{self, Vec2};

#[derive(Debug)]
pub struct MainAreaDna {
    dna: Arc<Mutex<DNAsequence>>,
    show_sequence: bool,
    show_map: bool,
}

impl MainAreaDna {
    pub fn new(dna: Arc<Mutex<DNAsequence>>) -> Self {
        Self {
            dna,
            show_sequence: true,
            show_map: true,
        }
    }

    pub fn dna(&self) -> &Arc<Mutex<DNAsequence>> {
        &self.dna
    }

    pub fn render(&mut self, ctx: &egui::Context, ui: &mut egui::Ui) {
        egui::TopBottomPanel::top("dna_top_buttons").show(ctx, |ui| {
            self.render_top_panel(ui);
        });

        if self.show_map {
            if self.show_sequence {
                let full_height = ui.available_rect_before_wrap().height();
                egui::TopBottomPanel::bottom("dna_sequence")
                    .resizable(true)
                    .max_height(full_height / 2.0)
                    .show(ctx, |ui| {
                        self.render_sequence(ui);
                    });
                egui::CentralPanel::default().show(ctx, |ui| {
                    self.render_middle(ui);
                });
            } else {
                egui::CentralPanel::default().show(ctx, |ui| {
                    self.render_dna_map(ui);
                });
            }
        } else {
            egui::CentralPanel::default().show(ctx, |ui| {
                self.render_sequence(ui);
            });
        }
    }

    pub fn render_top_panel(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            let button = egui::ImageButton::new(
                ICON_CIRCULAR_LINEAR
                    .clone()
                    .fit_to_fraction(Vec2::new(2.0, 2.0))
                    .rounding(5.0),
            );
            if ui.add(button).clicked() {
                let mut dna = self.dna.lock().unwrap();
                let is_circular = dna.is_circular();
                dna.set_circular(!is_circular);
            }

            let button = egui::ImageButton::new(ICON_SHOW_SEQUENCE.clone().rounding(5.0));
            if ui.add(button).clicked() {
                self.show_sequence = !self.show_sequence
            };

            let button = egui::ImageButton::new(ICON_SHOW_MAP.clone().rounding(5.0));
            if ui.add(button).clicked() {
                self.show_map = !self.show_map
            };
        });
    }

    pub fn render_sequence(&mut self, ui: &mut egui::Ui) {
        ui.label(self.dna.lock().unwrap().to_string());
    }

    pub fn render_features(&mut self, ui: &mut egui::Ui) {
        ui.heading(
            self.dna
                .lock()
                .unwrap()
                .name()
                .as_ref()
                .map(|s| s.as_str())
                .unwrap_or("<Unnamed DNA sequence>"),
        );
    }

    pub fn render_description(&mut self, ui: &mut egui::Ui) {
        ui.heading(self.dna.lock().unwrap().description().join("\n"));
    }

    pub fn render_dna_map(&mut self, ui: &mut egui::Ui) {
        if self.dna.lock().unwrap().is_circular() {
            ui.heading("Circular DNA");
        } else {
            ui.heading("Linear DNA");
        }
    }

    pub fn render_middle(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            ui.vertical(|ui| {
                self.render_features(ui);
                self.render_description(ui);
            });
            self.render_dna_map(ui);
        });
    }
}
