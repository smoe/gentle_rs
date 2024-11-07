use crate::{main_area::MainArea, window_dna::WindowDna};
use eframe::egui::{self, CollapsingHeader};
use std::sync::Arc;

#[derive(Debug, Default, Clone)]
pub struct LeftPanel {}

impl LeftPanel {
    pub fn render(&mut self, window: &mut WindowDna, ui: &mut egui::Ui) {
        egui::ScrollArea::both().show(ui, |ui| {
            ui.heading("Project");
            // egui::Separator::default().horizontal();
            ui.vertical(|ui| {
                CollapsingHeader::new("DNA")
                    .default_open(true)
                    .show(ui, |ui| {
                        ui.vertical(|ui| {
                            for (i, dna) in window.dna_sequence().clone().iter().enumerate() {
                                let name = dna
                                    .read()
                                    .expect("DNA lock poisoned")
                                    .name()
                                    .as_ref()
                                    .unwrap_or(&format!("DNA #{i}"))
                                    .to_string();
                                let selected = true;
                                let button = egui::Button::new(name).selected(selected);
                                if ui.add(button).clicked() {
                                    window.set_main_area(Some(MainArea::new_dna(dna.clone())));
                                }
                            }
                        });
                    });
            });
        });
    }
}
