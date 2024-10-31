use crate::{app::GENtleApp, main_area::MainArea};
use eframe::egui::{self, CollapsingHeader};
use std::sync::Arc;

#[derive(Debug, Default, Clone)]
pub struct LeftPanel {}

impl LeftPanel {
    pub fn render(&mut self, app: &mut GENtleApp, ui: &mut egui::Ui) {
        egui::ScrollArea::both().show(ui, |ui| {
            ui.heading("Project");
            // egui::Separator::default().horizontal();
            ui.vertical(|ui| {
                CollapsingHeader::new("DNA")
                    .default_open(true)
                    .show(ui, |ui| {
                        ui.vertical(|ui| {
                            for (i, dna) in app.dna_sequence().clone().iter().enumerate() {
                                let name = dna
                                    .lock()
                                    .expect("DNA lock poisoned")
                                    .name()
                                    .as_ref()
                                    .unwrap_or(&format!("DNA #{i}"))
                                    .to_string();
                                let selected = match &app.main_area() {
                                    Some(MainArea::Dna(ma)) => Arc::ptr_eq(dna, ma.dna()),
                                    _ => false,
                                };
                                let button = egui::Button::new(name).selected(selected);
                                if ui.add(button).clicked() {
                                    app.set_main_area(Some(MainArea::new_dna(dna.clone())));
                                }
                            }
                        });
                    });
            });
        });
    }
}
