use crate::{
    dna_sequence::DNAsequence,
    icons::{ICON_CIRCULAR_LINEAR, ICON_SHOW_MAP, ICON_SHOW_SEQUENCE},
    render_dna_circular::RenderDnaCircular,
};
use eframe::egui::{self, Frame, PointerState, Sense, Vec2};
use std::sync::{Arc, Mutex};

#[derive(Debug)]
pub struct MainAreaDna {
    dna: Arc<Mutex<DNAsequence>>,
    circular_map: Option<RenderDnaCircular>,
    show_sequence: bool,
    show_map: bool,
}

impl MainAreaDna {
    pub fn new(dna: Arc<Mutex<DNAsequence>>) -> Self {
        Self {
            dna,
            circular_map: None,
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
                let frame = Frame::default();
                egui::CentralPanel::default().frame(frame).show(ctx, |ui| {
                    self.render_middle(ctx, ui);
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
                let mut dna = self.dna.lock().expect("DNA lock poisoned");
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
        ui.label(self.dna.lock().expect("DNA lock poisoned").to_string());
    }

    pub fn render_features(&mut self, ui: &mut egui::Ui) {
        ui.heading(
            self.dna
                .lock()
                .expect("DNA lock poisoned")
                .name()
                .as_ref()
                .map(|s| s.as_str())
                .unwrap_or("<Unnamed DNA sequence>"),
        );
    }

    pub fn render_description(&mut self, ui: &mut egui::Ui) {
        let description = self
            .dna
            .lock()
            .expect("DNA lock poisoned")
            .description()
            .join("\n");
        // description += "\nTHIS IS THE END";
        ui.heading(description);
    }

    fn is_circular(&self) -> bool {
        self.dna.lock().expect("DNA lock poisoned").is_circular()
    }

    pub fn render_dna_map(&mut self, ui: &mut egui::Ui) {
        if self.is_circular() {
            // TODO remove linear renderer
            if self.circular_map.is_none() {
                self.circular_map = Some(RenderDnaCircular::new(self.dna.clone()));
            }
            if let Some(renderer) = &mut self.circular_map {
                renderer.render(ui);
            }
        } else {
            self.circular_map = None;
            ui.heading("Linear DNA");
        }
    }

    pub fn render_middle(&mut self, ctx: &egui::Context, ui: &mut egui::Ui) {
        ui.with_layout(egui::Layout::left_to_right(egui::Align::Center), |ui| {
            Frame::none().show(ui, |ui| {
                ui.vertical(|ui| {
                    self.render_features(ui);
                    self.render_description(ui);
                });
            });

            Frame::none()
                .fill(egui::Color32::LIGHT_BLUE)
                .show(ui, |ui| {
                    self.render_dna_map(ui);
                });
        });

        if ui.response().interact(Sense::click()).clicked() {
            let pointer_state: PointerState = ctx.input(|i| i.pointer.to_owned());
            if let Some(dna_map) = &self.circular_map {
                dna_map.on_click(pointer_state);
            }
            // TODO linear
        }
    }
}
