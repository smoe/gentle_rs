use crate::{
    dna_sequence::DNAsequence,
    icons::{ICON_CIRCULAR_LINEAR, ICON_SHOW_MAP, ICON_SHOW_SEQUENCE},
    render_dna::RenderDnaEnum,
};
use eframe::egui::{self, Frame, PointerState, Sense, Vec2};
use std::{
    collections::HashMap,
    sync::{Arc, RwLock},
};

#[derive(Debug)]
pub struct MainAreaDna {
    dna: Arc<RwLock<DNAsequence>>,
    map_dna: RenderDnaEnum,
    show_sequence: bool,
    show_map: bool,
}

impl MainAreaDna {
    pub fn new(dna: Arc<RwLock<DNAsequence>>) -> Self {
        Self {
            dna: dna.clone(),
            map_dna: RenderDnaEnum::new(dna),
            show_sequence: true,
            show_map: true,
        }
    }

    pub fn dna(&self) -> &Arc<RwLock<DNAsequence>> {
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
                let mut dna = self.dna.write().expect("DNA lock poisoned");
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
        ui.label(self.dna.read().expect("DNA lock poisoned").to_string());
    }

    fn get_selected_feature_id(&self) -> Option<usize> {
        self.map_dna.get_selected_feature_id()
    }

    pub fn render_features(&mut self, ui: &mut egui::Ui) {
        ui.heading(
            self.dna
                .read()
                .expect("DNA lock poisoned")
                .name()
                .as_ref()
                .map(|s| s.as_str())
                .unwrap_or("<Unnamed DNA sequence>"),
        );

        let selected_id = self.get_selected_feature_id();
        let typed_features = self
            .dna
            .read()
            .expect("DNA lock poisoned")
            .features()
            .iter()
            .enumerate()
            .map(|(id, feature)| {
                let kind = feature.kind.to_string();

                (kind, id)
            })
            .collect::<Vec<_>>();
        let mut grouped_features: HashMap<String, Vec<usize>> = HashMap::new();
        for (kind, id) in typed_features {
            grouped_features.entry(kind).or_default().push(id);
        }
        let mut group_keys = grouped_features.keys().collect::<Vec<_>>();
        group_keys.sort();
        for kind in group_keys {
            let ids = grouped_features.get(kind).unwrap();
            ui.collapsing(kind, |ui| {
                for id in ids {
                    let name = match &self
                        .dna
                        .read()
                        .expect("DNA lock poisoned")
                        .features()
                        .get(*id)
                    {
                        Some(feature) => RenderDnaEnum::feature_name(feature),
                        None => continue,
                    };
                    let selected = selected_id == Some(*id);
                    ui.horizontal(|ui| {
                        let button = egui::Button::new(name).selected(selected);
                        if ui.add(button).clicked() {
                            self.map_dna.select_feature(Some(*id));
                        }
                    });
                }
            });
        }
    }

    fn get_sequence_description(&self) -> String {
        self.dna
            .read()
            .expect("DNA lock poisoned")
            .description()
            .join("\n")
    }

    pub fn render_description(&mut self, ui: &mut egui::Ui) {
        egui::ScrollArea::vertical().show(ui, |ui| {
            ui.set_min_height(150.0);

            match self.get_selected_feature_id() {
                Some(id) => {
                    let feature = self
                        .dna
                        .read()
                        .expect("DNA lock poisoned")
                        .features()
                        .get(id)
                        .unwrap()
                        .to_owned(); // Temporary copy

                    let name = RenderDnaEnum::feature_name(&feature);
                    ui.heading(name);
                    let desc = &match feature.location.find_bounds() {
                        Ok((from, to)) => format!("{from}..{to}"),
                        Err(_) => String::new(),
                    };
                    ui.monospace(desc);
                }
                None => {
                    let description = self.get_sequence_description();
                    ui.heading(description);
                }
            };
        });
    }

    fn is_circular(&self) -> bool {
        self.dna.read().expect("DNA lock poisoned").is_circular()
    }

    pub fn render_dna_map(&mut self, ui: &mut egui::Ui) {
        if self.is_circular() != self.map_dna.is_circular() {
            self.map_dna = RenderDnaEnum::new(self.dna.clone());
        }
        self.map_dna.render(ui);
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

        if ui.response().clicked() {
            // .interact(Sense::click())
            let pointer_state: PointerState = ctx.input(|i| i.pointer.to_owned());
            println!("{pointer_state:?}");

            self.map_dna.on_click(pointer_state);
        }
    }
}
