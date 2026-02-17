use crate::{
    dna_display::DnaDisplay, dna_sequence::DNAsequence, icons::*, render_dna::RenderDna,
    render_sequence::RenderSequence,
    engine::{DisplayTarget, Engine, GentleEngine, Operation},
};
use eframe::egui::{self, Frame, PointerState, Vec2};
use std::{
    collections::HashMap,
    sync::{Arc, RwLock},
};

#[derive(Clone, Debug)]
pub struct MainAreaDna {
    dna: Arc<RwLock<DNAsequence>>,
    dna_display: Arc<RwLock<DnaDisplay>>,
    engine: Option<Arc<RwLock<GentleEngine>>>,
    seq_id: Option<String>,
    map_dna: RenderDna,
    map_sequence: RenderSequence,
    show_sequence: bool, // TODO move to DnaDisplay
    show_map: bool,      // TODO move to DnaDisplay
}

impl MainAreaDna {
    pub fn new(
        dna: DNAsequence,
        seq_id: Option<String>,
        engine: Option<Arc<RwLock<GentleEngine>>>,
    ) -> Self {
        let dna = Arc::new(RwLock::new(dna));
        let dna_display = Arc::new(RwLock::new(DnaDisplay::default()));
        let mut ret = Self {
            dna: dna.clone(),
            dna_display: dna_display.clone(),
            engine,
            seq_id,
            map_dna: RenderDna::new(dna.clone(), dna_display.clone()),
            map_sequence: RenderSequence::new_single_sequence(dna, dna_display),
            show_sequence: true,
            show_map: true,
        };
        ret.sync_from_engine_display();
        ret
    }

    pub fn dna(&self) -> &Arc<RwLock<DNAsequence>> {
        &self.dna
    }

    pub fn render(&mut self, ctx: &egui::Context, ui: &mut egui::Ui) {
        self.sync_from_engine_display();
        egui::TopBottomPanel::top("dna_top_buttons").show(ctx, |ui| {
            self.render_top_panel(ui);
        });

        if self.show_map {
            if self.show_sequence {
                let full_height = ui.available_rect_before_wrap().height();
                egui::TopBottomPanel::bottom("dna_sequence")
                    .resizable(true)
                    .default_height(full_height / 2.0)
                    .max_height(full_height / 2.0)
                    .min_height(full_height / 4.0)
                    .show(ctx, |ui| {
                        self.render_sequence(ui);
                    });
                let frame = Frame::default();
                egui::CentralPanel::default().frame(frame).show(ctx, |ui| {
                    self.render_middle(ctx, ui);
                });
            } else {
                egui::CentralPanel::default().show(ctx, |ui| {
                    self.update_dna_map();
                    ui.add(self.map_dna.to_owned());
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
            let response = ui
                .add(button)
                .on_hover_text("Toggle DNA map topology: circular <-> linear");
            if response.clicked() {
                if let (Some(engine), Some(seq_id)) = (&self.engine, &self.seq_id) {
                    let new_circular = !self.dna.read().expect("DNA lock poisoned").is_circular();
                    let op = Operation::SetTopology {
                        seq_id: seq_id.clone(),
                        circular: new_circular,
                    };
                    if engine.write().expect("Engine lock poisoned").apply(op).is_ok() {
                        if let Some(updated) = engine
                            .read()
                            .expect("Engine lock poisoned")
                            .state()
                            .sequences
                            .get(seq_id)
                            .cloned()
                        {
                            *self.dna.write().expect("DNA lock poisoned") = updated;
                        }
                    }
                } else {
                    let mut dna = self.dna.write().expect("DNA lock poisoned");
                    let is_circular = dna.is_circular();
                    dna.set_circular(!is_circular);
                }
            }

            let button =
                egui::ImageButton::new(ICON_SHOW_SEQUENCE.clone().rounding(5.0)).selected(self.show_sequence);
            let response = ui
                .add(button)
                .on_hover_text("Show or hide the sequence text panel");
            if response.clicked() {
                self.show_sequence = !self.show_sequence;
                self.set_display_visibility(DisplayTarget::SequencePanel, self.show_sequence);
            };

            let button =
                egui::ImageButton::new(ICON_SHOW_MAP.clone().rounding(5.0)).selected(self.show_map);
            let response = ui
                .add(button)
                .on_hover_text("Show or hide the DNA map panel");
            if response.clicked() {
                self.show_map = !self.show_map;
                self.set_display_visibility(DisplayTarget::MapPanel, self.show_map);
            };

            let features_active = self
                .dna_display
                .read()
                .expect("DNA display lock poisoned")
                .show_features();
            let button =
                egui::ImageButton::new(ICON_FEATURES.clone().rounding(5.0)).selected(features_active);
            let response = ui
                .add(button)
                .on_hover_text("Show or hide annotated sequence features");
            if response.clicked() {
                let visible = {
                    let display = self.dna_display.read().expect("DNA display lock poisoned");
                    !display.show_features()
                };
                self.dna_display
                    .write()
                    .expect("DNA display lock poisoned")
                    .set_show_features(visible);
                self.set_display_visibility(DisplayTarget::Features, visible);
            };

            let re_active = self
                .dna_display
                .read()
                .expect("DNA display lock poisoned")
                .show_restriction_enzyme_sites();
            let button = egui::ImageButton::new(ICON_RESTRICTION_ENZYMES.clone().rounding(5.0))
                .selected(re_active);
            let response = ui
                .add(button)
                .on_hover_text("Show or hide restriction enzyme cut sites");
            if response.clicked() {
                let visible = {
                    let display = self.dna_display.read().expect("DNA display lock poisoned");
                    !display.show_restriction_enzyme_sites()
                };
                self.dna_display
                    .write()
                    .expect("DNA display lock poisoned")
                    .set_show_restriction_enzyme_sites(visible);
                self.set_display_visibility(DisplayTarget::RestrictionEnzymes, visible);
            };

            let gc_active = self
                .dna_display
                .read()
                .expect("DNA display lock poisoned")
                .show_gc_contents();
            let button =
                egui::ImageButton::new(ICON_GC_CONTENT.clone().rounding(5.0)).selected(gc_active);
            let response = ui
                .add(button)
                .on_hover_text("Show or hide GC-content visualization");
            if response.clicked() {
                let visible = {
                    let display = self.dna_display.read().expect("DNA display lock poisoned");
                    !display.show_gc_contents()
                };
                self.dna_display
                    .write()
                    .expect("DNA display lock poisoned")
                    .set_show_gc_contents(visible);
                self.set_display_visibility(DisplayTarget::GcContents, visible);
            };

            let orf_active = self
                .dna_display
                .read()
                .expect("DNA display lock poisoned")
                .show_open_reading_frames();
            let button = egui::ImageButton::new(ICON_OPEN_READING_FRAMES.clone().rounding(5.0))
                .selected(orf_active);
            let response = ui
                .add(button)
                .on_hover_text("Show or hide predicted open reading frames (ORFs)");
            if response.clicked() {
                let visible = {
                    let display = self.dna_display.read().expect("DNA display lock poisoned");
                    !display.show_open_reading_frames()
                };
                self.dna_display
                    .write()
                    .expect("DNA display lock poisoned")
                    .set_show_open_reading_frames(visible);
                self.set_display_visibility(DisplayTarget::OpenReadingFrames, visible);
            };

            let methylation_active = self
                .dna_display
                .read()
                .expect("DNA display lock poisoned")
                .show_methylation_sites();
            let button = egui::ImageButton::new(ICON_METHYLATION_SITES.clone().rounding(5.0))
                .selected(methylation_active);
            let response = ui
                .add(button)
                .on_hover_text("Show or hide methylation-site markers");
            if response.clicked() {
                let visible = {
                    let display = self.dna_display.read().expect("DNA display lock poisoned");
                    !display.show_methylation_sites()
                };
                self.dna_display
                    .write()
                    .expect("DNA display lock poisoned")
                    .set_show_methylation_sites(visible);
                self.set_display_visibility(DisplayTarget::MethylationSites, visible);
            };
        });
    }

    fn set_display_visibility(&self, target: DisplayTarget, visible: bool) {
        if let Some(engine) = &self.engine {
            let _ = engine
                .write()
                .expect("Engine lock poisoned")
                .apply(Operation::SetDisplayVisibility { target, visible });
        }
    }

    fn sync_from_engine_display(&mut self) {
        let Some(engine) = &self.engine else {
            return;
        };
        let settings = engine.read().expect("Engine lock poisoned").state().display.clone();
        self.show_sequence = settings.show_sequence_panel;
        self.show_map = settings.show_map_panel;
        let mut display = self.dna_display.write().expect("DNA display lock poisoned");
        display.set_show_features(settings.show_features);
        display.set_show_restriction_enzyme_sites(settings.show_restriction_enzymes);
        display.set_show_gc_contents(settings.show_gc_contents);
        display.set_show_open_reading_frames(settings.show_open_reading_frames);
        display.set_show_methylation_sites(settings.show_methylation_sites);
    }

    pub fn render_sequence(&mut self, ui: &mut egui::Ui) {
        self.map_sequence.render(ui);
    }

    fn get_selected_feature_id(&self) -> Option<usize> {
        self.map_dna.get_selected_feature_id()
    }

    pub fn sequence_name(&self) -> Option<String> {
        self.dna
            .read()
            .expect("DNA lock poisoned")
            .name()
            .to_owned()
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
                        Some(feature) => RenderDna::feature_name(feature),
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

                    let name = RenderDna::feature_name(&feature);
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

    pub fn update_dna_map(&mut self) {
        if self.is_circular() != self.map_dna.is_circular() {
            self.map_dna = RenderDna::new(self.dna.clone(), self.dna_display.clone());
        }
    }

    pub fn render_middle(&mut self, ctx: &egui::Context, ui: &mut egui::Ui) {
        ui.with_layout(egui::Layout::left_to_right(egui::Align::Center), |ui| {
            self.update_dna_map();

            Frame::none().show(ui, |ui| {
                ui.vertical(|ui| {
                    self.render_features(ui);
                    self.render_description(ui);
                });
            });

            ui.separator();

            let response = ui.add(self.map_dna.to_owned());

            if response.hovered() {
                let pointer_state: PointerState = ctx.input(|i| i.pointer.to_owned());
                self.map_dna.on_hover(pointer_state);
            }

            if response.clicked() {
                let pointer_state: PointerState = ctx.input(|i| i.pointer.to_owned());
                self.map_dna.on_click(pointer_state);
            }

            if response.double_clicked() {
                let pointer_state: PointerState = ctx.input(|i| i.pointer.to_owned());
                self.map_dna.on_double_click(pointer_state);
            }
        });
    }
}
