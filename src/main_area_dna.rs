use crate::{
    dna_display::DnaDisplay,
    dna_sequence::DNAsequence,
    engine::{
        DisplayTarget, Engine, GentleEngine, LigationProtocol, Operation, PcrPrimerSpec,
        SnpMutationSpec,
    },
    icons::*,
    render_dna::RenderDna,
    render_sequence::RenderSequence,
};
use eframe::egui::{self, Frame, PointerState, Vec2};
use serde::{Deserialize, Serialize};
use std::{
    collections::HashMap,
    sync::{Arc, RwLock},
    time::Instant,
};

#[derive(Clone, Debug, Serialize, Deserialize)]
struct EngineOpsUiState {
    show_engine_ops: bool,
    digest_enzymes_text: String,
    digest_prefix_text: String,
    merge_inputs_text: String,
    ligation_inputs_text: String,
    ligation_output_prefix: String,
    ligation_protocol_sticky: bool,
    ligation_unique: bool,
    mw_inputs_text: String,
    mw_min_bp: String,
    mw_max_bp: String,
    mw_error: String,
    mw_unique: bool,
    pcr_forward: String,
    pcr_reverse: String,
    pcr_unique: bool,
    pcr_adv_forward: String,
    pcr_adv_reverse: String,
    pcr_adv_anneal_len: String,
    pcr_adv_max_mismatch: String,
    pcr_adv_3prime_exact: String,
    pcr_mut_position: String,
    pcr_mut_ref: String,
    pcr_mut_alt: String,
}

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
    show_engine_ops: bool,
    digest_enzymes_text: String,
    digest_prefix_text: String,
    merge_inputs_text: String,
    ligation_inputs_text: String,
    ligation_output_prefix: String,
    ligation_protocol_sticky: bool,
    ligation_unique: bool,
    mw_inputs_text: String,
    mw_min_bp: String,
    mw_max_bp: String,
    mw_error: String,
    mw_unique: bool,
    pcr_forward: String,
    pcr_reverse: String,
    pcr_unique: bool,
    pcr_adv_forward: String,
    pcr_adv_reverse: String,
    pcr_adv_anneal_len: String,
    pcr_adv_max_mismatch: String,
    pcr_adv_3prime_exact: String,
    pcr_mut_position: String,
    pcr_mut_ref: String,
    pcr_mut_alt: String,
    op_status: String,
    op_error_popup: Option<String>,
    last_created_seq_ids: Vec<String>,
    heartbeat_frame: u64,
    heartbeat_last_ms: Option<u128>,
    heartbeat_delta_ms: Option<u128>,
}

impl MainAreaDna {
    pub fn new(
        dna: DNAsequence,
        seq_id: Option<String>,
        engine: Option<Arc<RwLock<GentleEngine>>>,
    ) -> Self {
        let seq_id_for_defaults = seq_id.clone();
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
            show_engine_ops: false,
            digest_enzymes_text: "BamHI,EcoRI".to_string(),
            digest_prefix_text: "digest".to_string(),
            merge_inputs_text: seq_id_for_defaults.clone().unwrap_or_default(),
            ligation_inputs_text: seq_id_for_defaults.clone().unwrap_or_default(),
            ligation_output_prefix: "ligation".to_string(),
            ligation_protocol_sticky: true,
            ligation_unique: false,
            mw_inputs_text: seq_id_for_defaults.clone().unwrap_or_default(),
            mw_min_bp: "100".to_string(),
            mw_max_bp: "1000".to_string(),
            mw_error: "0.10".to_string(),
            mw_unique: false,
            pcr_forward: String::new(),
            pcr_reverse: String::new(),
            pcr_unique: false,
            pcr_adv_forward: String::new(),
            pcr_adv_reverse: String::new(),
            pcr_adv_anneal_len: "18".to_string(),
            pcr_adv_max_mismatch: "1".to_string(),
            pcr_adv_3prime_exact: "8".to_string(),
            pcr_mut_position: "0".to_string(),
            pcr_mut_ref: "A".to_string(),
            pcr_mut_alt: "G".to_string(),
            op_status: String::new(),
            op_error_popup: None,
            last_created_seq_ids: vec![],
            heartbeat_frame: 0,
            heartbeat_last_ms: None,
            heartbeat_delta_ms: None,
        };
        ret.sync_from_engine_display();
        ret.load_engine_ops_state();
        ret
    }

    pub fn dna(&self) -> &Arc<RwLock<DNAsequence>> {
        &self.dna
    }

    pub fn set_pool_context(&mut self, pool_seq_ids: Vec<String>) {
        self.last_created_seq_ids = pool_seq_ids;
        self.show_engine_ops = true;
        self.op_status = "Opened from lineage pool node".to_string();
    }

    pub fn render(&mut self, ctx: &egui::Context, ui: &mut egui::Ui) {
        self.update_heartbeat();
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
        self.render_error_popup(ctx);
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
                    if engine
                        .write()
                        .expect("Engine lock poisoned")
                        .apply(op)
                        .is_ok()
                    {
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

            let button = egui::ImageButton::new(ICON_SHOW_SEQUENCE.clone().rounding(5.0))
                .selected(self.show_sequence);
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
            let button = egui::ImageButton::new(ICON_FEATURES.clone().rounding(5.0))
                .selected(features_active);
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

            ui.separator();
            if ui
                .button("Rev")
                .on_hover_text("Create a reversed branch sequence")
                .clicked()
            {
                self.apply_sequence_derivation(Operation::Reverse {
                    input: self.seq_id.clone().unwrap_or_default(),
                    output_id: None,
                });
            }
            if ui
                .button("Comp")
                .on_hover_text("Create a complemented branch sequence")
                .clicked()
            {
                self.apply_sequence_derivation(Operation::Complement {
                    input: self.seq_id.clone().unwrap_or_default(),
                    output_id: None,
                });
            }
            if ui
                .button("RevComp")
                .on_hover_text("Create a reverse-complement branch sequence")
                .clicked()
            {
                self.apply_sequence_derivation(Operation::ReverseComplement {
                    input: self.seq_id.clone().unwrap_or_default(),
                    output_id: None,
                });
            }
            if ui
                .button("Branch")
                .on_hover_text("Create an unchanged branch copy for alternative workflows")
                .clicked()
            {
                self.apply_sequence_derivation(Operation::Branch {
                    input: self.seq_id.clone().unwrap_or_default(),
                    output_id: None,
                });
            }
            ui.separator();
            if ui
                .button("Engine Ops")
                .on_hover_text("Open strict engine operation controls")
                .clicked()
            {
                self.show_engine_ops = !self.show_engine_ops;
                self.save_engine_ops_state();
            }
            if !self.op_status.is_empty() {
                ui.monospace(&self.op_status);
            }
        });

        if self.show_engine_ops {
            ui.separator();
            ui.collapsing("Strict Engine Operations", |ui| {
                let heartbeat_now = self.heartbeat_last_ms.unwrap_or(0);
                let heartbeat_delta = self.heartbeat_delta_ms.unwrap_or(0);
                ui.monospace(format!(
                    "heartbeat frame={} now_ms={} delta_ms={}",
                    self.heartbeat_frame, heartbeat_now, heartbeat_delta
                ));
                ui.label("IDs are comma-separated sequence IDs.");
                let template_seq_id = self.seq_id.clone().unwrap_or_default();

                ui.horizontal(|ui| {
                    ui.label(format!("Digest template {}", template_seq_id));
                    ui.label("enzymes");
                    ui.text_edit_singleline(&mut self.digest_enzymes_text);
                    ui.label("prefix");
                    ui.text_edit_singleline(&mut self.digest_prefix_text);
                    if ui.button("Digest").clicked() {
                        let enzymes = self
                            .digest_enzymes_text
                            .split(',')
                            .map(|s| s.trim())
                            .filter(|s| !s.is_empty())
                            .map(|s| s.to_string())
                            .collect::<Vec<_>>();
                        if template_seq_id.is_empty() {
                            self.op_status = "No active template sequence".to_string();
                        } else if enzymes.is_empty() {
                            self.op_status = "Provide at least one enzyme".to_string();
                        } else {
                            self.apply_operation_with_feedback(Operation::Digest {
                                input: template_seq_id.clone(),
                                enzymes,
                                output_prefix: Some(self.digest_prefix_text.clone()),
                            });
                        }
                    }
                });

                ui.horizontal(|ui| {
                    ui.label("Merge inputs");
                    ui.text_edit_singleline(&mut self.merge_inputs_text);
                    if ui.button("Merge").clicked() {
                        let inputs = Self::parse_ids(&self.merge_inputs_text);
                        self.apply_operation_with_feedback(Operation::MergeContainers {
                            inputs,
                            output_prefix: Some("merged".to_string()),
                        });
                    }
                });

                ui.horizontal(|ui| {
                    ui.label("Ligation inputs");
                    ui.text_edit_singleline(&mut self.ligation_inputs_text);
                });
                ui.horizontal(|ui| {
                    ui.checkbox(&mut self.ligation_protocol_sticky, "Sticky (off = Blunt)");
                    ui.checkbox(&mut self.ligation_unique, "Unique");
                    ui.label("Prefix");
                    ui.text_edit_singleline(&mut self.ligation_output_prefix);
                    if ui.button("Ligate").clicked() {
                        let inputs = Self::parse_ids(&self.ligation_inputs_text);
                        let protocol = if self.ligation_protocol_sticky {
                            LigationProtocol::Sticky
                        } else {
                            LigationProtocol::Blunt
                        };
                        self.apply_operation_with_feedback(Operation::Ligation {
                            inputs,
                            circularize_if_possible: false,
                            output_id: None,
                            protocol,
                            output_prefix: Some(self.ligation_output_prefix.clone()),
                            unique: Some(self.ligation_unique),
                        });
                    }
                });

                ui.horizontal(|ui| {
                    ui.label("MW inputs");
                    ui.text_edit_singleline(&mut self.mw_inputs_text);
                });
                ui.horizontal(|ui| {
                    ui.label("min");
                    ui.text_edit_singleline(&mut self.mw_min_bp);
                    ui.label("max");
                    ui.text_edit_singleline(&mut self.mw_max_bp);
                    ui.label("error");
                    ui.text_edit_singleline(&mut self.mw_error);
                    ui.checkbox(&mut self.mw_unique, "Unique");
                    if ui.button("Filter MW").clicked() {
                        if let (Ok(min_bp), Ok(max_bp), Ok(error)) = (
                            self.mw_min_bp.parse::<usize>(),
                            self.mw_max_bp.parse::<usize>(),
                            self.mw_error.parse::<f64>(),
                        ) {
                            self.apply_operation_with_feedback(
                                Operation::FilterByMolecularWeight {
                                    inputs: Self::parse_ids(&self.mw_inputs_text),
                                    min_bp,
                                    max_bp,
                                    error,
                                    unique: self.mw_unique,
                                    output_prefix: Some("mw".to_string()),
                                },
                            );
                        } else {
                            self.op_status = "Invalid MW parameters (min/max/error)".to_string();
                        }
                    }
                });

                ui.horizontal(|ui| {
                    ui.label("PCR fwd");
                    ui.text_edit_singleline(&mut self.pcr_forward);
                    ui.label("rev");
                    ui.text_edit_singleline(&mut self.pcr_reverse);
                    ui.checkbox(&mut self.pcr_unique, "Unique");
                    if ui.button("PCR").clicked() {
                        let template = self.seq_id.clone().unwrap_or_default();
                        self.apply_operation_with_feedback(Operation::Pcr {
                            template,
                            forward_primer: self.pcr_forward.clone(),
                            reverse_primer: self.pcr_reverse.clone(),
                            output_id: None,
                            unique: Some(self.pcr_unique),
                        });
                    }
                });

                ui.horizontal(|ui| {
                    ui.label("PCR adv fwd");
                    ui.text_edit_singleline(&mut self.pcr_adv_forward);
                    ui.label("rev");
                    ui.text_edit_singleline(&mut self.pcr_adv_reverse);
                });
                ui.horizontal(|ui| {
                    ui.label("anneal");
                    ui.text_edit_singleline(&mut self.pcr_adv_anneal_len);
                    ui.label("mm");
                    ui.text_edit_singleline(&mut self.pcr_adv_max_mismatch);
                    ui.label("3' exact");
                    ui.text_edit_singleline(&mut self.pcr_adv_3prime_exact);
                    if ui.button("PCR Adv").clicked() {
                        if let (Ok(anneal_len), Ok(max_mismatches), Ok(exact3)) = (
                            self.pcr_adv_anneal_len.parse::<usize>(),
                            self.pcr_adv_max_mismatch.parse::<usize>(),
                            self.pcr_adv_3prime_exact.parse::<usize>(),
                        ) {
                            let template = self.seq_id.clone().unwrap_or_default();
                            let primer = |sequence: String| PcrPrimerSpec {
                                sequence,
                                anneal_len: Some(anneal_len),
                                max_mismatches: Some(max_mismatches),
                                require_3prime_exact_bases: Some(exact3),
                                library_mode: None,
                                max_variants: None,
                                sample_seed: None,
                            };
                            self.apply_operation_with_feedback(Operation::PcrAdvanced {
                                template,
                                forward_primer: primer(self.pcr_adv_forward.clone()),
                                reverse_primer: primer(self.pcr_adv_reverse.clone()),
                                output_id: None,
                                unique: Some(self.pcr_unique),
                            });
                        } else {
                            self.op_status = "Invalid advanced PCR parameters".to_string();
                        }
                    }
                });

                ui.horizontal(|ui| {
                    ui.label("SNP pos");
                    ui.text_edit_singleline(&mut self.pcr_mut_position);
                    ui.label("ref");
                    ui.text_edit_singleline(&mut self.pcr_mut_ref);
                    ui.label("alt");
                    ui.text_edit_singleline(&mut self.pcr_mut_alt);
                    if ui.button("PCR Mut").clicked() {
                        if let (Ok(pos), Ok(anneal_len), Ok(max_mismatches), Ok(exact3)) = (
                            self.pcr_mut_position.parse::<usize>(),
                            self.pcr_adv_anneal_len.parse::<usize>(),
                            self.pcr_adv_max_mismatch.parse::<usize>(),
                            self.pcr_adv_3prime_exact.parse::<usize>(),
                        ) {
                            let template = self.seq_id.clone().unwrap_or_default();
                            let primer = |sequence: String| PcrPrimerSpec {
                                sequence,
                                anneal_len: Some(anneal_len),
                                max_mismatches: Some(max_mismatches),
                                require_3prime_exact_bases: Some(exact3),
                                library_mode: None,
                                max_variants: None,
                                sample_seed: None,
                            };
                            self.apply_operation_with_feedback(Operation::PcrMutagenesis {
                                template,
                                forward_primer: primer(self.pcr_adv_forward.clone()),
                                reverse_primer: primer(self.pcr_adv_reverse.clone()),
                                mutations: vec![SnpMutationSpec {
                                    zero_based_position: pos,
                                    reference: self.pcr_mut_ref.clone(),
                                    alternate: self.pcr_mut_alt.clone(),
                                }],
                                output_id: None,
                                unique: Some(self.pcr_unique),
                                require_all_mutations: Some(true),
                            });
                        } else {
                            self.op_status = "Invalid mutagenesis parameters".to_string();
                        }
                    }
                });

                ui.separator();
                if ui
                    .button("Preset: Digest -> Merge -> Sticky Ligation")
                    .on_hover_text("Runs strict 3-step workflow from active template")
                    .clicked()
                {
                    self.run_preset_digest_merge_ligate();
                }

                if self.last_created_seq_ids.len() > 1 {
                    ui.separator();
                    self.render_pool_distribution(ui);
                }

                self.save_engine_ops_state();
            });
        }
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
        let settings = engine
            .read()
            .expect("Engine lock poisoned")
            .state()
            .display
            .clone();
        self.show_sequence = settings.show_sequence_panel;
        self.show_map = settings.show_map_panel;
        let mut display = self.dna_display.write().expect("DNA display lock poisoned");
        display.set_show_features(settings.show_features);
        display.set_show_restriction_enzyme_sites(settings.show_restriction_enzymes);
        display.set_show_gc_contents(settings.show_gc_contents);
        display.set_show_open_reading_frames(settings.show_open_reading_frames);
        display.set_show_methylation_sites(settings.show_methylation_sites);
    }

    fn apply_sequence_derivation(&mut self, op: Operation) {
        let Some(engine) = &self.engine else {
            return;
        };
        let Ok(result) = engine.write().expect("Engine lock poisoned").apply(op) else {
            return;
        };
        let Some(new_seq_id) = result.created_seq_ids.first() else {
            return;
        };
        let Some(new_dna) = engine
            .read()
            .expect("Engine lock poisoned")
            .state()
            .sequences
            .get(new_seq_id)
            .cloned()
        else {
            return;
        };

        self.seq_id = Some(new_seq_id.clone());
        *self.dna.write().expect("DNA lock poisoned") = new_dna;
    }

    fn apply_operation_with_feedback(&mut self, op: Operation) {
        let Some(engine) = &self.engine else {
            self.op_status = "No engine attached".to_string();
            return;
        };
        let started = Instant::now();
        match engine.write().expect("Engine lock poisoned").apply(op) {
            Ok(result) => {
                self.last_created_seq_ids = result.created_seq_ids.clone();
                // Only auto-switch current view if exactly one sequence is produced.
                // Multi-product ops (digest/ligation/merge) can trigger heavy redraw paths;
                // keep the current view stable and let the user select products explicitly.
                if result.created_seq_ids.len() == 1 {
                    if let Some(new_seq_id) = result.created_seq_ids.first() {
                        if let Some(new_dna) = engine
                            .read()
                            .expect("Engine lock poisoned")
                            .state()
                            .sequences
                            .get(new_seq_id)
                            .cloned()
                        {
                            self.seq_id = Some(new_seq_id.clone());
                            *self.dna.write().expect("DNA lock poisoned") = new_dna;
                        }
                    }
                }
                let elapsed_ms = started.elapsed().as_millis();
                let created = if result.created_seq_ids.is_empty() {
                    "-".to_string()
                } else {
                    result.created_seq_ids.join(", ")
                };
                let warnings = if result.warnings.is_empty() {
                    "-".to_string()
                } else {
                    result.warnings.join(" | ")
                };
                let messages = if result.messages.is_empty() {
                    "-".to_string()
                } else {
                    result.messages.join(" | ")
                };
                self.op_status = format!(
                    "ok in {} ms\ncreated: {}\ncounts: {} created, {} changed\nwarnings: {}\nmessages: {}",
                    elapsed_ms,
                    created,
                    result.created_seq_ids.len(),
                    result.changed_seq_ids.len(),
                    warnings,
                    messages
                );
                self.op_error_popup = None;
            }
            Err(e) => {
                self.last_created_seq_ids.clear();
                let elapsed_ms = started.elapsed().as_millis();
                self.op_status = format!("error after {} ms: {}", elapsed_ms, e.message);
                self.op_error_popup = Some(e.message);
            }
        }
    }

    fn run_preset_digest_merge_ligate(&mut self) {
        let Some(engine) = &self.engine else {
            self.op_status = "No engine attached".to_string();
            return;
        };
        let template = self.seq_id.clone().unwrap_or_default();
        if template.is_empty() {
            self.op_status = "No active template sequence".to_string();
            return;
        }
        let enzymes = self
            .digest_enzymes_text
            .split(',')
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .map(|s| s.to_string())
            .collect::<Vec<_>>();
        if enzymes.is_empty() {
            self.op_status = "Provide at least one enzyme".to_string();
            return;
        }

        let started = Instant::now();
        let mut guard = engine.write().expect("Engine lock poisoned");
        let digest = match guard.apply(Operation::Digest {
            input: template,
            enzymes,
            output_prefix: Some(self.digest_prefix_text.clone()),
        }) {
            Ok(v) => v,
            Err(e) => {
                self.op_status = format!("Preset failed at Digest: {}", e.message);
                self.op_error_popup = Some(e.message);
                return;
            }
        };
        let merge = match guard.apply(Operation::MergeContainers {
            inputs: digest.created_seq_ids.clone(),
            output_prefix: Some("preset_merged".to_string()),
        }) {
            Ok(v) => v,
            Err(e) => {
                self.op_status = format!("Preset failed at Merge: {}", e.message);
                self.op_error_popup = Some(e.message);
                return;
            }
        };
        let ligate = match guard.apply(Operation::Ligation {
            inputs: merge.created_seq_ids.clone(),
            circularize_if_possible: false,
            output_id: None,
            protocol: LigationProtocol::Sticky,
            output_prefix: Some(self.ligation_output_prefix.clone()),
            unique: Some(false),
        }) {
            Ok(v) => v,
            Err(e) => {
                self.op_status = format!("Preset failed at Ligation: {}", e.message);
                self.op_error_popup = Some(e.message);
                return;
            }
        };
        let new_seq = ligate.created_seq_ids.first().cloned();
        drop(guard);

        if let Some(new_seq_id) = new_seq {
            if let Some(new_dna) = engine
                .read()
                .expect("Engine lock poisoned")
                .state()
                .sequences
                .get(&new_seq_id)
                .cloned()
            {
                self.seq_id = Some(new_seq_id);
                *self.dna.write().expect("DNA lock poisoned") = new_dna;
            }
        }
        self.op_status = format!(
            "preset ok in {} ms: digest {} -> merge {} -> ligate {}",
            started.elapsed().as_millis(),
            digest.created_seq_ids.len(),
            merge.created_seq_ids.len(),
            ligate.created_seq_ids.len()
        );
        self.op_error_popup = None;
    }

    fn parse_ids(text: &str) -> Vec<String> {
        text.split(',')
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .map(|s| s.to_string())
            .collect()
    }

    fn engine_ops_state_key(&self) -> String {
        let seq = self.seq_id.clone().unwrap_or_else(|| "_global".to_string());
        format!("gui.engine_ops.{seq}")
    }

    fn current_engine_ops_state(&self) -> EngineOpsUiState {
        EngineOpsUiState {
            show_engine_ops: self.show_engine_ops,
            digest_enzymes_text: self.digest_enzymes_text.clone(),
            digest_prefix_text: self.digest_prefix_text.clone(),
            merge_inputs_text: self.merge_inputs_text.clone(),
            ligation_inputs_text: self.ligation_inputs_text.clone(),
            ligation_output_prefix: self.ligation_output_prefix.clone(),
            ligation_protocol_sticky: self.ligation_protocol_sticky,
            ligation_unique: self.ligation_unique,
            mw_inputs_text: self.mw_inputs_text.clone(),
            mw_min_bp: self.mw_min_bp.clone(),
            mw_max_bp: self.mw_max_bp.clone(),
            mw_error: self.mw_error.clone(),
            mw_unique: self.mw_unique,
            pcr_forward: self.pcr_forward.clone(),
            pcr_reverse: self.pcr_reverse.clone(),
            pcr_unique: self.pcr_unique,
            pcr_adv_forward: self.pcr_adv_forward.clone(),
            pcr_adv_reverse: self.pcr_adv_reverse.clone(),
            pcr_adv_anneal_len: self.pcr_adv_anneal_len.clone(),
            pcr_adv_max_mismatch: self.pcr_adv_max_mismatch.clone(),
            pcr_adv_3prime_exact: self.pcr_adv_3prime_exact.clone(),
            pcr_mut_position: self.pcr_mut_position.clone(),
            pcr_mut_ref: self.pcr_mut_ref.clone(),
            pcr_mut_alt: self.pcr_mut_alt.clone(),
        }
    }

    fn apply_engine_ops_state(&mut self, s: EngineOpsUiState) {
        self.show_engine_ops = s.show_engine_ops;
        self.digest_enzymes_text = s.digest_enzymes_text;
        self.digest_prefix_text = s.digest_prefix_text;
        self.merge_inputs_text = s.merge_inputs_text;
        self.ligation_inputs_text = s.ligation_inputs_text;
        self.ligation_output_prefix = s.ligation_output_prefix;
        self.ligation_protocol_sticky = s.ligation_protocol_sticky;
        self.ligation_unique = s.ligation_unique;
        self.mw_inputs_text = s.mw_inputs_text;
        self.mw_min_bp = s.mw_min_bp;
        self.mw_max_bp = s.mw_max_bp;
        self.mw_error = s.mw_error;
        self.mw_unique = s.mw_unique;
        self.pcr_forward = s.pcr_forward;
        self.pcr_reverse = s.pcr_reverse;
        self.pcr_unique = s.pcr_unique;
        self.pcr_adv_forward = s.pcr_adv_forward;
        self.pcr_adv_reverse = s.pcr_adv_reverse;
        self.pcr_adv_anneal_len = s.pcr_adv_anneal_len;
        self.pcr_adv_max_mismatch = s.pcr_adv_max_mismatch;
        self.pcr_adv_3prime_exact = s.pcr_adv_3prime_exact;
        self.pcr_mut_position = s.pcr_mut_position;
        self.pcr_mut_ref = s.pcr_mut_ref;
        self.pcr_mut_alt = s.pcr_mut_alt;
    }

    fn save_engine_ops_state(&self) {
        let Some(engine) = &self.engine else {
            return;
        };
        let value = match serde_json::to_value(self.current_engine_ops_state()) {
            Ok(v) => v,
            Err(_) => return,
        };
        let mut guard = engine.write().expect("Engine lock poisoned");
        guard
            .state_mut()
            .metadata
            .insert(self.engine_ops_state_key(), value);
    }

    fn load_engine_ops_state(&mut self) {
        let Some(engine) = &self.engine else {
            return;
        };
        let key = self.engine_ops_state_key();
        let value = engine
            .read()
            .expect("Engine lock poisoned")
            .state()
            .metadata
            .get(&key)
            .cloned();
        if let Some(value) = value {
            if let Ok(state) = serde_json::from_value::<EngineOpsUiState>(value) {
                self.apply_engine_ops_state(state);
            }
        }
    }

    fn update_heartbeat(&mut self) {
        self.heartbeat_frame = self.heartbeat_frame.saturating_add(1);
        let now_ms = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_millis())
            .unwrap_or(0);
        self.heartbeat_delta_ms = self
            .heartbeat_last_ms
            .and_then(|prev| now_ms.checked_sub(prev));
        self.heartbeat_last_ms = Some(now_ms);
    }

    fn render_error_popup(&mut self, ctx: &egui::Context) {
        let Some(message) = self.op_error_popup.clone() else {
            return;
        };
        let mut open = true;
        egui::Window::new("Operation Failed")
            .open(&mut open)
            .collapsible(false)
            .resizable(false)
            .show(ctx, |ui| {
                ui.label(message);
                if ui.button("Close").clicked() {
                    self.op_error_popup = None;
                }
            });
        if !open {
            self.op_error_popup = None;
        }
    }

    fn render_pool_distribution(&self, ui: &mut egui::Ui) {
        let Some(engine) = &self.engine else {
            return;
        };
        let state = engine.read().expect("Engine lock poisoned");
        let mut lengths: Vec<usize> = self
            .last_created_seq_ids
            .iter()
            .filter_map(|id| state.state().sequences.get(id))
            .map(|dna| dna.len())
            .collect();
        if lengths.is_empty() {
            return;
        }
        lengths.sort_unstable();
        let min_bp = *lengths.first().unwrap_or(&0);
        let max_bp = *lengths.last().unwrap_or(&0);
        let median_bp = lengths[lengths.len() / 2];
        ui.monospace(format!(
            "Pool distribution: n={} min={}bp median={}bp max={}bp",
            lengths.len(),
            min_bp,
            median_bp,
            max_bp
        ));

        let span = max_bp.saturating_sub(min_bp);
        let bucket_size = if span <= 500 {
            50
        } else if span <= 2_000 {
            100
        } else {
            250
        };

        let mut buckets: HashMap<usize, usize> = HashMap::new();
        for bp in lengths {
            let idx = (bp - min_bp) / bucket_size;
            *buckets.entry(idx).or_insert(0) += 1;
        }
        let mut bucket_keys: Vec<usize> = buckets.keys().cloned().collect();
        bucket_keys.sort_unstable();
        for idx in bucket_keys {
            let count = buckets.get(&idx).cloned().unwrap_or(0);
            let lo = min_bp + idx * bucket_size;
            let hi = lo + bucket_size.saturating_sub(1);
            let bar = "#".repeat(count.min(40));
            ui.monospace(format!("{:>6}-{:>6} bp | {:>3} {}", lo, hi, count, bar));
        }
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
