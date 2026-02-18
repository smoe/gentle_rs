use crate::{
    dna_display::{DnaDisplay, TfbsDisplayCriteria},
    dna_sequence::DNAsequence,
    engine::{
        AnchorBoundary, AnchorDirection, AnchoredRegionAnchor, DisplayTarget, Engine, EngineError,
        ErrorCode, ExportFormat, GentleEngine, LigationProtocol, OpResult, Operation,
        OperationProgress, PcrPrimerSpec, RenderSvgMode, SnpMutationSpec, TfThresholdOverride,
        TfbsProgress, Workflow,
    },
    icons::*,
    render_dna::RenderDna,
    render_sequence::RenderSequence,
    tf_motifs,
};
use eframe::egui::{self, Frame, PointerState, Vec2};
use serde::{Deserialize, Serialize};
use std::{
    collections::HashMap,
    sync::{
        mpsc::{self, Receiver, TryRecvError},
        Arc, Mutex, RwLock,
    },
    time::{Duration, Instant},
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
    #[serde(default)]
    extract_from: String,
    #[serde(default)]
    extract_to: String,
    #[serde(default)]
    extract_output_id: String,
    #[serde(default)]
    parameter_name: String,
    #[serde(default)]
    parameter_value_json: String,
    #[serde(default)]
    anchored_mode_feature: bool,
    #[serde(default)]
    anchored_position: String,
    #[serde(default)]
    anchored_feature_kind: String,
    #[serde(default)]
    anchored_feature_label: String,
    #[serde(default)]
    anchored_feature_boundary_start: bool,
    #[serde(default)]
    anchored_feature_occurrence: String,
    #[serde(default)]
    anchored_direction_upstream: bool,
    #[serde(default)]
    anchored_target_len: String,
    #[serde(default)]
    anchored_tolerance: String,
    #[serde(default)]
    anchored_required_re_sites: String,
    #[serde(default)]
    anchored_required_tf_motifs: String,
    #[serde(default)]
    anchored_forward_primer: String,
    #[serde(default)]
    anchored_reverse_primer: String,
    #[serde(default)]
    anchored_output_prefix: String,
    #[serde(default)]
    anchored_unique: bool,
    #[serde(default)]
    anchored_max_candidates: String,
    #[serde(default)]
    tfbs_motifs: String,
    #[serde(default)]
    tfbs_use_all_motifs: bool,
    #[serde(default)]
    tfbs_catalog_filter: String,
    #[serde(default)]
    tfbs_min_llr_bits: String,
    #[serde(default)]
    tfbs_min_llr_quantile: String,
    #[serde(default)]
    tfbs_per_tf_min_llr_bits: String,
    #[serde(default)]
    tfbs_per_tf_min_llr_quantile: String,
    #[serde(default = "default_true")]
    tfbs_clear_existing: bool,
    #[serde(default = "default_true")]
    tfbs_display_use_llr_bits: bool,
    #[serde(default = "default_zero_f64")]
    tfbs_display_min_llr_bits: f64,
    #[serde(default = "default_true")]
    tfbs_display_use_llr_quantile: bool,
    #[serde(default = "default_tfbs_quantile")]
    tfbs_display_min_llr_quantile: f64,
    #[serde(default)]
    tfbs_display_use_true_log_odds_bits: bool,
    #[serde(default = "default_zero_f64")]
    tfbs_display_min_true_log_odds_bits: f64,
    #[serde(default)]
    tfbs_display_use_true_log_odds_quantile: bool,
    #[serde(default = "default_tfbs_quantile")]
    tfbs_display_min_true_log_odds_quantile: f64,
    #[serde(default)]
    container_digest_id: String,
    #[serde(default)]
    container_merge_ids: String,
    #[serde(default)]
    container_ligation_id: String,
    #[serde(default)]
    container_ligation_output_prefix: String,
    #[serde(default = "default_true")]
    container_ligation_protocol_sticky: bool,
    #[serde(default)]
    container_ligation_circularize: bool,
    #[serde(default)]
    container_ligation_unique: bool,
    #[serde(default)]
    container_mw_id: String,
    #[serde(default)]
    container_mw_min_bp: String,
    #[serde(default)]
    container_mw_max_bp: String,
    #[serde(default)]
    container_mw_error: String,
    #[serde(default)]
    container_mw_unique: bool,
    #[serde(default)]
    container_mw_output_prefix: String,
    #[serde(default)]
    workflow_run_id: String,
    #[serde(default)]
    workflow_ops_json: String,
    #[serde(default)]
    export_pool_inputs_text: String,
    #[serde(default)]
    export_pool_id: String,
    #[serde(default)]
    export_pool_human_id: String,
}

fn default_true() -> bool {
    true
}

fn default_zero_f64() -> f64 {
    0.0
}

fn default_tfbs_quantile() -> f64 {
    0.95
}

#[derive(Clone, Debug)]
enum TfbsTaskMessage {
    Progress(TfbsProgress),
    Done(Result<OpResult, EngineError>),
}

#[derive(Clone, Debug)]
struct TfbsTask {
    started: Instant,
    motif_count: usize,
    receiver: Arc<Mutex<Receiver<TfbsTaskMessage>>>,
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
    extract_from: String,
    extract_to: String,
    extract_output_id: String,
    parameter_name: String,
    parameter_value_json: String,
    anchored_mode_feature: bool,
    anchored_position: String,
    anchored_feature_kind: String,
    anchored_feature_label: String,
    anchored_feature_boundary_start: bool,
    anchored_feature_occurrence: String,
    anchored_direction_upstream: bool,
    anchored_target_len: String,
    anchored_tolerance: String,
    anchored_required_re_sites: String,
    anchored_required_tf_motifs: String,
    anchored_forward_primer: String,
    anchored_reverse_primer: String,
    anchored_output_prefix: String,
    anchored_unique: bool,
    anchored_max_candidates: String,
    tfbs_motifs: String,
    tfbs_use_all_motifs: bool,
    tfbs_catalog_filter: String,
    tfbs_min_llr_bits: String,
    tfbs_min_llr_quantile: String,
    tfbs_per_tf_min_llr_bits: String,
    tfbs_per_tf_min_llr_quantile: String,
    tfbs_clear_existing: bool,
    tfbs_task: Option<TfbsTask>,
    tfbs_progress: Option<TfbsProgress>,
    op_status: String,
    op_error_popup: Option<String>,
    last_created_seq_ids: Vec<String>,
    container_digest_id: String,
    container_merge_ids: String,
    container_ligation_id: String,
    container_ligation_output_prefix: String,
    container_ligation_protocol_sticky: bool,
    container_ligation_circularize: bool,
    container_ligation_unique: bool,
    container_mw_id: String,
    container_mw_min_bp: String,
    container_mw_max_bp: String,
    container_mw_error: String,
    container_mw_unique: bool,
    container_mw_output_prefix: String,
    workflow_run_id: String,
    workflow_ops_json: String,
    export_pool_inputs_text: String,
    export_pool_id: String,
    export_pool_human_id: String,
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
            extract_from: "0".to_string(),
            extract_to: "0".to_string(),
            extract_output_id: String::new(),
            parameter_name: "max_fragments_per_container".to_string(),
            parameter_value_json: "80000".to_string(),
            anchored_mode_feature: true,
            anchored_position: "0".to_string(),
            anchored_feature_kind: "CDS".to_string(),
            anchored_feature_label: String::new(),
            anchored_feature_boundary_start: true,
            anchored_feature_occurrence: "0".to_string(),
            anchored_direction_upstream: true,
            anchored_target_len: "500".to_string(),
            anchored_tolerance: "100".to_string(),
            anchored_required_re_sites: String::new(),
            anchored_required_tf_motifs: String::new(),
            anchored_forward_primer: String::new(),
            anchored_reverse_primer: String::new(),
            anchored_output_prefix: "anchored".to_string(),
            anchored_unique: false,
            anchored_max_candidates: "20".to_string(),
            tfbs_motifs: String::new(),
            tfbs_use_all_motifs: false,
            tfbs_catalog_filter: String::new(),
            tfbs_min_llr_bits: "0.0".to_string(),
            tfbs_min_llr_quantile: "0.95".to_string(),
            tfbs_per_tf_min_llr_bits: String::new(),
            tfbs_per_tf_min_llr_quantile: String::new(),
            tfbs_clear_existing: true,
            tfbs_task: None,
            tfbs_progress: None,
            op_status: String::new(),
            op_error_popup: None,
            last_created_seq_ids: vec![],
            container_digest_id: String::new(),
            container_merge_ids: String::new(),
            container_ligation_id: String::new(),
            container_ligation_output_prefix: "ligation_container".to_string(),
            container_ligation_protocol_sticky: true,
            container_ligation_circularize: false,
            container_ligation_unique: false,
            container_mw_id: String::new(),
            container_mw_min_bp: "100".to_string(),
            container_mw_max_bp: "1000".to_string(),
            container_mw_error: "0.10".to_string(),
            container_mw_unique: false,
            container_mw_output_prefix: "mw_container".to_string(),
            workflow_run_id: "gui-workflow".to_string(),
            workflow_ops_json: "[]".to_string(),
            export_pool_inputs_text: seq_id_for_defaults.clone().unwrap_or_default(),
            export_pool_id: String::new(),
            export_pool_human_id: String::new(),
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

    fn latest_container_for_active_seq(&self) -> Option<String> {
        let seq_id = self.seq_id.as_ref()?;
        let engine = self.engine.as_ref()?;
        engine.read().ok().and_then(|guard| {
            guard
                .state()
                .container_state
                .seq_to_latest_container
                .get(seq_id)
                .cloned()
        })
    }

    fn prefill_container_ids(&mut self) {
        let Some(container_id) = self.latest_container_for_active_seq() else {
            return;
        };
        if self.container_digest_id.trim().is_empty() {
            self.container_digest_id = container_id.clone();
        }
        if self.container_ligation_id.trim().is_empty() {
            self.container_ligation_id = container_id.clone();
        }
        if self.container_mw_id.trim().is_empty() {
            self.container_mw_id = container_id;
        }
    }

    pub fn render(&mut self, ctx: &egui::Context, ui: &mut egui::Ui) {
        self.prefill_container_ids();
        self.poll_tfbs_task(ctx);
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

            let tfbs_active = self
                .dna_display
                .read()
                .expect("DNA display lock poisoned")
                .show_tfbs();
            let response = ui
                .selectable_label(tfbs_active, "TFBS")
                .on_hover_text("Show or hide computed TFBS features");
            if response.clicked() {
                let visible = {
                    let display = self.dna_display.read().expect("DNA display lock poisoned");
                    !display.show_tfbs()
                };
                self.dna_display
                    .write()
                    .expect("DNA display lock poisoned")
                    .set_show_tfbs(visible);
                self.set_display_visibility(DisplayTarget::Tfbs, visible);
            }

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
            if ui
                .button("Export Seq")
                .on_hover_text("Export active sequence as GenBank or FASTA via engine SaveFile")
                .clicked()
            {
                self.export_active_sequence();
            }
            if ui
                .button("Export SVG")
                .on_hover_text("Export active sequence map SVG via engine RenderSequenceSvg")
                .clicked()
            {
                self.export_active_sequence_svg();
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
        });
        if !self.op_status.is_empty() {
            ui.add(egui::Label::new(egui::RichText::new(&self.op_status).monospace()).wrap());
        }

        if self.show_engine_ops {
            ui.separator();
            ui.collapsing("Strict Engine Operations", |ui| {
                let resize_id = format!(
                    "engine_ops_resize_{}",
                    self.seq_id.as_deref().unwrap_or("_global")
                );
                egui::Resize::default()
                    .id_salt(resize_id)
                    .default_height(420.0)
                    .min_height(180.0)
                    .max_height(ui.available_height().max(220.0))
                    .resizable(egui::Vec2b::new(false, true))
                    .show(ui, |ui| {
                        egui::ScrollArea::vertical()
                            .auto_shrink([false, false])
                            .show(ui, |ui| {
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
                ui.label("Region extraction and engine settings");
                ui.horizontal(|ui| {
                    ui.label("Extract from");
                    ui.text_edit_singleline(&mut self.extract_from);
                    ui.label("to");
                    ui.text_edit_singleline(&mut self.extract_to);
                    ui.label("output_id");
                    ui.text_edit_singleline(&mut self.extract_output_id);
                    if ui.button("Extract Region").clicked() {
                        let template = self.seq_id.clone().unwrap_or_default();
                        if template.is_empty() {
                            self.op_status = "No active template sequence".to_string();
                        } else if let (Ok(from), Ok(to)) = (
                            self.extract_from.parse::<usize>(),
                            self.extract_to.parse::<usize>(),
                        ) {
                            self.apply_operation_with_feedback(Operation::ExtractRegion {
                                input: template,
                                from,
                                to,
                                output_id: if self.extract_output_id.trim().is_empty() {
                                    None
                                } else {
                                    Some(self.extract_output_id.trim().to_string())
                                },
                            });
                        } else {
                            self.op_status = "Invalid ExtractRegion bounds".to_string();
                        }
                    }
                });
                ui.horizontal(|ui| {
                    if ui.button("Recompute Features").clicked() {
                        let seq_id = self.seq_id.clone().unwrap_or_default();
                        if seq_id.is_empty() {
                            self.op_status = "No active template sequence".to_string();
                        } else {
                            self.apply_operation_with_feedback(Operation::RecomputeFeatures {
                                seq_id,
                            });
                        }
                    }
                    ui.label("SetParameter name");
                    ui.text_edit_singleline(&mut self.parameter_name);
                    ui.label("json");
                    ui.text_edit_singleline(&mut self.parameter_value_json);
                    if ui.button("Set Parameter").clicked() {
                        if self.parameter_name.trim().is_empty() {
                            self.op_status = "SetParameter name cannot be empty".to_string();
                        } else {
                            match serde_json::from_str::<serde_json::Value>(
                                self.parameter_value_json.trim(),
                            ) {
                                Ok(value) => {
                                    self.apply_operation_with_feedback(Operation::SetParameter {
                                        name: self.parameter_name.trim().to_string(),
                                        value,
                                    });
                                }
                                Err(e) => {
                                    self.op_status = format!("Invalid SetParameter JSON value: {e}");
                                }
                            }
                        }
                    }
                });

                ui.separator();
                ui.label("Container-first operations");
                ui.horizontal(|ui| {
                    ui.label("Digest container_id");
                    ui.text_edit_singleline(&mut self.container_digest_id);
                    ui.label("enzymes");
                    ui.text_edit_singleline(&mut self.digest_enzymes_text);
                    ui.label("prefix");
                    ui.text_edit_singleline(&mut self.digest_prefix_text);
                    if ui.button("Digest Container").clicked() {
                        let enzymes = self
                            .digest_enzymes_text
                            .split(',')
                            .map(|s| s.trim())
                            .filter(|s| !s.is_empty())
                            .map(|s| s.to_string())
                            .collect::<Vec<_>>();
                        if self.container_digest_id.trim().is_empty() {
                            self.op_status = "Provide container_id for DigestContainer".to_string();
                        } else if enzymes.is_empty() {
                            self.op_status = "Provide at least one enzyme".to_string();
                        } else {
                            self.apply_operation_with_feedback(Operation::DigestContainer {
                                container_id: self.container_digest_id.trim().to_string(),
                                enzymes,
                                output_prefix: Some(self.digest_prefix_text.clone()),
                            });
                        }
                    }
                });
                ui.horizontal(|ui| {
                    ui.label("Merge container_ids");
                    ui.text_edit_singleline(&mut self.container_merge_ids);
                    if ui.button("Merge ContainersById").clicked() {
                        let container_ids = Self::parse_ids(&self.container_merge_ids);
                        if container_ids.is_empty() {
                            self.op_status =
                                "Provide at least one container_id for merge".to_string();
                        } else {
                            self.apply_operation_with_feedback(Operation::MergeContainersById {
                                container_ids,
                                output_prefix: Some("merged_container".to_string()),
                            });
                        }
                    }
                });
                ui.horizontal(|ui| {
                    ui.label("Ligation container_id");
                    ui.text_edit_singleline(&mut self.container_ligation_id);
                    ui.checkbox(
                        &mut self.container_ligation_protocol_sticky,
                        "Sticky (off = Blunt)",
                    );
                    ui.checkbox(&mut self.container_ligation_circularize, "Circularize");
                    ui.checkbox(&mut self.container_ligation_unique, "Unique");
                });
                ui.horizontal(|ui| {
                    ui.label("prefix");
                    ui.text_edit_singleline(&mut self.container_ligation_output_prefix);
                    if ui.button("Ligate Container").clicked() {
                        if self.container_ligation_id.trim().is_empty() {
                            self.op_status =
                                "Provide container_id for LigationContainer".to_string();
                        } else {
                            self.apply_operation_with_feedback(Operation::LigationContainer {
                                container_id: self.container_ligation_id.trim().to_string(),
                                circularize_if_possible: self.container_ligation_circularize,
                                output_id: None,
                                protocol: if self.container_ligation_protocol_sticky {
                                    LigationProtocol::Sticky
                                } else {
                                    LigationProtocol::Blunt
                                },
                                output_prefix: if self.container_ligation_output_prefix.is_empty() {
                                    None
                                } else {
                                    Some(self.container_ligation_output_prefix.clone())
                                },
                                unique: Some(self.container_ligation_unique),
                            });
                        }
                    }
                });
                ui.horizontal(|ui| {
                    ui.label("MW container_id");
                    ui.text_edit_singleline(&mut self.container_mw_id);
                    ui.label("min");
                    ui.text_edit_singleline(&mut self.container_mw_min_bp);
                    ui.label("max");
                    ui.text_edit_singleline(&mut self.container_mw_max_bp);
                    ui.label("error");
                    ui.text_edit_singleline(&mut self.container_mw_error);
                    ui.checkbox(&mut self.container_mw_unique, "Unique");
                });
                ui.horizontal(|ui| {
                    ui.label("prefix");
                    ui.text_edit_singleline(&mut self.container_mw_output_prefix);
                    if ui.button("Filter Container MW").clicked() {
                        if self.container_mw_id.trim().is_empty() {
                            self.op_status =
                                "Provide container_id for FilterContainerByMolecularWeight"
                                    .to_string();
                        } else if let (Ok(min_bp), Ok(max_bp), Ok(error)) = (
                            self.container_mw_min_bp.parse::<usize>(),
                            self.container_mw_max_bp.parse::<usize>(),
                            self.container_mw_error.parse::<f64>(),
                        ) {
                            self.apply_operation_with_feedback(
                                Operation::FilterContainerByMolecularWeight {
                                    container_id: self.container_mw_id.trim().to_string(),
                                    min_bp,
                                    max_bp,
                                    error,
                                    unique: self.container_mw_unique,
                                    output_prefix: if self.container_mw_output_prefix.is_empty() {
                                        None
                                    } else {
                                        Some(self.container_mw_output_prefix.clone())
                                    },
                                },
                            );
                        } else {
                            self.op_status =
                                "Invalid container MW parameters (min/max/error)".to_string();
                        }
                    }
                });

                ui.separator();
                ui.label("Workflow runner (JSON array of operations)");
                ui.horizontal(|ui| {
                    ui.label("run_id");
                    ui.text_edit_singleline(&mut self.workflow_run_id);
                    if ui.button("Run Workflow").clicked() {
                        let ops = serde_json::from_str::<Vec<Operation>>(
                            self.workflow_ops_json.trim(),
                        );
                        match ops {
                            Ok(ops) => {
                                if ops.is_empty() {
                                    self.op_status = "Workflow has no operations".to_string();
                                } else {
                                    self.run_workflow_with_feedback(Workflow {
                                        run_id: if self.workflow_run_id.trim().is_empty() {
                                            "gui-workflow".to_string()
                                        } else {
                                            self.workflow_run_id.trim().to_string()
                                        },
                                        ops,
                                    });
                                }
                            }
                            Err(e) => {
                                self.op_status = format!("Invalid workflow JSON: {e}");
                            }
                        }
                    }
                });
                ui.add(
                    egui::TextEdit::multiline(&mut self.workflow_ops_json)
                        .desired_rows(6)
                        .hint_text(
                            r#"[{"Digest":{"input":"seq_1","enzymes":["EcoRI"],"output_prefix":"d"}}]"#,
                        ),
                );

                ui.separator();
                ui.label("Pool export");
                ui.horizontal(|ui| {
                    ui.label("inputs");
                    ui.text_edit_singleline(&mut self.export_pool_inputs_text);
                });
                ui.horizontal(|ui| {
                    ui.label("pool_id");
                    ui.text_edit_singleline(&mut self.export_pool_id);
                    ui.label("human_id");
                    ui.text_edit_singleline(&mut self.export_pool_human_id);
                    if ui.button("Export Pool").clicked() {
                        self.export_pool_to_file();
                    }
                });

                ui.separator();
                ui.label("Anchored region extraction (promoter-like)");
                if ui
                    .button("Preset: Promoter (CDS start, upstream 500±100)")
                    .on_hover_text(
                        "Fill anchored extraction with promoter defaults: CDS start anchor, upstream, 500 bp target length with 100 bp tolerance",
                    )
                    .clicked()
                {
                    self.anchored_mode_feature = true;
                    self.anchored_feature_kind = "CDS".to_string();
                    self.anchored_feature_label.clear();
                    self.anchored_feature_boundary_start = true;
                    self.anchored_feature_occurrence = "0".to_string();
                    self.anchored_direction_upstream = true;
                    self.anchored_target_len = "500".to_string();
                    self.anchored_tolerance = "100".to_string();
                    self.anchored_required_re_sites.clear();
                    self.anchored_required_tf_motifs.clear();
                    self.anchored_forward_primer.clear();
                    self.anchored_reverse_primer.clear();
                    self.anchored_output_prefix = "promoter".to_string();
                    self.anchored_unique = false;
                    self.anchored_max_candidates = "20".to_string();
                }
                ui.horizontal(|ui| {
                    ui.checkbox(
                        &mut self.anchored_mode_feature,
                        "Feature boundary anchor (off = absolute position)",
                    );
                    ui.checkbox(
                        &mut self.anchored_direction_upstream,
                        "Upstream (off = Downstream)",
                    );
                    ui.checkbox(&mut self.anchored_unique, "Unique");
                });
                if self.anchored_mode_feature {
                    ui.horizontal(|ui| {
                        ui.label("feature kind");
                        ui.text_edit_singleline(&mut self.anchored_feature_kind);
                        ui.label("feature label");
                        ui.text_edit_singleline(&mut self.anchored_feature_label);
                    });
                    ui.horizontal(|ui| {
                        ui.checkbox(&mut self.anchored_feature_boundary_start, "Boundary Start");
                        ui.label("occurrence");
                        ui.text_edit_singleline(&mut self.anchored_feature_occurrence);
                    });
                } else {
                    ui.horizontal(|ui| {
                        ui.label("anchor pos (0-based)");
                        ui.text_edit_singleline(&mut self.anchored_position);
                    });
                }
                ui.horizontal(|ui| {
                    ui.label("target len");
                    ui.text_edit_singleline(&mut self.anchored_target_len);
                    ui.label("tolerance");
                    ui.text_edit_singleline(&mut self.anchored_tolerance);
                    ui.label("max candidates");
                    ui.text_edit_singleline(&mut self.anchored_max_candidates);
                });
                ui.horizontal(|ui| {
                    ui.label("required RE sites");
                    ui.text_edit_singleline(&mut self.anchored_required_re_sites);
                });
                ui.horizontal(|ui| {
                    ui.label("required TF motifs");
                    ui.text_edit_singleline(&mut self.anchored_required_tf_motifs);
                });
                ui.horizontal(|ui| {
                    ui.label("fwd primer");
                    ui.text_edit_singleline(&mut self.anchored_forward_primer);
                    ui.label("rev primer");
                    ui.text_edit_singleline(&mut self.anchored_reverse_primer);
                });
                ui.horizontal(|ui| {
                    ui.label("prefix");
                    ui.text_edit_singleline(&mut self.anchored_output_prefix);
                    if ui.button("Extract Anchored").clicked() {
                        let template = self.seq_id.clone().unwrap_or_default();
                        if template.is_empty() {
                            self.op_status = "No active template sequence".to_string();
                        } else if let (Ok(target_length_bp), Ok(length_tolerance_bp)) = (
                            self.anchored_target_len.parse::<usize>(),
                            self.anchored_tolerance.parse::<usize>(),
                        ) {
                            let anchor_result: Result<AnchoredRegionAnchor, &'static str> = if self
                                .anchored_mode_feature
                            {
                                let occurrence =
                                    if self.anchored_feature_occurrence.trim().is_empty() {
                                        Ok(None)
                                    } else {
                                        self.anchored_feature_occurrence
                                            .parse::<usize>()
                                            .map(Some)
                                            .map_err(|_| "Invalid anchored occurrence")
                                    };
                                occurrence.map(|occurrence| AnchoredRegionAnchor::FeatureBoundary {
                                    feature_kind: if self.anchored_feature_kind.trim().is_empty() {
                                        None
                                    } else {
                                        Some(self.anchored_feature_kind.clone())
                                    },
                                    feature_label: if self.anchored_feature_label.trim().is_empty()
                                    {
                                        None
                                    } else {
                                        Some(self.anchored_feature_label.clone())
                                    },
                                    boundary: if self.anchored_feature_boundary_start {
                                        AnchorBoundary::Start
                                    } else {
                                        AnchorBoundary::End
                                    },
                                    occurrence,
                                })
                            } else {
                                self.anchored_position
                                    .parse::<usize>()
                                    .map(|zero_based| AnchoredRegionAnchor::Position { zero_based })
                                    .map_err(|_| "Invalid anchored position (0-based)")
                            };

                            let max_candidates_result: Result<Option<usize>, &'static str> =
                                if self.anchored_max_candidates.trim().is_empty() {
                                    Ok(None)
                                } else {
                                    self.anchored_max_candidates
                                        .parse::<usize>()
                                        .map(Some)
                                        .map_err(|_| "Invalid anchored max candidates")
                                };

                            match (anchor_result, max_candidates_result) {
                                (Ok(anchor), Ok(max_candidates)) => {
                                    self.apply_operation_with_feedback(
                                        Operation::ExtractAnchoredRegion {
                                            input: template,
                                            anchor,
                                            direction: if self.anchored_direction_upstream {
                                                AnchorDirection::Upstream
                                            } else {
                                                AnchorDirection::Downstream
                                            },
                                            target_length_bp,
                                            length_tolerance_bp,
                                            required_re_sites: Self::parse_ids(
                                                &self.anchored_required_re_sites,
                                            ),
                                            required_tf_motifs: Self::parse_ids(
                                                &self.anchored_required_tf_motifs,
                                            ),
                                            forward_primer: if self
                                                .anchored_forward_primer
                                                .trim()
                                                .is_empty()
                                            {
                                                None
                                            } else {
                                                Some(self.anchored_forward_primer.clone())
                                            },
                                            reverse_primer: if self
                                                .anchored_reverse_primer
                                                .trim()
                                                .is_empty()
                                            {
                                                None
                                            } else {
                                                Some(self.anchored_reverse_primer.clone())
                                            },
                                            output_prefix: if self
                                                .anchored_output_prefix
                                                .trim()
                                                .is_empty()
                                            {
                                                None
                                            } else {
                                                Some(self.anchored_output_prefix.clone())
                                            },
                                            unique: Some(self.anchored_unique),
                                            max_candidates,
                                        },
                                    );
                                }
                                (Err(e), _) | (_, Err(e)) => {
                                    self.op_status = e.to_string();
                                }
                            }
                        } else {
                            self.op_status =
                                "Invalid anchored parameters (target len / tolerance)".to_string();
                        }
                    }
                });

                ui.separator();
                ui.label("TFBS annotation (log-likelihood ratio)");
                ui.horizontal(|ui| {
                    ui.checkbox(&mut self.tfbs_use_all_motifs, "All known JASPAR motifs");
                    if self.tfbs_use_all_motifs {
                        let motif_count = tf_motifs::all_motif_ids().len();
                        ui.label(format!("({motif_count} motifs)"));
                    }
                });
                if !self.tfbs_use_all_motifs {
                    ui.horizontal(|ui| {
                        ui.label("Selected motifs (IDs/names/IUPAC)");
                        ui.text_edit_singleline(&mut self.tfbs_motifs);
                    });
                    ui.horizontal(|ui| {
                        ui.label("JASPAR filter");
                        ui.text_edit_singleline(&mut self.tfbs_catalog_filter);
                        if ui.small_button("Clear list").clicked() {
                            self.tfbs_motifs.clear();
                        }
                    });

                    let filter = self.tfbs_catalog_filter.trim().to_ascii_uppercase();
                    let motifs = tf_motifs::list_motif_summaries();
                    let mut shown = 0usize;
                    egui::ScrollArea::vertical()
                        .max_height(110.0)
                        .show(ui, |ui| {
                            for motif in motifs.iter() {
                                let id_match = motif.id.to_ascii_uppercase().contains(&filter);
                                let name_match = motif
                                    .name
                                    .as_deref()
                                    .unwrap_or("")
                                    .to_ascii_uppercase()
                                    .contains(&filter);
                                if !filter.is_empty() && !id_match && !name_match {
                                    continue;
                                }
                                shown += 1;
                                if shown > 40 {
                                    break;
                                }
                                ui.horizontal(|ui| {
                                    if ui.small_button("+").clicked() {
                                        self.add_tfbs_motif_selection(&motif.id);
                                    }
                                    if let Some(name) = &motif.name {
                                        ui.label(format!("{} ({})", motif.id, name));
                                    } else {
                                        ui.label(&motif.id);
                                    }
                                });
                            }
                        });
                    ui.small(format!("showing up to 40 matches, currently {shown}"));
                } else {
                    ui.small("Selection list is ignored while 'All known JASPAR motifs' is enabled.");
                }
                ui.horizontal(|ui| {
                    ui.label("min llr_bits");
                    ui.text_edit_singleline(&mut self.tfbs_min_llr_bits);
                    ui.label("min llr_quantile");
                    ui.text_edit_singleline(&mut self.tfbs_min_llr_quantile);
                    ui.checkbox(&mut self.tfbs_clear_existing, "Clear previous TFBS");
                });
                ui.horizontal(|ui| {
                    ui.label("per-TF min llr_bits (TF=VALUE,...)");
                    ui.text_edit_singleline(&mut self.tfbs_per_tf_min_llr_bits);
                });
                ui.horizontal(|ui| {
                    ui.label("per-TF min llr_quantile (TF=VALUE,...)");
                    ui.text_edit_singleline(&mut self.tfbs_per_tf_min_llr_quantile);
                });
                ui.separator();
                ui.label("TFBS display filter (checkbox = criterion enabled)");
                let mut tfbs_display = self
                    .dna_display
                    .read()
                    .expect("DNA display lock poisoned")
                    .tfbs_display_criteria();
                let before = tfbs_display;
                ui.horizontal(|ui| {
                    ui.checkbox(&mut tfbs_display.use_llr_bits, "llr_bits");
                    ui.add_enabled(
                        tfbs_display.use_llr_bits,
                        egui::DragValue::new(&mut tfbs_display.min_llr_bits).speed(0.1),
                    );
                });
                ui.horizontal(|ui| {
                    ui.checkbox(&mut tfbs_display.use_llr_quantile, "llr_quantile");
                    ui.add_enabled(
                        tfbs_display.use_llr_quantile,
                        egui::DragValue::new(&mut tfbs_display.min_llr_quantile)
                            .range(0.0..=1.0)
                            .speed(0.01),
                    );
                });
                ui.horizontal(|ui| {
                    ui.checkbox(
                        &mut tfbs_display.use_true_log_odds_bits,
                        "true_log_odds_bits",
                    );
                    ui.add_enabled(
                        tfbs_display.use_true_log_odds_bits,
                        egui::DragValue::new(&mut tfbs_display.min_true_log_odds_bits).speed(0.1),
                    );
                });
                ui.horizontal(|ui| {
                    ui.checkbox(
                        &mut tfbs_display.use_true_log_odds_quantile,
                        "true_log_odds_quantile",
                    );
                    ui.add_enabled(
                        tfbs_display.use_true_log_odds_quantile,
                        egui::DragValue::new(&mut tfbs_display.min_true_log_odds_quantile)
                            .range(0.0..=1.0)
                            .speed(0.01),
                    );
                });
                if !tfbs_display.use_llr_bits
                    && !tfbs_display.use_llr_quantile
                    && !tfbs_display.use_true_log_odds_bits
                    && !tfbs_display.use_true_log_odds_quantile
                {
                    ui.small("No display criterion enabled: all TFBS are shown.");
                }
                if tfbs_display != before {
                    self.dna_display
                        .write()
                        .expect("DNA display lock poisoned")
                        .set_tfbs_display_criteria(tfbs_display);
                    self.sync_tfbs_display_criteria_to_engine(tfbs_display);
                    self.save_engine_ops_state();
                }
                if let Some(task) = &self.tfbs_task {
                    ui.horizontal(|ui| {
                        ui.add(egui::Spinner::new());
                        if let Some(progress) = &self.tfbs_progress {
                            ui.label(format!(
                                "Annotating TFBS... motif {}/{} ({}) {:.1}% | total {:.1}% | {:.1}s elapsed",
                                progress.motif_index,
                                progress.motif_count,
                                progress.motif_id,
                                progress.motif_percent,
                                progress.total_percent,
                                task.started.elapsed().as_secs_f32()
                            ));
                        } else {
                            ui.label(format!(
                                "Annotating TFBS... {} motif(s), {:.1}s elapsed",
                                task.motif_count,
                                task.started.elapsed().as_secs_f32()
                            ));
                        }
                    });
                    if let Some(progress) = &self.tfbs_progress {
                        ui.add(
                            egui::ProgressBar::new((progress.total_percent / 100.0) as f32)
                                .show_percentage()
                                .text(format!(
                                    "Total TFs addressed: {}/{}",
                                    progress.motif_index, progress.motif_count
                                )),
                        );
                        ui.add(
                            egui::ProgressBar::new((progress.motif_percent / 100.0) as f32)
                                .show_percentage()
                                .text(format!("Current motif: {}", progress.motif_id)),
                        );
                    }
                }
                if ui
                    .add_enabled(self.tfbs_task.is_none(), egui::Button::new("Annotate TFBS"))
                    .clicked()
                {
                    self.op_status = "Validating TFBS inputs...".to_string();
                    let seq_id = self.seq_id.clone().unwrap_or_default();
                    if seq_id.is_empty() {
                        self.op_status = "No active sequence".to_string();
                    } else {
                        let motifs = self.collect_tfbs_motifs();
                        if motifs.is_empty() {
                            self.op_status = "Provide at least one TF motif, or enable all JASPAR motifs".to_string();
                        } else {
                            let min_llr_bits = if self.tfbs_min_llr_bits.trim().is_empty() {
                                Ok(None)
                            } else {
                                self.tfbs_min_llr_bits
                                    .trim()
                                    .parse::<f64>()
                                    .map(Some)
                                    .map_err(|_| "Invalid min llr_bits value".to_string())
                            };
                            let min_llr_quantile = if self.tfbs_min_llr_quantile.trim().is_empty() {
                                Ok(None)
                            } else {
                                self.tfbs_min_llr_quantile
                                    .trim()
                                    .parse::<f64>()
                                    .map(Some)
                                    .map_err(|_| "Invalid min llr_quantile value".to_string())
                            };
                            let per_bits =
                                Self::parse_tf_threshold_map(&self.tfbs_per_tf_min_llr_bits);
                            let per_quantiles =
                                Self::parse_tf_threshold_map(&self.tfbs_per_tf_min_llr_quantile);

                            match (min_llr_bits, min_llr_quantile, per_bits, per_quantiles) {
                                (
                                    Ok(min_llr_bits),
                                    Ok(min_llr_quantile),
                                    Ok(per_bits),
                                    Ok(per_quantiles),
                                ) => {
                                    let mut merged: HashMap<String, TfThresholdOverride> =
                                        HashMap::new();
                                    for (tf, v) in per_bits {
                                        merged.insert(
                                            tf.clone(),
                                            TfThresholdOverride {
                                                tf,
                                                min_llr_bits: Some(v),
                                                min_llr_quantile: None,
                                            },
                                        );
                                    }
                                    for (tf, v) in per_quantiles {
                                        let entry = merged
                                            .entry(tf.clone())
                                            .or_insert(TfThresholdOverride {
                                                tf: tf.clone(),
                                                min_llr_bits: None,
                                                min_llr_quantile: None,
                                            });
                                        entry.min_llr_quantile = Some(v);
                                    }

                                    let op = Operation::AnnotateTfbs {
                                        seq_id,
                                        motifs,
                                        min_llr_bits,
                                        min_llr_quantile,
                                        per_tf_thresholds: merged.into_values().collect(),
                                        clear_existing: Some(self.tfbs_clear_existing),
                                        max_hits: None,
                                    };
                                    self.start_tfbs_annotation(op);
                                }
                                (Err(e), _, _, _)
                                | (_, Err(e), _, _)
                                | (_, _, Err(e), _)
                                | (_, _, _, Err(e)) => {
                                    self.op_status = e;
                                }
                            }
                        }
                    }
                }

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

                            });
                    });
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

    fn sync_tfbs_display_criteria_to_engine(&self, criteria: TfbsDisplayCriteria) {
        let Some(engine) = &self.engine else {
            return;
        };
        let mut guard = engine.write().expect("Engine lock poisoned");
        let display = &mut guard.state_mut().display;
        display.tfbs_display_use_llr_bits = criteria.use_llr_bits;
        display.tfbs_display_min_llr_bits = criteria.min_llr_bits;
        display.tfbs_display_use_llr_quantile = criteria.use_llr_quantile;
        display.tfbs_display_min_llr_quantile = criteria.min_llr_quantile;
        display.tfbs_display_use_true_log_odds_bits = criteria.use_true_log_odds_bits;
        display.tfbs_display_min_true_log_odds_bits = criteria.min_true_log_odds_bits;
        display.tfbs_display_use_true_log_odds_quantile = criteria.use_true_log_odds_quantile;
        display.tfbs_display_min_true_log_odds_quantile = criteria.min_true_log_odds_quantile;
    }

    fn sync_from_engine_display(&mut self) {
        let Some(engine) = &self.engine else {
            return;
        };
        let Ok(guard) = engine.try_read() else {
            return;
        };
        let settings = guard.state().display.clone();
        drop(guard);
        self.show_sequence = settings.show_sequence_panel;
        self.show_map = settings.show_map_panel;
        let mut display = self.dna_display.write().expect("DNA display lock poisoned");
        display.set_show_features(settings.show_features);
        display.set_show_tfbs(settings.show_tfbs);
        display.set_tfbs_display_criteria(TfbsDisplayCriteria {
            use_llr_bits: settings.tfbs_display_use_llr_bits,
            min_llr_bits: settings.tfbs_display_min_llr_bits,
            use_llr_quantile: settings.tfbs_display_use_llr_quantile,
            min_llr_quantile: settings.tfbs_display_min_llr_quantile,
            use_true_log_odds_bits: settings.tfbs_display_use_true_log_odds_bits,
            min_true_log_odds_bits: settings.tfbs_display_min_true_log_odds_bits,
            use_true_log_odds_quantile: settings.tfbs_display_use_true_log_odds_quantile,
            min_true_log_odds_quantile: settings.tfbs_display_min_true_log_odds_quantile,
        });
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
        let Some(engine) = self.engine.clone() else {
            self.op_status = "No engine attached".to_string();
            return;
        };
        let started = Instant::now();
        let apply_result = {
            let mut guard = engine.write().expect("Engine lock poisoned");
            guard.apply(op)
        };
        match apply_result {
            Ok(result) => self.handle_operation_success(result, started),
            Err(e) => self.handle_operation_error(e, started),
        }
    }

    fn handle_operation_success(&mut self, result: OpResult, started: Instant) {
        self.last_created_seq_ids = result.created_seq_ids.clone();
        if !self.last_created_seq_ids.is_empty() {
            self.export_pool_inputs_text = self.last_created_seq_ids.join(", ");
        }
        if let Some(engine) = &self.engine {
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
            } else if let Some(current_seq_id) = self.seq_id.clone() {
                if result
                    .changed_seq_ids
                    .iter()
                    .any(|id| id == &current_seq_id)
                {
                    if let Some(updated_dna) = engine
                        .read()
                        .expect("Engine lock poisoned")
                        .state()
                        .sequences
                        .get(&current_seq_id)
                        .cloned()
                    {
                        *self.dna.write().expect("DNA lock poisoned") = updated_dna;
                    }
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

    fn handle_operation_error(&mut self, e: EngineError, started: Instant) {
        self.last_created_seq_ids.clear();
        let elapsed_ms = started.elapsed().as_millis();
        self.op_status = format!("error after {} ms: {}", elapsed_ms, e.message);
        self.op_error_popup = Some(e.message);
    }

    fn start_tfbs_annotation(&mut self, op: Operation) {
        let Some(engine) = self.engine.clone() else {
            self.op_status = "No engine attached".to_string();
            return;
        };
        if self.tfbs_task.is_some() {
            self.op_status = "TFBS annotation is already running".to_string();
            return;
        }
        let motif_count = match &op {
            Operation::AnnotateTfbs { motifs, .. } => motifs.len(),
            _ => 0,
        };
        let (tx, rx) = mpsc::channel::<TfbsTaskMessage>();
        let started = Instant::now();
        self.tfbs_progress = None;
        self.op_status = format!("TFBS annotation started for {motif_count} motif(s)");
        self.tfbs_task = Some(TfbsTask {
            started,
            motif_count,
            receiver: Arc::new(Mutex::new(rx)),
        });
        std::thread::spawn(move || {
            let tx_progress = tx.clone();
            let outcome = match engine.write() {
                Ok(mut guard) => {
                    let mut last_motif_index: Option<usize> = None;
                    let mut last_total_tenths: Option<i64> = None;
                    guard.apply_with_progress(op, move |OperationProgress::Tfbs(p)| {
                        let total_tenths = (p.total_percent * 10.0).floor() as i64;
                        let motif_changed = last_motif_index != Some(p.motif_index);
                        let progressed = last_total_tenths
                            .map(|prev| total_tenths > prev)
                            .unwrap_or(true);
                        let done = (p.total_percent - 100.0).abs() < f64::EPSILON;
                        if motif_changed || progressed || done {
                            last_motif_index = Some(p.motif_index);
                            last_total_tenths = Some(total_tenths);
                            let _ = tx_progress.send(TfbsTaskMessage::Progress(p));
                        }
                    })
                }
                Err(_) => Err(EngineError {
                    code: ErrorCode::Internal,
                    message: "Engine lock poisoned while running TFBS annotation".to_string(),
                }),
            };
            let _ = tx.send(TfbsTaskMessage::Done(outcome));
        });
    }

    fn poll_tfbs_task(&mut self, ctx: &egui::Context) {
        if self.tfbs_task.is_none() {
            return;
        }

        let mut done: Option<Result<OpResult, EngineError>> = None;
        if let Some(task) = &self.tfbs_task {
            ctx.request_repaint_after(Duration::from_millis(100));
            match task.receiver.lock() {
                Ok(rx) => {
                    const MAX_PROGRESS_MESSAGES_PER_TICK: usize = 256;
                    let mut processed_progress = 0usize;
                    loop {
                        match rx.try_recv() {
                            Ok(TfbsTaskMessage::Progress(progress)) => {
                                self.tfbs_progress = Some(progress);
                                processed_progress = processed_progress.saturating_add(1);
                                if processed_progress >= MAX_PROGRESS_MESSAGES_PER_TICK {
                                    break;
                                }
                            }
                            Ok(TfbsTaskMessage::Done(res)) => {
                                done = Some(res);
                                break;
                            }
                            Err(TryRecvError::Disconnected) => {
                                done = Some(Err(EngineError {
                                    code: ErrorCode::Internal,
                                    message: "TFBS worker disconnected unexpectedly".to_string(),
                                }));
                                break;
                            }
                            Err(TryRecvError::Empty) => {
                                break;
                            }
                        }
                    }
                    if done.is_none() {
                        if let Some(progress) = &self.tfbs_progress {
                            self.op_status = format!(
                                "TFBS annotation running: motif {}/{} ({}) {:.1}%, total {:.1}%, elapsed {:.1}s",
                                progress.motif_index,
                                progress.motif_count,
                                progress.motif_id,
                                progress.motif_percent,
                                progress.total_percent,
                                task.started.elapsed().as_secs_f32()
                            );
                        } else {
                            self.op_status = format!(
                                "TFBS annotation running... {} motif(s), {:.1}s elapsed",
                                task.motif_count,
                                task.started.elapsed().as_secs_f32()
                            );
                        }
                    }
                }
                Err(_) => {
                    done = Some(Err(EngineError {
                        code: ErrorCode::Internal,
                        message: "TFBS progress channel lock poisoned".to_string(),
                    }));
                }
            }
        }

        if let Some(done) = done {
            let started = self
                .tfbs_task
                .as_ref()
                .map(|t| t.started)
                .unwrap_or_else(Instant::now);
            self.tfbs_task = None;
            self.tfbs_progress = None;
            match done {
                Ok(result) => self.handle_operation_success(result, started),
                Err(e) => self.handle_operation_error(e, started),
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

    fn run_workflow_with_feedback(&mut self, workflow: Workflow) {
        let Some(engine) = self.engine.clone() else {
            self.op_status = "No engine attached".to_string();
            return;
        };
        let started = Instant::now();
        let workflow_run_id = workflow.run_id.clone();
        let apply_result = {
            let mut guard = engine.write().expect("Engine lock poisoned");
            guard.apply_workflow(workflow)
        };
        match apply_result {
            Ok(results) => {
                if results.is_empty() {
                    self.last_created_seq_ids.clear();
                    self.op_status = format!(
                        "workflow '{}' completed in {} ms with 0 operations",
                        workflow_run_id,
                        started.elapsed().as_millis()
                    );
                    self.op_error_popup = None;
                    return;
                }
                let op_count = results.len();
                let total_created: usize = results.iter().map(|r| r.created_seq_ids.len()).sum();
                let last_result = results.last().cloned().expect("non-empty workflow results");
                self.handle_operation_success(last_result, started);
                self.op_status = format!(
                    "workflow '{}' ok: {} op(s), {} created sequence(s)\n{}",
                    workflow_run_id, op_count, total_created, self.op_status
                );
            }
            Err(e) => self.handle_operation_error(e, started),
        }
    }

    fn infer_export_format(path: &str) -> ExportFormat {
        let lower = path.to_ascii_lowercase();
        if lower.ends_with(".fa") || lower.ends_with(".fasta") {
            ExportFormat::Fasta
        } else {
            ExportFormat::GenBank
        }
    }

    fn export_active_sequence(&mut self) {
        let Some(seq_id) = self.seq_id.clone() else {
            self.op_status = "No active sequence to export".to_string();
            return;
        };
        if self.engine.is_none() {
            self.op_status = "No engine attached".to_string();
            return;
        }

        let default_name = format!("{seq_id}.gb");
        let path = rfd::FileDialog::new()
            .set_file_name(&default_name)
            .add_filter("GenBank", &["gb", "gbk", "genbank"])
            .add_filter("FASTA", &["fa", "fasta"])
            .save_file();
        let Some(path) = path else {
            self.op_status = "Export canceled".to_string();
            return;
        };
        let path = path.display().to_string();
        let format = Self::infer_export_format(&path);
        self.apply_operation_with_feedback(Operation::SaveFile {
            seq_id,
            path,
            format,
        });
    }

    fn export_active_sequence_svg(&mut self) {
        let Some(seq_id) = self.seq_id.clone() else {
            self.op_status = "No active sequence to export".to_string();
            return;
        };
        if self.engine.is_none() {
            self.op_status = "No engine attached".to_string();
            return;
        }
        let default_name = format!("{seq_id}.svg");
        let path = rfd::FileDialog::new()
            .set_file_name(&default_name)
            .add_filter("SVG", &["svg"])
            .save_file();
        let Some(path) = path else {
            self.op_status = "Export canceled".to_string();
            return;
        };
        let path = path.display().to_string();
        let mode = if self.is_circular() {
            RenderSvgMode::Circular
        } else {
            RenderSvgMode::Linear
        };
        self.apply_operation_with_feedback(Operation::RenderSequenceSvg { seq_id, mode, path });
    }

    fn export_pool_to_file(&mut self) {
        if self.engine.is_none() {
            self.op_status = "No engine attached".to_string();
            return;
        }
        let inputs = Self::parse_ids(&self.export_pool_inputs_text);
        if inputs.is_empty() {
            self.op_status = "Provide at least one sequence id for ExportPool".to_string();
            return;
        }
        let path = rfd::FileDialog::new()
            .set_file_name("pool.json")
            .add_filter("JSON", &["json"])
            .save_file();
        let Some(path) = path else {
            self.op_status = "Export canceled".to_string();
            return;
        };
        let path = path.display().to_string();
        self.apply_operation_with_feedback(Operation::ExportPool {
            inputs,
            path,
            pool_id: if self.export_pool_id.trim().is_empty() {
                None
            } else {
                Some(self.export_pool_id.trim().to_string())
            },
            human_id: if self.export_pool_human_id.trim().is_empty() {
                None
            } else {
                Some(self.export_pool_human_id.trim().to_string())
            },
        });
    }

    fn parse_ids(text: &str) -> Vec<String> {
        text.split(',')
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .map(|s| s.to_string())
            .collect()
    }

    fn collect_tfbs_motifs(&self) -> Vec<String> {
        if self.tfbs_use_all_motifs {
            tf_motifs::all_motif_ids()
        } else {
            Self::parse_ids(&self.tfbs_motifs)
        }
    }

    fn add_tfbs_motif_selection(&mut self, motif_id: &str) {
        let mut motifs = Self::parse_ids(&self.tfbs_motifs);
        if motifs
            .iter()
            .any(|m| m.eq_ignore_ascii_case(motif_id.trim()))
        {
            return;
        }
        motifs.push(motif_id.trim().to_string());
        self.tfbs_motifs = motifs.join(",");
    }

    fn parse_tf_threshold_map(text: &str) -> Result<HashMap<String, f64>, String> {
        let mut out = HashMap::new();
        for chunk in text.split(|c| matches!(c, ',' | ';' | '\n' | '\r')) {
            let chunk = chunk.trim();
            if chunk.is_empty() {
                continue;
            }
            let Some((k, v)) = chunk.split_once('=') else {
                return Err(format!(
                    "Invalid per-TF threshold entry '{}'; expected TF=VALUE",
                    chunk
                ));
            };
            let key = k.trim().to_string();
            if key.is_empty() {
                return Err("Empty TF key in per-TF threshold map".to_string());
            }
            let value = v
                .trim()
                .parse::<f64>()
                .map_err(|_| format!("Invalid numeric threshold value '{}'", v.trim()))?;
            out.insert(key, value);
        }
        Ok(out)
    }

    fn engine_ops_state_key(&self) -> String {
        let seq = self.seq_id.clone().unwrap_or_else(|| "_global".to_string());
        format!("gui.engine_ops.{seq}")
    }

    fn current_engine_ops_state(&self) -> EngineOpsUiState {
        let tfbs_display = self
            .dna_display
            .read()
            .expect("DNA display lock poisoned")
            .tfbs_display_criteria();
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
            extract_from: self.extract_from.clone(),
            extract_to: self.extract_to.clone(),
            extract_output_id: self.extract_output_id.clone(),
            parameter_name: self.parameter_name.clone(),
            parameter_value_json: self.parameter_value_json.clone(),
            anchored_mode_feature: self.anchored_mode_feature,
            anchored_position: self.anchored_position.clone(),
            anchored_feature_kind: self.anchored_feature_kind.clone(),
            anchored_feature_label: self.anchored_feature_label.clone(),
            anchored_feature_boundary_start: self.anchored_feature_boundary_start,
            anchored_feature_occurrence: self.anchored_feature_occurrence.clone(),
            anchored_direction_upstream: self.anchored_direction_upstream,
            anchored_target_len: self.anchored_target_len.clone(),
            anchored_tolerance: self.anchored_tolerance.clone(),
            anchored_required_re_sites: self.anchored_required_re_sites.clone(),
            anchored_required_tf_motifs: self.anchored_required_tf_motifs.clone(),
            anchored_forward_primer: self.anchored_forward_primer.clone(),
            anchored_reverse_primer: self.anchored_reverse_primer.clone(),
            anchored_output_prefix: self.anchored_output_prefix.clone(),
            anchored_unique: self.anchored_unique,
            anchored_max_candidates: self.anchored_max_candidates.clone(),
            tfbs_motifs: self.tfbs_motifs.clone(),
            tfbs_use_all_motifs: self.tfbs_use_all_motifs,
            tfbs_catalog_filter: self.tfbs_catalog_filter.clone(),
            tfbs_min_llr_bits: self.tfbs_min_llr_bits.clone(),
            tfbs_min_llr_quantile: self.tfbs_min_llr_quantile.clone(),
            tfbs_per_tf_min_llr_bits: self.tfbs_per_tf_min_llr_bits.clone(),
            tfbs_per_tf_min_llr_quantile: self.tfbs_per_tf_min_llr_quantile.clone(),
            tfbs_clear_existing: self.tfbs_clear_existing,
            tfbs_display_use_llr_bits: tfbs_display.use_llr_bits,
            tfbs_display_min_llr_bits: tfbs_display.min_llr_bits,
            tfbs_display_use_llr_quantile: tfbs_display.use_llr_quantile,
            tfbs_display_min_llr_quantile: tfbs_display.min_llr_quantile,
            tfbs_display_use_true_log_odds_bits: tfbs_display.use_true_log_odds_bits,
            tfbs_display_min_true_log_odds_bits: tfbs_display.min_true_log_odds_bits,
            tfbs_display_use_true_log_odds_quantile: tfbs_display.use_true_log_odds_quantile,
            tfbs_display_min_true_log_odds_quantile: tfbs_display.min_true_log_odds_quantile,
            container_digest_id: self.container_digest_id.clone(),
            container_merge_ids: self.container_merge_ids.clone(),
            container_ligation_id: self.container_ligation_id.clone(),
            container_ligation_output_prefix: self.container_ligation_output_prefix.clone(),
            container_ligation_protocol_sticky: self.container_ligation_protocol_sticky,
            container_ligation_circularize: self.container_ligation_circularize,
            container_ligation_unique: self.container_ligation_unique,
            container_mw_id: self.container_mw_id.clone(),
            container_mw_min_bp: self.container_mw_min_bp.clone(),
            container_mw_max_bp: self.container_mw_max_bp.clone(),
            container_mw_error: self.container_mw_error.clone(),
            container_mw_unique: self.container_mw_unique,
            container_mw_output_prefix: self.container_mw_output_prefix.clone(),
            workflow_run_id: self.workflow_run_id.clone(),
            workflow_ops_json: self.workflow_ops_json.clone(),
            export_pool_inputs_text: self.export_pool_inputs_text.clone(),
            export_pool_id: self.export_pool_id.clone(),
            export_pool_human_id: self.export_pool_human_id.clone(),
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
        self.extract_from = s.extract_from;
        self.extract_to = s.extract_to;
        self.extract_output_id = s.extract_output_id;
        self.parameter_name = s.parameter_name;
        self.parameter_value_json = s.parameter_value_json;
        self.anchored_mode_feature = s.anchored_mode_feature;
        self.anchored_position = s.anchored_position;
        self.anchored_feature_kind = s.anchored_feature_kind;
        self.anchored_feature_label = s.anchored_feature_label;
        self.anchored_feature_boundary_start = s.anchored_feature_boundary_start;
        self.anchored_feature_occurrence = s.anchored_feature_occurrence;
        self.anchored_direction_upstream = s.anchored_direction_upstream;
        self.anchored_target_len = s.anchored_target_len;
        self.anchored_tolerance = s.anchored_tolerance;
        self.anchored_required_re_sites = s.anchored_required_re_sites;
        self.anchored_required_tf_motifs = s.anchored_required_tf_motifs;
        self.anchored_forward_primer = s.anchored_forward_primer;
        self.anchored_reverse_primer = s.anchored_reverse_primer;
        self.anchored_output_prefix = s.anchored_output_prefix;
        self.anchored_unique = s.anchored_unique;
        self.anchored_max_candidates = s.anchored_max_candidates;
        self.tfbs_motifs = s.tfbs_motifs;
        self.tfbs_use_all_motifs = s.tfbs_use_all_motifs;
        self.tfbs_catalog_filter = s.tfbs_catalog_filter;
        self.tfbs_min_llr_bits = s.tfbs_min_llr_bits;
        self.tfbs_min_llr_quantile = s.tfbs_min_llr_quantile;
        self.tfbs_per_tf_min_llr_bits = s.tfbs_per_tf_min_llr_bits;
        self.tfbs_per_tf_min_llr_quantile = s.tfbs_per_tf_min_llr_quantile;
        self.tfbs_clear_existing = s.tfbs_clear_existing;
        self.container_digest_id = s.container_digest_id;
        self.container_merge_ids = s.container_merge_ids;
        self.container_ligation_id = s.container_ligation_id;
        self.container_ligation_output_prefix = s.container_ligation_output_prefix;
        self.container_ligation_protocol_sticky = s.container_ligation_protocol_sticky;
        self.container_ligation_circularize = s.container_ligation_circularize;
        self.container_ligation_unique = s.container_ligation_unique;
        self.container_mw_id = s.container_mw_id;
        self.container_mw_min_bp = s.container_mw_min_bp;
        self.container_mw_max_bp = s.container_mw_max_bp;
        self.container_mw_error = s.container_mw_error;
        self.container_mw_unique = s.container_mw_unique;
        self.container_mw_output_prefix = s.container_mw_output_prefix;
        self.workflow_run_id = s.workflow_run_id;
        self.workflow_ops_json = s.workflow_ops_json;
        self.export_pool_inputs_text = s.export_pool_inputs_text;
        self.export_pool_id = s.export_pool_id;
        self.export_pool_human_id = s.export_pool_human_id;
        self.dna_display
            .write()
            .expect("DNA display lock poisoned")
            .set_tfbs_display_criteria(TfbsDisplayCriteria {
                use_llr_bits: s.tfbs_display_use_llr_bits,
                min_llr_bits: s.tfbs_display_min_llr_bits,
                use_llr_quantile: s.tfbs_display_use_llr_quantile,
                min_llr_quantile: s.tfbs_display_min_llr_quantile,
                use_true_log_odds_bits: s.tfbs_display_use_true_log_odds_bits,
                min_true_log_odds_bits: s.tfbs_display_min_true_log_odds_bits,
                use_true_log_odds_quantile: s.tfbs_display_use_true_log_odds_quantile,
                min_true_log_odds_quantile: s.tfbs_display_min_true_log_odds_quantile,
            });
        self.sync_tfbs_display_criteria_to_engine(TfbsDisplayCriteria {
            use_llr_bits: s.tfbs_display_use_llr_bits,
            min_llr_bits: s.tfbs_display_min_llr_bits,
            use_llr_quantile: s.tfbs_display_use_llr_quantile,
            min_llr_quantile: s.tfbs_display_min_llr_quantile,
            use_true_log_odds_bits: s.tfbs_display_use_true_log_odds_bits,
            min_true_log_odds_bits: s.tfbs_display_min_true_log_odds_bits,
            use_true_log_odds_quantile: s.tfbs_display_use_true_log_odds_quantile,
            min_true_log_odds_quantile: s.tfbs_display_min_true_log_odds_quantile,
        });
    }

    fn save_engine_ops_state(&self) {
        let Some(engine) = &self.engine else {
            return;
        };
        let key = self.engine_ops_state_key();
        let value = match serde_json::to_value(self.current_engine_ops_state()) {
            Ok(v) => v,
            Err(_) => return,
        };
        let mut guard = engine.write().expect("Engine lock poisoned");
        if guard.state().metadata.get(&key) == Some(&value) {
            return;
        }
        guard.state_mut().metadata.insert(key, value);
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
