use crate::{
    dna_display::{DnaDisplay, Selection, TfbsDisplayCriteria, VcfDisplayCriteria},
    dna_sequence::DNAsequence,
    engine::{
        AnchorBoundary, AnchorDirection, AnchoredRegionAnchor, CandidateFeatureStrandRelation,
        CandidateRecord, CandidateSetOperator, DisplayTarget, Engine, EngineError, ErrorCode,
        ExportFormat, GentleEngine, LigationProtocol, OpResult, Operation, OperationProgress,
        PcrPrimerSpec, RenderSvgMode, SnpMutationSpec, TfThresholdOverride, TfbsProgress, Workflow,
    },
    engine_shell::{
        execute_shell_command_with_options, parse_shell_line, shell_help_text, ShellCommand,
        ShellExecutionOptions,
    },
    feature_location::collect_location_ranges_usize,
    icons::*,
    open_reading_frame::OpenReadingFrame,
    pool_gel::build_pool_gel_layout,
    render_dna::RenderDna,
    render_sequence::RenderSequence,
    tf_motifs,
};
use eframe::egui::{self, Frame, PointerState, Vec2};
use serde::{Deserialize, Serialize};
use std::{
    cmp::Ordering,
    collections::{BTreeSet, HashMap},
    fs,
    path::PathBuf,
    sync::{
        mpsc::{self, Receiver, TryRecvError},
        Arc, Mutex, RwLock,
    },
    time::{Duration, Instant},
};

#[derive(Clone, Debug, Serialize, Deserialize)]
struct EngineOpsUiState {
    show_engine_ops: bool,
    #[serde(default)]
    show_shell: bool,
    #[serde(default)]
    shell_command_text: String,
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
    #[serde(default)]
    sq_inputs_text: String,
    #[serde(default)]
    sq_gc_min: String,
    #[serde(default)]
    sq_gc_max: String,
    #[serde(default)]
    sq_max_homopolymer_run: String,
    #[serde(default = "default_true")]
    sq_reject_ambiguous_bases: bool,
    #[serde(default = "default_true")]
    sq_avoid_u6_tttt: bool,
    #[serde(default)]
    sq_forbidden_motifs: String,
    #[serde(default)]
    sq_unique: bool,
    #[serde(default)]
    sq_output_prefix: String,
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
    candidate_set_name: String,
    #[serde(default)]
    candidate_source_seq_id: String,
    #[serde(default)]
    candidate_length_bp: String,
    #[serde(default)]
    candidate_step_bp: String,
    #[serde(default)]
    candidate_feature_kinds: String,
    #[serde(default)]
    candidate_feature_label_regex: String,
    #[serde(default)]
    candidate_feature_strand_relation: CandidateFeatureStrandRelation,
    #[serde(default)]
    candidate_max_distance_bp: String,
    #[serde(default)]
    candidate_limit: String,
    #[serde(default)]
    candidate_selected_set: String,
    #[serde(default)]
    candidate_page_limit: String,
    #[serde(default)]
    candidate_page_offset: String,
    #[serde(default)]
    candidate_sort_key: String,
    #[serde(default)]
    candidate_sort_desc: bool,
    #[serde(default)]
    candidate_score_metric: String,
    #[serde(default)]
    candidate_score_expression: String,
    #[serde(default)]
    candidate_distance_metric: String,
    #[serde(default)]
    candidate_distance_feature_kinds: String,
    #[serde(default)]
    candidate_distance_feature_label_regex: String,
    #[serde(default)]
    candidate_distance_feature_strand_relation: CandidateFeatureStrandRelation,
    #[serde(default)]
    candidate_filter_input_set: String,
    #[serde(default)]
    candidate_filter_output_set: String,
    #[serde(default)]
    candidate_filter_metric: String,
    #[serde(default)]
    candidate_filter_min: String,
    #[serde(default)]
    candidate_filter_max: String,
    #[serde(default)]
    candidate_filter_min_quantile: String,
    #[serde(default)]
    candidate_filter_max_quantile: String,
    #[serde(default)]
    candidate_setop_mode: String,
    #[serde(default)]
    candidate_setop_left: String,
    #[serde(default)]
    candidate_setop_right: String,
    #[serde(default)]
    candidate_setop_output: String,
    #[serde(default)]
    candidate_macro_script: String,
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
    #[serde(default = "default_true")]
    vcf_display_show_snp: bool,
    #[serde(default = "default_true")]
    vcf_display_show_ins: bool,
    #[serde(default = "default_true")]
    vcf_display_show_del: bool,
    #[serde(default = "default_true")]
    vcf_display_show_sv: bool,
    #[serde(default = "default_true")]
    vcf_display_show_other: bool,
    #[serde(default)]
    vcf_display_pass_only: bool,
    #[serde(default)]
    vcf_display_use_min_qual: bool,
    #[serde(default = "default_zero_f64")]
    vcf_display_min_qual: f64,
    #[serde(default)]
    vcf_display_use_max_qual: bool,
    #[serde(default = "default_zero_f64")]
    vcf_display_max_qual: f64,
    #[serde(default)]
    vcf_display_required_info_keys: String,
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
    #[serde(default)]
    pool_gel_ladders: String,
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

const TOP_PANEL_ICON_SIZE_PX: f32 = 20.0;
const UI_SIZE_MIN_PX: f32 = 1.0;
const UI_SIZE_MAX_PX: f32 = 4096.0;
const DECLUTTER_NOISE_SCORE_THRESHOLD: usize = 100;
const DECLUTTER_VISIBLE_FEATURE_THRESHOLD: usize = 70;
const POOL_GEL_LADDER_PRESETS: [(&str, &str); 5] = [
    ("Auto", ""),
    ("NEB 100bp + 1kb", "NEB 100bp DNA Ladder,NEB 1kb DNA Ladder"),
    (
        "MassRuler Low + High",
        "MassRuler Low Range,MassRuler High Range",
    ),
    ("GeneRuler 50bp + Mix", "GeneRuler 50bp,GeneRuler Mix"),
    (
        "GeneRuler 100bp+ + Mix",
        "GeneRuler 100bp DNA Ladder Plus,GeneRuler Mix",
    ),
];

#[derive(Clone, Copy, Debug)]
enum MapViewPreset {
    Anchored,
    Cloning,
    Annotation,
    Signal,
}

impl MapViewPreset {
    fn label(self) -> &'static str {
        match self {
            Self::Anchored => "Anchored",
            Self::Cloning => "Cloning",
            Self::Annotation => "Annotation",
            Self::Signal => "Signal",
        }
    }

    fn hover_text(self) -> &'static str {
        match self {
            Self::Anchored => {
                "Prioritize anchored-gene readability (core annotation lanes, low clutter overlays)"
            }
            Self::Cloning => "Prioritize cloning context (core features + restriction/GC overlays)",
            Self::Annotation => {
                "Prioritize annotation review (core features + TFBS, minimal auxiliary overlays)"
            }
            Self::Signal => {
                "Prioritize signal tracks and regulatory context with reduced coding-feature noise"
            }
        }
    }
}

#[derive(Clone, Debug, Default)]
struct LayerVisibilityCounts {
    feature_kind_counts: HashMap<String, usize>,
    regulatory_feature_count: usize,
    restriction_site_count: usize,
    gc_region_count: usize,
    orf_count: usize,
    methylation_site_count: usize,
}

impl LayerVisibilityCounts {
    fn kind_count(&self, kind: &str) -> usize {
        self.feature_kind_counts
            .get(&kind.trim().to_ascii_uppercase())
            .copied()
            .unwrap_or(0)
    }

    fn tfbs_count(&self) -> usize {
        self.kind_count("TFBS")
            + self.kind_count("TF_BINDING_SITE")
            + self.kind_count("PROTEIN_BIND")
    }

    fn total_feature_count(&self) -> usize {
        self.feature_kind_counts.values().sum()
    }
}

#[derive(Clone, Debug)]
struct DeclutterSnapshot {
    show_tfbs: bool,
    show_restriction_enzymes: bool,
    show_gc_contents: bool,
    show_open_reading_frames: bool,
    show_methylation_sites: bool,
    regulatory_tracks_near_baseline: bool,
    hidden_feature_kinds: BTreeSet<String>,
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
    show_shell: bool,
    shell_command_text: String,
    shell_preview_text: String,
    shell_output_text: String,
    rna_info_text: String,
    rna_svg_uri: String,
    rna_status: String,
    rna_preview_path: Option<PathBuf>,
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
    sq_inputs_text: String,
    sq_gc_min: String,
    sq_gc_max: String,
    sq_max_homopolymer_run: String,
    sq_reject_ambiguous_bases: bool,
    sq_avoid_u6_tttt: bool,
    sq_forbidden_motifs: String,
    sq_unique: bool,
    sq_output_prefix: String,
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
    candidate_set_name: String,
    candidate_source_seq_id: String,
    candidate_length_bp: String,
    candidate_step_bp: String,
    candidate_feature_kinds: String,
    candidate_feature_label_regex: String,
    candidate_feature_strand_relation: CandidateFeatureStrandRelation,
    candidate_max_distance_bp: String,
    candidate_limit: String,
    candidate_selected_set: String,
    candidate_page_limit: String,
    candidate_page_offset: String,
    candidate_sort_key: String,
    candidate_sort_desc: bool,
    candidate_score_metric: String,
    candidate_score_expression: String,
    candidate_distance_metric: String,
    candidate_distance_feature_kinds: String,
    candidate_distance_feature_label_regex: String,
    candidate_distance_feature_strand_relation: CandidateFeatureStrandRelation,
    candidate_filter_input_set: String,
    candidate_filter_output_set: String,
    candidate_filter_metric: String,
    candidate_filter_min: String,
    candidate_filter_max: String,
    candidate_filter_min_quantile: String,
    candidate_filter_max_quantile: String,
    candidate_setop_mode: String,
    candidate_setop_left: String,
    candidate_setop_right: String,
    candidate_setop_output: String,
    candidate_macro_script: String,
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
    vcf_display_required_info_keys: String,
    declutter_snapshot: Option<DeclutterSnapshot>,
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
    pool_gel_ladders: String,
    pending_feature_tree_scroll_to: Option<usize>,
    focused_feature_id: Option<usize>,
    description_cache_initialized: bool,
    description_cache_selected_id: Option<usize>,
    description_cache_seq_len: usize,
    description_cache_feature_count: usize,
    description_cache_title: String,
    description_cache_range: Option<String>,
}

impl MainAreaDna {
    fn sanitize_widget_size(size: Vec2) -> Vec2 {
        let sanitize = |value: f32| -> f32 {
            if value.is_finite() {
                value.clamp(UI_SIZE_MIN_PX, UI_SIZE_MAX_PX)
            } else {
                TOP_PANEL_ICON_SIZE_PX
            }
        };
        Vec2::new(sanitize(size.x), sanitize(size.y))
    }

    fn top_panel_icon_size(ui: &egui::Ui) -> Vec2 {
        let suggested = ui.spacing().interact_size.y;
        let edge = if suggested.is_finite() {
            suggested.clamp(UI_SIZE_MIN_PX, TOP_PANEL_ICON_SIZE_PX)
        } else {
            TOP_PANEL_ICON_SIZE_PX
        };
        Self::sanitize_widget_size(Vec2::new(edge, edge))
    }

    fn feature_kinds_for_toggle_buttons(&self) -> Vec<String> {
        let mut kinds = BTreeSet::new();
        if let Ok(dna) = self.dna.read() {
            for feature in dna.features() {
                let kind = feature.kind.to_string().trim().to_ascii_uppercase();
                if !kind.is_empty() && kind != "SOURCE" {
                    kinds.insert(kind);
                }
            }
        }
        kinds.into_iter().collect()
    }

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
            show_shell: false,
            shell_command_text: "state-summary".to_string(),
            shell_preview_text: String::new(),
            shell_output_text: String::new(),
            rna_info_text: String::new(),
            rna_svg_uri: String::new(),
            rna_status: String::new(),
            rna_preview_path: None,
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
            sq_inputs_text: seq_id_for_defaults.clone().unwrap_or_default(),
            sq_gc_min: "0.30".to_string(),
            sq_gc_max: "0.70".to_string(),
            sq_max_homopolymer_run: "4".to_string(),
            sq_reject_ambiguous_bases: true,
            sq_avoid_u6_tttt: true,
            sq_forbidden_motifs: String::new(),
            sq_unique: false,
            sq_output_prefix: "design".to_string(),
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
            candidate_set_name: "candidates".to_string(),
            candidate_source_seq_id: seq_id_for_defaults.clone().unwrap_or_default(),
            candidate_length_bp: "20".to_string(),
            candidate_step_bp: "1".to_string(),
            candidate_feature_kinds: String::new(),
            candidate_feature_label_regex: String::new(),
            candidate_feature_strand_relation: CandidateFeatureStrandRelation::Any,
            candidate_max_distance_bp: String::new(),
            candidate_limit: "5000".to_string(),
            candidate_selected_set: String::new(),
            candidate_page_limit: "50".to_string(),
            candidate_page_offset: "0".to_string(),
            candidate_sort_key: "start".to_string(),
            candidate_sort_desc: false,
            candidate_score_metric: "score".to_string(),
            candidate_score_expression: "gc_fraction".to_string(),
            candidate_distance_metric: "distance_to_gene_bp".to_string(),
            candidate_distance_feature_kinds: "gene".to_string(),
            candidate_distance_feature_label_regex: String::new(),
            candidate_distance_feature_strand_relation: CandidateFeatureStrandRelation::Any,
            candidate_filter_input_set: String::new(),
            candidate_filter_output_set: "filtered".to_string(),
            candidate_filter_metric: "gc_fraction".to_string(),
            candidate_filter_min: String::new(),
            candidate_filter_max: String::new(),
            candidate_filter_min_quantile: String::new(),
            candidate_filter_max_quantile: String::new(),
            candidate_setop_mode: "intersect".to_string(),
            candidate_setop_left: String::new(),
            candidate_setop_right: String::new(),
            candidate_setop_output: "set_op".to_string(),
            candidate_macro_script: "generate my_set seq_1 --length 20 --step 1 --limit 1000"
                .to_string(),
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
            vcf_display_required_info_keys: String::new(),
            declutter_snapshot: None,
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
            pool_gel_ladders: String::new(),
            pending_feature_tree_scroll_to: None,
            focused_feature_id: None,
            description_cache_initialized: false,
            description_cache_selected_id: None,
            description_cache_seq_len: 0,
            description_cache_feature_count: 0,
            description_cache_title: String::new(),
            description_cache_range: None,
        };
        ret.sync_from_engine_display();
        ret.load_engine_ops_state();
        ret
    }

    pub fn dna(&self) -> &Arc<RwLock<DNAsequence>> {
        &self.dna
    }

    pub fn set_pool_context(&mut self, pool_seq_ids: Vec<String>) {
        self.last_created_seq_ids = pool_seq_ids.clone();
        self.export_pool_inputs_text = pool_seq_ids.join(", ");
        self.show_engine_ops = true;
        self.op_status = "Opened from lineage pool node".to_string();
    }

    pub fn refresh_from_engine_settings(&mut self) {
        self.sync_from_engine_display();
        self.update_dna_map();
    }

    fn active_sequence_is_genome_anchored(&self, engine: &GentleEngine) -> bool {
        let Some(seq_id) = self.seq_id.as_deref() else {
            return false;
        };
        engine.describe_sequence_genome_anchor(seq_id).is_ok()
    }

    fn active_linear_viewport_range(&self) -> Option<(usize, usize)> {
        if self.is_circular() {
            return None;
        }
        let (start, span, sequence_length) = self.current_linear_viewport();
        if sequence_length == 0 || span == 0 {
            return None;
        }
        let end = start.saturating_add(span).min(sequence_length);
        Some((start, end))
    }

    fn ranges_overlap(
        a_start: usize,
        a_end_exclusive: usize,
        b_start: usize,
        b_end_exclusive: usize,
    ) -> bool {
        a_start < b_end_exclusive && a_end_exclusive > b_start
    }

    fn feature_overlaps_linear_viewport(
        feature: &gb_io::seq::Feature,
        sequence_length: usize,
        viewport_start: usize,
        viewport_end: usize,
    ) -> bool {
        let mut ranges: Vec<(usize, usize)> = vec![];
        collect_location_ranges_usize(&feature.location, &mut ranges);
        if ranges.is_empty() {
            if let Ok((raw_from, raw_to)) = feature.location.find_bounds() {
                if raw_from >= 0 && raw_to >= 0 {
                    let mut from = raw_from as usize;
                    let mut to = raw_to as usize;
                    if to < from {
                        std::mem::swap(&mut from, &mut to);
                    }
                    ranges.push((from, to));
                }
            }
        }
        ranges.into_iter().any(|(mut from, mut to)| {
            if to < from {
                std::mem::swap(&mut from, &mut to);
            }
            if from >= sequence_length {
                return false;
            }
            to = to.min(sequence_length);
            if to <= from {
                return false;
            }
            Self::ranges_overlap(from, to, viewport_start, viewport_end)
        })
    }

    fn orf_overlaps_linear_viewport(
        orf: &OpenReadingFrame,
        sequence_length: usize,
        viewport_start: usize,
        viewport_end: usize,
    ) -> bool {
        if sequence_length == 0 {
            return false;
        }
        let from = orf.from().max(0) as usize;
        let to = orf.to().max(0) as usize;
        if to >= from {
            let from = from.min(sequence_length);
            let to = to.min(sequence_length);
            if to <= from {
                return false;
            }
            return Self::ranges_overlap(from, to, viewport_start, viewport_end);
        }
        let wrapped_left = from.min(sequence_length);
        let wrapped_right = to.min(sequence_length);
        Self::ranges_overlap(wrapped_left, sequence_length, viewport_start, viewport_end)
            || Self::ranges_overlap(0, wrapped_right, viewport_start, viewport_end)
    }

    fn compute_layer_visibility_counts(&self) -> LayerVisibilityCounts {
        let viewport = self.active_linear_viewport_range();
        let (tfbs_display_criteria, vcf_display_criteria) = self
            .dna_display
            .read()
            .map(|display| {
                (
                    display.tfbs_display_criteria(),
                    display.vcf_display_criteria(),
                )
            })
            .unwrap_or((
                TfbsDisplayCriteria::default(),
                VcfDisplayCriteria::default(),
            ));

        let mut counts = LayerVisibilityCounts::default();
        if let Ok(dna) = self.dna.read() {
            let sequence_length = dna.len();
            for feature in dna.features() {
                if RenderDna::is_source_feature(feature) {
                    continue;
                }
                if let Some((start, end)) = viewport {
                    if !Self::feature_overlaps_linear_viewport(feature, sequence_length, start, end)
                    {
                        continue;
                    }
                }
                let kind = feature.kind.to_string().trim().to_ascii_uppercase();
                if kind.is_empty() {
                    continue;
                }
                if RenderDna::is_tfbs_feature(feature)
                    && !RenderDna::tfbs_feature_passes_display_filter(
                        feature,
                        tfbs_display_criteria,
                    )
                {
                    continue;
                }
                if RenderDna::is_vcf_track_feature(feature)
                    && !RenderDna::vcf_feature_passes_display_filter(feature, &vcf_display_criteria)
                {
                    continue;
                }
                *counts.feature_kind_counts.entry(kind).or_insert(0) += 1;
                if RenderDna::is_regulatory_feature(feature) {
                    counts.regulatory_feature_count += 1;
                }
            }

            counts.restriction_site_count = dna
                .restriction_enzyme_sites()
                .iter()
                .filter(|site| {
                    let offset = match usize::try_from(site.offset.max(0)) {
                        Ok(v) => v,
                        Err(_) => return false,
                    };
                    match viewport {
                        Some((start, end)) => offset >= start && offset < end,
                        None => true,
                    }
                })
                .count();
            counts.gc_region_count = dna
                .gc_content()
                .regions()
                .iter()
                .filter(|region| match viewport {
                    Some((start, end)) => {
                        Self::ranges_overlap(region.from(), region.to(), start, end)
                    }
                    None => true,
                })
                .count();
            counts.orf_count = dna
                .open_reading_frames()
                .iter()
                .filter(|orf| match viewport {
                    Some((start, end)) => {
                        Self::orf_overlaps_linear_viewport(orf, sequence_length, start, end)
                    }
                    None => true,
                })
                .count();
            counts.methylation_site_count = dna
                .methylation_sites()
                .sites()
                .iter()
                .filter(|site| match viewport {
                    Some((start, end)) => **site >= start && **site < end,
                    None => true,
                })
                .count();
        }
        counts
    }

    fn apply_display_preset_visibility(
        &mut self,
        show_features: bool,
        show_cds_features: bool,
        show_gene_features: bool,
        show_mrna_features: bool,
        show_tfbs: bool,
        show_restriction_enzymes: bool,
        show_gc_contents: bool,
        show_open_reading_frames: bool,
        show_methylation_sites: bool,
        regulatory_tracks_near_baseline: bool,
        hidden_feature_kinds: BTreeSet<String>,
    ) {
        if let Ok(mut display) = self.dna_display.write() {
            display.set_show_features(show_features);
            display.set_show_cds_features(show_cds_features);
            display.set_show_gene_features(show_gene_features);
            display.set_show_mrna_features(show_mrna_features);
            display.set_show_tfbs(show_tfbs);
            display.set_show_restriction_enzyme_sites(show_restriction_enzymes);
            display.set_show_gc_contents(show_gc_contents);
            display.set_show_open_reading_frames(show_open_reading_frames);
            display.set_show_methylation_sites(show_methylation_sites);
            display.set_regulatory_tracks_near_baseline(regulatory_tracks_near_baseline);
            display.set_hidden_feature_kinds(hidden_feature_kinds);
        }
        self.set_display_visibility(DisplayTarget::Features, show_features);
        self.set_display_visibility(DisplayTarget::CdsFeatures, show_cds_features);
        self.set_display_visibility(DisplayTarget::GeneFeatures, show_gene_features);
        self.set_display_visibility(DisplayTarget::MrnaFeatures, show_mrna_features);
        self.set_display_visibility(DisplayTarget::Tfbs, show_tfbs);
        self.set_display_visibility(DisplayTarget::RestrictionEnzymes, show_restriction_enzymes);
        self.set_display_visibility(DisplayTarget::GcContents, show_gc_contents);
        self.set_display_visibility(DisplayTarget::OpenReadingFrames, show_open_reading_frames);
        self.set_display_visibility(DisplayTarget::MethylationSites, show_methylation_sites);
        self.sync_regulatory_track_placement_to_engine(regulatory_tracks_near_baseline);
    }

    fn preset_keeps_feature_kind(preset: MapViewPreset, kind: &str) -> bool {
        let kind = kind.trim().to_ascii_uppercase();
        if kind.is_empty() {
            return true;
        }
        match preset {
            MapViewPreset::Anchored => {
                matches!(kind.as_str(), "CDS" | "GENE" | "MRNA")
                    || kind.contains("REGULATORY")
                    || matches!(
                        kind.as_str(),
                        "TRACK" | "TFBS" | "TF_BINDING_SITE" | "PROTEIN_BIND"
                    )
            }
            MapViewPreset::Cloning | MapViewPreset::Annotation => true,
            MapViewPreset::Signal => {
                kind == "GENE"
                    || kind.contains("REGULATORY")
                    || matches!(
                        kind.as_str(),
                        "TRACK" | "TFBS" | "TF_BINDING_SITE" | "PROTEIN_BIND"
                    )
            }
        }
    }

    fn hidden_feature_kinds_for_preset(&self, preset: MapViewPreset) -> BTreeSet<String> {
        self.feature_kinds_for_toggle_buttons()
            .into_iter()
            .filter(|kind| !Self::preset_keeps_feature_kind(preset, kind))
            .collect()
    }

    fn apply_map_view_preset(&mut self, preset: MapViewPreset) {
        let hidden_feature_kinds = self.hidden_feature_kinds_for_preset(preset);
        match preset {
            MapViewPreset::Anchored => self.apply_display_preset_visibility(
                true,
                true,
                true,
                true,
                false,
                false,
                true,
                false,
                false,
                false,
                hidden_feature_kinds,
            ),
            MapViewPreset::Cloning => self.apply_display_preset_visibility(
                true,
                true,
                true,
                true,
                false,
                true,
                true,
                false,
                false,
                false,
                hidden_feature_kinds,
            ),
            MapViewPreset::Annotation => self.apply_display_preset_visibility(
                true,
                true,
                true,
                true,
                true,
                false,
                false,
                false,
                false,
                false,
                hidden_feature_kinds,
            ),
            MapViewPreset::Signal => self.apply_display_preset_visibility(
                true,
                false,
                true,
                false,
                false,
                false,
                true,
                false,
                false,
                true,
                hidden_feature_kinds,
            ),
        }
        self.declutter_snapshot = None;
        self.op_status = format!("Applied {} view preset", preset.label());
    }

    fn is_low_value_feature_kind(kind: &str) -> bool {
        !matches!(kind, "CDS" | "GENE" | "MRNA")
    }

    fn visible_feature_noise_metrics(
        &self,
        counts: &LayerVisibilityCounts,
    ) -> (usize, usize, BTreeSet<String>) {
        let (
            show_features,
            show_cds,
            show_gene,
            show_mrna,
            show_tfbs,
            show_restriction_enzymes,
            show_open_reading_frames,
            show_methylation_sites,
            show_gc_contents,
            suppress_orfs_for_anchor,
            hidden_feature_kinds,
        ) = self
            .dna_display
            .read()
            .map(|display| {
                (
                    display.show_features(),
                    display.show_cds_features(),
                    display.show_gene_features(),
                    display.show_mrna_features(),
                    display.show_tfbs(),
                    display.show_restriction_enzyme_sites(),
                    display.show_open_reading_frames(),
                    display.show_methylation_sites(),
                    display.show_gc_contents(),
                    display.suppress_open_reading_frames_for_genome_anchor(),
                    display.hidden_feature_kinds().clone(),
                )
            })
            .unwrap_or((
                true,
                true,
                true,
                true,
                false,
                true,
                false,
                false,
                true,
                false,
                BTreeSet::new(),
            ));
        if !show_features {
            return (0, 0, hidden_feature_kinds);
        }

        let mut visible_feature_count = 0usize;
        let mut low_value_visible_count = 0usize;
        for (kind, count) in &counts.feature_kind_counts {
            if hidden_feature_kinds.contains(kind) {
                continue;
            }
            let visible = match kind.as_str() {
                "CDS" => show_cds,
                "GENE" => show_gene,
                "MRNA" => show_mrna,
                "TFBS" | "TF_BINDING_SITE" | "PROTEIN_BIND" => show_tfbs,
                _ => true,
            };
            if !visible {
                continue;
            }
            visible_feature_count = visible_feature_count.saturating_add(*count);
            if Self::is_low_value_feature_kind(kind) {
                low_value_visible_count = low_value_visible_count.saturating_add(*count);
            }
        }

        let mut noise_score = low_value_visible_count;
        if show_tfbs {
            noise_score = noise_score.saturating_add(counts.tfbs_count());
        }
        if show_restriction_enzymes {
            noise_score = noise_score.saturating_add(counts.restriction_site_count / 2);
        }
        if show_methylation_sites {
            noise_score = noise_score.saturating_add(counts.methylation_site_count / 3);
        }
        if show_open_reading_frames && !suppress_orfs_for_anchor {
            noise_score = noise_score.saturating_add(counts.orf_count / 2);
        }
        if show_gc_contents {
            noise_score = noise_score.saturating_add(counts.gc_region_count / 8);
        }
        (visible_feature_count, noise_score, hidden_feature_kinds)
    }

    fn apply_declutter_action(&mut self, counts: &LayerVisibilityCounts) {
        if let Some(snapshot) = self.declutter_snapshot.take() {
            if let Ok(mut display) = self.dna_display.write() {
                display.set_show_tfbs(snapshot.show_tfbs);
                display.set_show_restriction_enzyme_sites(snapshot.show_restriction_enzymes);
                display.set_show_gc_contents(snapshot.show_gc_contents);
                display.set_show_open_reading_frames(snapshot.show_open_reading_frames);
                display.set_show_methylation_sites(snapshot.show_methylation_sites);
                display
                    .set_regulatory_tracks_near_baseline(snapshot.regulatory_tracks_near_baseline);
                display.set_hidden_feature_kinds(snapshot.hidden_feature_kinds);
            }
            self.set_display_visibility(DisplayTarget::Tfbs, snapshot.show_tfbs);
            self.set_display_visibility(
                DisplayTarget::RestrictionEnzymes,
                snapshot.show_restriction_enzymes,
            );
            self.set_display_visibility(DisplayTarget::GcContents, snapshot.show_gc_contents);
            self.set_display_visibility(
                DisplayTarget::OpenReadingFrames,
                snapshot.show_open_reading_frames,
            );
            self.set_display_visibility(
                DisplayTarget::MethylationSites,
                snapshot.show_methylation_sites,
            );
            self.sync_regulatory_track_placement_to_engine(
                snapshot.regulatory_tracks_near_baseline,
            );
            self.op_status = "Restored view after declutter".to_string();
            return;
        }

        let (visible_feature_count, noise_score, current_hidden_feature_kinds) =
            self.visible_feature_noise_metrics(counts);
        let should_declutter = visible_feature_count >= DECLUTTER_VISIBLE_FEATURE_THRESHOLD
            || noise_score >= DECLUTTER_NOISE_SCORE_THRESHOLD;
        if !should_declutter {
            self.op_status = format!(
                "Declutter skipped: visible features={} noise score={} (below thresholds)",
                visible_feature_count, noise_score
            );
            return;
        }

        let snapshot = self
            .dna_display
            .read()
            .ok()
            .map(|display| DeclutterSnapshot {
                show_tfbs: display.show_tfbs(),
                show_restriction_enzymes: display.show_restriction_enzyme_sites(),
                show_gc_contents: display.show_gc_contents(),
                show_open_reading_frames: display.show_open_reading_frames(),
                show_methylation_sites: display.show_methylation_sites(),
                regulatory_tracks_near_baseline: display.regulatory_tracks_near_baseline(),
                hidden_feature_kinds: display.hidden_feature_kinds().clone(),
            });
        let Some(snapshot) = snapshot else {
            self.op_status = "Declutter unavailable: could not read display state".to_string();
            return;
        };

        let mut hidden_feature_kinds = current_hidden_feature_kinds;
        for kind in counts.feature_kind_counts.keys() {
            if Self::is_low_value_feature_kind(kind) {
                hidden_feature_kinds.insert(kind.clone());
            }
        }
        if let Ok(mut display) = self.dna_display.write() {
            display.set_show_tfbs(false);
            display.set_show_restriction_enzyme_sites(false);
            display.set_show_open_reading_frames(false);
            display.set_show_methylation_sites(false);
            display.set_regulatory_tracks_near_baseline(false);
            display.set_hidden_feature_kinds(hidden_feature_kinds);
        }
        self.set_display_visibility(DisplayTarget::Tfbs, false);
        self.set_display_visibility(DisplayTarget::RestrictionEnzymes, false);
        self.set_display_visibility(DisplayTarget::OpenReadingFrames, false);
        self.set_display_visibility(DisplayTarget::MethylationSites, false);
        self.sync_regulatory_track_placement_to_engine(false);
        self.declutter_snapshot = Some(snapshot);
        self.op_status = format!(
            "Declutter applied: visible features={} noise score={}",
            visible_feature_count, noise_score
        );
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
        ui.horizontal_wrapped(|ui| {
            let icon_size = Self::top_panel_icon_size(ui);
            let layer_counts = self.compute_layer_visibility_counts();
            let button = egui::ImageButton::new(
                ICON_CIRCULAR_LINEAR
                    .clone()
                    .fit_to_exact_size(icon_size)
                    .rounding(5.0),
            );
            let response = ui
                .add(button)
                .on_hover_text("Toggle DNA map topology: circular <-> linear");
            if response.clicked() {
                if let (Some(engine), Some(seq_id)) = (self.engine.clone(), self.seq_id.clone()) {
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
                            .get(&seq_id)
                            .cloned()
                        {
                            *self.dna.write().expect("DNA lock poisoned") = updated;
                            self.clear_feature_focus();
                        }
                    }
                } else {
                    let mut dna = self.dna.write().expect("DNA lock poisoned");
                    let is_circular = dna.is_circular();
                    dna.set_circular(!is_circular);
                }
            }

            let button = egui::ImageButton::new(
                ICON_SHOW_SEQUENCE
                    .clone()
                    .fit_to_exact_size(icon_size)
                    .rounding(5.0),
            )
            .selected(self.show_sequence);
            let response = ui
                .add(button)
                .on_hover_text("Show or hide the sequence text panel");
            if response.clicked() {
                self.show_sequence = !self.show_sequence;
                self.set_display_visibility(DisplayTarget::SequencePanel, self.show_sequence);
            };

            let button = egui::ImageButton::new(
                ICON_SHOW_MAP
                    .clone()
                    .fit_to_exact_size(icon_size)
                    .rounding(5.0),
            )
            .selected(self.show_map);
            let response = ui
                .add(button)
                .on_hover_text("Show or hide the DNA map panel");
            if response.clicked() {
                self.show_map = !self.show_map;
                self.set_display_visibility(DisplayTarget::MapPanel, self.show_map);
            };

            if !self.is_circular() {
                let (start_bp, span_bp, sequence_length) = self.current_linear_viewport();
                if sequence_length > 0 {
                    ui.separator();
                    if ui
                        .small_button("-")
                        .on_hover_text("Zoom out linear map")
                        .clicked()
                    {
                        self.zoom_linear_viewport_around(
                            start_bp.saturating_add(span_bp / 2),
                            false,
                        );
                    }
                    if ui
                        .small_button("+")
                        .on_hover_text("Zoom in linear map")
                        .clicked()
                    {
                        self.zoom_linear_viewport_around(
                            start_bp.saturating_add(span_bp / 2),
                            true,
                        );
                    }
                    if ui
                        .small_button("Fit")
                        .on_hover_text("Reset linear map to full sequence")
                        .clicked()
                    {
                        self.set_linear_viewport(0, sequence_length);
                    }
                    let mut pan_start = start_bp;
                    let max_start = sequence_length.saturating_sub(span_bp);
                    let pan = ui
                        .add_enabled(
                            max_start > 0,
                            egui::Slider::new(&mut pan_start, 0..=max_start).text("Pan"),
                        )
                        .on_hover_text("Pan linear map viewport left/right");
                    if pan.changed() {
                        self.set_linear_viewport(pan_start, span_bp);
                    }
                    let view_end = start_bp.saturating_add(span_bp).min(sequence_length);
                    ui.monospace(format!(
                        "view {}..{} ({} bp)",
                        start_bp.saturating_add(1),
                        view_end,
                        span_bp
                    ));
                }
            }

            ui.separator();
            for preset in [
                MapViewPreset::Anchored,
                MapViewPreset::Cloning,
                MapViewPreset::Annotation,
                MapViewPreset::Signal,
            ] {
                if ui
                    .small_button(preset.label())
                    .on_hover_text(preset.hover_text())
                    .clicked()
                {
                    self.apply_map_view_preset(preset);
                }
            }
            let declutter_label = if self.declutter_snapshot.is_some() {
                "Restore View"
            } else {
                "Declutter"
            };
            let visible_layer_total = layer_counts.total_feature_count();
            let declutter_hover = if self.declutter_snapshot.is_some() {
                "Restore the layer visibility state captured before the last declutter action"
                    .to_string()
            } else {
                format!(
                    "Temporarily hide low-value overlays/kinds when map noise is high ({} feature glyphs currently in view)",
                    visible_layer_total
                )
            };
            if ui
                .small_button(declutter_label)
                .on_hover_text(declutter_hover)
                .clicked()
            {
                self.apply_declutter_action(&layer_counts);
            }

            let cds_active = self
                .dna_display
                .read()
                .expect("DNA display lock poisoned")
                .show_cds_features();
            let cds_count = layer_counts.kind_count("CDS");
            let cds_response = ui
                .selectable_label(cds_active, format!("CDS ({cds_count})"))
                .on_hover_text(format!(
                    "Show or hide CDS features ({cds_count} in current view)"
                ));
            if cds_response.clicked() {
                let visible = {
                    let display = self.dna_display.read().expect("DNA display lock poisoned");
                    !display.show_cds_features()
                };
                self.dna_display
                    .write()
                    .expect("DNA display lock poisoned")
                    .set_show_cds_features(visible);
                self.set_display_visibility(DisplayTarget::CdsFeatures, visible);
            }

            let gene_active = self
                .dna_display
                .read()
                .expect("DNA display lock poisoned")
                .show_gene_features();
            let gene_count = layer_counts.kind_count("GENE");
            let gene_response = ui
                .selectable_label(gene_active, format!("Gene ({gene_count})"))
                .on_hover_text(format!(
                    "Show or hide gene features ({gene_count} in current view)"
                ));
            if gene_response.clicked() {
                let visible = {
                    let display = self.dna_display.read().expect("DNA display lock poisoned");
                    !display.show_gene_features()
                };
                self.dna_display
                    .write()
                    .expect("DNA display lock poisoned")
                    .set_show_gene_features(visible);
                self.set_display_visibility(DisplayTarget::GeneFeatures, visible);
            }

            let mrna_active = self
                .dna_display
                .read()
                .expect("DNA display lock poisoned")
                .show_mrna_features();
            let mrna_count = layer_counts.kind_count("MRNA");
            let mrna_response = ui
                .selectable_label(mrna_active, format!("mRNA ({mrna_count})"))
                .on_hover_text(format!(
                    "Show or hide mRNA features ({mrna_count} in current view)"
                ));
            if mrna_response.clicked() {
                let visible = {
                    let display = self.dna_display.read().expect("DNA display lock poisoned");
                    !display.show_mrna_features()
                };
                self.dna_display
                    .write()
                    .expect("DNA display lock poisoned")
                    .set_show_mrna_features(visible);
                self.set_display_visibility(DisplayTarget::MrnaFeatures, visible);
            }

            let tfbs_active = self
                .dna_display
                .read()
                .expect("DNA display lock poisoned")
                .show_tfbs();
            let tfbs_count = layer_counts.tfbs_count();
            let response = ui
                .selectable_label(tfbs_active, format!("TFBS ({tfbs_count})"))
                .on_hover_text(format!(
                    "Show or hide computed TFBS features ({tfbs_count} in current view)"
                ));
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

            let regulatory_near_baseline = self
                .dna_display
                .read()
                .expect("DNA display lock poisoned")
                .regulatory_tracks_near_baseline();
            let regulatory_label = if regulatory_near_baseline {
                "REG@DNA"
            } else {
                "REG@TOP"
            };
            let regulatory_count = layer_counts.regulatory_feature_count;
            let regulatory_tooltip = if regulatory_near_baseline {
                "Regulatory features are placed near the DNA/GC strip. Click to move them to dedicated top lanes."
            } else {
                "Regulatory features are placed in dedicated top lanes. Click to move them near the DNA/GC strip."
            };
            let response = ui
                .button(format!("{regulatory_label} ({regulatory_count})"))
                .on_hover_text(format!(
                    "{} ({} regulatory features in current view)",
                    regulatory_tooltip, regulatory_count
                ));
            if response.clicked() {
                let new_value = !regulatory_near_baseline;
                self.dna_display
                    .write()
                    .expect("DNA display lock poisoned")
                    .set_regulatory_tracks_near_baseline(new_value);
                self.sync_regulatory_track_placement_to_engine(new_value);
            }

            for kind in self
                .feature_kinds_for_toggle_buttons()
                .into_iter()
                .filter(|kind| !matches!(kind.as_str(), "CDS" | "GENE" | "MRNA" | "TFBS"))
            {
                let visible = self
                    .dna_display
                    .read()
                    .map(|display| display.feature_kind_visible(&kind))
                    .unwrap_or(true);
                let kind_count = layer_counts.kind_count(&kind);
                let response = ui
                    .selectable_label(visible, format!("{kind} ({kind_count})"))
                    .on_hover_text(format!(
                        "Show or hide {kind} features ({kind_count} in current view)"
                    ));
                if response.clicked() {
                    if let Ok(mut display) = self.dna_display.write() {
                        display.set_feature_kind_visible(&kind, !visible);
                    }
                }
            }

            let re_active = self
                .dna_display
                .read()
                .expect("DNA display lock poisoned")
                .show_restriction_enzyme_sites();
            let re_count = layer_counts.restriction_site_count;
            let button = egui::ImageButton::new(
                ICON_RESTRICTION_ENZYMES
                    .clone()
                    .fit_to_exact_size(icon_size)
                    .rounding(5.0),
            )
            .selected(re_active);
            let response = ui
                .add(button)
                .on_hover_text(format!(
                    "Show or hide restriction enzyme cut sites ({re_count} in current view)"
                ));
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
            ui.small(format!("{re_count}"));

            let gc_active = self
                .dna_display
                .read()
                .expect("DNA display lock poisoned")
                .show_gc_contents();
            let gc_count = layer_counts.gc_region_count;
            let button = egui::ImageButton::new(
                ICON_GC_CONTENT
                    .clone()
                    .fit_to_exact_size(icon_size)
                    .rounding(5.0),
            )
            .selected(gc_active);
            let response = ui
                .add(button)
                .on_hover_text(format!(
                    "Show or hide GC-content visualization ({gc_count} GC windows in current view)"
                ));
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
            ui.small(format!("{gc_count}"));

            let (orf_active, orf_suppressed_for_anchor) = {
                let display = self.dna_display.read().expect("DNA display lock poisoned");
                (
                    display.show_open_reading_frames_effective(),
                    display.suppress_open_reading_frames_for_genome_anchor(),
                )
            };
            let orf_count = layer_counts.orf_count;
            let button = egui::ImageButton::new(
                ICON_OPEN_READING_FRAMES
                    .clone()
                    .fit_to_exact_size(icon_size)
                    .rounding(5.0),
            )
            .selected(orf_active);
            let response = ui.add_enabled(!orf_suppressed_for_anchor, button).on_hover_text(
                if orf_suppressed_for_anchor {
                    format!(
                        "Predicted ORF overlays are hidden for genome-anchored sequences to keep gene tracks readable ({orf_count} predicted ORFs in current view)"
                    )
                } else {
                    format!(
                        "Show or hide predicted open reading frames (ORFs) ({orf_count} in current view)"
                    )
                },
            );
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
            ui.small(format!("{orf_count}"));

            let methylation_active = self
                .dna_display
                .read()
                .expect("DNA display lock poisoned")
                .show_methylation_sites();
            let methylation_count = layer_counts.methylation_site_count;
            let button = egui::ImageButton::new(
                ICON_METHYLATION_SITES
                    .clone()
                    .fit_to_exact_size(icon_size)
                    .rounding(5.0),
            )
            .selected(methylation_active);
            let response = ui
                .add(button)
                .on_hover_text(format!(
                    "Show or hide methylation-site markers ({methylation_count} in current view)"
                ));
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
            ui.small(format!("{methylation_count}"));

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
            if self.is_single_stranded_rna()
                && ui
                    .button("Export RNA SVG")
                    .on_hover_text("Export RNA secondary-structure SVG via rnapkin")
                    .clicked()
            {
                self.export_rna_structure_svg();
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
            if ui
                .button("Shell")
                .on_hover_text(
                    "Open GENtle shell (shared command parser/executor with gentle_cli shell)",
                )
                .clicked()
            {
                self.show_shell = !self.show_shell;
                self.save_engine_ops_state();
            }
        });
        if !self.op_status.is_empty() {
            ui.add(egui::Label::new(egui::RichText::new(&self.op_status).monospace()).wrap());
        }

        if self.is_single_stranded_rna() {
            ui.separator();
            self.render_rna_structure_panel(ui);
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
                        egui::CollapsingHeader::new("Core cloning operations")
                            .default_open(true)
                            .show(ui, |ui| {

                ui.horizontal(|ui| {
                    ui.label(format!("Digest template {}", template_seq_id));
                    ui.label("enzymes");
                    ui.text_edit_singleline(&mut self.digest_enzymes_text);
                    ui.label("prefix");
                    ui.text_edit_singleline(&mut self.digest_prefix_text);
                    if ui
                        .button("Digest")
                        .on_hover_text("Digest active template with listed enzymes")
                        .clicked()
                    {
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
                    if ui
                        .button("Merge")
                        .on_hover_text("Merge candidates from listed input sequence IDs")
                        .clicked()
                    {
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
                    if ui
                        .button("Ligate")
                        .on_hover_text("Ligate listed inputs using selected protocol")
                        .clicked()
                    {
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
                    if ui
                        .button("Filter MW")
                        .on_hover_text("Filter inputs by molecular-weight range in bp")
                        .clicked()
                    {
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
                    ui.label("Design inputs");
                    ui.text_edit_singleline(&mut self.sq_inputs_text);
                    ui.label("prefix");
                    ui.text_edit_singleline(&mut self.sq_output_prefix);
                });
                ui.horizontal(|ui| {
                    ui.label("GC min");
                    ui.text_edit_singleline(&mut self.sq_gc_min);
                    ui.label("GC max");
                    ui.text_edit_singleline(&mut self.sq_gc_max);
                    ui.label("max homopoly");
                    ui.text_edit_singleline(&mut self.sq_max_homopolymer_run);
                });
                ui.horizontal(|ui| {
                    ui.checkbox(&mut self.sq_reject_ambiguous_bases, "Reject ambiguous");
                    ui.checkbox(&mut self.sq_avoid_u6_tttt, "Avoid U6 TTTT");
                    ui.checkbox(&mut self.sq_unique, "Unique");
                    if ui
                        .button("Filter Design")
                        .on_hover_text(
                            "Design constraints: GC bounds, homopolymers, U6 TTTT avoidance, forbidden motifs",
                        )
                        .clicked()
                    {
                        let parsed_gc_min = if self.sq_gc_min.trim().is_empty() {
                            Ok(None)
                        } else {
                            self.sq_gc_min
                                .parse::<f64>()
                                .map(Some)
                                .map_err(|_| "Invalid design GC min".to_string())
                        };
                        let parsed_gc_max = if self.sq_gc_max.trim().is_empty() {
                            Ok(None)
                        } else {
                            self.sq_gc_max
                                .parse::<f64>()
                                .map(Some)
                                .map_err(|_| "Invalid design GC max".to_string())
                        };
                        let parsed_max_homopolymer_run =
                            if self.sq_max_homopolymer_run.trim().is_empty() {
                                Ok(None)
                            } else {
                                self.sq_max_homopolymer_run
                                    .parse::<usize>()
                                    .map(Some)
                                    .map_err(|_| "Invalid design max homopolymer run".to_string())
                            };
                        match (parsed_gc_min, parsed_gc_max, parsed_max_homopolymer_run) {
                            (Ok(gc_min), Ok(gc_max), Ok(max_homopolymer_run)) => {
                                let output_prefix = if self.sq_output_prefix.trim().is_empty() {
                                    None
                                } else {
                                    Some(self.sq_output_prefix.trim().to_string())
                                };
                                self.apply_operation_with_feedback(
                                    Operation::FilterByDesignConstraints {
                                        inputs: Self::parse_ids(&self.sq_inputs_text),
                                        gc_min,
                                        gc_max,
                                        max_homopolymer_run,
                                        reject_ambiguous_bases: Some(
                                            self.sq_reject_ambiguous_bases,
                                        ),
                                        avoid_u6_terminator_tttt: Some(self.sq_avoid_u6_tttt),
                                        forbidden_motifs: Self::parse_ids(
                                            &self.sq_forbidden_motifs,
                                        ),
                                        unique: self.sq_unique,
                                        output_prefix,
                                    },
                                );
                            }
                            (Err(e), _, _) | (_, Err(e), _) | (_, _, Err(e)) => {
                                self.op_status = e;
                            }
                        }
                    }
                });
                ui.horizontal(|ui| {
                    ui.label("Forbidden motifs");
                    ui.text_edit_singleline(&mut self.sq_forbidden_motifs);
                });

                ui.horizontal(|ui| {
                    ui.label("PCR fwd");
                    ui.text_edit_singleline(&mut self.pcr_forward);
                    ui.label("rev");
                    ui.text_edit_singleline(&mut self.pcr_reverse);
                    ui.checkbox(&mut self.pcr_unique, "Unique");
                    if ui
                        .button("PCR")
                        .on_hover_text("Run PCR on active template with provided primers")
                        .clicked()
                    {
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
                    if ui
                        .button("PCR Adv")
                        .on_hover_text("Run advanced PCR with anneal/mismatch constraints")
                        .clicked()
                    {
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
                    if ui
                        .button("PCR Mut")
                        .on_hover_text("Run SNP mutagenesis PCR on active template")
                        .clicked()
                    {
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
                            });

                egui::CollapsingHeader::new("Region extraction and engine settings")
                    .default_open(true)
                    .show(ui, |ui| {
                ui.horizontal(|ui| {
                    ui.label("Extract from");
                    ui.text_edit_singleline(&mut self.extract_from);
                    ui.label("to");
                    ui.text_edit_singleline(&mut self.extract_to);
                    ui.label("output_id");
                    ui.text_edit_singleline(&mut self.extract_output_id);
                    if ui
                        .button("Extract Region")
                        .on_hover_text("Extract [from,to) region from active template sequence")
                        .clicked()
                    {
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
                    if ui
                        .button("Recompute Features")
                        .on_hover_text("Recompute annotations/features for active sequence")
                        .clicked()
                    {
                        let seq_id = self.seq_id.clone().unwrap_or_default();
                        if seq_id.is_empty() {
                            self.op_status = "No active template sequence".to_string();
                        } else {
                            self.apply_operation_with_feedback(Operation::RecomputeFeatures {
                                seq_id,
                            });
                        }
                    }
                    ui.separator();
                    let mut feature_details_font_size = self.feature_details_font_size();
                    if ui
                        .add(
                            egui::Slider::new(&mut feature_details_font_size, 8.0..=24.0)
                                .step_by(0.5)
                                .text("Feature details font")
                                .suffix(" px"),
                        )
                        .on_hover_text("Adjust feature-tree/detail text size")
                        .changed()
                    {
                        self.dna_display
                            .write()
                            .expect("DNA display lock poisoned")
                            .set_feature_details_font_size(feature_details_font_size);
                        self.sync_feature_details_font_size_to_engine(feature_details_font_size);
                    }
                    if ui
                        .button("Reset Font")
                        .on_hover_text("Reset feature detail font size to default")
                        .clicked()
                    {
                        let default_size = 11.0;
                        self.dna_display
                            .write()
                            .expect("DNA display lock poisoned")
                            .set_feature_details_font_size(default_size);
                        self.sync_feature_details_font_size_to_engine(default_size);
                    }
                });

                ui.horizontal(|ui| {
                    ui.label("SetParameter name");
                    ui.text_edit_singleline(&mut self.parameter_name);
                    ui.label("json");
                    ui.text_edit_singleline(&mut self.parameter_value_json);
                    if ui
                        .button("Set Parameter")
                        .on_hover_text("Set strict-engine runtime parameter from JSON value")
                        .clicked()
                    {
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
                    });

                egui::CollapsingHeader::new("Candidate sets (scoring/filtering)")
                    .default_open(false)
                    .show(ui, |ui| {
                        self.render_candidate_sets_ops(ui);
                    });

                egui::CollapsingHeader::new("Container-first ops, workflow, and pool export")
                    .default_open(false)
                    .show(ui, |ui| {
                ui.horizontal(|ui| {
                    ui.label("Digest container_id");
                    ui.text_edit_singleline(&mut self.container_digest_id);
                    ui.label("enzymes");
                    ui.text_edit_singleline(&mut self.digest_enzymes_text);
                    ui.label("prefix");
                    ui.text_edit_singleline(&mut self.digest_prefix_text);
                    if ui
                        .button("Digest Container")
                        .on_hover_text("Digest members of an existing container by ID")
                        .clicked()
                    {
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
                    if ui
                        .button("Merge ContainersById")
                        .on_hover_text("Merge existing containers into one new container")
                        .clicked()
                    {
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
                    if ui
                        .button("Ligate Container")
                        .on_hover_text("Ligate members of selected container")
                        .clicked()
                    {
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
                    if ui
                        .button("Filter Container MW")
                        .on_hover_text("Filter one container by molecular-weight range in bp")
                        .clicked()
                    {
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
                    if ui
                        .button("Run Workflow")
                        .on_hover_text("Execute JSON workflow operation array")
                        .clicked()
                    {
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
                    if ui
                        .button("Export Pool")
                        .on_hover_text("Export current pool candidates as JSON file")
                        .clicked()
                    {
                        self.export_pool_to_file();
                    }
                });
                ui.horizontal(|ui| {
                    let selected_preset = Self::pool_gel_ladder_preset_label(&self.pool_gel_ladders);
                    egui::ComboBox::from_id_salt("pool_gel_ladder_preset")
                        .selected_text(selected_preset.clone())
                        .show_ui(ui, |ui| {
                            for (label, value) in POOL_GEL_LADDER_PRESETS {
                                if ui
                                    .selectable_label(selected_preset == label, label)
                                    .clicked()
                                {
                                    self.pool_gel_ladders = value.to_string();
                                }
                            }
                        });
                    ui.label("gel ladders");
                    ui.text_edit_singleline(&mut self.pool_gel_ladders)
                        .on_hover_text("Optional comma-separated ladder names; leave blank for auto");
                    if ui
                        .button("Export Pool Gel SVG")
                        .on_hover_text("Export pool gel preview with auto or selected ladders")
                        .clicked()
                    {
                        self.export_pool_gel_svg();
                    }
                });
                    });

                egui::CollapsingHeader::new("Anchored region extraction (promoter-like)")
                    .default_open(false)
                    .show(ui, |ui| {
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
                    if ui
                        .button("Extract Anchored")
                        .on_hover_text("Extract anchored region candidates from active sequence")
                        .clicked()
                    {
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
                    });

                egui::CollapsingHeader::new("TFBS annotation (log-likelihood ratio)")
                    .default_open(true)
                    .show(ui, |ui| {
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
                        if ui
                            .small_button("Clear list")
                            .on_hover_text("Clear selected motif list")
                            .clicked()
                        {
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
                                    if ui
                                        .small_button("+")
                                        .on_hover_text("Add motif to selected list")
                                        .clicked()
                                    {
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
                let llr_quantile_help = "Empirical quantile of the motif score among all scanned positions/windows for this transcription-factor motif on both strands of the current sequence. Quantiles are computed per motif (not across all motifs). 0.0 disables quantile filtering; 1.0 keeps only top-scoring hits.";
                ui.horizontal(|ui| {
                    ui.label("min llr_bits");
                    ui.text_edit_singleline(&mut self.tfbs_min_llr_bits);
                    ui.label("min llr_quantile")
                        .on_hover_text(llr_quantile_help);
                    ui.text_edit_singleline(&mut self.tfbs_min_llr_quantile);
                    ui.checkbox(&mut self.tfbs_clear_existing, "Clear previous TFBS");
                });
                ui.horizontal(|ui| {
                    ui.label("per-TF min llr_bits (TF=VALUE,...)");
                    ui.text_edit_singleline(&mut self.tfbs_per_tf_min_llr_bits);
                });
                ui.horizontal(|ui| {
                    ui.label("per-TF min llr_quantile (TF=VALUE,...)")
                        .on_hover_text(llr_quantile_help);
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
                    ui.checkbox(&mut tfbs_display.use_llr_quantile, "llr_quantile")
                        .on_hover_text(llr_quantile_help);
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
                    )
                    .on_hover_text(llr_quantile_help);
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
                ui.separator();
                ui.label("VCF display filter (applies to GUI + SVG export)");
                let mut vcf_display = self
                    .dna_display
                    .read()
                    .expect("DNA display lock poisoned")
                    .vcf_display_criteria();
                let before_vcf = vcf_display.clone();
                ui.horizontal(|ui| {
                    ui.checkbox(&mut vcf_display.show_snp, "SNP");
                    ui.checkbox(&mut vcf_display.show_ins, "INS");
                    ui.checkbox(&mut vcf_display.show_del, "DEL");
                    ui.checkbox(&mut vcf_display.show_sv, "SV");
                    ui.checkbox(&mut vcf_display.show_other, "OTHER");
                });
                ui.horizontal(|ui| {
                    ui.checkbox(&mut vcf_display.pass_only, "PASS only");
                    ui.small("Hide non-PASS calls when enabled.");
                });
                ui.horizontal(|ui| {
                    ui.checkbox(&mut vcf_display.use_min_qual, "min QUAL");
                    ui.add_enabled(
                        vcf_display.use_min_qual,
                        egui::DragValue::new(&mut vcf_display.min_qual).speed(0.5),
                    );
                    ui.checkbox(&mut vcf_display.use_max_qual, "max QUAL");
                    ui.add_enabled(
                        vcf_display.use_max_qual,
                        egui::DragValue::new(&mut vcf_display.max_qual).speed(0.5),
                    );
                });
                let mut vcf_keys_changed = false;
                ui.horizontal(|ui| {
                    ui.label("required INFO keys (comma-separated)");
                    if ui
                        .text_edit_singleline(&mut self.vcf_display_required_info_keys)
                        .changed()
                    {
                        vcf_keys_changed = true;
                    }
                });
                vcf_display.required_info_keys = Self::parse_ids(&self.vcf_display_required_info_keys);
                if !vcf_display.show_snp
                    && !vcf_display.show_ins
                    && !vcf_display.show_del
                    && !vcf_display.show_sv
                    && !vcf_display.show_other
                {
                    ui.small("No VCF class enabled: all VCF overlays are hidden.");
                }
                if vcf_display != before_vcf || vcf_keys_changed {
                    self.dna_display
                        .write()
                        .expect("DNA display lock poisoned")
                        .set_vcf_display_criteria(vcf_display.clone());
                    self.sync_vcf_display_criteria_to_engine(&vcf_display);
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
                    .on_hover_text("Run TFBS annotation with current motif and threshold settings")
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
                    ui.separator();
                    self.render_pool_gel_preview(ui);
                }

                            });
                    });
                self.save_engine_ops_state();
            });
        }

        if self.show_shell {
            ui.separator();
            self.render_shell_panel(ui);
        }
    }

    fn run_shell_command(&mut self) {
        let Some(engine) = self.engine.clone() else {
            self.op_status = "No engine attached".to_string();
            return;
        };
        let line = self.shell_command_text.trim().to_string();
        if line.is_empty() {
            self.op_status = "Shell command is empty".to_string();
            return;
        }
        let command = match parse_shell_line(&line) {
            Ok(cmd) => cmd,
            Err(e) => {
                self.op_status = format!("shell parse error: {e}");
                self.shell_preview_text.clear();
                return;
            }
        };
        self.shell_preview_text = command.preview();

        let outcome = {
            let mut guard = engine.write().expect("Engine lock poisoned");
            let options = ShellExecutionOptions::from_env();
            execute_shell_command_with_options(&mut guard, &command, &options)
        };
        match outcome {
            Ok(run) => {
                let output = serde_json::to_string_pretty(&run.output).unwrap_or_else(|e| {
                    format!("{{\"error\":\"Could not format shell output: {e}\"}}")
                });
                if !self.shell_output_text.is_empty() {
                    self.shell_output_text.push('\n');
                }
                self.shell_output_text
                    .push_str(&format!("$ {line}\n{output}\n"));
                self.op_status = format!("shell ok: {}", command.preview());
                self.sync_from_engine_display();
                self.save_engine_ops_state();
            }
            Err(e) => {
                if !self.shell_output_text.is_empty() {
                    self.shell_output_text.push('\n');
                }
                self.shell_output_text
                    .push_str(&format!("$ {line}\nERROR: {e}\n"));
                self.op_status = format!("shell error: {e}");
            }
        }
    }

    fn render_shell_panel(&mut self, ui: &mut egui::Ui) {
        ui.collapsing("GENtle Shell", |ui| {
            ui.set_width(ui.available_width());
            ui.label("Runs the same shared shell parser/executor as `gentle_cli shell`.");
            ui.monospace("Type `help` to list supported commands.");
            ui.small("Command details are also documented in Help -> CLI Manual.");

            let preview = parse_shell_line(self.shell_command_text.trim())
                .map(|cmd| cmd.preview())
                .unwrap_or_else(|e| format!("Parse error: {e}"));
            self.shell_preview_text = preview;

            let command_edit = ui.add_sized(
                [ui.available_width(), 0.0],
                egui::TextEdit::singleline(&mut self.shell_command_text),
            );
            if command_edit.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter)) {
                self.run_shell_command();
            }

            ui.horizontal(|ui| {
                if ui
                    .button("Run")
                    .on_hover_text("Execute current shell command")
                    .clicked()
                {
                    self.run_shell_command();
                }
                if ui
                    .button("Shell Help")
                    .on_hover_text("Print shell command reference")
                    .clicked()
                {
                    let help = shell_help_text();
                    if !self.shell_output_text.is_empty() {
                        self.shell_output_text.push('\n');
                    }
                    self.shell_output_text
                        .push_str(&format!("$ help\n{help}\n"));
                    self.shell_preview_text = "show shell command help".to_string();
                    self.op_status = "shell help shown".to_string();
                }
                if ui
                    .button("Engine Status")
                    .on_hover_text("Run state-summary shell command")
                    .clicked()
                {
                    self.shell_command_text = "state-summary".to_string();
                    self.run_shell_command();
                }
                if ui
                    .button("Clear Output")
                    .on_hover_text("Clear shell output pane")
                    .clicked()
                {
                    self.shell_output_text.clear();
                }
            });

            ui.monospace(format!("preview: {}", self.shell_preview_text));
            let output_height = (ui.available_height() * 0.6).max(180.0);
            egui::ScrollArea::vertical()
                .id_salt(format!(
                    "gentle_shell_output_{}",
                    self.seq_id.as_deref().unwrap_or("_global")
                ))
                .max_height(output_height)
                .auto_shrink([false, false])
                .stick_to_bottom(true)
                .show(ui, |ui| {
                    ui.add_sized(
                        [ui.available_width(), output_height],
                        egui::TextEdit::multiline(&mut self.shell_output_text)
                            .desired_width(f32::INFINITY)
                            .code_editor(),
                    );
                });
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

    fn sync_vcf_display_criteria_to_engine(&self, criteria: &VcfDisplayCriteria) {
        let Some(engine) = &self.engine else {
            return;
        };
        let mut guard = engine.write().expect("Engine lock poisoned");
        let display = &mut guard.state_mut().display;
        display.vcf_display_show_snp = criteria.show_snp;
        display.vcf_display_show_ins = criteria.show_ins;
        display.vcf_display_show_del = criteria.show_del;
        display.vcf_display_show_sv = criteria.show_sv;
        display.vcf_display_show_other = criteria.show_other;
        display.vcf_display_pass_only = criteria.pass_only;
        display.vcf_display_use_min_qual = criteria.use_min_qual;
        display.vcf_display_min_qual = criteria.min_qual;
        display.vcf_display_use_max_qual = criteria.use_max_qual;
        display.vcf_display_max_qual = criteria.max_qual;
        display.vcf_display_required_info_keys = criteria.required_info_keys.clone();
    }

    fn sync_regulatory_track_placement_to_engine(&self, near_baseline: bool) {
        let Some(engine) = &self.engine else {
            return;
        };
        let mut guard = engine.write().expect("Engine lock poisoned");
        guard.state_mut().display.regulatory_tracks_near_baseline = near_baseline;
    }

    fn sync_feature_details_font_size_to_engine(&self, value: f32) {
        let Some(engine) = &self.engine else {
            return;
        };
        let mut guard = engine.write().expect("Engine lock poisoned");
        guard.state_mut().display.feature_details_font_size = value.clamp(8.0, 24.0);
    }

    fn feature_details_font_size(&self) -> f32 {
        self.dna_display
            .read()
            .map(|display| display.feature_details_font_size())
            .unwrap_or(11.0)
    }

    fn sync_from_engine_display(&mut self) {
        let Some(engine) = &self.engine else {
            return;
        };
        let Ok(guard) = engine.try_read() else {
            return;
        };
        let settings = guard.state().display.clone();
        let suppress_orf_for_anchor = self.active_sequence_is_genome_anchored(&guard);
        drop(guard);
        self.show_sequence = settings.show_sequence_panel;
        self.show_map = settings.show_map_panel;
        let mut display = self.dna_display.write().expect("DNA display lock poisoned");
        display.set_show_features(settings.show_features);
        display.set_show_cds_features(settings.show_cds_features);
        display.set_show_gene_features(settings.show_gene_features);
        display.set_show_mrna_features(settings.show_mrna_features);
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
        display.set_vcf_display_criteria(VcfDisplayCriteria {
            show_snp: settings.vcf_display_show_snp,
            show_ins: settings.vcf_display_show_ins,
            show_del: settings.vcf_display_show_del,
            show_sv: settings.vcf_display_show_sv,
            show_other: settings.vcf_display_show_other,
            pass_only: settings.vcf_display_pass_only,
            use_min_qual: settings.vcf_display_use_min_qual,
            min_qual: settings.vcf_display_min_qual,
            use_max_qual: settings.vcf_display_use_max_qual,
            max_qual: settings.vcf_display_max_qual,
            required_info_keys: settings.vcf_display_required_info_keys.clone(),
        });
        self.vcf_display_required_info_keys = settings.vcf_display_required_info_keys.join(",");
        display.set_show_restriction_enzyme_sites(settings.show_restriction_enzymes);
        display.set_show_gc_contents(settings.show_gc_contents);
        display.set_show_open_reading_frames(settings.show_open_reading_frames);
        display.set_suppress_open_reading_frames_for_genome_anchor(suppress_orf_for_anchor);
        display.set_show_methylation_sites(settings.show_methylation_sites);
        display.set_regulatory_tracks_near_baseline(settings.regulatory_tracks_near_baseline);
        display.set_linear_viewport(settings.linear_view_start_bp, settings.linear_view_span_bp);
        display.set_feature_details_font_size(settings.feature_details_font_size);
    }

    fn current_linear_viewport(&self) -> (usize, usize, usize) {
        let sequence_length = self.dna.read().expect("DNA lock poisoned").len();
        if sequence_length == 0 {
            return (0, 0, 0);
        }
        let (start_bp, span_bp) = {
            let display = self.dna_display.read().expect("DNA display lock poisoned");
            (
                display.linear_view_start_bp(),
                display.linear_view_span_bp(),
            )
        };
        let span = if span_bp == 0 || span_bp > sequence_length {
            sequence_length
        } else {
            span_bp
        };
        let max_start = sequence_length.saturating_sub(span);
        let start = start_bp.min(max_start);
        (start, span, sequence_length)
    }

    fn set_linear_viewport(&self, start_bp: usize, span_bp: usize) {
        let sequence_length = self.dna.read().expect("DNA lock poisoned").len();
        if sequence_length == 0 {
            return;
        }
        let span = if span_bp == 0 {
            sequence_length
        } else {
            span_bp.clamp(1, sequence_length)
        };
        let max_start = sequence_length.saturating_sub(span);
        let start = start_bp.min(max_start);

        self.dna_display
            .write()
            .expect("DNA display lock poisoned")
            .set_linear_viewport(start, span);
        let Some(engine) = &self.engine else {
            return;
        };
        let mut guard = engine.write().expect("Engine lock poisoned");
        let display = &mut guard.state_mut().display;
        display.linear_view_start_bp = start;
        display.linear_view_span_bp = span;
    }

    fn zoom_linear_viewport_around(&self, center_bp: usize, zoom_in: bool) {
        let (start, span, sequence_length) = self.current_linear_viewport();
        if sequence_length == 0 {
            return;
        }
        let mut new_span = if zoom_in {
            (span / 2).max(1)
        } else {
            span.saturating_mul(2).min(sequence_length)
        };
        if new_span == 0 {
            new_span = sequence_length;
        }
        let center_bp = center_bp.min(sequence_length.saturating_sub(1));
        let half = new_span / 2;
        let proposed_start = center_bp.saturating_sub(half);
        let max_start = sequence_length.saturating_sub(new_span);
        let new_start = proposed_start.min(max_start);
        if new_start != start || new_span != span {
            self.set_linear_viewport(new_start, new_span);
        }
    }

    fn pan_linear_viewport(&self, delta_bp: isize) {
        let (start, span, sequence_length) = self.current_linear_viewport();
        if sequence_length == 0 || span >= sequence_length || delta_bp == 0 {
            return;
        }
        let max_start = sequence_length.saturating_sub(span) as isize;
        let new_start = (start as isize + delta_bp).clamp(0, max_start) as usize;
        if new_start != start {
            self.set_linear_viewport(new_start, span);
        }
    }

    fn feature_bounds(&self, feature_id: usize) -> Option<(usize, usize)> {
        let dna = self.dna.read().expect("DNA lock poisoned");
        let feature = dna.features().get(feature_id)?;
        let (from, to) = feature.location.find_bounds().ok()?;
        if from < 0 || to < 0 {
            return None;
        }
        let seq_len = dna.len();
        if seq_len == 0 {
            return None;
        }
        let start = (from as usize).min(seq_len.saturating_sub(1));
        let end = (to as usize).min(seq_len.saturating_sub(1));
        Some((start, end))
    }

    fn focus_feature(&mut self, feature_id: usize) {
        if self.focused_feature_id == Some(feature_id)
            && self.pending_feature_tree_scroll_to.is_none()
        {
            return;
        }
        self.focused_feature_id = Some(feature_id);
        self.map_dna.select_feature(Some(feature_id));
        self.pending_feature_tree_scroll_to = Some(feature_id);
        let Some((start, end)) = self.feature_bounds(feature_id) else {
            return;
        };
        let sequence_length = self.dna.read().expect("DNA lock poisoned").len();
        if sequence_length == 0 {
            return;
        }
        self.dna_display
            .write()
            .expect("DNA display lock poisoned")
            .select(Selection::new(start, end, sequence_length));
        self.map_sequence.request_scroll_to_selection();

        if !self.is_circular() {
            let center = start.saturating_add(end.saturating_sub(start) / 2);
            let (_, span, len) = self.current_linear_viewport();
            if len > 0 {
                let span = span.clamp(1, len);
                let max_start = len.saturating_sub(span);
                let new_start = center.saturating_sub(span / 2).min(max_start);
                self.set_linear_viewport(new_start, span);
            }
        }
    }

    fn clear_feature_focus(&mut self) {
        self.focused_feature_id = None;
        self.map_dna.select_feature(None);
        self.pending_feature_tree_scroll_to = None;
        self.dna_display
            .write()
            .expect("DNA display lock poisoned")
            .deselect();
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
        self.clear_feature_focus();
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
        if !result.created_seq_ids.is_empty() {
            self.last_created_seq_ids = result.created_seq_ids.clone();
            self.export_pool_inputs_text = self.last_created_seq_ids.join(", ");
        }
        if let Some(engine) = self.engine.clone() {
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
                        self.clear_feature_focus();
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
                        self.clear_feature_focus();
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
        let changed_note =
            if result.created_seq_ids.is_empty() && !result.changed_seq_ids.is_empty() {
                " (in-place update)"
            } else {
                ""
            };
        self.op_status = format!(
            "ok in {} ms\ncreated: {}\ncounts: {} created, {} changed{}\nwarnings: {}\nmessages: {}",
            elapsed_ms,
            created,
            result.created_seq_ids.len(),
            result.changed_seq_ids.len(),
            changed_note,
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
                    guard.apply_with_progress(op, move |progress| {
                        let OperationProgress::Tfbs(p) = progress else {
                            return true;
                        };
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
                        true
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
        let Some(engine) = self.engine.clone() else {
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
                self.clear_feature_focus();
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

    fn is_single_stranded_rna(&self) -> bool {
        self.dna
            .read()
            .ok()
            .and_then(|dna| dna.molecule_type().map(ToString::to_string))
            .map(|v| v.eq_ignore_ascii_case("RNA") || v.eq_ignore_ascii_case("ssRNA"))
            .unwrap_or(false)
    }

    fn refresh_rna_structure(&mut self) {
        if !self.is_single_stranded_rna() {
            self.rna_status = "RNA structure is only available for single-stranded RNA".to_string();
            return;
        }
        let Some(seq_id) = self.seq_id.clone() else {
            self.rna_status = "No active sequence id for RNA structure".to_string();
            return;
        };
        let Some(engine) = self.engine.clone() else {
            self.rna_status = "No engine attached".to_string();
            return;
        };

        let sanitized_seq_id: String = seq_id
            .chars()
            .map(|c| {
                if c.is_ascii_alphanumeric() || matches!(c, '-' | '_') {
                    c
                } else {
                    '_'
                }
            })
            .collect();
        let preview_path =
            std::env::temp_dir().join(format!("gentle_rna_preview_{sanitized_seq_id}.svg"));
        let preview_path_text = preview_path.display().to_string();

        let outcome = {
            let guard = engine.read().expect("Engine lock poisoned");
            let text_report = guard.inspect_rna_structure(&seq_id);
            let svg_report = guard.render_rna_structure_svg_to_path(&seq_id, &preview_path_text);
            (text_report, svg_report)
        };

        match outcome {
            (Ok(text_report), Ok(_svg_report)) => {
                let mut text = format!(
                    "Tool: {}\nExecutable: {}\nSequence length: {}\nCommand: rnapkin {}\n",
                    text_report.tool,
                    text_report.executable,
                    text_report.sequence_length,
                    text_report.command.join(" ")
                );
                if !text_report.stdout.trim().is_empty() {
                    text.push_str("\nstdout:\n");
                    text.push_str(&text_report.stdout);
                    if !text.ends_with('\n') {
                        text.push('\n');
                    }
                }
                if !text_report.stderr.trim().is_empty() {
                    text.push_str("\nstderr:\n");
                    text.push_str(&text_report.stderr);
                    if !text.ends_with('\n') {
                        text.push('\n');
                    }
                }
                self.rna_info_text = text;
                self.rna_preview_path = Some(preview_path.clone());
                self.rna_svg_uri = format!("file://{}", preview_path.display());
                self.rna_status = "RNA structure refreshed from rnapkin".to_string();
            }
            (Err(e), _) | (_, Err(e)) => {
                self.rna_status = format!("RNA structure error: {}", e);
                self.rna_info_text.clear();
                self.rna_svg_uri.clear();
                self.rna_preview_path = None;
            }
        }
    }

    fn export_rna_structure_svg(&mut self) {
        if !self.is_single_stranded_rna() {
            self.op_status =
                "RNA structure export is only available for single-stranded RNA".to_string();
            return;
        }
        let Some(seq_id) = self.seq_id.clone() else {
            self.op_status = "No active sequence to export".to_string();
            return;
        };
        if self.engine.is_none() {
            self.op_status = "No engine attached".to_string();
            return;
        }
        let default_name = format!("{seq_id}.rna.svg");
        let path = rfd::FileDialog::new()
            .set_file_name(&default_name)
            .add_filter("SVG", &["svg"])
            .save_file();
        let Some(path) = path else {
            self.op_status = "Export canceled".to_string();
            return;
        };
        let path = path.display().to_string();
        self.apply_operation_with_feedback(Operation::RenderRnaStructureSvg { seq_id, path });
    }

    fn render_rna_structure_panel(&mut self, ui: &mut egui::Ui) {
        if !self.is_single_stranded_rna() {
            return;
        }
        if self.rna_info_text.trim().is_empty()
            && self.rna_svg_uri.trim().is_empty()
            && self.rna_status.trim().is_empty()
        {
            self.refresh_rna_structure();
        }

        ui.collapsing("RNA Structure (rnapkin)", |ui| {
            ui.horizontal(|ui| {
                if ui
                    .button("Refresh RNA Structure")
                    .on_hover_text("Run rnapkin and refresh textual report + SVG preview")
                    .clicked()
                {
                    self.refresh_rna_structure();
                }
                if ui
                    .button("Export RNA SVG")
                    .on_hover_text("Export current RNA structure image as SVG via shared engine")
                    .clicked()
                {
                    self.export_rna_structure_svg();
                }
            });

            if !self.rna_status.trim().is_empty() {
                ui.monospace(self.rna_status.clone());
            }

            if !self.rna_svg_uri.trim().is_empty() {
                ui.add(
                    egui::Image::from_uri(self.rna_svg_uri.clone())
                        .max_width(ui.available_width())
                        .max_height(360.0)
                        .shrink_to_fit(),
                );
            }

            if !self.rna_info_text.trim().is_empty() {
                ui.label("Textual report");
                ui.add(
                    egui::TextEdit::multiline(&mut self.rna_info_text)
                        .desired_rows(10)
                        .code_editor()
                        .interactive(false),
                );
            }
        });
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

    fn export_pool_gel_svg(&mut self) {
        if self.engine.is_none() {
            self.op_status = "No engine attached".to_string();
            return;
        }
        let inputs = Self::parse_ids(&self.export_pool_inputs_text);
        if inputs.is_empty() {
            self.op_status =
                "Provide at least one sequence id in 'inputs' for pool gel export".to_string();
            return;
        }
        let path = rfd::FileDialog::new()
            .set_file_name("pool.gel.svg")
            .add_filter("SVG", &["svg"])
            .save_file();
        let Some(path) = path else {
            self.op_status = "Export canceled".to_string();
            return;
        };
        let path = path.display().to_string();
        let ladders = Self::parse_ids(&self.pool_gel_ladders);
        self.apply_operation_with_feedback(Operation::RenderPoolGelSvg {
            inputs,
            path,
            ladders: if ladders.is_empty() {
                None
            } else {
                Some(ladders)
            },
        });
    }

    fn parse_optional_usize_text(raw: &str, field_name: &str) -> Result<Option<usize>, String> {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            Ok(None)
        } else {
            trimmed
                .parse::<usize>()
                .map(Some)
                .map_err(|_| format!("Invalid {field_name}: expected an integer"))
        }
    }

    fn parse_optional_f64_text(raw: &str, field_name: &str) -> Result<Option<f64>, String> {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            Ok(None)
        } else {
            trimmed
                .parse::<f64>()
                .map(Some)
                .map_err(|_| format!("Invalid {field_name}: expected a number"))
        }
    }

    fn candidate_metric_lookup(record: &CandidateRecord, key: &str) -> Option<f64> {
        if let Some(value) = record.metrics.get(key).copied() {
            return Some(value);
        }
        let lower = key.to_ascii_lowercase();
        if lower == key {
            None
        } else {
            record.metrics.get(&lower).copied()
        }
    }

    fn compare_candidate_records(
        a: &CandidateRecord,
        b: &CandidateRecord,
        sort_key: &str,
    ) -> Ordering {
        let key = sort_key.trim().to_ascii_lowercase();
        match key.as_str() {
            "seq_id" => a
                .seq_id
                .cmp(&b.seq_id)
                .then(a.start_0based.cmp(&b.start_0based))
                .then(a.end_0based.cmp(&b.end_0based)),
            "start" | "start_0based" => a
                .start_0based
                .cmp(&b.start_0based)
                .then(a.end_0based.cmp(&b.end_0based))
                .then(a.seq_id.cmp(&b.seq_id)),
            "end" | "end_0based" => a
                .end_0based
                .cmp(&b.end_0based)
                .then(a.start_0based.cmp(&b.start_0based))
                .then(a.seq_id.cmp(&b.seq_id)),
            "length" | "length_bp" => a
                .end_0based
                .saturating_sub(a.start_0based)
                .cmp(&b.end_0based.saturating_sub(b.start_0based))
                .then(a.start_0based.cmp(&b.start_0based))
                .then(a.seq_id.cmp(&b.seq_id)),
            "sequence" => a
                .sequence
                .cmp(&b.sequence)
                .then(a.seq_id.cmp(&b.seq_id))
                .then(a.start_0based.cmp(&b.start_0based)),
            _ => {
                let a_value =
                    Self::candidate_metric_lookup(a, sort_key).unwrap_or(f64::NEG_INFINITY);
                let b_value =
                    Self::candidate_metric_lookup(b, sort_key).unwrap_or(f64::NEG_INFINITY);
                a_value
                    .partial_cmp(&b_value)
                    .unwrap_or(Ordering::Equal)
                    .then(a.seq_id.cmp(&b.seq_id))
                    .then(a.start_0based.cmp(&b.start_0based))
                    .then(a.end_0based.cmp(&b.end_0based))
            }
        }
    }

    fn sort_candidate_records(rows: &mut [CandidateRecord], sort_key: &str, descending: bool) {
        rows.sort_by(|a, b| Self::compare_candidate_records(a, b, sort_key));
        if descending {
            rows.reverse();
        }
    }

    fn render_candidate_sets_ops(&mut self, ui: &mut egui::Ui) {
        let Some(engine) = self.engine.clone() else {
            ui.label("No engine attached.");
            return;
        };
        if self.candidate_source_seq_id.trim().is_empty() {
            self.candidate_source_seq_id = self.seq_id.clone().unwrap_or_default();
        }

        let (set_summaries, sequence_ids) = {
            let guard = engine.read().expect("Engine lock poisoned");
            let mut set_summaries = guard.list_candidate_sets();
            set_summaries.sort_by(|a, b| a.name.cmp(&b.name));
            let mut sequence_ids = guard.state().sequences.keys().cloned().collect::<Vec<_>>();
            sequence_ids.sort_unstable();
            (set_summaries, sequence_ids)
        };
        let set_names: Vec<String> = set_summaries.iter().map(|s| s.name.clone()).collect();

        if !set_names.is_empty() {
            if self.candidate_selected_set.trim().is_empty()
                || !set_names.iter().any(|n| n == &self.candidate_selected_set)
            {
                self.candidate_selected_set = set_names[0].clone();
            }
        } else {
            self.candidate_selected_set.clear();
        }
        if self.candidate_filter_input_set.trim().is_empty()
            && !self.candidate_selected_set.is_empty()
        {
            self.candidate_filter_input_set = self.candidate_selected_set.clone();
        }

        ui.small(
            "Generate, score, filter, combine, inspect, export, and macro-run candidate sets.",
        );

        ui.separator();
        ui.label("Generate set");
        ui.horizontal(|ui| {
            ui.label("set");
            ui.text_edit_singleline(&mut self.candidate_set_name);
            ui.label("seq");
            ui.text_edit_singleline(&mut self.candidate_source_seq_id);
            if ui
                .small_button("Use active")
                .on_hover_text("Use currently active sequence ID as candidate source")
                .clicked()
            {
                self.candidate_source_seq_id = self.seq_id.clone().unwrap_or_default();
            }
        });
        if !sequence_ids.is_empty() {
            egui::ComboBox::from_id_salt(format!(
                "candidate_source_seq_combo_{}",
                self.seq_id.as_deref().unwrap_or("_global")
            ))
            .selected_text(if self.candidate_source_seq_id.trim().is_empty() {
                "choose source sequence"
            } else {
                self.candidate_source_seq_id.as_str()
            })
            .show_ui(ui, |ui| {
                for seq_id in &sequence_ids {
                    if ui
                        .selectable_label(self.candidate_source_seq_id == *seq_id, seq_id)
                        .clicked()
                    {
                        self.candidate_source_seq_id = seq_id.clone();
                    }
                }
            });
        }
        ui.horizontal(|ui| {
            ui.label("length");
            ui.add(egui::TextEdit::singleline(&mut self.candidate_length_bp).desired_width(64.0));
            ui.label("step");
            ui.add(egui::TextEdit::singleline(&mut self.candidate_step_bp).desired_width(64.0));
            ui.label("limit");
            ui.add(egui::TextEdit::singleline(&mut self.candidate_limit).desired_width(96.0));
            ui.label("max dist");
            ui.add(
                egui::TextEdit::singleline(&mut self.candidate_max_distance_bp).desired_width(96.0),
            );
        });
        ui.horizontal(|ui| {
            ui.label("feature kinds");
            ui.text_edit_singleline(&mut self.candidate_feature_kinds);
        });
        ui.horizontal(|ui| {
            ui.label("feature label regex");
            ui.text_edit_singleline(&mut self.candidate_feature_label_regex);
            ui.label("strand");
            egui::ComboBox::from_id_salt(format!(
                "candidate_generate_strand_relation_{}",
                self.seq_id.as_deref().unwrap_or("_global")
            ))
            .selected_text(self.candidate_feature_strand_relation.as_str())
            .show_ui(ui, |ui| {
                ui.selectable_value(
                    &mut self.candidate_feature_strand_relation,
                    CandidateFeatureStrandRelation::Any,
                    "any",
                );
                ui.selectable_value(
                    &mut self.candidate_feature_strand_relation,
                    CandidateFeatureStrandRelation::Same,
                    "same (+)",
                );
                ui.selectable_value(
                    &mut self.candidate_feature_strand_relation,
                    CandidateFeatureStrandRelation::Opposite,
                    "opposite (-)",
                );
            });
            if ui
                .button("Generate")
                .on_hover_text("Run GenerateCandidateSet")
                .clicked()
            {
                let set_name = self.candidate_set_name.trim().to_string();
                let seq_id = self.candidate_source_seq_id.trim().to_string();
                let length = self.candidate_length_bp.trim().parse::<usize>();
                let step = self.candidate_step_bp.trim().parse::<usize>();
                let limit =
                    Self::parse_optional_usize_text(&self.candidate_limit, "candidate limit");
                let max_distance = Self::parse_optional_usize_text(
                    &self.candidate_max_distance_bp,
                    "candidate max distance",
                );
                if set_name.is_empty() {
                    self.op_status = "Candidate set name cannot be empty".to_string();
                } else if seq_id.is_empty() {
                    self.op_status = "Candidate source sequence ID cannot be empty".to_string();
                } else if let (Ok(length_bp), Ok(step_bp), Ok(limit), Ok(max_distance_bp)) =
                    (length, step, limit, max_distance)
                {
                    self.apply_operation_with_feedback(Operation::GenerateCandidateSet {
                        set_name: set_name.clone(),
                        seq_id,
                        length_bp,
                        step_bp,
                        feature_kinds: Self::parse_ids(&self.candidate_feature_kinds),
                        feature_label_regex: if self.candidate_feature_label_regex.trim().is_empty()
                        {
                            None
                        } else {
                            Some(self.candidate_feature_label_regex.trim().to_string())
                        },
                        max_distance_bp,
                        feature_geometry_mode: None,
                        feature_boundary_mode: None,
                        feature_strand_relation: Some(self.candidate_feature_strand_relation),
                        limit,
                    });
                    self.candidate_selected_set = set_name.clone();
                    self.candidate_filter_input_set = set_name.clone();
                    if self.candidate_setop_left.trim().is_empty() {
                        self.candidate_setop_left = set_name;
                    }
                } else {
                    self.op_status = "Invalid candidate generation parameters".to_string();
                }
            }
        });

        ui.separator();
        ui.label("Score and transform sets");
        ui.horizontal(|ui| {
            ui.label("selected set");
            egui::ComboBox::from_id_salt(format!(
                "candidate_selected_set_combo_{}",
                self.seq_id.as_deref().unwrap_or("_global")
            ))
            .selected_text(if self.candidate_selected_set.trim().is_empty() {
                "none"
            } else {
                self.candidate_selected_set.as_str()
            })
            .show_ui(ui, |ui| {
                for name in &set_names {
                    if ui
                        .selectable_label(self.candidate_selected_set == *name, name)
                        .clicked()
                    {
                        self.candidate_selected_set = name.clone();
                    }
                }
            });
            if ui
                .button("Delete selected")
                .on_hover_text("Run DeleteCandidateSet for selected set")
                .clicked()
            {
                if self.candidate_selected_set.trim().is_empty() {
                    self.op_status = "Select a candidate set to delete".to_string();
                } else {
                    self.apply_operation_with_feedback(Operation::DeleteCandidateSet {
                        set_name: self.candidate_selected_set.trim().to_string(),
                    });
                }
            }
        });

        ui.horizontal(|ui| {
            ui.label("metric");
            ui.text_edit_singleline(&mut self.candidate_score_metric);
            ui.label("expr");
            ui.text_edit_singleline(&mut self.candidate_score_expression);
            if ui
                .button("Score expr")
                .on_hover_text("Run ScoreCandidateSetExpression on selected set")
                .clicked()
            {
                if self.candidate_selected_set.trim().is_empty() {
                    self.op_status = "Select a candidate set first".to_string();
                } else if self.candidate_score_metric.trim().is_empty() {
                    self.op_status = "Metric name cannot be empty".to_string();
                } else if self.candidate_score_expression.trim().is_empty() {
                    self.op_status = "Expression cannot be empty".to_string();
                } else {
                    self.apply_operation_with_feedback(Operation::ScoreCandidateSetExpression {
                        set_name: self.candidate_selected_set.trim().to_string(),
                        metric: self.candidate_score_metric.trim().to_string(),
                        expression: self.candidate_score_expression.trim().to_string(),
                    });
                }
            }
        });

        ui.horizontal(|ui| {
            ui.label("distance metric");
            ui.text_edit_singleline(&mut self.candidate_distance_metric);
            ui.label("feature kinds");
            ui.text_edit_singleline(&mut self.candidate_distance_feature_kinds);
        });
        ui.horizontal(|ui| {
            ui.label("feature label regex");
            ui.text_edit_singleline(&mut self.candidate_distance_feature_label_regex);
            ui.label("strand");
            egui::ComboBox::from_id_salt(format!(
                "candidate_distance_strand_relation_{}",
                self.seq_id.as_deref().unwrap_or("_global")
            ))
            .selected_text(self.candidate_distance_feature_strand_relation.as_str())
            .show_ui(ui, |ui| {
                ui.selectable_value(
                    &mut self.candidate_distance_feature_strand_relation,
                    CandidateFeatureStrandRelation::Any,
                    "any",
                );
                ui.selectable_value(
                    &mut self.candidate_distance_feature_strand_relation,
                    CandidateFeatureStrandRelation::Same,
                    "same (+)",
                );
                ui.selectable_value(
                    &mut self.candidate_distance_feature_strand_relation,
                    CandidateFeatureStrandRelation::Opposite,
                    "opposite (-)",
                );
            });
            if ui
                .button("Score distance")
                .on_hover_text("Run ScoreCandidateSetDistance on selected set")
                .clicked()
            {
                if self.candidate_selected_set.trim().is_empty() {
                    self.op_status = "Select a candidate set first".to_string();
                } else if self.candidate_distance_metric.trim().is_empty() {
                    self.op_status = "Distance metric name cannot be empty".to_string();
                } else {
                    self.apply_operation_with_feedback(Operation::ScoreCandidateSetDistance {
                        set_name: self.candidate_selected_set.trim().to_string(),
                        metric: self.candidate_distance_metric.trim().to_string(),
                        feature_kinds: Self::parse_ids(&self.candidate_distance_feature_kinds),
                        feature_label_regex: if self
                            .candidate_distance_feature_label_regex
                            .trim()
                            .is_empty()
                        {
                            None
                        } else {
                            Some(
                                self.candidate_distance_feature_label_regex
                                    .trim()
                                    .to_string(),
                            )
                        },
                        feature_geometry_mode: None,
                        feature_boundary_mode: None,
                        feature_strand_relation: Some(
                            self.candidate_distance_feature_strand_relation,
                        ),
                    });
                }
            }
        });

        ui.horizontal(|ui| {
            ui.label("filter input");
            ui.text_edit_singleline(&mut self.candidate_filter_input_set);
            ui.label("output");
            ui.text_edit_singleline(&mut self.candidate_filter_output_set);
            ui.label("metric");
            ui.text_edit_singleline(&mut self.candidate_filter_metric);
        });
        ui.horizontal(|ui| {
            ui.label("min");
            ui.add(egui::TextEdit::singleline(&mut self.candidate_filter_min).desired_width(72.0));
            ui.label("max");
            ui.add(egui::TextEdit::singleline(&mut self.candidate_filter_max).desired_width(72.0));
            ui.label("min q");
            ui.add(
                egui::TextEdit::singleline(&mut self.candidate_filter_min_quantile)
                    .desired_width(72.0),
            );
            ui.label("max q");
            ui.add(
                egui::TextEdit::singleline(&mut self.candidate_filter_max_quantile)
                    .desired_width(72.0),
            );
            if ui
                .button("Filter")
                .on_hover_text("Run FilterCandidateSet")
                .clicked()
            {
                let min = Self::parse_optional_f64_text(&self.candidate_filter_min, "filter min");
                let max = Self::parse_optional_f64_text(&self.candidate_filter_max, "filter max");
                let min_q = Self::parse_optional_f64_text(
                    &self.candidate_filter_min_quantile,
                    "filter min quantile",
                );
                let max_q = Self::parse_optional_f64_text(
                    &self.candidate_filter_max_quantile,
                    "filter max quantile",
                );
                if self.candidate_filter_input_set.trim().is_empty() {
                    self.op_status = "Filter input set cannot be empty".to_string();
                } else if self.candidate_filter_output_set.trim().is_empty() {
                    self.op_status = "Filter output set cannot be empty".to_string();
                } else if self.candidate_filter_metric.trim().is_empty() {
                    self.op_status = "Filter metric cannot be empty".to_string();
                } else if let (Ok(min), Ok(max), Ok(min_quantile), Ok(max_quantile)) =
                    (min, max, min_q, max_q)
                {
                    self.apply_operation_with_feedback(Operation::FilterCandidateSet {
                        input_set: self.candidate_filter_input_set.trim().to_string(),
                        output_set: self.candidate_filter_output_set.trim().to_string(),
                        metric: self.candidate_filter_metric.trim().to_string(),
                        min,
                        max,
                        min_quantile,
                        max_quantile,
                    });
                    self.candidate_selected_set =
                        self.candidate_filter_output_set.trim().to_string();
                } else {
                    self.op_status = "Invalid numeric filter threshold".to_string();
                }
            }
        });

        ui.horizontal(|ui| {
            ui.label("set-op");
            egui::ComboBox::from_id_salt(format!(
                "candidate_setop_mode_combo_{}",
                self.seq_id.as_deref().unwrap_or("_global")
            ))
            .selected_text(if self.candidate_setop_mode.trim().is_empty() {
                "intersect"
            } else {
                self.candidate_setop_mode.as_str()
            })
            .show_ui(ui, |ui| {
                for mode in ["union", "intersect", "subtract"] {
                    if ui
                        .selectable_label(self.candidate_setop_mode == mode, mode)
                        .clicked()
                    {
                        self.candidate_setop_mode = mode.to_string();
                    }
                }
            });
            ui.label("left");
            ui.text_edit_singleline(&mut self.candidate_setop_left);
            ui.label("right");
            ui.text_edit_singleline(&mut self.candidate_setop_right);
            ui.label("output");
            ui.text_edit_singleline(&mut self.candidate_setop_output);
            if ui
                .button("Apply set-op")
                .on_hover_text("Run CandidateSetOp")
                .clicked()
            {
                let op = match self
                    .candidate_setop_mode
                    .trim()
                    .to_ascii_lowercase()
                    .as_str()
                {
                    "union" => Some(CandidateSetOperator::Union),
                    "intersect" => Some(CandidateSetOperator::Intersect),
                    "subtract" => Some(CandidateSetOperator::Subtract),
                    _ => None,
                };
                if let Some(op) = op {
                    if self.candidate_setop_left.trim().is_empty()
                        || self.candidate_setop_right.trim().is_empty()
                        || self.candidate_setop_output.trim().is_empty()
                    {
                        self.op_status =
                            "Set operation requires left, right, and output set names".to_string();
                    } else {
                        self.apply_operation_with_feedback(Operation::CandidateSetOp {
                            op,
                            left_set: self.candidate_setop_left.trim().to_string(),
                            right_set: self.candidate_setop_right.trim().to_string(),
                            output_set: self.candidate_setop_output.trim().to_string(),
                        });
                        self.candidate_selected_set =
                            self.candidate_setop_output.trim().to_string();
                    }
                } else {
                    self.op_status =
                        "Set operation mode must be union/intersect/subtract".to_string();
                }
            }
        });

        ui.separator();
        ui.label("Inspect and export");
        ui.horizontal(|ui| {
            ui.label("page limit");
            ui.add(egui::TextEdit::singleline(&mut self.candidate_page_limit).desired_width(72.0));
            ui.label("offset");
            ui.add(egui::TextEdit::singleline(&mut self.candidate_page_offset).desired_width(96.0));
            ui.label("sort");
            ui.text_edit_singleline(&mut self.candidate_sort_key);
            ui.checkbox(&mut self.candidate_sort_desc, "descending");
        });

        if !self.candidate_selected_set.trim().is_empty() {
            let requested_limit = self
                .candidate_page_limit
                .trim()
                .parse::<usize>()
                .ok()
                .filter(|v| *v > 0)
                .unwrap_or(50);
            let requested_offset = self
                .candidate_page_offset
                .trim()
                .parse::<usize>()
                .ok()
                .unwrap_or(0);

            let page_result = {
                let guard = engine.read().expect("Engine lock poisoned");
                guard.inspect_candidate_set_page(
                    &self.candidate_selected_set,
                    requested_limit,
                    requested_offset,
                )
            };

            match page_result {
                Ok((mut page, total, clamped_offset)) => {
                    if self.candidate_page_offset.trim() != clamped_offset.to_string() {
                        self.candidate_page_offset = clamped_offset.to_string();
                    }
                    Self::sort_candidate_records(
                        &mut page.candidates,
                        &self.candidate_sort_key,
                        self.candidate_sort_desc,
                    );

                    let end_exclusive = clamped_offset.saturating_add(page.candidates.len());
                    ui.horizontal(|ui| {
                        if ui
                            .button("Prev")
                            .on_hover_text("Show previous candidate page")
                            .clicked()
                        {
                            let new_offset = clamped_offset.saturating_sub(requested_limit);
                            self.candidate_page_offset = new_offset.to_string();
                        }
                        if ui
                            .button("Next")
                            .on_hover_text("Show next candidate page")
                            .clicked()
                        {
                            let next = clamped_offset.saturating_add(requested_limit);
                            self.candidate_page_offset = next.min(total).to_string();
                        }
                        ui.monospace(format!(
                            "set '{}' rows {}..{} of {}",
                            self.candidate_selected_set,
                            if total == 0 {
                                0
                            } else {
                                clamped_offset.saturating_add(1)
                            },
                            end_exclusive,
                            total
                        ));
                    });

                    let mut metric_names = set_summaries
                        .iter()
                        .find(|s| s.name == self.candidate_selected_set)
                        .map(|s| s.metrics.clone())
                        .unwrap_or_default();
                    metric_names.sort();
                    if !metric_names.is_empty() {
                        ui.small(format!("metrics: {}", metric_names.join(", ")));
                    }

                    egui::ScrollArea::vertical()
                        .max_height(180.0)
                        .show(ui, |ui| {
                            egui::Grid::new(format!(
                                "candidate_rows_grid_{}",
                                self.seq_id.as_deref().unwrap_or("_global")
                            ))
                            .striped(true)
                            .show(ui, |ui| {
                                ui.strong("seq_id");
                                ui.strong("start");
                                ui.strong("end");
                                ui.strong("len");
                                ui.strong("sort value");
                                ui.end_row();
                                for row in &page.candidates {
                                    let length_bp = row.end_0based.saturating_sub(row.start_0based);
                                    let sort_value = match self
                                        .candidate_sort_key
                                        .trim()
                                        .to_ascii_lowercase()
                                        .as_str()
                                    {
                                        "seq_id" => row.seq_id.clone(),
                                        "sequence" => row.sequence.clone(),
                                        "start" | "start_0based" => row.start_0based.to_string(),
                                        "end" | "end_0based" => row.end_0based.to_string(),
                                        "length" | "length_bp" => length_bp.to_string(),
                                        _ => Self::candidate_metric_lookup(
                                            row,
                                            self.candidate_sort_key.trim(),
                                        )
                                        .map(|v| format!("{v:.6}"))
                                        .unwrap_or_else(|| "-".to_string()),
                                    };
                                    ui.monospace(&row.seq_id);
                                    ui.monospace(row.start_0based.to_string());
                                    ui.monospace(row.end_0based.to_string());
                                    ui.monospace(length_bp.to_string());
                                    ui.monospace(sort_value);
                                    ui.end_row();
                                }
                            });
                        });

                    if ui
                        .button("Export selected set as JSON")
                        .on_hover_text("Export all rows from selected set to a JSON file")
                        .clicked()
                    {
                        let total_rows = set_summaries
                            .iter()
                            .find(|s| s.name == self.candidate_selected_set)
                            .map(|s| s.candidate_count)
                            .unwrap_or(total);
                        let export_result = {
                            let guard = engine.read().expect("Engine lock poisoned");
                            guard.inspect_candidate_set_page(
                                &self.candidate_selected_set,
                                total_rows.max(1),
                                0,
                            )
                        };
                        match export_result {
                            Ok((set, _, _)) => {
                                let path = rfd::FileDialog::new()
                                    .set_file_name(&format!(
                                        "{}.candidate_set.json",
                                        self.candidate_selected_set.replace(' ', "_")
                                    ))
                                    .add_filter("JSON", &["json"])
                                    .save_file();
                                if let Some(path) = path {
                                    let payload = serde_json::json!({
                                        "schema": "gentle.candidate_set.export.v1",
                                        "exported_at_unix_ms": std::time::SystemTime::now()
                                            .duration_since(std::time::UNIX_EPOCH)
                                            .map(|d| d.as_millis())
                                            .unwrap_or(0),
                                        "set": set
                                    });
                                    match serde_json::to_string_pretty(&payload) {
                                        Ok(text) => {
                                            if let Err(e) = fs::write(&path, text) {
                                                self.op_status = format!(
                                                    "Could not export candidate set to '{}': {e}",
                                                    path.display()
                                                );
                                            } else {
                                                self.op_status = format!(
                                                    "Exported candidate set '{}' to '{}'",
                                                    self.candidate_selected_set,
                                                    path.display()
                                                );
                                            }
                                        }
                                        Err(e) => {
                                            self.op_status = format!(
                                                "Could not serialize candidate set export: {e}"
                                            );
                                        }
                                    }
                                }
                            }
                            Err(e) => {
                                self.op_status =
                                    format!("Could not load selected set for export: {e}");
                            }
                        }
                    }
                }
                Err(e) => {
                    ui.colored_label(
                        egui::Color32::LIGHT_RED,
                        format!("inspect failed: {}", e.message),
                    );
                }
            }
        } else {
            ui.small("No candidate sets available yet.");
        }

        ui.separator();
        ui.label("Macro");
        ui.add(
            egui::TextEdit::multiline(&mut self.candidate_macro_script)
                .desired_rows(3)
                .hint_text(
                    "generate setA seq_1 --length 20 --step 1; score setA gc_score gc_fraction",
                ),
        );
        if ui
            .button("Run candidates macro")
            .on_hover_text("Execute candidates macro statements (same behavior as CLI/shell)")
            .clicked()
        {
            let script = self.candidate_macro_script.trim().to_string();
            if script.is_empty() {
                self.op_status = "Candidates macro script is empty".to_string();
            } else {
                let command = ShellCommand::CandidatesMacro {
                    script: script.clone(),
                    transactional: false,
                };
                let outcome = {
                    let mut guard = engine.write().expect("Engine lock poisoned");
                    let options = ShellExecutionOptions::from_env();
                    execute_shell_command_with_options(&mut guard, &command, &options)
                };
                match outcome {
                    Ok(run) => {
                        let output =
                            serde_json::to_string_pretty(&run.output).unwrap_or_else(|e| {
                                format!("{{\"error\":\"Could not format macro output: {e}\"}}")
                            });
                        if !self.shell_output_text.is_empty() {
                            self.shell_output_text.push('\n');
                        }
                        self.shell_output_text
                            .push_str(&format!("$ candidates macro {script}\n{output}\n"));
                        self.op_status = format!(
                            "Candidates macro executed ({} statement(s))",
                            run.output
                                .get("executed")
                                .and_then(|v| v.as_u64())
                                .unwrap_or(0)
                        );
                        self.sync_from_engine_display();
                    }
                    Err(e) => {
                        self.op_status = format!("Candidates macro failed: {e}");
                    }
                }
            }
        }
    }

    fn parse_ids(text: &str) -> Vec<String> {
        text.split(',')
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .map(|s| s.to_string())
            .collect()
    }

    fn normalized_ids(text: &str) -> Vec<String> {
        let mut ids = Self::parse_ids(text)
            .into_iter()
            .map(|id| id.to_ascii_lowercase())
            .collect::<Vec<_>>();
        ids.sort_unstable();
        ids.dedup();
        ids
    }

    fn pool_gel_ladder_preset_label(text: &str) -> String {
        let current = Self::normalized_ids(text);
        for (label, preset) in POOL_GEL_LADDER_PRESETS {
            if current == Self::normalized_ids(preset) {
                return label.to_string();
            }
        }
        if current.is_empty() {
            "Auto".to_string()
        } else {
            "Custom".to_string()
        }
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
        let (tfbs_display, vcf_display) = {
            let display = self.dna_display.read().expect("DNA display lock poisoned");
            (
                display.tfbs_display_criteria(),
                display.vcf_display_criteria(),
            )
        };
        EngineOpsUiState {
            show_engine_ops: self.show_engine_ops,
            show_shell: self.show_shell,
            shell_command_text: self.shell_command_text.clone(),
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
            sq_inputs_text: self.sq_inputs_text.clone(),
            sq_gc_min: self.sq_gc_min.clone(),
            sq_gc_max: self.sq_gc_max.clone(),
            sq_max_homopolymer_run: self.sq_max_homopolymer_run.clone(),
            sq_reject_ambiguous_bases: self.sq_reject_ambiguous_bases,
            sq_avoid_u6_tttt: self.sq_avoid_u6_tttt,
            sq_forbidden_motifs: self.sq_forbidden_motifs.clone(),
            sq_unique: self.sq_unique,
            sq_output_prefix: self.sq_output_prefix.clone(),
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
            candidate_set_name: self.candidate_set_name.clone(),
            candidate_source_seq_id: self.candidate_source_seq_id.clone(),
            candidate_length_bp: self.candidate_length_bp.clone(),
            candidate_step_bp: self.candidate_step_bp.clone(),
            candidate_feature_kinds: self.candidate_feature_kinds.clone(),
            candidate_feature_label_regex: self.candidate_feature_label_regex.clone(),
            candidate_feature_strand_relation: self.candidate_feature_strand_relation,
            candidate_max_distance_bp: self.candidate_max_distance_bp.clone(),
            candidate_limit: self.candidate_limit.clone(),
            candidate_selected_set: self.candidate_selected_set.clone(),
            candidate_page_limit: self.candidate_page_limit.clone(),
            candidate_page_offset: self.candidate_page_offset.clone(),
            candidate_sort_key: self.candidate_sort_key.clone(),
            candidate_sort_desc: self.candidate_sort_desc,
            candidate_score_metric: self.candidate_score_metric.clone(),
            candidate_score_expression: self.candidate_score_expression.clone(),
            candidate_distance_metric: self.candidate_distance_metric.clone(),
            candidate_distance_feature_kinds: self.candidate_distance_feature_kinds.clone(),
            candidate_distance_feature_label_regex: self
                .candidate_distance_feature_label_regex
                .clone(),
            candidate_distance_feature_strand_relation: self
                .candidate_distance_feature_strand_relation,
            candidate_filter_input_set: self.candidate_filter_input_set.clone(),
            candidate_filter_output_set: self.candidate_filter_output_set.clone(),
            candidate_filter_metric: self.candidate_filter_metric.clone(),
            candidate_filter_min: self.candidate_filter_min.clone(),
            candidate_filter_max: self.candidate_filter_max.clone(),
            candidate_filter_min_quantile: self.candidate_filter_min_quantile.clone(),
            candidate_filter_max_quantile: self.candidate_filter_max_quantile.clone(),
            candidate_setop_mode: self.candidate_setop_mode.clone(),
            candidate_setop_left: self.candidate_setop_left.clone(),
            candidate_setop_right: self.candidate_setop_right.clone(),
            candidate_setop_output: self.candidate_setop_output.clone(),
            candidate_macro_script: self.candidate_macro_script.clone(),
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
            vcf_display_show_snp: vcf_display.show_snp,
            vcf_display_show_ins: vcf_display.show_ins,
            vcf_display_show_del: vcf_display.show_del,
            vcf_display_show_sv: vcf_display.show_sv,
            vcf_display_show_other: vcf_display.show_other,
            vcf_display_pass_only: vcf_display.pass_only,
            vcf_display_use_min_qual: vcf_display.use_min_qual,
            vcf_display_min_qual: vcf_display.min_qual,
            vcf_display_use_max_qual: vcf_display.use_max_qual,
            vcf_display_max_qual: vcf_display.max_qual,
            vcf_display_required_info_keys: vcf_display.required_info_keys.join(","),
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
            pool_gel_ladders: self.pool_gel_ladders.clone(),
        }
    }

    fn apply_engine_ops_state(&mut self, s: EngineOpsUiState) {
        self.show_engine_ops = s.show_engine_ops;
        self.show_shell = s.show_shell;
        self.shell_command_text = s.shell_command_text;
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
        self.sq_inputs_text = s.sq_inputs_text;
        self.sq_gc_min = s.sq_gc_min;
        self.sq_gc_max = s.sq_gc_max;
        self.sq_max_homopolymer_run = s.sq_max_homopolymer_run;
        self.sq_reject_ambiguous_bases = s.sq_reject_ambiguous_bases;
        self.sq_avoid_u6_tttt = s.sq_avoid_u6_tttt;
        self.sq_forbidden_motifs = s.sq_forbidden_motifs;
        self.sq_unique = s.sq_unique;
        self.sq_output_prefix = s.sq_output_prefix;
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
        self.candidate_set_name = s.candidate_set_name;
        self.candidate_source_seq_id = s.candidate_source_seq_id;
        self.candidate_length_bp = s.candidate_length_bp;
        self.candidate_step_bp = s.candidate_step_bp;
        self.candidate_feature_kinds = s.candidate_feature_kinds;
        self.candidate_feature_label_regex = s.candidate_feature_label_regex;
        self.candidate_feature_strand_relation = s.candidate_feature_strand_relation;
        self.candidate_max_distance_bp = s.candidate_max_distance_bp;
        self.candidate_limit = s.candidate_limit;
        self.candidate_selected_set = s.candidate_selected_set;
        self.candidate_page_limit = s.candidate_page_limit;
        self.candidate_page_offset = s.candidate_page_offset;
        self.candidate_sort_key = s.candidate_sort_key;
        self.candidate_sort_desc = s.candidate_sort_desc;
        self.candidate_score_metric = s.candidate_score_metric;
        self.candidate_score_expression = s.candidate_score_expression;
        self.candidate_distance_metric = s.candidate_distance_metric;
        self.candidate_distance_feature_kinds = s.candidate_distance_feature_kinds;
        self.candidate_distance_feature_label_regex = s.candidate_distance_feature_label_regex;
        self.candidate_distance_feature_strand_relation =
            s.candidate_distance_feature_strand_relation;
        self.candidate_filter_input_set = s.candidate_filter_input_set;
        self.candidate_filter_output_set = s.candidate_filter_output_set;
        self.candidate_filter_metric = s.candidate_filter_metric;
        self.candidate_filter_min = s.candidate_filter_min;
        self.candidate_filter_max = s.candidate_filter_max;
        self.candidate_filter_min_quantile = s.candidate_filter_min_quantile;
        self.candidate_filter_max_quantile = s.candidate_filter_max_quantile;
        self.candidate_setop_mode = s.candidate_setop_mode;
        self.candidate_setop_left = s.candidate_setop_left;
        self.candidate_setop_right = s.candidate_setop_right;
        self.candidate_setop_output = s.candidate_setop_output;
        self.candidate_macro_script = s.candidate_macro_script;
        self.tfbs_motifs = s.tfbs_motifs;
        self.tfbs_use_all_motifs = s.tfbs_use_all_motifs;
        self.tfbs_catalog_filter = s.tfbs_catalog_filter;
        self.tfbs_min_llr_bits = s.tfbs_min_llr_bits;
        self.tfbs_min_llr_quantile = s.tfbs_min_llr_quantile;
        self.tfbs_per_tf_min_llr_bits = s.tfbs_per_tf_min_llr_bits;
        self.tfbs_per_tf_min_llr_quantile = s.tfbs_per_tf_min_llr_quantile;
        self.tfbs_clear_existing = s.tfbs_clear_existing;
        self.vcf_display_required_info_keys = s.vcf_display_required_info_keys;
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
        self.pool_gel_ladders = s.pool_gel_ladders;
        let tfbs_criteria = TfbsDisplayCriteria {
            use_llr_bits: s.tfbs_display_use_llr_bits,
            min_llr_bits: s.tfbs_display_min_llr_bits,
            use_llr_quantile: s.tfbs_display_use_llr_quantile,
            min_llr_quantile: s.tfbs_display_min_llr_quantile,
            use_true_log_odds_bits: s.tfbs_display_use_true_log_odds_bits,
            min_true_log_odds_bits: s.tfbs_display_min_true_log_odds_bits,
            use_true_log_odds_quantile: s.tfbs_display_use_true_log_odds_quantile,
            min_true_log_odds_quantile: s.tfbs_display_min_true_log_odds_quantile,
        };
        let vcf_criteria = VcfDisplayCriteria {
            show_snp: s.vcf_display_show_snp,
            show_ins: s.vcf_display_show_ins,
            show_del: s.vcf_display_show_del,
            show_sv: s.vcf_display_show_sv,
            show_other: s.vcf_display_show_other,
            pass_only: s.vcf_display_pass_only,
            use_min_qual: s.vcf_display_use_min_qual,
            min_qual: s.vcf_display_min_qual,
            use_max_qual: s.vcf_display_use_max_qual,
            max_qual: s.vcf_display_max_qual,
            required_info_keys: Self::parse_ids(&self.vcf_display_required_info_keys),
        };
        let mut display = self.dna_display.write().expect("DNA display lock poisoned");
        display.set_tfbs_display_criteria(tfbs_criteria);
        display.set_vcf_display_criteria(vcf_criteria.clone());
        drop(display);
        self.sync_tfbs_display_criteria_to_engine(tfbs_criteria);
        self.sync_vcf_display_criteria_to_engine(&vcf_criteria);
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
                if ui
                    .button("Close")
                    .on_hover_text("Dismiss this error popup")
                    .clicked()
                {
                    self.op_error_popup = None;
                }
            });
        if !open {
            self.op_error_popup = None;
        }
    }

    fn current_pool_members(&self) -> Vec<(String, usize)> {
        let Some(engine) = &self.engine else {
            return vec![];
        };
        let state = engine.read().expect("Engine lock poisoned");
        self.last_created_seq_ids
            .iter()
            .filter_map(|id| {
                state
                    .state()
                    .sequences
                    .get(id)
                    .map(|dna| (id.clone(), dna.len()))
            })
            .collect()
    }

    fn render_pool_distribution(&self, ui: &mut egui::Ui) {
        let mut lengths: Vec<usize> = self
            .current_pool_members()
            .into_iter()
            .map(|(_, bp)| bp)
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

    fn render_pool_gel_preview(&self, ui: &mut egui::Ui) {
        let members = self.current_pool_members();
        if members.is_empty() {
            return;
        }
        let ladders = Self::parse_ids(&self.pool_gel_ladders);
        let layout = match build_pool_gel_layout(&members, &ladders) {
            Ok(layout) => layout,
            Err(e) => {
                ui.colored_label(egui::Color32::from_rgb(180, 40, 40), e);
                return;
            }
        };
        let ladder_text = if layout.selected_ladders.is_empty() {
            "auto".to_string()
        } else {
            layout.selected_ladders.join(" + ")
        };
        ui.monospace(format!(
            "Pool gel: ladders={} | range={}..{} bp",
            ladder_text, layout.range_min_bp, layout.range_max_bp
        ));
        let desired_h = 320.0;
        let desired_w = ui.available_width().max(460.0);
        let (rect, _response) =
            ui.allocate_exact_size(Vec2::new(desired_w, desired_h), egui::Sense::hover());
        let painter = ui.painter_at(rect);

        painter.rect_filled(rect, 8.0, egui::Color32::from_rgb(246, 248, 251));
        let gel_rect = egui::Rect::from_min_max(
            rect.min + Vec2::new(14.0, 18.0),
            rect.max - Vec2::new(180.0, 32.0),
        );
        painter.rect_filled(gel_rect, 8.0, egui::Color32::from_rgb(18, 20, 24));

        let lane_count = layout.lanes.len().max(1);
        let lane_gap = gel_rect.width() / (lane_count as f32 + 1.0);
        let mut ticks = BTreeSet::new();
        for lane in layout.lanes.iter().filter(|l| l.is_ladder) {
            for band in &lane.bands {
                ticks.insert(band.bp);
            }
        }
        if ticks.is_empty() {
            ticks.insert(layout.range_min_bp);
            ticks.insert(layout.range_max_bp);
        }
        let mut accepted_ticks: Vec<usize> = vec![];
        let mut last_y: Option<f32> = None;
        for bp in ticks.iter().rev() {
            let y = layout.y_for_bp(*bp, gel_rect.top(), gel_rect.bottom());
            if last_y.map(|v| (v - y).abs() >= 14.0).unwrap_or(true) {
                accepted_ticks.push(*bp);
                last_y = Some(y);
            }
            if accepted_ticks.len() >= 16 {
                break;
            }
        }
        for bp in accepted_ticks {
            let y = layout.y_for_bp(bp, gel_rect.top(), gel_rect.bottom());
            painter.line_segment(
                [
                    egui::pos2(gel_rect.left(), y),
                    egui::pos2(gel_rect.right(), y),
                ],
                egui::Stroke::new(1.0, egui::Color32::from_rgb(42, 46, 53)),
            );
            painter.text(
                egui::pos2(gel_rect.right() + 8.0, y),
                egui::Align2::LEFT_CENTER,
                format!("{bp} bp"),
                egui::FontId::monospace(10.0),
                egui::Color32::from_rgb(60, 72, 90),
            );
        }

        for (lane_idx, lane) in layout.lanes.iter().enumerate() {
            let x = gel_rect.left() + lane_gap * (lane_idx as f32 + 1.0);
            let lane_rect = egui::Rect::from_center_size(
                egui::pos2(x, gel_rect.center().y),
                Vec2::new(72.0, gel_rect.height() - 8.0),
            );
            painter.rect_filled(
                lane_rect,
                6.0,
                if lane.is_ladder {
                    egui::Color32::from_rgb(27, 33, 42)
                } else {
                    egui::Color32::from_rgb(30, 38, 50)
                },
            );
            for band in &lane.bands {
                let y = layout.y_for_bp(band.bp, gel_rect.top() + 4.0, gel_rect.bottom() - 4.0);
                let width = if lane.is_ladder {
                    30.0 + 18.0 * band.intensity
                } else {
                    36.0 + 26.0 * band.intensity
                };
                let band_rect = egui::Rect::from_center_size(
                    egui::pos2(x, y),
                    Vec2::new(width, if lane.is_ladder { 3.5 } else { 4.5 }),
                );
                painter.rect_filled(
                    band_rect,
                    2.0,
                    if lane.is_ladder {
                        egui::Color32::from_rgb(232, 236, 242)
                    } else {
                        egui::Color32::from_rgb(245, 158, 11)
                    },
                );
                if !lane.is_ladder {
                    let mut label = format!("{} bp", band.bp);
                    if band.count > 1 {
                        label.push_str(&format!(" (x{})", band.count));
                    }
                    painter.text(
                        egui::pos2(gel_rect.right() + 72.0, y),
                        egui::Align2::LEFT_CENTER,
                        label,
                        egui::FontId::monospace(10.0),
                        egui::Color32::from_rgb(17, 24, 39),
                    );
                }
            }
            painter.text(
                egui::pos2(x, gel_rect.bottom() + 14.0),
                egui::Align2::CENTER_TOP,
                lane.name.clone(),
                egui::FontId::monospace(10.5),
                egui::Color32::from_rgb(15, 23, 42),
            );
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
        struct FeatureTreeEntry {
            id: usize,
            button_label: String,
            hover_text: Option<String>,
            subgroup_key: Option<String>,
            subgroup_label: Option<String>,
        }

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
        let (
            show_cds_features,
            show_gene_features,
            show_mrna_features,
            show_tfbs,
            tfbs_display_criteria,
            vcf_display_criteria,
            hidden_feature_kinds,
            feature_details_font_size,
        ) = {
            let display = self.dna_display.read().expect("DNA display lock poisoned");
            (
                display.show_cds_features(),
                display.show_gene_features(),
                display.show_mrna_features(),
                display.show_tfbs(),
                display.tfbs_display_criteria(),
                display.vcf_display_criteria(),
                display.hidden_feature_kinds().clone(),
                display.feature_details_font_size(),
            )
        };
        let typed_features = self
            .dna
            .read()
            .expect("DNA lock poisoned")
            .features()
            .iter()
            .enumerate()
            .filter_map(|(id, feature)| {
                if RenderDna::is_source_feature(feature) {
                    return None;
                }
                if !RenderDna::feature_passes_kind_filter(
                    feature,
                    show_cds_features,
                    show_gene_features,
                    show_mrna_features,
                ) {
                    return None;
                }
                if hidden_feature_kinds.contains(&feature.kind.to_string().to_ascii_uppercase()) {
                    return None;
                }
                if RenderDna::is_tfbs_feature(feature) {
                    if !show_tfbs {
                        return None;
                    }
                    if !RenderDna::tfbs_feature_passes_display_filter(
                        feature,
                        tfbs_display_criteria,
                    ) {
                        return None;
                    }
                }
                if RenderDna::is_vcf_track_feature(feature)
                    && !RenderDna::vcf_feature_passes_display_filter(feature, &vcf_display_criteria)
                {
                    return None;
                }
                let (from, to) = feature.location.find_bounds().ok()?;
                if from < 0 || to < 0 {
                    return None;
                }
                let feature_label = {
                    let name = RenderDna::feature_name(feature);
                    if name.trim().is_empty() {
                        format!("{} #{}", feature.kind.to_string(), id + 1)
                    } else {
                        name
                    }
                };
                let subgroup_label = if RenderDna::is_track_feature(feature) {
                    RenderDna::track_group_label(feature)
                } else if RenderDna::is_tfbs_feature(feature) {
                    RenderDna::tfbs_group_label(feature)
                } else if RenderDna::is_restriction_site_feature(feature) {
                    RenderDna::restriction_site_group_label(feature)
                } else {
                    None
                };
                let subgroup_key = subgroup_label
                    .as_ref()
                    .map(|label| label.trim().to_ascii_uppercase());
                let kind_label = if RenderDna::is_track_feature(feature) {
                    "Tracks".to_string()
                } else {
                    feature.kind.to_string()
                };
                let range_label = RenderDna::feature_range_text(feature);
                let grouped_entry = subgroup_label.is_some();
                let button_label = if grouped_entry && !range_label.is_empty() {
                    range_label.clone()
                } else {
                    feature_label.clone()
                };
                let hover_text = if grouped_entry {
                    if range_label.is_empty() {
                        Some(feature_label.clone())
                    } else {
                        Some(format!("{feature_label} ({range_label})"))
                    }
                } else {
                    None
                };
                Some((
                    kind_label,
                    FeatureTreeEntry {
                        id,
                        button_label,
                        hover_text,
                        subgroup_key,
                        subgroup_label,
                    },
                ))
            })
            .collect::<Vec<_>>();
        let mut grouped_features: HashMap<String, Vec<FeatureTreeEntry>> = HashMap::new();
        let seq_key = self.seq_id.as_deref().unwrap_or("_global").to_string();
        for (kind, entry) in typed_features {
            grouped_features.entry(kind).or_default().push(entry);
        }
        let mut group_keys = grouped_features.keys().collect::<Vec<_>>();
        group_keys.sort();
        let mut clicked_feature: Option<usize> = None;
        let feature_font_size = feature_details_font_size;
        let kind_font_size = feature_font_size + 1.0;
        for kind in group_keys {
            let Some(entries) = grouped_features.get(kind) else {
                continue;
            };
            let has_selected = selected_id
                .map(|selected| entries.iter().any(|entry| entry.id == selected))
                .unwrap_or(false);
            egui::CollapsingHeader::new(
                egui::RichText::new(kind.as_str())
                    .size(kind_font_size)
                    .strong(),
            )
            .id_salt(format!("feature_kind_{seq_key}_{kind}"))
            .open(if has_selected { Some(true) } else { None })
            .show(ui, |ui| {
                let mut grouped_entries: HashMap<String, (String, Vec<&FeatureTreeEntry>)> =
                    HashMap::new();
                let mut ungrouped_entries: Vec<&FeatureTreeEntry> = Vec::new();
                for entry in entries {
                    if let (Some(subgroup_key), Some(subgroup_label)) =
                        (&entry.subgroup_key, &entry.subgroup_label)
                    {
                        grouped_entries
                            .entry(subgroup_key.clone())
                            .or_insert_with(|| (subgroup_label.clone(), Vec::new()))
                            .1
                            .push(entry);
                    } else {
                        ungrouped_entries.push(entry);
                    }
                }

                let mut render_entry = |ui: &mut egui::Ui, entry: &FeatureTreeEntry| {
                    let selected = selected_id == Some(entry.id);
                    ui.horizontal(|ui| {
                        let button = egui::Button::new(
                            egui::RichText::new(&entry.button_label).size(feature_font_size),
                        )
                        .selected(selected);
                        let mut response = ui.add(button);
                        if let Some(hover_text) = &entry.hover_text {
                            response = response.on_hover_text(hover_text);
                        }
                        if selected && self.pending_feature_tree_scroll_to == Some(entry.id) {
                            response.scroll_to_me(Some(egui::Align::Center));
                            self.pending_feature_tree_scroll_to = None;
                        }
                        if response.clicked() {
                            clicked_feature = Some(entry.id);
                        }
                    });
                };

                let mut subgroup_keys = grouped_entries.keys().cloned().collect::<Vec<_>>();
                subgroup_keys.sort_by(|left, right| {
                    let left_label = grouped_entries
                        .get(left)
                        .map(|(label, _)| label.as_str())
                        .unwrap_or("");
                    let right_label = grouped_entries
                        .get(right)
                        .map(|(label, _)| label.as_str())
                        .unwrap_or("");
                    left_label.cmp(right_label).then_with(|| left.cmp(right))
                });

                for subgroup_key in subgroup_keys {
                    let Some((subgroup_label, subgroup_entries)) =
                        grouped_entries.get(&subgroup_key)
                    else {
                        continue;
                    };
                    let subgroup_has_selected = selected_id
                        .map(|selected| subgroup_entries.iter().any(|entry| entry.id == selected))
                        .unwrap_or(false);
                    let subgroup_heading = format!("{subgroup_label} ({})", subgroup_entries.len());
                    egui::CollapsingHeader::new(
                        egui::RichText::new(subgroup_heading)
                            .size(feature_font_size)
                            .strong(),
                    )
                    .id_salt(format!(
                        "feature_kind_group_{seq_key}_{kind}_{}",
                        subgroup_key
                    ))
                    .open(if subgroup_has_selected {
                        Some(true)
                    } else {
                        None
                    })
                    .show(ui, |ui| {
                        for entry in subgroup_entries {
                            render_entry(ui, entry);
                        }
                    });
                }

                for entry in ungrouped_entries {
                    render_entry(ui, entry);
                }
            });
        }
        if let Some(id) = clicked_feature {
            self.focus_feature(id);
        }
    }

    fn refresh_description_cache(&mut self) {
        let selected_id = self.get_selected_feature_id();
        let mut clear_invalid_selection = false;
        {
            let dna = self.dna.read().expect("DNA lock poisoned");
            let seq_len = dna.len();
            let feature_count = dna.features().len();
            let needs_refresh = !self.description_cache_initialized
                || self.description_cache_selected_id != selected_id
                || self.description_cache_seq_len != seq_len
                || self.description_cache_feature_count != feature_count;
            if !needs_refresh {
                return;
            }

            self.description_cache_initialized = true;
            self.description_cache_selected_id = selected_id;
            self.description_cache_seq_len = seq_len;
            self.description_cache_feature_count = feature_count;

            if let Some(id) = selected_id {
                if let Some(feature) = dna.features().get(id) {
                    let label = RenderDna::feature_name(feature);
                    self.description_cache_title = if label.trim().is_empty() {
                        feature.kind.to_string()
                    } else {
                        label
                    };
                    self.description_cache_range = Some(RenderDna::feature_range_text(feature));
                } else {
                    clear_invalid_selection = true;
                    self.description_cache_selected_id = None;
                    self.description_cache_title = dna.description().join("\n");
                    self.description_cache_range = None;
                }
            } else {
                self.description_cache_title = dna.description().join("\n");
                self.description_cache_range = None;
            }
        }
        if clear_invalid_selection {
            self.clear_feature_focus();
        }
    }

    pub fn render_description(&mut self, ui: &mut egui::Ui) {
        egui::ScrollArea::vertical().show(ui, |ui| {
            ui.set_min_height(150.0);
            self.refresh_description_cache();
            ui.heading(&self.description_cache_title);
            if let Some(range) = &self.description_cache_range {
                ui.label(
                    egui::RichText::new(range)
                        .monospace()
                        .size(self.feature_details_font_size()),
                );
            }
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
                    ui.set_width(320.0);
                    let tree_height = (ui.available_height() * 0.55).max(180.0);
                    egui::ScrollArea::vertical()
                        .max_height(tree_height)
                        .show(ui, |ui| {
                            self.render_features(ui);
                        });
                    ui.separator();
                    self.render_description(ui);
                });
            });

            ui.separator();

            let response = ui.add(self.map_dna.to_owned());

            if response.hovered() {
                let pointer_state: PointerState = ctx.input(|i| i.pointer.to_owned());
                self.map_dna.on_hover(pointer_state);

                if let Some(feature_id) = self.map_dna.get_hovered_feature_id() {
                    let hover_text = self.dna.read().ok().and_then(|dna| {
                        let feature = dna.features().get(feature_id)?;
                        let name = {
                            let label = RenderDna::feature_name(feature);
                            if label.trim().is_empty() {
                                feature.kind.to_string()
                            } else {
                                label
                            }
                        };
                        let range = {
                            let text = RenderDna::feature_range_text(feature);
                            if text.is_empty() {
                                "-".to_string()
                            } else {
                                text
                            }
                        };
                        Some((name, feature.kind.to_string(), range))
                    });
                    if let Some((name, kind, range)) = hover_text {
                        let detail_font_size = self.feature_details_font_size();
                        response.clone().on_hover_ui_at_pointer(|ui| {
                            ui.strong(name);
                            ui.label(
                                egui::RichText::new(format!("{kind} | {range}"))
                                    .monospace()
                                    .size(detail_font_size),
                            );
                        });
                    }
                }

                if !self.is_circular() {
                    let scroll = ctx.input(|i| i.raw_scroll_delta);
                    if scroll.y.abs() > 0.0 {
                        let (start_bp, span_bp, sequence_length) = self.current_linear_viewport();
                        if sequence_length > 0 {
                            let center_bp = response
                                .hover_pos()
                                .map(|pos| {
                                    let frac = ((pos.x - response.rect.left())
                                        / response.rect.width().max(1.0))
                                    .clamp(0.0, 1.0);
                                    start_bp + (frac * span_bp as f32).floor() as usize
                                })
                                .unwrap_or_else(|| start_bp.saturating_add(span_bp / 2));
                            self.zoom_linear_viewport_around(center_bp, scroll.y > 0.0);
                        }
                    } else if scroll.x.abs() > 0.0 {
                        let (_, span_bp, _) = self.current_linear_viewport();
                        let delta_bp = ((-scroll.x / response.rect.width().max(1.0))
                            * span_bp as f32) as isize;
                        self.pan_linear_viewport(delta_bp);
                    }
                }
            }

            if response.clicked() {
                let previous_selected = self.map_dna.get_selected_feature_id();
                let pointer_state: PointerState = ctx.input(|i| i.pointer.to_owned());
                self.map_dna.on_click(pointer_state);
                let next_selected = self.map_dna.get_selected_feature_id();
                if next_selected != previous_selected {
                    if let Some(feature_id) = next_selected {
                        self.focus_feature(feature_id);
                    } else {
                        self.clear_feature_focus();
                    }
                }
            }

            if response.double_clicked() {
                let previous_selected = self.map_dna.get_selected_feature_id();
                let pointer_state: PointerState = ctx.input(|i| i.pointer.to_owned());
                self.map_dna.on_double_click(pointer_state);
                let next_selected = self.map_dna.get_selected_feature_id();
                if next_selected != previous_selected {
                    if let Some(feature_id) = next_selected {
                        self.focus_feature(feature_id);
                    } else {
                        self.clear_feature_focus();
                    }
                }
            }
        });
    }
}
