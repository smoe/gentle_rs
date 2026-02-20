use crate::{
    app::GENtleApp,
    dna_sequence::DNAsequence,
    enzymes::active_restriction_enzymes,
    genomes::{
        GenomeBlastReport, GenomeCatalog, GenomeGeneRecord, GenomeSourcePlan,
        PrepareGenomeProgress, PrepareGenomeReport, PreparedGenomeInspection,
        DEFAULT_GENOME_CATALOG_PATH, DEFAULT_HELPER_GENOME_CATALOG_PATH,
    },
    iupac_code::IupacCode,
    lineage_export::export_lineage_svg,
    methylation_sites::MethylationMode,
    pool_gel::{build_pool_gel_layout, export_pool_gel_svg},
    render_export::{export_circular_svg, export_linear_svg},
    restriction_enzyme::RestrictionEnzyme,
    rna_structure::{self, RnaStructureError, RnaStructureSvgReport, RnaStructureTextReport},
    tf_motifs, DNA_LADDERS, RNA_LADDERS,
};
use flate2::read::GzDecoder;
use regex::{Regex, RegexBuilder};
use serde::{Deserialize, Serialize};
use serde_json::json;
use std::{
    cmp::Ordering,
    collections::hash_map::DefaultHasher,
    collections::{BTreeSet, HashMap, HashSet},
    env,
    error::Error,
    fmt,
    fs::File,
    hash::{Hash, Hasher},
    io::{BufRead, BufReader, Write},
    path::Path,
    process::Command,
    time::Instant,
};
use tempfile::NamedTempFile;

pub type SeqId = String;
pub type OpId = String;
pub type RunId = String;
pub type NodeId = String;
pub type ContainerId = String;
const PROVENANCE_METADATA_KEY: &str = "provenance";
const GENOME_EXTRACTIONS_METADATA_KEY: &str = "genome_extractions";
pub const GENOME_TRACK_SUBSCRIPTIONS_METADATA_KEY: &str = "genome_bed_track_subscriptions";
const GENOME_TRACK_KNOWN_ANCHORS_METADATA_KEY: &str = "genome_track_autosync_known_anchors";
pub const CANDIDATE_SETS_METADATA_KEY: &str = "candidate_sets";
const CANDIDATE_SETS_SCHEMA: &str = "gentle.candidate_sets.v1";
const GENOME_BED_TRACK_GENERATED_TAG: &str = "genome_bed_track";
const GENOME_BIGWIG_TRACK_GENERATED_TAG: &str = "genome_bigwig_track";
const GENOME_VCF_TRACK_GENERATED_TAG: &str = "genome_vcf_track";
const BLAST_HIT_TRACK_GENERATED_TAG: &str = "blast_hit_track";
pub const DEFAULT_BIGWIG_TO_BEDGRAPH_BIN: &str = "bigWigToBedGraph";
pub const BIGWIG_TO_BEDGRAPH_ENV_BIN: &str = "GENTLE_BIGWIG_TO_BEDGRAPH_BIN";
const MAX_IMPORTED_SIGNAL_FEATURES: usize = 25_000;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum DisplayTarget {
    SequencePanel,
    MapPanel,
    Features,
    CdsFeatures,
    GeneFeatures,
    MrnaFeatures,
    Tfbs,
    RestrictionEnzymes,
    GcContents,
    OpenReadingFrames,
    MethylationSites,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct DisplaySettings {
    pub show_sequence_panel: bool,
    pub show_map_panel: bool,
    pub show_features: bool,
    pub show_cds_features: bool,
    pub show_gene_features: bool,
    pub show_mrna_features: bool,
    pub show_tfbs: bool,
    pub regulatory_tracks_near_baseline: bool,
    pub tfbs_display_use_llr_bits: bool,
    pub tfbs_display_min_llr_bits: f64,
    pub tfbs_display_use_llr_quantile: bool,
    pub tfbs_display_min_llr_quantile: f64,
    pub tfbs_display_use_true_log_odds_bits: bool,
    pub tfbs_display_min_true_log_odds_bits: f64,
    pub tfbs_display_use_true_log_odds_quantile: bool,
    pub tfbs_display_min_true_log_odds_quantile: f64,
    pub show_restriction_enzymes: bool,
    pub show_gc_contents: bool,
    pub show_open_reading_frames: bool,
    pub show_methylation_sites: bool,
    pub linear_view_start_bp: usize,
    pub linear_view_span_bp: usize,
    pub feature_details_font_size: f32,
}

impl Default for DisplaySettings {
    fn default() -> Self {
        Self {
            show_sequence_panel: true,
            show_map_panel: true,
            show_features: true,
            show_cds_features: true,
            show_gene_features: true,
            show_mrna_features: true,
            show_tfbs: false,
            regulatory_tracks_near_baseline: false,
            tfbs_display_use_llr_bits: true,
            tfbs_display_min_llr_bits: 0.0,
            tfbs_display_use_llr_quantile: true,
            tfbs_display_min_llr_quantile: 0.95,
            tfbs_display_use_true_log_odds_bits: false,
            tfbs_display_min_true_log_odds_bits: 0.0,
            tfbs_display_use_true_log_odds_quantile: false,
            tfbs_display_min_true_log_odds_quantile: 0.95,
            show_restriction_enzymes: true,
            show_gc_contents: true,
            show_open_reading_frames: false,
            show_methylation_sites: false,
            linear_view_start_bp: 0,
            linear_view_span_bp: 0,
            feature_details_font_size: 11.0,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum SequenceOrigin {
    ImportedGenomic,
    ImportedCdna,
    ImportedSynthetic,
    ImportedUnknown,
    Derived,
    InSilicoSelection,
    Branch,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LineageNode {
    pub node_id: NodeId,
    pub seq_id: SeqId,
    pub created_by_op: Option<OpId>,
    pub origin: SequenceOrigin,
    pub created_at_unix_ms: u128,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LineageEdge {
    pub from_node_id: NodeId,
    pub to_node_id: NodeId,
    pub op_id: OpId,
    pub run_id: RunId,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct LineageGraph {
    pub nodes: HashMap<NodeId, LineageNode>,
    pub seq_to_node: HashMap<SeqId, NodeId>,
    pub edges: Vec<LineageEdge>,
    pub next_node_counter: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ContainerKind {
    Singleton,
    Pool,
    Selection,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Container {
    pub container_id: ContainerId,
    pub kind: ContainerKind,
    pub name: Option<String>,
    pub members: Vec<SeqId>,
    pub created_by_op: Option<OpId>,
    pub created_at_unix_ms: u128,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ContainerState {
    pub containers: HashMap<ContainerId, Container>,
    pub seq_to_latest_container: HashMap<SeqId, ContainerId>,
    pub next_container_counter: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct EngineParameters {
    pub max_fragments_per_container: usize,
}

impl Default for EngineParameters {
    fn default() -> Self {
        Self {
            max_fragments_per_container: 80_000,
        }
    }
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ProjectState {
    pub sequences: HashMap<SeqId, DNAsequence>,
    pub metadata: HashMap<String, serde_json::Value>,
    #[serde(default)]
    pub display: DisplaySettings,
    #[serde(default)]
    pub lineage: LineageGraph,
    #[serde(default)]
    pub parameters: EngineParameters,
    #[serde(default)]
    pub container_state: ContainerState,
}

impl ProjectState {
    pub fn load_from_path(path: &str) -> Result<Self, EngineError> {
        let text = std::fs::read_to_string(path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not read state file '{path}': {e}"),
        })?;
        serde_json::from_str(&text).map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!("Could not parse state JSON '{path}': {e}"),
        })
    }

    pub fn save_to_path(&self, path: &str) -> Result<(), EngineError> {
        let text = serde_json::to_string_pretty(self).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize state: {e}"),
        })?;
        std::fs::write(path, text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write state file '{path}': {e}"),
        })
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ExportFormat {
    GenBank,
    Fasta,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum RenderSvgMode {
    Linear,
    Circular,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum PrimerLibraryMode {
    Enumerate,
    Sample,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PcrPrimerSpec {
    pub sequence: String,
    pub anneal_len: Option<usize>,
    pub max_mismatches: Option<usize>,
    pub require_3prime_exact_bases: Option<usize>,
    pub library_mode: Option<PrimerLibraryMode>,
    pub max_variants: Option<usize>,
    pub sample_seed: Option<u64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SnpMutationSpec {
    pub zero_based_position: usize,
    pub reference: String,
    pub alternate: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TfThresholdOverride {
    pub tf: String,
    pub min_llr_bits: Option<f64>,
    pub min_llr_quantile: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum LigationProtocol {
    Sticky,
    Blunt,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum AnchorBoundary {
    Start,
    End,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum AnchorDirection {
    Upstream,
    Downstream,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum AnchoredRegionAnchor {
    Position {
        zero_based: usize,
    },
    FeatureBoundary {
        feature_kind: Option<String>,
        feature_label: Option<String>,
        boundary: AnchorBoundary,
        occurrence: Option<usize>,
    },
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum GenomeTrackSource {
    #[default]
    Bed,
    BigWig,
    Vcf,
}

impl GenomeTrackSource {
    pub fn from_path(path: &str) -> Self {
        let lower = path.trim().to_ascii_lowercase();
        if lower.ends_with(".bw") || lower.ends_with(".bigwig") {
            Self::BigWig
        } else if lower.ends_with(".vcf") || lower.ends_with(".vcf.gz") {
            Self::Vcf
        } else {
            Self::Bed
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::Bed => "BED",
            Self::BigWig => "BigWig",
            Self::Vcf => "VCF",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default)]
pub struct GenomeTrackSubscription {
    pub source: GenomeTrackSource,
    pub path: String,
    pub track_name: Option<String>,
    pub min_score: Option<f64>,
    pub max_score: Option<f64>,
    pub clear_existing: bool,
}

impl Default for GenomeTrackSubscription {
    fn default() -> Self {
        Self {
            source: GenomeTrackSource::Bed,
            path: String::new(),
            track_name: None,
            min_score: None,
            max_score: None,
            clear_existing: false,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct GenomeTrackSyncReport {
    pub subscriptions_considered: usize,
    pub target_sequences: usize,
    pub applied_imports: usize,
    pub failed_imports: usize,
    pub warnings_count: usize,
    pub errors: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
struct CandidateStore {
    schema: String,
    updated_at_unix_ms: u128,
    sets: HashMap<String, CandidateSet>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct CandidateSet {
    pub name: String,
    pub created_at_unix_ms: u128,
    pub source_seq_ids: Vec<String>,
    pub candidates: Vec<CandidateRecord>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct CandidateRecord {
    pub seq_id: String,
    pub start_0based: usize,
    pub end_0based: usize,
    pub sequence: String,
    pub metrics: HashMap<String, f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CandidateSetSummary {
    pub name: String,
    pub created_at_unix_ms: u128,
    pub source_seq_ids: Vec<String>,
    pub candidate_count: usize,
    pub metrics: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CandidateMetricSummary {
    pub metric: String,
    pub present_in_candidates: usize,
    pub missing_in_candidates: usize,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CandidateSetOperator {
    Union,
    Intersect,
    Subtract,
}

impl CandidateSetOperator {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Union => "union",
            Self::Intersect => "intersect",
            Self::Subtract => "subtract",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlastHitFeatureInput {
    pub subject_id: String,
    pub query_start_1based: usize,
    pub query_end_1based: usize,
    pub subject_start_1based: usize,
    pub subject_end_1based: usize,
    pub identity_percent: f64,
    pub bit_score: f64,
    pub evalue: f64,
    pub query_coverage_percent: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Operation {
    LoadFile {
        path: String,
        as_id: Option<SeqId>,
    },
    SaveFile {
        seq_id: SeqId,
        path: String,
        format: ExportFormat,
    },
    RenderSequenceSvg {
        seq_id: SeqId,
        mode: RenderSvgMode,
        path: String,
    },
    RenderRnaStructureSvg {
        seq_id: SeqId,
        path: String,
    },
    RenderLineageSvg {
        path: String,
    },
    RenderPoolGelSvg {
        inputs: Vec<SeqId>,
        path: String,
        ladders: Option<Vec<String>>,
    },
    ExportDnaLadders {
        path: String,
        #[serde(default)]
        name_filter: Option<String>,
    },
    ExportRnaLadders {
        path: String,
        #[serde(default)]
        name_filter: Option<String>,
    },
    ExportPool {
        inputs: Vec<SeqId>,
        path: String,
        pool_id: Option<String>,
        human_id: Option<String>,
    },
    PrepareGenome {
        genome_id: String,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
    },
    ExtractGenomeRegion {
        genome_id: String,
        chromosome: String,
        start_1based: usize,
        end_1based: usize,
        output_id: Option<SeqId>,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
    },
    ExtractGenomeGene {
        genome_id: String,
        gene_query: String,
        occurrence: Option<usize>,
        output_id: Option<SeqId>,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
    },
    ImportGenomeBedTrack {
        seq_id: SeqId,
        path: String,
        track_name: Option<String>,
        min_score: Option<f64>,
        max_score: Option<f64>,
        clear_existing: Option<bool>,
    },
    ImportGenomeBigWigTrack {
        seq_id: SeqId,
        path: String,
        track_name: Option<String>,
        min_score: Option<f64>,
        max_score: Option<f64>,
        clear_existing: Option<bool>,
    },
    ImportGenomeVcfTrack {
        seq_id: SeqId,
        path: String,
        track_name: Option<String>,
        min_score: Option<f64>,
        max_score: Option<f64>,
        clear_existing: Option<bool>,
    },
    ImportBlastHitsTrack {
        seq_id: SeqId,
        hits: Vec<BlastHitFeatureInput>,
        track_name: Option<String>,
        clear_existing: Option<bool>,
    },
    DigestContainer {
        container_id: ContainerId,
        enzymes: Vec<String>,
        output_prefix: Option<String>,
    },
    MergeContainersById {
        container_ids: Vec<ContainerId>,
        output_prefix: Option<String>,
    },
    LigationContainer {
        container_id: ContainerId,
        circularize_if_possible: bool,
        output_id: Option<SeqId>,
        protocol: LigationProtocol,
        output_prefix: Option<String>,
        unique: Option<bool>,
    },
    FilterContainerByMolecularWeight {
        container_id: ContainerId,
        min_bp: usize,
        max_bp: usize,
        error: f64,
        unique: bool,
        output_prefix: Option<String>,
    },
    Digest {
        input: SeqId,
        enzymes: Vec<String>,
        output_prefix: Option<String>,
    },
    Ligation {
        inputs: Vec<SeqId>,
        circularize_if_possible: bool,
        output_id: Option<SeqId>,
        protocol: LigationProtocol,
        output_prefix: Option<String>,
        unique: Option<bool>,
    },
    MergeContainers {
        inputs: Vec<SeqId>,
        output_prefix: Option<String>,
    },
    Pcr {
        template: SeqId,
        forward_primer: String,
        reverse_primer: String,
        output_id: Option<SeqId>,
        unique: Option<bool>,
    },
    PcrAdvanced {
        template: SeqId,
        forward_primer: PcrPrimerSpec,
        reverse_primer: PcrPrimerSpec,
        output_id: Option<SeqId>,
        unique: Option<bool>,
    },
    PcrMutagenesis {
        template: SeqId,
        forward_primer: PcrPrimerSpec,
        reverse_primer: PcrPrimerSpec,
        mutations: Vec<SnpMutationSpec>,
        output_id: Option<SeqId>,
        unique: Option<bool>,
        require_all_mutations: Option<bool>,
    },
    ExtractRegion {
        input: SeqId,
        from: usize,
        to: usize,
        output_id: Option<SeqId>,
    },
    ExtractAnchoredRegion {
        input: SeqId,
        anchor: AnchoredRegionAnchor,
        direction: AnchorDirection,
        target_length_bp: usize,
        length_tolerance_bp: usize,
        required_re_sites: Vec<String>,
        required_tf_motifs: Vec<String>,
        forward_primer: Option<String>,
        reverse_primer: Option<String>,
        output_prefix: Option<String>,
        unique: Option<bool>,
        max_candidates: Option<usize>,
    },
    SelectCandidate {
        input: SeqId,
        criterion: String,
        output_id: Option<SeqId>,
    },
    FilterByMolecularWeight {
        inputs: Vec<SeqId>,
        min_bp: usize,
        max_bp: usize,
        error: f64,
        unique: bool,
        output_prefix: Option<String>,
    },
    #[serde(alias = "FilterBySequenceQuality")]
    FilterByDesignConstraints {
        inputs: Vec<SeqId>,
        gc_min: Option<f64>,
        gc_max: Option<f64>,
        max_homopolymer_run: Option<usize>,
        reject_ambiguous_bases: Option<bool>,
        avoid_u6_terminator_tttt: Option<bool>,
        #[serde(default)]
        forbidden_motifs: Vec<String>,
        unique: bool,
        output_prefix: Option<String>,
    },
    GenerateCandidateSet {
        set_name: String,
        seq_id: SeqId,
        length_bp: usize,
        step_bp: usize,
        #[serde(default)]
        feature_kinds: Vec<String>,
        feature_label_regex: Option<String>,
        max_distance_bp: Option<usize>,
        limit: Option<usize>,
    },
    DeleteCandidateSet {
        set_name: String,
    },
    ScoreCandidateSetExpression {
        set_name: String,
        metric: String,
        expression: String,
    },
    ScoreCandidateSetDistance {
        set_name: String,
        metric: String,
        #[serde(default)]
        feature_kinds: Vec<String>,
        feature_label_regex: Option<String>,
    },
    FilterCandidateSet {
        input_set: String,
        output_set: String,
        metric: String,
        min: Option<f64>,
        max: Option<f64>,
        min_quantile: Option<f64>,
        max_quantile: Option<f64>,
    },
    CandidateSetOp {
        op: CandidateSetOperator,
        left_set: String,
        right_set: String,
        output_set: String,
    },
    Reverse {
        input: SeqId,
        output_id: Option<SeqId>,
    },
    Complement {
        input: SeqId,
        output_id: Option<SeqId>,
    },
    ReverseComplement {
        input: SeqId,
        output_id: Option<SeqId>,
    },
    Branch {
        input: SeqId,
        output_id: Option<SeqId>,
    },
    SetDisplayVisibility {
        target: DisplayTarget,
        visible: bool,
    },
    SetLinearViewport {
        start_bp: usize,
        span_bp: usize,
    },
    SetTopology {
        seq_id: SeqId,
        circular: bool,
    },
    RecomputeFeatures {
        seq_id: SeqId,
    },
    SetParameter {
        name: String,
        value: serde_json::Value,
    },
    AnnotateTfbs {
        seq_id: SeqId,
        motifs: Vec<String>,
        min_llr_bits: Option<f64>,
        min_llr_quantile: Option<f64>,
        #[serde(default)]
        per_tf_thresholds: Vec<TfThresholdOverride>,
        clear_existing: Option<bool>,
        #[serde(default)]
        max_hits: Option<usize>,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct PoolEnd {
    end_type: String,
    forward_5: String,
    forward_3: String,
    reverse_5: String,
    reverse_3: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct PoolMember {
    seq_id: SeqId,
    human_id: String,
    name: Option<String>,
    sequence: String,
    length_bp: usize,
    topology: String,
    ends: PoolEnd,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct PoolExport {
    schema: String,
    pool_id: String,
    human_id: String,
    member_count: usize,
    members: Vec<PoolMember>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DnaLadderBandInfo {
    pub length_bp: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub relative_strength: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DnaLadderInfo {
    pub name: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub loading_hint: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub min_bp: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub max_bp: Option<usize>,
    pub band_count: usize,
    pub bands: Vec<DnaLadderBandInfo>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DnaLadderCatalog {
    pub schema: String,
    pub ladder_count: usize,
    pub ladders: Vec<DnaLadderInfo>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DnaLadderExportReport {
    pub path: String,
    pub ladder_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RnaLadderBandInfo {
    pub length_nt: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub relative_strength: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RnaLadderInfo {
    pub name: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub loading_hint: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub min_nt: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub max_nt: Option<usize>,
    pub band_count: usize,
    pub bands: Vec<RnaLadderBandInfo>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RnaLadderCatalog {
    pub schema: String,
    pub ladder_count: usize,
    pub ladders: Vec<RnaLadderInfo>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RnaLadderExportReport {
    pub path: String,
    pub ladder_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GenomeExtractionProvenance {
    pub seq_id: SeqId,
    pub recorded_at_unix_ms: u128,
    pub operation: String,
    pub genome_id: String,
    pub catalog_path: String,
    pub cache_dir: Option<String>,
    pub chromosome: Option<String>,
    pub start_1based: Option<usize>,
    pub end_1based: Option<usize>,
    pub gene_query: Option<String>,
    pub occurrence: Option<usize>,
    pub gene_id: Option<String>,
    pub gene_name: Option<String>,
    pub strand: Option<char>,
    pub anchor_strand: Option<char>,
    pub sequence_source_type: Option<String>,
    pub annotation_source_type: Option<String>,
    pub sequence_source: Option<String>,
    pub annotation_source: Option<String>,
    pub sequence_sha1: Option<String>,
    pub annotation_sha1: Option<String>,
}

#[derive(Debug, Clone)]
struct GenomeSequenceAnchor {
    genome_id: String,
    chromosome: String,
    start_1based: usize,
    end_1based: usize,
    strand: Option<char>,
}

#[derive(Debug, Clone)]
struct BedRecord {
    chromosome: String,
    start_0based: usize,
    end_0based: usize,
    name: Option<String>,
    score: Option<f64>,
    strand: Option<char>,
}

#[derive(Debug, Clone)]
struct VcfRecord {
    chromosome: String,
    pos_1based: usize,
    id: Option<String>,
    reference: String,
    alternates: Vec<String>,
    qual: Option<f64>,
    filter: Option<String>,
    info: Option<String>,
}

#[derive(Debug, Clone, Default)]
struct GenomeBedTrackImportReport {
    track_name: String,
    parsed_records: usize,
    imported_features: usize,
    skipped_records: usize,
    skipped_invalid: usize,
    skipped_wrong_chromosome: usize,
    skipped_non_overlap: usize,
    skipped_missing_score: usize,
    skipped_outside_score_range: usize,
    truncated_at_limit: bool,
    warnings: Vec<String>,
}

#[derive(Debug, Clone)]
struct FeatureInterval {
    kind_upper: String,
    labels_upper: Vec<String>,
    start_0based: usize,
    end_0based: usize,
}

#[derive(Debug, Clone, PartialEq)]
enum ExpressionToken {
    Number(f64),
    Ident(String),
    Plus,
    Minus,
    Star,
    Slash,
    LParen,
    RParen,
    Comma,
}

#[derive(Debug, Clone)]
enum MetricExpr {
    Number(f64),
    Variable(String),
    UnaryMinus(Box<MetricExpr>),
    Binary {
        op: ExpressionBinaryOp,
        left: Box<MetricExpr>,
        right: Box<MetricExpr>,
    },
    Function {
        name: String,
        args: Vec<MetricExpr>,
    },
}

#[derive(Debug, Clone, Copy)]
enum ExpressionBinaryOp {
    Add,
    Subtract,
    Multiply,
    Divide,
}

struct MetricExpressionParser {
    tokens: Vec<ExpressionToken>,
    index: usize,
}

impl MetricExpressionParser {
    fn new(tokens: Vec<ExpressionToken>) -> Self {
        Self { tokens, index: 0 }
    }

    fn parse(mut self) -> Result<MetricExpr, String> {
        let expr = self.parse_add_sub()?;
        if self.index != self.tokens.len() {
            return Err("Unexpected trailing expression tokens".to_string());
        }
        Ok(expr)
    }

    fn peek(&self) -> Option<&ExpressionToken> {
        self.tokens.get(self.index)
    }

    fn consume(&mut self) -> Option<ExpressionToken> {
        let token = self.tokens.get(self.index).cloned();
        if token.is_some() {
            self.index += 1;
        }
        token
    }

    fn parse_add_sub(&mut self) -> Result<MetricExpr, String> {
        let mut expr = self.parse_mul_div()?;
        loop {
            let op = match self.peek() {
                Some(ExpressionToken::Plus) => ExpressionBinaryOp::Add,
                Some(ExpressionToken::Minus) => ExpressionBinaryOp::Subtract,
                _ => break,
            };
            let _ = self.consume();
            let rhs = self.parse_mul_div()?;
            expr = MetricExpr::Binary {
                op,
                left: Box::new(expr),
                right: Box::new(rhs),
            };
        }
        Ok(expr)
    }

    fn parse_mul_div(&mut self) -> Result<MetricExpr, String> {
        let mut expr = self.parse_unary()?;
        loop {
            let op = match self.peek() {
                Some(ExpressionToken::Star) => ExpressionBinaryOp::Multiply,
                Some(ExpressionToken::Slash) => ExpressionBinaryOp::Divide,
                _ => break,
            };
            let _ = self.consume();
            let rhs = self.parse_unary()?;
            expr = MetricExpr::Binary {
                op,
                left: Box::new(expr),
                right: Box::new(rhs),
            };
        }
        Ok(expr)
    }

    fn parse_unary(&mut self) -> Result<MetricExpr, String> {
        if matches!(self.peek(), Some(ExpressionToken::Minus)) {
            let _ = self.consume();
            return Ok(MetricExpr::UnaryMinus(Box::new(self.parse_unary()?)));
        }
        self.parse_primary()
    }

    fn parse_primary(&mut self) -> Result<MetricExpr, String> {
        match self.consume() {
            Some(ExpressionToken::Number(value)) => Ok(MetricExpr::Number(value)),
            Some(ExpressionToken::Ident(name)) => {
                if matches!(self.peek(), Some(ExpressionToken::LParen)) {
                    let _ = self.consume();
                    let mut args = vec![];
                    if !matches!(self.peek(), Some(ExpressionToken::RParen)) {
                        loop {
                            args.push(self.parse_add_sub()?);
                            if matches!(self.peek(), Some(ExpressionToken::Comma)) {
                                let _ = self.consume();
                                continue;
                            }
                            break;
                        }
                    }
                    if !matches!(self.consume(), Some(ExpressionToken::RParen)) {
                        return Err(format!("Function '{}' is missing closing ')'", name));
                    }
                    Ok(MetricExpr::Function { name, args })
                } else {
                    Ok(MetricExpr::Variable(name))
                }
            }
            Some(ExpressionToken::LParen) => {
                let expr = self.parse_add_sub()?;
                if !matches!(self.consume(), Some(ExpressionToken::RParen)) {
                    return Err("Missing ')' in expression".to_string());
                }
                Ok(expr)
            }
            _ => Err("Unexpected token while parsing expression".to_string()),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Workflow {
    pub run_id: RunId,
    pub ops: Vec<Operation>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OpResult {
    pub op_id: OpId,
    pub created_seq_ids: Vec<SeqId>,
    pub changed_seq_ids: Vec<SeqId>,
    pub warnings: Vec<String>,
    pub messages: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TfbsProgress {
    pub seq_id: String,
    pub motif_id: String,
    pub motif_index: usize,
    pub motif_count: usize,
    pub scanned_steps: usize,
    pub total_steps: usize,
    pub motif_percent: f64,
    pub total_percent: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum OperationProgress {
    Tfbs(TfbsProgress),
    GenomePrepare(PrepareGenomeProgress),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OperationRecord {
    pub run_id: RunId,
    pub op: Operation,
    pub result: OpResult,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum ErrorCode {
    InvalidInput,
    NotFound,
    Unsupported,
    Io,
    Internal,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EngineError {
    pub code: ErrorCode,
    pub message: String,
}

impl fmt::Display for EngineError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}: {}", self.code, self.message)
    }
}

impl Error for EngineError {}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Capabilities {
    pub protocol_version: String,
    pub supported_operations: Vec<String>,
    pub supported_export_formats: Vec<String>,
    pub deterministic_operation_log: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EngineSequenceSummary {
    pub id: String,
    pub name: Option<String>,
    pub length: usize,
    pub circular: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EngineContainerSummary {
    pub id: String,
    pub kind: String,
    pub member_count: usize,
    pub members: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EngineStateSummary {
    pub sequence_count: usize,
    pub sequences: Vec<EngineSequenceSummary>,
    pub container_count: usize,
    pub containers: Vec<EngineContainerSummary>,
    pub display: DisplaySettings,
}

pub trait Engine {
    fn apply(&mut self, op: Operation) -> Result<OpResult, EngineError>;
    fn apply_workflow(&mut self, wf: Workflow) -> Result<Vec<OpResult>, EngineError>;
    fn snapshot(&self) -> &ProjectState;
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct GentleEngine {
    state: ProjectState,
    journal: Vec<OperationRecord>,
    op_counter: u64,
}

impl GentleEngine {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn from_state(state: ProjectState) -> Self {
        let mut ret = Self {
            state,
            ..Self::default()
        };
        ret.reconcile_lineage_nodes();
        ret.reconcile_containers();
        ret
    }

    pub fn state(&self) -> &ProjectState {
        &self.state
    }

    pub fn state_mut(&mut self) -> &mut ProjectState {
        &mut self.state
    }

    pub fn list_sequences_with_genome_anchor(&self) -> Vec<String> {
        let mut seq_ids: Vec<String> = self.state.sequences.keys().cloned().collect();
        seq_ids.sort_unstable();
        seq_ids
            .into_iter()
            .filter(|seq_id| self.latest_genome_anchor_for_seq(seq_id).is_ok())
            .collect()
    }

    pub fn describe_sequence_genome_anchor(&self, seq_id: &str) -> Result<String, EngineError> {
        let anchor = self.latest_genome_anchor_for_seq(seq_id)?;
        let strand = anchor.strand.unwrap_or('+');
        Ok(format!(
            "{}:{}-{} ({}, strand {})",
            anchor.chromosome, anchor.start_1based, anchor.end_1based, anchor.genome_id, strand
        ))
    }

    pub fn list_genome_track_subscriptions(&self) -> Vec<GenomeTrackSubscription> {
        Self::read_track_subscriptions_from_metadata(
            self.state
                .metadata
                .get(GENOME_TRACK_SUBSCRIPTIONS_METADATA_KEY),
        )
    }

    pub fn add_genome_track_subscription(
        &mut self,
        subscription: GenomeTrackSubscription,
    ) -> Result<bool, EngineError> {
        let normalized = Self::normalize_genome_track_subscription(subscription)?;
        let mut subscriptions = self.list_genome_track_subscriptions();
        if subscriptions.iter().any(|existing| existing == &normalized) {
            return Ok(false);
        }
        subscriptions.push(normalized);
        Self::sort_track_subscriptions(&mut subscriptions);
        self.write_track_subscriptions_to_metadata(&subscriptions)?;
        Ok(true)
    }

    pub fn remove_genome_track_subscription(
        &mut self,
        index: usize,
    ) -> Result<GenomeTrackSubscription, EngineError> {
        let mut subscriptions = self.list_genome_track_subscriptions();
        if index >= subscriptions.len() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Track subscription index {} is out of bounds (len={})",
                    index,
                    subscriptions.len()
                ),
            });
        }
        let removed = subscriptions.remove(index);
        self.write_track_subscriptions_to_metadata(&subscriptions)?;
        Ok(removed)
    }

    pub fn clear_genome_track_subscriptions(&mut self) {
        self.state
            .metadata
            .remove(GENOME_TRACK_SUBSCRIPTIONS_METADATA_KEY);
        self.state
            .metadata
            .remove(GENOME_TRACK_KNOWN_ANCHORS_METADATA_KEY);
    }

    pub fn import_genome_track_to_all_anchored(
        &mut self,
        subscription: GenomeTrackSubscription,
        track_subscription: bool,
    ) -> Result<GenomeTrackSyncReport, EngineError> {
        let normalized = Self::normalize_genome_track_subscription(subscription)?;
        let seq_ids = self.list_sequences_with_genome_anchor();
        if seq_ids.is_empty() {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message:
                    "No genome-anchored sequence is available. Extract a genome region or gene first."
                        .to_string(),
            });
        }
        let mut report = self.apply_track_subscription_to_seq_ids(&seq_ids, &normalized);
        report.subscriptions_considered = 1;
        report.target_sequences = seq_ids.len();
        if track_subscription {
            let _ = self.add_genome_track_subscription(normalized)?;
            let known: BTreeSet<String> = seq_ids.iter().cloned().collect();
            self.write_known_track_anchor_ids(&known);
        }
        Ok(report)
    }

    pub fn apply_tracked_genome_track_subscription(
        &mut self,
        index: usize,
    ) -> Result<GenomeTrackSyncReport, EngineError> {
        let subscriptions = self.list_genome_track_subscriptions();
        let Some(subscription) = subscriptions.get(index).cloned() else {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Track subscription index {} is out of bounds (len={})",
                    index,
                    subscriptions.len()
                ),
            });
        };
        let seq_ids = self.list_sequences_with_genome_anchor();
        if seq_ids.is_empty() {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message:
                    "No genome-anchored sequence is available. Extract a genome region or gene first."
                        .to_string(),
            });
        }
        let mut report = self.apply_track_subscription_to_seq_ids(&seq_ids, &subscription);
        report.subscriptions_considered = 1;
        report.target_sequences = seq_ids.len();
        let known: BTreeSet<String> = seq_ids.iter().cloned().collect();
        self.write_known_track_anchor_ids(&known);
        Ok(report)
    }

    pub fn sync_tracked_genome_track_subscriptions(
        &mut self,
        only_new_anchors: bool,
    ) -> Result<GenomeTrackSyncReport, EngineError> {
        let subscriptions = self.list_genome_track_subscriptions();
        let current_anchors = self.list_sequences_with_genome_anchor();
        let current_set: BTreeSet<String> = current_anchors.iter().cloned().collect();
        let target_seq_ids: Vec<String> = if only_new_anchors {
            let known = self.read_known_track_anchor_ids();
            current_set.difference(&known).cloned().collect()
        } else {
            current_anchors.clone()
        };

        if subscriptions.is_empty() || target_seq_ids.is_empty() {
            self.write_known_track_anchor_ids(&current_set);
            return Ok(GenomeTrackSyncReport {
                subscriptions_considered: subscriptions.len(),
                target_sequences: target_seq_ids.len(),
                ..GenomeTrackSyncReport::default()
            });
        }

        let mut report = GenomeTrackSyncReport {
            subscriptions_considered: subscriptions.len(),
            target_sequences: target_seq_ids.len(),
            ..GenomeTrackSyncReport::default()
        };
        for seq_id in &target_seq_ids {
            for subscription in &subscriptions {
                let op = Self::track_subscription_to_operation(seq_id, subscription);
                match self.apply(op) {
                    Ok(op_result) => {
                        report.applied_imports += 1;
                        report.warnings_count += op_result.warnings.len();
                    }
                    Err(e) => {
                        report.failed_imports += 1;
                        if report.errors.len() < 20 {
                            report.errors.push(format!(
                                "{} '{}' @ {}: {}",
                                subscription.source.label(),
                                subscription.path,
                                seq_id,
                                e.message
                            ));
                        }
                    }
                }
            }
        }
        self.write_known_track_anchor_ids(&current_set);
        Ok(report)
    }

    pub fn capabilities() -> Capabilities {
        Capabilities {
            protocol_version: "v1".to_string(),
            supported_operations: vec![
                "LoadFile".to_string(),
                "SaveFile".to_string(),
                "RenderSequenceSvg".to_string(),
                "RenderRnaStructureSvg".to_string(),
                "RenderLineageSvg".to_string(),
                "RenderPoolGelSvg".to_string(),
                "ExportDnaLadders".to_string(),
                "ExportRnaLadders".to_string(),
                "ExportPool".to_string(),
                "PrepareGenome".to_string(),
                "ExtractGenomeRegion".to_string(),
                "ExtractGenomeGene".to_string(),
                "ImportGenomeBedTrack".to_string(),
                "ImportGenomeBigWigTrack".to_string(),
                "ImportGenomeVcfTrack".to_string(),
                "ImportBlastHitsTrack".to_string(),
                "DigestContainer".to_string(),
                "MergeContainersById".to_string(),
                "LigationContainer".to_string(),
                "FilterContainerByMolecularWeight".to_string(),
                "Digest".to_string(),
                "MergeContainers".to_string(),
                "Ligation".to_string(),
                "Pcr".to_string(),
                "PcrAdvanced".to_string(),
                "PcrMutagenesis".to_string(),
                "ExtractRegion".to_string(),
                "ExtractAnchoredRegion".to_string(),
                "SelectCandidate".to_string(),
                "FilterByMolecularWeight".to_string(),
                "FilterByDesignConstraints".to_string(),
                "GenerateCandidateSet".to_string(),
                "DeleteCandidateSet".to_string(),
                "ScoreCandidateSetExpression".to_string(),
                "ScoreCandidateSetDistance".to_string(),
                "FilterCandidateSet".to_string(),
                "CandidateSetOp".to_string(),
                "Reverse".to_string(),
                "Complement".to_string(),
                "ReverseComplement".to_string(),
                "Branch".to_string(),
                "SetDisplayVisibility".to_string(),
                "SetLinearViewport".to_string(),
                "SetTopology".to_string(),
                "RecomputeFeatures".to_string(),
                "SetParameter".to_string(),
                "AnnotateTfbs".to_string(),
            ],
            supported_export_formats: vec!["GenBank".to_string(), "Fasta".to_string()],
            deterministic_operation_log: true,
        }
    }

    pub fn inspect_dna_ladders(name_filter: Option<&str>) -> DnaLadderCatalog {
        let filter = name_filter
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .map(|v| v.to_ascii_lowercase());

        let mut ladders: Vec<DnaLadderInfo> = vec![];
        for name in DNA_LADDERS.names_sorted() {
            if let Some(filter_text) = &filter {
                if !name.to_ascii_lowercase().contains(filter_text) {
                    continue;
                }
            }
            let Some(ladder) = DNA_LADDERS.get(&name) else {
                continue;
            };
            let bands = ladder
                .bands()
                .iter()
                .map(|band| DnaLadderBandInfo {
                    length_bp: band.length_bp(),
                    relative_strength: band.relative_strength,
                })
                .collect::<Vec<_>>();
            ladders.push(DnaLadderInfo {
                name,
                loading_hint: ladder.loading_hint(),
                min_bp: ladder.min_bp(),
                max_bp: ladder.max_bp(),
                band_count: bands.len(),
                bands,
            });
        }

        DnaLadderCatalog {
            schema: "gentle.dna_ladders.v1".to_string(),
            ladder_count: ladders.len(),
            ladders,
        }
    }

    pub fn export_dna_ladders(
        path: &str,
        name_filter: Option<&str>,
    ) -> Result<DnaLadderExportReport, EngineError> {
        let catalog = Self::inspect_dna_ladders(name_filter);
        let text = serde_json::to_string_pretty(&catalog).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize DNA ladders JSON: {e}"),
        })?;
        std::fs::write(path, text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write DNA ladders file '{path}': {e}"),
        })?;
        Ok(DnaLadderExportReport {
            path: path.to_string(),
            ladder_count: catalog.ladder_count,
        })
    }

    pub fn inspect_rna_ladders(name_filter: Option<&str>) -> RnaLadderCatalog {
        let filter = name_filter
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .map(|v| v.to_ascii_lowercase());

        let mut ladders: Vec<RnaLadderInfo> = vec![];
        for name in RNA_LADDERS.names_sorted() {
            if let Some(filter_text) = &filter {
                if !name.to_ascii_lowercase().contains(filter_text) {
                    continue;
                }
            }
            let Some(ladder) = RNA_LADDERS.get(&name) else {
                continue;
            };
            let bands = ladder
                .bands()
                .iter()
                .map(|band| RnaLadderBandInfo {
                    length_nt: band.length_nt(),
                    relative_strength: band.relative_strength,
                })
                .collect::<Vec<_>>();
            ladders.push(RnaLadderInfo {
                name,
                loading_hint: ladder.loading_hint(),
                min_nt: ladder.min_nt(),
                max_nt: ladder.max_nt(),
                band_count: bands.len(),
                bands,
            });
        }

        RnaLadderCatalog {
            schema: "gentle.rna_ladders.v1".to_string(),
            ladder_count: ladders.len(),
            ladders,
        }
    }

    pub fn export_rna_ladders(
        path: &str,
        name_filter: Option<&str>,
    ) -> Result<RnaLadderExportReport, EngineError> {
        let catalog = Self::inspect_rna_ladders(name_filter);
        let text = serde_json::to_string_pretty(&catalog).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize RNA ladders JSON: {e}"),
        })?;
        std::fs::write(path, text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write RNA ladders file '{path}': {e}"),
        })?;
        Ok(RnaLadderExportReport {
            path: path.to_string(),
            ladder_count: catalog.ladder_count,
        })
    }

    fn map_rna_structure_error(err: RnaStructureError) -> EngineError {
        match err {
            RnaStructureError::UnsupportedBiotype { .. } | RnaStructureError::EmptySequence => {
                EngineError {
                    code: ErrorCode::InvalidInput,
                    message: err.to_string(),
                }
            }
            RnaStructureError::ToolNotFound { .. } => EngineError {
                code: ErrorCode::Unsupported,
                message: err.to_string(),
            },
            RnaStructureError::ToolFailed { .. } => EngineError {
                code: ErrorCode::Internal,
                message: err.to_string(),
            },
            RnaStructureError::Io { .. } => EngineError {
                code: ErrorCode::Io,
                message: err.to_string(),
            },
        }
    }

    pub fn inspect_rna_structure(
        &self,
        seq_id: &str,
    ) -> Result<RnaStructureTextReport, EngineError> {
        let dna = self
            .state
            .sequences
            .get(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{seq_id}' not found"),
            })?;
        rna_structure::inspect_text(dna).map_err(Self::map_rna_structure_error)
    }

    pub fn render_rna_structure_svg_to_path(
        &self,
        seq_id: &str,
        path: &str,
    ) -> Result<RnaStructureSvgReport, EngineError> {
        let dna = self
            .state
            .sequences
            .get(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{seq_id}' not found"),
            })?;
        rna_structure::render_svg(dna, path).map_err(Self::map_rna_structure_error)
    }

    fn open_reference_genome_catalog(
        catalog_path: Option<&str>,
    ) -> Result<(GenomeCatalog, String), EngineError> {
        let catalog_path = catalog_path
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .unwrap_or(DEFAULT_GENOME_CATALOG_PATH)
            .to_string();
        let catalog = GenomeCatalog::from_json_file(&catalog_path).map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!("Could not open genome catalog '{catalog_path}': {e}"),
        })?;
        Ok((catalog, catalog_path))
    }

    pub fn list_reference_genomes(catalog_path: Option<&str>) -> Result<Vec<String>, EngineError> {
        let (catalog, _) = Self::open_reference_genome_catalog(catalog_path)?;
        Ok(catalog.list_genomes())
    }

    pub fn describe_reference_genome_sources(
        catalog_path: Option<&str>,
        genome_id: &str,
        cache_dir: Option<&str>,
    ) -> Result<GenomeSourcePlan, EngineError> {
        let (catalog, catalog_path) = Self::open_reference_genome_catalog(catalog_path)?;
        catalog
            .source_plan(
                genome_id,
                cache_dir.map(str::trim).filter(|v| !v.is_empty()),
            )
            .map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not resolve source plan for genome '{}' in catalog '{}': {}",
                    genome_id, catalog_path, e
                ),
            })
    }

    pub fn list_helper_genomes(catalog_path: Option<&str>) -> Result<Vec<String>, EngineError> {
        let chosen = catalog_path
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .unwrap_or(DEFAULT_HELPER_GENOME_CATALOG_PATH);
        Self::list_reference_genomes(Some(chosen))
    }

    pub fn describe_helper_genome_sources(
        genome_id: &str,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
    ) -> Result<GenomeSourcePlan, EngineError> {
        let chosen = catalog_path
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .unwrap_or(DEFAULT_HELPER_GENOME_CATALOG_PATH);
        Self::describe_reference_genome_sources(Some(chosen), genome_id, cache_dir)
    }

    pub fn is_reference_genome_prepared(
        catalog_path: Option<&str>,
        genome_id: &str,
        cache_dir: Option<&str>,
    ) -> Result<bool, EngineError> {
        let (catalog, _) = Self::open_reference_genome_catalog(catalog_path)?;
        catalog
            .is_prepared(
                genome_id,
                cache_dir.map(str::trim).filter(|v| !v.is_empty()),
            )
            .map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not check prepared status for genome '{}': {}",
                    genome_id, e
                ),
            })
    }

    pub fn is_helper_genome_prepared(
        genome_id: &str,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
    ) -> Result<bool, EngineError> {
        let chosen = catalog_path
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .unwrap_or(DEFAULT_HELPER_GENOME_CATALOG_PATH);
        Self::is_reference_genome_prepared(Some(chosen), genome_id, cache_dir)
    }

    pub fn list_reference_genome_genes(
        catalog_path: Option<&str>,
        genome_id: &str,
        cache_dir: Option<&str>,
    ) -> Result<Vec<GenomeGeneRecord>, EngineError> {
        let (catalog, catalog_path) = Self::open_reference_genome_catalog(catalog_path)?;
        catalog
            .list_gene_regions(
                genome_id,
                cache_dir.map(str::trim).filter(|v| !v.is_empty()),
            )
            .map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not list genes for genome '{}' in catalog '{}': {}",
                    genome_id, catalog_path, e
                ),
            })
    }

    pub fn list_helper_genome_features(
        genome_id: &str,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
    ) -> Result<Vec<GenomeGeneRecord>, EngineError> {
        let chosen = catalog_path
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .unwrap_or(DEFAULT_HELPER_GENOME_CATALOG_PATH);
        Self::list_reference_genome_genes(Some(chosen), genome_id, cache_dir)
    }

    pub fn blast_reference_genome(
        catalog_path: Option<&str>,
        genome_id: &str,
        query_sequence: &str,
        max_hits: usize,
        task: Option<&str>,
        cache_dir: Option<&str>,
    ) -> Result<GenomeBlastReport, EngineError> {
        let (catalog, catalog_path) = Self::open_reference_genome_catalog(catalog_path)?;
        catalog
            .blast_sequence_with_cache(
                genome_id,
                query_sequence,
                max_hits,
                task,
                cache_dir.map(str::trim).filter(|v| !v.is_empty()),
            )
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not run BLAST search against genome '{}' in catalog '{}': {}",
                    genome_id, catalog_path, e
                ),
            })
    }

    pub fn blast_helper_genome(
        genome_id: &str,
        query_sequence: &str,
        max_hits: usize,
        task: Option<&str>,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
    ) -> Result<GenomeBlastReport, EngineError> {
        let chosen = catalog_path
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .unwrap_or(DEFAULT_HELPER_GENOME_CATALOG_PATH);
        Self::blast_reference_genome(
            Some(chosen),
            genome_id,
            query_sequence,
            max_hits,
            task,
            cache_dir,
        )
    }

    pub fn prepare_reference_genome_once(
        genome_id: &str,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
        on_progress: &mut dyn FnMut(PrepareGenomeProgress),
    ) -> Result<PrepareGenomeReport, EngineError> {
        let (catalog, _) = Self::open_reference_genome_catalog(catalog_path)?;
        catalog
            .prepare_genome_once_with_progress(
                genome_id,
                cache_dir.map(str::trim).filter(|v| !v.is_empty()),
                on_progress,
            )
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not prepare genome '{genome_id}': {e}"),
            })
    }

    pub fn prepare_helper_genome_once(
        genome_id: &str,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
        on_progress: &mut dyn FnMut(PrepareGenomeProgress),
    ) -> Result<PrepareGenomeReport, EngineError> {
        let chosen = catalog_path
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .unwrap_or(DEFAULT_HELPER_GENOME_CATALOG_PATH);
        Self::prepare_reference_genome_once(genome_id, Some(chosen), cache_dir, on_progress)
    }

    pub fn format_prepare_genome_message(
        genome_id: &str,
        cache_dir: Option<&str>,
        report: &PrepareGenomeReport,
    ) -> String {
        let status = if report.reused_existing {
            "reused existing local cache"
        } else {
            "downloaded and installed"
        };
        let sequence_type = report.sequence_source_type.as_deref().unwrap_or("unknown");
        let annotation_type = report
            .annotation_source_type
            .as_deref()
            .unwrap_or("unknown");
        let blast_status = if report.blast_index_ready {
            format!(
                "ready ({})",
                report
                    .blast_db_prefix
                    .as_deref()
                    .unwrap_or("unknown BLAST DB prefix")
            )
        } else {
            "not available".to_string()
        };
        format!(
            "Prepared genome '{}' ({status}). cache='{}' sequence='{}' [{}], annotation='{}' [{}], blast_index={}",
            genome_id,
            cache_dir
                .map(str::trim)
                .filter(|v| !v.is_empty())
                .unwrap_or("catalog/default"),
            report.sequence_path,
            sequence_type,
            report.annotation_path,
            annotation_type,
            blast_status
        )
    }

    pub fn operation_log(&self) -> &[OperationRecord] {
        &self.journal
    }

    pub fn apply_with_progress<F>(
        &mut self,
        op: Operation,
        mut on_progress: F,
    ) -> Result<OpResult, EngineError>
    where
        F: FnMut(OperationProgress),
    {
        let run_id = "interactive".to_string();
        let result = self.apply_internal(op.clone(), &run_id, &mut on_progress)?;
        self.journal.push(OperationRecord {
            run_id,
            op,
            result: result.clone(),
        });
        Ok(result)
    }

    pub fn apply_workflow_with_progress<F>(
        &mut self,
        wf: Workflow,
        mut on_progress: F,
    ) -> Result<Vec<OpResult>, EngineError>
    where
        F: FnMut(OperationProgress),
    {
        let mut results = Vec::new();
        for op in &wf.ops {
            let result = self.apply_internal(op.clone(), &wf.run_id, &mut on_progress)?;
            self.journal.push(OperationRecord {
                run_id: wf.run_id.clone(),
                op: op.clone(),
                result: result.clone(),
            });
            results.push(result);
        }
        Ok(results)
    }

    fn next_op_id(&mut self) -> OpId {
        self.op_counter += 1;
        format!("op-{}", self.op_counter)
    }

    fn derive_seq_id(path: &str) -> SeqId {
        Path::new(path)
            .file_stem()
            .map(|s| s.to_string_lossy().to_string())
            .filter(|s| !s.is_empty())
            .unwrap_or_else(|| "sequence".to_string())
    }

    fn normalize_id_token(raw: &str) -> String {
        let mut out = String::new();
        for c in raw.chars() {
            if c.is_ascii_alphanumeric() {
                out.push(c.to_ascii_lowercase());
            } else if matches!(c, '_' | '-' | '.') {
                out.push('_');
            }
        }
        let trimmed = out.trim_matches('_');
        if trimmed.is_empty() {
            "region".to_string()
        } else {
            trimmed.to_string()
        }
    }

    fn genome_gene_matches_exact(record: &GenomeGeneRecord, query: &str) -> bool {
        record
            .gene_name
            .as_ref()
            .map(|v| v.eq_ignore_ascii_case(query))
            .unwrap_or(false)
            || record
                .gene_id
                .as_ref()
                .map(|v| v.eq_ignore_ascii_case(query))
                .unwrap_or(false)
    }

    fn genome_gene_matches_contains(record: &GenomeGeneRecord, query_lower: &str) -> bool {
        record
            .gene_name
            .as_ref()
            .map(|v| v.to_ascii_lowercase().contains(query_lower))
            .unwrap_or(false)
            || record
                .gene_id
                .as_ref()
                .map(|v| v.to_ascii_lowercase().contains(query_lower))
                .unwrap_or(false)
    }

    fn genome_gene_display_label(record: &GenomeGeneRecord) -> String {
        let label = record
            .gene_name
            .as_ref()
            .or(record.gene_id.as_ref())
            .cloned()
            .unwrap_or_else(|| "unnamed_gene".to_string());
        format!(
            "{}:{}-{} ({})",
            record.chromosome, record.start_1based, record.end_1based, label
        )
    }

    fn import_genome_slice_sequence(
        &mut self,
        result: &mut OpResult,
        sequence: String,
        default_id: String,
    ) -> Result<SeqId, EngineError> {
        let mut dna = DNAsequence::from_sequence(&sequence).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not construct DNA sequence from genome slice: {e}"),
        })?;
        Self::prepare_sequence(&mut dna);
        let seq_id = self.unique_seq_id(&default_id);
        self.state.sequences.insert(seq_id.clone(), dna);
        self.add_lineage_node(
            &seq_id,
            SequenceOrigin::ImportedGenomic,
            Some(&result.op_id),
        );
        result.created_seq_ids.push(seq_id.clone());
        Ok(seq_id)
    }

    fn now_unix_ms() -> u128 {
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_millis())
            .unwrap_or(0)
    }

    fn append_genome_extraction_provenance(&mut self, entry: GenomeExtractionProvenance) {
        let mut provenance = self
            .state
            .metadata
            .get(PROVENANCE_METADATA_KEY)
            .cloned()
            .unwrap_or_else(|| json!({}));
        if !provenance.is_object() {
            provenance = json!({});
        }
        if let Some(obj) = provenance.as_object_mut() {
            let mut records: Vec<GenomeExtractionProvenance> = obj
                .get(GENOME_EXTRACTIONS_METADATA_KEY)
                .cloned()
                .and_then(|v| serde_json::from_value(v).ok())
                .unwrap_or_default();
            records.push(entry);
            if let Ok(records_value) = serde_json::to_value(records) {
                obj.insert(GENOME_EXTRACTIONS_METADATA_KEY.to_string(), records_value);
            }
            obj.insert("updated_at_unix_ms".to_string(), json!(Self::now_unix_ms()));
            self.state
                .metadata
                .insert(PROVENANCE_METADATA_KEY.to_string(), provenance);
        }
    }

    fn genome_source_snapshot(
        source_plan: Option<&GenomeSourcePlan>,
        inspection: Option<&PreparedGenomeInspection>,
    ) -> (
        Option<String>,
        Option<String>,
        Option<String>,
        Option<String>,
        Option<String>,
        Option<String>,
    ) {
        let sequence_source_type = source_plan
            .map(|p| p.sequence_source_type.clone())
            .or_else(|| inspection.map(|i| i.sequence_source_type.clone()));
        let annotation_source_type = source_plan
            .map(|p| p.annotation_source_type.clone())
            .or_else(|| inspection.map(|i| i.annotation_source_type.clone()));
        let sequence_source = source_plan
            .map(|p| p.sequence_source.clone())
            .or_else(|| inspection.map(|i| i.sequence_source.clone()));
        let annotation_source = source_plan
            .map(|p| p.annotation_source.clone())
            .or_else(|| inspection.map(|i| i.annotation_source.clone()));
        let sequence_sha1 = inspection.and_then(|i| i.sequence_sha1.clone());
        let annotation_sha1 = inspection.and_then(|i| i.annotation_sha1.clone());
        (
            sequence_source_type,
            annotation_source_type,
            sequence_source,
            annotation_source,
            sequence_sha1,
            annotation_sha1,
        )
    }

    fn normalize_genome_track_subscription(
        mut subscription: GenomeTrackSubscription,
    ) -> Result<GenomeTrackSubscription, EngineError> {
        subscription.path = subscription.path.trim().to_string();
        if subscription.path.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Track path cannot be empty".to_string(),
            });
        }
        if subscription
            .min_score
            .zip(subscription.max_score)
            .map(|(min, max)| min > max)
            .unwrap_or(false)
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Track subscription requires min_score <= max_score".to_string(),
            });
        }
        subscription.track_name = subscription
            .track_name
            .map(|v| v.trim().to_string())
            .filter(|v| !v.is_empty());
        Ok(subscription)
    }

    fn sort_track_subscriptions(subscriptions: &mut Vec<GenomeTrackSubscription>) {
        subscriptions.sort_by(|a, b| {
            a.source
                .label()
                .cmp(b.source.label())
                .then(a.path.cmp(&b.path))
                .then(
                    a.track_name
                        .clone()
                        .unwrap_or_default()
                        .cmp(&b.track_name.clone().unwrap_or_default()),
                )
        });
    }

    fn read_track_subscriptions_from_metadata(
        value: Option<&serde_json::Value>,
    ) -> Vec<GenomeTrackSubscription> {
        let mut subscriptions = value
            .cloned()
            .and_then(|v| serde_json::from_value::<Vec<GenomeTrackSubscription>>(v).ok())
            .unwrap_or_default();
        subscriptions.retain(|subscription| !subscription.path.trim().is_empty());
        Self::sort_track_subscriptions(&mut subscriptions);
        subscriptions
    }

    fn write_track_subscriptions_to_metadata(
        &mut self,
        subscriptions: &[GenomeTrackSubscription],
    ) -> Result<(), EngineError> {
        if subscriptions.is_empty() {
            self.state
                .metadata
                .remove(GENOME_TRACK_SUBSCRIPTIONS_METADATA_KEY);
            return Ok(());
        }
        let value = serde_json::to_value(subscriptions).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize genome track subscriptions: {e}"),
        })?;
        self.state
            .metadata
            .insert(GENOME_TRACK_SUBSCRIPTIONS_METADATA_KEY.to_string(), value);
        Ok(())
    }

    fn read_known_track_anchor_ids(&self) -> BTreeSet<String> {
        self.state
            .metadata
            .get(GENOME_TRACK_KNOWN_ANCHORS_METADATA_KEY)
            .and_then(|v| v.as_array())
            .map(|values| {
                values
                    .iter()
                    .filter_map(|v| v.as_str())
                    .map(str::trim)
                    .filter(|v| !v.is_empty())
                    .map(str::to_string)
                    .collect::<BTreeSet<_>>()
            })
            .unwrap_or_default()
    }

    fn write_known_track_anchor_ids(&mut self, anchors: &BTreeSet<String>) {
        if anchors.is_empty() {
            self.state
                .metadata
                .remove(GENOME_TRACK_KNOWN_ANCHORS_METADATA_KEY);
            return;
        }
        let values = anchors
            .iter()
            .map(|v| serde_json::Value::String(v.clone()))
            .collect::<Vec<_>>();
        let new_value = serde_json::Value::Array(values);
        if self
            .state
            .metadata
            .get(GENOME_TRACK_KNOWN_ANCHORS_METADATA_KEY)
            == Some(&new_value)
        {
            return;
        }
        self.state.metadata.insert(
            GENOME_TRACK_KNOWN_ANCHORS_METADATA_KEY.to_string(),
            new_value,
        );
    }

    fn track_subscription_to_operation(
        seq_id: &str,
        subscription: &GenomeTrackSubscription,
    ) -> Operation {
        match subscription.source {
            GenomeTrackSource::Bed => Operation::ImportGenomeBedTrack {
                seq_id: seq_id.to_string(),
                path: subscription.path.clone(),
                track_name: subscription.track_name.clone(),
                min_score: subscription.min_score,
                max_score: subscription.max_score,
                clear_existing: Some(subscription.clear_existing),
            },
            GenomeTrackSource::BigWig => Operation::ImportGenomeBigWigTrack {
                seq_id: seq_id.to_string(),
                path: subscription.path.clone(),
                track_name: subscription.track_name.clone(),
                min_score: subscription.min_score,
                max_score: subscription.max_score,
                clear_existing: Some(subscription.clear_existing),
            },
            GenomeTrackSource::Vcf => Operation::ImportGenomeVcfTrack {
                seq_id: seq_id.to_string(),
                path: subscription.path.clone(),
                track_name: subscription.track_name.clone(),
                min_score: subscription.min_score,
                max_score: subscription.max_score,
                clear_existing: Some(subscription.clear_existing),
            },
        }
    }

    fn apply_track_subscription_to_seq_ids(
        &mut self,
        seq_ids: &[String],
        subscription: &GenomeTrackSubscription,
    ) -> GenomeTrackSyncReport {
        let mut report = GenomeTrackSyncReport {
            subscriptions_considered: 1,
            target_sequences: seq_ids.len(),
            ..GenomeTrackSyncReport::default()
        };
        for seq_id in seq_ids {
            let op = Self::track_subscription_to_operation(seq_id, subscription);
            match self.apply(op) {
                Ok(op_result) => {
                    report.applied_imports += 1;
                    report.warnings_count += op_result.warnings.len();
                }
                Err(e) => {
                    report.failed_imports += 1;
                    if report.errors.len() < 20 {
                        report.errors.push(format!(
                            "{} '{}' @ {}: {}",
                            subscription.source.label(),
                            subscription.path,
                            seq_id,
                            e.message
                        ));
                    }
                }
            }
        }
        report
    }

    fn normalize_candidate_set_name(raw: &str) -> Result<String, EngineError> {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Candidate set name cannot be empty".to_string(),
            });
        }
        Ok(trimmed.to_string())
    }

    fn normalize_metric_name(raw: &str) -> String {
        let mut out = String::with_capacity(raw.len());
        for ch in raw.chars() {
            if ch.is_ascii_alphanumeric() {
                out.push(ch.to_ascii_lowercase());
            } else {
                out.push('_');
            }
        }
        let trimmed = out.trim_matches('_');
        if trimmed.is_empty() {
            "metric".to_string()
        } else {
            trimmed.to_string()
        }
    }

    fn read_candidate_store_from_metadata(value: Option<&serde_json::Value>) -> CandidateStore {
        let mut store = value
            .cloned()
            .and_then(|v| serde_json::from_value::<CandidateStore>(v).ok())
            .unwrap_or_default();
        if store.schema.trim().is_empty() {
            store.schema = CANDIDATE_SETS_SCHEMA.to_string();
        }
        store
    }

    fn read_candidate_store(&self) -> CandidateStore {
        Self::read_candidate_store_from_metadata(
            self.state.metadata.get(CANDIDATE_SETS_METADATA_KEY),
        )
    }

    fn write_candidate_store(&mut self, mut store: CandidateStore) -> Result<(), EngineError> {
        if store.sets.is_empty() {
            self.state.metadata.remove(CANDIDATE_SETS_METADATA_KEY);
            return Ok(());
        }
        store.schema = CANDIDATE_SETS_SCHEMA.to_string();
        store.updated_at_unix_ms = Self::now_unix_ms();
        let value = serde_json::to_value(store).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize candidate-set metadata: {e}"),
        })?;
        self.state
            .metadata
            .insert(CANDIDATE_SETS_METADATA_KEY.to_string(), value);
        Ok(())
    }

    fn metric_names_for_candidate_set(set: &CandidateSet) -> Vec<String> {
        let mut names = BTreeSet::new();
        for candidate in &set.candidates {
            for metric in candidate.metrics.keys() {
                names.insert(metric.to_string());
            }
        }
        names.into_iter().collect()
    }

    fn candidate_key(record: &CandidateRecord) -> String {
        format!(
            "{}:{}:{}:{}",
            record.seq_id, record.start_0based, record.end_0based, record.sequence
        )
    }

    pub fn list_candidate_sets(&self) -> Vec<CandidateSetSummary> {
        let store = self.read_candidate_store();
        let mut set_names: Vec<String> = store.sets.keys().cloned().collect();
        set_names.sort();
        set_names
            .iter()
            .filter_map(|name| store.sets.get(name))
            .map(|set| CandidateSetSummary {
                name: set.name.clone(),
                created_at_unix_ms: set.created_at_unix_ms,
                source_seq_ids: set.source_seq_ids.clone(),
                candidate_count: set.candidates.len(),
                metrics: Self::metric_names_for_candidate_set(set),
            })
            .collect()
    }

    pub fn inspect_candidate_set_page(
        &self,
        set_name: &str,
        limit: usize,
        offset: usize,
    ) -> Result<(CandidateSet, usize, usize), EngineError> {
        if limit == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Candidate page limit must be >= 1".to_string(),
            });
        }
        let set_name = Self::normalize_candidate_set_name(set_name)?;
        let store = self.read_candidate_store();
        let mut set = store
            .sets
            .get(&set_name)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Candidate set '{}' not found", set_name),
            })?;
        let total = set.candidates.len();
        let clamped_offset = offset.min(total);
        set.candidates = set
            .candidates
            .into_iter()
            .skip(clamped_offset)
            .take(limit)
            .collect();
        Ok((set, total, clamped_offset))
    }

    pub fn list_candidate_set_metrics(
        &self,
        set_name: &str,
    ) -> Result<Vec<CandidateMetricSummary>, EngineError> {
        let set_name = Self::normalize_candidate_set_name(set_name)?;
        let store = self.read_candidate_store();
        let set = store.sets.get(&set_name).ok_or_else(|| EngineError {
            code: ErrorCode::NotFound,
            message: format!("Candidate set '{}' not found", set_name),
        })?;
        let mut counts: HashMap<String, usize> = HashMap::new();
        for candidate in &set.candidates {
            for metric in candidate.metrics.keys() {
                *counts.entry(metric.to_string()).or_insert(0) += 1;
            }
        }
        let mut summaries = counts
            .into_iter()
            .map(|(metric, present_in_candidates)| CandidateMetricSummary {
                metric,
                present_in_candidates,
                missing_in_candidates: set.candidates.len().saturating_sub(present_in_candidates),
            })
            .collect::<Vec<_>>();
        summaries.sort_by(|a, b| a.metric.cmp(&b.metric));
        Ok(summaries)
    }

    fn sequence_base_counts(sequence: &[u8]) -> (usize, usize, usize, usize, usize) {
        let mut a = 0usize;
        let mut c = 0usize;
        let mut g = 0usize;
        let mut t = 0usize;
        for b in sequence {
            match b.to_ascii_uppercase() {
                b'A' => a += 1,
                b'C' => c += 1,
                b'G' => g += 1,
                b'T' => t += 1,
                _ => {}
            }
        }
        let canonical = a + c + g + t;
        (a, c, g, t, canonical)
    }

    fn approximate_molecular_weight_da(sequence: &[u8]) -> f64 {
        let mut total = 0.0;
        for b in sequence {
            total += match b.to_ascii_uppercase() {
                b'A' => 313.21,
                b'C' => 289.18,
                b'G' => 329.21,
                b'T' => 304.20,
                _ => 0.0,
            };
        }
        total
    }

    fn feature_labels_upper(feature: &gb_io::seq::Feature) -> Vec<String> {
        let mut labels = vec![];
        for key in [
            "label",
            "gene",
            "gene_name",
            "locus_tag",
            "product",
            "standard_name",
            "note",
        ] {
            for value in feature.qualifier_values(key.into()) {
                let trimmed = value.trim();
                if !trimmed.is_empty() {
                    labels.push(trimmed.to_ascii_uppercase());
                }
            }
        }
        labels
    }

    fn collect_feature_intervals(dna: &DNAsequence) -> Vec<FeatureInterval> {
        let mut out = vec![];
        for feature in dna.features() {
            if feature.kind.to_string().eq_ignore_ascii_case("SOURCE") {
                continue;
            }
            let Ok((from, to)) = feature.location.find_bounds() else {
                continue;
            };
            if from < 0 || to < 0 {
                continue;
            }
            let start = from.min(to) as usize;
            let end_inclusive = from.max(to) as usize;
            if dna.is_empty() || start >= dna.len() {
                continue;
            }
            let end_0based = end_inclusive.saturating_add(1).min(dna.len());
            if end_0based <= start {
                continue;
            }
            out.push(FeatureInterval {
                kind_upper: feature.kind.to_string().to_ascii_uppercase(),
                labels_upper: Self::feature_labels_upper(feature),
                start_0based: start,
                end_0based,
            });
        }
        out
    }

    fn compile_optional_regex(
        pattern: &Option<String>,
        option_name: &str,
    ) -> Result<Option<Regex>, EngineError> {
        let Some(raw) = pattern else {
            return Ok(None);
        };
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return Ok(None);
        }
        RegexBuilder::new(trimmed)
            .case_insensitive(true)
            .build()
            .map(Some)
            .map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Invalid {option_name} regex '{}': {}", trimmed, e),
            })
    }

    fn feature_matches_filter(
        feature: &FeatureInterval,
        kind_filter_upper: &[String],
        label_regex: Option<&Regex>,
    ) -> bool {
        let kind_ok = kind_filter_upper.is_empty()
            || kind_filter_upper
                .iter()
                .any(|kind| feature.kind_upper == *kind);
        if !kind_ok {
            return false;
        }
        if let Some(regex) = label_regex {
            feature
                .labels_upper
                .iter()
                .any(|label| regex.is_match(label))
        } else {
            true
        }
    }

    fn interval_distance(
        left_start: usize,
        left_end: usize,
        right_start: usize,
        right_end: usize,
    ) -> usize {
        if left_end <= right_start {
            right_start.saturating_sub(left_end)
        } else if right_end <= left_start {
            left_start.saturating_sub(right_end)
        } else {
            0
        }
    }

    fn nearest_feature_distance(
        candidate_start: usize,
        candidate_end: usize,
        features: &[FeatureInterval],
        kind_filter_upper: &[String],
        label_regex: Option<&Regex>,
    ) -> Option<usize> {
        features
            .iter()
            .filter(|feature| Self::feature_matches_filter(feature, kind_filter_upper, label_regex))
            .map(|feature| {
                Self::interval_distance(
                    candidate_start,
                    candidate_end,
                    feature.start_0based,
                    feature.end_0based,
                )
            })
            .min()
    }

    fn compute_candidate_metrics(
        sequence: &[u8],
        start_0based: usize,
        end_0based: usize,
        source_len: usize,
        nearest_feature_distance_bp: Option<usize>,
    ) -> HashMap<String, f64> {
        let mut metrics = HashMap::new();
        let length_bp = end_0based.saturating_sub(start_0based);
        let (count_a, count_c, count_g, count_t, canonical) = Self::sequence_base_counts(sequence);
        metrics.insert("length_bp".to_string(), length_bp as f64);
        metrics.insert(
            "molecular_weight_da".to_string(),
            Self::approximate_molecular_weight_da(sequence),
        );
        metrics.insert("distance_to_seq_start_bp".to_string(), start_0based as f64);
        metrics.insert(
            "distance_to_seq_end_bp".to_string(),
            source_len.saturating_sub(end_0based) as f64,
        );
        metrics.insert("count_a".to_string(), count_a as f64);
        metrics.insert("count_c".to_string(), count_c as f64);
        metrics.insert("count_g".to_string(), count_g as f64);
        metrics.insert("count_t".to_string(), count_t as f64);
        if canonical > 0 {
            let canonical_f = canonical as f64;
            metrics.insert(
                "gc_fraction".to_string(),
                (count_g + count_c) as f64 / canonical_f,
            );
            metrics.insert(
                "at_fraction".to_string(),
                (count_a + count_t) as f64 / canonical_f,
            );
            metrics.insert("a_fraction".to_string(), count_a as f64 / canonical_f);
            metrics.insert("c_fraction".to_string(), count_c as f64 / canonical_f);
            metrics.insert("g_fraction".to_string(), count_g as f64 / canonical_f);
            metrics.insert("t_fraction".to_string(), count_t as f64 / canonical_f);
        }
        if count_t > 0 {
            metrics.insert("at_ratio".to_string(), count_a as f64 / count_t as f64);
        }
        if count_c > 0 {
            metrics.insert("gc_ratio".to_string(), count_g as f64 / count_c as f64);
        }
        if let Some(distance) = nearest_feature_distance_bp {
            metrics.insert(
                "distance_to_nearest_feature_bp".to_string(),
                distance as f64,
            );
        }
        metrics.retain(|_, value| value.is_finite());
        metrics
    }

    fn quantile_threshold(values: &[f64], quantile: f64) -> Option<f64> {
        if values.is_empty() || !quantile.is_finite() || !(0.0..=1.0).contains(&quantile) {
            return None;
        }
        let mut sorted = values.to_vec();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
        let max_idx = sorted.len().saturating_sub(1);
        let idx = ((max_idx as f64) * quantile).round() as usize;
        sorted.get(idx.min(max_idx)).copied()
    }

    fn tokenize_expression(expression: &str) -> Result<Vec<ExpressionToken>, EngineError> {
        let bytes = expression.as_bytes();
        let mut idx = 0usize;
        let mut tokens = vec![];
        while idx < bytes.len() {
            let b = bytes[idx];
            if b.is_ascii_whitespace() {
                idx += 1;
                continue;
            }
            match b {
                b'+' => {
                    tokens.push(ExpressionToken::Plus);
                    idx += 1;
                }
                b'-' => {
                    tokens.push(ExpressionToken::Minus);
                    idx += 1;
                }
                b'*' => {
                    tokens.push(ExpressionToken::Star);
                    idx += 1;
                }
                b'/' => {
                    tokens.push(ExpressionToken::Slash);
                    idx += 1;
                }
                b'(' => {
                    tokens.push(ExpressionToken::LParen);
                    idx += 1;
                }
                b')' => {
                    tokens.push(ExpressionToken::RParen);
                    idx += 1;
                }
                b',' => {
                    tokens.push(ExpressionToken::Comma);
                    idx += 1;
                }
                _ if b.is_ascii_digit() || b == b'.' => {
                    let start = idx;
                    idx += 1;
                    while idx < bytes.len()
                        && (bytes[idx].is_ascii_digit()
                            || bytes[idx] == b'.'
                            || matches!(bytes[idx], b'e' | b'E' | b'+' | b'-'))
                    {
                        if matches!(bytes[idx], b'+' | b'-') {
                            let prev = bytes[idx.saturating_sub(1)];
                            if !matches!(prev, b'e' | b'E') {
                                break;
                            }
                        }
                        idx += 1;
                    }
                    let raw = &expression[start..idx];
                    let value = raw.parse::<f64>().map_err(|e| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("Invalid numeric literal '{}': {}", raw, e),
                    })?;
                    tokens.push(ExpressionToken::Number(value));
                }
                _ if b.is_ascii_alphabetic() || b == b'_' => {
                    let start = idx;
                    idx += 1;
                    while idx < bytes.len()
                        && (bytes[idx].is_ascii_alphanumeric() || bytes[idx] == b'_')
                    {
                        idx += 1;
                    }
                    tokens.push(ExpressionToken::Ident(expression[start..idx].to_string()));
                }
                _ => {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("Unsupported character '{}' in expression", b as char),
                    })
                }
            }
        }
        if tokens.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Expression is empty".to_string(),
            });
        }
        Ok(tokens)
    }

    fn parse_metric_expression(expression: &str) -> Result<MetricExpr, EngineError> {
        let tokens = Self::tokenize_expression(expression)?;
        MetricExpressionParser::new(tokens)
            .parse()
            .map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Invalid score expression: {e}"),
            })
    }

    fn evaluate_metric_expression(
        expr: &MetricExpr,
        metrics: &HashMap<String, f64>,
    ) -> Result<f64, EngineError> {
        match expr {
            MetricExpr::Number(value) => Ok(*value),
            MetricExpr::Variable(name) => {
                if let Some(value) = metrics.get(name) {
                    return Ok(*value);
                }
                let normalized = Self::normalize_metric_name(name);
                metrics
                    .get(&normalized)
                    .copied()
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("Expression references unknown metric '{}'", name),
                    })
            }
            MetricExpr::UnaryMinus(inner) => Ok(-Self::evaluate_metric_expression(inner, metrics)?),
            MetricExpr::Binary { op, left, right } => {
                let lhs = Self::evaluate_metric_expression(left, metrics)?;
                let rhs = Self::evaluate_metric_expression(right, metrics)?;
                let value = match op {
                    ExpressionBinaryOp::Add => lhs + rhs,
                    ExpressionBinaryOp::Subtract => lhs - rhs,
                    ExpressionBinaryOp::Multiply => lhs * rhs,
                    ExpressionBinaryOp::Divide => {
                        if rhs.abs() < f64::EPSILON {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "Expression division by zero".to_string(),
                            });
                        }
                        lhs / rhs
                    }
                };
                if value.is_finite() {
                    Ok(value)
                } else {
                    Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "Expression produced a non-finite value".to_string(),
                    })
                }
            }
            MetricExpr::Function { name, args } => {
                let normalized = name.trim().to_ascii_lowercase();
                let value = match normalized.as_str() {
                    "abs" => {
                        if args.len() != 1 {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "Function abs() expects exactly 1 argument".to_string(),
                            });
                        }
                        Self::evaluate_metric_expression(&args[0], metrics)?.abs()
                    }
                    "sqrt" => {
                        if args.len() != 1 {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "Function sqrt() expects exactly 1 argument".to_string(),
                            });
                        }
                        let x = Self::evaluate_metric_expression(&args[0], metrics)?;
                        if x < 0.0 {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "Function sqrt() requires non-negative input".to_string(),
                            });
                        }
                        x.sqrt()
                    }
                    "log" => {
                        if args.len() != 1 {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "Function log() expects exactly 1 argument".to_string(),
                            });
                        }
                        let x = Self::evaluate_metric_expression(&args[0], metrics)?;
                        if x <= 0.0 {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "Function log() requires positive input".to_string(),
                            });
                        }
                        x.ln()
                    }
                    "exp" => {
                        if args.len() != 1 {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "Function exp() expects exactly 1 argument".to_string(),
                            });
                        }
                        Self::evaluate_metric_expression(&args[0], metrics)?.exp()
                    }
                    "min" => {
                        if args.len() != 2 {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "Function min() expects exactly 2 arguments".to_string(),
                            });
                        }
                        let a = Self::evaluate_metric_expression(&args[0], metrics)?;
                        let b = Self::evaluate_metric_expression(&args[1], metrics)?;
                        a.min(b)
                    }
                    "max" => {
                        if args.len() != 2 {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "Function max() expects exactly 2 arguments".to_string(),
                            });
                        }
                        let a = Self::evaluate_metric_expression(&args[0], metrics)?;
                        let b = Self::evaluate_metric_expression(&args[1], metrics)?;
                        a.max(b)
                    }
                    "pow" => {
                        if args.len() != 2 {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "Function pow() expects exactly 2 arguments".to_string(),
                            });
                        }
                        let a = Self::evaluate_metric_expression(&args[0], metrics)?;
                        let b = Self::evaluate_metric_expression(&args[1], metrics)?;
                        a.powf(b)
                    }
                    "clamp" => {
                        if args.len() != 3 {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "Function clamp() expects exactly 3 arguments".to_string(),
                            });
                        }
                        let x = Self::evaluate_metric_expression(&args[0], metrics)?;
                        let lo = Self::evaluate_metric_expression(&args[1], metrics)?;
                        let hi = Self::evaluate_metric_expression(&args[2], metrics)?;
                        if lo > hi {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "Function clamp() requires lo <= hi".to_string(),
                            });
                        }
                        x.max(lo).min(hi)
                    }
                    _ => {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!("Unknown expression function '{}'", name),
                        })
                    }
                };
                if value.is_finite() {
                    Ok(value)
                } else {
                    Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("Function '{}' produced a non-finite value", name),
                    })
                }
            }
        }
    }

    fn op_generate_candidate_set(
        &mut self,
        set_name: String,
        seq_id: SeqId,
        length_bp: usize,
        step_bp: usize,
        feature_kinds: Vec<String>,
        feature_label_regex: Option<String>,
        max_distance_bp: Option<usize>,
        limit: Option<usize>,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let set_name = Self::normalize_candidate_set_name(&set_name)?;
        if length_bp == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "GenerateCandidateSet requires length_bp >= 1".to_string(),
            });
        }
        if step_bp == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "GenerateCandidateSet requires step_bp >= 1".to_string(),
            });
        }
        let limit = limit.unwrap_or(self.max_fragments_per_container());
        if limit == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "GenerateCandidateSet requires limit >= 1".to_string(),
            });
        }
        let dna = self
            .state
            .sequences
            .get(&seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{}' not found", seq_id),
            })?;
        if dna.len() < length_bp {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "GenerateCandidateSet length_bp {} exceeds sequence '{}' length {}",
                    length_bp,
                    seq_id,
                    dna.len()
                ),
            });
        }
        let label_regex =
            Self::compile_optional_regex(&feature_label_regex, "feature_label_regex")?;
        let mut kind_filter_upper = feature_kinds
            .iter()
            .map(|kind| kind.trim().to_ascii_uppercase())
            .filter(|kind| !kind.is_empty())
            .collect::<Vec<_>>();
        kind_filter_upper.sort();
        kind_filter_upper.dedup();
        let feature_intervals = Self::collect_feature_intervals(dna);
        let matching_feature_count = feature_intervals
            .iter()
            .filter(|feature| {
                Self::feature_matches_filter(feature, &kind_filter_upper, label_regex.as_ref())
            })
            .count();
        if (!kind_filter_upper.is_empty() || label_regex.is_some() || max_distance_bp.is_some())
            && matching_feature_count == 0
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "No features matched feature filters while generating candidate set"
                    .to_string(),
            });
        }

        let mut candidates = vec![];
        let mut considered = 0usize;
        let mut truncated = false;
        let upper = dna.len().saturating_sub(length_bp);
        let mut start = 0usize;
        while start <= upper {
            let end = start + length_bp;
            considered += 1;
            let distance_any =
                Self::nearest_feature_distance(start, end, &feature_intervals, &[], None);
            let distance_filtered = Self::nearest_feature_distance(
                start,
                end,
                &feature_intervals,
                &kind_filter_upper,
                label_regex.as_ref(),
            );
            let selected_distance = if !kind_filter_upper.is_empty() || label_regex.is_some() {
                distance_filtered
            } else {
                distance_any
            };
            if let Some(max_distance) = max_distance_bp {
                let Some(distance) = selected_distance else {
                    start = start.saturating_add(step_bp);
                    continue;
                };
                if distance > max_distance {
                    start = start.saturating_add(step_bp);
                    continue;
                }
            }

            let Some(fragment) = dna.get_range_safe(start..end) else {
                start = start.saturating_add(step_bp);
                continue;
            };
            let sequence = String::from_utf8_lossy(&fragment).to_string();
            let mut metrics =
                Self::compute_candidate_metrics(&fragment, start, end, dna.len(), distance_any);
            if let Some(distance) = selected_distance {
                metrics.insert(
                    "distance_to_filtered_feature_bp".to_string(),
                    distance as f64,
                );
            }
            candidates.push(CandidateRecord {
                seq_id: seq_id.clone(),
                start_0based: start,
                end_0based: end,
                sequence,
                metrics,
            });
            if candidates.len() >= limit {
                truncated = true;
                break;
            }
            start = start.saturating_add(step_bp);
        }
        if candidates.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "No candidates matched generation constraints".to_string(),
            });
        }
        let generated = candidates.len();
        let mut store = self.read_candidate_store();
        let replaced_existing = store
            .sets
            .insert(
                set_name.clone(),
                CandidateSet {
                    name: set_name.clone(),
                    created_at_unix_ms: Self::now_unix_ms(),
                    source_seq_ids: vec![seq_id.clone()],
                    candidates,
                },
            )
            .is_some();
        let metric_names = store
            .sets
            .get(&set_name)
            .map(Self::metric_names_for_candidate_set)
            .unwrap_or_default();
        self.write_candidate_store(store)?;
        result.messages.push(format!(
            "Generated candidate set '{}' from '{}' ({} candidates, {} windows considered)",
            set_name, seq_id, generated, considered
        ));
        if truncated {
            result.warnings.push(format!(
                "Candidate generation for '{}' was truncated at limit={}",
                set_name, limit
            ));
        }
        if replaced_existing {
            result.warnings.push(format!(
                "Candidate set '{}' replaced existing set",
                set_name
            ));
        }
        if !metric_names.is_empty() {
            result.messages.push(format!(
                "Candidate set '{}' metrics: {}",
                set_name,
                metric_names.join(", ")
            ));
        }
        if matching_feature_count > 0 {
            result.messages.push(format!(
                "Candidate set '{}' matching feature count: {}",
                set_name, matching_feature_count
            ));
        }
        Ok(())
    }

    fn op_delete_candidate_set(
        &mut self,
        set_name: String,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let set_name = Self::normalize_candidate_set_name(&set_name)?;
        let mut store = self.read_candidate_store();
        let removed = store.sets.remove(&set_name).is_some();
        self.write_candidate_store(store)?;
        if removed {
            result
                .messages
                .push(format!("Deleted candidate set '{}'", set_name));
        } else {
            result
                .warnings
                .push(format!("Candidate set '{}' was not present", set_name));
        }
        Ok(())
    }

    fn op_score_candidate_set_expression(
        &mut self,
        set_name: String,
        metric: String,
        expression: String,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let set_name = Self::normalize_candidate_set_name(&set_name)?;
        let metric_name = Self::normalize_metric_name(&metric);
        let expr = Self::parse_metric_expression(&expression)?;
        let mut store = self.read_candidate_store();
        let set = store.sets.get_mut(&set_name).ok_or_else(|| EngineError {
            code: ErrorCode::NotFound,
            message: format!("Candidate set '{}' not found", set_name),
        })?;
        if set.candidates.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Candidate set '{}' is empty", set_name),
            });
        }
        let mut values = Vec::with_capacity(set.candidates.len());
        for (idx, candidate) in set.candidates.iter().enumerate() {
            let value =
                Self::evaluate_metric_expression(&expr, &candidate.metrics).map_err(|e| {
                    EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Could not evaluate expression for candidate {} in '{}': {}",
                            idx, set_name, e.message
                        ),
                    }
                })?;
            values.push(value);
        }
        for (candidate, value) in set.candidates.iter_mut().zip(values.iter()) {
            candidate.metrics.insert(metric_name.clone(), *value);
        }
        let min_value = values.iter().copied().fold(f64::INFINITY, f64::min);
        let max_value = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        self.write_candidate_store(store)?;
        result.messages.push(format!(
            "Scored candidate set '{}' with metric '{}' from expression '{}'",
            set_name, metric_name, expression
        ));
        result.messages.push(format!(
            "Metric '{}' range in '{}': [{:.6}, {:.6}]",
            metric_name, set_name, min_value, max_value
        ));
        Ok(())
    }

    fn op_score_candidate_set_distance(
        &mut self,
        set_name: String,
        metric: String,
        feature_kinds: Vec<String>,
        feature_label_regex: Option<String>,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let set_name = Self::normalize_candidate_set_name(&set_name)?;
        let metric_name = Self::normalize_metric_name(&metric);
        let label_regex =
            Self::compile_optional_regex(&feature_label_regex, "feature_label_regex")?;
        let mut kind_filter_upper = feature_kinds
            .iter()
            .map(|kind| kind.trim().to_ascii_uppercase())
            .filter(|kind| !kind.is_empty())
            .collect::<Vec<_>>();
        kind_filter_upper.sort();
        kind_filter_upper.dedup();

        let mut store = self.read_candidate_store();
        let set = store.sets.get_mut(&set_name).ok_or_else(|| EngineError {
            code: ErrorCode::NotFound,
            message: format!("Candidate set '{}' not found", set_name),
        })?;
        if set.candidates.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Candidate set '{}' is empty", set_name),
            });
        }

        let mut feature_cache: HashMap<String, Vec<FeatureInterval>> = HashMap::new();
        for seq_id in set
            .candidates
            .iter()
            .map(|c| c.seq_id.clone())
            .collect::<BTreeSet<_>>()
        {
            let Some(dna) = self.state.sequences.get(&seq_id) else {
                return Err(EngineError {
                    code: ErrorCode::NotFound,
                    message: format!(
                        "Candidate set '{}' references missing sequence '{}'",
                        set_name, seq_id
                    ),
                });
            };
            feature_cache.insert(seq_id.clone(), Self::collect_feature_intervals(dna));
        }

        let mut values = Vec::with_capacity(set.candidates.len());
        for (idx, candidate) in set.candidates.iter().enumerate() {
            let features = feature_cache
                .get(&candidate.seq_id)
                .ok_or_else(|| EngineError {
                    code: ErrorCode::Internal,
                    message: format!(
                        "Missing feature cache for sequence '{}' while scoring candidate set '{}'",
                        candidate.seq_id, set_name
                    ),
                })?;
            let distance = Self::nearest_feature_distance(
                candidate.start_0based,
                candidate.end_0based,
                features,
                &kind_filter_upper,
                label_regex.as_ref(),
            )
            .ok_or_else(|| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "No matching features found for candidate {} (seq='{}')",
                    idx, candidate.seq_id
                ),
            })?;
            values.push(distance as f64);
        }

        for (candidate, value) in set.candidates.iter_mut().zip(values.iter()) {
            candidate.metrics.insert(metric_name.clone(), *value);
        }
        let min_value = values.iter().copied().fold(f64::INFINITY, f64::min);
        let max_value = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        self.write_candidate_store(store)?;
        result.messages.push(format!(
            "Scored candidate set '{}' with distance metric '{}'",
            set_name, metric_name
        ));
        result.messages.push(format!(
            "Metric '{}' range in '{}': [{:.6}, {:.6}]",
            metric_name, set_name, min_value, max_value
        ));
        Ok(())
    }

    fn op_filter_candidate_set(
        &mut self,
        input_set: String,
        output_set: String,
        metric: String,
        min: Option<f64>,
        max: Option<f64>,
        min_quantile: Option<f64>,
        max_quantile: Option<f64>,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let input_set = Self::normalize_candidate_set_name(&input_set)?;
        let output_set = Self::normalize_candidate_set_name(&output_set)?;
        let metric_name = Self::normalize_metric_name(&metric);
        if min.zip(max).map(|(a, b)| a > b).unwrap_or(false) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "FilterCandidateSet requires min <= max".to_string(),
            });
        }
        if min_quantile
            .zip(max_quantile)
            .map(|(a, b)| a > b)
            .unwrap_or(false)
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "FilterCandidateSet requires min_quantile <= max_quantile".to_string(),
            });
        }
        for (name, value) in [
            ("min_quantile", min_quantile),
            ("max_quantile", max_quantile),
        ] {
            if let Some(q) = value {
                if !(0.0..=1.0).contains(&q) {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("FilterCandidateSet {} must be between 0 and 1", name),
                    });
                }
            }
        }
        if min.is_none() && max.is_none() && min_quantile.is_none() && max_quantile.is_none() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "FilterCandidateSet requires at least one bound or quantile".to_string(),
            });
        }

        let mut store = self.read_candidate_store();
        let input = store
            .sets
            .get(&input_set)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Candidate set '{}' not found", input_set),
            })?;
        if input.candidates.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Candidate set '{}' is empty", input_set),
            });
        }

        let mut metric_values = Vec::with_capacity(input.candidates.len());
        for (idx, candidate) in input.candidates.iter().enumerate() {
            let value = candidate
                .metrics
                .get(&metric_name)
                .copied()
                .ok_or_else(|| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Candidate {} in '{}' is missing metric '{}'",
                        idx, input_set, metric_name
                    ),
                })?;
            if !value.is_finite() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Candidate {} in '{}' has non-finite metric '{}'",
                        idx, input_set, metric_name
                    ),
                });
            }
            metric_values.push(value);
        }

        let min_q_threshold =
            min_quantile.and_then(|q| Self::quantile_threshold(&metric_values, q));
        let max_q_threshold =
            max_quantile.and_then(|q| Self::quantile_threshold(&metric_values, q));
        let lower_bound = match (min, min_q_threshold) {
            (Some(a), Some(b)) => Some(a.max(b)),
            (Some(a), None) => Some(a),
            (None, Some(b)) => Some(b),
            (None, None) => None,
        };
        let upper_bound = match (max, max_q_threshold) {
            (Some(a), Some(b)) => Some(a.min(b)),
            (Some(a), None) => Some(a),
            (None, Some(b)) => Some(b),
            (None, None) => None,
        };

        let filtered_candidates = input
            .candidates
            .into_iter()
            .filter(|candidate| {
                let Some(value) = candidate.metrics.get(&metric_name).copied() else {
                    return false;
                };
                if let Some(lo) = lower_bound {
                    if value < lo {
                        return false;
                    }
                }
                if let Some(hi) = upper_bound {
                    if value > hi {
                        return false;
                    }
                }
                true
            })
            .collect::<Vec<_>>();
        let kept = filtered_candidates.len();
        let dropped = metric_values.len().saturating_sub(kept);
        let replaced_existing = store
            .sets
            .insert(
                output_set.clone(),
                CandidateSet {
                    name: output_set.clone(),
                    created_at_unix_ms: Self::now_unix_ms(),
                    source_seq_ids: input.source_seq_ids.clone(),
                    candidates: filtered_candidates,
                },
            )
            .is_some();
        self.write_candidate_store(store)?;
        result.messages.push(format!(
            "Filtered candidate set '{}' by metric '{}' into '{}' (kept {}, dropped {})",
            input_set, metric_name, output_set, kept, dropped
        ));
        if replaced_existing {
            result.warnings.push(format!(
                "FilterCandidateSet output '{}' replaced existing candidate set",
                output_set
            ));
        }
        Ok(())
    }

    fn op_candidate_set_op(
        &mut self,
        op: CandidateSetOperator,
        left_set: String,
        right_set: String,
        output_set: String,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let left_set = Self::normalize_candidate_set_name(&left_set)?;
        let right_set = Self::normalize_candidate_set_name(&right_set)?;
        let output_set = Self::normalize_candidate_set_name(&output_set)?;
        let mut store = self.read_candidate_store();
        let left = store
            .sets
            .get(&left_set)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Candidate set '{}' not found", left_set),
            })?;
        let right = store
            .sets
            .get(&right_set)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Candidate set '{}' not found", right_set),
            })?;
        let mut output_candidates = vec![];
        match op {
            CandidateSetOperator::Union => {
                let mut seen = HashSet::new();
                for candidate in left.candidates.iter().chain(right.candidates.iter()) {
                    let key = Self::candidate_key(candidate);
                    if seen.insert(key) {
                        output_candidates.push(candidate.clone());
                    }
                }
            }
            CandidateSetOperator::Intersect => {
                let right_keys = right
                    .candidates
                    .iter()
                    .map(Self::candidate_key)
                    .collect::<HashSet<_>>();
                for candidate in &left.candidates {
                    if right_keys.contains(&Self::candidate_key(candidate)) {
                        output_candidates.push(candidate.clone());
                    }
                }
            }
            CandidateSetOperator::Subtract => {
                let right_keys = right
                    .candidates
                    .iter()
                    .map(Self::candidate_key)
                    .collect::<HashSet<_>>();
                for candidate in &left.candidates {
                    if !right_keys.contains(&Self::candidate_key(candidate)) {
                        output_candidates.push(candidate.clone());
                    }
                }
            }
        }
        let mut source_seq_ids = left.source_seq_ids.clone();
        source_seq_ids.extend(right.source_seq_ids.clone());
        source_seq_ids.sort();
        source_seq_ids.dedup();
        let output_count = output_candidates.len();
        let replaced_existing = store
            .sets
            .insert(
                output_set.clone(),
                CandidateSet {
                    name: output_set.clone(),
                    created_at_unix_ms: Self::now_unix_ms(),
                    source_seq_ids,
                    candidates: output_candidates,
                },
            )
            .is_some();
        self.write_candidate_store(store)?;
        result.messages.push(format!(
            "Candidate set {} '{}' and '{}' into '{}' ({} candidates)",
            op.as_str(),
            left_set,
            right_set,
            output_set,
            output_count
        ));
        if replaced_existing {
            result.warnings.push(format!(
                "CandidateSetOp output '{}' replaced existing candidate set",
                output_set
            ));
        }
        Ok(())
    }

    fn latest_genome_anchor_for_seq(
        &self,
        seq_id: &str,
    ) -> Result<GenomeSequenceAnchor, EngineError> {
        let Some(provenance) = self.state.metadata.get(PROVENANCE_METADATA_KEY) else {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "Sequence '{seq_id}' has no genome anchor provenance (missing metadata '{}')",
                    PROVENANCE_METADATA_KEY
                ),
            });
        };
        let Some(entries) = provenance
            .get(GENOME_EXTRACTIONS_METADATA_KEY)
            .and_then(|v| v.as_array())
        else {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "Sequence '{seq_id}' has no genome anchor provenance (missing '{}')",
                    GENOME_EXTRACTIONS_METADATA_KEY
                ),
            });
        };

        let mut latest: Option<GenomeSequenceAnchor> = None;
        let mut latest_recorded_at = 0u128;

        for entry in entries {
            let Some(entry_seq_id) = entry.get("seq_id").and_then(|v| v.as_str()) else {
                continue;
            };
            if entry_seq_id != seq_id {
                continue;
            }
            let Some(chromosome) = entry.get("chromosome").and_then(|v| v.as_str()) else {
                continue;
            };
            let Some(start_1based) = entry.get("start_1based").and_then(|v| v.as_u64()) else {
                continue;
            };
            let Some(end_1based) = entry.get("end_1based").and_then(|v| v.as_u64()) else {
                continue;
            };
            if start_1based == 0 || end_1based < start_1based {
                continue;
            }
            let anchor_strand = entry
                .get("anchor_strand")
                .and_then(|v| v.as_str())
                .and_then(|v| v.trim().chars().next())
                .filter(|c| matches!(c, '+' | '-'));
            let recorded_at = entry
                .get("recorded_at_unix_ms")
                .and_then(|v| v.as_u64())
                .map(|v| v as u128)
                .unwrap_or(0);
            if latest.is_none() || recorded_at >= latest_recorded_at {
                latest_recorded_at = recorded_at;
                latest = Some(GenomeSequenceAnchor {
                    genome_id: entry
                        .get("genome_id")
                        .and_then(|v| v.as_str())
                        .unwrap_or_default()
                        .to_string(),
                    chromosome: chromosome.to_string(),
                    start_1based: start_1based as usize,
                    end_1based: end_1based as usize,
                    strand: anchor_strand,
                });
            }
        }

        latest.ok_or_else(|| EngineError {
            code: ErrorCode::NotFound,
            message: format!(
                "Sequence '{seq_id}' is not anchored to a genome interval; run ExtractGenomeRegion/ExtractGenomeGene first"
            ),
        })
    }

    fn chromosomes_match(left: &str, right: &str) -> bool {
        let normalize = |raw: &str| -> String {
            let trimmed = raw.trim();
            let lower = trimmed.to_ascii_lowercase();
            lower
                .strip_prefix("chr")
                .map(|v| v.to_string())
                .unwrap_or(lower)
        };
        normalize(left) == normalize(right)
    }

    fn is_generated_genome_bed_feature(feature: &gb_io::seq::Feature) -> bool {
        feature
            .qualifier_values("gentle_generated".into())
            .any(|v| v.eq_ignore_ascii_case(GENOME_BED_TRACK_GENERATED_TAG))
    }

    fn is_generated_genome_bigwig_feature(feature: &gb_io::seq::Feature) -> bool {
        feature
            .qualifier_values("gentle_generated".into())
            .any(|v| v.eq_ignore_ascii_case(GENOME_BIGWIG_TRACK_GENERATED_TAG))
    }

    fn is_generated_genome_vcf_feature(feature: &gb_io::seq::Feature) -> bool {
        feature
            .qualifier_values("gentle_generated".into())
            .any(|v| v.eq_ignore_ascii_case(GENOME_VCF_TRACK_GENERATED_TAG))
    }

    fn is_generated_blast_hit_feature(feature: &gb_io::seq::Feature) -> bool {
        feature
            .qualifier_values("gentle_generated".into())
            .any(|v| v.eq_ignore_ascii_case(BLAST_HIT_TRACK_GENERATED_TAG))
    }

    fn is_generated_genome_signal_feature(feature: &gb_io::seq::Feature) -> bool {
        Self::is_generated_genome_bed_feature(feature)
            || Self::is_generated_genome_bigwig_feature(feature)
            || Self::is_generated_genome_vcf_feature(feature)
    }

    fn remove_generated_genome_signal_features(features: &mut Vec<gb_io::seq::Feature>) {
        features.retain(|f| !Self::is_generated_genome_signal_feature(f));
    }

    fn remove_generated_blast_hit_features(features: &mut Vec<gb_io::seq::Feature>) {
        features.retain(|f| !Self::is_generated_blast_hit_feature(f));
    }

    fn open_text_reader(path: &str) -> Result<Box<dyn BufRead>, EngineError> {
        let file = std::fs::File::open(path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not open track file '{path}': {e}"),
        })?;
        let lower = path.to_ascii_lowercase();
        if lower.ends_with(".gz") {
            let decoder = GzDecoder::new(BufReader::new(file));
            Ok(Box::new(BufReader::new(decoder)))
        } else {
            Ok(Box::new(BufReader::new(file)))
        }
    }

    fn default_track_name(path: &str) -> String {
        let file_name = Path::new(path)
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("track");
        let without_gz = file_name
            .strip_suffix(".gz")
            .or_else(|| file_name.strip_suffix(".GZ"))
            .unwrap_or(file_name);
        let mut base = without_gz.to_string();
        let lower = base.to_ascii_lowercase();
        for suffix in [".bed", ".bigwig", ".bw", ".vcf"] {
            if lower.ends_with(suffix) && base.len() > suffix.len() {
                base.truncate(base.len() - suffix.len());
                break;
            }
        }
        if base.trim().is_empty() {
            "track".to_string()
        } else {
            base
        }
    }

    fn parse_bed_record(line: &str) -> Result<BedRecord, String> {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 3 {
            return Err("BED record needs at least 3 columns (chrom, start, end)".to_string());
        }
        let chromosome = fields[0].trim();
        if chromosome.is_empty() {
            return Err("BED chromosome is empty".to_string());
        }
        let start_0based = fields[1]
            .parse::<usize>()
            .map_err(|e| format!("Invalid BED start '{}': {e}", fields[1]))?;
        let end_0based = fields[2]
            .parse::<usize>()
            .map_err(|e| format!("Invalid BED end '{}': {e}", fields[2]))?;
        if end_0based <= start_0based {
            return Err(format!(
                "Invalid BED interval {}:{}-{} (end must be > start)",
                chromosome, start_0based, end_0based
            ));
        }
        let name = fields
            .get(3)
            .map(|v| v.trim())
            .filter(|v| !v.is_empty() && *v != ".")
            .map(|v| v.to_string());
        let score = fields
            .get(4)
            .map(|v| v.trim())
            .filter(|v| !v.is_empty() && *v != ".")
            .map(|v| {
                v.parse::<f64>()
                    .map_err(|e| format!("Invalid BED score '{}': {e}", v))
            })
            .transpose()?;
        let strand = fields
            .get(5)
            .and_then(|v| v.trim().chars().next())
            .filter(|c| matches!(c, '+' | '-'));

        Ok(BedRecord {
            chromosome: chromosome.to_string(),
            start_0based,
            end_0based,
            name,
            score,
            strand,
        })
    }

    fn parse_bedgraph_record(line: &str) -> Result<BedRecord, String> {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 4 {
            return Err(
                "bedGraph record needs at least 4 columns (chrom, start, end, value)".to_string(),
            );
        }
        let chromosome = fields[0].trim();
        if chromosome.is_empty() {
            return Err("bedGraph chromosome is empty".to_string());
        }
        let start_0based = fields[1]
            .parse::<usize>()
            .map_err(|e| format!("Invalid bedGraph start '{}': {e}", fields[1]))?;
        let end_0based = fields[2]
            .parse::<usize>()
            .map_err(|e| format!("Invalid bedGraph end '{}': {e}", fields[2]))?;
        if end_0based <= start_0based {
            return Err(format!(
                "Invalid bedGraph interval {}:{}-{} (end must be > start)",
                chromosome, start_0based, end_0based
            ));
        }
        let value = fields[3].parse::<f64>().map_err(|e| {
            format!(
                "Invalid bedGraph value '{}' (expected numeric): {e}",
                fields[3]
            )
        })?;
        Ok(BedRecord {
            chromosome: chromosome.to_string(),
            start_0based,
            end_0based,
            name: None,
            score: Some(value),
            strand: None,
        })
    }

    fn parse_vcf_record(line: &str) -> Result<VcfRecord, String> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 8 {
            return Err(
                "VCF record needs at least 8 columns (CHROM POS ID REF ALT QUAL FILTER INFO)"
                    .to_string(),
            );
        }
        let chromosome = fields[0].trim();
        if chromosome.is_empty() {
            return Err("VCF chromosome is empty".to_string());
        }
        let pos_1based = fields[1]
            .trim()
            .parse::<usize>()
            .map_err(|e| format!("Invalid VCF POS '{}': {e}", fields[1].trim()))?;
        if pos_1based == 0 {
            return Err("VCF POS must be >= 1".to_string());
        }
        let id = fields
            .get(2)
            .map(|v| v.trim())
            .filter(|v| !v.is_empty() && *v != ".")
            .map(|v| v.to_string());
        let reference = fields
            .get(3)
            .map(|v| v.trim())
            .filter(|v| !v.is_empty() && *v != ".")
            .ok_or_else(|| "VCF REF must be non-empty".to_string())?
            .to_string();
        let alternates = fields
            .get(4)
            .map(|v| v.trim())
            .ok_or_else(|| "VCF ALT must be present".to_string())?
            .split(',')
            .map(str::trim)
            .filter(|v| !v.is_empty() && *v != ".")
            .map(|v| v.to_string())
            .collect::<Vec<_>>();
        if alternates.is_empty() {
            return Err("VCF ALT must contain at least one allele".to_string());
        }
        let qual = fields
            .get(5)
            .map(|v| v.trim())
            .filter(|v| !v.is_empty() && *v != ".")
            .map(|v| {
                v.parse::<f64>()
                    .map_err(|e| format!("Invalid VCF QUAL '{}': {e}", v))
            })
            .transpose()?;
        let filter = fields
            .get(6)
            .map(|v| v.trim())
            .filter(|v| !v.is_empty() && *v != ".")
            .map(|v| v.to_string());
        let info = fields
            .get(7)
            .map(|v| v.trim())
            .filter(|v| !v.is_empty() && *v != ".")
            .map(|v| v.to_string());

        Ok(VcfRecord {
            chromosome: chromosome.to_string(),
            pos_1based,
            id,
            reference,
            alternates,
            qual,
            filter,
            info,
        })
    }

    fn build_genome_signal_feature(
        record: &BedRecord,
        track_name: &str,
        path: &str,
        local_start_0based: usize,
        local_end_0based_exclusive: usize,
        local_strand: Option<char>,
        generated_tag: &str,
        source_label: &str,
    ) -> gb_io::seq::Feature {
        let label = record
            .name
            .as_ref()
            .filter(|v| !v.trim().is_empty())
            .cloned()
            .unwrap_or_else(|| {
                format!(
                    "{}:{}:{}-{}",
                    track_name, record.chromosome, record.start_0based, record.end_0based
                )
            });
        let mut qualifiers = vec![
            ("label".into(), Some(label.clone())),
            (
                "note".into(),
                Some(format!(
                    "{} track '{}' from '{}' [{}:{}-{}]",
                    source_label,
                    track_name,
                    path,
                    record.chromosome,
                    record.start_0based,
                    record.end_0based
                )),
            ),
            ("gentle_track_source".into(), Some(source_label.to_string())),
            ("gentle_generated".into(), Some(generated_tag.to_string())),
            ("gentle_track_name".into(), Some(track_name.to_string())),
            ("gentle_track_file".into(), Some(path.to_string())),
            ("chromosome".into(), Some(record.chromosome.clone())),
            (
                "bed_start_0based".into(),
                Some(record.start_0based.to_string()),
            ),
            ("bed_end_0based".into(), Some(record.end_0based.to_string())),
        ];
        if let Some(score) = record.score {
            qualifiers.push(("score".into(), Some(format!("{score:.6}"))));
        }
        if let Some(strand) = record.strand {
            qualifiers.push(("bed_strand".into(), Some(strand.to_string())));
        }
        if let Some(strand) = local_strand {
            qualifiers.push(("strand".into(), Some(strand.to_string())));
        }

        let base_location = gb_io::seq::Location::simple_range(
            local_start_0based as i64,
            local_end_0based_exclusive as i64,
        );
        let location = if local_strand == Some('-') {
            gb_io::seq::Location::Complement(Box::new(base_location))
        } else {
            base_location
        };

        gb_io::seq::Feature {
            kind: gb_io::FeatureKind::from("track"),
            location,
            qualifiers,
        }
    }

    fn build_genome_bed_feature(
        record: &BedRecord,
        track_name: &str,
        path: &str,
        local_start_0based: usize,
        local_end_0based_exclusive: usize,
        local_strand: Option<char>,
    ) -> gb_io::seq::Feature {
        Self::build_genome_signal_feature(
            record,
            track_name,
            path,
            local_start_0based,
            local_end_0based_exclusive,
            local_strand,
            GENOME_BED_TRACK_GENERATED_TAG,
            "BED",
        )
    }

    fn build_genome_vcf_feature(
        record: &VcfRecord,
        alt: &str,
        track_name: &str,
        path: &str,
        local_start_0based: usize,
        local_end_0based_exclusive: usize,
    ) -> gb_io::seq::Feature {
        let label = record.id.clone().unwrap_or_else(|| {
            format!(
                "{}:{}:{} {}>{}",
                track_name, record.chromosome, record.pos_1based, record.reference, alt
            )
        });
        let mut qualifiers = vec![
            ("label".into(), Some(label.clone())),
            (
                "note".into(),
                Some(format!(
                    "VCF track '{}' from '{}' [{}:{} {}>{}]",
                    track_name, path, record.chromosome, record.pos_1based, record.reference, alt
                )),
            ),
            ("gentle_track_source".into(), Some("VCF".to_string())),
            (
                "gentle_generated".into(),
                Some(GENOME_VCF_TRACK_GENERATED_TAG.to_string()),
            ),
            ("gentle_track_name".into(), Some(track_name.to_string())),
            ("gentle_track_file".into(), Some(path.to_string())),
            ("chromosome".into(), Some(record.chromosome.clone())),
            ("vcf_pos_1based".into(), Some(record.pos_1based.to_string())),
            ("vcf_ref".into(), Some(record.reference.clone())),
            ("vcf_alt".into(), Some(alt.to_string())),
        ];
        if let Some(id) = &record.id {
            qualifiers.push(("vcf_id".into(), Some(id.clone())));
        }
        if let Some(qual) = record.qual {
            qualifiers.push(("score".into(), Some(format!("{qual:.6}"))));
            qualifiers.push(("vcf_qual".into(), Some(format!("{qual:.6}"))));
        }
        if let Some(filter) = &record.filter {
            qualifiers.push(("vcf_filter".into(), Some(filter.clone())));
        }
        if let Some(info) = &record.info {
            qualifiers.push(("vcf_info".into(), Some(info.clone())));
        }

        gb_io::seq::Feature {
            kind: gb_io::FeatureKind::from("track"),
            location: gb_io::seq::Location::simple_range(
                local_start_0based as i64,
                local_end_0based_exclusive as i64,
            ),
            qualifiers,
        }
    }

    fn build_blast_hit_feature(
        hit: &BlastHitFeatureInput,
        track_name: &str,
        local_start_0based: usize,
        local_end_0based_exclusive: usize,
        local_strand: Option<char>,
    ) -> gb_io::seq::Feature {
        let label = format!(
            "{}:{}:{}-{}",
            track_name, hit.subject_id, hit.query_start_1based, hit.query_end_1based
        );
        let mut qualifiers = vec![
            ("label".into(), Some(label)),
            (
                "note".into(),
                Some(format!(
                    "BLAST hit '{}' query={}..{} subject={}..{} identity={:.2}% bitscore={:.2} evalue={:.3e}",
                    hit.subject_id,
                    hit.query_start_1based,
                    hit.query_end_1based,
                    hit.subject_start_1based,
                    hit.subject_end_1based,
                    hit.identity_percent,
                    hit.bit_score,
                    hit.evalue
                )),
            ),
            ("gentle_track_source".into(), Some("BLAST".to_string())),
            (
                "gentle_generated".into(),
                Some(BLAST_HIT_TRACK_GENERATED_TAG.to_string()),
            ),
            ("gentle_track_name".into(), Some(track_name.to_string())),
            ("blast_subject_id".into(), Some(hit.subject_id.clone())),
            (
                "blast_query_start_1based".into(),
                Some(hit.query_start_1based.to_string()),
            ),
            (
                "blast_query_end_1based".into(),
                Some(hit.query_end_1based.to_string()),
            ),
            (
                "blast_subject_start_1based".into(),
                Some(hit.subject_start_1based.to_string()),
            ),
            (
                "blast_subject_end_1based".into(),
                Some(hit.subject_end_1based.to_string()),
            ),
            (
                "blast_identity_percent".into(),
                Some(format!("{:.6}", hit.identity_percent)),
            ),
            ("score".into(), Some(format!("{:.6}", hit.bit_score))),
            (
                "blast_bit_score".into(),
                Some(format!("{:.6}", hit.bit_score)),
            ),
            ("blast_evalue".into(), Some(format!("{:.3e}", hit.evalue))),
        ];
        if let Some(qcov) = hit.query_coverage_percent {
            qualifiers.push(("blast_qcov_percent".into(), Some(format!("{qcov:.4}"))));
        }
        if let Some(strand) = local_strand {
            qualifiers.push(("strand".into(), Some(strand.to_string())));
        }

        let base_location = gb_io::seq::Location::simple_range(
            local_start_0based as i64,
            local_end_0based_exclusive as i64,
        );
        let location = if local_strand == Some('-') {
            gb_io::seq::Location::Complement(Box::new(base_location))
        } else {
            base_location
        };

        gb_io::seq::Feature {
            kind: gb_io::FeatureKind::from("track"),
            location,
            qualifiers,
        }
    }

    fn import_genome_bed_track(
        dna: &mut DNAsequence,
        anchor: &GenomeSequenceAnchor,
        path: &str,
        track_name: Option<&str>,
        min_score: Option<f64>,
        max_score: Option<f64>,
        clear_existing: bool,
    ) -> Result<GenomeBedTrackImportReport, EngineError> {
        let selected_track_name = track_name
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .map(|v| v.to_string())
            .unwrap_or_else(|| Self::default_track_name(path));

        if clear_existing {
            Self::remove_generated_genome_signal_features(dna.features_mut());
        }

        let mut report = GenomeBedTrackImportReport {
            track_name: selected_track_name.clone(),
            ..Default::default()
        };
        let mut reader = Self::open_text_reader(path)?;
        let mut line = String::new();
        let mut line_no = 0usize;

        while {
            line.clear();
            reader.read_line(&mut line).map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not read BED file '{path}': {e}"),
            })? > 0
        } {
            line_no += 1;
            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }
            if trimmed.starts_with('#')
                || trimmed.to_ascii_lowercase().starts_with("track ")
                || trimmed.to_ascii_lowercase().starts_with("browser ")
            {
                continue;
            }
            report.parsed_records += 1;

            let record = match Self::parse_bed_record(trimmed) {
                Ok(v) => v,
                Err(e) => {
                    report.skipped_records += 1;
                    report.skipped_invalid += 1;
                    if report.warnings.len() < 20 {
                        report
                            .warnings
                            .push(format!("BED line {} skipped: {}", line_no, e));
                    }
                    continue;
                }
            };

            if !Self::chromosomes_match(&record.chromosome, &anchor.chromosome) {
                report.skipped_records += 1;
                report.skipped_wrong_chromosome += 1;
                continue;
            }

            let bed_start_1based = record.start_0based.saturating_add(1);
            let bed_end_1based = record.end_0based;
            if bed_end_1based < anchor.start_1based || bed_start_1based > anchor.end_1based {
                report.skipped_records += 1;
                report.skipped_non_overlap += 1;
                continue;
            }
            let overlap_start_1based = bed_start_1based.max(anchor.start_1based);
            let overlap_end_1based = bed_end_1based.min(anchor.end_1based);
            if overlap_end_1based < overlap_start_1based {
                report.skipped_records += 1;
                report.skipped_non_overlap += 1;
                continue;
            }

            if min_score.is_some() || max_score.is_some() {
                let Some(score) = record.score else {
                    report.skipped_records += 1;
                    report.skipped_missing_score += 1;
                    continue;
                };
                if min_score.map(|v| score < v).unwrap_or(false)
                    || max_score.map(|v| score > v).unwrap_or(false)
                {
                    report.skipped_records += 1;
                    report.skipped_outside_score_range += 1;
                    continue;
                }
            }

            let (local_start_0based, local_end_0based_exclusive) = if anchor.strand == Some('-') {
                (
                    anchor.end_1based.saturating_sub(overlap_end_1based),
                    anchor.end_1based.saturating_sub(overlap_start_1based) + 1,
                )
            } else {
                (
                    overlap_start_1based - anchor.start_1based,
                    overlap_end_1based - anchor.start_1based + 1,
                )
            };
            let local_strand = match (record.strand, anchor.strand) {
                (Some('+'), Some('-')) => Some('-'),
                (Some('-'), Some('-')) => Some('+'),
                (Some(strand), _) => Some(strand),
                _ => None,
            };
            if report.imported_features >= MAX_IMPORTED_SIGNAL_FEATURES {
                report.truncated_at_limit = true;
                break;
            }
            let feature = Self::build_genome_bed_feature(
                &record,
                &selected_track_name,
                path,
                local_start_0based,
                local_end_0based_exclusive,
                local_strand,
            );
            dna.features_mut().push(feature);
            report.imported_features += 1;
        }

        Ok(report)
    }

    fn resolve_bigwig_to_bedgraph_executable() -> String {
        env::var(BIGWIG_TO_BEDGRAPH_ENV_BIN)
            .ok()
            .map(|v| v.trim().to_string())
            .filter(|v| !v.is_empty())
            .unwrap_or_else(|| DEFAULT_BIGWIG_TO_BEDGRAPH_BIN.to_string())
    }

    fn convert_bigwig_to_bedgraph(path: &str) -> Result<NamedTempFile, EngineError> {
        let executable = Self::resolve_bigwig_to_bedgraph_executable();
        let output = NamedTempFile::new().map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not create temporary bedGraph file: {e}"),
        })?;
        let output_path = output.path().to_path_buf();
        let command_output = Command::new(&executable)
            .arg(path)
            .arg(&output_path)
            .output()
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not execute '{}' for BigWig conversion (set {} to override): {}",
                    executable, BIGWIG_TO_BEDGRAPH_ENV_BIN, e
                ),
            })?;
        if !command_output.status.success() {
            let stderr = String::from_utf8_lossy(&command_output.stderr)
                .trim()
                .to_string();
            let stdout = String::from_utf8_lossy(&command_output.stdout)
                .trim()
                .to_string();
            let detail = if !stderr.is_empty() {
                stderr
            } else if !stdout.is_empty() {
                stdout
            } else {
                format!("exit status {}", command_output.status)
            };
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "BigWig conversion failed for '{}' via '{}' (set {} to override): {}",
                    path, executable, BIGWIG_TO_BEDGRAPH_ENV_BIN, detail
                ),
            });
        }
        Ok(output)
    }

    fn import_genome_bigwig_track(
        dna: &mut DNAsequence,
        anchor: &GenomeSequenceAnchor,
        path: &str,
        track_name: Option<&str>,
        min_score: Option<f64>,
        max_score: Option<f64>,
        clear_existing: bool,
    ) -> Result<GenomeBedTrackImportReport, EngineError> {
        let selected_track_name = track_name
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .map(|v| v.to_string())
            .unwrap_or_else(|| Self::default_track_name(path));

        if clear_existing {
            Self::remove_generated_genome_signal_features(dna.features_mut());
        }

        let mut report = GenomeBedTrackImportReport {
            track_name: selected_track_name.clone(),
            ..Default::default()
        };
        let bedgraph_file = Self::convert_bigwig_to_bedgraph(path)?;
        let bedgraph_path = bedgraph_file.path().to_string_lossy().to_string();
        let mut reader = Self::open_text_reader(&bedgraph_path)?;
        let mut line = String::new();
        let mut line_no = 0usize;

        while {
            line.clear();
            reader.read_line(&mut line).map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not read converted bedGraph for '{path}': {e}"),
            })? > 0
        } {
            line_no += 1;
            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }
            if trimmed.starts_with('#')
                || trimmed.to_ascii_lowercase().starts_with("track ")
                || trimmed.to_ascii_lowercase().starts_with("browser ")
            {
                continue;
            }
            report.parsed_records += 1;

            let record = match Self::parse_bedgraph_record(trimmed) {
                Ok(v) => v,
                Err(e) => {
                    report.skipped_records += 1;
                    report.skipped_invalid += 1;
                    if report.warnings.len() < 20 {
                        report
                            .warnings
                            .push(format!("bedGraph line {} skipped: {}", line_no, e));
                    }
                    continue;
                }
            };

            if !Self::chromosomes_match(&record.chromosome, &anchor.chromosome) {
                report.skipped_records += 1;
                report.skipped_wrong_chromosome += 1;
                continue;
            }

            let bed_start_1based = record.start_0based.saturating_add(1);
            let bed_end_1based = record.end_0based;
            if bed_end_1based < anchor.start_1based || bed_start_1based > anchor.end_1based {
                report.skipped_records += 1;
                report.skipped_non_overlap += 1;
                continue;
            }
            let overlap_start_1based = bed_start_1based.max(anchor.start_1based);
            let overlap_end_1based = bed_end_1based.min(anchor.end_1based);
            if overlap_end_1based < overlap_start_1based {
                report.skipped_records += 1;
                report.skipped_non_overlap += 1;
                continue;
            }

            if min_score.is_some() || max_score.is_some() {
                let Some(score) = record.score else {
                    report.skipped_records += 1;
                    report.skipped_missing_score += 1;
                    continue;
                };
                if min_score.map(|v| score < v).unwrap_or(false)
                    || max_score.map(|v| score > v).unwrap_or(false)
                {
                    report.skipped_records += 1;
                    report.skipped_outside_score_range += 1;
                    continue;
                }
            }

            let (local_start_0based, local_end_0based_exclusive) = if anchor.strand == Some('-') {
                (
                    anchor.end_1based.saturating_sub(overlap_end_1based),
                    anchor.end_1based.saturating_sub(overlap_start_1based) + 1,
                )
            } else {
                (
                    overlap_start_1based - anchor.start_1based,
                    overlap_end_1based - anchor.start_1based + 1,
                )
            };
            if report.imported_features >= MAX_IMPORTED_SIGNAL_FEATURES {
                report.truncated_at_limit = true;
                break;
            }
            let feature = Self::build_genome_signal_feature(
                &record,
                &selected_track_name,
                path,
                local_start_0based,
                local_end_0based_exclusive,
                None,
                GENOME_BIGWIG_TRACK_GENERATED_TAG,
                "BigWig",
            );
            dna.features_mut().push(feature);
            report.imported_features += 1;
        }

        Ok(report)
    }

    fn import_genome_vcf_track(
        dna: &mut DNAsequence,
        anchor: &GenomeSequenceAnchor,
        path: &str,
        track_name: Option<&str>,
        min_score: Option<f64>,
        max_score: Option<f64>,
        clear_existing: bool,
    ) -> Result<GenomeBedTrackImportReport, EngineError> {
        let selected_track_name = track_name
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .map(|v| v.to_string())
            .unwrap_or_else(|| Self::default_track_name(path));

        if clear_existing {
            Self::remove_generated_genome_signal_features(dna.features_mut());
        }

        let mut report = GenomeBedTrackImportReport {
            track_name: selected_track_name.clone(),
            ..Default::default()
        };
        let mut reader = Self::open_text_reader(path)?;
        let mut line = String::new();
        let mut line_no = 0usize;

        while {
            line.clear();
            reader.read_line(&mut line).map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not read VCF file '{path}': {e}"),
            })? > 0
        } {
            line_no += 1;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            report.parsed_records += 1;

            let record = match Self::parse_vcf_record(trimmed) {
                Ok(v) => v,
                Err(e) => {
                    report.skipped_records += 1;
                    report.skipped_invalid += 1;
                    if report.warnings.len() < 20 {
                        report
                            .warnings
                            .push(format!("VCF line {} skipped: {}", line_no, e));
                    }
                    continue;
                }
            };

            if !Self::chromosomes_match(&record.chromosome, &anchor.chromosome) {
                report.skipped_records += 1;
                report.skipped_wrong_chromosome += 1;
                continue;
            }

            let ref_len = record.reference.len().max(1);
            let variant_start_1based = record.pos_1based;
            let variant_end_1based = variant_start_1based.saturating_add(ref_len.saturating_sub(1));
            if variant_end_1based < anchor.start_1based || variant_start_1based > anchor.end_1based
            {
                report.skipped_records += 1;
                report.skipped_non_overlap += 1;
                continue;
            }

            if min_score.is_some() || max_score.is_some() {
                let Some(score) = record.qual else {
                    report.skipped_records += 1;
                    report.skipped_missing_score += 1;
                    continue;
                };
                if min_score.map(|v| score < v).unwrap_or(false)
                    || max_score.map(|v| score > v).unwrap_or(false)
                {
                    report.skipped_records += 1;
                    report.skipped_outside_score_range += 1;
                    continue;
                }
            }

            let overlap_start_1based = variant_start_1based.max(anchor.start_1based);
            let overlap_end_1based = variant_end_1based.min(anchor.end_1based);
            if overlap_end_1based < overlap_start_1based {
                report.skipped_records += 1;
                report.skipped_non_overlap += 1;
                continue;
            }

            let (local_start_0based, local_end_0based_exclusive) = if anchor.strand == Some('-') {
                (
                    anchor.end_1based.saturating_sub(overlap_end_1based),
                    anchor.end_1based.saturating_sub(overlap_start_1based) + 1,
                )
            } else {
                (
                    overlap_start_1based - anchor.start_1based,
                    overlap_end_1based - anchor.start_1based + 1,
                )
            };

            let mut stop = false;
            for alt in &record.alternates {
                if report.imported_features >= MAX_IMPORTED_SIGNAL_FEATURES {
                    report.truncated_at_limit = true;
                    stop = true;
                    break;
                }
                let feature = Self::build_genome_vcf_feature(
                    &record,
                    alt,
                    &selected_track_name,
                    path,
                    local_start_0based,
                    local_end_0based_exclusive,
                );
                dna.features_mut().push(feature);
                report.imported_features += 1;
            }
            if stop {
                break;
            }
        }

        Ok(report)
    }

    fn is_genbank_like_path(path: &str) -> bool {
        let lower = path.to_ascii_lowercase();
        lower.ends_with(".gb")
            || lower.ends_with(".gbk")
            || lower.ends_with(".genbank")
            || lower.ends_with(".gbff")
    }

    fn parse_first_usize_tokens(raw: &str, max_items: usize) -> Vec<usize> {
        if max_items == 0 {
            return vec![];
        }
        let mut out = Vec::with_capacity(max_items);
        let mut current = String::new();
        for ch in raw.chars() {
            if ch.is_ascii_digit() || ch == ',' {
                current.push(ch);
                continue;
            }
            if !current.is_empty() {
                if let Ok(value) = current.replace(',', "").parse::<usize>() {
                    out.push(value);
                    if out.len() >= max_items {
                        return out;
                    }
                }
                current.clear();
            }
        }
        if !current.is_empty() && out.len() < max_items {
            if let Ok(value) = current.replace(',', "").parse::<usize>() {
                out.push(value);
            }
        }
        out
    }

    fn parse_genbank_header_field_block(text: &str, key: &str) -> Option<String> {
        let mut block = String::new();
        let mut in_block = false;
        for line in text.lines() {
            if line.starts_with("ORIGIN") || line.starts_with("//") {
                break;
            }
            if let Some(rest) = line.strip_prefix(key) {
                in_block = true;
                let trimmed = rest.trim();
                if !trimmed.is_empty() {
                    block.push_str(trimmed);
                }
                continue;
            }
            if in_block {
                if line.starts_with("            ") {
                    let trimmed = line.trim();
                    if !trimmed.is_empty() {
                        if !block.is_empty() {
                            block.push(' ');
                        }
                        block.push_str(trimmed);
                    }
                    continue;
                }
                break;
            }
        }
        if block.is_empty() {
            None
        } else {
            Some(block)
        }
    }

    fn parse_genbank_accession_region(path: &str) -> Option<(String, usize, usize, String, char)> {
        if !Self::is_genbank_like_path(path) {
            return None;
        }
        let text = std::fs::read_to_string(path).ok()?;
        let definition =
            Self::parse_genbank_header_field_block(&text, "DEFINITION").unwrap_or_default();
        let accession_block = Self::parse_genbank_header_field_block(&text, "ACCESSION")?;
        if accession_block.is_empty() {
            return None;
        }
        let accession = accession_block.split_whitespace().next()?.to_string();
        let lower_block = accession_block.to_ascii_lowercase();
        let region_pos = lower_block.find("region:")?;
        let region_spec = &accession_block[region_pos + "region:".len()..];
        let lower_region_spec = region_spec.to_ascii_lowercase();
        let anchor_strand = if lower_region_spec.contains("complement(") {
            '-'
        } else {
            '+'
        };
        if lower_region_spec.contains("join(") || lower_region_spec.contains("order(") {
            return None;
        }
        let numbers = Self::parse_first_usize_tokens(region_spec, 2);
        if numbers.len() < 2 {
            return None;
        }
        let mut start_1based = numbers[0];
        let mut end_1based = numbers[1];
        if start_1based == 0 || end_1based == 0 {
            return None;
        }
        if end_1based < start_1based {
            std::mem::swap(&mut start_1based, &mut end_1based);
        }
        Some((
            accession,
            start_1based,
            end_1based,
            definition,
            anchor_strand,
        ))
    }

    fn parse_chromosome_from_definition(definition: &str) -> Option<String> {
        let lower = definition.to_ascii_lowercase();
        let marker = "chromosome ";
        let marker_pos = lower.find(marker)?;
        let raw_tail = &definition[marker_pos + marker.len()..];
        let token = raw_tail
            .chars()
            .take_while(|c| c.is_ascii_alphanumeric() || matches!(c, '_' | '-' | '.'))
            .collect::<String>();
        if token.is_empty() {
            None
        } else {
            Some(token)
        }
    }

    fn parse_genome_id_from_definition(definition: &str) -> Option<String> {
        let lower = definition.to_ascii_lowercase();
        let marker_pos = lower.find("primary assembly")?;
        let prefix = definition[..marker_pos]
            .trim()
            .trim_end_matches(|c: char| c.is_ascii_whitespace() || c == ',' || c == '.');
        let candidate = prefix.rsplit(',').next()?.trim();
        if candidate.is_empty() {
            None
        } else {
            Some(candidate.to_string())
        }
    }

    fn infer_imported_genbank_anchor(
        path: &str,
        dna: &DNAsequence,
    ) -> Option<GenomeSequenceAnchor> {
        let (accession, start_1based, end_1based, definition, anchor_strand) =
            Self::parse_genbank_accession_region(path)?;
        let region_len = end_1based - start_1based + 1;
        if dna.len() > 0 && dna.len() != region_len {
            return None;
        }
        let chromosome = Self::parse_chromosome_from_definition(&definition)
            .unwrap_or_else(|| accession.clone());
        let genome_id =
            Self::parse_genome_id_from_definition(&definition).unwrap_or_else(|| accession.clone());
        Some(GenomeSequenceAnchor {
            genome_id,
            chromosome,
            start_1based,
            end_1based,
            strand: Some(anchor_strand),
        })
    }

    fn classify_import_origin(path: &str, dna: &DNAsequence) -> SequenceOrigin {
        let lower = path.to_ascii_lowercase();
        if lower.ends_with(".fa")
            || lower.ends_with(".fasta")
            || lower.ends_with(".fna")
            || lower.ends_with(".ffn")
            || lower.ends_with(".faa")
        {
            // Requested policy: treat all FASTA imports as synthetic for now.
            return SequenceOrigin::ImportedSynthetic;
        }

        // For GenBank/EMBL-like records, use SOURCE/mol_type metadata if present.
        for feature in dna.features() {
            if feature.kind.to_string().to_ascii_uppercase() != "SOURCE" {
                continue;
            }
            for key in ["mol_type", "molecule_type"] {
                if let Some(value) = feature.qualifier_values(key.into()).next() {
                    let v = value.to_ascii_lowercase();
                    if v.contains("synthetic") {
                        return SequenceOrigin::ImportedSynthetic;
                    }
                    if v.contains("cdna") || v.contains("mrna") || v.contains("transcript") {
                        return SequenceOrigin::ImportedCdna;
                    }
                    if v.contains("genomic") {
                        return SequenceOrigin::ImportedGenomic;
                    }
                }
            }
        }

        SequenceOrigin::ImportedUnknown
    }

    fn add_lineage_node(
        &mut self,
        seq_id: &str,
        origin: SequenceOrigin,
        created_by_op: Option<&str>,
    ) -> NodeId {
        self.state.lineage.next_node_counter += 1;
        let node_id = format!("n-{}", self.state.lineage.next_node_counter);
        let node = LineageNode {
            node_id: node_id.clone(),
            seq_id: seq_id.to_string(),
            created_by_op: created_by_op.map(|s| s.to_string()),
            origin,
            created_at_unix_ms: Self::now_unix_ms(),
        };
        self.state
            .lineage
            .seq_to_node
            .insert(seq_id.to_string(), node_id.clone());
        self.state.lineage.nodes.insert(node_id.clone(), node);
        node_id
    }

    fn ensure_lineage_node(&mut self, seq_id: &str) -> NodeId {
        if let Some(node_id) = self.state.lineage.seq_to_node.get(seq_id) {
            return node_id.clone();
        }
        self.add_lineage_node(seq_id, SequenceOrigin::ImportedUnknown, None)
    }

    fn reconcile_lineage_nodes(&mut self) {
        let seq_ids: Vec<String> = self.state.sequences.keys().cloned().collect();
        for seq_id in seq_ids {
            let _ = self.ensure_lineage_node(&seq_id);
        }
    }

    fn next_container_id(&mut self) -> ContainerId {
        loop {
            self.state.container_state.next_container_counter += 1;
            let id = format!(
                "container-{}",
                self.state.container_state.next_container_counter
            );
            if !self.state.container_state.containers.contains_key(&id) {
                return id;
            }
        }
    }

    fn add_container(
        &mut self,
        members: &[SeqId],
        kind: ContainerKind,
        name: Option<String>,
        created_by_op: Option<&str>,
    ) -> Option<ContainerId> {
        if members.is_empty() {
            return None;
        }
        let container_id = self.next_container_id();
        let container = Container {
            container_id: container_id.clone(),
            kind,
            name,
            members: members.to_vec(),
            created_by_op: created_by_op.map(ToString::to_string),
            created_at_unix_ms: Self::now_unix_ms(),
        };
        self.state
            .container_state
            .containers
            .insert(container_id.clone(), container);
        for seq_id in members {
            self.state
                .container_state
                .seq_to_latest_container
                .insert(seq_id.clone(), container_id.clone());
        }
        Some(container_id)
    }

    fn reconcile_containers(&mut self) {
        let seq_ids: Vec<SeqId> = self.state.sequences.keys().cloned().collect();
        for seq_id in &seq_ids {
            if !self
                .state
                .container_state
                .seq_to_latest_container
                .contains_key(seq_id)
            {
                let _ = self.add_container(
                    std::slice::from_ref(seq_id),
                    ContainerKind::Singleton,
                    Some(format!("Imported sequence {seq_id}")),
                    None,
                );
            }
        }
    }

    fn add_container_from_result(&mut self, op: &Operation, result: &OpResult) {
        if result.created_seq_ids.is_empty() {
            return;
        }
        let kind = if matches!(op, Operation::SelectCandidate { .. }) {
            ContainerKind::Selection
        } else if result.created_seq_ids.len() > 1 {
            ContainerKind::Pool
        } else {
            ContainerKind::Singleton
        };
        let name = match op {
            Operation::LoadFile { .. } => Some("Imported sequence".to_string()),
            Operation::Digest { .. } => Some("Digest products".to_string()),
            Operation::DigestContainer { .. } => Some("Digest products".to_string()),
            Operation::MergeContainers { .. } => Some("Merged container".to_string()),
            Operation::MergeContainersById { .. } => Some("Merged container".to_string()),
            Operation::Ligation { .. } => Some("Ligation products".to_string()),
            Operation::LigationContainer { .. } => Some("Ligation products".to_string()),
            Operation::Pcr { .. }
            | Operation::PcrAdvanced { .. }
            | Operation::PcrMutagenesis { .. } => Some("PCR products".to_string()),
            Operation::ExtractRegion { .. } => Some("Extracted region".to_string()),
            Operation::ExtractAnchoredRegion { .. } => Some("Extracted region".to_string()),
            Operation::ExtractGenomeRegion { .. } => Some("Extracted genome region".to_string()),
            Operation::ExtractGenomeGene { .. } => Some("Extracted genome gene".to_string()),
            Operation::ImportGenomeBedTrack { .. } => Some("Imported BED track".to_string()),
            Operation::ImportGenomeBigWigTrack { .. } => Some("Imported BigWig track".to_string()),
            Operation::ImportGenomeVcfTrack { .. } => Some("Imported VCF track".to_string()),
            Operation::ImportBlastHitsTrack { .. } => Some("Imported BLAST hit track".to_string()),
            Operation::SelectCandidate { .. } => Some("Selected candidate".to_string()),
            Operation::FilterByMolecularWeight { .. } => {
                Some("Molecular-weight filtered".to_string())
            }
            Operation::FilterByDesignConstraints { .. } => {
                Some("Design constraints filtered".to_string())
            }
            Operation::FilterContainerByMolecularWeight { .. } => {
                Some("Molecular-weight filtered".to_string())
            }
            Operation::Reverse { .. }
            | Operation::Complement { .. }
            | Operation::ReverseComplement { .. }
            | Operation::Branch { .. } => Some("Derived sequence".to_string()),
            _ => None,
        };
        let _ = self.add_container(&result.created_seq_ids, kind, name, Some(&result.op_id));
    }

    fn add_lineage_edges(
        &mut self,
        parent_seq_ids: &[SeqId],
        created_seq_ids: &[SeqId],
        op_id: &str,
        run_id: &str,
    ) {
        if parent_seq_ids.is_empty() || created_seq_ids.is_empty() {
            return;
        }
        let parent_nodes: Vec<NodeId> = parent_seq_ids
            .iter()
            .map(|seq_id| self.ensure_lineage_node(seq_id))
            .collect();
        let child_nodes: Vec<NodeId> = created_seq_ids
            .iter()
            .map(|seq_id| self.ensure_lineage_node(seq_id))
            .collect();
        for from_node_id in &parent_nodes {
            for to_node_id in &child_nodes {
                self.state.lineage.edges.push(LineageEdge {
                    from_node_id: from_node_id.clone(),
                    to_node_id: to_node_id.clone(),
                    op_id: op_id.to_string(),
                    run_id: run_id.to_string(),
                });
            }
        }
    }

    fn unique_seq_id(&self, base: &str) -> SeqId {
        if !self.state.sequences.contains_key(base) {
            return base.to_string();
        }
        let mut i = 2usize;
        loop {
            let candidate = format!("{base}_{i}");
            if !self.state.sequences.contains_key(&candidate) {
                return candidate;
            }
            i += 1;
        }
    }

    fn prepare_sequence(dna: &mut DNAsequence) {
        Self::prepare_sequence_light(dna);
        dna.update_computed_features();
    }

    fn prepare_sequence_light(dna: &mut DNAsequence) {
        *dna.restriction_enzymes_mut() = active_restriction_enzymes();
        dna.set_max_restriction_enzyme_sites(Some(2));
        dna.set_methylation_mode(MethylationMode::both());
    }

    fn digest_with_guard(
        dna: &DNAsequence,
        enzymes: Vec<RestrictionEnzyme>,
        max_fragments: usize,
    ) -> Result<Vec<DNAsequence>, EngineError> {
        let mut fragments = vec![dna.clone()];
        for enzyme in &enzymes {
            println!("Digesting with enzyme: {}", enzyme.name);
            let mut seen_states: HashSet<u64> = HashSet::new();
            let mut rounds: usize = 0;
            let mut last_fragment_count = fragments.len();
            let enzyme_started = Instant::now();
            // Conservative guard against non-converging digest loops.
            let max_rounds = max_fragments.min(1_024).max(64);
            loop {
                rounds += 1;
                if enzyme_started.elapsed().as_millis() > 750 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Digest timed out for enzyme '{}' (>{} ms)",
                            enzyme.name, 750
                        ),
                    });
                }
                if rounds > max_rounds {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Digest did not converge for enzyme '{}' within {} rounds",
                            enzyme.name, max_rounds
                        ),
                    });
                }

                // Detect cyclic/non-converging fragment states.
                let mut state_hasher = DefaultHasher::new();
                fragments.len().hash(&mut state_hasher);
                for seq in &fragments {
                    seq.get_forward_string().hash(&mut state_hasher);
                    seq.overhang().forward_3.hash(&mut state_hasher);
                    seq.overhang().forward_5.hash(&mut state_hasher);
                    seq.overhang().reverse_3.hash(&mut state_hasher);
                    seq.overhang().reverse_5.hash(&mut state_hasher);
                }
                let state_sig = state_hasher.finish();
                if !seen_states.insert(state_sig) {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Digest entered a repeated state for enzyme '{}'",
                            enzyme.name
                        ),
                    });
                }

                let mut found_one = false;
                let mut new_fragments: Vec<DNAsequence> = vec![];
                for seq in fragments.drain(..) {
                    if let Some(site) = enzyme.get_sites(&seq, None).first() {
                        let split = seq.split_at_restriction_enzyme_site(site);
                        found_one = true;
                        new_fragments.extend(split);
                    } else {
                        new_fragments.push(seq);
                    }

                    if new_fragments.len() > max_fragments {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!(
                                "Digest produced more than max_fragments_per_container={}",
                                max_fragments
                            ),
                        });
                    }
                }
                fragments = new_fragments;
                let current_count = fragments.len();

                if !found_one {
                    println!(
                        "Digest enzyme '{}' completed in {} round(s), fragments: {}",
                        enzyme.name, rounds, current_count
                    );
                    break;
                }
                // For linear-digest progression, total fragment count should increase.
                // If not, we are likely cutting in a pathological cycle.
                if current_count <= last_fragment_count && rounds > 1 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Digest stalled for enzyme '{}' (no fragment-count progress)",
                            enzyme.name
                        ),
                    });
                }
                last_fragment_count = current_count;
            }
        }
        Ok(fragments)
    }

    fn max_fragments_per_container(&self) -> usize {
        self.state.parameters.max_fragments_per_container
    }

    fn container_members(&self, container_id: &str) -> Result<Vec<SeqId>, EngineError> {
        let container = self
            .state
            .container_state
            .containers
            .get(container_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Container '{container_id}' not found"),
            })?;
        if container.members.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Container '{container_id}' has no members"),
            });
        }
        for seq_id in &container.members {
            if !self.state.sequences.contains_key(seq_id) {
                return Err(EngineError {
                    code: ErrorCode::NotFound,
                    message: format!(
                        "Container '{container_id}' references unknown sequence '{seq_id}'"
                    ),
                });
            }
        }
        Ok(container.members.clone())
    }

    fn flatten_container_members(
        &self,
        container_ids: &[ContainerId],
    ) -> Result<Vec<SeqId>, EngineError> {
        if container_ids.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "At least one container id is required".to_string(),
            });
        }
        let mut members = Vec::new();
        for container_id in container_ids {
            members.extend(self.container_members(container_id)?);
            if members.len() > self.max_fragments_per_container() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Container merge input count exceeds max_fragments_per_container={}",
                        self.max_fragments_per_container()
                    ),
                });
            }
        }
        if members.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "No sequences found in selected containers".to_string(),
            });
        }
        Ok(members)
    }

    fn resolve_enzymes(
        &self,
        enzymes: &[String],
    ) -> Result<(Vec<RestrictionEnzyme>, Vec<String>), EngineError> {
        if enzymes.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Digest requires at least one enzyme".to_string(),
            });
        }

        let mut by_name: HashMap<String, RestrictionEnzyme> = HashMap::new();
        for enzyme in active_restriction_enzymes() {
            by_name.entry(enzyme.name.clone()).or_insert(enzyme);
        }
        for dna in self.state.sequences.values() {
            for enzyme in dna.restriction_enzymes() {
                by_name
                    .entry(enzyme.name.clone())
                    .or_insert_with(|| enzyme.clone());
            }
        }

        let mut found = Vec::new();
        let mut missing = Vec::new();
        let mut seen_names: HashSet<String> = HashSet::new();
        for name in enzymes {
            if let Some(enzyme) = by_name.get(name) {
                if seen_names.insert(name.clone()) {
                    found.push(enzyme.clone());
                }
            } else if !missing.contains(name) {
                missing.push(name.clone());
            }
        }

        if found.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "None of the requested enzymes are known: {}",
                    enzymes.join(",")
                ),
            });
        }
        Ok((found, missing))
    }

    pub fn summarize_state(&self) -> EngineStateSummary {
        let mut sequences: Vec<EngineSequenceSummary> = self
            .state
            .sequences
            .iter()
            .map(|(id, dna)| EngineSequenceSummary {
                id: id.to_string(),
                name: dna.name().clone(),
                length: dna.len(),
                circular: dna.is_circular(),
            })
            .collect();
        sequences.sort_by(|a, b| a.id.cmp(&b.id));

        let mut containers: Vec<EngineContainerSummary> = self
            .state
            .container_state
            .containers
            .iter()
            .map(|(id, c)| EngineContainerSummary {
                id: id.to_string(),
                kind: format!("{:?}", c.kind),
                member_count: c.members.len(),
                members: c.members.clone(),
            })
            .collect();
        containers.sort_by(|a, b| a.id.cmp(&b.id));

        EngineStateSummary {
            sequence_count: sequences.len(),
            sequences,
            container_count: containers.len(),
            containers,
            display: self.state.display.clone(),
        }
    }

    fn canonical_fasta_molecule(raw: Option<&str>) -> &'static str {
        let normalized = raw
            .unwrap_or("dsdna")
            .trim()
            .to_ascii_lowercase()
            .replace(['_', '-'], "")
            .replace(' ', "");
        match normalized.as_str() {
            "ssdna" | "singlestrandeddna" | "single" | "ssdnaoligo" | "ss" => "ssdna",
            "rna" | "ssrna" | "singlestrandedrna" | "transcript" | "mrna" | "cdna" => "rna",
            _ => "dsdna",
        }
    }

    fn fasta_metadata_tokens(dna: &DNAsequence) -> Vec<String> {
        let molecule = Self::canonical_fasta_molecule(dna.molecule_type());
        let mut tokens = vec![
            format!("molecule={molecule}"),
            format!(
                "topology={}",
                if dna.is_circular() {
                    "circular"
                } else {
                    "linear"
                }
            ),
        ];

        if molecule == "dsdna" {
            let overhang = dna.overhang();
            if !overhang.forward_5.is_empty() {
                tokens.push(format!("f5={}", Self::overhang_text(&overhang.forward_5)));
            }
            if !overhang.forward_3.is_empty() {
                tokens.push(format!("f3={}", Self::overhang_text(&overhang.forward_3)));
            }
            if !overhang.reverse_5.is_empty() {
                tokens.push(format!("r5={}", Self::overhang_text(&overhang.reverse_5)));
            }
            if !overhang.reverse_3.is_empty() {
                tokens.push(format!("r3={}", Self::overhang_text(&overhang.reverse_3)));
            }
        }

        tokens
    }

    fn save_as_fasta(seq_id: &str, dna: &DNAsequence, path: &str) -> Result<(), EngineError> {
        let mut file = File::create(path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not create FASTA file '{path}': {e}"),
        })?;

        let header = dna
            .name()
            .clone()
            .unwrap_or_else(|| seq_id.to_string())
            .replace(' ', "_");
        let seq = dna.get_forward_string();

        let mut header_line = format!(">{header}");
        let metadata = Self::fasta_metadata_tokens(dna);
        if !metadata.is_empty() {
            header_line.push(' ');
            header_line.push_str(&metadata.join(" "));
        }

        writeln!(file, "{header_line}").map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write FASTA header to '{path}': {e}"),
        })?;

        for chunk in seq.as_bytes().chunks(80) {
            file.write_all(chunk).map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not write FASTA sequence to '{path}': {e}"),
            })?;
            file.write_all(b"\n").map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not write FASTA newline to '{path}': {e}"),
            })?;
        }

        Ok(())
    }

    fn overhang_text(v: &[u8]) -> String {
        String::from_utf8_lossy(v).to_string()
    }

    fn infer_pool_end_type(end: &PoolEnd) -> String {
        let any = !end.forward_5.is_empty()
            || !end.forward_3.is_empty()
            || !end.reverse_5.is_empty()
            || !end.reverse_3.is_empty();
        if !any {
            return "blunt".to_string();
        }
        if !end.forward_5.is_empty() || !end.reverse_5.is_empty() {
            return "sticky_5p".to_string();
        }
        if !end.forward_3.is_empty() || !end.reverse_3.is_empty() {
            return "sticky_3p".to_string();
        }
        "unknown".to_string()
    }

    fn default_pool_human_id(inputs: &[SeqId]) -> String {
        if inputs.len() <= 4 {
            format!("Pool({})", inputs.join(", "))
        } else {
            format!("Pool({} members, first: {})", inputs.len(), inputs[0])
        }
    }

    fn export_pool_file(
        &self,
        inputs: &[SeqId],
        path: &str,
        pool_id: Option<String>,
        human_id: Option<String>,
    ) -> Result<(String, String, usize), EngineError> {
        if inputs.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "ExportPool requires at least one input sequence id".to_string(),
            });
        }
        let mut members: Vec<PoolMember> = Vec::with_capacity(inputs.len());
        for seq_id in inputs {
            let dna = self
                .state
                .sequences
                .get(seq_id)
                .ok_or_else(|| EngineError {
                    code: ErrorCode::NotFound,
                    message: format!("Sequence '{seq_id}' not found"),
                })?;
            let mut end = PoolEnd {
                end_type: String::new(),
                forward_5: Self::overhang_text(&dna.overhang().forward_5),
                forward_3: Self::overhang_text(&dna.overhang().forward_3),
                reverse_5: Self::overhang_text(&dna.overhang().reverse_5),
                reverse_3: Self::overhang_text(&dna.overhang().reverse_3),
            };
            end.end_type = Self::infer_pool_end_type(&end);
            members.push(PoolMember {
                seq_id: seq_id.clone(),
                human_id: dna.name().clone().unwrap_or_else(|| seq_id.clone()),
                name: dna.name().clone(),
                sequence: dna.get_forward_string(),
                length_bp: dna.len(),
                topology: if dna.is_circular() {
                    "circular".to_string()
                } else {
                    "linear".to_string()
                },
                ends: end,
            });
        }

        let pool_id = pool_id.unwrap_or_else(|| "pool_export".to_string());
        let human_id = human_id.unwrap_or_else(|| Self::default_pool_human_id(inputs));
        let export = PoolExport {
            schema: "gentle.pool.v1".to_string(),
            pool_id: pool_id.clone(),
            human_id: human_id.clone(),
            member_count: members.len(),
            members,
        };
        let text = serde_json::to_string_pretty(&export).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize pool JSON: {e}"),
        })?;
        std::fs::write(path, text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write pool file '{path}': {e}"),
        })?;
        Ok((pool_id, human_id, export.member_count))
    }

    fn reverse_complement(seq: &str) -> String {
        seq.as_bytes()
            .iter()
            .rev()
            .map(|c| IupacCode::letter_complement(*c))
            .map(char::from)
            .collect()
    }

    fn reverse_complement_bytes(seq: &[u8]) -> Vec<u8> {
        seq.iter()
            .rev()
            .map(|c| IupacCode::letter_complement(*c))
            .collect()
    }

    fn complement_bytes(seq: &[u8]) -> Vec<u8> {
        seq.iter()
            .map(|c| IupacCode::letter_complement(*c))
            .collect()
    }

    fn normalize_dna_text(seq: &str) -> String {
        let cleaned = DNAsequence::validate_dna_sequence(seq.as_bytes());
        String::from_utf8_lossy(&cleaned).to_string()
    }

    fn normalize_iupac_text(seq: &str) -> Result<String, EngineError> {
        let upper = seq.trim().to_ascii_uppercase();
        if upper.is_empty() {
            return Ok(upper);
        }
        for b in upper.as_bytes() {
            if !IupacCode::is_valid_letter(*b) {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("Invalid IUPAC nucleotide '{}'", *b as char),
                });
            }
        }
        Ok(upper)
    }

    fn resolve_tf_motif_or_iupac(token: &str) -> Result<String, EngineError> {
        let normalized = token.trim().to_ascii_uppercase();
        if normalized.is_empty() {
            return Ok(normalized);
        }
        if normalized
            .as_bytes()
            .iter()
            .all(|b| IupacCode::is_valid_letter(*b))
        {
            return Self::normalize_iupac_text(&normalized);
        }
        if let Some(motif) = tf_motifs::resolve_motif_definition(&normalized) {
            return Self::normalize_iupac_text(&motif.consensus_iupac);
        }
        Err(EngineError {
            code: ErrorCode::NotFound,
            message: format!(
                "TF motif '{}' was neither valid IUPAC text nor found in local motif registry",
                token
            ),
        })
    }

    fn iupac_letter_counts(letter: u8) -> [f64; 4] {
        match letter.to_ascii_uppercase() {
            b'A' => [1.0, 0.0, 0.0, 0.0],
            b'C' => [0.0, 1.0, 0.0, 0.0],
            b'G' => [0.0, 0.0, 1.0, 0.0],
            b'T' | b'U' => [0.0, 0.0, 0.0, 1.0],
            b'M' => [1.0, 1.0, 0.0, 0.0], // A/C
            b'R' => [1.0, 0.0, 1.0, 0.0], // A/G
            b'W' => [1.0, 0.0, 0.0, 1.0], // A/T
            b'S' => [0.0, 1.0, 1.0, 0.0], // C/G
            b'Y' => [0.0, 1.0, 0.0, 1.0], // C/T
            b'K' => [0.0, 0.0, 1.0, 1.0], // G/T
            b'V' => [1.0, 1.0, 1.0, 0.0], // A/C/G
            b'H' => [1.0, 1.0, 0.0, 1.0], // A/C/T
            b'D' => [1.0, 0.0, 1.0, 1.0], // A/G/T
            b'B' => [0.0, 1.0, 1.0, 1.0], // C/G/T
            _ => [1.0, 1.0, 1.0, 1.0],    // N/unknown
        }
    }

    fn matrix_from_iupac(consensus: &str) -> Vec<[f64; 4]> {
        consensus
            .as_bytes()
            .iter()
            .map(|b| Self::iupac_letter_counts(*b))
            .collect()
    }

    fn resolve_tf_motif_for_scoring(
        token: &str,
    ) -> Result<(String, Option<String>, String, Vec<[f64; 4]>), EngineError> {
        let normalized = token.trim().to_ascii_uppercase();
        if normalized.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Empty TF motif token".to_string(),
            });
        }
        if normalized
            .as_bytes()
            .iter()
            .all(|b| IupacCode::is_valid_letter(*b))
        {
            let consensus = Self::normalize_iupac_text(&normalized)?;
            let matrix = Self::matrix_from_iupac(&consensus);
            return Ok((consensus.clone(), None, consensus, matrix));
        }
        if let Some(motif) = tf_motifs::resolve_motif_definition(&normalized) {
            let consensus = Self::normalize_iupac_text(&motif.consensus_iupac)?;
            return Ok((motif.id, motif.name, consensus, motif.matrix_counts));
        }
        Err(EngineError {
            code: ErrorCode::NotFound,
            message: format!(
                "TF motif '{}' was neither valid IUPAC text nor found in local motif registry",
                token
            ),
        })
    }

    fn validate_tf_thresholds(min_llr_quantile: f64) -> Result<(), EngineError> {
        if !(0.0..=1.0).contains(&min_llr_quantile) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "TFBS min_llr_quantile must be between 0.0 and 1.0, got {}",
                    min_llr_quantile
                ),
            });
        }
        Ok(())
    }

    fn format_tf_threshold_summary(min_llr_bits: f64, min_llr_quantile: f64) -> String {
        let mut parts: Vec<String> = vec![];
        if min_llr_bits.is_finite() {
            parts.push(format!("min_llr_bits={min_llr_bits}"));
        }
        if min_llr_quantile > 0.0 {
            parts.push(format!("min_llr_quantile={min_llr_quantile}"));
        }
        if parts.is_empty() {
            String::new()
        } else {
            format!(" with {}", parts.join(", "))
        }
    }

    fn smooth_probability_matrix(matrix_counts: &[[f64; 4]]) -> Vec<[f64; 4]> {
        if matrix_counts.is_empty() {
            return vec![];
        }
        let max_col_sum = matrix_counts
            .iter()
            .map(|c| c.iter().sum::<f64>())
            .fold(0.0_f64, f64::max);
        let baseline = max_col_sum.max(1.0);

        let mut out = Vec::with_capacity(matrix_counts.len());
        for col in matrix_counts {
            let mut adjusted = *col;
            let col_sum = adjusted.iter().sum::<f64>();

            // Missing observations are distributed uniformly to match the
            // highest-supported column count in this motif.
            if baseline > col_sum {
                let add = (baseline - col_sum) / 4.0;
                for v in &mut adjusted {
                    *v += add;
                }
            }

            let epsilon = baseline * 1e-9;
            for v in &mut adjusted {
                *v += epsilon;
            }

            let total = adjusted.iter().sum::<f64>();
            let mut p_col = [0.0_f64; 4];
            for i in 0..4 {
                let p = (adjusted[i] / total).clamp(f64::MIN_POSITIVE, 1.0 - f64::EPSILON);
                p_col[i] = p;
            }
            out.push(p_col);
        }
        out
    }

    fn prepare_scoring_matrices(matrix_counts: &[[f64; 4]]) -> (Vec<[f64; 4]>, Vec<[f64; 4]>) {
        let probabilities = Self::smooth_probability_matrix(matrix_counts);
        let background = [0.25_f64, 0.25_f64, 0.25_f64, 0.25_f64];
        let mut llr = Vec::with_capacity(probabilities.len());
        let mut true_log_odds = Vec::with_capacity(probabilities.len());

        for col in probabilities {
            let mut llr_col = [0.0_f64; 4];
            let mut lor_col = [0.0_f64; 4];
            for i in 0..4 {
                let p = col[i];
                let q = background[i];
                llr_col[i] = (p / q).log2();
                let odds_p = p / (1.0 - p);
                let odds_q = q / (1.0 - q);
                lor_col[i] = (odds_p / odds_q).log2();
            }
            llr.push(llr_col);
            true_log_odds.push(lor_col);
        }
        (llr, true_log_odds)
    }

    fn base_to_idx(base: u8) -> Option<usize> {
        match base.to_ascii_uppercase() {
            b'A' => Some(0),
            b'C' => Some(1),
            b'G' => Some(2),
            b'T' => Some(3),
            _ => None,
        }
    }

    fn score_matrix_window(window: &[u8], score_matrix: &[[f64; 4]]) -> Option<f64> {
        if window.len() != score_matrix.len() {
            return None;
        }
        let mut score = 0.0_f64;
        for (idx, base) in window.iter().enumerate() {
            let b = Self::base_to_idx(*base)?;
            score += score_matrix[idx][b];
        }
        Some(score)
    }

    fn empirical_quantile(sorted_scores: &[f64], score: f64) -> f64 {
        if sorted_scores.is_empty() {
            return 0.0;
        }
        let mut lo = 0usize;
        let mut hi = sorted_scores.len();
        while lo < hi {
            let mid = (lo + hi) / 2;
            if sorted_scores[mid] <= score {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        lo as f64 / sorted_scores.len() as f64
    }

    fn scan_tf_scores(
        sequence: &[u8],
        llr_matrix: &[[f64; 4]],
        true_log_odds_matrix: &[[f64; 4]],
        mut on_progress: impl FnMut(usize, usize),
    ) -> Vec<(usize, bool, f64, f64, f64, f64)> {
        if llr_matrix.is_empty()
            || sequence.len() < llr_matrix.len()
            || llr_matrix.len() != true_log_odds_matrix.len()
        {
            return vec![];
        }
        let mut raw_hits = Vec::new();
        let mut all_llr_scores = Vec::new();
        let mut all_true_log_odds_scores = Vec::new();
        let len = llr_matrix.len();
        let windows = sequence.len().saturating_sub(len).saturating_add(1);
        let total_steps = windows.saturating_mul(2);
        let progress_stride = (total_steps / 200).max(1);
        let mut scanned_steps = 0usize;
        on_progress(scanned_steps, total_steps);
        for start in 0..=(sequence.len() - len) {
            let window = &sequence[start..start + len];
            if let (Some(llr), Some(true_log_odds)) = (
                Self::score_matrix_window(window, llr_matrix),
                Self::score_matrix_window(window, true_log_odds_matrix),
            ) {
                all_llr_scores.push(llr);
                all_true_log_odds_scores.push(true_log_odds);
                raw_hits.push((start, false, llr, true_log_odds));
            }
            scanned_steps = scanned_steps.saturating_add(1);
            if scanned_steps % progress_stride == 0 || scanned_steps == total_steps {
                on_progress(scanned_steps, total_steps);
            }
            let rc_window = Self::reverse_complement_bytes(window);
            if let (Some(llr), Some(true_log_odds)) = (
                Self::score_matrix_window(&rc_window, llr_matrix),
                Self::score_matrix_window(&rc_window, true_log_odds_matrix),
            ) {
                all_llr_scores.push(llr);
                all_true_log_odds_scores.push(true_log_odds);
                raw_hits.push((start, true, llr, true_log_odds));
            }
            scanned_steps = scanned_steps.saturating_add(1);
            if scanned_steps % progress_stride == 0 || scanned_steps == total_steps {
                on_progress(scanned_steps, total_steps);
            }
        }
        if scanned_steps != total_steps {
            on_progress(total_steps, total_steps);
        }
        all_llr_scores.sort_by(|a, b| a.total_cmp(b));
        all_true_log_odds_scores.sort_by(|a, b| a.total_cmp(b));
        raw_hits
            .into_iter()
            .map(|(start, reverse, llr_bits, true_log_odds_bits)| {
                (
                    start,
                    reverse,
                    llr_bits,
                    Self::empirical_quantile(&all_llr_scores, llr_bits),
                    true_log_odds_bits,
                    Self::empirical_quantile(&all_true_log_odds_scores, true_log_odds_bits),
                )
            })
            .collect()
    }

    fn is_generated_tfbs_feature(feature: &gb_io::seq::Feature) -> bool {
        feature
            .qualifier_values("gentle_generated".into())
            .any(|v| v.eq_ignore_ascii_case("tfbs"))
    }

    fn remove_generated_tfbs_features(features: &mut Vec<gb_io::seq::Feature>) {
        features.retain(|f| !Self::is_generated_tfbs_feature(f));
    }

    fn build_tfbs_feature(
        start: usize,
        end: usize,
        reverse: bool,
        motif_len: usize,
        tf_id: &str,
        tf_name: Option<&str>,
        llr_bits: f64,
        llr_quantile: f64,
        true_log_odds_bits: f64,
        true_log_odds_quantile: f64,
    ) -> gb_io::seq::Feature {
        let base_location = gb_io::seq::Location::simple_range(start as i64, end as i64);
        let location = if reverse {
            gb_io::seq::Location::Complement(Box::new(base_location))
        } else {
            base_location
        };
        let mut qualifiers = vec![
            ("label".into(), Some(format!("TFBS {tf_id}"))),
            ("tf_id".into(), Some(tf_id.to_string())),
            ("motif_length_bp".into(), Some(motif_len.to_string())),
            ("llr_bits".into(), Some(format!("{llr_bits:.6}"))),
            ("llr_quantile".into(), Some(format!("{llr_quantile:.6}"))),
            (
                "true_log_odds_bits".into(),
                Some(format!("{true_log_odds_bits:.6}")),
            ),
            (
                "true_log_odds_quantile".into(),
                Some(format!("{true_log_odds_quantile:.6}")),
            ),
            (
                "log_odds_ratio_bits".into(),
                Some(format!("{true_log_odds_bits:.6}")),
            ),
            (
                "log_odds_ratio_quantile".into(),
                Some(format!("{true_log_odds_quantile:.6}")),
            ),
            (
                "quantile_scope".into(),
                Some("per_motif_windows_both_strands".to_string()),
            ),
            (
                "note".into(),
                Some(format!(
                    "tf_id={tf_id}; motif_length_bp={motif_len}; llr_bits={llr_bits:.4}; llr_quantile={llr_quantile:.4}; true_log_odds_bits={true_log_odds_bits:.4}; true_log_odds_quantile={true_log_odds_quantile:.4}"
                )),
            ),
            ("gentle_generated".into(), Some("tfbs".to_string())),
        ];
        if let Some(name) = tf_name {
            if !name.trim().is_empty() {
                qualifiers.push(("bound_moiety".into(), Some(name.trim().to_string())));
            }
        }
        gb_io::seq::Feature {
            kind: gb_io::FeatureKind::from("TFBS"),
            location,
            qualifiers,
        }
    }

    fn iupac_letter_complement(letter: u8) -> Option<u8> {
        match letter.to_ascii_uppercase() {
            b'A' => Some(b'T'),
            b'C' => Some(b'G'),
            b'G' => Some(b'C'),
            b'T' | b'U' => Some(b'A'),
            b'W' => Some(b'W'),
            b'S' => Some(b'S'),
            b'M' => Some(b'K'),
            b'K' => Some(b'M'),
            b'R' => Some(b'Y'),
            b'Y' => Some(b'R'),
            b'B' => Some(b'V'),
            b'D' => Some(b'H'),
            b'H' => Some(b'D'),
            b'V' => Some(b'B'),
            b'N' => Some(b'N'),
            _ => None,
        }
    }

    fn reverse_complement_iupac(seq: &str) -> Result<String, EngineError> {
        let seq = Self::normalize_iupac_text(seq)?;
        let mut out = String::with_capacity(seq.len());
        for b in seq.as_bytes().iter().rev() {
            let c = Self::iupac_letter_complement(*b).ok_or_else(|| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Invalid IUPAC nucleotide '{}'", *b as char),
            })?;
            out.push(c as char);
        }
        Ok(out)
    }

    fn iupac_match_at(sequence: &[u8], pattern: &[u8], start: usize) -> bool {
        if pattern.is_empty() {
            return true;
        }
        if start > sequence.len() || start + pattern.len() > sequence.len() {
            return false;
        }
        for i in 0..pattern.len() {
            let s = IupacCode::from_letter(sequence[start + i]);
            let p = IupacCode::from_letter(pattern[i]);
            if s.is_empty() || p.is_empty() || s.subset(p).is_empty() {
                return false;
            }
        }
        true
    }

    fn contains_iupac_pattern(sequence: &[u8], pattern: &[u8]) -> bool {
        if pattern.is_empty() {
            return true;
        }
        if sequence.len() < pattern.len() {
            return false;
        }
        for start in 0..=(sequence.len() - pattern.len()) {
            if Self::iupac_match_at(sequence, pattern, start) {
                return true;
            }
        }
        false
    }

    fn contains_motif_any_strand(sequence: &[u8], motif: &str) -> Result<bool, EngineError> {
        let motif = Self::normalize_iupac_text(motif)?;
        if motif.is_empty() {
            return Ok(true);
        }
        let motif_bytes = motif.as_bytes();
        if Self::contains_iupac_pattern(sequence, motif_bytes) {
            return Ok(true);
        }
        let motif_rc = Self::reverse_complement_iupac(&motif)?;
        Ok(Self::contains_iupac_pattern(sequence, motif_rc.as_bytes()))
    }

    fn normalized_sequence_for_quality(dna: &DNAsequence) -> Vec<u8> {
        dna.get_forward_string()
            .as_bytes()
            .iter()
            .filter(|b| !b.is_ascii_whitespace())
            .map(|b| match b.to_ascii_uppercase() {
                b'U' => b'T',
                other => other,
            })
            .collect()
    }

    fn sequence_gc_fraction(sequence: &[u8]) -> Option<f64> {
        let mut canonical = 0usize;
        let mut gc = 0usize;
        for b in sequence {
            match b.to_ascii_uppercase() {
                b'G' | b'C' => {
                    canonical += 1;
                    gc += 1;
                }
                b'A' | b'T' => {
                    canonical += 1;
                }
                _ => {}
            }
        }
        if canonical == 0 {
            None
        } else {
            Some(gc as f64 / canonical as f64)
        }
    }

    fn max_homopolymer_run(sequence: &[u8]) -> usize {
        let mut best = 0usize;
        let mut current = 0usize;
        let mut prev = 0u8;
        for b in sequence {
            let base = b.to_ascii_uppercase();
            if !matches!(base, b'A' | b'C' | b'G' | b'T') {
                current = 0;
                prev = 0;
                continue;
            }
            if base == prev {
                current += 1;
            } else {
                current = 1;
                prev = base;
            }
            if current > best {
                best = current;
            }
        }
        best
    }

    fn contains_u6_terminator_t4(sequence: &[u8]) -> bool {
        sequence.windows(4).any(|w| w == b"TTTT")
    }

    fn has_ambiguous_bases(sequence: &[u8]) -> bool {
        sequence
            .iter()
            .any(|b| !matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
    }

    fn feature_labels(feature: &gb_io::seq::Feature) -> Vec<String> {
        let mut labels = Vec::new();
        for key in [
            "label",
            "gene",
            "locus_tag",
            "product",
            "standard_name",
            "note",
        ] {
            for value in feature.qualifier_values(key.into()) {
                let v = value.trim();
                if !v.is_empty() {
                    labels.push(v.to_string());
                }
            }
        }
        labels
    }

    fn resolve_anchor_position(
        dna: &DNAsequence,
        anchor: &AnchoredRegionAnchor,
    ) -> Result<usize, EngineError> {
        match anchor {
            AnchoredRegionAnchor::Position { zero_based } => {
                if *zero_based > dna.len() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Anchor position {} is out of bounds for sequence length {}",
                            zero_based,
                            dna.len()
                        ),
                    });
                }
                Ok(*zero_based)
            }
            AnchoredRegionAnchor::FeatureBoundary {
                feature_kind,
                feature_label,
                boundary,
                occurrence,
            } => {
                let kind_filter = feature_kind.as_ref().map(|s| s.to_ascii_uppercase());
                let label_filter = feature_label.as_ref().map(|s| s.to_ascii_uppercase());
                let mut matches = Vec::new();
                for feature in dna.features() {
                    if feature.kind.to_string().eq_ignore_ascii_case("SOURCE") {
                        continue;
                    }
                    if let Some(expected_kind) = &kind_filter {
                        if feature.kind.to_string().to_ascii_uppercase() != *expected_kind {
                            continue;
                        }
                    }
                    if let Some(expected_label) = &label_filter {
                        let labels = Self::feature_labels(feature);
                        let found_label = labels.iter().any(|label| {
                            let upper = label.to_ascii_uppercase();
                            upper == *expected_label || upper.contains(expected_label)
                        });
                        if !found_label {
                            continue;
                        }
                    }
                    let Ok((from, to)) = feature.location.find_bounds() else {
                        continue;
                    };
                    if from < 0 || to < 0 {
                        continue;
                    }
                    let pos = match boundary {
                        AnchorBoundary::Start => from as usize,
                        AnchorBoundary::End => to as usize,
                    };
                    if pos <= dna.len() {
                        matches.push(pos);
                    }
                }
                if matches.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::NotFound,
                        message: "No feature matched anchored-region anchor".to_string(),
                    });
                }
                matches.sort_unstable();
                let idx = occurrence.unwrap_or(0);
                matches.get(idx).cloned().ok_or_else(|| EngineError {
                    code: ErrorCode::NotFound,
                    message: format!(
                        "Feature anchor occurrence {} was requested, but only {} match(es) found",
                        idx,
                        matches.len()
                    ),
                })
            }
        }
    }

    fn anchored_range(
        anchor_pos: usize,
        len: usize,
        direction: &AnchorDirection,
        seq_len: usize,
        circular: bool,
    ) -> Option<(usize, usize)> {
        if len == 0 || seq_len == 0 {
            return None;
        }
        match direction {
            AnchorDirection::Upstream => {
                if !circular {
                    if anchor_pos > seq_len || anchor_pos < len {
                        return None;
                    }
                    Some((anchor_pos - len, anchor_pos))
                } else {
                    if anchor_pos > seq_len || len > seq_len {
                        return None;
                    }
                    let start = if anchor_pos >= len {
                        anchor_pos - len
                    } else {
                        seq_len - (len - anchor_pos)
                    };
                    Some((start, anchor_pos))
                }
            }
            AnchorDirection::Downstream => {
                if !circular {
                    if anchor_pos > seq_len || anchor_pos + len > seq_len {
                        return None;
                    }
                    Some((anchor_pos, anchor_pos + len))
                } else {
                    if anchor_pos > seq_len || len > seq_len {
                        return None;
                    }
                    Some((anchor_pos, anchor_pos + len))
                }
            }
        }
    }

    fn primer_options(seq: &str) -> Result<Vec<Vec<u8>>, EngineError> {
        let mut ret = Vec::with_capacity(seq.len());
        for b in seq.as_bytes() {
            let opts = IupacCode::from_letter(*b).to_vec();
            if opts.is_empty() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("Invalid primer base '{}'", *b as char),
                });
            }
            ret.push(opts);
        }
        Ok(ret)
    }

    fn total_primer_variants(options: &[Vec<u8>]) -> usize {
        options
            .iter()
            .fold(1usize, |acc, v| acc.saturating_mul(v.len().max(1)))
    }

    fn primer_variant_by_index(options: &[Vec<u8>], mut idx: usize) -> String {
        if options.is_empty() {
            return String::new();
        }
        let mut out: Vec<u8> = vec![b'A'; options.len()];
        for pos in (0..options.len()).rev() {
            let radix = options[pos].len();
            let choice = idx % radix;
            out[pos] = options[pos][choice];
            idx /= radix;
        }
        String::from_utf8(out).unwrap_or_default()
    }

    fn expand_primer_variants(
        spec: &PcrPrimerSpec,
        cap: usize,
    ) -> Result<Vec<String>, EngineError> {
        let normalized = Self::normalize_iupac_text(&spec.sequence)?;
        if normalized.is_empty() {
            return Ok(vec![]);
        }
        let options = Self::primer_options(&normalized)?;
        let total = Self::total_primer_variants(&options);
        if total == 0 {
            return Ok(vec![]);
        }

        let max_variants = spec.max_variants.unwrap_or(total).min(cap).max(1);
        let mode = spec
            .library_mode
            .clone()
            .unwrap_or(PrimerLibraryMode::Enumerate);
        match mode {
            PrimerLibraryMode::Enumerate => {
                if total > max_variants {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Primer variant space ({total}) exceeds max_variants ({max_variants}); use Sample mode or raise limits"
                        ),
                    });
                }
                let mut ret = Vec::with_capacity(total);
                for idx in 0..total {
                    ret.push(Self::primer_variant_by_index(&options, idx));
                }
                Ok(ret)
            }
            PrimerLibraryMode::Sample => {
                let target = max_variants.min(total);
                let mut chosen: Vec<usize> = Vec::with_capacity(target);
                let mut seen: HashSet<usize> = HashSet::with_capacity(target * 2);
                let mut state = spec.sample_seed.unwrap_or(0x9E3779B97F4A7C15);

                while chosen.len() < target {
                    state = state
                        .wrapping_mul(6364136223846793005)
                        .wrapping_add(1442695040888963407);
                    let idx = (state as usize) % total;
                    if seen.insert(idx) {
                        chosen.push(idx);
                    }
                    if seen.len() == total {
                        break;
                    }
                }
                chosen.sort_unstable();
                let mut ret = Vec::with_capacity(chosen.len());
                for idx in chosen {
                    ret.push(Self::primer_variant_by_index(&options, idx));
                }
                Ok(ret)
            }
        }
    }

    fn find_subsequence(haystack: &[u8], needle: &[u8], start: usize) -> Option<usize> {
        if needle.is_empty() || haystack.len() < needle.len() || start >= haystack.len() {
            return None;
        }
        let end = haystack.len() - needle.len();
        (start..=end).find(|idx| &haystack[*idx..*idx + needle.len()] == needle)
    }

    fn find_all_subsequences(haystack: &[u8], needle: &[u8]) -> Vec<usize> {
        let mut ret = vec![];
        if needle.is_empty() || haystack.len() < needle.len() {
            return ret;
        }
        let mut start = 0usize;
        while let Some(pos) = Self::find_subsequence(haystack, needle, start) {
            ret.push(pos);
            start = pos + 1;
            if start >= haystack.len() {
                break;
            }
        }
        ret
    }

    fn right_end_overhangs(dna: &DNAsequence) -> Vec<Vec<u8>> {
        let mut ret = Vec::new();
        if !dna.overhang().forward_3.is_empty() {
            ret.push(dna.overhang().forward_3.clone());
        }
        if !dna.overhang().reverse_5.is_empty() {
            ret.push(dna.overhang().reverse_5.clone());
        }
        ret
    }

    fn left_end_overhangs(dna: &DNAsequence) -> Vec<Vec<u8>> {
        let mut ret = Vec::new();
        if !dna.overhang().forward_5.is_empty() {
            ret.push(dna.overhang().forward_5.clone());
        }
        if !dna.overhang().reverse_3.is_empty() {
            ret.push(dna.overhang().reverse_3.clone());
        }
        ret
    }

    fn right_end_is_blunt(dna: &DNAsequence) -> bool {
        dna.overhang().forward_3.is_empty() && dna.overhang().reverse_5.is_empty()
    }

    fn left_end_is_blunt(dna: &DNAsequence) -> bool {
        dna.overhang().forward_5.is_empty() && dna.overhang().reverse_3.is_empty()
    }

    fn sticky_compatible(left: &DNAsequence, right: &DNAsequence) -> bool {
        let right_opts = Self::right_end_overhangs(left);
        let left_opts = Self::left_end_overhangs(right);
        if right_opts.is_empty() || left_opts.is_empty() {
            return false;
        }
        for r in &right_opts {
            let rc_r = Self::reverse_complement_bytes(r);
            let c_r = Self::complement_bytes(r);
            if left_opts.iter().any(|l| *l == rc_r || *l == c_r) {
                return true;
            }
        }
        false
    }

    fn find_anneal_sites(
        template: &[u8],
        anneal: &[u8],
        max_mismatches: usize,
        require_3prime_exact_bases: usize,
        three_prime_is_window_end: bool,
    ) -> Vec<usize> {
        let mut ret = vec![];
        if anneal.is_empty() || template.len() < anneal.len() {
            return ret;
        }
        let window_len = anneal.len();
        for start in 0..=(template.len() - window_len) {
            let window = &template[start..start + window_len];
            let mismatches = window
                .iter()
                .zip(anneal.iter())
                .filter(|(a, b)| a != b)
                .count();
            if mismatches > max_mismatches {
                continue;
            }
            if require_3prime_exact_bases > 0 {
                if require_3prime_exact_bases > window_len {
                    continue;
                }
                let exact_ok = if three_prime_is_window_end {
                    let from = window_len - require_3prime_exact_bases;
                    window[from..] == anneal[from..]
                } else {
                    window[..require_3prime_exact_bases] == anneal[..require_3prime_exact_bases]
                };
                if !exact_ok {
                    continue;
                }
            }
            ret.push(start);
        }
        ret
    }

    fn apply_internal(
        &mut self,
        op: Operation,
        run_id: &str,
        on_progress: &mut dyn FnMut(OperationProgress),
    ) -> Result<OpResult, EngineError> {
        self.reconcile_lineage_nodes();
        self.reconcile_containers();
        let op_for_containers = op.clone();
        let op = match op {
            Operation::MergeContainersById {
                container_ids,
                output_prefix,
            } => Operation::MergeContainers {
                inputs: self.flatten_container_members(&container_ids)?,
                output_prefix,
            },
            Operation::LigationContainer {
                container_id,
                circularize_if_possible,
                output_id,
                protocol,
                output_prefix,
                unique,
            } => Operation::Ligation {
                inputs: self.container_members(&container_id)?,
                circularize_if_possible,
                output_id,
                protocol,
                output_prefix,
                unique,
            },
            Operation::FilterContainerByMolecularWeight {
                container_id,
                min_bp,
                max_bp,
                error,
                unique,
                output_prefix,
            } => Operation::FilterByMolecularWeight {
                inputs: self.container_members(&container_id)?,
                min_bp,
                max_bp,
                error,
                unique,
                output_prefix,
            },
            other => other,
        };
        let op_id = self.next_op_id();
        let mut parent_seq_ids: Vec<SeqId> = vec![];
        let mut result = OpResult {
            op_id,
            created_seq_ids: vec![],
            changed_seq_ids: vec![],
            warnings: vec![],
            messages: vec![],
        };

        match op {
            Operation::LoadFile { path, as_id } => {
                let mut dna = GENtleApp::load_from_file(&path).map_err(|e| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("Could not load sequence file '{path}': {e}"),
                })?;
                Self::prepare_sequence(&mut dna);

                let base = as_id.unwrap_or_else(|| Self::derive_seq_id(&path));
                let seq_id = self.unique_seq_id(&base);
                self.state.sequences.insert(seq_id.clone(), dna);
                self.add_lineage_node(
                    &seq_id,
                    Self::classify_import_origin(
                        &path,
                        self.state
                            .sequences
                            .get(&seq_id)
                            .expect("sequence just inserted"),
                    ),
                    Some(&result.op_id),
                );
                let imported_anchor = self
                    .state
                    .sequences
                    .get(&seq_id)
                    .and_then(|loaded| Self::infer_imported_genbank_anchor(&path, loaded));
                if let Some(anchor) = imported_anchor {
                    self.append_genome_extraction_provenance(GenomeExtractionProvenance {
                        seq_id: seq_id.clone(),
                        recorded_at_unix_ms: Self::now_unix_ms(),
                        operation: "LoadFileGenBankRegion".to_string(),
                        genome_id: anchor.genome_id.clone(),
                        catalog_path: path.clone(),
                        cache_dir: None,
                        chromosome: Some(anchor.chromosome.clone()),
                        start_1based: Some(anchor.start_1based),
                        end_1based: Some(anchor.end_1based),
                        gene_query: None,
                        occurrence: None,
                        gene_id: None,
                        gene_name: None,
                        strand: None,
                        anchor_strand: anchor.strand,
                        sequence_source_type: Some("genbank_file".to_string()),
                        annotation_source_type: Some("genbank_file".to_string()),
                        sequence_source: Some(path.clone()),
                        annotation_source: Some(path.clone()),
                        sequence_sha1: None,
                        annotation_sha1: None,
                    });
                    let strand = anchor.strand.unwrap_or('+');
                    result.messages.push(format!(
                        "Detected GenBank genome anchor for '{}': {}:{}-{} ({}, strand {})",
                        seq_id,
                        anchor.chromosome,
                        anchor.start_1based,
                        anchor.end_1based,
                        anchor.genome_id,
                        strand
                    ));
                }
                result.created_seq_ids.push(seq_id.clone());
                result
                    .messages
                    .push(format!("Loaded '{path}' as '{seq_id}'"));
            }
            Operation::SaveFile {
                seq_id,
                path,
                format,
            } => {
                let dna = self
                    .state
                    .sequences
                    .get(&seq_id)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{seq_id}' not found"),
                    })?;

                match format {
                    ExportFormat::GenBank => {
                        dna.write_genbank_file(&path).map_err(|e| EngineError {
                            code: ErrorCode::Io,
                            message: format!("Could not write GenBank file '{path}': {e}"),
                        })?;
                    }
                    ExportFormat::Fasta => Self::save_as_fasta(&seq_id, dna, &path)?,
                }

                result.changed_seq_ids.push(seq_id.clone());
                result
                    .messages
                    .push(format!("Wrote '{seq_id}' to '{path}'"));
            }
            Operation::RenderSequenceSvg { seq_id, mode, path } => {
                let dna = self
                    .state
                    .sequences
                    .get(&seq_id)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{seq_id}' not found"),
                    })?;
                let svg = match mode {
                    RenderSvgMode::Linear => export_linear_svg(dna, &self.state.display),
                    RenderSvgMode::Circular => export_circular_svg(dna, &self.state.display),
                };
                std::fs::write(&path, svg).map_err(|e| EngineError {
                    code: ErrorCode::Io,
                    message: format!("Could not write SVG output '{path}': {e}"),
                })?;
                result.messages.push(format!(
                    "Wrote {:?} SVG for '{}' to '{}'",
                    mode, seq_id, path
                ));
            }
            Operation::RenderRnaStructureSvg { seq_id, path } => {
                let report = self.render_rna_structure_svg_to_path(&seq_id, &path)?;
                result.messages.push(format!(
                    "Wrote RNA structure SVG for '{}' to '{}' using {}",
                    seq_id, path, report.tool
                ));
            }
            Operation::RenderLineageSvg { path } => {
                let svg = export_lineage_svg(&self.state);
                std::fs::write(&path, svg).map_err(|e| EngineError {
                    code: ErrorCode::Io,
                    message: format!("Could not write SVG output '{path}': {e}"),
                })?;
                result
                    .messages
                    .push(format!("Wrote lineage SVG to '{}'", path));
            }
            Operation::RenderPoolGelSvg {
                inputs,
                path,
                ladders,
            } => {
                if inputs.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "RenderPoolGelSvg requires at least one input sequence id"
                            .to_string(),
                    });
                }
                let mut members: Vec<(String, usize)> = Vec::with_capacity(inputs.len());
                for seq_id in &inputs {
                    let dna = self
                        .state
                        .sequences
                        .get(seq_id)
                        .ok_or_else(|| EngineError {
                            code: ErrorCode::NotFound,
                            message: format!("Sequence '{seq_id}' not found"),
                        })?;
                    members.push((seq_id.clone(), dna.len()));
                }
                let ladder_names = ladders
                    .unwrap_or_default()
                    .into_iter()
                    .map(|v| v.trim().to_string())
                    .filter(|v| !v.is_empty())
                    .collect::<Vec<_>>();
                let layout =
                    build_pool_gel_layout(&members, &ladder_names).map_err(|e| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: e,
                    })?;
                let svg = export_pool_gel_svg(&layout);
                std::fs::write(&path, svg).map_err(|e| EngineError {
                    code: ErrorCode::Io,
                    message: format!("Could not write SVG output '{path}': {e}"),
                })?;
                let ladders_used = if layout.selected_ladders.is_empty() {
                    "auto".to_string()
                } else {
                    layout.selected_ladders.join(" + ")
                };
                result.messages.push(format!(
                    "Wrote pool gel SVG for {} sequence(s) to '{}' (ladders: {})",
                    members.len(),
                    path,
                    ladders_used
                ));
            }
            Operation::ExportDnaLadders { path, name_filter } => {
                let report = Self::export_dna_ladders(&path, name_filter.as_deref())?;
                let filter_text = name_filter
                    .as_deref()
                    .map(str::trim)
                    .filter(|v| !v.is_empty())
                    .unwrap_or("-");
                result.messages.push(format!(
                    "Wrote DNA ladders catalog ({} ladder(s), filter={}) to '{}'",
                    report.ladder_count, filter_text, report.path
                ));
            }
            Operation::ExportRnaLadders { path, name_filter } => {
                let report = Self::export_rna_ladders(&path, name_filter.as_deref())?;
                let filter_text = name_filter
                    .as_deref()
                    .map(str::trim)
                    .filter(|v| !v.is_empty())
                    .unwrap_or("-");
                result.messages.push(format!(
                    "Wrote RNA ladders catalog ({} ladder(s), filter={}) to '{}'",
                    report.ladder_count, filter_text, report.path
                ));
            }
            Operation::ExportPool {
                inputs,
                path,
                pool_id,
                human_id,
            } => {
                let (pool_id, human_id, count) =
                    self.export_pool_file(&inputs, &path, pool_id, human_id)?;
                result.messages.push(format!(
                    "Wrote pool export '{}' ({}) with {} members to '{}'",
                    pool_id, human_id, count, path
                ));
            }
            Operation::PrepareGenome {
                genome_id,
                catalog_path,
                cache_dir,
            } => {
                let report = Self::prepare_reference_genome_once(
                    &genome_id,
                    catalog_path.as_deref(),
                    cache_dir.as_deref(),
                    &mut |p| on_progress(OperationProgress::GenomePrepare(p)),
                )?;
                result.messages.push(Self::format_prepare_genome_message(
                    &genome_id,
                    cache_dir.as_deref(),
                    &report,
                ));
                if !report.warnings.is_empty() {
                    result.warnings.extend(report.warnings.clone());
                }
            }
            Operation::ExtractGenomeRegion {
                genome_id,
                chromosome,
                start_1based,
                end_1based,
                output_id,
                catalog_path,
                cache_dir,
            } => {
                let catalog_path =
                    catalog_path.unwrap_or_else(|| DEFAULT_GENOME_CATALOG_PATH.to_string());
                let catalog =
                    GenomeCatalog::from_json_file(&catalog_path).map_err(|e| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("Could not open genome catalog '{catalog_path}': {e}"),
                    })?;
                let sequence = catalog
                    .get_sequence_region_with_cache(
                        &genome_id,
                        &chromosome,
                        start_1based,
                        end_1based,
                        cache_dir.as_deref(),
                    )
                    .map_err(|e| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!(
                            "Could not load genome region {}:{}-{} from '{}': {}",
                            chromosome, start_1based, end_1based, genome_id, e
                        ),
                    })?;
                let default_id = format!(
                    "{}_{}_{}_{}",
                    Self::normalize_id_token(&genome_id),
                    Self::normalize_id_token(&chromosome),
                    start_1based,
                    end_1based
                );
                let base = output_id.unwrap_or(default_id);
                let seq_id = self.import_genome_slice_sequence(&mut result, sequence, base)?;
                let source_plan = catalog.source_plan(&genome_id, cache_dir.as_deref()).ok();
                let inspection = catalog
                    .inspect_prepared_genome(&genome_id, cache_dir.as_deref())
                    .ok()
                    .flatten();
                let (
                    sequence_source_type,
                    annotation_source_type,
                    sequence_source,
                    annotation_source,
                    sequence_sha1,
                    annotation_sha1,
                ) = Self::genome_source_snapshot(source_plan.as_ref(), inspection.as_ref());
                self.append_genome_extraction_provenance(GenomeExtractionProvenance {
                    seq_id: seq_id.clone(),
                    recorded_at_unix_ms: Self::now_unix_ms(),
                    operation: "ExtractGenomeRegion".to_string(),
                    genome_id: genome_id.clone(),
                    catalog_path: catalog_path.clone(),
                    cache_dir: cache_dir.clone(),
                    chromosome: Some(chromosome.clone()),
                    start_1based: Some(start_1based),
                    end_1based: Some(end_1based),
                    gene_query: None,
                    occurrence: None,
                    gene_id: None,
                    gene_name: None,
                    strand: None,
                    anchor_strand: Some('+'),
                    sequence_source_type,
                    annotation_source_type,
                    sequence_source,
                    annotation_source,
                    sequence_sha1,
                    annotation_sha1,
                });
                result.messages.push(format!(
                    "Extracted genome region {}:{}-{} from '{}' as '{}'",
                    chromosome, start_1based, end_1based, genome_id, seq_id
                ));
            }
            Operation::ExtractGenomeGene {
                genome_id,
                gene_query,
                occurrence,
                output_id,
                catalog_path,
                cache_dir,
            } => {
                let query = gene_query.trim();
                if query.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "Gene query cannot be empty".to_string(),
                    });
                }
                let catalog_path =
                    catalog_path.unwrap_or_else(|| DEFAULT_GENOME_CATALOG_PATH.to_string());
                let catalog =
                    GenomeCatalog::from_json_file(&catalog_path).map_err(|e| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("Could not open genome catalog '{catalog_path}': {e}"),
                    })?;
                let genes = catalog
                    .list_gene_regions(&genome_id, cache_dir.as_deref())
                    .map_err(|e| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!(
                            "Could not load gene index for genome '{}': {}",
                            genome_id, e
                        ),
                    })?;
                let mut exact_matches: Vec<&GenomeGeneRecord> = genes
                    .iter()
                    .filter(|record| Self::genome_gene_matches_exact(record, query))
                    .collect();
                let used_fuzzy = if exact_matches.is_empty() {
                    let query_lower = query.to_ascii_lowercase();
                    exact_matches = genes
                        .iter()
                        .filter(|record| Self::genome_gene_matches_contains(record, &query_lower))
                        .collect();
                    true
                } else {
                    false
                };
                if exact_matches.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("No genes in '{}' match query '{}'", genome_id, query),
                    });
                }
                let requested_occurrence = occurrence;
                let occurrence = requested_occurrence.unwrap_or(1);
                if occurrence == 0 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "Gene occurrence must be >= 1".to_string(),
                    });
                }
                if exact_matches.len() > 1 && requested_occurrence.is_none() {
                    result.warnings.push(format!(
                        "Gene query '{}' matched {} records in '{}'; using first match (set occurrence for another match).",
                        query,
                        exact_matches.len(),
                        genome_id
                    ));
                }
                let Some(selected_gene) = exact_matches.get(occurrence - 1) else {
                    return Err(EngineError {
                        code: ErrorCode::NotFound,
                        message: format!(
                            "Gene query '{}' matched {} records, but occurrence {} was requested",
                            query,
                            exact_matches.len(),
                            occurrence
                        ),
                    });
                };
                let sequence = catalog
                    .get_sequence_region_with_cache(
                        &genome_id,
                        &selected_gene.chromosome,
                        selected_gene.start_1based,
                        selected_gene.end_1based,
                        cache_dir.as_deref(),
                    )
                    .map_err(|e| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!(
                            "Could not load gene region {}:{}-{} from '{}': {}",
                            selected_gene.chromosome,
                            selected_gene.start_1based,
                            selected_gene.end_1based,
                            genome_id,
                            e
                        ),
                    })?;
                let gene_label = selected_gene
                    .gene_name
                    .as_ref()
                    .or(selected_gene.gene_id.as_ref())
                    .cloned()
                    .unwrap_or_else(|| query.to_string());
                let default_id = format!(
                    "{}_{}_{}_{}",
                    Self::normalize_id_token(&genome_id),
                    Self::normalize_id_token(&gene_label),
                    selected_gene.start_1based,
                    selected_gene.end_1based
                );
                let base = output_id.unwrap_or(default_id);
                let seq_id = self.import_genome_slice_sequence(&mut result, sequence, base)?;
                let source_plan = catalog.source_plan(&genome_id, cache_dir.as_deref()).ok();
                let inspection = catalog
                    .inspect_prepared_genome(&genome_id, cache_dir.as_deref())
                    .ok()
                    .flatten();
                let (
                    sequence_source_type,
                    annotation_source_type,
                    sequence_source,
                    annotation_source,
                    sequence_sha1,
                    annotation_sha1,
                ) = Self::genome_source_snapshot(source_plan.as_ref(), inspection.as_ref());
                self.append_genome_extraction_provenance(GenomeExtractionProvenance {
                    seq_id: seq_id.clone(),
                    recorded_at_unix_ms: Self::now_unix_ms(),
                    operation: "ExtractGenomeGene".to_string(),
                    genome_id: genome_id.clone(),
                    catalog_path: catalog_path.clone(),
                    cache_dir: cache_dir.clone(),
                    chromosome: Some(selected_gene.chromosome.clone()),
                    start_1based: Some(selected_gene.start_1based),
                    end_1based: Some(selected_gene.end_1based),
                    gene_query: Some(query.to_string()),
                    occurrence: Some(occurrence),
                    gene_id: selected_gene.gene_id.clone(),
                    gene_name: selected_gene.gene_name.clone(),
                    strand: selected_gene.strand,
                    anchor_strand: Some('+'),
                    sequence_source_type,
                    annotation_source_type,
                    sequence_source,
                    annotation_source,
                    sequence_sha1,
                    annotation_sha1,
                });
                let match_mode = if used_fuzzy { "fuzzy" } else { "exact" };
                result.messages.push(format!(
                    "Extracted genome gene '{}' [{} match, occurrence {}] as '{}' from '{}' ({})",
                    query,
                    match_mode,
                    occurrence,
                    seq_id,
                    genome_id,
                    Self::genome_gene_display_label(selected_gene)
                ));
            }
            Operation::ImportGenomeBedTrack {
                seq_id,
                path,
                track_name,
                min_score,
                max_score,
                clear_existing,
            } => {
                if path.trim().is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "ImportGenomeBedTrack requires a non-empty BED path".to_string(),
                    });
                }
                if min_score
                    .zip(max_score)
                    .map(|(min, max)| min > max)
                    .unwrap_or(false)
                {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "ImportGenomeBedTrack requires min_score <= max_score".to_string(),
                    });
                }

                let anchor = self.latest_genome_anchor_for_seq(&seq_id)?;
                let _ = self.ensure_lineage_node(&seq_id);
                let dna = self
                    .state
                    .sequences
                    .get_mut(&seq_id)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{seq_id}' not found"),
                    })?;

                let report = Self::import_genome_bed_track(
                    dna,
                    &anchor,
                    &path,
                    track_name.as_deref(),
                    min_score,
                    max_score,
                    clear_existing.unwrap_or(false),
                )?;

                result.changed_seq_ids.push(seq_id.clone());
                result.warnings.extend(report.warnings);
                let anchor_strand = anchor.strand.unwrap_or('+');
                result.messages.push(format!(
                    "Imported {} BED feature(s) into '{}' from '{}' as track '{}' (anchor={} {}:{}-{} strand {}, parsed={}, skipped={})",
                    report.imported_features,
                    seq_id,
                    path,
                    report.track_name,
                    anchor.genome_id,
                    anchor.chromosome,
                    anchor.start_1based,
                    anchor.end_1based,
                    anchor_strand,
                    report.parsed_records,
                    report.skipped_records
                ));
                if report.skipped_missing_score > 0 {
                    result.warnings.push(format!(
                        "{} BED record(s) were skipped because score filters were set but the BED score column was missing",
                        report.skipped_missing_score
                    ));
                }
                if report.skipped_outside_score_range > 0 {
                    result.messages.push(format!(
                        "{} BED record(s) were outside score filter bounds",
                        report.skipped_outside_score_range
                    ));
                }
                if report.truncated_at_limit {
                    result.warnings.push(format!(
                        "BED import was truncated after {} features (limit={})",
                        report.imported_features, MAX_IMPORTED_SIGNAL_FEATURES
                    ));
                }
            }
            Operation::ImportGenomeBigWigTrack {
                seq_id,
                path,
                track_name,
                min_score,
                max_score,
                clear_existing,
            } => {
                if path.trim().is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "ImportGenomeBigWigTrack requires a non-empty BigWig path"
                            .to_string(),
                    });
                }
                if min_score
                    .zip(max_score)
                    .map(|(min, max)| min > max)
                    .unwrap_or(false)
                {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "ImportGenomeBigWigTrack requires min_score <= max_score"
                            .to_string(),
                    });
                }

                let anchor = self.latest_genome_anchor_for_seq(&seq_id)?;
                let _ = self.ensure_lineage_node(&seq_id);
                let dna = self
                    .state
                    .sequences
                    .get_mut(&seq_id)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{seq_id}' not found"),
                    })?;

                let report = Self::import_genome_bigwig_track(
                    dna,
                    &anchor,
                    &path,
                    track_name.as_deref(),
                    min_score,
                    max_score,
                    clear_existing.unwrap_or(false),
                )?;

                result.changed_seq_ids.push(seq_id.clone());
                result.warnings.extend(report.warnings);
                let anchor_strand = anchor.strand.unwrap_or('+');
                result.messages.push(format!(
                    "Imported {} BigWig feature(s) into '{}' from '{}' as track '{}' (anchor={} {}:{}-{} strand {}, parsed={}, skipped={})",
                    report.imported_features,
                    seq_id,
                    path,
                    report.track_name,
                    anchor.genome_id,
                    anchor.chromosome,
                    anchor.start_1based,
                    anchor.end_1based,
                    anchor_strand,
                    report.parsed_records,
                    report.skipped_records
                ));
                if report.skipped_missing_score > 0 {
                    result.warnings.push(format!(
                        "{} converted bedGraph record(s) were skipped because score filters were set but no value was available",
                        report.skipped_missing_score
                    ));
                }
                if report.skipped_outside_score_range > 0 {
                    result.messages.push(format!(
                        "{} converted bedGraph record(s) were outside score filter bounds",
                        report.skipped_outside_score_range
                    ));
                }
                if report.truncated_at_limit {
                    result.warnings.push(format!(
                        "BigWig import was truncated after {} features (limit={})",
                        report.imported_features, MAX_IMPORTED_SIGNAL_FEATURES
                    ));
                }
            }
            Operation::ImportGenomeVcfTrack {
                seq_id,
                path,
                track_name,
                min_score,
                max_score,
                clear_existing,
            } => {
                if path.trim().is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "ImportGenomeVcfTrack requires a non-empty VCF path".to_string(),
                    });
                }
                if min_score
                    .zip(max_score)
                    .map(|(min, max)| min > max)
                    .unwrap_or(false)
                {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "ImportGenomeVcfTrack requires min_score <= max_score".to_string(),
                    });
                }

                let anchor = self.latest_genome_anchor_for_seq(&seq_id)?;
                let _ = self.ensure_lineage_node(&seq_id);
                let dna = self
                    .state
                    .sequences
                    .get_mut(&seq_id)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{seq_id}' not found"),
                    })?;

                let report = Self::import_genome_vcf_track(
                    dna,
                    &anchor,
                    &path,
                    track_name.as_deref(),
                    min_score,
                    max_score,
                    clear_existing.unwrap_or(false),
                )?;

                result.changed_seq_ids.push(seq_id.clone());
                result.warnings.extend(report.warnings);
                let anchor_strand = anchor.strand.unwrap_or('+');
                result.messages.push(format!(
                    "Imported {} VCF feature(s) into '{}' from '{}' as track '{}' (anchor={} {}:{}-{} strand {}, parsed={}, skipped={})",
                    report.imported_features,
                    seq_id,
                    path,
                    report.track_name,
                    anchor.genome_id,
                    anchor.chromosome,
                    anchor.start_1based,
                    anchor.end_1based,
                    anchor_strand,
                    report.parsed_records,
                    report.skipped_records
                ));
                if report.skipped_missing_score > 0 {
                    result.warnings.push(format!(
                        "{} VCF record(s) were skipped because QUAL-based score filters were set but QUAL was missing",
                        report.skipped_missing_score
                    ));
                }
                if report.skipped_outside_score_range > 0 {
                    result.messages.push(format!(
                        "{} VCF record(s) were outside QUAL score filter bounds",
                        report.skipped_outside_score_range
                    ));
                }
                if report.truncated_at_limit {
                    result.warnings.push(format!(
                        "VCF import was truncated after {} features (limit={})",
                        report.imported_features, MAX_IMPORTED_SIGNAL_FEATURES
                    ));
                }
            }
            Operation::ImportBlastHitsTrack {
                seq_id,
                hits,
                track_name,
                clear_existing,
            } => {
                if hits.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "ImportBlastHitsTrack requires at least one hit".to_string(),
                    });
                }
                let selected_track_name = track_name
                    .as_deref()
                    .map(str::trim)
                    .filter(|v| !v.is_empty())
                    .unwrap_or("blast_hits")
                    .to_string();
                let clear_existing = clear_existing.unwrap_or(false);

                let _ = self.ensure_lineage_node(&seq_id);
                let dna = self
                    .state
                    .sequences
                    .get_mut(&seq_id)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{seq_id}' not found"),
                    })?;
                let seq_len = dna.len();
                if seq_len == 0 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Sequence '{}' is empty and cannot receive BLAST hit features",
                            seq_id
                        ),
                    });
                }

                let mut removed_count = 0usize;
                if clear_existing {
                    let before = dna.features().len();
                    Self::remove_generated_blast_hit_features(dna.features_mut());
                    removed_count = before.saturating_sub(dna.features().len());
                }

                let mut imported_count = 0usize;
                let mut skipped_count = 0usize;
                for (idx, hit) in hits.iter().enumerate() {
                    let start = hit.query_start_1based.min(hit.query_end_1based);
                    let end = hit.query_start_1based.max(hit.query_end_1based);
                    if start == 0 || end == 0 {
                        skipped_count += 1;
                        if result.warnings.len() < 20 {
                            result.warnings.push(format!(
                                "BLAST hit {} skipped because query coordinates are invalid: {}..{}",
                                idx + 1,
                                hit.query_start_1based,
                                hit.query_end_1based
                            ));
                        }
                        continue;
                    }
                    if start > seq_len {
                        skipped_count += 1;
                        if result.warnings.len() < 20 {
                            result.warnings.push(format!(
                                "BLAST hit {} skipped because query start {} is outside sequence length {}",
                                idx + 1,
                                start,
                                seq_len
                            ));
                        }
                        continue;
                    }
                    let end_clamped = end.min(seq_len);
                    if end_clamped < start {
                        skipped_count += 1;
                        continue;
                    }
                    if end_clamped < end && result.warnings.len() < 20 {
                        result.warnings.push(format!(
                            "BLAST hit {} query range {}..{} was clamped to {}..{} for sequence length {}",
                            idx + 1,
                            start,
                            end,
                            start,
                            end_clamped,
                            seq_len
                        ));
                    }
                    let local_start_0based = start.saturating_sub(1);
                    let local_end_0based_exclusive = end_clamped;
                    let local_strand =
                        if hit.subject_start_1based == 0 || hit.subject_end_1based == 0 {
                            None
                        } else if hit.subject_start_1based <= hit.subject_end_1based {
                            Some('+')
                        } else {
                            Some('-')
                        };
                    let feature = Self::build_blast_hit_feature(
                        hit,
                        &selected_track_name,
                        local_start_0based,
                        local_end_0based_exclusive,
                        local_strand,
                    );
                    dna.features_mut().push(feature);
                    imported_count += 1;
                }

                if imported_count > 0 || removed_count > 0 {
                    result.changed_seq_ids.push(seq_id.clone());
                }
                result.messages.push(format!(
                    "Imported {} BLAST hit feature(s) into '{}' as track '{}' (input_hits={}, skipped={}, cleared_existing={})",
                    imported_count,
                    seq_id,
                    selected_track_name,
                    hits.len(),
                    skipped_count,
                    removed_count
                ));
            }
            Operation::Digest {
                input,
                enzymes,
                output_prefix,
            } => {
                parent_seq_ids.push(input.clone());
                let dna = self
                    .state
                    .sequences
                    .get(&input)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{input}' not found"),
                    })?
                    .clone();

                let (found, missing) = self.resolve_enzymes(&enzymes)?;
                if !missing.is_empty() {
                    result
                        .warnings
                        .push(format!("Unknown enzymes ignored: {}", missing.join(",")));
                }

                let fragments =
                    Self::digest_with_guard(&dna, found, self.max_fragments_per_container())?;
                let prefix = output_prefix.unwrap_or_else(|| format!("{input}_digest"));

                for (i, mut fragment) in fragments.into_iter().enumerate() {
                    // Keep digest interactive by deferring expensive feature recomputation.
                    Self::prepare_sequence_light(&mut fragment);
                    let candidate = format!("{}_{}", prefix, i + 1);
                    let seq_id = self.unique_seq_id(&candidate);
                    self.state.sequences.insert(seq_id.clone(), fragment);
                    self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
                    result.created_seq_ids.push(seq_id);
                }

                result.messages.push(format!(
                    "Digest created {} fragment(s); feature recomputation deferred",
                    result.created_seq_ids.len()
                ));
            }
            Operation::DigestContainer {
                container_id,
                enzymes,
                output_prefix,
            } => {
                let inputs = self.container_members(&container_id)?;
                parent_seq_ids.extend(inputs.clone());
                let (found, missing) = self.resolve_enzymes(&enzymes)?;
                if !missing.is_empty() {
                    result
                        .warnings
                        .push(format!("Unknown enzymes ignored: {}", missing.join(",")));
                }
                let prefix = output_prefix.unwrap_or_else(|| format!("{container_id}_digest"));
                for input in inputs {
                    let dna = self
                        .state
                        .sequences
                        .get(&input)
                        .ok_or_else(|| EngineError {
                            code: ErrorCode::NotFound,
                            message: format!("Sequence '{input}' not found"),
                        })?
                        .clone();
                    let fragments = Self::digest_with_guard(
                        &dna,
                        found.clone(),
                        self.max_fragments_per_container(),
                    )?;
                    for (i, mut fragment) in fragments.into_iter().enumerate() {
                        // Keep digest interactive by deferring expensive feature recomputation.
                        Self::prepare_sequence_light(&mut fragment);
                        let candidate = format!("{}_{}_{}", prefix, input, i + 1);
                        let seq_id = self.unique_seq_id(&candidate);
                        self.state.sequences.insert(seq_id.clone(), fragment);
                        self.add_lineage_node(
                            &seq_id,
                            SequenceOrigin::Derived,
                            Some(&result.op_id),
                        );
                        result.created_seq_ids.push(seq_id);
                        if result.created_seq_ids.len() > self.max_fragments_per_container() {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: format!(
                                    "Digest produced more than max_fragments_per_container={}",
                                    self.max_fragments_per_container()
                                ),
                            });
                        }
                    }
                }

                result.messages.push(format!(
                    "Digest container '{}' created {} fragment(s); feature recomputation deferred",
                    container_id,
                    result.created_seq_ids.len()
                ));
            }
            Operation::MergeContainersById { .. }
            | Operation::LigationContainer { .. }
            | Operation::FilterContainerByMolecularWeight { .. } => {
                unreachable!("container operation variants are normalized before execution")
            }
            Operation::MergeContainers {
                inputs,
                output_prefix,
            } => {
                if inputs.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "MergeContainers requires at least one input sequence".to_string(),
                    });
                }
                if inputs.len() > self.max_fragments_per_container() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "MergeContainers input count {} exceeds max_fragments_per_container={}",
                            inputs.len(),
                            self.max_fragments_per_container()
                        ),
                    });
                }
                parent_seq_ids.extend(inputs.clone());
                let prefix = output_prefix.unwrap_or_else(|| "merged".to_string());
                for (i, input) in inputs.iter().enumerate() {
                    let dna = self
                        .state
                        .sequences
                        .get(input)
                        .ok_or_else(|| EngineError {
                            code: ErrorCode::NotFound,
                            message: format!("Sequence '{input}' not found"),
                        })?
                        .clone();
                    let seq_id = self.unique_seq_id(&format!("{}_{}", prefix, i + 1));
                    self.state.sequences.insert(seq_id.clone(), dna);
                    self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
                    result.created_seq_ids.push(seq_id);
                }
                result.messages.push(format!(
                    "Merged {} input sequence(s) into container prefix '{}'",
                    inputs.len(),
                    prefix
                ));
            }
            Operation::Ligation {
                inputs,
                circularize_if_possible,
                output_id,
                protocol,
                output_prefix,
                unique,
            } => {
                parent_seq_ids.extend(inputs.clone());
                if inputs.len() < 2 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "Ligation requires at least two input sequences".to_string(),
                    });
                }
                if 1 > self.max_fragments_per_container() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Ligation product count exceeds max_fragments_per_container={}",
                            self.max_fragments_per_container()
                        ),
                    });
                }
                let mut accepted: Vec<(String, String, String)> = vec![];
                for (i, left_id) in inputs.iter().enumerate() {
                    for (j, right_id) in inputs.iter().enumerate() {
                        if i == j {
                            continue;
                        }
                        let left =
                            self.state
                                .sequences
                                .get(left_id)
                                .ok_or_else(|| EngineError {
                                    code: ErrorCode::NotFound,
                                    message: format!("Sequence '{left_id}' not found"),
                                })?;
                        let right =
                            self.state
                                .sequences
                                .get(right_id)
                                .ok_or_else(|| EngineError {
                                    code: ErrorCode::NotFound,
                                    message: format!("Sequence '{right_id}' not found"),
                                })?;

                        let ok = match protocol {
                            LigationProtocol::Sticky => Self::sticky_compatible(left, right),
                            LigationProtocol::Blunt => {
                                Self::right_end_is_blunt(left) && Self::left_end_is_blunt(right)
                            }
                        };
                        if !ok {
                            continue;
                        }

                        let product = format!(
                            "{}{}",
                            left.get_forward_string(),
                            right.get_forward_string()
                        );
                        accepted.push((left_id.clone(), right_id.clone(), product));
                        if accepted.len() > self.max_fragments_per_container() {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: format!(
                                    "Ligation produced more than max_fragments_per_container={}",
                                    self.max_fragments_per_container()
                                ),
                            });
                        }
                    }
                }

                if accepted.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "No ligation products found for protocol '{:?}'",
                            protocol
                        ),
                    });
                }

                let require_unique = unique.unwrap_or(false);
                if require_unique && accepted.len() != 1 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Ligation unique=true requires exactly one product, found {}",
                            accepted.len()
                        ),
                    });
                }
                if output_id.is_some() && accepted.len() != 1 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "Ligation output_id can only be used when exactly one product is produced"
                            .to_string(),
                    });
                }

                let prefix = output_prefix.unwrap_or_else(|| "ligation".to_string());
                for (idx, (left_id, right_id, merged)) in accepted.into_iter().enumerate() {
                    let mut product =
                        DNAsequence::from_sequence(&merged).map_err(|e| EngineError {
                            code: ErrorCode::Internal,
                            message: format!("Could not create ligation product: {e}"),
                        })?;
                    product.set_circular(circularize_if_possible);
                    Self::prepare_sequence(&mut product);

                    let seq_id = if idx == 0 {
                        if let Some(id) = output_id.clone() {
                            self.unique_seq_id(&id)
                        } else {
                            self.unique_seq_id(&format!("{}_{}", prefix, idx + 1))
                        }
                    } else {
                        self.unique_seq_id(&format!("{}_{}", prefix, idx + 1))
                    };
                    self.state.sequences.insert(seq_id.clone(), product);
                    self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
                    result.created_seq_ids.push(seq_id.clone());
                    result.messages.push(format!(
                        "Ligation product '{}' from '{}' + '{}'",
                        seq_id, left_id, right_id
                    ));
                }
            }
            Operation::Pcr {
                template,
                forward_primer,
                reverse_primer,
                output_id,
                unique,
            } => {
                parent_seq_ids.push(template.clone());
                let dna = self
                    .state
                    .sequences
                    .get(&template)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{template}' not found"),
                    })?
                    .clone();

                if dna.is_circular() {
                    return Err(EngineError {
                        code: ErrorCode::Unsupported,
                        message: "PCR on circular templates is not implemented yet".to_string(),
                    });
                }

                let template_seq = dna.get_forward_string().to_ascii_uppercase();
                let template_bytes = template_seq.as_bytes();
                let fwd = Self::normalize_dna_text(&forward_primer);
                let rev = Self::normalize_dna_text(&reverse_primer);
                if fwd.is_empty() || rev.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "PCR primers must not be empty".to_string(),
                    });
                }

                let rev_binding = Self::reverse_complement(&rev);
                let fwd_sites = Self::find_all_subsequences(template_bytes, fwd.as_bytes());
                if fwd_sites.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "Forward primer not found on template".to_string(),
                    });
                }
                let rev_sites = Self::find_all_subsequences(template_bytes, rev_binding.as_bytes());
                if rev_sites.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "Reverse primer binding site not found on template".to_string(),
                    });
                }

                let mut amplicon_ranges: Vec<(usize, usize)> = vec![];
                for fwd_pos in &fwd_sites {
                    for rev_pos in &rev_sites {
                        if *rev_pos < *fwd_pos {
                            continue;
                        }
                        let amplicon_end = rev_pos + rev_binding.len();
                        if amplicon_end <= template_seq.len() {
                            amplicon_ranges.push((*fwd_pos, amplicon_end));
                        }
                    }
                }
                amplicon_ranges.sort_unstable();
                amplicon_ranges.dedup();

                if amplicon_ranges.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "No valid forward/reverse primer pair produced an amplicon"
                            .to_string(),
                    });
                }
                if amplicon_ranges.len() > self.max_fragments_per_container() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "PCR produced {} amplicons, exceeding max_fragments_per_container={}",
                            amplicon_ranges.len(),
                            self.max_fragments_per_container()
                        ),
                    });
                }

                let require_unique = unique.unwrap_or(false);
                if require_unique && amplicon_ranges.len() != 1 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "PCR unique=true requires exactly one amplicon, found {}",
                            amplicon_ranges.len()
                        ),
                    });
                }
                if output_id.is_some() && amplicon_ranges.len() != 1 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message:
                            "PCR output_id can only be used when exactly one amplicon is produced"
                                .to_string(),
                    });
                }

                let default_base = format!("{template}_pcr");
                for (i, (start, end)) in amplicon_ranges.iter().enumerate() {
                    let amplicon = &template_seq[*start..*end];
                    let mut pcr_product =
                        DNAsequence::from_sequence(amplicon).map_err(|e| EngineError {
                            code: ErrorCode::Internal,
                            message: format!("Could not create PCR product: {e}"),
                        })?;
                    pcr_product.set_circular(false);
                    Self::prepare_sequence(&mut pcr_product);

                    let requested = if i == 0 { output_id.clone() } else { None };
                    let seq_id = if let Some(id) = requested {
                        self.unique_seq_id(&id)
                    } else if amplicon_ranges.len() == 1 {
                        self.unique_seq_id(&default_base)
                    } else {
                        self.unique_seq_id(&format!("{}_{}", default_base, i + 1))
                    };
                    self.state.sequences.insert(seq_id.clone(), pcr_product);
                    self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
                    result.created_seq_ids.push(seq_id.clone());
                    result.messages.push(format!(
                        "PCR product '{}' created: {}..{} (len {})",
                        seq_id,
                        start,
                        end,
                        end - start
                    ));
                }
            }
            Operation::PcrAdvanced {
                template,
                forward_primer,
                reverse_primer,
                output_id,
                unique,
            } => {
                parent_seq_ids.push(template.clone());
                let dna = self
                    .state
                    .sequences
                    .get(&template)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{template}' not found"),
                    })?
                    .clone();

                if dna.is_circular() {
                    return Err(EngineError {
                        code: ErrorCode::Unsupported,
                        message: "PCR on circular templates is not implemented yet".to_string(),
                    });
                }

                let template_seq = dna.get_forward_string().to_ascii_uppercase();
                let template_bytes = template_seq.as_bytes();

                let fwd_variants = Self::expand_primer_variants(
                    &forward_primer,
                    self.max_fragments_per_container(),
                )?;
                let rev_variants = Self::expand_primer_variants(
                    &reverse_primer,
                    self.max_fragments_per_container(),
                )?;
                if fwd_variants.is_empty() || rev_variants.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "PCR primers must not be empty".to_string(),
                    });
                }

                let mut candidates: Vec<(usize, usize, String)> = vec![];
                let mut seen_amplicons: HashSet<String> = HashSet::new();
                for fwd_full in &fwd_variants {
                    let fwd_anneal_len = forward_primer.anneal_len.unwrap_or(fwd_full.len());
                    if fwd_anneal_len == 0 || fwd_anneal_len > fwd_full.len() {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: "Invalid anneal_len in PCR primer spec".to_string(),
                        });
                    }
                    let fwd_anneal = &fwd_full[fwd_full.len() - fwd_anneal_len..];
                    let fwd_sites = Self::find_anneal_sites(
                        template_bytes,
                        fwd_anneal.as_bytes(),
                        forward_primer.max_mismatches.unwrap_or(0),
                        forward_primer.require_3prime_exact_bases.unwrap_or(0),
                        true,
                    );
                    if fwd_sites.is_empty() {
                        continue;
                    }

                    for rev_full in &rev_variants {
                        let rev_anneal_len = reverse_primer.anneal_len.unwrap_or(rev_full.len());
                        if rev_anneal_len == 0 || rev_anneal_len > rev_full.len() {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "Invalid anneal_len in PCR primer spec".to_string(),
                            });
                        }
                        let rev_anneal = &rev_full[rev_full.len() - rev_anneal_len..];
                        let rev_binding = Self::reverse_complement(rev_anneal);
                        let rev_sites = Self::find_anneal_sites(
                            template_bytes,
                            rev_binding.as_bytes(),
                            reverse_primer.max_mismatches.unwrap_or(0),
                            reverse_primer.require_3prime_exact_bases.unwrap_or(0),
                            false,
                        );
                        if rev_sites.is_empty() {
                            continue;
                        }
                        let rev_full_rc = Self::reverse_complement(rev_full);
                        for fwd_pos in &fwd_sites {
                            for rev_pos in &rev_sites {
                                let interior_start = *fwd_pos + fwd_anneal_len;
                                let interior_end = *rev_pos;
                                if interior_start > interior_end {
                                    continue;
                                }
                                let interior = &template_seq[interior_start..interior_end];
                                let amplicon = format!("{fwd_full}{interior}{rev_full_rc}");
                                if seen_amplicons.insert(amplicon.clone()) {
                                    candidates.push((*fwd_pos, *rev_pos, amplicon));
                                    if candidates.len() > self.max_fragments_per_container() {
                                        return Err(EngineError {
                                            code: ErrorCode::InvalidInput,
                                            message: format!(
                                                "PCR produced more than max_fragments_per_container={}",
                                                self.max_fragments_per_container()
                                            ),
                                        });
                                    }
                                }
                            }
                        }
                    }
                }
                if candidates.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "No valid advanced PCR amplicon could be formed".to_string(),
                    });
                }

                let require_unique = unique.unwrap_or(false);
                if require_unique && candidates.len() != 1 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "PCR unique=true requires exactly one amplicon, found {}",
                            candidates.len()
                        ),
                    });
                }
                if output_id.is_some() && candidates.len() != 1 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message:
                            "PCR output_id can only be used when exactly one amplicon is produced"
                                .to_string(),
                    });
                }

                let default_base = format!("{template}_pcr");
                for (i, (fwd_pos, rev_pos, amplicon)) in candidates.into_iter().enumerate() {
                    let mut pcr_product =
                        DNAsequence::from_sequence(&amplicon).map_err(|e| EngineError {
                            code: ErrorCode::Internal,
                            message: format!("Could not create PCR product: {e}"),
                        })?;
                    pcr_product.set_circular(false);
                    Self::prepare_sequence(&mut pcr_product);

                    let requested = if i == 0 { output_id.clone() } else { None };
                    let seq_id = if let Some(id) = requested {
                        self.unique_seq_id(&id)
                    } else if seen_amplicons.len() == 1 {
                        self.unique_seq_id(&default_base)
                    } else {
                        self.unique_seq_id(&format!("{}_{}", default_base, i + 1))
                    };
                    self.state.sequences.insert(seq_id.clone(), pcr_product);
                    self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
                    result.created_seq_ids.push(seq_id.clone());
                    result.messages.push(format!(
                        "Advanced PCR product '{}' created from fwd@{} rev@{} (len {})",
                        seq_id,
                        fwd_pos,
                        rev_pos,
                        amplicon.len()
                    ));
                }
            }
            Operation::PcrMutagenesis {
                template,
                forward_primer,
                reverse_primer,
                mutations,
                output_id,
                unique,
                require_all_mutations,
            } => {
                parent_seq_ids.push(template.clone());
                let dna = self
                    .state
                    .sequences
                    .get(&template)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{template}' not found"),
                    })?
                    .clone();

                if dna.is_circular() {
                    return Err(EngineError {
                        code: ErrorCode::Unsupported,
                        message: "PCR on circular templates is not implemented yet".to_string(),
                    });
                }
                if mutations.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "PcrMutagenesis requires at least one mutation".to_string(),
                    });
                }

                let template_seq = dna.get_forward_string().to_ascii_uppercase();
                let template_bytes = template_seq.as_bytes();

                let mut normalized_mutations: Vec<(usize, u8, u8)> = vec![];
                for m in &mutations {
                    if m.zero_based_position >= template_bytes.len() {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!(
                                "Mutation position {} is out of bounds for template length {}",
                                m.zero_based_position,
                                template_bytes.len()
                            ),
                        });
                    }
                    let ref_nt = Self::normalize_dna_text(&m.reference).to_ascii_uppercase();
                    let alt_nt = Self::normalize_dna_text(&m.alternate).to_ascii_uppercase();
                    if ref_nt.len() != 1 || alt_nt.len() != 1 {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: "Mutation reference/alternate must be single nucleotides"
                                .to_string(),
                        });
                    }
                    let ref_b = ref_nt.as_bytes()[0];
                    let alt_b = alt_nt.as_bytes()[0];
                    if template_bytes[m.zero_based_position] != ref_b {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!(
                                "Mutation reference mismatch at position {}: template has '{}', expected '{}'",
                                m.zero_based_position,
                                template_bytes[m.zero_based_position] as char,
                                ref_b as char
                            ),
                        });
                    }
                    normalized_mutations.push((m.zero_based_position, ref_b, alt_b));
                }

                let fwd_variants = Self::expand_primer_variants(
                    &forward_primer,
                    self.max_fragments_per_container(),
                )?;
                let rev_variants = Self::expand_primer_variants(
                    &reverse_primer,
                    self.max_fragments_per_container(),
                )?;
                if fwd_variants.is_empty() || rev_variants.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "PCR primers must not be empty".to_string(),
                    });
                }

                let require_all = require_all_mutations.unwrap_or(true);
                let mut selected: Vec<((usize, usize), String)> = vec![];
                let mut seen_amplicons: HashSet<String> = HashSet::new();
                for fwd_full in &fwd_variants {
                    let fwd_anneal_len = forward_primer.anneal_len.unwrap_or(fwd_full.len());
                    if fwd_anneal_len == 0 || fwd_anneal_len > fwd_full.len() {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: "Invalid anneal_len in PCR primer spec".to_string(),
                        });
                    }
                    let fwd_anneal = &fwd_full[fwd_full.len() - fwd_anneal_len..];
                    let fwd_sites = Self::find_anneal_sites(
                        template_bytes,
                        fwd_anneal.as_bytes(),
                        forward_primer.max_mismatches.unwrap_or(0),
                        forward_primer.require_3prime_exact_bases.unwrap_or(0),
                        true,
                    );
                    if fwd_sites.is_empty() {
                        continue;
                    }

                    for rev_full in &rev_variants {
                        let rev_anneal_len = reverse_primer.anneal_len.unwrap_or(rev_full.len());
                        if rev_anneal_len == 0 || rev_anneal_len > rev_full.len() {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "Invalid anneal_len in PCR primer spec".to_string(),
                            });
                        }
                        let rev_anneal = &rev_full[rev_full.len() - rev_anneal_len..];
                        let rev_binding = Self::reverse_complement(rev_anneal);
                        let rev_sites = Self::find_anneal_sites(
                            template_bytes,
                            rev_binding.as_bytes(),
                            reverse_primer.max_mismatches.unwrap_or(0),
                            reverse_primer.require_3prime_exact_bases.unwrap_or(0),
                            false,
                        );
                        if rev_sites.is_empty() {
                            continue;
                        }
                        let rev_full_rc = Self::reverse_complement(rev_full);

                        for fwd_pos in &fwd_sites {
                            for rev_pos in &rev_sites {
                                let interior_start = *fwd_pos + fwd_anneal_len;
                                let interior_end = *rev_pos;
                                if interior_start > interior_end {
                                    continue;
                                }
                                let interior = &template_seq[interior_start..interior_end];
                                let amplicon = format!("{fwd_full}{interior}{rev_full_rc}");
                                let mut introduced_count = 0usize;
                                let mut valid = true;
                                for (pos, _ref_b, alt_b) in &normalized_mutations {
                                    if *pos < *fwd_pos || *pos >= (*rev_pos + rev_anneal_len) {
                                        if require_all {
                                            valid = false;
                                            break;
                                        }
                                        continue;
                                    }

                                    let observed = if *pos < (*fwd_pos + fwd_anneal_len) {
                                        let offset =
                                            fwd_full.len() - fwd_anneal_len + (*pos - *fwd_pos);
                                        fwd_full.as_bytes()[offset]
                                    } else if *pos < *rev_pos {
                                        template_bytes[*pos]
                                    } else {
                                        let offset = *pos - *rev_pos;
                                        rev_full_rc.as_bytes()[offset]
                                    };

                                    if observed == *alt_b {
                                        introduced_count += 1;
                                    } else if require_all {
                                        valid = false;
                                        break;
                                    }
                                }

                                let keep = if require_all {
                                    valid && introduced_count == normalized_mutations.len()
                                } else {
                                    introduced_count > 0
                                };
                                if keep && seen_amplicons.insert(amplicon.clone()) {
                                    selected.push(((*fwd_pos, *rev_pos), amplicon));
                                    if selected.len() > self.max_fragments_per_container() {
                                        return Err(EngineError {
                                            code: ErrorCode::InvalidInput,
                                            message: format!(
                                                "Mutagenesis PCR produced more than max_fragments_per_container={}",
                                                self.max_fragments_per_container()
                                            ),
                                        });
                                    }
                                }
                            }
                        }
                    }
                }

                if selected.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: if require_all {
                            "No amplicon introduced all requested mutations".to_string()
                        } else {
                            "No amplicon introduced any requested mutation".to_string()
                        },
                    });
                }
                if selected.len() > self.max_fragments_per_container() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Mutagenesis PCR produced {} amplicons, exceeding max_fragments_per_container={}",
                            selected.len(),
                            self.max_fragments_per_container()
                        ),
                    });
                }

                let require_unique = unique.unwrap_or(false);
                if require_unique && selected.len() != 1 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "PCR unique=true requires exactly one amplicon, found {}",
                            selected.len()
                        ),
                    });
                }
                if output_id.is_some() && selected.len() != 1 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message:
                            "PCR output_id can only be used when exactly one amplicon is produced"
                                .to_string(),
                    });
                }

                let default_base = format!("{template}_pcr_mut");
                for (i, ((fwd_pos, rev_pos), amplicon)) in selected.into_iter().enumerate() {
                    let mut pcr_product =
                        DNAsequence::from_sequence(&amplicon).map_err(|e| EngineError {
                            code: ErrorCode::Internal,
                            message: format!("Could not create PCR product: {e}"),
                        })?;
                    pcr_product.set_circular(false);
                    Self::prepare_sequence(&mut pcr_product);

                    let requested = if i == 0 { output_id.clone() } else { None };
                    let seq_id = if let Some(id) = requested {
                        self.unique_seq_id(&id)
                    } else if seen_amplicons.len() == 1 {
                        self.unique_seq_id(&default_base)
                    } else {
                        self.unique_seq_id(&format!("{}_{}", default_base, i + 1))
                    };
                    self.state.sequences.insert(seq_id.clone(), pcr_product);
                    self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
                    result.created_seq_ids.push(seq_id.clone());
                    result.messages.push(format!(
                        "Mutagenesis PCR product '{}' created from fwd@{} rev@{}",
                        seq_id, fwd_pos, rev_pos
                    ));
                }
            }
            Operation::ExtractRegion {
                input,
                from,
                to,
                output_id,
            } => {
                parent_seq_ids.push(input.clone());
                let dna = self
                    .state
                    .sequences
                    .get(&input)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{input}' not found"),
                    })?
                    .clone();

                if from == to {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "ExtractRegion requires from != to".to_string(),
                    });
                }
                let fragment = dna.get_range_safe(from..to).ok_or_else(|| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Could not extract region {}..{} from sequence '{}'",
                        from, to, input
                    ),
                })?;
                let text = String::from_utf8_lossy(&fragment).to_string();
                let mut out = DNAsequence::from_sequence(&text).map_err(|e| EngineError {
                    code: ErrorCode::Internal,
                    message: format!("Could not create extracted sequence: {e}"),
                })?;
                out.set_circular(false);
                Self::prepare_sequence(&mut out);

                let base = output_id.unwrap_or_else(|| format!("{input}_region"));
                let seq_id = self.unique_seq_id(&base);
                self.state.sequences.insert(seq_id.clone(), out);
                self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
                result.created_seq_ids.push(seq_id.clone());
                result.messages.push(format!(
                    "Extracted region {}..{} from '{}' into '{}'",
                    from, to, input, seq_id
                ));
            }
            Operation::ExtractAnchoredRegion {
                input,
                anchor,
                direction,
                target_length_bp,
                length_tolerance_bp,
                required_re_sites,
                required_tf_motifs,
                forward_primer,
                reverse_primer,
                output_prefix,
                unique,
                max_candidates,
            } => {
                parent_seq_ids.push(input.clone());
                let dna = self
                    .state
                    .sequences
                    .get(&input)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{input}' not found"),
                    })?
                    .clone();

                if target_length_bp == 0 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "ExtractAnchoredRegion requires target_length_bp >= 1".to_string(),
                    });
                }

                let anchor_pos = Self::resolve_anchor_position(&dna, &anchor)?;
                let min_len = target_length_bp.saturating_sub(length_tolerance_bp).max(1);
                let max_len = target_length_bp
                    .saturating_add(length_tolerance_bp)
                    .max(min_len);

                let forward_primer = match forward_primer {
                    Some(p) => {
                        let v = Self::normalize_iupac_text(&p)?;
                        if v.is_empty() {
                            None
                        } else {
                            Some(v)
                        }
                    }
                    None => None,
                };
                let reverse_primer = match reverse_primer {
                    Some(p) => {
                        let v = Self::normalize_iupac_text(&p)?;
                        if v.is_empty() {
                            None
                        } else {
                            Some(v)
                        }
                    }
                    None => None,
                };
                let reverse_primer_rc = reverse_primer
                    .as_ref()
                    .map(|p| Self::reverse_complement_iupac(p))
                    .transpose()?;

                let mut tf_motifs = Vec::new();
                for motif in required_tf_motifs {
                    let m = Self::resolve_tf_motif_or_iupac(&motif)?;
                    if !m.is_empty() {
                        tf_motifs.push(m);
                    }
                }

                let (required_enzymes, missing_enzymes) = if required_re_sites.is_empty() {
                    (Vec::new(), Vec::new())
                } else {
                    self.resolve_enzymes(&required_re_sites)?
                };
                if !missing_enzymes.is_empty() {
                    result.warnings.push(format!(
                        "Unknown anchored-region enzymes ignored: {}",
                        missing_enzymes.join(",")
                    ));
                }

                #[derive(Clone)]
                struct Candidate {
                    start: usize,
                    end: usize,
                    len: usize,
                    sequence: String,
                    score: usize,
                }

                let mut candidates = Vec::new();
                for len in min_len..=max_len {
                    let Some((start, end)) = Self::anchored_range(
                        anchor_pos,
                        len,
                        &direction,
                        dna.len(),
                        dna.is_circular(),
                    ) else {
                        continue;
                    };
                    let Some(fragment) = dna.get_range_safe(start..end) else {
                        continue;
                    };
                    if fragment.is_empty() {
                        continue;
                    }
                    let fragment_text = String::from_utf8_lossy(&fragment).to_string();
                    let fragment_bytes = fragment_text.as_bytes();

                    if let Some(fwd) = &forward_primer {
                        if !Self::iupac_match_at(fragment_bytes, fwd.as_bytes(), 0) {
                            continue;
                        }
                    }
                    if let Some(rev_rc) = &reverse_primer_rc {
                        if rev_rc.len() > fragment_bytes.len() {
                            continue;
                        }
                        let start_idx = fragment_bytes.len() - rev_rc.len();
                        if !Self::iupac_match_at(fragment_bytes, rev_rc.as_bytes(), start_idx) {
                            continue;
                        }
                    }

                    let mut motif_ok = true;
                    for motif in &tf_motifs {
                        if !Self::contains_motif_any_strand(fragment_bytes, motif)? {
                            motif_ok = false;
                            break;
                        }
                    }
                    if !motif_ok {
                        continue;
                    }

                    if !required_enzymes.is_empty() {
                        let frag_dna = DNAsequence::from_sequence(&fragment_text).map_err(|e| {
                            EngineError {
                                code: ErrorCode::Internal,
                                message: format!(
                                    "Could not evaluate anchored-region enzyme constraints: {e}"
                                ),
                            }
                        })?;
                        let mut enzymes_ok = true;
                        for enzyme in &required_enzymes {
                            if enzyme.get_sites(&frag_dna, None).is_empty() {
                                enzymes_ok = false;
                                break;
                            }
                        }
                        if !enzymes_ok {
                            continue;
                        }
                    }

                    candidates.push(Candidate {
                        start,
                        end,
                        len,
                        sequence: fragment_text,
                        score: len.abs_diff(target_length_bp),
                    });
                }

                if candidates.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message:
                            "No anchored-region candidate satisfied the configured constraints"
                                .to_string(),
                    });
                }

                candidates.sort_by(|a, b| {
                    a.score
                        .cmp(&b.score)
                        .then(a.start.cmp(&b.start))
                        .then(a.end.cmp(&b.end))
                });

                let require_unique = unique.unwrap_or(false);
                if require_unique && candidates.len() != 1 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "ExtractAnchoredRegion unique=true requires exactly one candidate, found {}",
                            candidates.len()
                        ),
                    });
                }

                let limit = max_candidates
                    .unwrap_or(candidates.len())
                    .min(self.max_fragments_per_container());
                if limit == 0 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "ExtractAnchoredRegion max_candidates must be >= 1".to_string(),
                    });
                }
                if candidates.len() > limit {
                    result.warnings.push(format!(
                        "Anchored-region candidates truncated from {} to {} by max_candidates/max_fragments_per_container",
                        candidates.len(),
                        limit
                    ));
                    candidates.truncate(limit);
                }

                let prefix = output_prefix.unwrap_or_else(|| format!("{input}_anchored"));
                for (idx, cand) in candidates.into_iter().enumerate() {
                    let mut out =
                        DNAsequence::from_sequence(&cand.sequence).map_err(|e| EngineError {
                            code: ErrorCode::Internal,
                            message: format!("Could not create anchored-region sequence: {e}"),
                        })?;
                    out.set_circular(false);
                    Self::prepare_sequence(&mut out);
                    let seq_id = self.unique_seq_id(&format!("{}_{}", prefix, idx + 1));
                    self.state.sequences.insert(seq_id.clone(), out);
                    self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
                    result.created_seq_ids.push(seq_id.clone());
                    result.messages.push(format!(
                        "Anchored-region candidate '{}' [{}..{}, {} bp] score={}",
                        seq_id, cand.start, cand.end, cand.len, cand.score
                    ));
                }
            }
            Operation::SelectCandidate {
                input,
                criterion,
                output_id,
            } => {
                parent_seq_ids.push(input.clone());
                let dna = self
                    .state
                    .sequences
                    .get(&input)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{input}' not found"),
                    })?
                    .clone();

                let base = output_id.unwrap_or_else(|| format!("{input}_selected"));
                let seq_id = self.unique_seq_id(&base);
                self.state.sequences.insert(seq_id.clone(), dna);
                self.add_lineage_node(
                    &seq_id,
                    SequenceOrigin::InSilicoSelection,
                    Some(&result.op_id),
                );
                result.created_seq_ids.push(seq_id.clone());
                result.warnings.push(
                    "Selection operation is in-silico and may not directly correspond to a unique wet-lab product"
                        .to_string(),
                );
                result.messages.push(format!(
                    "Selected candidate '{}' from '{}' using criterion '{}'",
                    seq_id, input, criterion
                ));
            }
            Operation::FilterByMolecularWeight {
                inputs,
                min_bp,
                max_bp,
                error,
                unique,
                output_prefix,
            } => {
                if inputs.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "FilterByMolecularWeight requires at least one input sequence"
                            .to_string(),
                    });
                }
                if min_bp > max_bp {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("min_bp ({min_bp}) must be <= max_bp ({max_bp})"),
                    });
                }
                if !(0.0..=1.0).contains(&error) {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "error must be between 0.0 and 1.0".to_string(),
                    });
                }

                let min_allowed = ((min_bp as f64) * (1.0 - error)).floor() as usize;
                let max_allowed = ((max_bp as f64) * (1.0 + error)).ceil() as usize;

                let mut matches: Vec<(SeqId, DNAsequence)> = vec![];
                for input in &inputs {
                    parent_seq_ids.push(input.clone());
                    let dna = self
                        .state
                        .sequences
                        .get(input)
                        .ok_or_else(|| EngineError {
                            code: ErrorCode::NotFound,
                            message: format!("Sequence '{input}' not found"),
                        })?
                        .clone();
                    let bp = dna.len();
                    if bp >= min_allowed && bp <= max_allowed {
                        matches.push((input.clone(), dna));
                    }
                }

                if matches.len() > self.max_fragments_per_container() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "FilterByMolecularWeight produced {} candidates, exceeding max_fragments_per_container={}",
                            matches.len(),
                            self.max_fragments_per_container()
                        ),
                    });
                }

                if unique && matches.len() != 1 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "unique=true requires exactly one match, found {}",
                            matches.len()
                        ),
                    });
                }

                let prefix = output_prefix.unwrap_or_else(|| "mw_filter".to_string());
                for (i, (_source_id, mut dna)) in matches.into_iter().enumerate() {
                    Self::prepare_sequence(&mut dna);
                    let candidate = format!("{}_{}", prefix, i + 1);
                    let seq_id = self.unique_seq_id(&candidate);
                    self.state.sequences.insert(seq_id.clone(), dna);
                    self.add_lineage_node(
                        &seq_id,
                        SequenceOrigin::InSilicoSelection,
                        Some(&result.op_id),
                    );
                    result.created_seq_ids.push(seq_id);
                }

                result.messages.push(format!(
                    "Molecular-weight filter kept {} sequence(s) in effective range {}-{} bp (requested {}-{} bp, error={:.3})",
                    result.created_seq_ids.len(),
                    min_allowed,
                    max_allowed,
                    min_bp,
                    max_bp,
                    error
                ));
            }
            Operation::FilterByDesignConstraints {
                inputs,
                gc_min,
                gc_max,
                max_homopolymer_run,
                reject_ambiguous_bases,
                avoid_u6_terminator_tttt,
                forbidden_motifs,
                unique,
                output_prefix,
            } => {
                if inputs.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "FilterByDesignConstraints requires at least one input sequence"
                            .to_string(),
                    });
                }

                if let Some(min) = gc_min {
                    if !(0.0..=1.0).contains(&min) {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!("gc_min ({min}) must be between 0.0 and 1.0"),
                        });
                    }
                }
                if let Some(max) = gc_max {
                    if !(0.0..=1.0).contains(&max) {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!("gc_max ({max}) must be between 0.0 and 1.0"),
                        });
                    }
                }
                if let (Some(min), Some(max)) = (gc_min, gc_max) {
                    if min > max {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!("gc_min ({min}) must be <= gc_max ({max})"),
                        });
                    }
                }
                if let Some(max_run) = max_homopolymer_run {
                    if max_run == 0 {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: "max_homopolymer_run must be >= 1".to_string(),
                        });
                    }
                }

                let reject_ambiguous_bases = reject_ambiguous_bases.unwrap_or(true);
                let avoid_u6_terminator_tttt = avoid_u6_terminator_tttt.unwrap_or(true);
                let mut forbidden_motifs_normalized: Vec<String> = vec![];
                for motif in forbidden_motifs {
                    let normalized = Self::normalize_iupac_text(&motif)?;
                    if !normalized.is_empty() {
                        forbidden_motifs_normalized.push(normalized);
                    }
                }

                let mut matches: Vec<(SeqId, DNAsequence)> = vec![];
                let mut rejected = 0usize;
                let mut rejection_warnings_left = 32usize;
                let input_count = inputs.len();

                for input in &inputs {
                    parent_seq_ids.push(input.clone());
                    let dna = self
                        .state
                        .sequences
                        .get(input)
                        .ok_or_else(|| EngineError {
                            code: ErrorCode::NotFound,
                            message: format!("Sequence '{input}' not found"),
                        })?
                        .clone();
                    let sequence = Self::normalized_sequence_for_quality(&dna);
                    let mut reasons: Vec<String> = vec![];

                    if reject_ambiguous_bases && Self::has_ambiguous_bases(&sequence) {
                        reasons.push("contains_ambiguous_base".to_string());
                    }

                    if gc_min.is_some() || gc_max.is_some() {
                        match Self::sequence_gc_fraction(&sequence) {
                            Some(gc) => {
                                if let Some(min) = gc_min {
                                    if gc < min {
                                        reasons.push(format!("gc_too_low({gc:.3}<{min:.3})"));
                                    }
                                }
                                if let Some(max) = gc_max {
                                    if gc > max {
                                        reasons.push(format!("gc_too_high({gc:.3}>{max:.3})"));
                                    }
                                }
                            }
                            None => reasons.push("gc_not_computable".to_string()),
                        }
                    }

                    if let Some(max_run) = max_homopolymer_run {
                        let observed = Self::max_homopolymer_run(&sequence);
                        if observed > max_run {
                            reasons.push(format!("homopolymer_run_exceeded({observed}>{max_run})"));
                        }
                    }

                    if avoid_u6_terminator_tttt && Self::contains_u6_terminator_t4(&sequence) {
                        reasons.push("u6_terminator_t4".to_string());
                    }

                    for motif in &forbidden_motifs_normalized {
                        if Self::contains_motif_any_strand(&sequence, motif)? {
                            reasons.push(format!("forbidden_motif_present({motif})"));
                        }
                    }

                    if reasons.is_empty() {
                        matches.push((input.clone(), dna));
                    } else {
                        rejected += 1;
                        if rejection_warnings_left > 0 {
                            result.warnings.push(format!(
                                "Sequence '{}' rejected by design constraints: {}",
                                input,
                                reasons.join(", ")
                            ));
                            rejection_warnings_left -= 1;
                        }
                    }
                }

                if matches.len() > self.max_fragments_per_container() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "FilterByDesignConstraints produced {} candidates, exceeding max_fragments_per_container={}",
                            matches.len(),
                            self.max_fragments_per_container()
                        ),
                    });
                }

                if unique && matches.len() != 1 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "unique=true requires exactly one match, found {}",
                            matches.len()
                        ),
                    });
                }

                let prefix = output_prefix.unwrap_or_else(|| "design_filter".to_string());
                for (i, (_source_id, mut dna)) in matches.into_iter().enumerate() {
                    Self::prepare_sequence(&mut dna);
                    let candidate = format!("{}_{}", prefix, i + 1);
                    let seq_id = self.unique_seq_id(&candidate);
                    self.state.sequences.insert(seq_id.clone(), dna);
                    self.add_lineage_node(
                        &seq_id,
                        SequenceOrigin::InSilicoSelection,
                        Some(&result.op_id),
                    );
                    result.created_seq_ids.push(seq_id);
                }

                if rejected > 32 {
                    result.warnings.push(format!(
                        "{} additional sequence(s) were rejected (warning output truncated)",
                        rejected - 32
                    ));
                }

                result.messages.push(format!(
                    "Design-constraint filter kept {} of {} sequence(s)",
                    result.created_seq_ids.len(),
                    input_count
                ));
                result.messages.push(format!(
                    "Applied filters: gc_min={:?}, gc_max={:?}, max_homopolymer_run={:?}, reject_ambiguous_bases={}, avoid_u6_terminator_tttt={}, forbidden_motifs={}",
                    gc_min,
                    gc_max,
                    max_homopolymer_run,
                    reject_ambiguous_bases,
                    avoid_u6_terminator_tttt,
                    forbidden_motifs_normalized.len()
                ));
            }
            Operation::GenerateCandidateSet {
                set_name,
                seq_id,
                length_bp,
                step_bp,
                feature_kinds,
                feature_label_regex,
                max_distance_bp,
                limit,
            } => {
                self.op_generate_candidate_set(
                    set_name,
                    seq_id,
                    length_bp,
                    step_bp,
                    feature_kinds,
                    feature_label_regex,
                    max_distance_bp,
                    limit,
                    &mut result,
                )?;
            }
            Operation::DeleteCandidateSet { set_name } => {
                self.op_delete_candidate_set(set_name, &mut result)?;
            }
            Operation::ScoreCandidateSetExpression {
                set_name,
                metric,
                expression,
            } => {
                self.op_score_candidate_set_expression(set_name, metric, expression, &mut result)?;
            }
            Operation::ScoreCandidateSetDistance {
                set_name,
                metric,
                feature_kinds,
                feature_label_regex,
            } => {
                self.op_score_candidate_set_distance(
                    set_name,
                    metric,
                    feature_kinds,
                    feature_label_regex,
                    &mut result,
                )?;
            }
            Operation::FilterCandidateSet {
                input_set,
                output_set,
                metric,
                min,
                max,
                min_quantile,
                max_quantile,
            } => {
                self.op_filter_candidate_set(
                    input_set,
                    output_set,
                    metric,
                    min,
                    max,
                    min_quantile,
                    max_quantile,
                    &mut result,
                )?;
            }
            Operation::CandidateSetOp {
                op,
                left_set,
                right_set,
                output_set,
            } => {
                self.op_candidate_set_op(op, left_set, right_set, output_set, &mut result)?;
            }
            Operation::Reverse { input, output_id } => {
                parent_seq_ids.push(input.clone());
                let dna = self
                    .state
                    .sequences
                    .get(&input)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{input}' not found"),
                    })?
                    .clone();
                let mut text = dna.get_forward_string();
                text = text.chars().rev().collect();
                let mut out = DNAsequence::from_sequence(&text).map_err(|e| EngineError {
                    code: ErrorCode::Internal,
                    message: format!("Could not create reverse sequence: {e}"),
                })?;
                out.set_circular(dna.is_circular());
                Self::prepare_sequence(&mut out);

                let base = output_id.unwrap_or_else(|| format!("{input}_rev"));
                let seq_id = self.unique_seq_id(&base);
                self.state.sequences.insert(seq_id.clone(), out);
                self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
                result.created_seq_ids.push(seq_id.clone());
                result.messages.push(format!(
                    "Created reverse sequence '{}' from '{}'",
                    seq_id, input
                ));
            }
            Operation::Complement { input, output_id } => {
                parent_seq_ids.push(input.clone());
                let dna = self
                    .state
                    .sequences
                    .get(&input)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{input}' not found"),
                    })?
                    .clone();
                let text: String = dna
                    .get_forward_string()
                    .as_bytes()
                    .iter()
                    .map(|c| IupacCode::letter_complement(*c))
                    .map(char::from)
                    .collect();
                let mut out = DNAsequence::from_sequence(&text).map_err(|e| EngineError {
                    code: ErrorCode::Internal,
                    message: format!("Could not create complement sequence: {e}"),
                })?;
                out.set_circular(dna.is_circular());
                Self::prepare_sequence(&mut out);

                let base = output_id.unwrap_or_else(|| format!("{input}_comp"));
                let seq_id = self.unique_seq_id(&base);
                self.state.sequences.insert(seq_id.clone(), out);
                self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
                result.created_seq_ids.push(seq_id.clone());
                result.messages.push(format!(
                    "Created complement sequence '{}' from '{}'",
                    seq_id, input
                ));
            }
            Operation::ReverseComplement { input, output_id } => {
                parent_seq_ids.push(input.clone());
                let dna = self
                    .state
                    .sequences
                    .get(&input)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{input}' not found"),
                    })?
                    .clone();
                let text = Self::reverse_complement(&dna.get_forward_string());
                let mut out = DNAsequence::from_sequence(&text).map_err(|e| EngineError {
                    code: ErrorCode::Internal,
                    message: format!("Could not create reverse-complement sequence: {e}"),
                })?;
                out.set_circular(dna.is_circular());
                Self::prepare_sequence(&mut out);

                let base = output_id.unwrap_or_else(|| format!("{input}_revcomp"));
                let seq_id = self.unique_seq_id(&base);
                self.state.sequences.insert(seq_id.clone(), out);
                self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
                result.created_seq_ids.push(seq_id.clone());
                result.messages.push(format!(
                    "Created reverse-complement sequence '{}' from '{}'",
                    seq_id, input
                ));
            }
            Operation::Branch { input, output_id } => {
                parent_seq_ids.push(input.clone());
                let dna = self
                    .state
                    .sequences
                    .get(&input)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{input}' not found"),
                    })?
                    .clone();

                let base = output_id.unwrap_or_else(|| format!("{input}_branch"));
                let seq_id = self.unique_seq_id(&base);
                self.state.sequences.insert(seq_id.clone(), dna);
                self.add_lineage_node(&seq_id, SequenceOrigin::Branch, Some(&result.op_id));
                result.created_seq_ids.push(seq_id.clone());
                result
                    .messages
                    .push(format!("Branched '{}' into '{}'", input, seq_id));
            }
            Operation::SetDisplayVisibility { target, visible } => {
                let (name, slot): (&str, &mut bool) = match target {
                    DisplayTarget::SequencePanel => (
                        "sequence_panel",
                        &mut self.state.display.show_sequence_panel,
                    ),
                    DisplayTarget::MapPanel => {
                        ("map_panel", &mut self.state.display.show_map_panel)
                    }
                    DisplayTarget::Features => ("features", &mut self.state.display.show_features),
                    DisplayTarget::CdsFeatures => {
                        ("cds_features", &mut self.state.display.show_cds_features)
                    }
                    DisplayTarget::GeneFeatures => {
                        ("gene_features", &mut self.state.display.show_gene_features)
                    }
                    DisplayTarget::MrnaFeatures => {
                        ("mrna_features", &mut self.state.display.show_mrna_features)
                    }
                    DisplayTarget::Tfbs => ("tfbs", &mut self.state.display.show_tfbs),
                    DisplayTarget::RestrictionEnzymes => (
                        "restriction_enzymes",
                        &mut self.state.display.show_restriction_enzymes,
                    ),
                    DisplayTarget::GcContents => {
                        ("gc_contents", &mut self.state.display.show_gc_contents)
                    }
                    DisplayTarget::OpenReadingFrames => (
                        "open_reading_frames",
                        &mut self.state.display.show_open_reading_frames,
                    ),
                    DisplayTarget::MethylationSites => (
                        "methylation_sites",
                        &mut self.state.display.show_methylation_sites,
                    ),
                };
                *slot = visible;
                result
                    .messages
                    .push(format!("Set display target '{name}' to {visible}"));
            }
            Operation::SetLinearViewport { start_bp, span_bp } => {
                self.state.display.linear_view_start_bp = start_bp;
                self.state.display.linear_view_span_bp = span_bp;
                result.messages.push(format!(
                    "Set linear viewport start_bp={start_bp}, span_bp={span_bp}"
                ));
            }
            Operation::SetTopology { seq_id, circular } => {
                let _ = self.ensure_lineage_node(&seq_id);
                let dna = self
                    .state
                    .sequences
                    .get_mut(&seq_id)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{seq_id}' not found"),
                    })?;
                dna.set_circular(circular);
                dna.update_computed_features();
                result.changed_seq_ids.push(seq_id.clone());
                result.messages.push(format!(
                    "Set topology of '{seq_id}' to {}",
                    if circular { "circular" } else { "linear" }
                ));
            }
            Operation::RecomputeFeatures { seq_id } => {
                let _ = self.ensure_lineage_node(&seq_id);
                let dna = self
                    .state
                    .sequences
                    .get_mut(&seq_id)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{seq_id}' not found"),
                    })?;
                dna.update_computed_features();
                result.changed_seq_ids.push(seq_id.clone());
                result
                    .messages
                    .push(format!("Recomputed features for '{seq_id}'"));
            }
            Operation::AnnotateTfbs {
                seq_id,
                motifs,
                min_llr_bits,
                min_llr_quantile,
                per_tf_thresholds,
                clear_existing,
                max_hits,
            } => {
                const DEFAULT_MAX_TFBS_HITS: usize = 500;
                let motifs = if motifs.len() == 1
                    && matches!(motifs[0].trim().to_ascii_uppercase().as_str(), "ALL" | "*")
                {
                    tf_motifs::all_motif_ids()
                } else {
                    motifs
                };

                if motifs.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "AnnotateTfbs requires at least one motif".to_string(),
                    });
                }
                let default_min_llr_bits = min_llr_bits.unwrap_or(f64::NEG_INFINITY);
                let default_min_llr_quantile = min_llr_quantile.unwrap_or(0.0);
                Self::validate_tf_thresholds(default_min_llr_quantile)?;

                let mut override_map: HashMap<String, (Option<f64>, Option<f64>)> = HashMap::new();
                for o in &per_tf_thresholds {
                    let key = o.tf.trim().to_ascii_uppercase();
                    if key.is_empty() {
                        continue;
                    }
                    if let Some(q) = o.min_llr_quantile {
                        Self::validate_tf_thresholds(q)?;
                    }
                    override_map.insert(key, (o.min_llr_bits, o.min_llr_quantile));
                }

                let _ = self.ensure_lineage_node(&seq_id);
                let dna = self
                    .state
                    .sequences
                    .get_mut(&seq_id)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{seq_id}' not found"),
                    })?;
                let seq_text = dna.get_forward_string();
                let seq_bytes = seq_text.as_bytes();

                if clear_existing.unwrap_or(true) {
                    Self::remove_generated_tfbs_features(dna.features_mut());
                }

                let mut added = 0usize;
                let max_hits = match max_hits {
                    Some(0) => None,
                    Some(v) => Some(v),
                    None => Some(DEFAULT_MAX_TFBS_HITS),
                };
                let motif_count = motifs.len();
                let mut motifs_scanned = 0usize;
                let mut cap_reached = false;
                'motif_loop: for (motif_idx, token) in motifs.into_iter().enumerate() {
                    let token_key = token.trim().to_ascii_uppercase();
                    if token_key.is_empty() {
                        continue;
                    }
                    let (tf_id, tf_name, _consensus, matrix_counts) =
                        Self::resolve_tf_motif_for_scoring(&token)?;
                    let (llr_matrix, true_log_odds_matrix) =
                        Self::prepare_scoring_matrices(&matrix_counts);
                    if llr_matrix.is_empty() || llr_matrix.len() > seq_bytes.len() {
                        result.warnings.push(format!(
                            "TF '{}' skipped: motif length {} exceeds sequence length {}",
                            tf_id,
                            llr_matrix.len(),
                            seq_bytes.len()
                        ));
                        continue;
                    }

                    let mut eff_bits = default_min_llr_bits;
                    let mut eff_quantile = default_min_llr_quantile;
                    let id_key = tf_id.to_ascii_uppercase();
                    let name_key = tf_name
                        .as_ref()
                        .map(|n| n.trim().to_ascii_uppercase())
                        .unwrap_or_default();
                    for key in [token_key.as_str(), id_key.as_str(), name_key.as_str()] {
                        if key.is_empty() {
                            continue;
                        }
                        if let Some((b, q)) = override_map.get(key) {
                            if let Some(v) = b {
                                eff_bits = *v;
                            }
                            if let Some(v) = q {
                                eff_quantile = *v;
                            }
                            break;
                        }
                    }

                    motifs_scanned += 1;
                    let seq_id_for_progress = seq_id.clone();
                    let tf_id_for_progress = tf_id.clone();
                    let motif_index = motif_idx + 1;
                    let hits = Self::scan_tf_scores(
                        seq_bytes,
                        &llr_matrix,
                        &true_log_odds_matrix,
                        |scanned_steps, total_steps| {
                            let motif_fraction = if total_steps == 0 {
                                1.0
                            } else {
                                (scanned_steps as f64 / total_steps as f64).clamp(0.0, 1.0)
                            };
                            let total_fraction = if motif_count == 0 {
                                1.0
                            } else {
                                ((motif_index - 1) as f64 + motif_fraction) / motif_count as f64
                            }
                            .clamp(0.0, 1.0);
                            on_progress(OperationProgress::Tfbs(TfbsProgress {
                                seq_id: seq_id_for_progress.clone(),
                                motif_id: tf_id_for_progress.clone(),
                                motif_index,
                                motif_count,
                                scanned_steps,
                                total_steps,
                                motif_percent: motif_fraction * 100.0,
                                total_percent: total_fraction * 100.0,
                            }));
                        },
                    );
                    let mut kept = 0usize;
                    for (
                        start,
                        reverse,
                        llr_bits,
                        llr_quantile,
                        true_log_odds_bits,
                        true_log_odds_quantile,
                    ) in hits
                    {
                        if llr_bits < eff_bits || llr_quantile < eff_quantile {
                            continue;
                        }
                        let end = start + llr_matrix.len();
                        dna.features_mut().push(Self::build_tfbs_feature(
                            start,
                            end,
                            reverse,
                            llr_matrix.len(),
                            &tf_id,
                            tf_name.as_deref(),
                            llr_bits,
                            llr_quantile,
                            true_log_odds_bits,
                            true_log_odds_quantile,
                        ));
                        kept += 1;
                        added += 1;
                        if let Some(limit) = max_hits {
                            if added >= limit {
                                cap_reached = true;
                                break;
                            }
                        }
                    }
                    result.messages.push(format!(
                        "TF '{}' annotated {} hit(s){}",
                        tf_id,
                        kept,
                        Self::format_tf_threshold_summary(eff_bits, eff_quantile)
                    ));
                    on_progress(OperationProgress::Tfbs(TfbsProgress {
                        seq_id: seq_id.clone(),
                        motif_id: tf_id,
                        motif_index,
                        motif_count,
                        scanned_steps: 1,
                        total_steps: 1,
                        motif_percent: 100.0,
                        total_percent: (motif_index as f64 / motif_count.max(1) as f64) * 100.0,
                    }));
                    if cap_reached {
                        if let Some(limit) = max_hits {
                            result.warnings.push(format!(
                                "TFBS hit cap ({limit}) reached after scanning {motifs_scanned}/{motif_count} motif(s); skipping remaining motif scans"
                            ));
                        }
                        break 'motif_loop;
                    }
                }
                result.messages.push(format!(
                    "TFBS motif scan coverage: {motifs_scanned}/{motif_count} motif(s)"
                ));

                result.changed_seq_ids.push(seq_id.clone());
                result.messages.push(format!(
                    "Annotated {} TFBS feature(s) on '{}'",
                    added, seq_id
                ));
            }
            Operation::SetParameter { name, value } => {
                match name.as_str() {
                    "max_fragments_per_container" => {
                        let raw = value.as_u64().ok_or_else(|| EngineError {
                            code: ErrorCode::InvalidInput,
                            message: "SetParameter max_fragments_per_container requires a positive integer".to_string(),
                        })?;
                        if raw == 0 {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "max_fragments_per_container must be >= 1".to_string(),
                            });
                        }
                        self.state.parameters.max_fragments_per_container = raw as usize;
                        result.messages.push(format!(
                            "Set parameter '{}' to {}",
                            name, self.state.parameters.max_fragments_per_container
                        ));
                    }
                    "feature_details_font_size" | "feature_detail_font_size" => {
                        let raw = value.as_f64().ok_or_else(|| EngineError {
                            code: ErrorCode::InvalidInput,
                            message: "SetParameter feature_details_font_size requires a number"
                                .to_string(),
                        })?;
                        if !raw.is_finite() {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "feature_details_font_size must be a finite number"
                                    .to_string(),
                            });
                        }
                        if !(8.0..=24.0).contains(&raw) {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "feature_details_font_size must be between 8.0 and 24.0"
                                    .to_string(),
                            });
                        }
                        self.state.display.feature_details_font_size = raw as f32;
                        result.messages.push(format!(
                            "Set parameter 'feature_details_font_size' to {:.2}",
                            self.state.display.feature_details_font_size
                        ));
                    }
                    _ => {
                        return Err(EngineError {
                            code: ErrorCode::Unsupported,
                            message: format!("Unknown parameter '{}'", name),
                        });
                    }
                }
            }
        }

        self.add_lineage_edges(
            &parent_seq_ids,
            &result.created_seq_ids,
            &result.op_id,
            run_id,
        );
        self.add_container_from_result(&op_for_containers, &result);

        Ok(result)
    }
}

impl Engine for GentleEngine {
    fn apply(&mut self, op: Operation) -> Result<OpResult, EngineError> {
        let run_id = "interactive".to_string();
        let mut noop = |_p: OperationProgress| {};
        let result = self.apply_internal(op.clone(), &run_id, &mut noop)?;
        self.journal.push(OperationRecord {
            run_id,
            op,
            result: result.clone(),
        });
        Ok(result)
    }

    fn apply_workflow(&mut self, wf: Workflow) -> Result<Vec<OpResult>, EngineError> {
        let mut results = Vec::new();
        for op in &wf.ops {
            let mut noop = |_p: OperationProgress| {};
            let result = self.apply_internal(op.clone(), &wf.run_id, &mut noop)?;
            self.journal.push(OperationRecord {
                run_id: wf.run_id.clone(),
                op: op.clone(),
                result: result.clone(),
            });
            results.push(result);
        }
        Ok(results)
    }

    fn snapshot(&self) -> &ProjectState {
        &self.state
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bio::io::fasta;
    use flate2::{write::GzEncoder, Compression};
    use std::env;
    use std::fs;
    use std::io::Write;
    #[cfg(unix)]
    use std::os::unix::fs::PermissionsExt;
    use std::path::Path;
    use tempfile::tempdir;

    fn seq(s: &str) -> DNAsequence {
        DNAsequence::from_sequence(s).unwrap()
    }

    fn synth_oligo(desc: &str, sequence: &[u8]) -> DNAsequence {
        let record = fasta::Record::with_attrs("synthetic", Some(desc), sequence);
        DNAsequence::from_fasta_record(&record)
    }

    fn assert_fasta_roundtrip(expected: DNAsequence) {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("synth".to_string(), expected.clone());
        let mut engine = GentleEngine::from_state(state);
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = tmp.path().with_extension("fa");
        let path_text = path.display().to_string();

        engine
            .apply(Operation::SaveFile {
                seq_id: "synth".to_string(),
                path: path_text.clone(),
                format: ExportFormat::Fasta,
            })
            .unwrap();
        engine
            .apply(Operation::LoadFile {
                path: path_text,
                as_id: Some("roundtrip".to_string()),
            })
            .unwrap();

        let actual = engine
            .state()
            .sequences
            .get("roundtrip")
            .expect("roundtrip sequence should exist");
        assert_eq!(actual.assert_sequence_equality(&expected), Ok(()));
    }

    fn write_gzip(path: &Path, text: &str) {
        let file = std::fs::File::create(path).unwrap();
        let mut encoder = GzEncoder::new(file, Compression::default());
        encoder.write_all(text.as_bytes()).unwrap();
        encoder.finish().unwrap();
    }

    #[cfg(unix)]
    fn install_fake_rnapkin(path: &Path) -> String {
        let script_path = path.join("fake_rnapkin.sh");
        let script = r#"#!/bin/sh
if [ "$1" = "-v" ] && [ "$2" = "-p" ]; then
  seq="$3"
  echo "rnapkin textual report"
  echo "sequence_length=${#seq}"
  echo "points:"
  echo "0.0 0.0"
  echo "1.0 1.0"
  exit 0
fi

if [ "$#" -eq 2 ]; then
  seq="$1"
  out="$2"
  cat > "$out" <<EOF
<svg xmlns="http://www.w3.org/2000/svg" width="160" height="80"><text x="10" y="40">$seq</text></svg>
EOF
  echo "wrote $out" >&2
  exit 0
fi

echo "unexpected args: $@" >&2
exit 2
"#;
        std::fs::write(&script_path, script).expect("write fake rnapkin");
        let mut perms = std::fs::metadata(&script_path)
            .expect("metadata fake rnapkin")
            .permissions();
        perms.set_mode(0o755);
        std::fs::set_permissions(&script_path, perms).expect("chmod fake rnapkin");
        script_path.display().to_string()
    }

    #[cfg(unix)]
    fn install_fake_bigwig_to_bedgraph(path: &Path, bedgraph_source: &Path) -> String {
        let script_path = path.join("fake_bigwig_to_bedgraph.sh");
        let script = format!(
            "#!/bin/sh\nif [ \"$#\" -ne 2 ]; then\n  echo \"expected INPUT.bw OUTPUT.bedGraph\" >&2\n  exit 2\nfi\ncp \"{}\" \"$2\"\n",
            bedgraph_source.display()
        );
        std::fs::write(&script_path, script).expect("write fake bigwig converter");
        let mut perms = std::fs::metadata(&script_path)
            .expect("metadata fake bigwig converter")
            .permissions();
        perms.set_mode(0o755);
        std::fs::set_permissions(&script_path, perms).expect("chmod fake bigwig converter");
        script_path.display().to_string()
    }

    fn file_url(path: &Path) -> String {
        format!("file://{}", path.display())
    }

    struct EnvVarGuard {
        key: &'static str,
        previous: Option<String>,
    }

    impl EnvVarGuard {
        fn set(key: &'static str, value: &str) -> Self {
            let previous = env::var(key).ok();
            env::set_var(key, value);
            Self { key, previous }
        }
    }

    impl Drop for EnvVarGuard {
        fn drop(&mut self) {
            match &self.previous {
                Some(value) => env::set_var(self.key, value),
                None => env::remove_var(self.key),
            }
        }
    }

    #[test]
    fn test_set_display_visibility() {
        let mut engine = GentleEngine::new();
        assert!(!engine.state().display.show_tfbs);
        assert!(engine.state().display.show_cds_features);
        let res = engine
            .apply(Operation::SetDisplayVisibility {
                target: DisplayTarget::Features,
                visible: false,
            })
            .unwrap();
        assert!(res.messages.iter().any(|m| m.contains("features")));
        assert!(!engine.state().display.show_features);
        let res_tfbs = engine
            .apply(Operation::SetDisplayVisibility {
                target: DisplayTarget::Tfbs,
                visible: true,
            })
            .unwrap();
        assert!(res_tfbs.messages.iter().any(|m| m.contains("tfbs")));
        assert!(engine.state().display.show_tfbs);

        let res_cds = engine
            .apply(Operation::SetDisplayVisibility {
                target: DisplayTarget::CdsFeatures,
                visible: false,
            })
            .unwrap();
        assert!(res_cds.messages.iter().any(|m| m.contains("cds_features")));
        assert!(!engine.state().display.show_cds_features);
    }

    #[test]
    fn test_set_linear_viewport() {
        let mut engine = GentleEngine::new();
        let res = engine
            .apply(Operation::SetLinearViewport {
                start_bp: 123,
                span_bp: 456,
            })
            .unwrap();
        assert!(res.messages.iter().any(|m| m.contains("linear viewport")));
        assert_eq!(engine.state().display.linear_view_start_bp, 123);
        assert_eq!(engine.state().display.linear_view_span_bp, 456);
    }

    #[test]
    fn test_extract_region() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("x".to_string(), seq(&"ATGC".repeat(100)));
        let mut engine = GentleEngine::from_state(state);
        let res = engine
            .apply(Operation::ExtractRegion {
                input: "x".to_string(),
                from: 2,
                to: 7,
                output_id: Some("part".to_string()),
            })
            .unwrap();
        assert_eq!(res.created_seq_ids, vec!["part".to_string()]);
        assert_eq!(
            engine
                .state()
                .sequences
                .get("part")
                .unwrap()
                .get_forward_string(),
            "GCATG"
        );
    }

    #[test]
    fn test_extract_anchored_region_with_constraints_unique() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("tpl".to_string(), seq("TTTGAATTCACGTACGTACGTTAAACCC"));
        let mut engine = GentleEngine::from_state(state);
        let res = engine
            .apply(Operation::ExtractAnchoredRegion {
                input: "tpl".to_string(),
                anchor: AnchoredRegionAnchor::Position { zero_based: 21 },
                direction: AnchorDirection::Upstream,
                target_length_bp: 18,
                length_tolerance_bp: 2,
                required_re_sites: vec!["EcoRI".to_string()],
                required_tf_motifs: vec!["CACGTA".to_string()],
                forward_primer: Some("GAATTC".to_string()),
                reverse_primer: Some("ACGT".to_string()),
                output_prefix: Some("prom".to_string()),
                unique: Some(true),
                max_candidates: Some(5),
            })
            .unwrap();
        assert_eq!(res.created_seq_ids, vec!["prom_1".to_string()]);
        assert_eq!(
            engine
                .state()
                .sequences
                .get("prom_1")
                .unwrap()
                .get_forward_string(),
            "GAATTCACGTACGTACGT"
        );
    }

    #[test]
    fn test_extract_anchored_region_from_feature_boundary() {
        let mut engine = GentleEngine::new();
        engine
            .apply(Operation::LoadFile {
                path: "test_files/pGEX-3X.gb".to_string(),
                as_id: Some("pgex".to_string()),
            })
            .unwrap();
        let res = engine
            .apply(Operation::ExtractAnchoredRegion {
                input: "pgex".to_string(),
                anchor: AnchoredRegionAnchor::FeatureBoundary {
                    feature_kind: Some("CDS".to_string()),
                    feature_label: None,
                    boundary: AnchorBoundary::Start,
                    occurrence: Some(0),
                },
                direction: AnchorDirection::Upstream,
                target_length_bp: 60,
                length_tolerance_bp: 0,
                required_re_sites: vec![],
                required_tf_motifs: vec![],
                forward_primer: None,
                reverse_primer: None,
                output_prefix: Some("cds_up".to_string()),
                unique: Some(true),
                max_candidates: Some(1),
            })
            .unwrap();
        assert_eq!(res.created_seq_ids, vec!["cds_up_1".to_string()]);
        assert_eq!(engine.state().sequences.get("cds_up_1").unwrap().len(), 60);
    }

    #[test]
    fn test_ligation_simple_concatenation() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("a".to_string(), seq(&"ATGC".repeat(40)));
        state
            .sequences
            .insert("b".to_string(), seq(&"TTAA".repeat(40)));
        let mut engine = GentleEngine::from_state(state);
        let res = engine
            .apply(Operation::Ligation {
                inputs: vec!["a".to_string(), "b".to_string()],
                circularize_if_possible: false,
                output_id: None,
                protocol: LigationProtocol::Blunt,
                output_prefix: Some("ab".to_string()),
                unique: None,
            })
            .unwrap();
        assert_eq!(res.created_seq_ids.len(), 2);
        let out = engine.state().sequences.get("ab_1").unwrap();
        assert_eq!(
            out.get_forward_string(),
            format!("{}{}", "ATGC".repeat(40), "TTAA".repeat(40))
        );
        assert!(!out.is_circular());
    }

    #[test]
    fn test_merge_containers_creates_pool_copies() {
        let mut state = ProjectState::default();
        state.sequences.insert("a".to_string(), seq("ATGC"));
        state.sequences.insert("b".to_string(), seq("TTAA"));
        let mut engine = GentleEngine::from_state(state);
        let res = engine
            .apply(Operation::MergeContainers {
                inputs: vec!["a".to_string(), "b".to_string()],
                output_prefix: Some("pool".to_string()),
            })
            .unwrap();
        assert_eq!(
            res.created_seq_ids,
            vec!["pool_1".to_string(), "pool_2".to_string()]
        );
        assert_eq!(
            engine
                .state()
                .sequences
                .get("pool_1")
                .unwrap()
                .get_forward_string(),
            "ATGC"
        );
        assert_eq!(
            engine
                .state()
                .sequences
                .get("pool_2")
                .unwrap()
                .get_forward_string(),
            "TTAA"
        );
    }

    #[test]
    fn test_merge_containers_respects_max_fragments() {
        let mut state = ProjectState::default();
        state.parameters.max_fragments_per_container = 1;
        state.sequences.insert("a".to_string(), seq("ATGC"));
        state.sequences.insert("b".to_string(), seq("TTAA"));
        let mut engine = GentleEngine::from_state(state);
        let err = engine
            .apply(Operation::MergeContainers {
                inputs: vec!["a".to_string(), "b".to_string()],
                output_prefix: Some("pool".to_string()),
            })
            .unwrap_err();
        assert!(err.message.contains("max_fragments_per_container"));
    }

    #[test]
    fn test_ligation_protocol_blunt_enumerates_ordered_pairs() {
        let mut state = ProjectState::default();
        state.sequences.insert("a".to_string(), seq("ATGC"));
        state.sequences.insert("b".to_string(), seq("TTAA"));
        let mut engine = GentleEngine::from_state(state);
        let res = engine
            .apply(Operation::Ligation {
                inputs: vec!["a".to_string(), "b".to_string()],
                circularize_if_possible: false,
                output_id: None,
                protocol: LigationProtocol::Blunt,
                output_prefix: Some("lig".to_string()),
                unique: Some(false),
            })
            .unwrap();
        assert_eq!(res.created_seq_ids.len(), 2);
        assert_eq!(
            engine
                .state()
                .sequences
                .get("lig_1")
                .unwrap()
                .get_forward_string(),
            "ATGCTTAA"
        );
        assert_eq!(
            engine
                .state()
                .sequences
                .get("lig_2")
                .unwrap()
                .get_forward_string(),
            "TTAAATGC"
        );
    }

    #[test]
    fn test_ligation_protocol_sticky_uses_overhang_compatibility() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("x".to_string(), seq("ATGGATCCGCATGGATCCGCATGGATCCGC"));
        let mut engine = GentleEngine::from_state(state);
        let digest_res = engine
            .apply(Operation::Digest {
                input: "x".to_string(),
                enzymes: vec!["BamHI".to_string()],
                output_prefix: Some("frag".to_string()),
            })
            .unwrap();
        assert!(digest_res.created_seq_ids.len() >= 2);
        let a = digest_res.created_seq_ids[0].clone();
        let b = digest_res.created_seq_ids[1].clone();

        let lig_res = engine
            .apply(Operation::Ligation {
                inputs: vec![a, b],
                circularize_if_possible: false,
                output_id: None,
                protocol: LigationProtocol::Sticky,
                output_prefix: Some("st".to_string()),
                unique: Some(false),
            })
            .unwrap();
        assert!(!lig_res.created_seq_ids.is_empty());
    }

    #[test]
    fn test_workflow_digest_merge_ligation_is_deterministic() {
        let mut base = ProjectState::default();
        base.sequences
            .insert("x".to_string(), seq("ATGGATCCGCATGGATCCGCATGGATCCGC"));

        let run_once = |state: ProjectState| {
            let mut engine = GentleEngine::from_state(state);
            let digest = engine
                .apply(Operation::Digest {
                    input: "x".to_string(),
                    enzymes: vec!["BamHI".to_string(), "EcoRI".to_string()],
                    output_prefix: Some("d".to_string()),
                })
                .unwrap();
            let merge = engine
                .apply(Operation::MergeContainers {
                    inputs: digest.created_seq_ids.clone(),
                    output_prefix: Some("m".to_string()),
                })
                .unwrap();
            let lig = engine
                .apply(Operation::Ligation {
                    inputs: merge.created_seq_ids.clone(),
                    circularize_if_possible: false,
                    output_id: None,
                    protocol: LigationProtocol::Sticky,
                    output_prefix: Some("lig".to_string()),
                    unique: Some(false),
                })
                .unwrap();
            lig.created_seq_ids
        };

        let a = run_once(base.clone());
        let b = run_once(base.clone());
        assert_eq!(a, b);
        assert!(!a.is_empty());
        assert_eq!(a.first().unwrap(), "lig_1");
    }

    #[test]
    fn test_lineage_extract_creates_parent_child_edge() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("x".to_string(), seq(&"ATGC".repeat(100)));
        let mut engine = GentleEngine::from_state(state);
        let _ = engine
            .apply(Operation::ExtractRegion {
                input: "x".to_string(),
                from: 2,
                to: 7,
                output_id: Some("part".to_string()),
            })
            .unwrap();

        let lineage = &engine.state().lineage;
        let x_node = lineage.seq_to_node.get("x").unwrap();
        let part_node = lineage.seq_to_node.get("part").unwrap();
        assert!(lineage
            .edges
            .iter()
            .any(|e| e.from_node_id == *x_node && e.to_node_id == *part_node));
    }

    #[test]
    fn test_lineage_ligation_has_two_parents() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("a".to_string(), seq(&"ATGC".repeat(40)));
        state
            .sequences
            .insert("b".to_string(), seq(&"TTAA".repeat(40)));
        let mut engine = GentleEngine::from_state(state);
        let res = engine
            .apply(Operation::Ligation {
                inputs: vec!["a".to_string(), "b".to_string()],
                circularize_if_possible: false,
                output_id: None,
                protocol: LigationProtocol::Blunt,
                output_prefix: Some("ab".to_string()),
                unique: None,
            })
            .unwrap();

        let lineage = &engine.state().lineage;
        let a_node = lineage.seq_to_node.get("a").unwrap();
        let b_node = lineage.seq_to_node.get("b").unwrap();
        let ab_node = lineage.seq_to_node.get(&res.created_seq_ids[0]).unwrap();
        assert!(lineage
            .edges
            .iter()
            .any(|e| e.from_node_id == *a_node && e.to_node_id == *ab_node));
        assert!(lineage
            .edges
            .iter()
            .any(|e| e.from_node_id == *b_node && e.to_node_id == *ab_node));
    }

    #[test]
    fn test_select_candidate_creates_in_silico_node() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("frag".to_string(), seq(&"ATGC".repeat(40)));
        let mut engine = GentleEngine::from_state(state);
        let res = engine
            .apply(Operation::SelectCandidate {
                input: "frag".to_string(),
                criterion: "band_size_range:150-170bp".to_string(),
                output_id: Some("picked".to_string()),
            })
            .unwrap();

        assert_eq!(res.created_seq_ids, vec!["picked".to_string()]);
        assert!(!res.warnings.is_empty());
        let lineage = &engine.state().lineage;
        let picked_node = lineage.seq_to_node.get("picked").unwrap();
        let picked_origin = &lineage.nodes.get(picked_node).unwrap().origin;
        assert!(matches!(picked_origin, SequenceOrigin::InSilicoSelection));
    }

    #[test]
    fn test_reverse_complement_reverse_complement_and_branch() {
        let mut state = ProjectState::default();
        state.sequences.insert("s".to_string(), seq("ATGCCA"));
        let mut engine = GentleEngine::from_state(state);

        let res_rev = engine
            .apply(Operation::Reverse {
                input: "s".to_string(),
                output_id: Some("s_rev".to_string()),
            })
            .unwrap();
        assert_eq!(res_rev.created_seq_ids, vec!["s_rev".to_string()]);
        assert_eq!(
            engine
                .state()
                .sequences
                .get("s_rev")
                .unwrap()
                .get_forward_string(),
            "ACCGTA"
        );

        let res_comp = engine
            .apply(Operation::Complement {
                input: "s".to_string(),
                output_id: Some("s_comp".to_string()),
            })
            .unwrap();
        assert_eq!(res_comp.created_seq_ids, vec!["s_comp".to_string()]);
        assert_eq!(
            engine
                .state()
                .sequences
                .get("s_comp")
                .unwrap()
                .get_forward_string(),
            "TACGGT"
        );

        let res_rc = engine
            .apply(Operation::ReverseComplement {
                input: "s".to_string(),
                output_id: Some("s_rc".to_string()),
            })
            .unwrap();
        assert_eq!(res_rc.created_seq_ids, vec!["s_rc".to_string()]);
        assert_eq!(
            engine
                .state()
                .sequences
                .get("s_rc")
                .unwrap()
                .get_forward_string(),
            "TGGCAT"
        );

        let res_split = engine
            .apply(Operation::Branch {
                input: "s".to_string(),
                output_id: Some("s_branch".to_string()),
            })
            .unwrap();
        assert_eq!(res_split.created_seq_ids, vec!["s_branch".to_string()]);
        assert_eq!(
            engine
                .state()
                .sequences
                .get("s_branch")
                .unwrap()
                .get_forward_string(),
            "ATGCCA"
        );

        let lineage = &engine.state().lineage;
        let s_node = lineage.seq_to_node.get("s").unwrap();
        for derived in ["s_rev", "s_comp", "s_rc", "s_branch"] {
            let dnode = lineage.seq_to_node.get(derived).unwrap();
            assert!(lineage
                .edges
                .iter()
                .any(|e| e.from_node_id == *s_node && e.to_node_id == *dnode));
        }
    }

    #[test]
    fn test_default_max_fragments_per_container_is_80000() {
        let state = ProjectState::default();
        assert_eq!(state.parameters.max_fragments_per_container, 80_000);
    }

    #[test]
    fn test_set_parameter_max_fragments_per_container() {
        let mut engine = GentleEngine::new();
        let res = engine
            .apply(Operation::SetParameter {
                name: "max_fragments_per_container".to_string(),
                value: serde_json::json!(1234),
            })
            .unwrap();
        assert!(res
            .messages
            .iter()
            .any(|m| m.contains("max_fragments_per_container")));
        assert_eq!(engine.state().parameters.max_fragments_per_container, 1234);
    }

    #[test]
    fn test_set_parameter_feature_details_font_size() {
        let mut engine = GentleEngine::new();
        let res = engine
            .apply(Operation::SetParameter {
                name: "feature_details_font_size".to_string(),
                value: serde_json::json!(9.5),
            })
            .unwrap();
        assert!(res
            .messages
            .iter()
            .any(|m| m.contains("feature_details_font_size")));
        assert!((engine.state().display.feature_details_font_size - 9.5).abs() < f32::EPSILON);
    }

    #[test]
    fn test_digest_respects_max_fragments_per_container() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("x".to_string(), seq("ATGGATCCGCATGGATCCGCATGGATCCGC"));
        state.parameters.max_fragments_per_container = 2;
        let mut engine = GentleEngine::from_state(state);
        let err = engine
            .apply(Operation::Digest {
                input: "x".to_string(),
                enzymes: vec!["BamHI".to_string()],
                output_prefix: Some("frag".to_string()),
            })
            .unwrap_err();
        assert!(err.message.contains("max_fragments_per_container"));
    }

    #[test]
    fn test_filter_by_molecular_weight_with_error_range() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("a".to_string(), seq(&"A".repeat(100)));
        state
            .sequences
            .insert("b".to_string(), seq(&"A".repeat(180)));
        state
            .sequences
            .insert("c".to_string(), seq(&"A".repeat(260)));
        let mut engine = GentleEngine::from_state(state);

        let res = engine
            .apply(Operation::FilterByMolecularWeight {
                inputs: vec!["a".to_string(), "b".to_string(), "c".to_string()],
                min_bp: 150,
                max_bp: 200,
                error: 0.10,
                unique: false,
                output_prefix: Some("mw".to_string()),
            })
            .unwrap();

        assert_eq!(res.created_seq_ids.len(), 1);
        let out = engine.state().sequences.get("mw_1").unwrap();
        assert_eq!(out.len(), 180);
    }

    #[test]
    fn test_filter_by_molecular_weight_unique_fails_on_multiple_matches() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("a".to_string(), seq(&"A".repeat(100)));
        state
            .sequences
            .insert("b".to_string(), seq(&"A".repeat(105)));
        state
            .sequences
            .insert("c".to_string(), seq(&"A".repeat(200)));
        let mut engine = GentleEngine::from_state(state);

        let err = engine
            .apply(Operation::FilterByMolecularWeight {
                inputs: vec!["a".to_string(), "b".to_string(), "c".to_string()],
                min_bp: 95,
                max_bp: 105,
                error: 0.10,
                unique: true,
                output_prefix: Some("mw".to_string()),
            })
            .unwrap_err();

        assert!(err.message.contains("exactly one match"));
    }

    #[test]
    fn test_filter_by_design_constraints_practical_filters() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("good".to_string(), seq("GACTGACTGACTGACTGACT"));
        state
            .sequences
            .insert("low_gc".to_string(), seq("ATATATATATATATATATAT"));
        state
            .sequences
            .insert("homopoly".to_string(), seq("GACAAAAAGACTGACTGACT"));
        state
            .sequences
            .insert("u6_t4".to_string(), seq("GACTTTTGACTGACTGACT"));
        state
            .sequences
            .insert("amb".to_string(), seq("GACTNNACTGACTGACTGAC"));
        let mut engine = GentleEngine::from_state(state);

        let res = engine
            .apply(Operation::FilterByDesignConstraints {
                inputs: vec![
                    "good".to_string(),
                    "low_gc".to_string(),
                    "homopoly".to_string(),
                    "u6_t4".to_string(),
                    "amb".to_string(),
                ],
                gc_min: Some(0.30),
                gc_max: Some(0.70),
                max_homopolymer_run: Some(4),
                reject_ambiguous_bases: Some(true),
                avoid_u6_terminator_tttt: Some(true),
                forbidden_motifs: vec![],
                unique: false,
                output_prefix: Some("design".to_string()),
            })
            .unwrap();

        assert_eq!(res.created_seq_ids.len(), 1);
        assert_eq!(res.created_seq_ids[0], "design_1".to_string());
        assert!(engine.state().sequences.contains_key("design_1"));
    }

    #[test]
    fn test_filter_by_design_constraints_unique_fails_on_multiple_matches() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("a".to_string(), seq("GACTGACTGACTGACTGACT"));
        state
            .sequences
            .insert("b".to_string(), seq("GACCGACTGACTGACTGACC"));
        let mut engine = GentleEngine::from_state(state);

        let err = engine
            .apply(Operation::FilterByDesignConstraints {
                inputs: vec!["a".to_string(), "b".to_string()],
                gc_min: Some(0.20),
                gc_max: Some(0.80),
                max_homopolymer_run: Some(6),
                reject_ambiguous_bases: Some(true),
                avoid_u6_terminator_tttt: Some(true),
                forbidden_motifs: vec![],
                unique: true,
                output_prefix: Some("design".to_string()),
            })
            .unwrap_err();

        assert!(err.message.contains("exactly one match"));
    }

    #[test]
    fn test_filter_by_design_constraints_accepts_legacy_operation_name() {
        let json = r#"{
            "FilterBySequenceQuality": {
                "inputs": ["a"],
                "gc_min": 0.30,
                "gc_max": 0.70,
                "max_homopolymer_run": 4,
                "reject_ambiguous_bases": true,
                "avoid_u6_terminator_tttt": true,
                "forbidden_motifs": [],
                "unique": false,
                "output_prefix": "legacy"
            }
        }"#;
        let op: Operation = serde_json::from_str(json).expect("legacy op json parses");
        match op {
            Operation::FilterByDesignConstraints {
                inputs,
                gc_min,
                gc_max,
                max_homopolymer_run,
                unique,
                output_prefix,
                ..
            } => {
                assert_eq!(inputs, vec!["a".to_string()]);
                assert_eq!(gc_min, Some(0.30));
                assert_eq!(gc_max, Some(0.70));
                assert_eq!(max_homopolymer_run, Some(4));
                assert!(!unique);
                assert_eq!(output_prefix, Some("legacy".to_string()));
            }
            other => panic!("unexpected operation variant: {:?}", other),
        }
    }

    #[test]
    fn test_pcr_single_amplicon() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("tpl".to_string(), seq("ATGAAACCCGGGTTT"));
        let mut engine = GentleEngine::from_state(state);
        let res = engine
            .apply(Operation::Pcr {
                template: "tpl".to_string(),
                forward_primer: "ATGAAA".to_string(),
                reverse_primer: "AAACCC".to_string(),
                output_id: Some("amp".to_string()),
                unique: Some(true),
            })
            .unwrap();
        assert_eq!(res.created_seq_ids, vec!["amp".to_string()]);
        assert_eq!(
            engine
                .state()
                .sequences
                .get("amp")
                .unwrap()
                .get_forward_string(),
            "ATGAAACCCGGGTTT"
        );
    }

    #[test]
    fn test_pcr_unique_fails_on_multiple_amplicons() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("tpl".to_string(), seq("AAAACCCCGGGGAAAACCCCGGGG"));
        let mut engine = GentleEngine::from_state(state);
        let err = engine
            .apply(Operation::Pcr {
                template: "tpl".to_string(),
                forward_primer: "AAAA".to_string(),
                reverse_primer: "CCCC".to_string(),
                output_id: None,
                unique: Some(true),
            })
            .unwrap_err();
        assert!(err.message.contains("unique=true"));
    }

    #[test]
    fn test_pcr_multiple_amplicons_without_unique() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("tpl".to_string(), seq("AAAACCCCGGGGAAAACCCCGGGG"));
        let mut engine = GentleEngine::from_state(state);
        let res = engine
            .apply(Operation::Pcr {
                template: "tpl".to_string(),
                forward_primer: "AAAA".to_string(),
                reverse_primer: "CCCC".to_string(),
                output_id: None,
                unique: Some(false),
            })
            .unwrap();
        assert_eq!(res.created_seq_ids.len(), 3);
    }

    #[test]
    fn test_pcr_advanced_inserts_5prime_tail() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("tpl".to_string(), seq("ATGAAACCCGGGTTT"));
        let mut engine = GentleEngine::from_state(state);
        let res = engine
            .apply(Operation::PcrAdvanced {
                template: "tpl".to_string(),
                forward_primer: PcrPrimerSpec {
                    sequence: "GGATCCATGAAA".to_string(),
                    anneal_len: Some(6),
                    max_mismatches: Some(0),
                    require_3prime_exact_bases: Some(4),
                    library_mode: None,
                    max_variants: None,
                    sample_seed: None,
                },
                reverse_primer: PcrPrimerSpec {
                    sequence: "AAACCC".to_string(),
                    anneal_len: Some(6),
                    max_mismatches: Some(0),
                    require_3prime_exact_bases: Some(4),
                    library_mode: None,
                    max_variants: None,
                    sample_seed: None,
                },
                output_id: Some("amp_adv".to_string()),
                unique: Some(true),
            })
            .unwrap();
        assert_eq!(res.created_seq_ids, vec!["amp_adv".to_string()]);
        assert_eq!(
            engine
                .state()
                .sequences
                .get("amp_adv")
                .unwrap()
                .get_forward_string(),
            "GGATCCATGAAACCCGGGTTT"
        );
    }

    #[test]
    fn test_pcr_advanced_allows_partial_match_and_introduces_mutation() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("tpl".to_string(), seq("ATGAAACCCGGGTTT"));
        let mut engine = GentleEngine::from_state(state);
        let res = engine
            .apply(Operation::PcrAdvanced {
                template: "tpl".to_string(),
                forward_primer: PcrPrimerSpec {
                    sequence: "ATCAAA".to_string(),
                    anneal_len: Some(6),
                    max_mismatches: Some(1),
                    require_3prime_exact_bases: Some(3),
                    library_mode: None,
                    max_variants: None,
                    sample_seed: None,
                },
                reverse_primer: PcrPrimerSpec {
                    sequence: "AAACCC".to_string(),
                    anneal_len: Some(6),
                    max_mismatches: Some(0),
                    require_3prime_exact_bases: Some(4),
                    library_mode: None,
                    max_variants: None,
                    sample_seed: None,
                },
                output_id: Some("amp_mut".to_string()),
                unique: Some(true),
            })
            .unwrap();
        assert_eq!(res.created_seq_ids, vec!["amp_mut".to_string()]);
        assert_eq!(
            engine
                .state()
                .sequences
                .get("amp_mut")
                .unwrap()
                .get_forward_string(),
            "ATCAAACCCGGGTTT"
        );
    }

    #[test]
    fn test_pcr_mutagenesis_single_snp_success() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("tpl".to_string(), seq("ATGAAACCCGGGTTT"));
        let mut engine = GentleEngine::from_state(state);
        let res = engine
            .apply(Operation::PcrMutagenesis {
                template: "tpl".to_string(),
                forward_primer: PcrPrimerSpec {
                    sequence: "ATCAAA".to_string(),
                    anneal_len: Some(6),
                    max_mismatches: Some(1),
                    require_3prime_exact_bases: Some(3),
                    library_mode: None,
                    max_variants: None,
                    sample_seed: None,
                },
                reverse_primer: PcrPrimerSpec {
                    sequence: "AAACCC".to_string(),
                    anneal_len: Some(6),
                    max_mismatches: Some(0),
                    require_3prime_exact_bases: Some(4),
                    library_mode: None,
                    max_variants: None,
                    sample_seed: None,
                },
                mutations: vec![SnpMutationSpec {
                    zero_based_position: 2,
                    reference: "G".to_string(),
                    alternate: "C".to_string(),
                }],
                output_id: Some("mut1".to_string()),
                unique: Some(true),
                require_all_mutations: Some(true),
            })
            .unwrap();
        assert_eq!(res.created_seq_ids, vec!["mut1".to_string()]);
        assert_eq!(
            engine
                .state()
                .sequences
                .get("mut1")
                .unwrap()
                .get_forward_string(),
            "ATCAAACCCGGGTTT"
        );
    }

    #[test]
    fn test_pcr_mutagenesis_fails_when_requested_snp_not_introduced() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("tpl".to_string(), seq("ATGAAACCCGGGTTT"));
        let mut engine = GentleEngine::from_state(state);
        let err = engine
            .apply(Operation::PcrMutagenesis {
                template: "tpl".to_string(),
                forward_primer: PcrPrimerSpec {
                    sequence: "ATCAAA".to_string(),
                    anneal_len: Some(6),
                    max_mismatches: Some(1),
                    require_3prime_exact_bases: Some(3),
                    library_mode: None,
                    max_variants: None,
                    sample_seed: None,
                },
                reverse_primer: PcrPrimerSpec {
                    sequence: "AAACCC".to_string(),
                    anneal_len: Some(6),
                    max_mismatches: Some(0),
                    require_3prime_exact_bases: Some(4),
                    library_mode: None,
                    max_variants: None,
                    sample_seed: None,
                },
                mutations: vec![SnpMutationSpec {
                    zero_based_position: 2,
                    reference: "G".to_string(),
                    alternate: "T".to_string(),
                }],
                output_id: Some("mut_fail".to_string()),
                unique: Some(true),
                require_all_mutations: Some(true),
            })
            .unwrap_err();
        assert!(err
            .message
            .contains("No amplicon introduced all requested mutations"));
    }

    #[test]
    fn test_pcr_advanced_degenerate_primer_enumerate() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("tpl".to_string(), seq("ATGAAACCCGGGTTT"));
        let mut engine = GentleEngine::from_state(state);
        let res = engine
            .apply(Operation::PcrAdvanced {
                template: "tpl".to_string(),
                forward_primer: PcrPrimerSpec {
                    sequence: "ATNAAA".to_string(),
                    anneal_len: Some(6),
                    max_mismatches: Some(1),
                    require_3prime_exact_bases: Some(3),
                    library_mode: Some(PrimerLibraryMode::Enumerate),
                    max_variants: Some(4),
                    sample_seed: None,
                },
                reverse_primer: PcrPrimerSpec {
                    sequence: "AAACCC".to_string(),
                    anneal_len: Some(6),
                    max_mismatches: Some(0),
                    require_3prime_exact_bases: Some(4),
                    library_mode: None,
                    max_variants: None,
                    sample_seed: None,
                },
                output_id: None,
                unique: Some(false),
            })
            .unwrap();
        assert_eq!(res.created_seq_ids.len(), 4);
    }

    #[test]
    fn test_pcr_advanced_degenerate_primer_sample_mode() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("tpl".to_string(), seq("ATGAAACCCGGGTTT"));
        let mut engine = GentleEngine::from_state(state);
        let res = engine
            .apply(Operation::PcrAdvanced {
                template: "tpl".to_string(),
                forward_primer: PcrPrimerSpec {
                    sequence: "ATNAAA".to_string(),
                    anneal_len: Some(6),
                    max_mismatches: Some(1),
                    require_3prime_exact_bases: Some(3),
                    library_mode: Some(PrimerLibraryMode::Sample),
                    max_variants: Some(2),
                    sample_seed: Some(42),
                },
                reverse_primer: PcrPrimerSpec {
                    sequence: "AAACCC".to_string(),
                    anneal_len: Some(6),
                    max_mismatches: Some(0),
                    require_3prime_exact_bases: Some(4),
                    library_mode: None,
                    max_variants: None,
                    sample_seed: None,
                },
                output_id: None,
                unique: Some(false),
            })
            .unwrap();
        assert_eq!(res.created_seq_ids.len(), 2);
    }

    #[test]
    fn test_load_file_operation() {
        let mut engine = GentleEngine::new();
        let res = engine
            .apply(Operation::LoadFile {
                path: "test_files/pGEX_3X.fa".to_string(),
                as_id: Some("pgex".to_string()),
            })
            .unwrap();
        assert_eq!(res.created_seq_ids, vec!["pgex".to_string()]);
        assert!(engine.state().sequences.contains_key("pgex"));
    }

    #[test]
    fn test_load_file_operation_genbank_region_anchor_enables_bed_import() {
        let mut engine = GentleEngine::new();
        let res = engine
            .apply(Operation::LoadFile {
                path: "test_files/tp73.ncbi.gb".to_string(),
                as_id: Some("tp73".to_string()),
            })
            .unwrap();
        assert_eq!(res.created_seq_ids, vec!["tp73".to_string()]);
        assert!(engine
            .list_sequences_with_genome_anchor()
            .iter()
            .any(|seq_id| seq_id == "tp73"));

        let td = tempdir().unwrap();
        let bed_path = td.path().join("tp73_anchor_test.bed");
        std::fs::write(&bed_path, "chr1\t3652515\t3652525\tpeak1\t100\t+\n").unwrap();

        let import_res = engine
            .apply(Operation::ImportGenomeBedTrack {
                seq_id: "tp73".to_string(),
                path: bed_path.display().to_string(),
                track_name: Some("anchor-test".to_string()),
                min_score: None,
                max_score: None,
                clear_existing: Some(true),
            })
            .unwrap();
        assert!(import_res
            .changed_seq_ids
            .iter()
            .any(|seq_id| seq_id == "tp73"));

        let tp73 = engine
            .state()
            .sequences
            .get("tp73")
            .expect("tp73 should exist");
        assert!(tp73
            .features()
            .iter()
            .any(GentleEngine::is_generated_genome_bed_feature));
    }

    #[test]
    fn test_parse_genbank_accession_region_supports_complement() {
        let td = tempdir().unwrap();
        let gb_path = td.path().join("complement_header.gb");
        let gb_text = "\
LOCUS       TEST000001              11 bp    DNA     linear   CON 01-JAN-2000
DEFINITION  Homo sapiens chromosome 1, GRCh38.p14 Primary Assembly.
ACCESSION   NC_000001 REGION: complement(3652516..3652526)
ORIGIN
        1 atgcgatgcga
//
";
        std::fs::write(&gb_path, gb_text).unwrap();
        let parsed =
            GentleEngine::parse_genbank_accession_region(gb_path.to_string_lossy().as_ref())
                .expect("region should parse");
        assert_eq!(parsed.0, "NC_000001");
        assert_eq!(parsed.1, 3652516);
        assert_eq!(parsed.2, 3652526);
        assert_eq!(parsed.4, '-');
    }

    #[test]
    fn test_import_genome_bed_track_remaps_for_reverse_anchor() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("rev".to_string(), seq("ACGTACGTACG"));
        state.metadata.insert(
            PROVENANCE_METADATA_KEY.to_string(),
            serde_json::json!({
                GENOME_EXTRACTIONS_METADATA_KEY: [
                    {
                        "seq_id": "rev",
                        "recorded_at_unix_ms": 1,
                        "operation": "LoadFileGenBankRegion",
                        "genome_id": "GRCh38.p14",
                        "catalog_path": "synthetic",
                        "cache_dir": null,
                        "chromosome": "1",
                        "start_1based": 100,
                        "end_1based": 110,
                        "gene_query": null,
                        "occurrence": null,
                        "gene_id": null,
                        "gene_name": null,
                        "strand": null,
                        "anchor_strand": "-",
                        "sequence_source_type": "genbank_file",
                        "annotation_source_type": "genbank_file",
                        "sequence_source": "synthetic",
                        "annotation_source": "synthetic",
                        "sequence_sha1": null,
                        "annotation_sha1": null
                    }
                ]
            }),
        );
        let mut engine = GentleEngine::from_state(state);

        let td = tempdir().unwrap();
        let bed_path = td.path().join("reverse_anchor.bed");
        std::fs::write(&bed_path, "chr1\t99\t101\tpeak1\t100\t+\n").unwrap();

        let result = engine
            .apply(Operation::ImportGenomeBedTrack {
                seq_id: "rev".to_string(),
                path: bed_path.display().to_string(),
                track_name: Some("rev-track".to_string()),
                min_score: None,
                max_score: None,
                clear_existing: Some(true),
            })
            .unwrap();
        assert!(result.changed_seq_ids.iter().any(|id| id == "rev"));

        let dna = engine
            .state()
            .sequences
            .get("rev")
            .expect("rev sequence should exist");
        let feature = dna
            .features()
            .iter()
            .find(|f| GentleEngine::is_generated_genome_bed_feature(f))
            .expect("generated BED feature should exist");
        assert!(crate::feature_location::feature_is_reverse(feature));
        let ranges = crate::feature_location::feature_ranges_sorted_i64(feature);
        assert_eq!(ranges, vec![(9, 11)]);
        assert_eq!(
            feature.qualifier_values("bed_strand".into()).next(),
            Some("+")
        );
        assert_eq!(feature.qualifier_values("strand".into()).next(), Some("-"));
    }

    #[test]
    fn test_load_file_operation_fasta_synthetic_oligo_metadata() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("oligo.fa");
        std::fs::write(&path, ">oligo1 molecule=ssdna\nATGCATGC\n").unwrap();

        let mut engine = GentleEngine::new();
        let res = engine
            .apply(Operation::LoadFile {
                path: path.display().to_string(),
                as_id: Some("oligo".to_string()),
            })
            .unwrap();
        assert_eq!(res.created_seq_ids, vec!["oligo".to_string()]);

        let dna = engine
            .state()
            .sequences
            .get("oligo")
            .expect("oligo sequence should exist");
        assert_eq!(dna.molecule_type(), Some("ssDNA"));
        assert!(dna.overhang().is_blunt());
    }

    #[test]
    fn test_load_file_operation_fasta_dsdna_with_overhang_metadata() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("sticky.fa");
        std::fs::write(&path, ">sticky molecule=dsdna f5=GATC r5=CTAG\nATGCATGC\n").unwrap();

        let mut engine = GentleEngine::new();
        engine
            .apply(Operation::LoadFile {
                path: path.display().to_string(),
                as_id: Some("sticky".to_string()),
            })
            .unwrap();

        let dna = engine
            .state()
            .sequences
            .get("sticky")
            .expect("sticky sequence should exist");
        assert_eq!(dna.molecule_type(), Some("dsDNA"));
        assert_eq!(dna.overhang().forward_5, b"GATC".to_vec());
        assert_eq!(dna.overhang().reverse_5, b"CTAG".to_vec());
    }

    #[test]
    fn test_fasta_roundtrip_synthetic_dsdna_blunt() {
        let expected = synth_oligo("molecule=dsdna topology=linear", b"ATGCATGC");
        assert_fasta_roundtrip(expected);
    }

    #[test]
    fn test_fasta_roundtrip_synthetic_ssdna() {
        let expected = synth_oligo("molecule=ssdna topology=linear", b"ATGCATGC");
        assert_fasta_roundtrip(expected);
    }

    #[test]
    fn test_fasta_roundtrip_synthetic_rna() {
        let expected = synth_oligo("molecule=rna topology=linear", b"AUGTT");
        assert_fasta_roundtrip(expected);
    }

    #[test]
    fn test_fasta_roundtrip_synthetic_dsdna_with_overhangs() {
        let expected = synth_oligo(
            "molecule=dsdna f5=GATC r5=CTAG topology=linear",
            b"ATGCATGC",
        );
        assert_fasta_roundtrip(expected);
    }

    #[cfg(unix)]
    #[test]
    fn test_inspect_rna_structure_returns_rnapkin_textual_output() {
        let td = tempdir().unwrap();
        let fake_rnapkin = install_fake_rnapkin(td.path());
        let _bin_guard = EnvVarGuard::set("GENTLE_RNAPKIN_BIN", &fake_rnapkin);

        let mut state = ProjectState::default();
        state
            .sequences
            .insert("rna".to_string(), synth_oligo("molecule=rna", b"AUGCAU"));
        let engine = GentleEngine::from_state(state);
        let report = engine.inspect_rna_structure("rna").unwrap();

        assert_eq!(report.tool, "rnapkin");
        assert!(report.stdout.contains("rnapkin textual report"));
        assert!(report.stdout.contains("points:"));
    }

    #[cfg(unix)]
    #[test]
    fn test_render_rna_structure_svg_operation() {
        let td = tempdir().unwrap();
        let fake_rnapkin = install_fake_rnapkin(td.path());
        let _bin_guard = EnvVarGuard::set("GENTLE_RNAPKIN_BIN", &fake_rnapkin);

        let mut state = ProjectState::default();
        state
            .sequences
            .insert("rna".to_string(), synth_oligo("molecule=rna", b"AUGCAU"));
        let mut engine = GentleEngine::from_state(state);
        let output = td.path().join("rna.structure.svg");
        let output_text = output.display().to_string();

        let res = engine
            .apply(Operation::RenderRnaStructureSvg {
                seq_id: "rna".to_string(),
                path: output_text.clone(),
            })
            .unwrap();

        assert!(res.messages.iter().any(|m| m.contains("RNA structure SVG")));
        let svg = std::fs::read_to_string(output_text).unwrap();
        assert!(svg.contains("<svg"));
    }

    #[test]
    fn test_render_rna_structure_svg_requires_rna_biotype() {
        let mut state = ProjectState::default();
        state.sequences.insert("dna".to_string(), seq("ATGCATGC"));
        let mut engine = GentleEngine::from_state(state);
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let output = tmp.path().with_extension("svg");
        let err = engine
            .apply(Operation::RenderRnaStructureSvg {
                seq_id: "dna".to_string(),
                path: output.display().to_string(),
            })
            .unwrap_err();
        assert!(matches!(err.code, ErrorCode::InvalidInput));
        assert!(err.message.to_ascii_lowercase().contains("rna"));
    }

    #[test]
    fn test_save_file_operation_genbank() {
        let mut state = ProjectState::default();
        state.sequences.insert("s".to_string(), seq("ATGCCA"));
        let mut engine = GentleEngine::from_state(state);
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = tmp.path().with_extension("gb");
        let path_text = path.display().to_string();
        let res = engine
            .apply(Operation::SaveFile {
                seq_id: "s".to_string(),
                path: path_text.clone(),
                format: ExportFormat::GenBank,
            })
            .unwrap();
        assert!(res.changed_seq_ids.contains(&"s".to_string()));
        let text = std::fs::read_to_string(path_text).unwrap();
        assert!(text.contains("LOCUS"));
    }

    #[test]
    fn test_set_topology_operation() {
        let mut state = ProjectState::default();
        state.sequences.insert("s".to_string(), seq("ATGCCA"));
        let mut engine = GentleEngine::from_state(state);
        let res = engine
            .apply(Operation::SetTopology {
                seq_id: "s".to_string(),
                circular: true,
            })
            .unwrap();
        assert!(res.changed_seq_ids.contains(&"s".to_string()));
        assert!(engine.state().sequences.get("s").unwrap().is_circular());
    }

    #[test]
    fn test_recompute_features_operation() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("s".to_string(), seq("ATGAAACCCGGGTTT"));
        let mut engine = GentleEngine::from_state(state);
        let res = engine
            .apply(Operation::RecomputeFeatures {
                seq_id: "s".to_string(),
            })
            .unwrap();
        assert!(res.changed_seq_ids.contains(&"s".to_string()));
    }

    #[test]
    fn test_prepare_scoring_matrices_avoid_negative_infinity() {
        let matrix = vec![[10.0, 0.0, 0.0, 0.0], [5.0, 0.0, 0.0, 0.0]];
        let (llr, true_log_odds) = GentleEngine::prepare_scoring_matrices(&matrix);
        assert_eq!(llr.len(), 2);
        for col in llr {
            for v in col {
                assert!(v.is_finite());
            }
        }
        assert_eq!(true_log_odds.len(), 2);
        for col in true_log_odds {
            for v in col {
                assert!(v.is_finite());
            }
        }
    }

    #[test]
    fn test_annotate_tfbs_adds_scored_features() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("s".to_string(), seq("TTTACGTAAACGTGGG"));
        let mut engine = GentleEngine::from_state(state);
        let res = engine
            .apply(Operation::AnnotateTfbs {
                seq_id: "s".to_string(),
                motifs: vec!["ACGT".to_string()],
                min_llr_bits: Some(0.0),
                min_llr_quantile: Some(0.0),
                per_tf_thresholds: vec![],
                clear_existing: Some(true),
                max_hits: None,
            })
            .unwrap();
        assert!(res.changed_seq_ids.contains(&"s".to_string()));
        let dna = engine.state().sequences.get("s").unwrap();
        let tfbs_features: Vec<_> = dna
            .features()
            .iter()
            .filter(|f| {
                f.qualifier_values("gentle_generated".into())
                    .any(|v| v.eq_ignore_ascii_case("tfbs"))
            })
            .collect();
        assert!(!tfbs_features.is_empty());
        assert!(tfbs_features.iter().all(|f| {
            let motif_len = f
                .qualifier_values("motif_length_bp".into())
                .next()
                .and_then(|v| v.parse::<usize>().ok())
                .unwrap_or(0);
            let span_matches = f
                .location
                .find_bounds()
                .ok()
                .map(|(from, to)| (to - from) as usize == motif_len)
                .unwrap_or(false);
            span_matches
                && motif_len == 4
                && f.qualifier_values("llr_bits".into()).next().is_some()
                && f.qualifier_values("llr_quantile".into()).next().is_some()
                && f.qualifier_values("true_log_odds_bits".into())
                    .next()
                    .is_some()
                && f.qualifier_values("true_log_odds_quantile".into())
                    .next()
                    .is_some()
                && f.qualifier_values("quantile_scope".into()).next().is_some()
        }));
    }

    #[test]
    fn test_annotate_tfbs_progress_reaches_completion() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("s".to_string(), seq("TTTACGTAAACGTGGG"));
        let mut engine = GentleEngine::from_state(state);
        let mut progress_events: Vec<TfbsProgress> = vec![];
        let res = engine
            .apply_with_progress(
                Operation::AnnotateTfbs {
                    seq_id: "s".to_string(),
                    motifs: vec!["ACGT".to_string(), "TATAAA".to_string()],
                    min_llr_bits: Some(0.0),
                    min_llr_quantile: Some(0.0),
                    per_tf_thresholds: vec![],
                    clear_existing: Some(true),
                    max_hits: None,
                },
                |progress| {
                    if let OperationProgress::Tfbs(p) = progress {
                        progress_events.push(p);
                    }
                },
            )
            .unwrap();
        assert!(res.changed_seq_ids.contains(&"s".to_string()));
        assert!(!progress_events.is_empty());
        let last = progress_events.last().unwrap();
        assert_eq!(last.motif_index, last.motif_count);
        assert!((last.motif_percent - 100.0).abs() < f64::EPSILON);
        assert!((last.total_percent - 100.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_annotate_tfbs_per_tf_override_changes_quantile_threshold() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("s".to_string(), seq("TTTACGTAAACGTGGG"));
        let mut engine = GentleEngine::from_state(state);
        let strict = engine
            .apply(Operation::AnnotateTfbs {
                seq_id: "s".to_string(),
                motifs: vec!["ACGT".to_string()],
                min_llr_bits: Some(0.0),
                min_llr_quantile: Some(1.0),
                per_tf_thresholds: vec![],
                clear_existing: Some(true),
                max_hits: None,
            })
            .unwrap();
        let strict_count = engine
            .state()
            .sequences
            .get("s")
            .unwrap()
            .features()
            .iter()
            .filter(|f| {
                f.qualifier_values("gentle_generated".into())
                    .any(|v| v.eq_ignore_ascii_case("tfbs"))
            })
            .count();
        assert!(strict.changed_seq_ids.contains(&"s".to_string()));

        let relaxed = engine
            .apply(Operation::AnnotateTfbs {
                seq_id: "s".to_string(),
                motifs: vec!["ACGT".to_string()],
                min_llr_bits: Some(0.0),
                min_llr_quantile: Some(1.0),
                per_tf_thresholds: vec![TfThresholdOverride {
                    tf: "ACGT".to_string(),
                    min_llr_bits: None,
                    min_llr_quantile: Some(0.0),
                }],
                clear_existing: Some(true),
                max_hits: None,
            })
            .unwrap();
        let relaxed_count = engine
            .state()
            .sequences
            .get("s")
            .unwrap()
            .features()
            .iter()
            .filter(|f| {
                f.qualifier_values("gentle_generated".into())
                    .any(|v| v.eq_ignore_ascii_case("tfbs"))
            })
            .count();
        assert!(relaxed.changed_seq_ids.contains(&"s".to_string()));
        assert!(relaxed_count >= strict_count);
    }

    #[test]
    fn test_annotate_tfbs_max_hits_cap_and_unlimited() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("s".to_string(), seq("ACGTACGTACGTACGTACGT"));
        let mut engine = GentleEngine::from_state(state);

        let capped = engine
            .apply(Operation::AnnotateTfbs {
                seq_id: "s".to_string(),
                motifs: vec!["ACGT".to_string()],
                min_llr_bits: Some(0.0),
                min_llr_quantile: Some(0.0),
                per_tf_thresholds: vec![],
                clear_existing: Some(true),
                max_hits: Some(2),
            })
            .unwrap();
        assert!(capped.changed_seq_ids.contains(&"s".to_string()));
        let capped_count = engine
            .state()
            .sequences
            .get("s")
            .unwrap()
            .features()
            .iter()
            .filter(|f| {
                f.qualifier_values("gentle_generated".into())
                    .any(|v| v.eq_ignore_ascii_case("tfbs"))
            })
            .count();
        assert_eq!(capped_count, 2);

        let unlimited = engine
            .apply(Operation::AnnotateTfbs {
                seq_id: "s".to_string(),
                motifs: vec!["ACGT".to_string()],
                min_llr_bits: Some(0.0),
                min_llr_quantile: Some(0.0),
                per_tf_thresholds: vec![],
                clear_existing: Some(true),
                max_hits: Some(0),
            })
            .unwrap();
        assert!(unlimited.changed_seq_ids.contains(&"s".to_string()));
        let unlimited_count = engine
            .state()
            .sequences
            .get("s")
            .unwrap()
            .features()
            .iter()
            .filter(|f| {
                f.qualifier_values("gentle_generated".into())
                    .any(|v| v.eq_ignore_ascii_case("tfbs"))
            })
            .count();
        assert!(unlimited_count > capped_count);
    }

    #[test]
    fn test_render_sequence_svg_operation() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("s".to_string(), seq(&"ATGC".repeat(40)));
        let mut engine = GentleEngine::from_state(state);
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = tmp.path().with_extension("svg");
        let path_text = path.display().to_string();
        let res = engine
            .apply(Operation::RenderSequenceSvg {
                seq_id: "s".to_string(),
                mode: RenderSvgMode::Linear,
                path: path_text.clone(),
            })
            .unwrap();
        assert!(res.messages.iter().any(|m| m.contains("SVG")));
        let text = std::fs::read_to_string(path_text).unwrap();
        assert!(text.contains("<svg"));
    }

    #[test]
    fn test_render_lineage_svg_operation() {
        let mut state = ProjectState::default();
        state.sequences.insert("s".to_string(), seq("ATGCCA"));
        let mut engine = GentleEngine::from_state(state);
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = tmp.path().with_extension("svg");
        let path_text = path.display().to_string();
        let res = engine
            .apply(Operation::RenderLineageSvg {
                path: path_text.clone(),
            })
            .unwrap();
        assert!(res.messages.iter().any(|m| m.contains("lineage SVG")));
        let text = std::fs::read_to_string(path_text).unwrap();
        assert!(text.contains("<svg"));
    }

    #[test]
    fn test_render_pool_gel_svg_operation() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("a".to_string(), seq(&"ATGC".repeat(80)));
        state
            .sequences
            .insert("b".to_string(), seq(&"ATGC".repeat(150)));
        state
            .sequences
            .insert("c".to_string(), seq(&"ATGC".repeat(260)));
        let mut engine = GentleEngine::from_state(state);
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = tmp.path().with_extension("pool.gel.svg");
        let path_text = path.display().to_string();
        let res = engine
            .apply(Operation::RenderPoolGelSvg {
                inputs: vec!["a".to_string(), "b".to_string(), "c".to_string()],
                path: path_text.clone(),
                ladders: None,
            })
            .unwrap();
        assert!(res.messages.iter().any(|m| m.contains("pool gel SVG")));
        let text = std::fs::read_to_string(path_text).unwrap();
        assert!(text.contains("<svg"));
        assert!(text.contains("Pool Gel Preview"));
    }

    #[test]
    fn test_render_pool_gel_svg_operation_missing_input_fails() {
        let mut state = ProjectState::default();
        state.sequences.insert("a".to_string(), seq("ATGC"));
        let mut engine = GentleEngine::from_state(state);
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = tmp.path().with_extension("pool.gel.svg");
        let path_text = path.display().to_string();
        let err = engine
            .apply(Operation::RenderPoolGelSvg {
                inputs: vec!["missing".to_string()],
                path: path_text,
                ladders: None,
            })
            .unwrap_err();
        assert!(err.message.contains("not found"));
    }

    #[test]
    fn test_export_pool_operation() {
        let mut state = ProjectState::default();
        state.sequences.insert("a".to_string(), seq("ATGC"));
        state.sequences.insert("b".to_string(), seq("TTAA"));
        let mut engine = GentleEngine::from_state(state);
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = tmp.path().with_extension("pool.gentle.json");
        let path_text = path.display().to_string();
        let res = engine
            .apply(Operation::ExportPool {
                inputs: vec!["a".to_string(), "b".to_string()],
                path: path_text.clone(),
                pool_id: Some("pool_1".to_string()),
                human_id: Some("test pool".to_string()),
            })
            .unwrap();
        assert!(res.messages.iter().any(|m| m.contains("Wrote pool export")));
        let text = std::fs::read_to_string(path_text).unwrap();
        let v: serde_json::Value = serde_json::from_str(&text).unwrap();
        assert_eq!(v["schema"], "gentle.pool.v1");
        assert_eq!(v["pool_id"], "pool_1");
        assert_eq!(v["human_id"], "test pool");
        assert_eq!(v["member_count"], 2);
    }

    #[test]
    fn test_inspect_dna_ladders() {
        let catalog = GentleEngine::inspect_dna_ladders(None);
        assert_eq!(catalog.schema, "gentle.dna_ladders.v1");
        assert!(catalog.ladder_count > 0);
        assert_eq!(catalog.ladder_count, catalog.ladders.len());
        assert!(catalog
            .ladders
            .iter()
            .any(|ladder| ladder.name == "NEB 100bp DNA Ladder"));
    }

    #[test]
    fn test_export_dna_ladders_operation() {
        let mut engine = GentleEngine::new();
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = tmp.path().with_extension("dna.ladders.json");
        let path_text = path.display().to_string();
        let res = engine
            .apply(Operation::ExportDnaLadders {
                path: path_text.clone(),
                name_filter: Some("NEB".to_string()),
            })
            .unwrap();
        assert!(res
            .messages
            .iter()
            .any(|m| m.contains("DNA ladders catalog")));
        let text = std::fs::read_to_string(path_text).unwrap();
        let catalog: DnaLadderCatalog = serde_json::from_str(&text).unwrap();
        assert_eq!(catalog.schema, "gentle.dna_ladders.v1");
        assert!(catalog.ladder_count > 0);
        assert!(catalog
            .ladders
            .iter()
            .all(|ladder| ladder.name.to_ascii_lowercase().contains("neb")));
    }

    #[test]
    fn test_inspect_rna_ladders() {
        let catalog = GentleEngine::inspect_rna_ladders(None);
        assert_eq!(catalog.schema, "gentle.rna_ladders.v1");
        assert!(catalog.ladder_count > 0);
        assert_eq!(catalog.ladder_count, catalog.ladders.len());
        assert!(catalog
            .ladders
            .iter()
            .any(|ladder| ladder.name.contains("RNA")));
    }

    #[test]
    fn test_export_rna_ladders_operation() {
        let mut engine = GentleEngine::new();
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = tmp.path().with_extension("rna.ladders.json");
        let path_text = path.display().to_string();
        let res = engine
            .apply(Operation::ExportRnaLadders {
                path: path_text.clone(),
                name_filter: Some("NEB".to_string()),
            })
            .unwrap();
        assert!(res
            .messages
            .iter()
            .any(|m| m.contains("RNA ladders catalog")));
        let text = std::fs::read_to_string(path_text).unwrap();
        let catalog: RnaLadderCatalog = serde_json::from_str(&text).unwrap();
        assert_eq!(catalog.schema, "gentle.rna_ladders.v1");
        assert!(catalog.ladder_count > 0);
        assert!(catalog
            .ladders
            .iter()
            .all(|ladder| ladder.name.to_ascii_lowercase().contains("neb")));
    }

    #[test]
    fn test_save_file_operation_fasta() {
        let mut state = ProjectState::default();
        state.sequences.insert("s".to_string(), seq("ATGCCA"));
        let mut engine = GentleEngine::from_state(state);
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = tmp.path().with_extension("fa");
        let path_text = path.display().to_string();
        let res = engine
            .apply(Operation::SaveFile {
                seq_id: "s".to_string(),
                path: path_text.clone(),
                format: ExportFormat::Fasta,
            })
            .unwrap();
        assert!(res.changed_seq_ids.contains(&"s".to_string()));
        let text = std::fs::read_to_string(path_text).unwrap();
        assert!(text.starts_with(">"));
        assert!(text.contains("ATGCCA"));
    }

    #[test]
    fn test_save_file_operation_fasta_includes_synthetic_metadata() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "s".to_string(),
            synth_oligo("molecule=ssdna topology=linear", b"ATGCCA"),
        );
        let mut engine = GentleEngine::from_state(state);
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = tmp.path().with_extension("fa");
        let path_text = path.display().to_string();
        engine
            .apply(Operation::SaveFile {
                seq_id: "s".to_string(),
                path: path_text.clone(),
                format: ExportFormat::Fasta,
            })
            .unwrap();
        let text = std::fs::read_to_string(path_text).unwrap();
        assert!(text.contains("molecule=ssdna"));
        assert!(text.contains("topology=linear"));
    }

    #[test]
    fn test_render_sequence_svg_operation_circular() {
        let mut state = ProjectState::default();
        let mut s = seq(&"ATGC".repeat(40));
        s.set_circular(true);
        state.sequences.insert("s".to_string(), s);
        let mut engine = GentleEngine::from_state(state);
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = tmp.path().with_extension("svg");
        let path_text = path.display().to_string();
        let res = engine
            .apply(Operation::RenderSequenceSvg {
                seq_id: "s".to_string(),
                mode: RenderSvgMode::Circular,
                path: path_text.clone(),
            })
            .unwrap();
        assert!(res.messages.iter().any(|m| m.contains("SVG")));
        let text = std::fs::read_to_string(path_text).unwrap();
        assert!(text.contains("<svg"));
    }

    #[test]
    fn test_export_pool_operation_defaults() {
        let mut state = ProjectState::default();
        state.sequences.insert("a".to_string(), seq("ATGC"));
        let mut engine = GentleEngine::from_state(state);
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = tmp.path().with_extension("pool.gentle.json");
        let path_text = path.display().to_string();
        engine
            .apply(Operation::ExportPool {
                inputs: vec!["a".to_string()],
                path: path_text.clone(),
                pool_id: None,
                human_id: None,
            })
            .unwrap();
        let text = std::fs::read_to_string(path_text).unwrap();
        let v: serde_json::Value = serde_json::from_str(&text).unwrap();
        assert_eq!(v["pool_id"], "pool_export");
        assert!(v["human_id"].as_str().unwrap().starts_with("Pool("));
    }

    #[test]
    fn test_export_pool_operation_empty_inputs_fails() {
        let mut engine = GentleEngine::new();
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = tmp.path().with_extension("pool.gentle.json");
        let path_text = path.display().to_string();
        let err = engine
            .apply(Operation::ExportPool {
                inputs: vec![],
                path: path_text,
                pool_id: None,
                human_id: None,
            })
            .unwrap_err();
        assert!(err.message.contains("at least one input"));
    }

    #[test]
    fn test_export_pool_operation_missing_sequence_fails() {
        let mut engine = GentleEngine::new();
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = tmp.path().with_extension("pool.gentle.json");
        let path_text = path.display().to_string();
        let err = engine
            .apply(Operation::ExportPool {
                inputs: vec!["missing".to_string()],
                path: path_text,
                pool_id: None,
                human_id: None,
            })
            .unwrap_err();
        assert!(err.message.contains("not found"));
    }

    #[test]
    fn test_set_parameter_unknown_name_fails() {
        let mut engine = GentleEngine::new();
        let err = engine
            .apply(Operation::SetParameter {
                name: "unknown_param".to_string(),
                value: serde_json::json!(1),
            })
            .unwrap_err();
        assert!(err.message.contains("Unknown parameter"));
    }

    #[test]
    fn test_set_parameter_invalid_type_fails() {
        let mut engine = GentleEngine::new();
        let err = engine
            .apply(Operation::SetParameter {
                name: "max_fragments_per_container".to_string(),
                value: serde_json::json!("not-a-number"),
            })
            .unwrap_err();
        assert!(err.message.contains("requires a positive integer"));
    }

    #[test]
    fn test_set_parameter_zero_fails() {
        let mut engine = GentleEngine::new();
        let err = engine
            .apply(Operation::SetParameter {
                name: "max_fragments_per_container".to_string(),
                value: serde_json::json!(0),
            })
            .unwrap_err();
        assert!(err.message.contains("must be >= 1"));
    }

    #[test]
    fn test_set_parameter_feature_details_font_size_invalid_type_fails() {
        let mut engine = GentleEngine::new();
        let err = engine
            .apply(Operation::SetParameter {
                name: "feature_details_font_size".to_string(),
                value: serde_json::json!("small"),
            })
            .unwrap_err();
        assert!(err
            .message
            .contains("feature_details_font_size requires a number"));
    }

    #[test]
    fn test_set_parameter_feature_details_font_size_out_of_range_fails() {
        let mut engine = GentleEngine::new();
        let err = engine
            .apply(Operation::SetParameter {
                name: "feature_details_font_size".to_string(),
                value: serde_json::json!(3.0),
            })
            .unwrap_err();
        assert!(err
            .message
            .contains("feature_details_font_size must be between 8.0 and 24.0"));
    }

    #[test]
    fn test_containers_created_on_load_and_digest() {
        let mut engine = GentleEngine::new();
        let load = engine
            .apply(Operation::LoadFile {
                path: "test_files/pGEX_3X.fa".to_string(),
                as_id: Some("pgex".to_string()),
            })
            .unwrap();
        let seq_id = load.created_seq_ids.first().unwrap().clone();
        let latest = engine
            .state()
            .container_state
            .seq_to_latest_container
            .get(&seq_id)
            .cloned();
        assert!(latest.is_some());
        let digest = engine
            .apply(Operation::Digest {
                input: "pgex".to_string(),
                enzymes: vec!["BamHI".to_string(), "EcoRI".to_string()],
                output_prefix: Some("frag".to_string()),
            })
            .unwrap();
        assert!(digest.created_seq_ids.len() >= 2);
        let container_id = engine
            .state()
            .container_state
            .seq_to_latest_container
            .get(digest.created_seq_ids.first().unwrap())
            .unwrap();
        let container = engine
            .state()
            .container_state
            .containers
            .get(container_id)
            .unwrap();
        assert!(matches!(container.kind, ContainerKind::Pool));
        assert_eq!(container.members.len(), digest.created_seq_ids.len());
    }

    #[test]
    fn test_select_candidate_creates_selection_container_kind() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("frag".to_string(), seq(&"ATGC".repeat(40)));
        let mut engine = GentleEngine::from_state(state);
        let res = engine
            .apply(Operation::SelectCandidate {
                input: "frag".to_string(),
                criterion: "manual".to_string(),
                output_id: Some("picked".to_string()),
            })
            .unwrap();
        let picked = res.created_seq_ids.first().unwrap();
        let cid = engine
            .state()
            .container_state
            .seq_to_latest_container
            .get(picked)
            .unwrap();
        let c = engine.state().container_state.containers.get(cid).unwrap();
        assert!(matches!(c.kind, ContainerKind::Selection));
    }

    #[test]
    fn test_container_operations_map_to_core_ops() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("a".to_string(), seq(&"ATGC".repeat(40)));
        state
            .sequences
            .insert("b".to_string(), seq(&"TTAA".repeat(40)));
        let mut engine = GentleEngine::from_state(state);
        let merge = engine
            .apply(Operation::MergeContainersById {
                container_ids: vec![
                    engine
                        .state()
                        .container_state
                        .seq_to_latest_container
                        .get("a")
                        .unwrap()
                        .clone(),
                    engine
                        .state()
                        .container_state
                        .seq_to_latest_container
                        .get("b")
                        .unwrap()
                        .clone(),
                ],
                output_prefix: Some("pool".to_string()),
            })
            .unwrap();
        assert_eq!(merge.created_seq_ids.len(), 2);

        let pool_container = engine
            .state()
            .container_state
            .seq_to_latest_container
            .get(merge.created_seq_ids.first().unwrap())
            .unwrap()
            .clone();
        let lig = engine
            .apply(Operation::LigationContainer {
                container_id: pool_container.clone(),
                circularize_if_possible: false,
                output_id: None,
                protocol: LigationProtocol::Blunt,
                output_prefix: Some("lig".to_string()),
                unique: None,
            })
            .unwrap();
        assert!(!lig.created_seq_ids.is_empty());

        let filtered = engine
            .apply(Operation::FilterContainerByMolecularWeight {
                container_id: pool_container,
                min_bp: 120,
                max_bp: 400,
                error: 0.0,
                unique: false,
                output_prefix: Some("mw".to_string()),
            })
            .unwrap();
        assert!(!filtered.created_seq_ids.is_empty());
    }

    #[test]
    fn test_digest_container_and_state_summary_include_containers() {
        let mut engine = GentleEngine::new();
        engine
            .apply(Operation::LoadFile {
                path: "test_files/pGEX_3X.fa".to_string(),
                as_id: Some("pgex".to_string()),
            })
            .unwrap();
        let cid = engine
            .state()
            .container_state
            .seq_to_latest_container
            .get("pgex")
            .unwrap()
            .clone();
        let res = engine
            .apply(Operation::DigestContainer {
                container_id: cid,
                enzymes: vec!["BamHI".to_string(), "EcoRI".to_string()],
                output_prefix: Some("dig".to_string()),
            })
            .unwrap();
        assert!(!res.created_seq_ids.is_empty());
        let summary = engine.summarize_state();
        assert!(summary.container_count > 0);
        assert!(!summary.containers.is_empty());
    }

    #[test]
    fn test_prepare_genome_and_extract_region_operations() {
        let td = tempdir().unwrap();
        let root = td.path();
        let fasta_gz = root.join("toy.fa.gz");
        let ann_gz = root.join("toy.gtf.gz");
        write_gzip(&fasta_gz, ">chr1\nACGT\nACGT\nACGT\n");
        write_gzip(
            &ann_gz,
            "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\";\n",
        );
        let cache_dir = root.join("cache");
        let catalog_path = root.join("catalog.json");
        let catalog_json = format!(
            r#"{{
  "ToyGenome": {{
    "description": "toy genome",
    "sequence_remote": "{}",
    "annotations_remote": "{}",
    "cache_dir": "{}"
  }}
}}"#,
            file_url(&fasta_gz),
            file_url(&ann_gz),
            cache_dir.display()
        );
        fs::write(&catalog_path, catalog_json).unwrap();

        let mut engine = GentleEngine::new();
        let prep = engine
            .apply(Operation::PrepareGenome {
                genome_id: "ToyGenome".to_string(),
                catalog_path: Some(catalog_path.to_string_lossy().to_string()),
                cache_dir: None,
            })
            .unwrap();
        assert!(prep
            .messages
            .iter()
            .any(|m| m.contains("Prepared genome 'ToyGenome'")));
        let catalog_path_str = catalog_path.to_string_lossy().to_string();
        let catalog_names = GentleEngine::list_reference_genomes(Some(&catalog_path_str)).unwrap();
        assert!(catalog_names.contains(&"ToyGenome".to_string()));
        let prepared =
            GentleEngine::is_reference_genome_prepared(Some(&catalog_path_str), "ToyGenome", None)
                .unwrap();
        assert!(prepared);
        let listed_genes =
            GentleEngine::list_reference_genome_genes(Some(&catalog_path_str), "ToyGenome", None)
                .unwrap();
        assert_eq!(listed_genes.len(), 1);

        let extract = engine
            .apply(Operation::ExtractGenomeRegion {
                genome_id: "ToyGenome".to_string(),
                chromosome: "chr1".to_string(),
                start_1based: 3,
                end_1based: 10,
                output_id: Some("toy_slice".to_string()),
                catalog_path: Some(catalog_path_str.clone()),
                cache_dir: None,
            })
            .unwrap();
        assert_eq!(extract.created_seq_ids, vec!["toy_slice".to_string()]);
        let loaded = engine.state().sequences.get("toy_slice").unwrap();
        assert_eq!(loaded.get_forward_string(), "GTACGTAC");

        let extract_gene = engine
            .apply(Operation::ExtractGenomeGene {
                genome_id: "ToyGenome".to_string(),
                gene_query: "MYGENE".to_string(),
                occurrence: None,
                output_id: Some("toy_gene".to_string()),
                catalog_path: Some(catalog_path_str),
                cache_dir: None,
            })
            .unwrap();
        assert_eq!(extract_gene.created_seq_ids, vec!["toy_gene".to_string()]);
        let loaded_gene = engine.state().sequences.get("toy_gene").unwrap();
        assert_eq!(loaded_gene.get_forward_string(), "ACGTACGTACGT");
        let provenance = engine
            .state()
            .metadata
            .get("provenance")
            .and_then(|v| v.as_object())
            .expect("provenance metadata object");
        let extractions = provenance
            .get("genome_extractions")
            .and_then(|v| v.as_array())
            .expect("genome_extractions array");
        assert_eq!(extractions.len(), 2);
        assert!(extractions.iter().any(|entry| {
            entry
                .get("seq_id")
                .and_then(|v| v.as_str())
                .map(|v| v == "toy_slice")
                .unwrap_or(false)
                && entry
                    .get("operation")
                    .and_then(|v| v.as_str())
                    .map(|v| v == "ExtractGenomeRegion")
                    .unwrap_or(false)
        }));
        assert!(extractions.iter().any(|entry| {
            entry
                .get("seq_id")
                .and_then(|v| v.as_str())
                .map(|v| v == "toy_gene")
                .unwrap_or(false)
                && entry
                    .get("operation")
                    .and_then(|v| v.as_str())
                    .map(|v| v == "ExtractGenomeGene")
                    .unwrap_or(false)
                && entry
                    .get("sequence_sha1")
                    .and_then(|v| v.as_str())
                    .map(|v| !v.is_empty())
                    .unwrap_or(false)
        }));
    }

    #[test]
    fn test_import_genome_bed_track_supports_plain_and_gzip() {
        let td = tempdir().unwrap();
        let root = td.path();
        let fasta_gz = root.join("toy.fa.gz");
        let ann_gz = root.join("toy.gtf.gz");
        write_gzip(&fasta_gz, ">chr1\nACGT\nACGT\nACGT\n");
        write_gzip(
            &ann_gz,
            "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\";\n",
        );
        let cache_dir = root.join("cache");
        let catalog_path = root.join("catalog.json");
        let catalog_json = format!(
            r#"{{
  "ToyGenome": {{
    "description": "toy genome",
    "sequence_remote": "{}",
    "annotations_remote": "{}",
    "cache_dir": "{}"
  }}
}}"#,
            file_url(&fasta_gz),
            file_url(&ann_gz),
            cache_dir.display()
        );
        fs::write(&catalog_path, catalog_json).unwrap();
        let catalog_path_str = catalog_path.to_string_lossy().to_string();

        let mut engine = GentleEngine::new();
        engine
            .apply(Operation::PrepareGenome {
                genome_id: "ToyGenome".to_string(),
                catalog_path: Some(catalog_path_str.clone()),
                cache_dir: None,
            })
            .unwrap();
        engine
            .apply(Operation::ExtractGenomeRegion {
                genome_id: "ToyGenome".to_string(),
                chromosome: "chr1".to_string(),
                start_1based: 3,
                end_1based: 10,
                output_id: Some("toy_slice".to_string()),
                catalog_path: Some(catalog_path_str),
                cache_dir: None,
            })
            .unwrap();

        let bed_path = root.join("signals.bed");
        fs::write(
            &bed_path,
            "track name=toy\nchr1\t1\t4\tpeak_a\t42\t+\nchr1\t5\t12\tpeak_b\t7\t-\nchr2\t1\t4\twrong_chr\t50\t+\nchr1\tbad\t9\tbroken\n",
        )
        .unwrap();
        let plain = engine
            .apply(Operation::ImportGenomeBedTrack {
                seq_id: "toy_slice".to_string(),
                path: bed_path.to_string_lossy().to_string(),
                track_name: Some("chipseq_plain".to_string()),
                min_score: None,
                max_score: None,
                clear_existing: Some(true),
            })
            .unwrap();
        assert!(plain.changed_seq_ids.contains(&"toy_slice".to_string()));
        assert!(plain.warnings.iter().any(|w| w.contains("BED line")));

        let dna_plain = engine.state().sequences.get("toy_slice").unwrap();
        let generated_plain: Vec<_> = dna_plain
            .features()
            .iter()
            .filter(|f| {
                f.qualifier_values("gentle_generated".into())
                    .any(|v| v.eq_ignore_ascii_case(GENOME_BED_TRACK_GENERATED_TAG))
            })
            .collect();
        assert_eq!(generated_plain.len(), 2);
        assert!(generated_plain.iter().any(|f| {
            f.qualifier_values("label".into())
                .next()
                .map(|v| v.contains("peak_a"))
                .unwrap_or(false)
        }));

        let bed_gz = root.join("signals.bed.gz");
        write_gzip(
            &bed_gz,
            "chr1\t1\t4\tpeak_a\t42\t+\nchr1\t5\t12\tpeak_b\t.\t-\n",
        );
        let gz = engine
            .apply(Operation::ImportGenomeBedTrack {
                seq_id: "toy_slice".to_string(),
                path: bed_gz.to_string_lossy().to_string(),
                track_name: Some("chipseq_gz".to_string()),
                min_score: Some(10.0),
                max_score: None,
                clear_existing: Some(true),
            })
            .unwrap();
        assert!(gz.changed_seq_ids.contains(&"toy_slice".to_string()));
        assert!(gz.warnings.iter().any(|w| w.contains("score column")));

        let dna_gz = engine.state().sequences.get("toy_slice").unwrap();
        let generated_gz: Vec<_> = dna_gz
            .features()
            .iter()
            .filter(|f| {
                f.qualifier_values("gentle_generated".into())
                    .any(|v| v.eq_ignore_ascii_case(GENOME_BED_TRACK_GENERATED_TAG))
            })
            .collect();
        assert_eq!(generated_gz.len(), 1);
        assert!(generated_gz[0]
            .qualifier_values("label".into())
            .next()
            .map(|v| v.contains("peak_a"))
            .unwrap_or(false));
    }

    #[cfg(unix)]
    #[test]
    fn test_import_genome_bigwig_track_uses_converter_and_filters_scores() {
        let td = tempdir().unwrap();
        let root = td.path();
        let fasta_gz = root.join("toy.fa.gz");
        let ann_gz = root.join("toy.gtf.gz");
        write_gzip(&fasta_gz, ">chr1\nACGT\nACGT\nACGT\n");
        write_gzip(
            &ann_gz,
            "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\";\n",
        );
        let cache_dir = root.join("cache");
        let catalog_path = root.join("catalog.json");
        let catalog_json = format!(
            r#"{{
  "ToyGenome": {{
    "description": "toy genome",
    "sequence_remote": "{}",
    "annotations_remote": "{}",
    "cache_dir": "{}"
  }}
}}"#,
            file_url(&fasta_gz),
            file_url(&ann_gz),
            cache_dir.display()
        );
        fs::write(&catalog_path, catalog_json).unwrap();
        let catalog_path_str = catalog_path.to_string_lossy().to_string();

        let mut engine = GentleEngine::new();
        engine
            .apply(Operation::PrepareGenome {
                genome_id: "ToyGenome".to_string(),
                catalog_path: Some(catalog_path_str.clone()),
                cache_dir: None,
            })
            .unwrap();
        engine
            .apply(Operation::ExtractGenomeRegion {
                genome_id: "ToyGenome".to_string(),
                chromosome: "chr1".to_string(),
                start_1based: 3,
                end_1based: 10,
                output_id: Some("toy_slice".to_string()),
                catalog_path: Some(catalog_path_str),
                cache_dir: None,
            })
            .unwrap();

        let bed_path = root.join("baseline.bed");
        fs::write(&bed_path, "chr1\t1\t4\tpeak_a\t42\t+\n").unwrap();
        engine
            .apply(Operation::ImportGenomeBedTrack {
                seq_id: "toy_slice".to_string(),
                path: bed_path.to_string_lossy().to_string(),
                track_name: Some("baseline".to_string()),
                min_score: None,
                max_score: None,
                clear_existing: Some(true),
            })
            .unwrap();

        let converted_bedgraph = root.join("signals.from_bigwig.bedgraph");
        fs::write(
            &converted_bedgraph,
            "chr1\t1\t4\t0.5\nchr1\t5\t12\t2.0\nchr2\t1\t4\t7.0\nchr1\tbad\t9\tbroken\n",
        )
        .unwrap();
        let fake_converter = install_fake_bigwig_to_bedgraph(root, &converted_bedgraph);
        let _converter_guard = EnvVarGuard::set(BIGWIG_TO_BEDGRAPH_ENV_BIN, &fake_converter);

        let fake_bigwig = root.join("signals.bw");
        fs::write(&fake_bigwig, "placeholder").unwrap();
        let result = engine
            .apply(Operation::ImportGenomeBigWigTrack {
                seq_id: "toy_slice".to_string(),
                path: fake_bigwig.to_string_lossy().to_string(),
                track_name: Some("signal".to_string()),
                min_score: Some(1.0),
                max_score: Some(3.0),
                clear_existing: Some(true),
            })
            .unwrap();
        assert!(result.changed_seq_ids.contains(&"toy_slice".to_string()));
        assert!(result.warnings.iter().any(|w| w.contains("bedGraph line")));

        let dna = engine.state().sequences.get("toy_slice").unwrap();
        let bed_features: Vec<_> = dna
            .features()
            .iter()
            .filter(|f| {
                f.qualifier_values("gentle_generated".into())
                    .any(|v| v.eq_ignore_ascii_case(GENOME_BED_TRACK_GENERATED_TAG))
            })
            .collect();
        assert_eq!(bed_features.len(), 0);
        let bigwig_features: Vec<_> = dna
            .features()
            .iter()
            .filter(|f| GentleEngine::is_generated_genome_bigwig_feature(f))
            .collect();
        assert_eq!(bigwig_features.len(), 1);
        assert_eq!(
            bigwig_features[0]
                .qualifier_values("score".into())
                .next()
                .unwrap_or_default(),
            "2.000000"
        );
    }

    #[test]
    fn test_import_genome_vcf_track_supports_multiallelic_and_qual_filters() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("toy_slice".to_string(), seq("ACGTACGTACGTACGTACGT"));
        state.metadata.insert(
            PROVENANCE_METADATA_KEY.to_string(),
            serde_json::json!({
                GENOME_EXTRACTIONS_METADATA_KEY: [
                    {
                        "seq_id": "toy_slice",
                        "recorded_at_unix_ms": 1,
                        "operation": "ExtractGenomeRegion",
                        "genome_id": "ToyGenome",
                        "catalog_path": "synthetic",
                        "cache_dir": null,
                        "chromosome": "chr1",
                        "start_1based": 100,
                        "end_1based": 119,
                        "gene_query": null,
                        "occurrence": null,
                        "gene_id": null,
                        "gene_name": null,
                        "strand": null,
                        "anchor_strand": "+",
                        "sequence_source_type": "synthetic",
                        "annotation_source_type": "synthetic",
                        "sequence_source": "synthetic",
                        "annotation_source": "synthetic",
                        "sequence_sha1": null,
                        "annotation_sha1": null
                    }
                ]
            }),
        );
        let mut engine = GentleEngine::from_state(state);

        let td = tempdir().unwrap();
        let vcf_path = td.path().join("variants.vcf");
        fs::write(
            &vcf_path,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t101\trs1\tA\tG\t50\tPASS\tAC=1\nchr1\t105\t.\tA\tT,<DEL>\t.\tq10\tDP=5\nchr2\t101\trs2\tA\tC\t60\tPASS\tAC=1\nchr1\tbad\tbroken\n",
        )
        .unwrap();

        let plain = engine
            .apply(Operation::ImportGenomeVcfTrack {
                seq_id: "toy_slice".to_string(),
                path: vcf_path.to_string_lossy().to_string(),
                track_name: Some("variants".to_string()),
                min_score: None,
                max_score: None,
                clear_existing: Some(true),
            })
            .unwrap();
        assert!(plain.changed_seq_ids.contains(&"toy_slice".to_string()));
        assert!(plain.warnings.iter().any(|w| w.contains("VCF line")));

        let dna_plain = engine.state().sequences.get("toy_slice").unwrap();
        let plain_features: Vec<_> = dna_plain
            .features()
            .iter()
            .filter(|f| GentleEngine::is_generated_genome_vcf_feature(f))
            .collect();
        assert_eq!(plain_features.len(), 3);
        assert!(plain_features
            .iter()
            .any(|f| { f.qualifier_values("vcf_alt".into()).any(|v| v == "<DEL>") }));

        let vcf_gz = td.path().join("variants.vcf.gz");
        write_gzip(
            &vcf_gz,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t110\trs3\tC\tT\t25\tPASS\tAC=1\nchr1\t111\trs4\tC\tG\t.\tPASS\tAC=1\n",
        );
        let filtered = engine
            .apply(Operation::ImportGenomeVcfTrack {
                seq_id: "toy_slice".to_string(),
                path: vcf_gz.to_string_lossy().to_string(),
                track_name: Some("variants_gz".to_string()),
                min_score: Some(20.0),
                max_score: Some(30.0),
                clear_existing: Some(true),
            })
            .unwrap();
        assert!(filtered.changed_seq_ids.contains(&"toy_slice".to_string()));
        assert!(filtered
            .warnings
            .iter()
            .any(|w| w.contains("QUAL-based score filters")));

        let dna_filtered = engine.state().sequences.get("toy_slice").unwrap();
        let filtered_features: Vec<_> = dna_filtered
            .features()
            .iter()
            .filter(|f| GentleEngine::is_generated_genome_vcf_feature(f))
            .collect();
        assert_eq!(filtered_features.len(), 1);
        assert_eq!(
            filtered_features[0]
                .qualifier_values("vcf_id".into())
                .next()
                .unwrap_or_default(),
            "rs3"
        );
    }

    #[test]
    fn test_prepare_helper_genome_via_genbank_accession_and_extract() {
        let _guard = crate::genomes::genbank_env_lock().lock().unwrap();
        let td = tempdir().unwrap();
        let root = td.path();

        let mock_dir = root.join("mock");
        fs::create_dir_all(&mock_dir).unwrap();
        fs::copy("test_files/pGEX_3X.fa", mock_dir.join("L09137.fasta")).unwrap();
        fs::copy("test_files/pGEX-3X.gb", mock_dir.join("L09137.gbwithparts")).unwrap();

        let efetch_template = format!("file://{}/{{accession}}.{{rettype}}", mock_dir.display());
        let _efetch_env = EnvVarGuard::set("GENTLE_NCBI_EFETCH_URL", &efetch_template);

        let cache_dir = root.join("cache");
        let catalog_path = root.join("helper_catalog.json");
        let catalog_json = format!(
            r#"{{
  "Helper pUC19": {{
    "description": "helper vector from GenBank accession",
    "genbank_accession": "L09137",
    "cache_dir": "{}"
  }}
}}"#,
            cache_dir.display()
        );
        fs::write(&catalog_path, catalog_json).unwrap();
        let catalog_path_str = catalog_path.to_string_lossy().to_string();

        let mut engine = GentleEngine::new();
        let prep = engine
            .apply(Operation::PrepareGenome {
                genome_id: "Helper pUC19".to_string(),
                catalog_path: Some(catalog_path_str.clone()),
                cache_dir: None,
            })
            .unwrap();
        assert!(prep
            .messages
            .iter()
            .any(|m| m.contains("[genbank_accession]")));

        let plan = GentleEngine::describe_reference_genome_sources(
            Some(&catalog_path_str),
            "Helper pUC19",
            None,
        )
        .unwrap();
        assert_eq!(plan.sequence_source_type, "genbank_accession");
        assert_eq!(plan.annotation_source_type, "genbank_accession");

        let genes = GentleEngine::list_reference_genome_genes(
            Some(&catalog_path_str),
            "Helper pUC19",
            None,
        )
        .unwrap();
        assert!(!genes.is_empty());
        assert!(genes.iter().any(|g| {
            g.gene_name
                .as_ref()
                .map(|v| v.eq_ignore_ascii_case("bla"))
                .unwrap_or(false)
                || g.gene_id
                    .as_ref()
                    .map(|v| v.eq_ignore_ascii_case("bla"))
                    .unwrap_or(false)
        }));

        let extract_gene = engine
            .apply(Operation::ExtractGenomeGene {
                genome_id: "Helper pUC19".to_string(),
                gene_query: "bla".to_string(),
                occurrence: Some(1),
                output_id: Some("helper_bla".to_string()),
                catalog_path: Some(catalog_path_str.clone()),
                cache_dir: None,
            })
            .unwrap();
        assert_eq!(extract_gene.created_seq_ids, vec!["helper_bla".to_string()]);
        let seq = engine.state().sequences.get("helper_bla").unwrap();
        assert!(seq.len() > 0);

        let extract_region = engine
            .apply(Operation::ExtractGenomeRegion {
                genome_id: "Helper pUC19".to_string(),
                chromosome: "U13852.1".to_string(),
                start_1based: 1,
                end_1based: 40,
                output_id: Some("helper_head".to_string()),
                catalog_path: Some(catalog_path_str),
                cache_dir: None,
            })
            .unwrap();
        assert_eq!(
            extract_region.created_seq_ids,
            vec!["helper_head".to_string()]
        );
    }

    #[test]
    fn test_sync_tracked_genome_track_subscriptions_only_applies_to_new_anchors() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("anch1".to_string(), seq("ACGTACGTACGT"));
        state.metadata.insert(
            PROVENANCE_METADATA_KEY.to_string(),
            serde_json::json!({
                GENOME_EXTRACTIONS_METADATA_KEY: [
                    {
                        "seq_id": "anch1",
                        "recorded_at_unix_ms": 1,
                        "operation": "ExtractGenomeRegion",
                        "genome_id": "ToyGenome",
                        "catalog_path": "synthetic",
                        "cache_dir": null,
                        "chromosome": "chr1",
                        "start_1based": 1,
                        "end_1based": 12,
                        "gene_query": null,
                        "occurrence": null,
                        "gene_id": null,
                        "gene_name": null,
                        "strand": null,
                        "anchor_strand": "+",
                        "sequence_source_type": null,
                        "annotation_source_type": null,
                        "sequence_source": null,
                        "annotation_source": null,
                        "sequence_sha1": null,
                        "annotation_sha1": null
                    }
                ]
            }),
        );
        let mut engine = GentleEngine::from_state(state);

        let td = tempdir().unwrap();
        let bed_path = td.path().join("tracked.bed");
        std::fs::write(&bed_path, "chr1\t0\t4\tpeak1\t100\t+\n").unwrap();

        let inserted = engine
            .add_genome_track_subscription(GenomeTrackSubscription {
                source: GenomeTrackSource::Bed,
                path: bed_path.to_string_lossy().to_string(),
                track_name: Some("tracked".to_string()),
                min_score: None,
                max_score: None,
                clear_existing: true,
            })
            .unwrap();
        assert!(inserted);

        let first = engine
            .sync_tracked_genome_track_subscriptions(false)
            .unwrap();
        assert_eq!(first.target_sequences, 1);
        assert_eq!(first.applied_imports, 1);
        assert_eq!(first.failed_imports, 0);

        let anch1_features_after_first = engine
            .state()
            .sequences
            .get("anch1")
            .unwrap()
            .features()
            .iter()
            .filter(|f| GentleEngine::is_generated_genome_bed_feature(f))
            .count();
        assert_eq!(anch1_features_after_first, 1);

        engine
            .state_mut()
            .sequences
            .insert("anch2".to_string(), seq("ACGTACGTACGT"));
        let provenance = engine
            .state_mut()
            .metadata
            .get_mut(PROVENANCE_METADATA_KEY)
            .and_then(|v| v.as_object_mut())
            .unwrap();
        let records = provenance
            .get_mut(GENOME_EXTRACTIONS_METADATA_KEY)
            .and_then(|v| v.as_array_mut())
            .unwrap();
        records.push(serde_json::json!({
            "seq_id": "anch2",
            "recorded_at_unix_ms": 2,
            "operation": "ExtractGenomeRegion",
            "genome_id": "ToyGenome",
            "catalog_path": "synthetic",
            "cache_dir": null,
            "chromosome": "chr1",
            "start_1based": 1,
            "end_1based": 12,
            "gene_query": null,
            "occurrence": null,
            "gene_id": null,
            "gene_name": null,
            "strand": null,
            "anchor_strand": "+",
            "sequence_source_type": null,
            "annotation_source_type": null,
            "sequence_source": null,
            "annotation_source": null,
            "sequence_sha1": null,
            "annotation_sha1": null
        }));

        let second = engine
            .sync_tracked_genome_track_subscriptions(true)
            .unwrap();
        assert_eq!(second.target_sequences, 1);
        assert_eq!(second.applied_imports, 1);
        assert_eq!(second.failed_imports, 0);

        let anch1_features_after_second = engine
            .state()
            .sequences
            .get("anch1")
            .unwrap()
            .features()
            .iter()
            .filter(|f| GentleEngine::is_generated_genome_bed_feature(f))
            .count();
        assert_eq!(anch1_features_after_second, 1);
        let anch2_features = engine
            .state()
            .sequences
            .get("anch2")
            .unwrap()
            .features()
            .iter()
            .filter(|f| GentleEngine::is_generated_genome_bed_feature(f))
            .count();
        assert_eq!(anch2_features, 1);
    }

    #[test]
    fn test_sync_tracked_genome_track_subscriptions_does_not_persist_empty_known_anchor_key() {
        let mut engine = GentleEngine::new();
        let report = engine
            .sync_tracked_genome_track_subscriptions(true)
            .expect("sync should succeed for empty state");
        assert_eq!(report.subscriptions_considered, 0);
        assert_eq!(report.target_sequences, 0);
        assert_eq!(report.applied_imports, 0);
        assert!(!engine
            .state()
            .metadata
            .contains_key(GENOME_TRACK_KNOWN_ANCHORS_METADATA_KEY));
    }

    #[test]
    fn test_import_blast_hits_track_operation_adds_features_and_clears_previous() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("query".to_string(), seq("ACGTACGTACGTACGTACGTACGTACGT"));
        let mut engine = GentleEngine::from_state(state);

        let first = engine
            .apply(Operation::ImportBlastHitsTrack {
                seq_id: "query".to_string(),
                hits: vec![
                    BlastHitFeatureInput {
                        subject_id: "chr1".to_string(),
                        query_start_1based: 1,
                        query_end_1based: 8,
                        subject_start_1based: 100,
                        subject_end_1based: 107,
                        identity_percent: 99.5,
                        bit_score: 42.0,
                        evalue: 1e-8,
                        query_coverage_percent: Some(100.0),
                    },
                    BlastHitFeatureInput {
                        subject_id: "chr2".to_string(),
                        query_start_1based: 20,
                        query_end_1based: 40,
                        subject_start_1based: 500,
                        subject_end_1based: 480,
                        identity_percent: 95.0,
                        bit_score: 33.0,
                        evalue: 1e-4,
                        query_coverage_percent: Some(75.0),
                    },
                ],
                track_name: Some("blast_hits_demo".to_string()),
                clear_existing: Some(true),
            })
            .unwrap();
        assert!(first.changed_seq_ids.contains(&"query".to_string()));
        assert!(first.warnings.iter().any(|w| w.contains("was clamped")));

        let dna = engine.state().sequences.get("query").unwrap();
        let blast_features: Vec<_> = dna
            .features()
            .iter()
            .filter(|f| GentleEngine::is_generated_blast_hit_feature(f))
            .collect();
        assert_eq!(blast_features.len(), 2);
        assert!(blast_features.iter().any(|f| {
            f.qualifier_values("blast_subject_id".into())
                .any(|v| v == "chr1")
        }));
        assert!(blast_features.iter().any(|f| {
            f.qualifier_values("blast_subject_id".into())
                .any(|v| v == "chr2")
        }));

        let second = engine
            .apply(Operation::ImportBlastHitsTrack {
                seq_id: "query".to_string(),
                hits: vec![BlastHitFeatureInput {
                    subject_id: "chr3".to_string(),
                    query_start_1based: 5,
                    query_end_1based: 12,
                    subject_start_1based: 1000,
                    subject_end_1based: 1007,
                    identity_percent: 98.0,
                    bit_score: 50.0,
                    evalue: 1e-12,
                    query_coverage_percent: Some(90.0),
                }],
                track_name: Some("blast_hits_demo".to_string()),
                clear_existing: Some(true),
            })
            .unwrap();
        assert!(second.changed_seq_ids.contains(&"query".to_string()));

        let dna = engine.state().sequences.get("query").unwrap();
        let blast_features_after_clear: Vec<_> = dna
            .features()
            .iter()
            .filter(|f| GentleEngine::is_generated_blast_hit_feature(f))
            .collect();
        assert_eq!(blast_features_after_clear.len(), 1);
        assert_eq!(
            blast_features_after_clear[0]
                .qualifier_values("blast_subject_id".into())
                .next()
                .unwrap_or_default(),
            "chr3"
        );
    }
}
