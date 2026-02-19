use crate::{
    app::GENtleApp,
    dna_sequence::DNAsequence,
    enzymes::active_restriction_enzymes,
    genomes::{
        GenomeCatalog, GenomeGeneRecord, GenomeSourcePlan, PrepareGenomeProgress,
        PrepareGenomeReport, PreparedGenomeInspection, DEFAULT_GENOME_CATALOG_PATH,
        DEFAULT_HELPER_GENOME_CATALOG_PATH,
    },
    iupac_code::IupacCode,
    lineage_export::export_lineage_svg,
    methylation_sites::MethylationMode,
    pool_gel::{build_pool_gel_layout, export_pool_gel_svg},
    render_export::{export_circular_svg, export_linear_svg},
    restriction_enzyme::RestrictionEnzyme,
    tf_motifs, DNA_LADDERS,
};
use serde::{Deserialize, Serialize};
use serde_json::json;
use std::{
    collections::hash_map::DefaultHasher,
    collections::{HashMap, HashSet},
    error::Error,
    fmt,
    fs::File,
    hash::{Hash, Hasher},
    io::Write,
    path::Path,
    time::Instant,
};

pub type SeqId = String;
pub type OpId = String;
pub type RunId = String;
pub type NodeId = String;
pub type ContainerId = String;
const PROVENANCE_METADATA_KEY: &str = "provenance";
const GENOME_EXTRACTIONS_METADATA_KEY: &str = "genome_extractions";

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
            show_open_reading_frames: true,
            show_methylation_sites: false,
            linear_view_start_bp: 0,
            linear_view_span_bp: 0,
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
    pub sequence_source_type: Option<String>,
    pub annotation_source_type: Option<String>,
    pub sequence_source: Option<String>,
    pub annotation_source: Option<String>,
    pub sequence_sha1: Option<String>,
    pub annotation_sha1: Option<String>,
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

    pub fn capabilities() -> Capabilities {
        Capabilities {
            protocol_version: "v1".to_string(),
            supported_operations: vec![
                "LoadFile".to_string(),
                "SaveFile".to_string(),
                "RenderSequenceSvg".to_string(),
                "RenderLineageSvg".to_string(),
                "RenderPoolGelSvg".to_string(),
                "ExportDnaLadders".to_string(),
                "ExportPool".to_string(),
                "PrepareGenome".to_string(),
                "ExtractGenomeRegion".to_string(),
                "ExtractGenomeGene".to_string(),
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
                    length_bp: band.length_bp,
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
        format!(
            "Prepared genome '{}' ({status}). cache='{}' sequence='{}' [{}], annotation='{}' [{}]",
            genome_id,
            cache_dir
                .map(str::trim)
                .filter(|v| !v.is_empty())
                .unwrap_or("catalog/default"),
            report.sequence_path,
            sequence_type,
            report.annotation_path,
            annotation_type
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
            Operation::SelectCandidate { .. } => Some("Selected candidate".to_string()),
            Operation::FilterByMolecularWeight { .. } => {
                Some("Molecular-weight filtered".to_string())
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
                "note".into(),
                Some(format!(
                    "tf_id={tf_id}; llr_bits={llr_bits:.4}; llr_quantile={llr_quantile:.4}; true_log_odds_bits={true_log_odds_bits:.4}; true_log_odds_quantile={true_log_odds_quantile:.4}"
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
                let default_min_llr_quantile = min_llr_quantile.unwrap_or(0.95);
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
                        "TF '{}' annotated {} hit(s) with min_llr_bits={} and min_llr_quantile={}",
                        tf_id, kept, eff_bits, eff_quantile
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
                                "TFBS hit cap ({limit}) reached; skipping remaining motif scans"
                            ));
                        }
                        break 'motif_loop;
                    }
                }

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
            f.qualifier_values("llr_bits".into()).next().is_some()
                && f.qualifier_values("llr_quantile".into()).next().is_some()
                && f.qualifier_values("true_log_odds_bits".into())
                    .next()
                    .is_some()
                && f.qualifier_values("true_log_odds_quantile".into())
                    .next()
                    .is_some()
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
}
