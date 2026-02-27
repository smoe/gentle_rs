//! Shared deterministic operation engine and project-state contract.
//!
//! This module is the single execution layer for GENtle biology/workflow
//! behavior. GUI, CLI, JavaScript, and Lua adapters are expected to translate
//! user intent into operations defined here instead of re-implementing domain
//! logic in frontend code.
//!
//! Core responsibilities:
//! - Define persistent state contracts (`ProjectState`, lineage, containers,
//!   display settings, and metadata schemas).
//! - Define typed operation/workflow payloads and structured results/warnings.
//! - Execute cloning/editing/analysis/render-export operations through one
//!   deterministic path.
//! - Maintain provenance and operation identifiers used by replay, audit, and
//!   adapter-equivalent automation.
//!
//! Invariants:
//! - Behavioral parity across adapters: same operation input yields equivalent
//!   state/result output regardless of entry point.
//! - Operation semantics are explicit and machine-readable; hidden adapter state
//!   must not change core biological outcomes.
//! - Long-running operations use cooperative progress callbacks and explicit
//!   cancellation signaling where supported.

use crate::{
    DNA_LADDERS, RNA_LADDERS,
    app::GENtleApp,
    dna_sequence::DNAsequence,
    enzymes::active_restriction_enzymes,
    feature_location::{collect_location_ranges_usize, feature_is_reverse},
    genomes::{
        DEFAULT_GENOME_CATALOG_PATH, DEFAULT_HELPER_GENOME_CATALOG_PATH, GenomeBlastReport,
        GenomeCatalog, GenomeGeneRecord, GenomeSourcePlan, PrepareGenomeProgress,
        PrepareGenomeReport, PreparedGenomeInspection, is_prepare_cancelled_error,
    },
    iupac_code::IupacCode,
    lineage_export::export_lineage_svg,
    methylation_sites::MethylationMode,
    pool_gel::{GelSampleInput, build_serial_gel_layout, export_pool_gel_svg},
    render_export::{export_circular_svg, export_linear_svg},
    render_feature_expert::render_feature_expert_svg,
    restriction_enzyme::{RestrictionEnzyme, RestrictionEnzymeKey},
    rna_structure::{self, RnaStructureError, RnaStructureSvgReport, RnaStructureTextReport},
    tf_motifs,
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
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    process::Command,
    time::{Duration, Instant},
};
use tempfile::NamedTempFile;

pub type SeqId = String;
pub type OpId = String;
pub type RunId = String;
pub type NodeId = String;
pub type ContainerId = String;
pub use crate::feature_expert::{
    FeatureExpertTarget, FeatureExpertView, RESTRICTION_EXPERT_INSTRUCTION,
    RestrictionSiteExpertView, SPLICING_EXPERT_INSTRUCTION, SplicingBoundaryMarker,
    SplicingEventSummary, SplicingExonSummary, SplicingExpertView, SplicingJunctionArc,
    SplicingMatrixRow, SplicingRange, SplicingTranscriptLane, TFBS_EXPERT_INSTRUCTION,
    TfbsExpertColumn, TfbsExpertView,
};
const PROVENANCE_METADATA_KEY: &str = "provenance";
const GENOME_EXTRACTIONS_METADATA_KEY: &str = "genome_extractions";
pub const GENOME_TRACK_SUBSCRIPTIONS_METADATA_KEY: &str = "genome_bed_track_subscriptions";
const GENOME_TRACK_KNOWN_ANCHORS_METADATA_KEY: &str = "genome_track_autosync_known_anchors";
pub const CANDIDATE_SETS_METADATA_KEY: &str = "candidate_sets";
const CANDIDATE_SETS_SCHEMA: &str = "gentle.candidate_sets.v1";
const CANDIDATE_SETS_REF_SCHEMA: &str = "gentle.candidate_sets.ref.v1";
const CANDIDATE_SETS_DISK_INDEX_SCHEMA: &str = "gentle.candidate_sets.disk_index.v1";
const CANDIDATE_SETS_LOAD_WARNING_METADATA_KEY: &str = "candidate_sets_load_warning";
const CANDIDATE_STORE_STRICT_LOAD_ENV: &str = "GENTLE_CANDIDATE_STORE_STRICT_LOAD";
pub const GUIDE_DESIGN_METADATA_KEY: &str = "guide_design";
const GUIDE_DESIGN_SCHEMA: &str = "gentle.guide_design.v1";
pub const WORKFLOW_MACRO_TEMPLATES_METADATA_KEY: &str = "workflow_macro_templates";
const WORKFLOW_MACRO_TEMPLATES_SCHEMA: &str = "gentle.workflow_macro_templates.v1";
pub const CLONING_MACRO_TEMPLATE_SCHEMA: &str = "gentle.cloning_macro_template.v1";
pub const CANDIDATE_MACRO_TEMPLATES_METADATA_KEY: &str = "candidate_macro_templates";
const CANDIDATE_MACRO_TEMPLATES_SCHEMA: &str = "gentle.candidate_macro_templates.v1";
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

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Hash, Default)]
#[serde(rename_all = "snake_case")]
pub enum LinearSequenceLetterLayoutMode {
    #[default]
    ContinuousHelical,
    Condensed10Row,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct DisplaySettings {
    pub show_sequence_panel: bool,
    pub auto_hide_sequence_panel_when_linear_bases_visible: bool,
    pub show_map_panel: bool,
    pub show_features: bool,
    pub show_cds_features: bool,
    pub show_gene_features: bool,
    pub show_mrna_features: bool,
    pub show_tfbs: bool,
    pub regulatory_tracks_near_baseline: bool,
    pub regulatory_feature_max_view_span_bp: usize,
    pub tfbs_display_use_llr_bits: bool,
    pub tfbs_display_min_llr_bits: f64,
    pub tfbs_display_use_llr_quantile: bool,
    pub tfbs_display_min_llr_quantile: f64,
    pub tfbs_display_use_true_log_odds_bits: bool,
    pub tfbs_display_min_true_log_odds_bits: f64,
    pub tfbs_display_use_true_log_odds_quantile: bool,
    pub tfbs_display_min_true_log_odds_quantile: f64,
    pub vcf_display_show_snp: bool,
    pub vcf_display_show_ins: bool,
    pub vcf_display_show_del: bool,
    pub vcf_display_show_sv: bool,
    pub vcf_display_show_other: bool,
    pub vcf_display_pass_only: bool,
    pub vcf_display_use_min_qual: bool,
    pub vcf_display_min_qual: f64,
    pub vcf_display_use_max_qual: bool,
    pub vcf_display_max_qual: f64,
    #[serde(default)]
    pub vcf_display_required_info_keys: Vec<String>,
    pub show_restriction_enzymes: bool,
    pub show_gc_contents: bool,
    #[serde(default = "DisplaySettings::default_gc_content_bin_size_bp")]
    pub gc_content_bin_size_bp: usize,
    pub show_open_reading_frames: bool,
    pub show_methylation_sites: bool,
    pub linear_view_start_bp: usize,
    pub linear_view_span_bp: usize,
    pub linear_sequence_base_text_max_view_span_bp: usize,
    pub linear_sequence_helical_letters_enabled: bool,
    pub linear_sequence_helical_max_view_span_bp: usize,
    #[serde(default = "DisplaySettings::default_linear_sequence_condensed_max_view_span_bp")]
    pub linear_sequence_condensed_max_view_span_bp: usize,
    #[serde(default)]
    pub linear_sequence_letter_layout_mode: LinearSequenceLetterLayoutMode,
    pub linear_sequence_helical_phase_offset_bp: usize,
    pub linear_show_double_strand_bases: bool,
    pub linear_hide_backbone_when_sequence_bases_visible: bool,
    pub linear_reverse_strand_use_upside_down_letters: bool,
    pub feature_details_font_size: f32,
    pub linear_external_feature_label_font_size: f32,
    pub linear_external_feature_label_background_opacity: f32,
}

impl DisplaySettings {
    pub const fn default_gc_content_bin_size_bp() -> usize {
        100
    }

    pub const fn default_linear_sequence_condensed_max_view_span_bp() -> usize {
        1500
    }
}

impl Default for DisplaySettings {
    fn default() -> Self {
        Self {
            show_sequence_panel: true,
            auto_hide_sequence_panel_when_linear_bases_visible: false,
            show_map_panel: true,
            show_features: true,
            show_cds_features: true,
            show_gene_features: true,
            show_mrna_features: true,
            show_tfbs: false,
            regulatory_tracks_near_baseline: false,
            regulatory_feature_max_view_span_bp: 50_000,
            tfbs_display_use_llr_bits: true,
            tfbs_display_min_llr_bits: 0.0,
            tfbs_display_use_llr_quantile: true,
            tfbs_display_min_llr_quantile: 0.95,
            tfbs_display_use_true_log_odds_bits: false,
            tfbs_display_min_true_log_odds_bits: 0.0,
            tfbs_display_use_true_log_odds_quantile: false,
            tfbs_display_min_true_log_odds_quantile: 0.95,
            vcf_display_show_snp: true,
            vcf_display_show_ins: true,
            vcf_display_show_del: true,
            vcf_display_show_sv: true,
            vcf_display_show_other: true,
            vcf_display_pass_only: false,
            vcf_display_use_min_qual: false,
            vcf_display_min_qual: 0.0,
            vcf_display_use_max_qual: false,
            vcf_display_max_qual: 0.0,
            vcf_display_required_info_keys: vec![],
            show_restriction_enzymes: true,
            show_gc_contents: true,
            gc_content_bin_size_bp: Self::default_gc_content_bin_size_bp(),
            show_open_reading_frames: false,
            show_methylation_sites: false,
            linear_view_start_bp: 0,
            linear_view_span_bp: 0,
            linear_sequence_base_text_max_view_span_bp: 500,
            linear_sequence_helical_letters_enabled: false,
            linear_sequence_helical_max_view_span_bp: 2000,
            linear_sequence_condensed_max_view_span_bp:
                Self::default_linear_sequence_condensed_max_view_span_bp(),
            linear_sequence_letter_layout_mode: LinearSequenceLetterLayoutMode::ContinuousHelical,
            linear_sequence_helical_phase_offset_bp: 0,
            linear_show_double_strand_bases: true,
            linear_hide_backbone_when_sequence_bases_visible: false,
            linear_reverse_strand_use_upside_down_letters: true,
            feature_details_font_size: 8.25,
            linear_external_feature_label_font_size: 11.0,
            linear_external_feature_label_background_opacity: 0.9,
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

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum ArrangementMode {
    Serial,
    Plate,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Arrangement {
    pub arrangement_id: String,
    pub mode: ArrangementMode,
    pub name: Option<String>,
    pub lane_container_ids: Vec<ContainerId>,
    #[serde(default)]
    pub ladders: Vec<String>,
    pub created_by_op: Option<OpId>,
    pub created_at_unix_ms: u128,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ContainerState {
    pub containers: HashMap<ContainerId, Container>,
    pub arrangements: HashMap<String, Arrangement>,
    pub seq_to_latest_container: HashMap<SeqId, ContainerId>,
    pub next_container_counter: u64,
    pub next_arrangement_counter: u64,
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

#[derive(Debug)]
enum CandidateSidecarTransaction {
    Noop,
    Replace {
        final_dir: PathBuf,
        staging_dir: PathBuf,
    },
    Remove {
        final_dir: PathBuf,
    },
}

#[derive(Debug)]
enum CandidateSidecarCommitted {
    Noop,
    Replaced {
        final_dir: PathBuf,
        backup_dir: Option<PathBuf>,
    },
    Removed {
        final_dir: PathBuf,
        backup_dir: Option<PathBuf>,
    },
}

impl ProjectState {
    pub fn load_from_path(path: &str) -> Result<Self, EngineError> {
        let text = std::fs::read_to_string(path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not read state file '{path}': {e}"),
        })?;
        let mut state: Self = serde_json::from_str(&text).map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!("Could not parse state JSON '{path}': {e}"),
        })?;
        if let Err(err) = state.hydrate_candidate_store_from_external_ref(path) {
            if Self::strict_candidate_store_load_enabled() {
                return Err(err);
            }
            let warning = err.message;
            state.metadata.remove(CANDIDATE_SETS_METADATA_KEY);
            state.metadata.insert(
                CANDIDATE_SETS_LOAD_WARNING_METADATA_KEY.to_string(),
                json!({
                    "message": warning,
                    "strict_env": CANDIDATE_STORE_STRICT_LOAD_ENV,
                }),
            );
        } else {
            state
                .metadata
                .remove(CANDIDATE_SETS_LOAD_WARNING_METADATA_KEY);
        }
        Ok(state)
    }

    pub fn save_to_path(&self, path: &str) -> Result<(), EngineError> {
        let mut state_for_disk = self.clone();
        let project_path = Path::new(path);
        let sidecar_tx =
            state_for_disk.prepare_candidate_store_sidecar_transaction(project_path)?;
        let text = serde_json::to_string_pretty(&state_for_disk).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize state: {e}"),
        })?;
        let committed = Self::commit_candidate_store_transaction(sidecar_tx)?;
        if let Err(write_err) = Self::write_text_file_atomically(project_path, &text) {
            if let Err(rollback_err) = Self::rollback_candidate_store_transaction(committed) {
                return Err(EngineError {
                    code: ErrorCode::Io,
                    message: format!(
                        "{}; candidate-sidecar rollback also failed: {}",
                        write_err.message, rollback_err.message
                    ),
                });
            }
            return Err(write_err);
        }
        Self::finalize_candidate_store_transaction(committed)
    }

    fn candidate_store_sidecar_basename(project_path: &Path) -> String {
        let fallback = "project.gentle.json".to_string();
        let base = project_path
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or(&fallback);
        format!("{base}.candidates")
    }

    fn candidate_store_sidecar_dir(project_path: &Path) -> PathBuf {
        let parent = project_path
            .parent()
            .map(Path::to_path_buf)
            .unwrap_or_else(|| PathBuf::from("."));
        parent.join(Self::candidate_store_sidecar_basename(project_path))
    }

    fn candidate_store_sidecar_index_rel(project_path: &Path) -> PathBuf {
        PathBuf::from(Self::candidate_store_sidecar_basename(project_path)).join("index.json")
    }

    fn strict_candidate_store_load_enabled() -> bool {
        matches!(
            env::var(CANDIDATE_STORE_STRICT_LOAD_ENV)
                .unwrap_or_default()
                .trim()
                .to_ascii_lowercase()
                .as_str(),
            "1" | "true" | "yes" | "on"
        )
    }

    fn write_text_file_atomically(path: &Path, text: &str) -> Result<(), EngineError> {
        let parent = path
            .parent()
            .map(Path::to_path_buf)
            .unwrap_or_else(|| PathBuf::from("."));
        std::fs::create_dir_all(&parent).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not create parent directory '{}' for state save: {e}",
                parent.display()
            ),
        })?;
        let mut tmp = NamedTempFile::new_in(&parent).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not create temporary state file in '{}': {e}",
                parent.display()
            ),
        })?;
        tmp.write_all(text.as_bytes()).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not write temporary state file for '{}': {e}",
                path.display()
            ),
        })?;
        tmp.flush().map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not flush temporary state file for '{}': {e}",
                path.display()
            ),
        })?;
        tmp.as_file().sync_all().map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not sync temporary state file for '{}': {e}",
                path.display()
            ),
        })?;
        tmp.persist(path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not replace state file '{}': {}",
                path.display(),
                e.error
            ),
        })?;
        Ok(())
    }

    fn candidate_sidecar_staging_dir(final_dir: &Path) -> Result<PathBuf, EngineError> {
        let parent = final_dir
            .parent()
            .map(Path::to_path_buf)
            .unwrap_or_else(|| PathBuf::from("."));
        let base = final_dir
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("candidate_sidecar");
        let nonce = format!("{}-{}", std::process::id(), GentleEngine::now_unix_ms());
        for attempt in 0..64 {
            let suffix = if attempt == 0 {
                nonce.clone()
            } else {
                format!("{nonce}-{attempt}")
            };
            let candidate = parent.join(format!(".{base}.staging.{suffix}"));
            match std::fs::create_dir(&candidate) {
                Ok(_) => return Ok(candidate),
                Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(e) => {
                    return Err(EngineError {
                        code: ErrorCode::Io,
                        message: format!(
                            "Could not create candidate-sidecar staging directory '{}': {e}",
                            candidate.display()
                        ),
                    });
                }
            }
        }
        Err(EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not create unique candidate-sidecar staging directory under '{}'",
                parent.display()
            ),
        })
    }

    fn candidate_sidecar_backup_dir(final_dir: &Path) -> Result<PathBuf, EngineError> {
        let parent = final_dir
            .parent()
            .map(Path::to_path_buf)
            .unwrap_or_else(|| PathBuf::from("."));
        let base = final_dir
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("candidate_sidecar");
        let nonce = format!("{}-{}", std::process::id(), GentleEngine::now_unix_ms());
        for attempt in 0..64 {
            let suffix = if attempt == 0 {
                nonce.clone()
            } else {
                format!("{nonce}-{attempt}")
            };
            let candidate = parent.join(format!(".{base}.backup.{suffix}"));
            if !candidate.exists() {
                return Ok(candidate);
            }
        }
        Err(EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not allocate candidate-sidecar backup path under '{}'",
                parent.display()
            ),
        })
    }

    fn sanitize_candidate_set_file_stem(name: &str) -> String {
        let mut out = String::with_capacity(name.len());
        for ch in name.chars() {
            if ch.is_ascii_alphanumeric() {
                out.push(ch.to_ascii_lowercase());
            } else if ch.is_whitespace() || matches!(ch, '-' | '_' | '.') {
                out.push('_');
            }
        }
        let trimmed = out.trim_matches('_');
        if trimmed.is_empty() {
            "set".to_string()
        } else {
            trimmed.chars().take(48).collect()
        }
    }

    fn parse_inline_candidate_store_value(value: &serde_json::Value) -> Option<CandidateStore> {
        let schema = value
            .get("schema")
            .and_then(|v| v.as_str())
            .unwrap_or_default();
        let has_sets = value.get("sets").map(|v| v.is_object()).unwrap_or(false);
        if schema == CANDIDATE_SETS_REF_SCHEMA {
            return None;
        }
        if !has_sets && schema != CANDIDATE_SETS_SCHEMA {
            return None;
        }
        let mut store: CandidateStore = serde_json::from_value(value.clone()).ok()?;
        if store.schema.trim().is_empty() {
            store.schema = CANDIDATE_SETS_SCHEMA.to_string();
        }
        Some(store)
    }

    fn candidate_set_metrics_for_disk(set: &CandidateSet) -> Vec<String> {
        let mut names = BTreeSet::new();
        for candidate in &set.candidates {
            for metric in candidate.metrics.keys() {
                names.insert(metric.to_string());
            }
        }
        names.into_iter().collect()
    }

    fn load_candidate_store_from_ref(
        project_path: &Path,
        reference: &CandidateStoreReference,
    ) -> Result<CandidateStore, EngineError> {
        let index_path = if Path::new(&reference.index_path).is_absolute() {
            PathBuf::from(&reference.index_path)
        } else {
            project_path
                .parent()
                .map(Path::to_path_buf)
                .unwrap_or_else(|| PathBuf::from("."))
                .join(&reference.index_path)
        };
        let index_text = std::fs::read_to_string(&index_path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not read candidate-store index '{}': {e}",
                index_path.display()
            ),
        })?;
        let index: CandidateStoreDiskIndex =
            serde_json::from_str(&index_text).map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not parse candidate-store index '{}': {e}",
                    index_path.display()
                ),
            })?;
        if !index.schema.trim().is_empty() && index.schema != CANDIDATE_SETS_DISK_INDEX_SCHEMA {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Unsupported candidate-store index schema '{}' in '{}'",
                    index.schema,
                    index_path.display()
                ),
            });
        }
        let index_dir = index_path
            .parent()
            .map(Path::to_path_buf)
            .unwrap_or_else(|| PathBuf::from("."));
        let mut sets: HashMap<String, CandidateSet> = HashMap::new();
        for entry in index.sets {
            let records_path = if Path::new(&entry.records_path).is_absolute() {
                PathBuf::from(&entry.records_path)
            } else {
                index_dir.join(&entry.records_path)
            };
            let file = File::open(&records_path).map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not read candidate records '{}' for set '{}': {e}",
                    records_path.display(),
                    entry.name
                ),
            })?;
            let reader = BufReader::new(file);
            let mut candidates: Vec<CandidateRecord> = vec![];
            for (line_no, line_result) in reader.lines().enumerate() {
                let line = line_result.map_err(|e| EngineError {
                    code: ErrorCode::Io,
                    message: format!(
                        "Could not read candidate record line {} from '{}': {e}",
                        line_no + 1,
                        records_path.display()
                    ),
                })?;
                let trimmed = line.trim();
                if trimmed.is_empty() {
                    continue;
                }
                let candidate =
                    serde_json::from_str::<CandidateRecord>(trimmed).map_err(|e| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Invalid candidate record JSON at '{}':{}: {}",
                            records_path.display(),
                            line_no + 1,
                            e
                        ),
                    })?;
                candidates.push(candidate);
            }
            sets.insert(
                entry.name.clone(),
                CandidateSet {
                    name: entry.name,
                    created_at_unix_ms: entry.created_at_unix_ms,
                    source_seq_ids: entry.source_seq_ids,
                    candidates,
                },
            );
        }
        Ok(CandidateStore {
            schema: CANDIDATE_SETS_SCHEMA.to_string(),
            updated_at_unix_ms: index.updated_at_unix_ms,
            sets,
        })
    }

    fn hydrate_candidate_store_from_external_ref(
        &mut self,
        project_path: &str,
    ) -> Result<(), EngineError> {
        let Some(value) = self.metadata.get(CANDIDATE_SETS_METADATA_KEY).cloned() else {
            return Ok(());
        };
        if Self::parse_inline_candidate_store_value(&value).is_some() {
            return Ok(());
        }
        let reference: CandidateStoreReference =
            serde_json::from_value(value.clone()).map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not parse candidate-store reference metadata '{}': {e}",
                    CANDIDATE_SETS_METADATA_KEY
                ),
            })?;
        if reference.schema != CANDIDATE_SETS_REF_SCHEMA {
            return Ok(());
        }
        let store = Self::load_candidate_store_from_ref(Path::new(project_path), &reference)?;
        let store_value = serde_json::to_value(store).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not hydrate candidate-store metadata: {e}"),
        })?;
        self.metadata
            .insert(CANDIDATE_SETS_METADATA_KEY.to_string(), store_value);
        Ok(())
    }

    fn load_candidate_store_for_externalization(
        &self,
        project_path: &Path,
        value: serde_json::Value,
    ) -> Result<Option<CandidateStore>, EngineError> {
        if let Some(store) = Self::parse_inline_candidate_store_value(&value) {
            return Ok(Some(store));
        }
        let reference: CandidateStoreReference =
            serde_json::from_value(value).map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not parse candidate-store reference metadata '{}': {e}",
                    CANDIDATE_SETS_METADATA_KEY
                ),
            })?;
        if reference.schema != CANDIDATE_SETS_REF_SCHEMA {
            return Ok(None);
        }
        let store = Self::load_candidate_store_from_ref(project_path, &reference)?;
        Ok(Some(store))
    }

    fn prepare_candidate_store_sidecar_transaction(
        &mut self,
        project_path: &Path,
    ) -> Result<CandidateSidecarTransaction, EngineError> {
        let sidecar_dir = Self::candidate_store_sidecar_dir(project_path);
        let Some(value) = self.metadata.get(CANDIDATE_SETS_METADATA_KEY).cloned() else {
            if sidecar_dir.exists() {
                return Ok(CandidateSidecarTransaction::Remove {
                    final_dir: sidecar_dir,
                });
            }
            return Ok(CandidateSidecarTransaction::Noop);
        };
        let Some(mut store) = self.load_candidate_store_for_externalization(project_path, value)?
        else {
            return Ok(CandidateSidecarTransaction::Noop);
        };
        if store.sets.is_empty() {
            self.metadata.remove(CANDIDATE_SETS_METADATA_KEY);
            if sidecar_dir.exists() {
                return Ok(CandidateSidecarTransaction::Remove {
                    final_dir: sidecar_dir,
                });
            }
            return Ok(CandidateSidecarTransaction::Noop);
        }
        if store.schema.trim().is_empty() {
            store.schema = CANDIDATE_SETS_SCHEMA.to_string();
        }
        let staging_dir = Self::candidate_sidecar_staging_dir(&sidecar_dir)?;

        let mut set_names: Vec<String> = store.sets.keys().cloned().collect();
        set_names.sort_unstable();
        let mut index_entries: Vec<CandidateStoreDiskSetIndexEntry> = vec![];
        for (idx, set_name) in set_names.iter().enumerate() {
            let Some(set) = store.sets.get(set_name) else {
                continue;
            };
            let mut hasher = DefaultHasher::new();
            set_name.hash(&mut hasher);
            let stem = Self::sanitize_candidate_set_file_stem(set_name);
            let records_filename = format!("{:03}_{}_{}.jsonl", idx + 1, stem, hasher.finish());
            let records_path = staging_dir.join(&records_filename);
            let records_file = File::create(&records_path).map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not create candidate-set records file '{}': {e}",
                    records_path.display()
                ),
            })?;
            let mut writer = BufWriter::new(records_file);
            for candidate in &set.candidates {
                serde_json::to_writer(&mut writer, candidate).map_err(|e| EngineError {
                    code: ErrorCode::Internal,
                    message: format!(
                        "Could not serialize candidate record for set '{}': {e}",
                        set_name
                    ),
                })?;
                writer.write_all(b"\n").map_err(|e| EngineError {
                    code: ErrorCode::Io,
                    message: format!(
                        "Could not write candidate record for set '{}': {e}",
                        set_name
                    ),
                })?;
            }
            writer.flush().map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not flush candidate-set records file '{}': {e}",
                    records_path.display()
                ),
            })?;

            index_entries.push(CandidateStoreDiskSetIndexEntry {
                name: set.name.clone(),
                created_at_unix_ms: set.created_at_unix_ms,
                source_seq_ids: set.source_seq_ids.clone(),
                candidate_count: set.candidates.len(),
                metrics: Self::candidate_set_metrics_for_disk(set),
                records_path: records_filename,
            });
        }

        let index = CandidateStoreDiskIndex {
            schema: CANDIDATE_SETS_DISK_INDEX_SCHEMA.to_string(),
            updated_at_unix_ms: store.updated_at_unix_ms,
            set_count: index_entries.len(),
            sets: index_entries,
        };
        let index_path = staging_dir.join("index.json");
        let index_text = serde_json::to_string_pretty(&index).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize candidate-store index: {e}"),
        })?;
        std::fs::write(&index_path, index_text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not write candidate-store index '{}': {e}",
                index_path.display()
            ),
        })?;

        let reference = CandidateStoreReference {
            schema: CANDIDATE_SETS_REF_SCHEMA.to_string(),
            storage: "jsonl_indexed".to_string(),
            index_path: Self::candidate_store_sidecar_index_rel(project_path)
                .to_string_lossy()
                .to_string(),
            set_count: index.set_count,
            updated_at_unix_ms: index.updated_at_unix_ms,
        };
        let reference_value = serde_json::to_value(reference).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize candidate-store reference metadata: {e}"),
        })?;
        self.metadata
            .insert(CANDIDATE_SETS_METADATA_KEY.to_string(), reference_value);
        Ok(CandidateSidecarTransaction::Replace {
            final_dir: sidecar_dir,
            staging_dir,
        })
    }

    fn commit_candidate_store_transaction(
        tx: CandidateSidecarTransaction,
    ) -> Result<CandidateSidecarCommitted, EngineError> {
        match tx {
            CandidateSidecarTransaction::Noop => Ok(CandidateSidecarCommitted::Noop),
            CandidateSidecarTransaction::Replace {
                final_dir,
                staging_dir,
            } => {
                let mut backup_dir: Option<PathBuf> = None;
                if final_dir.exists() {
                    let backup = Self::candidate_sidecar_backup_dir(&final_dir)?;
                    std::fs::rename(&final_dir, &backup).map_err(|e| EngineError {
                        code: ErrorCode::Io,
                        message: format!(
                            "Could not move existing candidate-sidecar directory '{}' to backup '{}': {e}",
                            final_dir.display(),
                            backup.display()
                        ),
                    })?;
                    backup_dir = Some(backup);
                }
                if let Err(e) = std::fs::rename(&staging_dir, &final_dir) {
                    if let Some(backup) = &backup_dir {
                        let _ = std::fs::rename(backup, &final_dir);
                    }
                    let _ = std::fs::remove_dir_all(&staging_dir);
                    return Err(EngineError {
                        code: ErrorCode::Io,
                        message: format!(
                            "Could not activate candidate-sidecar directory '{}' from staging '{}': {e}",
                            final_dir.display(),
                            staging_dir.display()
                        ),
                    });
                }
                Ok(CandidateSidecarCommitted::Replaced {
                    final_dir,
                    backup_dir,
                })
            }
            CandidateSidecarTransaction::Remove { final_dir } => {
                if !final_dir.exists() {
                    return Ok(CandidateSidecarCommitted::Noop);
                }
                let backup = Self::candidate_sidecar_backup_dir(&final_dir)?;
                std::fs::rename(&final_dir, &backup).map_err(|e| EngineError {
                    code: ErrorCode::Io,
                    message: format!(
                        "Could not stage candidate-sidecar removal '{}' -> '{}': {e}",
                        final_dir.display(),
                        backup.display()
                    ),
                })?;
                Ok(CandidateSidecarCommitted::Removed {
                    final_dir,
                    backup_dir: Some(backup),
                })
            }
        }
    }

    fn rollback_candidate_store_transaction(
        committed: CandidateSidecarCommitted,
    ) -> Result<(), EngineError> {
        match committed {
            CandidateSidecarCommitted::Noop => Ok(()),
            CandidateSidecarCommitted::Replaced {
                final_dir,
                backup_dir,
            } => {
                if let Some(backup) = backup_dir {
                    if final_dir.exists() {
                        std::fs::remove_dir_all(&final_dir).map_err(|e| EngineError {
                            code: ErrorCode::Io,
                            message: format!(
                                "Could not remove partially committed candidate-sidecar directory '{}': {e}",
                                final_dir.display()
                            ),
                        })?;
                    }
                    std::fs::rename(&backup, &final_dir).map_err(|e| EngineError {
                        code: ErrorCode::Io,
                        message: format!(
                            "Could not restore candidate-sidecar backup '{}' -> '{}': {e}",
                            backup.display(),
                            final_dir.display()
                        ),
                    })?;
                } else if final_dir.exists() {
                    std::fs::remove_dir_all(&final_dir).map_err(|e| EngineError {
                        code: ErrorCode::Io,
                        message: format!(
                            "Could not remove new candidate-sidecar directory '{}' during rollback: {e}",
                            final_dir.display()
                        ),
                    })?;
                }
                Ok(())
            }
            CandidateSidecarCommitted::Removed {
                final_dir,
                backup_dir,
            } => {
                if let Some(backup) = backup_dir {
                    std::fs::rename(&backup, &final_dir).map_err(|e| EngineError {
                        code: ErrorCode::Io,
                        message: format!(
                            "Could not restore candidate-sidecar directory '{}' -> '{}': {e}",
                            backup.display(),
                            final_dir.display()
                        ),
                    })?;
                }
                Ok(())
            }
        }
    }

    fn finalize_candidate_store_transaction(
        committed: CandidateSidecarCommitted,
    ) -> Result<(), EngineError> {
        match committed {
            CandidateSidecarCommitted::Noop => Ok(()),
            CandidateSidecarCommitted::Replaced { backup_dir, .. }
            | CandidateSidecarCommitted::Removed { backup_dir, .. } => {
                if let Some(backup) = backup_dir {
                    std::fs::remove_dir_all(&backup).map_err(|e| EngineError {
                        code: ErrorCode::Io,
                        message: format!(
                            "Could not remove candidate-sidecar backup directory '{}': {e}",
                            backup.display()
                        ),
                    })?;
                }
                Ok(())
            }
        }
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

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub enum AnchorBoundary {
    Start,
    End,
    Middle,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum AnchorDirection {
    Upstream,
    Downstream,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub enum SequenceAnchor {
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

// Backward-compatible alias kept while adapter/docs migrate to SequenceAnchor.
pub type AnchoredRegionAnchor = SequenceAnchor;

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
struct CandidateStoreReference {
    schema: String,
    storage: String,
    index_path: String,
    set_count: usize,
    updated_at_unix_ms: u128,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
struct CandidateStoreDiskIndex {
    schema: String,
    updated_at_unix_ms: u128,
    set_count: usize,
    sets: Vec<CandidateStoreDiskSetIndexEntry>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
struct CandidateStoreDiskSetIndexEntry {
    name: String,
    created_at_unix_ms: u128,
    source_seq_ids: Vec<String>,
    candidate_count: usize,
    metrics: Vec<String>,
    records_path: String,
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

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum CandidateFeatureGeometryMode {
    #[default]
    FeatureSpan,
    FeatureParts,
    FeatureBoundaries,
}

impl CandidateFeatureGeometryMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::FeatureSpan => "feature_span",
            Self::FeatureParts => "feature_parts",
            Self::FeatureBoundaries => "feature_boundaries",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum CandidateFeatureBoundaryMode {
    #[default]
    Any,
    FivePrime,
    ThreePrime,
    Start,
    End,
}

impl CandidateFeatureBoundaryMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Any => "any",
            Self::FivePrime => "five_prime",
            Self::ThreePrime => "three_prime",
            Self::Start => "start",
            Self::End => "end",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum CandidateFeatureStrandRelation {
    #[default]
    Any,
    Same,
    Opposite,
}

impl CandidateFeatureStrandRelation {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Any => "any",
            Self::Same => "same",
            Self::Opposite => "opposite",
        }
    }
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

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum CandidateObjectiveDirection {
    #[default]
    Maximize,
    Minimize,
}

impl CandidateObjectiveDirection {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Maximize => "maximize",
            Self::Minimize => "minimize",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CandidateObjectiveSpec {
    pub metric: String,
    #[serde(default)]
    pub direction: CandidateObjectiveDirection,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CandidateWeightedObjectiveTerm {
    pub metric: String,
    pub weight: f64,
    #[serde(default)]
    pub direction: CandidateObjectiveDirection,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum CandidateTieBreakPolicy {
    #[default]
    SeqStartEnd,
    SeqEndStart,
    LengthAscending,
    LengthDescending,
    SequenceLexicographic,
}

impl CandidateTieBreakPolicy {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::SeqStartEnd => "seq_start_end",
            Self::SeqEndStart => "seq_end_start",
            Self::LengthAscending => "length_ascending",
            Self::LengthDescending => "length_descending",
            Self::SequenceLexicographic => "sequence_lexicographic",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum GuideU6TerminatorWindow {
    SpacerOnly,
    #[default]
    SpacerPlusTail,
}

impl GuideU6TerminatorWindow {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::SpacerOnly => "spacer_only",
            Self::SpacerPlusTail => "spacer_plus_tail",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GuideOligoExportFormat {
    CsvTable,
    PlateCsv,
    Fasta,
}

impl GuideOligoExportFormat {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::CsvTable => "csv_table",
            Self::PlateCsv => "plate_csv",
            Self::Fasta => "fasta",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum GuideOligoPlateFormat {
    #[default]
    Plate96,
    Plate384,
}

impl GuideOligoPlateFormat {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Plate96 => "96",
            Self::Plate384 => "384",
        }
    }

    fn dimensions(self) -> (usize, usize) {
        match self {
            Self::Plate96 => (8, 12),
            Self::Plate384 => (16, 24),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
struct GuideDesignStore {
    schema: String,
    updated_at_unix_ms: u128,
    guide_sets: HashMap<String, GuideSet>,
    practical_filter_reports: HashMap<String, GuidePracticalFilterReport>,
    oligo_sets: HashMap<String, GuideOligoSet>,
    latest_oligo_set_by_guide_set: HashMap<String, String>,
    audit_log: Vec<GuideDesignAuditEntry>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GuideSet {
    pub guide_set_id: String,
    pub guides: Vec<GuideCandidate>,
    pub created_at_unix_ms: u128,
    pub updated_at_unix_ms: u128,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GuideCandidate {
    pub guide_id: String,
    pub seq_id: String,
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub strand: String,
    pub protospacer: String,
    pub pam: String,
    pub nuclease: String,
    pub cut_offset_from_protospacer_start: usize,
    pub rank: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GuideSetSummary {
    pub guide_set_id: String,
    pub guide_count: usize,
    pub created_at_unix_ms: u128,
    pub updated_at_unix_ms: u128,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct GuidePracticalFilterConfig {
    pub gc_min: Option<f64>,
    pub gc_max: Option<f64>,
    pub max_homopolymer_run: Option<usize>,
    #[serde(default)]
    pub max_homopolymer_run_per_base: HashMap<String, usize>,
    pub reject_ambiguous_bases: bool,
    pub avoid_u6_terminator_tttt: bool,
    pub u6_terminator_window: GuideU6TerminatorWindow,
    pub max_dinucleotide_repeat_units: Option<usize>,
    #[serde(default)]
    pub forbidden_motifs: Vec<String>,
    pub required_5prime_base: Option<String>,
    pub allow_5prime_g_extension: bool,
}

impl Default for GuidePracticalFilterConfig {
    fn default() -> Self {
        Self {
            gc_min: None,
            gc_max: None,
            max_homopolymer_run: None,
            max_homopolymer_run_per_base: HashMap::new(),
            reject_ambiguous_bases: true,
            avoid_u6_terminator_tttt: true,
            u6_terminator_window: GuideU6TerminatorWindow::SpacerPlusTail,
            max_dinucleotide_repeat_units: None,
            forbidden_motifs: vec![],
            required_5prime_base: None,
            allow_5prime_g_extension: false,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GuideFilterReason {
    pub code: String,
    pub message: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GuidePracticalFilterResult {
    pub guide_id: String,
    pub passed: bool,
    #[serde(default)]
    pub reasons: Vec<GuideFilterReason>,
    #[serde(default)]
    pub warnings: Vec<String>,
    #[serde(default)]
    pub metrics: HashMap<String, f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GuidePracticalFilterReport {
    pub guide_set_id: String,
    pub generated_at_unix_ms: u128,
    pub config: GuidePracticalFilterConfig,
    pub passed_count: usize,
    pub rejected_count: usize,
    #[serde(default)]
    pub results: Vec<GuidePracticalFilterResult>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GuideOligoTemplateSpec {
    pub template_id: String,
    pub description: String,
    pub forward_prefix: String,
    pub forward_suffix: String,
    pub reverse_prefix: String,
    pub reverse_suffix: String,
    pub reverse_uses_reverse_complement_of_spacer: bool,
    pub uppercase_output: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GuideOligoRecord {
    pub guide_id: String,
    pub rank: Option<usize>,
    pub forward_oligo: String,
    pub reverse_oligo: String,
    #[serde(default)]
    pub notes: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GuideOligoSet {
    pub oligo_set_id: String,
    pub guide_set_id: String,
    pub generated_at_unix_ms: u128,
    pub template: GuideOligoTemplateSpec,
    pub apply_5prime_g_extension: bool,
    #[serde(default)]
    pub records: Vec<GuideOligoRecord>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GuideOligoExportReport {
    pub guide_set_id: String,
    pub oligo_set_id: String,
    pub format: String,
    pub path: String,
    pub exported_records: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GuideProtocolExportReport {
    pub guide_set_id: String,
    pub oligo_set_id: String,
    pub path: String,
    pub guide_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
struct GuideDesignAuditEntry {
    unix_ms: u128,
    operation: String,
    guide_set_id: String,
    payload: serde_json::Value,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
struct WorkflowMacroTemplateStore {
    schema: String,
    updated_at_unix_ms: u128,
    templates: HashMap<String, WorkflowMacroTemplate>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct WorkflowMacroTemplate {
    pub name: String,
    pub description: Option<String>,
    pub details_url: Option<String>,
    pub parameters: Vec<WorkflowMacroTemplateParam>,
    #[serde(default = "default_cloning_macro_template_schema")]
    pub template_schema: String,
    pub script: String,
    pub created_at_unix_ms: u128,
    pub updated_at_unix_ms: u128,
}

fn default_cloning_macro_template_schema() -> String {
    CLONING_MACRO_TEMPLATE_SCHEMA.to_string()
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct WorkflowMacroTemplateParam {
    pub name: String,
    pub default_value: Option<String>,
    pub required: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WorkflowMacroTemplateSummary {
    pub name: String,
    pub description: Option<String>,
    pub details_url: Option<String>,
    pub parameter_count: usize,
    pub created_at_unix_ms: u128,
    pub updated_at_unix_ms: u128,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
struct CandidateMacroTemplateStore {
    schema: String,
    updated_at_unix_ms: u128,
    templates: HashMap<String, CandidateMacroTemplate>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct CandidateMacroTemplate {
    pub name: String,
    pub description: Option<String>,
    pub details_url: Option<String>,
    pub parameters: Vec<CandidateMacroTemplateParam>,
    pub script: String,
    pub created_at_unix_ms: u128,
    pub updated_at_unix_ms: u128,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct CandidateMacroTemplateParam {
    pub name: String,
    pub default_value: Option<String>,
    pub required: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CandidateMacroTemplateSummary {
    pub name: String,
    pub description: Option<String>,
    pub details_url: Option<String>,
    pub parameter_count: usize,
    pub created_at_unix_ms: u128,
    pub updated_at_unix_ms: u128,
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

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GenomeAnchorSide {
    FivePrime,
    ThreePrime,
}

impl GenomeAnchorSide {
    fn as_str(self) -> &'static str {
        match self {
            Self::FivePrime => "5prime",
            Self::ThreePrime => "3prime",
        }
    }
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
    RenderFeatureExpertSvg {
        seq_id: SeqId,
        target: FeatureExpertTarget,
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
        #[serde(default)]
        container_ids: Option<Vec<ContainerId>>,
        #[serde(default)]
        arrangement_id: Option<String>,
    },
    CreateArrangementSerial {
        container_ids: Vec<ContainerId>,
        arrangement_id: Option<String>,
        name: Option<String>,
        #[serde(default)]
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
        #[serde(default)]
        timeout_seconds: Option<u64>,
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
    ExtendGenomeAnchor {
        seq_id: SeqId,
        side: GenomeAnchorSide,
        length_bp: usize,
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
        anchor: SequenceAnchor,
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
    GenerateCandidateSetBetweenAnchors {
        set_name: String,
        seq_id: SeqId,
        anchor_a: SequenceAnchor,
        anchor_b: SequenceAnchor,
        length_bp: usize,
        step_bp: usize,
        limit: Option<usize>,
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
        #[serde(default)]
        feature_geometry_mode: Option<CandidateFeatureGeometryMode>,
        #[serde(default)]
        feature_boundary_mode: Option<CandidateFeatureBoundaryMode>,
        #[serde(default)]
        feature_strand_relation: Option<CandidateFeatureStrandRelation>,
        limit: Option<usize>,
    },
    DeleteCandidateSet {
        set_name: String,
    },
    UpsertGuideSet {
        guide_set_id: String,
        guides: Vec<GuideCandidate>,
    },
    DeleteGuideSet {
        guide_set_id: String,
    },
    FilterGuidesPractical {
        guide_set_id: String,
        #[serde(default)]
        config: GuidePracticalFilterConfig,
        #[serde(default)]
        output_guide_set_id: Option<String>,
    },
    GenerateGuideOligos {
        guide_set_id: String,
        template_id: String,
        #[serde(default)]
        apply_5prime_g_extension: Option<bool>,
        #[serde(default)]
        output_oligo_set_id: Option<String>,
        #[serde(default)]
        passed_only: Option<bool>,
    },
    ExportGuideOligos {
        guide_set_id: String,
        #[serde(default)]
        oligo_set_id: Option<String>,
        format: GuideOligoExportFormat,
        path: String,
        #[serde(default)]
        plate_format: Option<GuideOligoPlateFormat>,
    },
    ExportGuideProtocolText {
        guide_set_id: String,
        #[serde(default)]
        oligo_set_id: Option<String>,
        path: String,
        #[serde(default)]
        include_qc_checklist: Option<bool>,
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
        #[serde(default)]
        feature_geometry_mode: Option<CandidateFeatureGeometryMode>,
        #[serde(default)]
        feature_boundary_mode: Option<CandidateFeatureBoundaryMode>,
        #[serde(default)]
        feature_strand_relation: Option<CandidateFeatureStrandRelation>,
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
    ScoreCandidateSetWeightedObjective {
        set_name: String,
        metric: String,
        #[serde(default)]
        objectives: Vec<CandidateWeightedObjectiveTerm>,
        #[serde(default)]
        normalize_metrics: Option<bool>,
    },
    TopKCandidateSet {
        input_set: String,
        output_set: String,
        metric: String,
        k: usize,
        #[serde(default)]
        direction: Option<CandidateObjectiveDirection>,
        #[serde(default)]
        tie_break: Option<CandidateTieBreakPolicy>,
    },
    ParetoFrontierCandidateSet {
        input_set: String,
        output_set: String,
        #[serde(default)]
        objectives: Vec<CandidateObjectiveSpec>,
        #[serde(default)]
        max_candidates: Option<usize>,
        #[serde(default)]
        tie_break: Option<CandidateTieBreakPolicy>,
    },
    UpsertWorkflowMacroTemplate {
        name: String,
        description: Option<String>,
        #[serde(default)]
        details_url: Option<String>,
        #[serde(default)]
        parameters: Vec<WorkflowMacroTemplateParam>,
        script: String,
    },
    DeleteWorkflowMacroTemplate {
        name: String,
    },
    UpsertCandidateMacroTemplate {
        name: String,
        description: Option<String>,
        #[serde(default)]
        details_url: Option<String>,
        #[serde(default)]
        parameters: Vec<CandidateMacroTemplateParam>,
        script: String,
    },
    DeleteCandidateMacroTemplate {
        name: String,
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

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SequenceGenomeAnchorSummary {
    pub seq_id: String,
    pub genome_id: String,
    pub chromosome: String,
    pub start_1based: usize,
    pub end_1based: usize,
    pub strand: Option<char>,
}

#[derive(Debug, Clone)]
struct GenomeSequenceAnchor {
    genome_id: String,
    chromosome: String,
    start_1based: usize,
    end_1based: usize,
    strand: Option<char>,
    catalog_path: Option<String>,
    cache_dir: Option<String>,
}

#[derive(Debug, Clone, Default)]
struct EngineHistoryCheckpoint {
    state: ProjectState,
    journal: Vec<OperationRecord>,
    op_counter: u64,
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
    format: Option<String>,
    sample_columns: Vec<String>,
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
    cancelled: bool,
    warnings: Vec<String>,
    skipped_wrong_chromosome_examples: Vec<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum VcfVariantClass {
    Snp,
    Ins,
    Del,
    Sv,
    Other,
}

impl VcfVariantClass {
    fn as_str(self) -> &'static str {
        match self {
            Self::Snp => "SNP",
            Self::Ins => "INS",
            Self::Del => "DEL",
            Self::Sv => "SV",
            Self::Other => "OTHER",
        }
    }
}

#[derive(Debug, Clone, Default)]
struct VcfAltGenotypeSummary {
    carriers: usize,
    het: usize,
    hom_alt: usize,
    mixed_alt: usize,
    haploid_alt: usize,
    phased_carriers: usize,
    unphased_carriers: usize,
    carrier_samples: Vec<String>,
}

#[derive(Debug, Clone)]
struct FeatureDistanceTarget {
    feature_index: usize,
    kind_upper: String,
    labels_upper: Vec<String>,
    start_0based: usize,
    end_0based: usize,
    strand: Option<char>,
}

#[derive(Debug, Clone, Copy)]
struct FeatureLocationSegment {
    start_0based: usize,
    end_0based: usize,
    strand: Option<char>,
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
pub struct GenomeTrackImportProgress {
    pub seq_id: String,
    pub source: String,
    pub path: String,
    pub parsed_records: usize,
    pub imported_features: usize,
    pub skipped_records: usize,
    pub done: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum OperationProgress {
    Tfbs(TfbsProgress),
    GenomePrepare(PrepareGenomeProgress),
    GenomeTrackImport(GenomeTrackImportProgress),
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
pub struct EngineArrangementSummary {
    pub id: String,
    pub mode: String,
    pub lane_count: usize,
    pub lane_container_ids: Vec<String>,
    pub ladders: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EngineStateSummary {
    pub sequence_count: usize,
    pub sequences: Vec<EngineSequenceSummary>,
    pub container_count: usize,
    pub containers: Vec<EngineContainerSummary>,
    pub arrangement_count: usize,
    pub arrangements: Vec<EngineArrangementSummary>,
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
    #[serde(skip, default)]
    undo_stack: Vec<EngineHistoryCheckpoint>,
    #[serde(skip, default)]
    redo_stack: Vec<EngineHistoryCheckpoint>,
    #[serde(skip, default = "GentleEngine::default_history_limit")]
    history_limit: usize,
}

impl GentleEngine {
    fn default_history_limit() -> usize {
        256
    }

    pub fn new() -> Self {
        let mut ret = Self::default();
        if ret.history_limit == 0 {
            ret.history_limit = Self::default_history_limit();
        }
        ret
    }

    pub fn from_state(state: ProjectState) -> Self {
        let mut ret = Self {
            state,
            ..Self::default()
        };
        if ret.history_limit == 0 {
            ret.history_limit = Self::default_history_limit();
        }
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

    pub fn sequence_genome_anchor_summary(
        &self,
        seq_id: &str,
    ) -> Result<SequenceGenomeAnchorSummary, EngineError> {
        let anchor = self.latest_genome_anchor_for_seq(seq_id)?;
        Ok(SequenceGenomeAnchorSummary {
            seq_id: seq_id.to_string(),
            genome_id: anchor.genome_id,
            chromosome: anchor.chromosome,
            start_1based: anchor.start_1based,
            end_1based: anchor.end_1based,
            strand: anchor.strand,
        })
    }

    pub fn list_sequence_genome_anchor_summaries(&self) -> Vec<SequenceGenomeAnchorSummary> {
        let mut seq_ids: Vec<String> = self.state.sequences.keys().cloned().collect();
        seq_ids.sort_unstable();
        seq_ids
            .into_iter()
            .filter_map(|seq_id| self.sequence_genome_anchor_summary(&seq_id).ok())
            .collect()
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
                "RenderFeatureExpertSvg".to_string(),
                "RenderRnaStructureSvg".to_string(),
                "RenderLineageSvg".to_string(),
                "RenderPoolGelSvg".to_string(),
                "CreateArrangementSerial".to_string(),
                "ExportDnaLadders".to_string(),
                "ExportRnaLadders".to_string(),
                "ExportPool".to_string(),
                "PrepareGenome".to_string(),
                "ExtractGenomeRegion".to_string(),
                "ExtractGenomeGene".to_string(),
                "ExtendGenomeAnchor".to_string(),
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
                "GenerateCandidateSetBetweenAnchors".to_string(),
                "DeleteCandidateSet".to_string(),
                "UpsertGuideSet".to_string(),
                "DeleteGuideSet".to_string(),
                "FilterGuidesPractical".to_string(),
                "GenerateGuideOligos".to_string(),
                "ExportGuideOligos".to_string(),
                "ExportGuideProtocolText".to_string(),
                "ScoreCandidateSetExpression".to_string(),
                "ScoreCandidateSetDistance".to_string(),
                "FilterCandidateSet".to_string(),
                "CandidateSetOp".to_string(),
                "ScoreCandidateSetWeightedObjective".to_string(),
                "TopKCandidateSet".to_string(),
                "ParetoFrontierCandidateSet".to_string(),
                "UpsertWorkflowMacroTemplate".to_string(),
                "DeleteWorkflowMacroTemplate".to_string(),
                "UpsertCandidateMacroTemplate".to_string(),
                "DeleteCandidateMacroTemplate".to_string(),
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

    fn feature_qualifier_text(feature: &gb_io::seq::Feature, key: &str) -> Option<String> {
        feature
            .qualifier_values(key.into())
            .map(|value| value.split_whitespace().collect::<Vec<_>>().join(" "))
            .map(|value| value.trim().to_string())
            .find(|value| !value.is_empty())
    }

    fn feature_qualifier_f64(feature: &gb_io::seq::Feature, key: &str) -> Option<f64> {
        feature
            .qualifier_values(key.into())
            .next()
            .and_then(|v| v.trim().parse::<f64>().ok())
    }

    fn first_nonempty_feature_qualifier(
        feature: &gb_io::seq::Feature,
        keys: &[&str],
    ) -> Option<String> {
        for key in keys {
            if let Some(value) = Self::feature_qualifier_text(feature, key) {
                return Some(value);
            }
        }
        None
    }

    fn is_tfbs_feature(feature: &gb_io::seq::Feature) -> bool {
        matches!(
            feature.kind.to_string().to_ascii_uppercase().as_str(),
            "TFBS" | "TF_BINDING_SITE" | "PROTEIN_BIND"
        )
    }

    fn feature_display_label(feature: &gb_io::seq::Feature, feature_id: usize) -> String {
        let fallback = format!("{} #{}", feature.kind.to_string(), feature_id + 1);
        let label = Self::first_nonempty_feature_qualifier(
            feature,
            &[
                "label",
                "name",
                "standard_name",
                "gene",
                "protein_id",
                "product",
                "region_name",
                "bound_moiety",
            ],
        )
        .unwrap_or(fallback);
        let trimmed = label.trim();
        if trimmed.is_empty() {
            format!("{} #{}", feature.kind.to_string(), feature_id + 1)
        } else {
            trimmed.to_string()
        }
    }

    fn normalize_column_frequencies(counts: [f64; 4]) -> [f64; 4] {
        let total = counts
            .iter()
            .copied()
            .map(|v| if v.is_finite() && v > 0.0 { v } else { 0.0 })
            .sum::<f64>();
        if total <= 0.0 {
            [0.25_f64, 0.25_f64, 0.25_f64, 0.25_f64]
        } else {
            [
                (counts[0].max(0.0) / total).clamp(0.0, 1.0),
                (counts[1].max(0.0) / total).clamp(0.0, 1.0),
                (counts[2].max(0.0) / total).clamp(0.0, 1.0),
                (counts[3].max(0.0) / total).clamp(0.0, 1.0),
            ]
        }
    }

    fn information_content_bits(frequencies: &[f64; 4]) -> f64 {
        let entropy = frequencies
            .iter()
            .copied()
            .filter(|p| *p > 0.0)
            .map(|p| -p * p.log2())
            .sum::<f64>();
        (2.0 - entropy).clamp(0.0, 2.0)
    }

    fn resolve_tfbs_scoring_motif(
        feature: &gb_io::seq::Feature,
    ) -> Result<(String, Option<String>, Vec<[f64; 4]>), EngineError> {
        let label = Self::feature_qualifier_text(feature, "label");
        let mut tokens = vec![];
        if let Some(tf_id) = Self::feature_qualifier_text(feature, "tf_id") {
            tokens.push(tf_id);
        }
        if let Some(bound) = Self::first_nonempty_feature_qualifier(
            feature,
            &["bound_moiety", "standard_name", "gene", "name"],
        ) {
            tokens.push(bound);
        }
        if let Some(raw_label) = label {
            let trimmed = raw_label.trim().to_string();
            if !trimmed.is_empty() {
                tokens.push(trimmed.clone());
                let label_upper = trimmed.to_ascii_uppercase();
                if let Some(stripped) = label_upper.strip_prefix("TFBS ") {
                    tokens.push(stripped.trim().to_string());
                }
            }
        }

        let mut last_error: Option<EngineError> = None;
        for token in tokens {
            if token.trim().is_empty() {
                continue;
            }
            match Self::resolve_tf_motif_for_scoring(&token) {
                Ok((tf_id, tf_name, _, matrix_counts)) => {
                    return Ok((tf_id, tf_name, matrix_counts));
                }
                Err(err) => {
                    last_error = Some(err);
                }
            }
        }
        Err(last_error.unwrap_or(EngineError {
            code: ErrorCode::InvalidInput,
            message:
                "Could not resolve TF motif for feature (missing tf_id/label resolvable in motif registry)"
                    .to_string(),
        }))
    }

    fn build_tfbs_expert_view(
        &self,
        seq_id: &str,
        feature_id: usize,
    ) -> Result<TfbsExpertView, EngineError> {
        let dna = self
            .state
            .sequences
            .get(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{seq_id}' not found"),
            })?;
        let feature = dna.features().get(feature_id).ok_or_else(|| EngineError {
            code: ErrorCode::NotFound,
            message: format!(
                "Feature id '{}' was not found in sequence '{}'",
                feature_id, seq_id
            ),
        })?;
        if !Self::is_tfbs_feature(feature) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Feature '{}' in '{}' is not a TFBS/protein-bind feature",
                    feature_id, seq_id
                ),
            });
        }
        let (from, _to) = feature.location.find_bounds().map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!("Could not read TFBS feature range: {e}"),
        })?;
        if from < 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "TFBS feature has negative range bounds".to_string(),
            });
        }

        let (tf_id, tf_name, matrix_counts) = Self::resolve_tfbs_scoring_motif(feature)?;
        if matrix_counts.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Resolved motif '{tf_id}' has empty matrix"),
            });
        }
        let motif_length = matrix_counts.len();
        let start = from as usize;
        let end = start.saturating_add(motif_length);
        let mut matched_bytes = dna.get_range_safe(start..end).ok_or_else(|| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!(
                "TFBS window {}..{} for feature '{}' exceeds sequence bounds",
                start + 1,
                end,
                feature_id
            ),
        })?;
        let is_reverse = feature_is_reverse(feature);
        if is_reverse {
            matched_bytes = Self::reverse_complement_bytes(&matched_bytes);
        }
        let matched_sequence = String::from_utf8_lossy(&matched_bytes).to_ascii_uppercase();

        let (llr_matrix, true_log_odds_matrix) = Self::prepare_scoring_matrices(&matrix_counts);
        let mut columns: Vec<TfbsExpertColumn> = Vec::with_capacity(motif_length);
        for idx in 0..motif_length {
            let counts = matrix_counts[idx];
            let frequencies = Self::normalize_column_frequencies(counts);
            let information_content_bits = Self::information_content_bits(&frequencies);
            let match_base = matched_bytes
                .get(idx)
                .map(|b| (*b as char).to_ascii_uppercase());
            let match_idx = matched_bytes.get(idx).and_then(|b| Self::base_to_idx(*b));
            let (match_frequency, llr_bits, true_log_odds_bits) = if let Some(base_idx) = match_idx
            {
                (
                    Some(frequencies[base_idx]),
                    Some(llr_matrix[idx][base_idx]),
                    Some(true_log_odds_matrix[idx][base_idx]),
                )
            } else {
                (None, None, None)
            };
            columns.push(TfbsExpertColumn {
                index_1based: idx + 1,
                counts,
                frequencies,
                information_content_bits,
                match_base,
                match_frequency,
                llr_bits,
                true_log_odds_bits,
                llr_rank_desc: None,
            });
        }

        let mut ranked = columns
            .iter()
            .enumerate()
            .filter_map(|(idx, col)| col.llr_bits.map(|score| (idx, score)))
            .collect::<Vec<_>>();
        ranked.sort_by(|a, b| b.1.total_cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
        for (rank_idx, (col_idx, _)) in ranked.iter().enumerate() {
            if let Some(col) = columns.get_mut(*col_idx) {
                col.llr_rank_desc = Some(rank_idx + 1);
            }
        }

        let llr_from_columns = columns.iter().filter_map(|col| col.llr_bits).sum::<f64>();
        let true_log_odds_from_columns = columns
            .iter()
            .filter_map(|col| col.true_log_odds_bits)
            .sum::<f64>();
        let llr_total_bits =
            Self::feature_qualifier_f64(feature, "llr_bits").or(Some(llr_from_columns));
        let llr_quantile = Self::feature_qualifier_f64(feature, "llr_quantile");
        let true_log_odds_total_bits = Self::feature_qualifier_f64(feature, "true_log_odds_bits")
            .or_else(|| Self::feature_qualifier_f64(feature, "log_odds_ratio_bits"))
            .or(Some(true_log_odds_from_columns));
        let true_log_odds_quantile = Self::feature_qualifier_f64(feature, "true_log_odds_quantile")
            .or_else(|| Self::feature_qualifier_f64(feature, "log_odds_ratio_quantile"));

        Ok(TfbsExpertView {
            seq_id: seq_id.to_string(),
            feature_id,
            feature_label: Self::feature_display_label(feature, feature_id),
            tf_id,
            tf_name,
            strand: if is_reverse {
                "-".to_string()
            } else {
                "+".to_string()
            },
            start_1based: start + 1,
            end_1based: end,
            motif_length,
            matched_sequence,
            llr_total_bits,
            llr_quantile,
            true_log_odds_total_bits,
            true_log_odds_quantile,
            instruction: TFBS_EXPERT_INSTRUCTION.to_string(),
            columns,
        })
    }

    fn resolve_restriction_site_target(
        &self,
        seq_id: &str,
        cut_pos_1based: usize,
        enzyme: Option<&str>,
        recognition_start_1based: Option<usize>,
        recognition_end_1based: Option<usize>,
    ) -> Result<(RestrictionEnzymeKey, Vec<String>, Option<String>), EngineError> {
        if cut_pos_1based == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Restriction-site cut_pos_1based must be >= 1".to_string(),
            });
        }
        let dna = self
            .state
            .sequences
            .get(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{seq_id}' not found"),
            })?;
        let enzyme_filter = enzyme
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .map(|v| v.to_ascii_uppercase());
        let mut candidates: Vec<(RestrictionEnzymeKey, Vec<String>)> = Vec::new();
        for (key, names) in dna.restriction_enzyme_groups() {
            if key.pos() < 0 || (key.pos() as usize + 1) != cut_pos_1based {
                continue;
            }
            if let Some(start) = recognition_start_1based {
                if key.from() < 0 || (key.from() as usize + 1) != start {
                    continue;
                }
            }
            if let Some(end) = recognition_end_1based {
                if key.to() < 0 || key.to() as usize != end {
                    continue;
                }
            }
            if let Some(filter) = &enzyme_filter {
                if !names
                    .iter()
                    .any(|name| name.to_ascii_uppercase() == *filter)
                {
                    continue;
                }
            }
            candidates.push((key.clone(), names.clone()));
        }
        if candidates.is_empty() {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "No restriction-site match found for cut position {} in '{}'",
                    cut_pos_1based, seq_id
                ),
            });
        }
        if candidates.len() > 1 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Restriction-site target is ambiguous for cut position {} in '{}'; provide enzyme and/or recognition start/end",
                    cut_pos_1based, seq_id
                ),
            });
        }
        let (key, mut names) = candidates.into_iter().next().unwrap_or_default();
        names.sort_by(|a, b| a.to_ascii_uppercase().cmp(&b.to_ascii_uppercase()));
        names.dedup_by(|left, right| left.eq_ignore_ascii_case(right));
        let selected_enzyme = if let Some(filter) = enzyme_filter {
            names
                .iter()
                .find(|name| name.to_ascii_uppercase() == filter)
                .cloned()
        } else {
            names.first().cloned()
        };
        Ok((key, names, selected_enzyme))
    }

    fn complement_iupac_text(seq: &str) -> String {
        seq.as_bytes()
            .iter()
            .map(|b| {
                Self::iupac_letter_complement(*b)
                    .unwrap_or(b'N')
                    .to_ascii_uppercase() as char
            })
            .collect()
    }

    fn build_restriction_site_expert_view(
        &self,
        seq_id: &str,
        cut_pos_1based: usize,
        enzyme: Option<&str>,
        recognition_start_1based: Option<usize>,
        recognition_end_1based: Option<usize>,
    ) -> Result<RestrictionSiteExpertView, EngineError> {
        let (key, enzyme_names, selected_enzyme) = self.resolve_restriction_site_target(
            seq_id,
            cut_pos_1based,
            enzyme,
            recognition_start_1based,
            recognition_end_1based,
        )?;
        let dna = self
            .state
            .sequences
            .get(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{seq_id}' not found"),
            })?;
        let start_1based = key.from().max(0) as usize + 1;
        let end_1based = key.to().max(0) as usize;
        let start0 = start_1based.saturating_sub(1);
        let end0 = end_1based;
        let mut site_sequence = dna
            .get_range_safe(start0..end0)
            .map(|bytes| String::from_utf8_lossy(&bytes).to_ascii_uppercase())
            .unwrap_or_default();
        let selected_enzyme_def = selected_enzyme.as_ref().and_then(|name| {
            dna.restriction_enzymes()
                .iter()
                .find(|enzyme| enzyme.name.eq_ignore_ascii_case(name))
        });
        let recognition_iupac = selected_enzyme_def.map(|enzyme| enzyme.sequence.clone());
        if site_sequence.is_empty() {
            if let Some(iupac) = &recognition_iupac {
                site_sequence = iupac.to_ascii_uppercase();
            }
        }
        let site_sequence_complement = Self::complement_iupac_text(&site_sequence);
        let max_cut = site_sequence.chars().count();
        let cut_index_0based = key.cut_size().max(0) as usize;
        let cut_index_0based = cut_index_0based.min(max_cut);

        Ok(RestrictionSiteExpertView {
            seq_id: seq_id.to_string(),
            cut_pos_1based,
            recognition_start_1based: start_1based,
            recognition_end_1based: end_1based,
            cut_index_0based,
            number_of_cuts_for_enzyme: key.number_of_cuts(),
            selected_enzyme,
            enzyme_names,
            recognition_iupac,
            site_sequence,
            site_sequence_complement,
            overlap_bp: selected_enzyme_def.map(|enzyme| enzyme.overlap),
            instruction: RESTRICTION_EXPERT_INSTRUCTION.to_string(),
        })
    }

    fn is_mrna_feature(feature: &gb_io::seq::Feature) -> bool {
        feature.kind.to_string().eq_ignore_ascii_case("mRNA")
    }

    fn is_exon_feature(feature: &gb_io::seq::Feature) -> bool {
        feature.kind.to_string().eq_ignore_ascii_case("exon")
    }

    fn splicing_group_label(feature: &gb_io::seq::Feature, fallback: usize) -> String {
        Self::first_nonempty_feature_qualifier(
            feature,
            &[
                "gene",
                "gene_id",
                "locus_tag",
                "standard_name",
                "gene_synonym",
            ],
        )
        .unwrap_or_else(|| format!("mRNA-group-{}", fallback + 1))
    }

    fn feature_transcript_id(feature: &gb_io::seq::Feature, fallback: usize) -> String {
        Self::first_nonempty_feature_qualifier(
            feature,
            &[
                "transcript_id",
                "standard_name",
                "product",
                "label",
                "name",
                "gene",
            ],
        )
        .unwrap_or_else(|| format!("transcript-{}", fallback + 1))
    }

    fn range_vec_to_splicing(mut ranges: Vec<(usize, usize)>) -> Vec<SplicingRange> {
        ranges.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        ranges
            .into_iter()
            .filter(|(start, end)| end > start)
            .map(|(start, end)| SplicingRange {
                start_1based: start + 1,
                end_1based: end,
            })
            .collect()
    }

    fn sequence_slice_upper(dna: &DNAsequence, start: usize, end: usize) -> String {
        if end <= start {
            return String::new();
        }
        dna.get_range_safe(start..end)
            .map(|bytes| String::from_utf8_lossy(&bytes).to_ascii_uppercase())
            .unwrap_or_default()
    }

    fn splice_boundary_markers_for_introns(
        dna: &DNAsequence,
        transcript_feature_id: usize,
        transcript_id: &str,
        is_reverse: bool,
        introns: &[(usize, usize)],
    ) -> Vec<SplicingBoundaryMarker> {
        let mut out = Vec::new();
        for (start, end) in introns {
            if *end <= *start || end - start < 2 {
                continue;
            }
            if is_reverse {
                let donor_raw = Self::sequence_slice_upper(dna, end.saturating_sub(2), *end);
                let acceptor_raw = Self::sequence_slice_upper(dna, *start, (start + 2).min(*end));
                let donor_rc = Self::reverse_complement_bytes(donor_raw.as_bytes());
                let acceptor_rc = Self::reverse_complement_bytes(acceptor_raw.as_bytes());
                let donor_motif = String::from_utf8_lossy(&donor_rc).to_ascii_uppercase();
                let acceptor_motif = String::from_utf8_lossy(&acceptor_rc).to_ascii_uppercase();
                out.push(SplicingBoundaryMarker {
                    transcript_feature_id,
                    transcript_id: transcript_id.to_string(),
                    side: "donor".to_string(),
                    position_1based: *end,
                    canonical: donor_motif == "GT",
                    motif_2bp: donor_motif,
                });
                out.push(SplicingBoundaryMarker {
                    transcript_feature_id,
                    transcript_id: transcript_id.to_string(),
                    side: "acceptor".to_string(),
                    position_1based: start + 1,
                    canonical: acceptor_motif == "AG",
                    motif_2bp: acceptor_motif,
                });
            } else {
                let donor_motif = Self::sequence_slice_upper(dna, *start, (start + 2).min(*end));
                let acceptor_motif = Self::sequence_slice_upper(dna, end.saturating_sub(2), *end);
                out.push(SplicingBoundaryMarker {
                    transcript_feature_id,
                    transcript_id: transcript_id.to_string(),
                    side: "donor".to_string(),
                    position_1based: start + 1,
                    canonical: donor_motif == "GT",
                    motif_2bp: donor_motif,
                });
                out.push(SplicingBoundaryMarker {
                    transcript_feature_id,
                    transcript_id: transcript_id.to_string(),
                    side: "acceptor".to_string(),
                    position_1based: *end,
                    canonical: acceptor_motif == "AG",
                    motif_2bp: acceptor_motif,
                });
            }
        }
        out
    }

    fn build_splicing_expert_view(
        &self,
        seq_id: &str,
        feature_id: usize,
    ) -> Result<SplicingExpertView, EngineError> {
        let dna = self
            .state
            .sequences
            .get(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{seq_id}' not found"),
            })?;
        let features = dna.features();
        let target_feature = features.get(feature_id).ok_or_else(|| EngineError {
            code: ErrorCode::NotFound,
            message: format!(
                "Feature id '{}' was not found in sequence '{}'",
                feature_id, seq_id
            ),
        })?;
        if !Self::is_mrna_feature(target_feature) && !Self::is_exon_feature(target_feature) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Feature '{}' in '{}' is not mRNA/exon and cannot seed a splicing view",
                    feature_id, seq_id
                ),
            });
        }
        let target_group = Self::splicing_group_label(target_feature, feature_id);
        let target_is_reverse = feature_is_reverse(target_feature);

        #[derive(Clone)]
        struct TranscriptWork {
            feature_id: usize,
            transcript_id: String,
            label: String,
            is_reverse: bool,
            exon_ranges: Vec<(usize, usize)>,
            intron_ranges: Vec<(usize, usize)>,
        }

        let mut transcripts: Vec<TranscriptWork> = Vec::new();
        for (idx, feature) in features.iter().enumerate() {
            if !Self::is_mrna_feature(feature) {
                continue;
            }
            if Self::splicing_group_label(feature, idx) != target_group {
                continue;
            }
            if feature_is_reverse(feature) != target_is_reverse {
                continue;
            }
            let mut exon_ranges = vec![];
            collect_location_ranges_usize(&feature.location, &mut exon_ranges);
            if exon_ranges.is_empty() {
                let (from, to) = feature.location.find_bounds().map_err(|e| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("Could not parse transcript range: {e}"),
                })?;
                if from >= 0 && to >= 0 {
                    exon_ranges.push((from as usize, to as usize));
                }
            }
            exon_ranges.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
            exon_ranges.retain(|(start, end)| end > start);
            if exon_ranges.is_empty() {
                continue;
            }
            let mut introns: Vec<(usize, usize)> = vec![];
            for pair in exon_ranges.windows(2) {
                let left = pair[0];
                let right = pair[1];
                if right.0 > left.1 {
                    introns.push((left.1, right.0));
                }
            }
            transcripts.push(TranscriptWork {
                feature_id: idx,
                transcript_id: Self::feature_transcript_id(feature, idx),
                label: Self::feature_display_label(feature, idx),
                is_reverse: feature_is_reverse(feature),
                exon_ranges,
                intron_ranges: introns,
            });
        }

        if transcripts.is_empty() {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "No mRNA transcripts found for splicing group '{}' in '{}'",
                    target_group, seq_id
                ),
            });
        }

        transcripts.sort_by(|a, b| a.transcript_id.cmp(&b.transcript_id));

        let mut exon_support: HashMap<(usize, usize), HashSet<usize>> = HashMap::new();
        let mut junction_support: HashMap<(usize, usize), HashSet<usize>> = HashMap::new();
        for transcript in &transcripts {
            for exon in &transcript.exon_ranges {
                exon_support
                    .entry(*exon)
                    .or_default()
                    .insert(transcript.feature_id);
            }
            for intron in &transcript.intron_ranges {
                junction_support
                    .entry(*intron)
                    .or_default()
                    .insert(transcript.feature_id);
            }
        }

        let mut unique_exons: Vec<(usize, usize)> = exon_support.keys().copied().collect();
        unique_exons.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        let region_start = unique_exons
            .iter()
            .map(|(start, _)| *start)
            .min()
            .unwrap_or(0);
        let region_end = unique_exons
            .iter()
            .map(|(_, end)| *end)
            .max()
            .unwrap_or(region_start + 1);

        let transcript_count = transcripts.len();
        let unique_exon_summaries: Vec<SplicingExonSummary> = unique_exons
            .iter()
            .map(|(start, end)| {
                let support = exon_support
                    .get(&(*start, *end))
                    .map(|set| set.len())
                    .unwrap_or(0);
                SplicingExonSummary {
                    start_1based: start + 1,
                    end_1based: *end,
                    support_transcript_count: support,
                    constitutive: support == transcript_count,
                }
            })
            .collect();

        let transcript_lanes: Vec<SplicingTranscriptLane> = transcripts
            .iter()
            .map(|transcript| SplicingTranscriptLane {
                transcript_feature_id: transcript.feature_id,
                transcript_id: transcript.transcript_id.clone(),
                label: transcript.label.clone(),
                strand: if transcript.is_reverse {
                    "-".to_string()
                } else {
                    "+".to_string()
                },
                exons: Self::range_vec_to_splicing(transcript.exon_ranges.clone()),
                introns: Self::range_vec_to_splicing(transcript.intron_ranges.clone()),
                has_target_feature: transcript.feature_id == feature_id,
            })
            .collect();

        let matrix_rows: Vec<SplicingMatrixRow> = transcripts
            .iter()
            .map(|transcript| {
                let set = transcript
                    .exon_ranges
                    .iter()
                    .copied()
                    .collect::<HashSet<(usize, usize)>>();
                SplicingMatrixRow {
                    transcript_feature_id: transcript.feature_id,
                    transcript_id: transcript.transcript_id.clone(),
                    label: transcript.label.clone(),
                    exon_presence: unique_exons.iter().map(|exon| set.contains(exon)).collect(),
                }
            })
            .collect();

        let mut boundaries = Vec::new();
        for transcript in &transcripts {
            boundaries.extend(Self::splice_boundary_markers_for_introns(
                dna,
                transcript.feature_id,
                &transcript.transcript_id,
                transcript.is_reverse,
                &transcript.intron_ranges,
            ));
        }
        boundaries.sort_by(|a, b| {
            a.position_1based
                .cmp(&b.position_1based)
                .then_with(|| a.side.cmp(&b.side))
                .then_with(|| a.transcript_feature_id.cmp(&b.transcript_feature_id))
        });

        let mut junctions: Vec<SplicingJunctionArc> = junction_support
            .into_iter()
            .map(|((donor0, acceptor0), ids)| {
                let mut transcript_feature_ids = ids.into_iter().collect::<Vec<_>>();
                transcript_feature_ids.sort_unstable();
                SplicingJunctionArc {
                    donor_1based: donor0 + 1,
                    acceptor_1based: acceptor0,
                    support_transcript_count: transcript_feature_ids.len(),
                    transcript_feature_ids,
                }
            })
            .collect();
        junctions.sort_by(|a, b| {
            a.donor_1based
                .cmp(&b.donor_1based)
                .then(a.acceptor_1based.cmp(&b.acceptor_1based))
        });

        let mut exon_skip_details = Vec::new();
        let mut intron_retention_details = Vec::new();
        let mut alt5_details = Vec::new();
        let mut alt3_details = Vec::new();
        let mut mex_details = Vec::new();

        let exon_set_by_transcript: HashMap<usize, HashSet<(usize, usize)>> = transcripts
            .iter()
            .map(|transcript| {
                (
                    transcript.feature_id,
                    transcript.exon_ranges.iter().copied().collect(),
                )
            })
            .collect();

        let mut alt5_groups = 0usize;
        let mut alt3_groups = 0usize;
        let mut grouped_by_start: HashMap<usize, Vec<(usize, usize)>> = HashMap::new();
        let mut grouped_by_end: HashMap<usize, Vec<(usize, usize)>> = HashMap::new();
        for exon in &unique_exons {
            grouped_by_start.entry(exon.0).or_default().push(*exon);
            grouped_by_end.entry(exon.1).or_default().push(*exon);
        }
        if target_is_reverse {
            for (start, group) in grouped_by_start {
                let distinct = group
                    .iter()
                    .map(|(_, end)| *end)
                    .collect::<HashSet<_>>()
                    .len();
                if distinct > 1 {
                    alt5_groups += 1;
                    alt5_details.push(format!(
                        "shared start={} with {} alternative 5' exon end(s)",
                        start + 1,
                        distinct
                    ));
                }
            }
            for (end, group) in grouped_by_end {
                let distinct = group
                    .iter()
                    .map(|(start, _)| *start)
                    .collect::<HashSet<_>>()
                    .len();
                if distinct > 1 {
                    alt3_groups += 1;
                    alt3_details.push(format!(
                        "shared end={} with {} alternative 3' exon start(s)",
                        end, distinct
                    ));
                }
            }
        } else {
            for (start, group) in grouped_by_start {
                let distinct = group
                    .iter()
                    .map(|(_, end)| *end)
                    .collect::<HashSet<_>>()
                    .len();
                if distinct > 1 {
                    alt3_groups += 1;
                    alt3_details.push(format!(
                        "shared start={} with {} alternative 3' exon end(s)",
                        start + 1,
                        distinct
                    ));
                }
            }
            for (end, group) in grouped_by_end {
                let distinct = group
                    .iter()
                    .map(|(start, _)| *start)
                    .collect::<HashSet<_>>()
                    .len();
                if distinct > 1 {
                    alt5_groups += 1;
                    alt5_details.push(format!(
                        "shared end={} with {} alternative 5' exon start(s)",
                        end, distinct
                    ));
                }
            }
        }

        let mut exon_skip_count = 0usize;
        for transcript in &transcripts {
            for (intron_start, intron_end) in &transcript.intron_ranges {
                let skipped = unique_exons
                    .iter()
                    .filter(|(exon_start, exon_end)| {
                        *exon_start >= *intron_start
                            && *exon_end <= *intron_end
                            && !transcript
                                .exon_ranges
                                .iter()
                                .any(|own| own.0 == *exon_start && own.1 == *exon_end)
                    })
                    .count();
                if skipped > 0 {
                    exon_skip_count += skipped;
                    exon_skip_details.push(format!(
                        "{} junction {}..{} skips {} exon(s)",
                        transcript.transcript_id,
                        intron_start + 1,
                        intron_end,
                        skipped
                    ));
                }
            }
        }

        let mut intron_retention_count = 0usize;
        for transcript in &transcripts {
            for (intron_start, intron_end) in &transcript.intron_ranges {
                let retained_by = transcripts.iter().find_map(|other| {
                    if other.feature_id == transcript.feature_id {
                        return None;
                    }
                    if other.exon_ranges.iter().any(|(exon_start, exon_end)| {
                        exon_start <= intron_start && exon_end >= intron_end
                    }) {
                        Some(other.transcript_id.clone())
                    } else {
                        None
                    }
                });
                if let Some(other_id) = retained_by {
                    intron_retention_count += 1;
                    intron_retention_details.push(format!(
                        "{} intron {}..{} retained in {}",
                        transcript.transcript_id,
                        intron_start + 1,
                        intron_end,
                        other_id
                    ));
                }
            }
        }

        let mut mex_count = 0usize;
        for left_idx in 0..unique_exons.len() {
            for right_idx in (left_idx + 1)..unique_exons.len() {
                let left = unique_exons[left_idx];
                let right = unique_exons[right_idx];
                if left.1 <= right.0 || right.1 <= left.0 {
                    continue;
                }
                let left_support = exon_support.get(&left).cloned().unwrap_or_default();
                let right_support = exon_support.get(&right).cloned().unwrap_or_default();
                if left_support.is_empty()
                    || right_support.is_empty()
                    || !left_support.is_disjoint(&right_support)
                {
                    continue;
                }
                let co_occurs = transcripts.iter().any(|transcript| {
                    let set = exon_set_by_transcript
                        .get(&transcript.feature_id)
                        .cloned()
                        .unwrap_or_default();
                    set.contains(&left) && set.contains(&right)
                });
                if !co_occurs {
                    mex_count += 1;
                    mex_details.push(format!(
                        "{}..{} and {}..{} appear mutually exclusive",
                        left.0 + 1,
                        left.1,
                        right.0 + 1,
                        right.1
                    ));
                }
            }
        }

        let alt_exon_count = unique_exon_summaries
            .iter()
            .filter(|exon| !exon.constitutive)
            .count();

        let events = vec![
            SplicingEventSummary {
                event_type: "alternative_exon".to_string(),
                count: alt_exon_count,
                details: unique_exon_summaries
                    .iter()
                    .filter(|exon| !exon.constitutive)
                    .take(6)
                    .map(|exon| {
                        format!(
                            "{}..{} support={}/{}",
                            exon.start_1based,
                            exon.end_1based,
                            exon.support_transcript_count,
                            transcript_count
                        )
                    })
                    .collect(),
            },
            SplicingEventSummary {
                event_type: "alternative_5prime".to_string(),
                count: alt5_groups,
                details: alt5_details.into_iter().take(6).collect(),
            },
            SplicingEventSummary {
                event_type: "alternative_3prime".to_string(),
                count: alt3_groups,
                details: alt3_details.into_iter().take(6).collect(),
            },
            SplicingEventSummary {
                event_type: "exon_skipping".to_string(),
                count: exon_skip_count,
                details: exon_skip_details.into_iter().take(6).collect(),
            },
            SplicingEventSummary {
                event_type: "intron_retention".to_string(),
                count: intron_retention_count,
                details: intron_retention_details.into_iter().take(6).collect(),
            },
            SplicingEventSummary {
                event_type: "mutually_exclusive_exons".to_string(),
                count: mex_count,
                details: mex_details.into_iter().take(6).collect(),
            },
        ];

        Ok(SplicingExpertView {
            seq_id: seq_id.to_string(),
            target_feature_id: feature_id,
            group_label: target_group,
            strand: if target_is_reverse {
                "-".to_string()
            } else {
                "+".to_string()
            },
            region_start_1based: region_start + 1,
            region_end_1based: region_end,
            transcript_count,
            unique_exon_count: unique_exons.len(),
            instruction: SPLICING_EXPERT_INSTRUCTION.to_string(),
            transcripts: transcript_lanes,
            unique_exons: unique_exon_summaries,
            matrix_rows,
            boundaries,
            junctions,
            events,
        })
    }

    pub fn inspect_feature_expert(
        &self,
        seq_id: &str,
        target: &FeatureExpertTarget,
    ) -> Result<FeatureExpertView, EngineError> {
        match target {
            FeatureExpertTarget::TfbsFeature { feature_id } => self
                .build_tfbs_expert_view(seq_id, *feature_id)
                .map(FeatureExpertView::Tfbs),
            FeatureExpertTarget::RestrictionSite {
                cut_pos_1based,
                enzyme,
                recognition_start_1based,
                recognition_end_1based,
            } => self
                .build_restriction_site_expert_view(
                    seq_id,
                    *cut_pos_1based,
                    enzyme.as_deref(),
                    *recognition_start_1based,
                    *recognition_end_1based,
                )
                .map(FeatureExpertView::RestrictionSite),
            FeatureExpertTarget::SplicingFeature { feature_id } => self
                .build_splicing_expert_view(seq_id, *feature_id)
                .map(FeatureExpertView::Splicing),
        }
    }

    pub fn render_feature_expert_svg_to_path(
        &self,
        seq_id: &str,
        target: &FeatureExpertTarget,
        path: &str,
    ) -> Result<FeatureExpertView, EngineError> {
        let view = self.inspect_feature_expert(seq_id, target)?;
        let svg = render_feature_expert_svg(&view);
        std::fs::write(path, svg).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write feature-expert SVG to '{path}': {e}"),
        })?;
        Ok(view)
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
        timeout_seconds: Option<u64>,
        on_progress: &mut dyn FnMut(PrepareGenomeProgress) -> bool,
    ) -> Result<PrepareGenomeReport, EngineError> {
        let (catalog, _) = Self::open_reference_genome_catalog(catalog_path)?;
        let timeout = timeout_seconds.map(Duration::from_secs);
        let started = Instant::now();
        let mut timed_out = false;
        let mut guarded_progress = |progress: PrepareGenomeProgress| -> bool {
            if let Some(limit) = timeout {
                if started.elapsed() >= limit {
                    timed_out = true;
                    return false;
                }
            }
            on_progress(progress)
        };
        catalog
            .prepare_genome_once_with_progress(
                genome_id,
                cache_dir.map(str::trim).filter(|v| !v.is_empty()),
                &mut guarded_progress,
            )
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: if is_prepare_cancelled_error(&e) {
                    if timed_out {
                        format!(
                            "Genome preparation timed out for '{genome_id}' after {} second(s)",
                            timeout_seconds.unwrap_or(0)
                        )
                    } else {
                        format!("Genome preparation cancelled for '{genome_id}'")
                    }
                } else {
                    format!("Could not prepare genome '{genome_id}': {e}")
                },
            })
    }

    pub fn prepare_helper_genome_once(
        genome_id: &str,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
        timeout_seconds: Option<u64>,
        on_progress: &mut dyn FnMut(PrepareGenomeProgress) -> bool,
    ) -> Result<PrepareGenomeReport, EngineError> {
        let chosen = catalog_path
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .unwrap_or(DEFAULT_HELPER_GENOME_CATALOG_PATH);
        Self::prepare_reference_genome_once(
            genome_id,
            Some(chosen),
            cache_dir,
            timeout_seconds,
            on_progress,
        )
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

    fn history_limit_or_default(&self) -> usize {
        if self.history_limit == 0 {
            Self::default_history_limit()
        } else {
            self.history_limit
        }
    }

    fn capture_history_checkpoint(&self) -> EngineHistoryCheckpoint {
        EngineHistoryCheckpoint {
            state: self.state.clone(),
            journal: self.journal.clone(),
            op_counter: self.op_counter,
        }
    }

    fn restore_history_checkpoint(&mut self, checkpoint: EngineHistoryCheckpoint) {
        self.state = checkpoint.state;
        self.journal = checkpoint.journal;
        self.op_counter = checkpoint.op_counter;
        self.reconcile_lineage_nodes();
        self.reconcile_containers();
    }

    fn op_records_history_checkpoint(op: &Operation) -> bool {
        !matches!(
            op,
            Operation::SaveFile { .. }
                | Operation::RenderSequenceSvg { .. }
                | Operation::RenderFeatureExpertSvg { .. }
                | Operation::RenderRnaStructureSvg { .. }
                | Operation::RenderLineageSvg { .. }
                | Operation::RenderPoolGelSvg { .. }
                | Operation::ExportDnaLadders { .. }
                | Operation::ExportRnaLadders { .. }
                | Operation::ExportPool { .. }
        )
    }

    fn maybe_capture_checkpoint(&self, op: &Operation) -> Option<EngineHistoryCheckpoint> {
        Self::op_records_history_checkpoint(op).then(|| self.capture_history_checkpoint())
    }

    fn push_undo_checkpoint(&mut self, checkpoint: EngineHistoryCheckpoint) {
        self.undo_stack.push(checkpoint);
        let limit = self.history_limit_or_default();
        if self.undo_stack.len() > limit {
            let drain_len = self.undo_stack.len() - limit;
            self.undo_stack.drain(0..drain_len);
        }
        self.redo_stack.clear();
    }

    pub fn undo_available(&self) -> usize {
        self.undo_stack.len()
    }

    pub fn redo_available(&self) -> usize {
        self.redo_stack.len()
    }

    pub fn undo_last_operation(&mut self) -> Result<(), EngineError> {
        let Some(previous) = self.undo_stack.pop() else {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: "No operation to undo".to_string(),
            });
        };
        let current = self.capture_history_checkpoint();
        self.redo_stack.push(current);
        let limit = self.history_limit_or_default();
        if self.redo_stack.len() > limit {
            let drain_len = self.redo_stack.len() - limit;
            self.redo_stack.drain(0..drain_len);
        }
        self.restore_history_checkpoint(previous);
        Ok(())
    }

    pub fn redo_last_operation(&mut self) -> Result<(), EngineError> {
        let Some(next) = self.redo_stack.pop() else {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: "No operation to redo".to_string(),
            });
        };
        let current = self.capture_history_checkpoint();
        self.undo_stack.push(current);
        let limit = self.history_limit_or_default();
        if self.undo_stack.len() > limit {
            let drain_len = self.undo_stack.len() - limit;
            self.undo_stack.drain(0..drain_len);
        }
        self.restore_history_checkpoint(next);
        Ok(())
    }

    pub fn apply_with_progress<F>(
        &mut self,
        op: Operation,
        mut on_progress: F,
    ) -> Result<OpResult, EngineError>
    where
        F: FnMut(OperationProgress) -> bool,
    {
        let run_id = "interactive".to_string();
        let checkpoint = self.maybe_capture_checkpoint(&op);
        let result = self.apply_internal(op.clone(), &run_id, &mut on_progress)?;
        self.journal.push(OperationRecord {
            run_id,
            op,
            result: result.clone(),
        });
        if let Some(checkpoint) = checkpoint {
            self.push_undo_checkpoint(checkpoint);
        }
        Ok(result)
    }

    pub fn apply_workflow_with_progress<F>(
        &mut self,
        wf: Workflow,
        mut on_progress: F,
    ) -> Result<Vec<OpResult>, EngineError>
    where
        F: FnMut(OperationProgress) -> bool,
    {
        let mut results = Vec::new();
        for op in &wf.ops {
            let checkpoint = self.maybe_capture_checkpoint(op);
            let result = self.apply_internal(op.clone(), &wf.run_id, &mut on_progress)?;
            self.journal.push(OperationRecord {
                run_id: wf.run_id.clone(),
                op: op.clone(),
                result: result.clone(),
            });
            if let Some(checkpoint) = checkpoint {
                self.push_undo_checkpoint(checkpoint);
            }
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

    fn read_guide_design_store_from_metadata(
        value: Option<&serde_json::Value>,
    ) -> GuideDesignStore {
        let mut store = value
            .cloned()
            .and_then(|v| serde_json::from_value::<GuideDesignStore>(v).ok())
            .unwrap_or_default();
        if store.schema.trim().is_empty() {
            store.schema = GUIDE_DESIGN_SCHEMA.to_string();
        }
        store
    }

    fn read_guide_design_store(&self) -> GuideDesignStore {
        Self::read_guide_design_store_from_metadata(
            self.state.metadata.get(GUIDE_DESIGN_METADATA_KEY),
        )
    }

    fn write_guide_design_store(&mut self, mut store: GuideDesignStore) -> Result<(), EngineError> {
        if store.guide_sets.is_empty()
            && store.practical_filter_reports.is_empty()
            && store.oligo_sets.is_empty()
            && store.latest_oligo_set_by_guide_set.is_empty()
            && store.audit_log.is_empty()
        {
            self.state.metadata.remove(GUIDE_DESIGN_METADATA_KEY);
            return Ok(());
        }
        store.schema = GUIDE_DESIGN_SCHEMA.to_string();
        store.updated_at_unix_ms = Self::now_unix_ms();
        let value = serde_json::to_value(store).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize guide-design metadata: {e}"),
        })?;
        self.state
            .metadata
            .insert(GUIDE_DESIGN_METADATA_KEY.to_string(), value);
        Ok(())
    }

    fn normalize_guide_set_id(raw: &str) -> Result<String, EngineError> {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "guide_set_id cannot be empty".to_string(),
            });
        }
        Ok(trimmed.to_string())
    }

    fn normalize_guide_id(raw: &str, index: usize) -> String {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            format!("g_{:04}", index + 1)
        } else {
            trimmed.to_string()
        }
    }

    fn normalize_guide_strand(raw: &str) -> Result<String, EngineError> {
        match raw.trim().to_ascii_lowercase().as_str() {
            "+" | "plus" | "forward" | "fwd" => Ok("+".to_string()),
            "-" | "minus" | "reverse" | "rev" => Ok("-".to_string()),
            other => Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Unsupported guide strand '{}'; expected '+' or '-'", other),
            }),
        }
    }

    fn normalize_guide_candidate(
        guide: GuideCandidate,
        index: usize,
    ) -> Result<GuideCandidate, EngineError> {
        if guide.end_0based_exclusive <= guide.start_0based {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Guide {} has invalid genomic interval {}..{}",
                    index + 1,
                    guide.start_0based,
                    guide.end_0based_exclusive
                ),
            });
        }
        if guide.seq_id.trim().is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Guide {} has empty seq_id", index + 1),
            });
        }
        if matches!(guide.rank, Some(0)) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Guide {} rank must be >= 1", index + 1),
            });
        }
        let protospacer = Self::normalize_iupac_text(&guide.protospacer)?;
        if protospacer.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Guide {} has empty protospacer", index + 1),
            });
        }
        let pam = Self::normalize_iupac_text(&guide.pam)?;
        if pam.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Guide {} has empty PAM", index + 1),
            });
        }
        if guide.cut_offset_from_protospacer_start >= protospacer.len() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Guide {} cut_offset_from_protospacer_start ({}) must be < protospacer length ({})",
                    index + 1,
                    guide.cut_offset_from_protospacer_start,
                    protospacer.len()
                ),
            });
        }
        Ok(GuideCandidate {
            guide_id: Self::normalize_guide_id(&guide.guide_id, index),
            seq_id: guide.seq_id.trim().to_string(),
            start_0based: guide.start_0based,
            end_0based_exclusive: guide.end_0based_exclusive,
            strand: Self::normalize_guide_strand(&guide.strand)?,
            protospacer,
            pam,
            nuclease: if guide.nuclease.trim().is_empty() {
                "SpCas9".to_string()
            } else {
                guide.nuclease.trim().to_string()
            },
            cut_offset_from_protospacer_start: guide.cut_offset_from_protospacer_start,
            rank: guide.rank,
        })
    }

    fn normalize_guide_candidates(
        guides: Vec<GuideCandidate>,
    ) -> Result<Vec<GuideCandidate>, EngineError> {
        if guides.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "guide set requires at least one guide".to_string(),
            });
        }
        let mut normalized = Vec::with_capacity(guides.len());
        let mut seen_ids = HashSet::new();
        for (idx, guide) in guides.into_iter().enumerate() {
            let guide = Self::normalize_guide_candidate(guide, idx)?;
            if !seen_ids.insert(guide.guide_id.clone()) {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("Duplicate guide_id '{}' in guide set", guide.guide_id),
                });
            }
            normalized.push(guide);
        }
        normalized.sort_by(|a, b| {
            a.rank
                .unwrap_or(usize::MAX)
                .cmp(&b.rank.unwrap_or(usize::MAX))
                .then(a.guide_id.cmp(&b.guide_id))
        });
        Ok(normalized)
    }

    fn max_homopolymer_run_for_base(sequence: &[u8], base: u8) -> usize {
        let base = base.to_ascii_uppercase();
        let mut best = 0usize;
        let mut current = 0usize;
        for b in sequence {
            if b.to_ascii_uppercase() == base {
                current += 1;
                best = best.max(current);
            } else {
                current = 0;
            }
        }
        best
    }

    fn max_dinucleotide_repeat_units(sequence: &[u8]) -> usize {
        if sequence.len() < 2 {
            return 0;
        }
        let canonical = |b: u8| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T');
        let mut best = 1usize;
        for start in 0..(sequence.len() - 1) {
            let b0 = sequence[start].to_ascii_uppercase();
            let b1 = sequence[start + 1].to_ascii_uppercase();
            if !canonical(b0) || !canonical(b1) {
                continue;
            }
            let mut units = 1usize;
            let mut idx = start + 2;
            while idx + 1 < sequence.len()
                && sequence[idx].to_ascii_uppercase() == b0
                && sequence[idx + 1].to_ascii_uppercase() == b1
            {
                units += 1;
                idx += 2;
            }
            best = best.max(units);
        }
        best
    }

    fn normalize_practical_filter_config(
        mut config: GuidePracticalFilterConfig,
    ) -> Result<GuidePracticalFilterConfig, EngineError> {
        if let Some(min) = config.gc_min {
            if !(0.0..=1.0).contains(&min) {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("gc_min ({min}) must be between 0.0 and 1.0"),
                });
            }
        }
        if let Some(max) = config.gc_max {
            if !(0.0..=1.0).contains(&max) {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("gc_max ({max}) must be between 0.0 and 1.0"),
                });
            }
        }
        if let (Some(min), Some(max)) = (config.gc_min, config.gc_max) {
            if min > max {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("gc_min ({min}) must be <= gc_max ({max})"),
                });
            }
        }
        if let Some(max_run) = config.max_homopolymer_run {
            if max_run == 0 {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: "max_homopolymer_run must be >= 1".to_string(),
                });
            }
        }
        if let Some(max_repeat) = config.max_dinucleotide_repeat_units {
            if max_repeat == 0 {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: "max_dinucleotide_repeat_units must be >= 1".to_string(),
                });
            }
        }

        let mut normalized_per_base = HashMap::new();
        for (base, value) in &config.max_homopolymer_run_per_base {
            let key = base.trim().to_ascii_uppercase();
            if !matches!(key.as_str(), "A" | "C" | "G" | "T" | "U") {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Unsupported max_homopolymer_run_per_base key '{}'; expected A/C/G/T",
                        base
                    ),
                });
            }
            if *value == 0 {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("max_homopolymer_run_per_base for '{}' must be >= 1", base),
                });
            }
            let key = if key == "U" {
                "T".to_string()
            } else {
                key.to_string()
            };
            normalized_per_base.insert(key, *value);
        }
        config.max_homopolymer_run_per_base = normalized_per_base;

        if let Some(required) = config.required_5prime_base.as_ref() {
            let normalized = Self::normalize_iupac_text(required)?;
            if normalized.len() != 1
                || !matches!(normalized.as_bytes()[0], b'A' | b'C' | b'G' | b'T' | b'U')
            {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: "required_5prime_base must be one canonical nucleotide (A/C/G/T)"
                        .to_string(),
                });
            }
            config.required_5prime_base = Some(if normalized == "U" {
                "T".to_string()
            } else {
                normalized
            });
        }

        let mut motifs = vec![];
        for motif in &config.forbidden_motifs {
            let normalized = Self::normalize_iupac_text(motif)?;
            if !normalized.is_empty() {
                motifs.push(normalized);
            }
        }
        motifs.sort();
        motifs.dedup();
        config.forbidden_motifs = motifs;
        Ok(config)
    }

    pub fn list_guide_sets(&self) -> Vec<GuideSetSummary> {
        let store = self.read_guide_design_store();
        let mut names = store.guide_sets.keys().cloned().collect::<Vec<_>>();
        names.sort();
        names
            .into_iter()
            .filter_map(|name| store.guide_sets.get(&name).cloned())
            .map(|set| GuideSetSummary {
                guide_set_id: set.guide_set_id,
                guide_count: set.guides.len(),
                created_at_unix_ms: set.created_at_unix_ms,
                updated_at_unix_ms: set.updated_at_unix_ms,
            })
            .collect()
    }

    pub fn inspect_guide_set_page(
        &self,
        guide_set_id: &str,
        limit: usize,
        offset: usize,
    ) -> Result<(GuideSet, usize, usize), EngineError> {
        if limit == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Guide page limit must be >= 1".to_string(),
            });
        }
        let guide_set_id = Self::normalize_guide_set_id(guide_set_id)?;
        let store = self.read_guide_design_store();
        let mut set = store
            .guide_sets
            .get(&guide_set_id)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Guide set '{}' not found", guide_set_id),
            })?;
        let total = set.guides.len();
        let clamped_offset = offset.min(total);
        set.guides = set
            .guides
            .into_iter()
            .skip(clamped_offset)
            .take(limit)
            .collect();
        Ok((set, total, clamped_offset))
    }

    pub fn get_guide_practical_filter_report(
        &self,
        guide_set_id: &str,
    ) -> Result<GuidePracticalFilterReport, EngineError> {
        let guide_set_id = Self::normalize_guide_set_id(guide_set_id)?;
        let store = self.read_guide_design_store();
        store
            .practical_filter_reports
            .get(&guide_set_id)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "No practical guide filter report found for '{}'",
                    guide_set_id
                ),
            })
    }

    pub fn list_guide_oligo_sets(&self, guide_set_id: Option<&str>) -> Vec<GuideOligoSet> {
        let store = self.read_guide_design_store();
        let filter = guide_set_id.map(|v| v.trim().to_string());
        let mut ids = store.oligo_sets.keys().cloned().collect::<Vec<_>>();
        ids.sort();
        ids.into_iter()
            .filter_map(|id| store.oligo_sets.get(&id).cloned())
            .filter(|set| {
                filter
                    .as_ref()
                    .map(|f| set.guide_set_id == *f)
                    .unwrap_or(true)
            })
            .collect()
    }

    pub fn get_guide_oligo_set(&self, oligo_set_id: &str) -> Result<GuideOligoSet, EngineError> {
        let id = oligo_set_id.trim();
        if id.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "oligo_set_id cannot be empty".to_string(),
            });
        }
        let store = self.read_guide_design_store();
        store
            .oligo_sets
            .get(id)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Guide oligo set '{}' not found", id),
            })
    }

    fn built_in_guide_oligo_template(template_id: &str) -> Option<GuideOligoTemplateSpec> {
        match template_id.trim() {
            "lenti_bsmbi_u6_default" => Some(GuideOligoTemplateSpec {
                template_id: "lenti_bsmbi_u6_default".to_string(),
                description: "U6 sgRNA cloning oligos with BsmBI overhangs".to_string(),
                forward_prefix: "CACC".to_string(),
                forward_suffix: "".to_string(),
                reverse_prefix: "AAAC".to_string(),
                reverse_suffix: "C".to_string(),
                reverse_uses_reverse_complement_of_spacer: true,
                uppercase_output: true,
            }),
            "plain_forward_reverse" => Some(GuideOligoTemplateSpec {
                template_id: "plain_forward_reverse".to_string(),
                description: "Raw spacer and reverse-complement spacer".to_string(),
                forward_prefix: "".to_string(),
                forward_suffix: "".to_string(),
                reverse_prefix: "".to_string(),
                reverse_suffix: "".to_string(),
                reverse_uses_reverse_complement_of_spacer: true,
                uppercase_output: true,
            }),
            _ => None,
        }
    }

    fn normalize_oligo_set_id(raw: &str) -> Result<String, EngineError> {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "oligo_set_id cannot be empty".to_string(),
            });
        }
        Ok(trimmed.to_string())
    }

    fn unique_oligo_set_id(store: &GuideDesignStore, base: &str) -> String {
        if !store.oligo_sets.contains_key(base) {
            return base.to_string();
        }
        let mut idx = 2usize;
        loop {
            let candidate = format!("{base}_{idx}");
            if !store.oligo_sets.contains_key(&candidate) {
                return candidate;
            }
            idx += 1;
        }
    }

    fn resolve_oligo_set_for_export(
        store: &GuideDesignStore,
        guide_set_id: &str,
        requested_oligo_set_id: Option<&str>,
    ) -> Result<GuideOligoSet, EngineError> {
        if let Some(id) = requested_oligo_set_id {
            let id = id.trim();
            if id.is_empty() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: "oligo_set_id cannot be empty".to_string(),
                });
            }
            let set = store.oligo_sets.get(id).ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Guide oligo set '{}' not found", id),
            })?;
            if set.guide_set_id != guide_set_id {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Guide oligo set '{}' belongs to '{}' (expected '{}')",
                        id, set.guide_set_id, guide_set_id
                    ),
                });
            }
            return Ok(set.clone());
        }

        let latest = store
            .latest_oligo_set_by_guide_set
            .get(guide_set_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "No oligo set is registered for guide set '{}'; run GenerateGuideOligos first",
                    guide_set_id
                ),
            })?;
        store
            .oligo_sets
            .get(latest)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "Guide oligo set '{}' referenced by '{}' is missing",
                    latest, guide_set_id
                ),
            })
    }

    fn normalize_workflow_macro_template_name(raw: &str) -> Result<String, EngineError> {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Workflow macro template name cannot be empty".to_string(),
            });
        }
        Ok(trimmed.to_string())
    }

    fn normalize_workflow_macro_param_name(raw: &str) -> Result<String, EngineError> {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Workflow macro parameter name cannot be empty".to_string(),
            });
        }
        let valid = trimmed.chars().enumerate().all(|(idx, ch)| match idx {
            0 => ch.is_ascii_alphabetic() || ch == '_',
            _ => ch.is_ascii_alphanumeric() || ch == '_',
        });
        if !valid {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Invalid workflow macro parameter name '{}' (expected [A-Za-z_][A-Za-z0-9_]*)",
                    trimmed
                ),
            });
        }
        Ok(trimmed.to_string())
    }

    fn normalize_workflow_macro_template_details_url(
        raw: Option<String>,
    ) -> Result<Option<String>, EngineError> {
        let Some(raw) = raw else {
            return Ok(None);
        };
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return Ok(None);
        }
        let lower = trimmed.to_ascii_lowercase();
        if !(lower.starts_with("https://") || lower.starts_with("http://")) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Workflow macro template details_url '{}' must start with http:// or https://",
                    trimmed
                ),
            });
        }
        Ok(Some(trimmed.to_string()))
    }

    fn read_workflow_macro_template_store_from_metadata(
        value: Option<&serde_json::Value>,
    ) -> WorkflowMacroTemplateStore {
        let mut store = value
            .cloned()
            .and_then(|v| serde_json::from_value::<WorkflowMacroTemplateStore>(v).ok())
            .unwrap_or_default();
        if store.schema.trim().is_empty() {
            store.schema = WORKFLOW_MACRO_TEMPLATES_SCHEMA.to_string();
        }
        store
    }

    fn read_workflow_macro_template_store(&self) -> WorkflowMacroTemplateStore {
        Self::read_workflow_macro_template_store_from_metadata(
            self.state
                .metadata
                .get(WORKFLOW_MACRO_TEMPLATES_METADATA_KEY),
        )
    }

    fn write_workflow_macro_template_store(
        &mut self,
        mut store: WorkflowMacroTemplateStore,
    ) -> Result<(), EngineError> {
        if store.templates.is_empty() {
            self.state
                .metadata
                .remove(WORKFLOW_MACRO_TEMPLATES_METADATA_KEY);
            return Ok(());
        }
        store.schema = WORKFLOW_MACRO_TEMPLATES_SCHEMA.to_string();
        store.updated_at_unix_ms = Self::now_unix_ms();
        let value = serde_json::to_value(store).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize workflow macro template metadata: {e}"),
        })?;
        self.state
            .metadata
            .insert(WORKFLOW_MACRO_TEMPLATES_METADATA_KEY.to_string(), value);
        Ok(())
    }

    pub fn list_workflow_macro_templates(&self) -> Vec<WorkflowMacroTemplateSummary> {
        let store = self.read_workflow_macro_template_store();
        let mut names = store.templates.keys().cloned().collect::<Vec<_>>();
        names.sort();
        names
            .into_iter()
            .filter_map(|name| store.templates.get(&name))
            .map(|template| WorkflowMacroTemplateSummary {
                name: template.name.clone(),
                description: template.description.clone(),
                details_url: template.details_url.clone(),
                parameter_count: template.parameters.len(),
                created_at_unix_ms: template.created_at_unix_ms,
                updated_at_unix_ms: template.updated_at_unix_ms,
            })
            .collect()
    }

    pub fn get_workflow_macro_template(
        &self,
        name: &str,
    ) -> Result<WorkflowMacroTemplate, EngineError> {
        let name = Self::normalize_workflow_macro_template_name(name)?;
        let store = self.read_workflow_macro_template_store();
        store
            .templates
            .get(&name)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Workflow macro template '{}' not found", name),
            })
    }

    pub fn render_workflow_macro_template_script(
        &self,
        name: &str,
        bindings: &HashMap<String, String>,
    ) -> Result<String, EngineError> {
        let template = self.get_workflow_macro_template(name)?;
        let declared = template
            .parameters
            .iter()
            .map(|p| p.name.clone())
            .collect::<HashSet<_>>();
        for key in bindings.keys() {
            if !declared.contains(key) {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Workflow macro template '{}' does not define parameter '{}'",
                        template.name, key
                    ),
                });
            }
        }

        let mut resolved: HashMap<String, String> = HashMap::new();
        for param in &template.parameters {
            if let Some(value) = bindings.get(&param.name) {
                resolved.insert(param.name.clone(), value.clone());
            } else if let Some(default_value) = &param.default_value {
                resolved.insert(param.name.clone(), default_value.clone());
            } else if param.required {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Workflow macro template '{}' is missing required parameter '{}'",
                        template.name, param.name
                    ),
                });
            }
        }

        let placeholder_regex =
            Regex::new(r"\$\{([A-Za-z_][A-Za-z0-9_]*)\}").map_err(|e| EngineError {
                code: ErrorCode::Internal,
                message: format!("Could not compile workflow macro placeholder regex: {e}"),
            })?;
        let mut missing: Vec<String> = vec![];
        for captures in placeholder_regex.captures_iter(&template.script) {
            if let Some(name) = captures.get(1).map(|m| m.as_str()) {
                if !declared.contains(name) {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Workflow macro template '{}' references undeclared parameter '{}'",
                            template.name, name
                        ),
                    });
                }
                if !resolved.contains_key(name) {
                    missing.push(name.to_string());
                }
            }
        }
        if !missing.is_empty() {
            missing.sort();
            missing.dedup();
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Workflow macro template '{}' is missing parameter bindings for: {}",
                    template.name,
                    missing.join(", ")
                ),
            });
        }

        let rendered = placeholder_regex
            .replace_all(&template.script, |captures: &regex::Captures<'_>| {
                let key = captures.get(1).map(|m| m.as_str()).unwrap_or_default();
                resolved.get(key).cloned().unwrap_or_default()
            })
            .to_string();
        Ok(rendered)
    }

    fn normalize_candidate_macro_template_name(raw: &str) -> Result<String, EngineError> {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Candidate macro template name cannot be empty".to_string(),
            });
        }
        Ok(trimmed.to_string())
    }

    fn normalize_candidate_macro_param_name(raw: &str) -> Result<String, EngineError> {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Candidate macro parameter name cannot be empty".to_string(),
            });
        }
        let valid = trimmed.chars().enumerate().all(|(idx, ch)| match idx {
            0 => ch.is_ascii_alphabetic() || ch == '_',
            _ => ch.is_ascii_alphanumeric() || ch == '_',
        });
        if !valid {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Invalid candidate macro parameter name '{}' (expected [A-Za-z_][A-Za-z0-9_]*)",
                    trimmed
                ),
            });
        }
        Ok(trimmed.to_string())
    }

    fn normalize_candidate_macro_template_details_url(
        raw: Option<String>,
    ) -> Result<Option<String>, EngineError> {
        let Some(raw) = raw else {
            return Ok(None);
        };
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return Ok(None);
        }
        let lower = trimmed.to_ascii_lowercase();
        if !(lower.starts_with("https://") || lower.starts_with("http://")) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Candidate macro template details_url '{}' must start with http:// or https://",
                    trimmed
                ),
            });
        }
        Ok(Some(trimmed.to_string()))
    }

    fn read_candidate_macro_template_store_from_metadata(
        value: Option<&serde_json::Value>,
    ) -> CandidateMacroTemplateStore {
        let mut store = value
            .cloned()
            .and_then(|v| serde_json::from_value::<CandidateMacroTemplateStore>(v).ok())
            .unwrap_or_default();
        if store.schema.trim().is_empty() {
            store.schema = CANDIDATE_MACRO_TEMPLATES_SCHEMA.to_string();
        }
        store
    }

    fn read_candidate_macro_template_store(&self) -> CandidateMacroTemplateStore {
        Self::read_candidate_macro_template_store_from_metadata(
            self.state
                .metadata
                .get(CANDIDATE_MACRO_TEMPLATES_METADATA_KEY),
        )
    }

    fn write_candidate_macro_template_store(
        &mut self,
        mut store: CandidateMacroTemplateStore,
    ) -> Result<(), EngineError> {
        if store.templates.is_empty() {
            self.state
                .metadata
                .remove(CANDIDATE_MACRO_TEMPLATES_METADATA_KEY);
            return Ok(());
        }
        store.schema = CANDIDATE_MACRO_TEMPLATES_SCHEMA.to_string();
        store.updated_at_unix_ms = Self::now_unix_ms();
        let value = serde_json::to_value(store).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize candidate macro template metadata: {e}"),
        })?;
        self.state
            .metadata
            .insert(CANDIDATE_MACRO_TEMPLATES_METADATA_KEY.to_string(), value);
        Ok(())
    }

    pub fn list_candidate_macro_templates(&self) -> Vec<CandidateMacroTemplateSummary> {
        let store = self.read_candidate_macro_template_store();
        let mut names = store.templates.keys().cloned().collect::<Vec<_>>();
        names.sort();
        names
            .into_iter()
            .filter_map(|name| store.templates.get(&name))
            .map(|template| CandidateMacroTemplateSummary {
                name: template.name.clone(),
                description: template.description.clone(),
                details_url: template.details_url.clone(),
                parameter_count: template.parameters.len(),
                created_at_unix_ms: template.created_at_unix_ms,
                updated_at_unix_ms: template.updated_at_unix_ms,
            })
            .collect()
    }

    pub fn get_candidate_macro_template(
        &self,
        name: &str,
    ) -> Result<CandidateMacroTemplate, EngineError> {
        let name = Self::normalize_candidate_macro_template_name(name)?;
        let store = self.read_candidate_macro_template_store();
        store
            .templates
            .get(&name)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Candidate macro template '{}' not found", name),
            })
    }

    pub fn render_candidate_macro_template_script(
        &self,
        name: &str,
        bindings: &HashMap<String, String>,
    ) -> Result<String, EngineError> {
        let template = self.get_candidate_macro_template(name)?;
        let declared = template
            .parameters
            .iter()
            .map(|p| p.name.clone())
            .collect::<HashSet<_>>();
        for key in bindings.keys() {
            if !declared.contains(key) {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Candidate macro template '{}' does not define parameter '{}'",
                        template.name, key
                    ),
                });
            }
        }

        let mut resolved: HashMap<String, String> = HashMap::new();
        for param in &template.parameters {
            if let Some(value) = bindings.get(&param.name) {
                resolved.insert(param.name.clone(), value.clone());
            } else if let Some(default_value) = &param.default_value {
                resolved.insert(param.name.clone(), default_value.clone());
            } else if param.required {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Candidate macro template '{}' is missing required parameter '{}'",
                        template.name, param.name
                    ),
                });
            }
        }

        let placeholder_regex =
            Regex::new(r"\$\{([A-Za-z_][A-Za-z0-9_]*)\}").map_err(|e| EngineError {
                code: ErrorCode::Internal,
                message: format!("Could not compile candidate macro placeholder regex: {e}"),
            })?;
        let mut missing: Vec<String> = vec![];
        for captures in placeholder_regex.captures_iter(&template.script) {
            if let Some(name) = captures.get(1).map(|m| m.as_str()) {
                if !declared.contains(name) {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Candidate macro template '{}' references undeclared parameter '{}'",
                            template.name, name
                        ),
                    });
                }
                if !resolved.contains_key(name) {
                    missing.push(name.to_string());
                }
            }
        }
        if !missing.is_empty() {
            missing.sort();
            missing.dedup();
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Candidate macro template '{}' is missing parameter bindings for: {}",
                    template.name,
                    missing.join(", ")
                ),
            });
        }

        let rendered = placeholder_regex
            .replace_all(&template.script, |captures: &regex::Captures<'_>| {
                let key = captures.get(1).map(|m| m.as_str()).unwrap_or_default();
                resolved.get(key).cloned().unwrap_or_default()
            })
            .to_string();
        Ok(rendered)
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

    fn collect_location_segments(
        location: &gb_io::seq::Location,
        seq_len: usize,
        reverse: bool,
        out: &mut Vec<FeatureLocationSegment>,
    ) {
        use gb_io::seq::Location;
        match location {
            Location::Range((a, _), (b, _)) => {
                let raw_start = (*a).min(*b);
                let raw_end = (*a).max(*b);
                let raw_end_exclusive = if raw_end == raw_start {
                    raw_end.saturating_add(1)
                } else {
                    raw_end
                };
                if raw_end_exclusive <= 0 || raw_start >= seq_len as i64 {
                    return;
                }
                let start = raw_start.max(0).min(seq_len as i64) as usize;
                let end_0based = raw_end_exclusive.max(0).min(seq_len as i64) as usize;
                if end_0based <= start {
                    return;
                }
                out.push(FeatureLocationSegment {
                    start_0based: start,
                    end_0based,
                    strand: Some(if reverse { '-' } else { '+' }),
                });
            }
            Location::Between(a, b) => {
                let raw_start = (*a).min(*b);
                let raw_end_exclusive = (*a).max(*b).saturating_add(1);
                if raw_end_exclusive <= 0 || raw_start >= seq_len as i64 {
                    return;
                }
                let start = raw_start.max(0).min(seq_len as i64) as usize;
                let end_0based = raw_end_exclusive.max(0).min(seq_len as i64) as usize;
                if end_0based <= start {
                    return;
                }
                out.push(FeatureLocationSegment {
                    start_0based: start,
                    end_0based,
                    strand: Some(if reverse { '-' } else { '+' }),
                });
            }
            Location::Complement(inner) => {
                Self::collect_location_segments(inner, seq_len, !reverse, out);
            }
            Location::Join(parts)
            | Location::Order(parts)
            | Location::Bond(parts)
            | Location::OneOf(parts) => {
                for part in parts {
                    Self::collect_location_segments(part, seq_len, reverse, out);
                }
            }
            Location::External(_, maybe_inner) => {
                if let Some(inner) = maybe_inner.as_deref() {
                    Self::collect_location_segments(inner, seq_len, reverse, out);
                }
            }
            Location::Gap(_) => {}
        }
    }

    fn feature_segment_boundary_positions(
        segment: &FeatureLocationSegment,
        boundary_mode: CandidateFeatureBoundaryMode,
    ) -> Vec<usize> {
        if segment.end_0based <= segment.start_0based {
            return vec![];
        }
        let start = segment.start_0based;
        let end = segment.end_0based.saturating_sub(1);
        let mut out = vec![];
        match boundary_mode {
            CandidateFeatureBoundaryMode::Any => {
                out.push(start);
                out.push(end);
            }
            CandidateFeatureBoundaryMode::Start => out.push(start),
            CandidateFeatureBoundaryMode::End => out.push(end),
            CandidateFeatureBoundaryMode::FivePrime => match segment.strand {
                Some('+') => out.push(start),
                Some('-') => out.push(end),
                _ => {
                    out.push(start);
                    out.push(end);
                }
            },
            CandidateFeatureBoundaryMode::ThreePrime => match segment.strand {
                Some('+') => out.push(end),
                Some('-') => out.push(start),
                _ => {
                    out.push(start);
                    out.push(end);
                }
            },
        }
        out.sort_unstable();
        out.dedup();
        out
    }

    fn collect_feature_distance_targets(
        dna: &DNAsequence,
        geometry_mode: CandidateFeatureGeometryMode,
        boundary_mode: CandidateFeatureBoundaryMode,
    ) -> Vec<FeatureDistanceTarget> {
        let mut out = vec![];
        for (feature_index, feature) in dna.features().iter().enumerate() {
            if feature.kind.to_string().eq_ignore_ascii_case("SOURCE") {
                continue;
            }
            let kind_upper = feature.kind.to_string().to_ascii_uppercase();
            let labels_upper = Self::feature_labels_upper(feature);
            let mut segments = vec![];
            Self::collect_location_segments(&feature.location, dna.len(), false, &mut segments);
            if segments.is_empty() {
                let Ok((from, to)) = feature.location.find_bounds() else {
                    continue;
                };
                if from < 0 || to <= 0 {
                    continue;
                }
                let start = from.min(to).max(0).min(dna.len() as i64) as usize;
                let raw_end = from.max(to);
                let raw_end_exclusive = if raw_end == from.min(to) {
                    raw_end.saturating_add(1)
                } else {
                    raw_end
                };
                let end_0based = raw_end_exclusive.max(0).min(dna.len() as i64) as usize;
                if end_0based <= start {
                    continue;
                }
                segments.push(FeatureLocationSegment {
                    start_0based: start,
                    end_0based,
                    strand: None,
                });
            }
            if segments.is_empty() {
                continue;
            }
            let feature_strand = {
                let mut strand = None;
                let mut ambiguous = false;
                for segment in &segments {
                    match (strand, segment.strand) {
                        (_, None) => {
                            ambiguous = true;
                            break;
                        }
                        (None, Some(s)) => strand = Some(s),
                        (Some(prev), Some(s)) if prev == s => {}
                        (Some(_), Some(_)) => {
                            ambiguous = true;
                            break;
                        }
                    }
                }
                if ambiguous { None } else { strand }
            };
            match geometry_mode {
                CandidateFeatureGeometryMode::FeatureSpan => {
                    let start = segments
                        .iter()
                        .map(|segment| segment.start_0based)
                        .min()
                        .unwrap_or(0);
                    let end_0based = segments
                        .iter()
                        .map(|segment| segment.end_0based)
                        .max()
                        .unwrap_or(start);
                    if end_0based <= start {
                        continue;
                    }
                    out.push(FeatureDistanceTarget {
                        feature_index,
                        kind_upper: kind_upper.clone(),
                        labels_upper: labels_upper.clone(),
                        start_0based: start,
                        end_0based,
                        strand: feature_strand,
                    });
                }
                CandidateFeatureGeometryMode::FeatureParts => {
                    for segment in segments {
                        if segment.end_0based <= segment.start_0based {
                            continue;
                        }
                        out.push(FeatureDistanceTarget {
                            feature_index,
                            kind_upper: kind_upper.clone(),
                            labels_upper: labels_upper.clone(),
                            start_0based: segment.start_0based,
                            end_0based: segment.end_0based,
                            strand: segment.strand,
                        });
                    }
                }
                CandidateFeatureGeometryMode::FeatureBoundaries => {
                    for segment in &segments {
                        for point in
                            Self::feature_segment_boundary_positions(segment, boundary_mode)
                        {
                            if point >= dna.len() {
                                continue;
                            }
                            out.push(FeatureDistanceTarget {
                                feature_index,
                                kind_upper: kind_upper.clone(),
                                labels_upper: labels_upper.clone(),
                                start_0based: point,
                                end_0based: point.saturating_add(1).min(dna.len()),
                                strand: segment.strand,
                            });
                        }
                    }
                }
            }
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
        feature: &FeatureDistanceTarget,
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

    fn feature_matches_strand_relation(
        feature: &FeatureDistanceTarget,
        strand_relation: CandidateFeatureStrandRelation,
    ) -> bool {
        match strand_relation {
            CandidateFeatureStrandRelation::Any => true,
            CandidateFeatureStrandRelation::Same => feature.strand == Some('+'),
            CandidateFeatureStrandRelation::Opposite => feature.strand == Some('-'),
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

    fn boundary_distance(
        candidate_start: usize,
        candidate_end: usize,
        boundary_pos: usize,
    ) -> usize {
        if boundary_pos < candidate_start {
            candidate_start.saturating_sub(boundary_pos)
        } else if boundary_pos > candidate_end {
            boundary_pos.saturating_sub(candidate_end)
        } else {
            0
        }
    }

    fn nearest_feature_distance(
        candidate_start: usize,
        candidate_end: usize,
        features: &[FeatureDistanceTarget],
        kind_filter_upper: &[String],
        label_regex: Option<&Regex>,
        strand_relation: CandidateFeatureStrandRelation,
    ) -> Option<usize> {
        features
            .iter()
            .filter(|feature| Self::feature_matches_filter(feature, kind_filter_upper, label_regex))
            .filter(|feature| Self::feature_matches_strand_relation(feature, strand_relation))
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

    fn matching_feature_count(
        features: &[FeatureDistanceTarget],
        kind_filter_upper: &[String],
        label_regex: Option<&Regex>,
        strand_relation: CandidateFeatureStrandRelation,
    ) -> usize {
        features
            .iter()
            .filter(|feature| Self::feature_matches_filter(feature, kind_filter_upper, label_regex))
            .filter(|feature| Self::feature_matches_strand_relation(feature, strand_relation))
            .map(|feature| feature.feature_index)
            .collect::<HashSet<_>>()
            .len()
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
                    });
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
                        });
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
        feature_geometry_mode: Option<CandidateFeatureGeometryMode>,
        feature_boundary_mode: Option<CandidateFeatureBoundaryMode>,
        feature_strand_relation: Option<CandidateFeatureStrandRelation>,
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
        let feature_geometry_mode = feature_geometry_mode.unwrap_or_default();
        let requested_boundary_mode = feature_boundary_mode.unwrap_or_default();
        let feature_strand_relation = feature_strand_relation.unwrap_or_default();
        let effective_boundary_mode =
            if feature_geometry_mode == CandidateFeatureGeometryMode::FeatureBoundaries {
                requested_boundary_mode
            } else {
                CandidateFeatureBoundaryMode::Any
            };
        if feature_geometry_mode != CandidateFeatureGeometryMode::FeatureBoundaries
            && feature_boundary_mode.is_some()
        {
            result.warnings.push(
                "feature_boundary_mode is ignored unless feature_geometry_mode=feature_boundaries"
                    .to_string(),
            );
        }
        let feature_targets = Self::collect_feature_distance_targets(
            dna,
            feature_geometry_mode,
            effective_boundary_mode,
        );
        let matching_feature_count = Self::matching_feature_count(
            &feature_targets,
            &kind_filter_upper,
            label_regex.as_ref(),
            feature_strand_relation,
        );
        let has_feature_filter = !kind_filter_upper.is_empty()
            || label_regex.is_some()
            || max_distance_bp.is_some()
            || feature_strand_relation != CandidateFeatureStrandRelation::Any;
        if has_feature_filter && matching_feature_count == 0 {
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
            let distance_any = Self::nearest_feature_distance(
                start,
                end,
                &feature_targets,
                &[],
                None,
                CandidateFeatureStrandRelation::Any,
            );
            let distance_filtered = Self::nearest_feature_distance(
                start,
                end,
                &feature_targets,
                &kind_filter_upper,
                label_regex.as_ref(),
                feature_strand_relation,
            );
            let selected_distance = if !kind_filter_upper.is_empty()
                || label_regex.is_some()
                || feature_strand_relation != CandidateFeatureStrandRelation::Any
            {
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
        if !feature_targets.is_empty() {
            result.messages.push(format!(
                "Candidate set '{}' feature distance mode: geometry='{}', boundary='{}', strand_relation='{}'",
                set_name,
                feature_geometry_mode.as_str(),
                effective_boundary_mode.as_str(),
                feature_strand_relation.as_str()
            ));
        }
        Ok(())
    }

    fn op_generate_candidate_set_between_anchors(
        &mut self,
        set_name: String,
        seq_id: SeqId,
        anchor_a: SequenceAnchor,
        anchor_b: SequenceAnchor,
        length_bp: usize,
        step_bp: usize,
        limit: Option<usize>,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let set_name = Self::normalize_candidate_set_name(&set_name)?;
        if length_bp == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "GenerateCandidateSetBetweenAnchors requires length_bp >= 1".to_string(),
            });
        }
        if step_bp == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "GenerateCandidateSetBetweenAnchors requires step_bp >= 1".to_string(),
            });
        }
        let limit = limit.unwrap_or(self.max_fragments_per_container());
        if limit == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "GenerateCandidateSetBetweenAnchors requires limit >= 1".to_string(),
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
        let is_circular = dna.is_circular();

        let anchor_a_pos = Self::resolve_sequence_anchor_position(dna, &anchor_a, "anchor_a")?;
        let anchor_b_pos = Self::resolve_sequence_anchor_position(dna, &anchor_b, "anchor_b")?;
        if anchor_a_pos == anchor_b_pos {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "GenerateCandidateSetBetweenAnchors requires distinct anchor positions"
                    .to_string(),
            });
        }

        let left = anchor_a_pos.min(anchor_b_pos);
        let right = anchor_a_pos.max(anchor_b_pos);
        let span = right.saturating_sub(left);
        if span < length_bp {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Anchor span {} bp is shorter than requested candidate length {} bp",
                    span, length_bp
                ),
            });
        }

        let mut candidates = vec![];
        let mut considered = 0usize;
        let mut truncated = false;
        let upper = right.saturating_sub(length_bp);
        let mut start = left;
        while start <= upper {
            let end = start + length_bp;
            considered += 1;
            let Some(fragment) = dna.get_range_safe(start..end) else {
                start = start.saturating_add(step_bp);
                continue;
            };
            let sequence = String::from_utf8_lossy(&fragment).to_string();
            let mut metrics =
                Self::compute_candidate_metrics(&fragment, start, end, dna.len(), None);
            let dist_a = Self::boundary_distance(start, end, anchor_a_pos);
            let dist_b = Self::boundary_distance(start, end, anchor_b_pos);
            metrics.insert("anchor_a_position_bp".to_string(), anchor_a_pos as f64);
            metrics.insert("anchor_b_position_bp".to_string(), anchor_b_pos as f64);
            metrics.insert("distance_to_anchor_a_bp".to_string(), dist_a as f64);
            metrics.insert("distance_to_anchor_b_bp".to_string(), dist_b as f64);
            metrics.insert(
                "distance_to_nearest_anchor_bp".to_string(),
                dist_a.min(dist_b) as f64,
            );
            metrics.insert("anchor_interval_start_bp".to_string(), left as f64);
            metrics.insert("anchor_interval_end_bp".to_string(), right as f64);
            metrics.insert("anchor_interval_span_bp".to_string(), span as f64);
            metrics.insert(
                "anchor_order_sign".to_string(),
                if anchor_a_pos <= anchor_b_pos {
                    1.0
                } else {
                    -1.0
                },
            );
            metrics.insert(
                "candidate_offset_from_anchor_a_bp".to_string(),
                if anchor_a_pos <= anchor_b_pos {
                    start.saturating_sub(anchor_a_pos) as f64
                } else {
                    anchor_a_pos.saturating_sub(end) as f64
                },
            );
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
                message: "No candidates generated between the two anchors".to_string(),
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
            "Generated candidate set '{}' between anchors on '{}' ({} candidates, {} windows considered, span={} bp)",
            set_name, seq_id, generated, considered, span
        ));
        result.messages.push(format!(
            "Anchor positions on '{}': anchor_a={}, anchor_b={}",
            seq_id, anchor_a_pos, anchor_b_pos
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
        if is_circular {
            result.warnings.push(
                "GenerateCandidateSetBetweenAnchors currently uses linearized anchor interval semantics on circular sequences"
                    .to_string(),
            );
        }
        Ok(())
    }

    fn append_guide_design_audit(
        store: &mut GuideDesignStore,
        operation: &str,
        guide_set_id: &str,
        payload: serde_json::Value,
    ) {
        store.audit_log.push(GuideDesignAuditEntry {
            unix_ms: Self::now_unix_ms(),
            operation: operation.to_string(),
            guide_set_id: guide_set_id.to_string(),
            payload,
        });
        if store.audit_log.len() > 512 {
            let drain = store.audit_log.len() - 512;
            store.audit_log.drain(0..drain);
        }
    }

    fn op_upsert_guide_set(
        &mut self,
        guide_set_id: String,
        guides: Vec<GuideCandidate>,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let guide_set_id = Self::normalize_guide_set_id(&guide_set_id)?;
        let guides = Self::normalize_guide_candidates(guides)?;
        let mut store = self.read_guide_design_store();
        let now = Self::now_unix_ms();
        let created_at_unix_ms = store
            .guide_sets
            .get(&guide_set_id)
            .map(|set| set.created_at_unix_ms)
            .unwrap_or(now);
        let replaced_existing = store
            .guide_sets
            .insert(
                guide_set_id.clone(),
                GuideSet {
                    guide_set_id: guide_set_id.clone(),
                    guides: guides.clone(),
                    created_at_unix_ms,
                    updated_at_unix_ms: now,
                },
            )
            .is_some();
        Self::append_guide_design_audit(
            &mut store,
            "UpsertGuideSet",
            &guide_set_id,
            json!({
                "guide_count": guides.len(),
                "replaced_existing": replaced_existing
            }),
        );
        self.write_guide_design_store(store)?;
        result.messages.push(format!(
            "Upserted guide set '{}' with {} guide(s)",
            guide_set_id,
            guides.len()
        ));
        if replaced_existing {
            result.warnings.push(format!(
                "Guide set '{}' replaced existing content",
                guide_set_id
            ));
        }
        Ok(())
    }

    fn op_delete_guide_set(
        &mut self,
        guide_set_id: String,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let guide_set_id = Self::normalize_guide_set_id(&guide_set_id)?;
        let mut store = self.read_guide_design_store();
        let removed = store.guide_sets.remove(&guide_set_id).is_some();
        let _ = store.practical_filter_reports.remove(&guide_set_id);
        let oligo_ids = store
            .oligo_sets
            .iter()
            .filter(|(_, set)| set.guide_set_id == guide_set_id)
            .map(|(id, _)| id.clone())
            .collect::<Vec<_>>();
        for oligo_id in &oligo_ids {
            let _ = store.oligo_sets.remove(oligo_id);
        }
        let _ = store.latest_oligo_set_by_guide_set.remove(&guide_set_id);
        if removed {
            Self::append_guide_design_audit(
                &mut store,
                "DeleteGuideSet",
                &guide_set_id,
                json!({
                    "removed_oligo_sets": oligo_ids.len()
                }),
            );
        }
        self.write_guide_design_store(store)?;
        if removed {
            result.messages.push(format!(
                "Deleted guide set '{}' and {} associated oligo set(s)",
                guide_set_id,
                oligo_ids.len()
            ));
        } else {
            result
                .warnings
                .push(format!("Guide set '{}' was not present", guide_set_id));
        }
        Ok(())
    }

    fn op_filter_guides_practical(
        &mut self,
        guide_set_id: String,
        config: GuidePracticalFilterConfig,
        output_guide_set_id: Option<String>,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let guide_set_id = Self::normalize_guide_set_id(&guide_set_id)?;
        let output_guide_set_id = output_guide_set_id
            .as_deref()
            .map(Self::normalize_guide_set_id)
            .transpose()?;
        let config = Self::normalize_practical_filter_config(config)?;
        let mut store = self.read_guide_design_store();
        let guide_set = store
            .guide_sets
            .get(&guide_set_id)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Guide set '{}' not found", guide_set_id),
            })?;
        if guide_set.guides.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Guide set '{}' is empty", guide_set_id),
            });
        }

        let mut passed_guides: Vec<GuideCandidate> = Vec::new();
        let mut rows: Vec<GuidePracticalFilterResult> = Vec::with_capacity(guide_set.guides.len());

        for guide in &guide_set.guides {
            let spacer = guide
                .protospacer
                .as_bytes()
                .iter()
                .map(|b| match b.to_ascii_uppercase() {
                    b'U' => b'T',
                    other => other,
                })
                .collect::<Vec<_>>();
            let tail = guide
                .pam
                .as_bytes()
                .iter()
                .map(|b| match b.to_ascii_uppercase() {
                    b'U' => b'T',
                    other => other,
                })
                .collect::<Vec<_>>();

            let mut row = GuidePracticalFilterResult {
                guide_id: guide.guide_id.clone(),
                passed: true,
                reasons: vec![],
                warnings: vec![],
                metrics: HashMap::new(),
            };

            if config.reject_ambiguous_bases && Self::has_ambiguous_bases(&spacer) {
                row.reasons.push(GuideFilterReason {
                    code: "contains_ambiguous_base".to_string(),
                    message: "Protospacer contains non-ACGT bases".to_string(),
                });
            }

            if config.gc_min.is_some() || config.gc_max.is_some() {
                if let Some(gc) = Self::sequence_gc_fraction(&spacer) {
                    row.metrics.insert("gc_fraction".to_string(), gc);
                    if let Some(min) = config.gc_min {
                        if gc < min {
                            row.reasons.push(GuideFilterReason {
                                code: "gc_too_low".to_string(),
                                message: format!(
                                    "GC fraction {:.3} is below minimum {:.3}",
                                    gc, min
                                ),
                            });
                        }
                    }
                    if let Some(max) = config.gc_max {
                        if gc > max {
                            row.reasons.push(GuideFilterReason {
                                code: "gc_too_high".to_string(),
                                message: format!(
                                    "GC fraction {:.3} is above maximum {:.3}",
                                    gc, max
                                ),
                            });
                        }
                    }
                }
            }

            let max_run = Self::max_homopolymer_run(&spacer);
            row.metrics
                .insert("max_homopolymer_run".to_string(), max_run as f64);
            if let Some(limit) = config.max_homopolymer_run {
                if max_run > limit {
                    row.reasons.push(GuideFilterReason {
                        code: "homopolymer_run_exceeded".to_string(),
                        message: format!(
                            "Max homopolymer run {} exceeds configured limit {}",
                            max_run, limit
                        ),
                    });
                }
            }
            for (base, limit) in &config.max_homopolymer_run_per_base {
                let probe = base.as_bytes()[0];
                let observed = Self::max_homopolymer_run_for_base(&spacer, probe);
                row.metrics.insert(
                    format!("max_homopolymer_run_{}", base.to_ascii_lowercase()),
                    observed as f64,
                );
                if observed > *limit {
                    row.reasons.push(GuideFilterReason {
                        code: "homopolymer_run_exceeded".to_string(),
                        message: format!(
                            "Homopolymer run for base '{}' ({}) exceeds configured limit {}",
                            base, observed, limit
                        ),
                    });
                }
            }

            if config.avoid_u6_terminator_tttt {
                let mut window = spacer.clone();
                if config.u6_terminator_window == GuideU6TerminatorWindow::SpacerPlusTail {
                    window.extend_from_slice(&tail);
                }
                if Self::contains_u6_terminator_t4(&window) {
                    row.reasons.push(GuideFilterReason {
                        code: "u6_terminator_t4".to_string(),
                        message: "Contains TTTT in configured U6 terminator scan window"
                            .to_string(),
                    });
                }
            }

            let max_repeat = Self::max_dinucleotide_repeat_units(&spacer);
            row.metrics.insert(
                "max_dinucleotide_repeat_units".to_string(),
                max_repeat as f64,
            );
            if let Some(limit) = config.max_dinucleotide_repeat_units {
                if max_repeat > limit {
                    row.reasons.push(GuideFilterReason {
                        code: "dinucleotide_repeat_exceeded".to_string(),
                        message: format!(
                            "Max dinucleotide repeat units {} exceeds configured limit {}",
                            max_repeat, limit
                        ),
                    });
                }
            }

            for motif in &config.forbidden_motifs {
                if Self::contains_motif_any_strand(&spacer, motif)? {
                    row.reasons.push(GuideFilterReason {
                        code: "forbidden_motif_present".to_string(),
                        message: format!("Protospacer contains forbidden motif '{}'", motif),
                    });
                }
            }

            if let Some(required_base) = config.required_5prime_base.as_ref() {
                if let Some(actual) = spacer.first().map(|b| (*b as char).to_string()) {
                    if actual != *required_base {
                        if config.allow_5prime_g_extension && required_base == "G" {
                            row.warnings.push(
                                "Missing required 5' G but can be rescued by 5' G extension"
                                    .to_string(),
                            );
                        } else {
                            row.reasons.push(GuideFilterReason {
                                code: "required_5prime_base_missing".to_string(),
                                message: format!(
                                    "Protospacer 5' base '{}' does not match required '{}'",
                                    actual, required_base
                                ),
                            });
                        }
                    }
                }
            }

            row.passed = row.reasons.is_empty();
            if row.passed {
                passed_guides.push(guide.clone());
            }
            rows.push(row);
        }

        let report = GuidePracticalFilterReport {
            guide_set_id: guide_set_id.clone(),
            generated_at_unix_ms: Self::now_unix_ms(),
            config: config.clone(),
            passed_count: passed_guides.len(),
            rejected_count: rows.len().saturating_sub(passed_guides.len()),
            results: rows,
        };
        store
            .practical_filter_reports
            .insert(guide_set_id.clone(), report.clone());

        if let Some(output_guide_set_id) = output_guide_set_id {
            let now = Self::now_unix_ms();
            let created_at_unix_ms = store
                .guide_sets
                .get(&output_guide_set_id)
                .map(|set| set.created_at_unix_ms)
                .unwrap_or(now);
            let replaced = store
                .guide_sets
                .insert(
                    output_guide_set_id.clone(),
                    GuideSet {
                        guide_set_id: output_guide_set_id.clone(),
                        guides: passed_guides.clone(),
                        created_at_unix_ms,
                        updated_at_unix_ms: now,
                    },
                )
                .is_some();
            result.messages.push(format!(
                "Wrote passed guides to set '{}' ({} guide(s))",
                output_guide_set_id,
                passed_guides.len()
            ));
            if replaced {
                result.warnings.push(format!(
                    "Guide set '{}' replaced existing content",
                    output_guide_set_id
                ));
            }
        }

        Self::append_guide_design_audit(
            &mut store,
            "FilterGuidesPractical",
            &guide_set_id,
            json!({
                "passed_count": report.passed_count,
                "rejected_count": report.rejected_count
            }),
        );
        self.write_guide_design_store(store)?;
        result.messages.push(format!(
            "Filtered guide set '{}' (passed {}, rejected {})",
            guide_set_id, report.passed_count, report.rejected_count
        ));
        Ok(())
    }

    fn op_generate_guide_oligos(
        &mut self,
        guide_set_id: String,
        template_id: String,
        apply_5prime_g_extension: Option<bool>,
        output_oligo_set_id: Option<String>,
        passed_only: Option<bool>,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let guide_set_id = Self::normalize_guide_set_id(&guide_set_id)?;
        let template = Self::built_in_guide_oligo_template(&template_id).ok_or_else(|| {
            EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "Unknown guide oligo template '{}'; supported templates: lenti_bsmbi_u6_default, plain_forward_reverse",
                    template_id
                ),
            }
        })?;
        let apply_5prime_g_extension = apply_5prime_g_extension.unwrap_or(false);
        let passed_only = passed_only.unwrap_or(false);

        let mut store = self.read_guide_design_store();
        let guide_set = store
            .guide_sets
            .get(&guide_set_id)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Guide set '{}' not found", guide_set_id),
            })?;
        if guide_set.guides.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Guide set '{}' is empty", guide_set_id),
            });
        }

        let passed_lookup = if passed_only {
            let report = store
                .practical_filter_reports
                .get(&guide_set_id)
                .ok_or_else(|| EngineError {
                    code: ErrorCode::NotFound,
                    message: format!(
                        "No practical filter report found for '{}' (required because passed_only=true)",
                        guide_set_id
                    ),
                })?;
            Some(
                report
                    .results
                    .iter()
                    .filter(|r| r.passed)
                    .map(|r| r.guide_id.clone())
                    .collect::<HashSet<_>>(),
            )
        } else {
            None
        };

        let mut guides = guide_set.guides.clone();
        guides.sort_by(|a, b| {
            a.rank
                .unwrap_or(usize::MAX)
                .cmp(&b.rank.unwrap_or(usize::MAX))
                .then(a.guide_id.cmp(&b.guide_id))
        });
        let mut records = vec![];
        for guide in &guides {
            if let Some(lookup) = &passed_lookup {
                if !lookup.contains(&guide.guide_id) {
                    continue;
                }
            }
            let mut notes = vec![];
            let mut spacer = guide
                .protospacer
                .as_bytes()
                .iter()
                .map(|b| match b.to_ascii_uppercase() {
                    b'U' => b'T',
                    other => other,
                } as char)
                .collect::<String>();
            if apply_5prime_g_extension && !spacer.starts_with('G') {
                spacer = format!("G{spacer}");
                notes.push("5' G extension applied".to_string());
            }
            let mut forward = format!(
                "{}{}{}",
                template.forward_prefix, spacer, template.forward_suffix
            );
            let reverse_spacer = if template.reverse_uses_reverse_complement_of_spacer {
                Self::reverse_complement(&spacer)
            } else {
                spacer
            };
            let mut reverse = format!(
                "{}{}{}",
                template.reverse_prefix, reverse_spacer, template.reverse_suffix
            );
            if template.uppercase_output {
                forward = forward.to_ascii_uppercase();
                reverse = reverse.to_ascii_uppercase();
            }
            records.push(GuideOligoRecord {
                guide_id: guide.guide_id.clone(),
                rank: guide.rank,
                forward_oligo: forward,
                reverse_oligo: reverse,
                notes,
            });
        }
        if records.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "No guides selected for oligo generation in set '{}'",
                    guide_set_id
                ),
            });
        }

        let now = Self::now_unix_ms();
        let requested = output_oligo_set_id
            .as_deref()
            .map(Self::normalize_oligo_set_id)
            .transpose()?;
        let default_id = format!("{}_{}_{}", guide_set_id, template.template_id, now);
        let oligo_set_id =
            Self::unique_oligo_set_id(&store, requested.as_deref().unwrap_or(default_id.as_str()));
        store
            .latest_oligo_set_by_guide_set
            .insert(guide_set_id.clone(), oligo_set_id.clone());
        store.oligo_sets.insert(
            oligo_set_id.clone(),
            GuideOligoSet {
                oligo_set_id: oligo_set_id.clone(),
                guide_set_id: guide_set_id.clone(),
                generated_at_unix_ms: now,
                template: template.clone(),
                apply_5prime_g_extension,
                records: records.clone(),
            },
        );
        Self::append_guide_design_audit(
            &mut store,
            "GenerateGuideOligos",
            &guide_set_id,
            json!({
                "oligo_set_id": oligo_set_id,
                "record_count": records.len(),
                "template_id": template.template_id,
                "passed_only": passed_only,
                "apply_5prime_g_extension": apply_5prime_g_extension
            }),
        );
        self.write_guide_design_store(store)?;
        result.messages.push(format!(
            "Generated oligos for guide set '{}' ({} guide(s), template '{}')",
            guide_set_id,
            records.len(),
            template.template_id
        ));
        if passed_only {
            result
                .messages
                .push("Oligo generation used only practical-filter passing guides".to_string());
        }
        Ok(())
    }

    fn csv_escape(value: &str) -> String {
        if value.contains(',') || value.contains('"') || value.contains('\n') {
            format!("\"{}\"", value.replace('"', "\"\""))
        } else {
            value.to_string()
        }
    }

    fn ensure_output_parent_dir(path: &str) -> Result<(), EngineError> {
        let parent = Path::new(path)
            .parent()
            .map(Path::to_path_buf)
            .unwrap_or_else(|| PathBuf::from("."));
        std::fs::create_dir_all(&parent).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not create output directory '{}' for export '{}': {e}",
                parent.display(),
                path
            ),
        })
    }

    fn op_export_guide_oligos(
        &mut self,
        guide_set_id: String,
        oligo_set_id: Option<String>,
        format: GuideOligoExportFormat,
        path: String,
        plate_format: Option<GuideOligoPlateFormat>,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let guide_set_id = Self::normalize_guide_set_id(&guide_set_id)?;
        if path.trim().is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "ExportGuideOligos requires non-empty path".to_string(),
            });
        }
        let mut store = self.read_guide_design_store();
        if !store.guide_sets.contains_key(&guide_set_id) {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!("Guide set '{}' not found", guide_set_id),
            });
        }
        let oligo_set =
            Self::resolve_oligo_set_for_export(&store, &guide_set_id, oligo_set_id.as_deref())?;

        let text = match format {
            GuideOligoExportFormat::CsvTable => {
                let mut rows = vec!["guide_id,rank,forward_oligo,reverse_oligo,notes".to_string()];
                for record in &oligo_set.records {
                    let rank = record.rank.map(|v| v.to_string()).unwrap_or_default();
                    let notes = record.notes.join("; ");
                    rows.push(format!(
                        "{},{},{},{},{}",
                        Self::csv_escape(&record.guide_id),
                        Self::csv_escape(&rank),
                        Self::csv_escape(&record.forward_oligo),
                        Self::csv_escape(&record.reverse_oligo),
                        Self::csv_escape(&notes),
                    ));
                }
                rows.join("\n")
            }
            GuideOligoExportFormat::PlateCsv => {
                let plate_format = plate_format.unwrap_or_default();
                let (rows_per_plate, cols_per_plate) = plate_format.dimensions();
                let capacity = rows_per_plate * cols_per_plate;
                let mut rows =
                    vec!["plate,well,guide_id,rank,forward_oligo,reverse_oligo,notes".to_string()];
                for (idx, record) in oligo_set.records.iter().enumerate() {
                    let plate_index = idx / capacity + 1;
                    let within_plate = idx % capacity;
                    let row = within_plate / cols_per_plate;
                    let col = within_plate % cols_per_plate + 1;
                    let row_char = (b'A' + row as u8) as char;
                    let well = format!("{row_char}{:02}", col);
                    let rank = record.rank.map(|v| v.to_string()).unwrap_or_default();
                    let notes = record.notes.join("; ");
                    rows.push(format!(
                        "{},{},{},{},{},{},{}",
                        plate_index,
                        Self::csv_escape(&well),
                        Self::csv_escape(&record.guide_id),
                        Self::csv_escape(&rank),
                        Self::csv_escape(&record.forward_oligo),
                        Self::csv_escape(&record.reverse_oligo),
                        Self::csv_escape(&notes),
                    ));
                }
                rows.join("\n")
            }
            GuideOligoExportFormat::Fasta => {
                let mut out = String::new();
                for record in &oligo_set.records {
                    out.push_str(&format!(
                        ">{}|forward\n{}\n",
                        record.guide_id, record.forward_oligo
                    ));
                    out.push_str(&format!(
                        ">{}|reverse\n{}\n",
                        record.guide_id, record.reverse_oligo
                    ));
                }
                out
            }
        };

        Self::ensure_output_parent_dir(&path)?;
        std::fs::write(&path, text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write oligo export '{}': {e}", path),
        })?;

        let report = GuideOligoExportReport {
            guide_set_id: guide_set_id.clone(),
            oligo_set_id: oligo_set.oligo_set_id.clone(),
            format: format.as_str().to_string(),
            path: path.clone(),
            exported_records: oligo_set.records.len(),
        };
        Self::append_guide_design_audit(
            &mut store,
            "ExportGuideOligos",
            &guide_set_id,
            serde_json::to_value(&report).unwrap_or_else(|_| json!({ "path": path })),
        );
        self.write_guide_design_store(store)?;
        result.messages.push(format!(
            "Exported {} oligo records for '{}' to '{}' as {}",
            report.exported_records, guide_set_id, report.path, report.format
        ));
        Ok(())
    }

    fn op_export_guide_protocol_text(
        &mut self,
        guide_set_id: String,
        oligo_set_id: Option<String>,
        path: String,
        include_qc_checklist: Option<bool>,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let guide_set_id = Self::normalize_guide_set_id(&guide_set_id)?;
        if path.trim().is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "ExportGuideProtocolText requires non-empty path".to_string(),
            });
        }
        let mut store = self.read_guide_design_store();
        if !store.guide_sets.contains_key(&guide_set_id) {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!("Guide set '{}' not found", guide_set_id),
            });
        }
        let oligo_set =
            Self::resolve_oligo_set_for_export(&store, &guide_set_id, oligo_set_id.as_deref())?;
        let include_qc = include_qc_checklist.unwrap_or(true);

        let mut text = String::new();
        text.push_str("GENtle Guide Oligo Protocol\n");
        text.push_str("==========================\n\n");
        text.push_str(&format!("Guide set: {}\n", guide_set_id));
        text.push_str(&format!("Oligo set: {}\n", oligo_set.oligo_set_id));
        text.push_str(&format!(
            "Template: {} ({})\n",
            oligo_set.template.template_id, oligo_set.template.description
        ));
        text.push_str(&format!(
            "Generated guides: {}\n\n",
            oligo_set.records.len()
        ));
        text.push_str("Suggested steps:\n");
        text.push_str("1. Prepare oligo stocks according to ordering sheet.\n");
        text.push_str("2. Anneal forward/reverse oligo pairs.\n");
        text.push_str("3. Clone into vector backbone using template-defined overhang strategy.\n");
        text.push_str("4. Transform, recover colonies, and verify insertion by sequencing.\n\n");
        text.push_str("Guide oligos:\n");
        for record in &oligo_set.records {
            let rank = record
                .rank
                .map(|v| format!("rank={v}"))
                .unwrap_or_else(|| "rank=-".to_string());
            text.push_str(&format!(
                "- {} ({})\n  forward: {}\n  reverse: {}\n",
                record.guide_id, rank, record.forward_oligo, record.reverse_oligo
            ));
            if !record.notes.is_empty() {
                text.push_str(&format!("  notes: {}\n", record.notes.join("; ")));
            }
        }
        if include_qc {
            text.push_str("\nQC checklist:\n");
            text.push_str("- Confirm oligo lengths and overhang sequences.\n");
            text.push_str(
                "- Confirm no guide contains forbidden motifs for your cloning strategy.\n",
            );
            text.push_str("- Confirm expected insert size by colony PCR or digest.\n");
            text.push_str("- Confirm sequence identity by Sanger/NGS.\n");
        }

        Self::ensure_output_parent_dir(&path)?;
        std::fs::write(&path, text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write protocol export '{}': {e}", path),
        })?;

        let report = GuideProtocolExportReport {
            guide_set_id: guide_set_id.clone(),
            oligo_set_id: oligo_set.oligo_set_id.clone(),
            path: path.clone(),
            guide_count: oligo_set.records.len(),
        };
        Self::append_guide_design_audit(
            &mut store,
            "ExportGuideProtocolText",
            &guide_set_id,
            serde_json::to_value(&report).unwrap_or_else(|_| json!({ "path": path })),
        );
        self.write_guide_design_store(store)?;
        result.messages.push(format!(
            "Exported guide protocol text for '{}' to '{}'",
            guide_set_id, report.path
        ));
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
        feature_geometry_mode: Option<CandidateFeatureGeometryMode>,
        feature_boundary_mode: Option<CandidateFeatureBoundaryMode>,
        feature_strand_relation: Option<CandidateFeatureStrandRelation>,
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
        let feature_geometry_mode = feature_geometry_mode.unwrap_or_default();
        let requested_boundary_mode = feature_boundary_mode.unwrap_or_default();
        let feature_strand_relation = feature_strand_relation.unwrap_or_default();
        let effective_boundary_mode =
            if feature_geometry_mode == CandidateFeatureGeometryMode::FeatureBoundaries {
                requested_boundary_mode
            } else {
                CandidateFeatureBoundaryMode::Any
            };
        if feature_geometry_mode != CandidateFeatureGeometryMode::FeatureBoundaries
            && feature_boundary_mode.is_some()
        {
            result.warnings.push(
                "feature_boundary_mode is ignored unless feature_geometry_mode=feature_boundaries"
                    .to_string(),
            );
        }

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

        let mut feature_cache: HashMap<String, Vec<FeatureDistanceTarget>> = HashMap::new();
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
            feature_cache.insert(
                seq_id.clone(),
                Self::collect_feature_distance_targets(
                    dna,
                    feature_geometry_mode,
                    effective_boundary_mode,
                ),
            );
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
                feature_strand_relation,
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
        result.messages.push(format!(
            "Distance scoring mode for '{}': geometry='{}', boundary='{}', strand_relation='{}'",
            set_name,
            feature_geometry_mode.as_str(),
            effective_boundary_mode.as_str(),
            feature_strand_relation.as_str()
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

    fn compare_candidates_by_tie_break(
        a: &CandidateRecord,
        b: &CandidateRecord,
        policy: CandidateTieBreakPolicy,
    ) -> Ordering {
        let a_len = a.end_0based.saturating_sub(a.start_0based);
        let b_len = b.end_0based.saturating_sub(b.start_0based);
        match policy {
            CandidateTieBreakPolicy::SeqStartEnd => a
                .seq_id
                .cmp(&b.seq_id)
                .then(a.start_0based.cmp(&b.start_0based))
                .then(a.end_0based.cmp(&b.end_0based))
                .then(a.sequence.cmp(&b.sequence)),
            CandidateTieBreakPolicy::SeqEndStart => a
                .seq_id
                .cmp(&b.seq_id)
                .then(a.end_0based.cmp(&b.end_0based))
                .then(a.start_0based.cmp(&b.start_0based))
                .then(a.sequence.cmp(&b.sequence)),
            CandidateTieBreakPolicy::LengthAscending => a_len
                .cmp(&b_len)
                .then(a.seq_id.cmp(&b.seq_id))
                .then(a.start_0based.cmp(&b.start_0based))
                .then(a.end_0based.cmp(&b.end_0based))
                .then(a.sequence.cmp(&b.sequence)),
            CandidateTieBreakPolicy::LengthDescending => b_len
                .cmp(&a_len)
                .then(a.seq_id.cmp(&b.seq_id))
                .then(a.start_0based.cmp(&b.start_0based))
                .then(a.end_0based.cmp(&b.end_0based))
                .then(a.sequence.cmp(&b.sequence)),
            CandidateTieBreakPolicy::SequenceLexicographic => a
                .sequence
                .cmp(&b.sequence)
                .then(a.seq_id.cmp(&b.seq_id))
                .then(a.start_0based.cmp(&b.start_0based))
                .then(a.end_0based.cmp(&b.end_0based)),
        }
    }

    fn op_score_candidate_set_weighted_objective(
        &mut self,
        set_name: String,
        metric: String,
        objectives: Vec<CandidateWeightedObjectiveTerm>,
        normalize_metrics: Option<bool>,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let set_name = Self::normalize_candidate_set_name(&set_name)?;
        let metric_name = Self::normalize_metric_name(&metric);
        if objectives.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "ScoreCandidateSetWeightedObjective requires at least one objective term"
                    .to_string(),
            });
        }
        let normalize_metrics = normalize_metrics.unwrap_or(true);

        let mut compiled = Vec::with_capacity(objectives.len());
        for term in objectives {
            let metric = Self::normalize_metric_name(&term.metric);
            if !term.weight.is_finite() || term.weight <= 0.0 {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Invalid weighted objective term for metric '{}': weight must be finite and > 0",
                        metric
                    ),
                });
            }
            compiled.push((metric, term.weight, term.direction));
        }

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

        let mut bounds = Vec::with_capacity(compiled.len());
        for (metric_name, _, _) in &compiled {
            let mut min_value = f64::INFINITY;
            let mut max_value = f64::NEG_INFINITY;
            for (idx, candidate) in set.candidates.iter().enumerate() {
                let value =
                    candidate
                        .metrics
                        .get(metric_name)
                        .copied()
                        .ok_or_else(|| EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!(
                                "Candidate {} in '{}' is missing metric '{}'",
                                idx, set_name, metric_name
                            ),
                        })?;
                if !value.is_finite() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Candidate {} in '{}' has non-finite metric '{}'",
                            idx, set_name, metric_name
                        ),
                    });
                }
                min_value = min_value.min(value);
                max_value = max_value.max(value);
            }
            bounds.push((min_value, max_value));
        }

        let mut combined_values = Vec::with_capacity(set.candidates.len());
        for candidate in &mut set.candidates {
            let mut combined = 0.0f64;
            for ((metric_name, weight, direction), (min_value, max_value)) in
                compiled.iter().zip(bounds.iter())
            {
                let raw_value =
                    candidate
                        .metrics
                        .get(metric_name)
                        .copied()
                        .ok_or_else(|| EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!(
                                "Candidate in '{}' is missing metric '{}'",
                                set_name, metric_name
                            ),
                        })?;
                let objective_value = if normalize_metrics {
                    let scaled = if *max_value > *min_value {
                        (raw_value - *min_value) / (*max_value - *min_value)
                    } else {
                        0.5
                    };
                    match direction {
                        CandidateObjectiveDirection::Maximize => scaled,
                        CandidateObjectiveDirection::Minimize => 1.0 - scaled,
                    }
                } else {
                    match direction {
                        CandidateObjectiveDirection::Maximize => raw_value,
                        CandidateObjectiveDirection::Minimize => -raw_value,
                    }
                };
                combined += *weight * objective_value;
            }
            candidate.metrics.insert(metric_name.clone(), combined);
            combined_values.push(combined);
        }

        let min_combined = combined_values
            .iter()
            .copied()
            .fold(f64::INFINITY, f64::min);
        let max_combined = combined_values
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        self.write_candidate_store(store)?;
        result.messages.push(format!(
            "Scored candidate set '{}' with weighted objective metric '{}'",
            set_name, metric_name
        ));
        result.messages.push(format!(
            "Weighted objective mode for '{}': normalize_metrics={}",
            set_name, normalize_metrics
        ));
        result.messages.push(format!(
            "Metric '{}' range in '{}': [{:.6}, {:.6}]",
            metric_name, set_name, min_combined, max_combined
        ));
        Ok(())
    }

    fn op_top_k_candidate_set(
        &mut self,
        input_set: String,
        output_set: String,
        metric: String,
        k: usize,
        direction: Option<CandidateObjectiveDirection>,
        tie_break: Option<CandidateTieBreakPolicy>,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let input_set = Self::normalize_candidate_set_name(&input_set)?;
        let output_set = Self::normalize_candidate_set_name(&output_set)?;
        let metric_name = Self::normalize_metric_name(&metric);
        if k == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "TopKCandidateSet requires k >= 1".to_string(),
            });
        }
        let direction = direction.unwrap_or_default();
        let tie_break = tie_break.unwrap_or_default();

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

        let mut scored = Vec::with_capacity(input.candidates.len());
        for (idx, candidate) in input.candidates.into_iter().enumerate() {
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
            scored.push((candidate, value));
        }
        scored.sort_by(|(left, left_value), (right, right_value)| {
            let primary = match direction {
                CandidateObjectiveDirection::Maximize => right_value
                    .partial_cmp(left_value)
                    .unwrap_or(Ordering::Equal),
                CandidateObjectiveDirection::Minimize => left_value
                    .partial_cmp(right_value)
                    .unwrap_or(Ordering::Equal),
            };
            if primary == Ordering::Equal {
                Self::compare_candidates_by_tie_break(left, right, tie_break)
            } else {
                primary
            }
        });

        let selected = scored
            .into_iter()
            .take(k)
            .map(|(candidate, _)| candidate)
            .collect::<Vec<_>>();
        let selected_count = selected.len();
        let replaced_existing = store
            .sets
            .insert(
                output_set.clone(),
                CandidateSet {
                    name: output_set.clone(),
                    created_at_unix_ms: Self::now_unix_ms(),
                    source_seq_ids: input.source_seq_ids,
                    candidates: selected,
                },
            )
            .is_some();
        self.write_candidate_store(store)?;
        result.messages.push(format!(
            "Selected top {} candidate(s) from '{}' into '{}' by metric '{}' ({}, tie_break={})",
            selected_count,
            input_set,
            output_set,
            metric_name,
            direction.as_str(),
            tie_break.as_str()
        ));
        if replaced_existing {
            result.warnings.push(format!(
                "TopKCandidateSet output '{}' replaced existing candidate set",
                output_set
            ));
        }
        Ok(())
    }

    fn candidate_dominates(
        left_values: &[f64],
        right_values: &[f64],
        objectives: &[CandidateObjectiveSpec],
    ) -> bool {
        let mut strictly_better = false;
        for (idx, objective) in objectives.iter().enumerate() {
            let left = left_values.get(idx).copied().unwrap_or(0.0);
            let right = right_values.get(idx).copied().unwrap_or(0.0);
            match objective.direction {
                CandidateObjectiveDirection::Maximize => {
                    if left < right {
                        return false;
                    }
                    if left > right {
                        strictly_better = true;
                    }
                }
                CandidateObjectiveDirection::Minimize => {
                    if left > right {
                        return false;
                    }
                    if left < right {
                        strictly_better = true;
                    }
                }
            }
        }
        strictly_better
    }

    fn op_pareto_frontier_candidate_set(
        &mut self,
        input_set: String,
        output_set: String,
        objectives: Vec<CandidateObjectiveSpec>,
        max_candidates: Option<usize>,
        tie_break: Option<CandidateTieBreakPolicy>,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let input_set = Self::normalize_candidate_set_name(&input_set)?;
        let output_set = Self::normalize_candidate_set_name(&output_set)?;
        if objectives.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "ParetoFrontierCandidateSet requires at least one objective".to_string(),
            });
        }
        if max_candidates == Some(0) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "ParetoFrontierCandidateSet max_candidates must be >= 1".to_string(),
            });
        }
        let tie_break = tie_break.unwrap_or_default();

        let mut compiled = Vec::with_capacity(objectives.len());
        for objective in objectives {
            compiled.push(CandidateObjectiveSpec {
                metric: Self::normalize_metric_name(&objective.metric),
                direction: objective.direction,
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

        let mut objective_values: Vec<Vec<f64>> = vec![];
        for (candidate_idx, candidate) in input.candidates.iter().enumerate() {
            let mut row = Vec::with_capacity(compiled.len());
            for objective in &compiled {
                let value = candidate
                    .metrics
                    .get(&objective.metric)
                    .copied()
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Candidate {} in '{}' is missing metric '{}'",
                            candidate_idx, input_set, objective.metric
                        ),
                    })?;
                if !value.is_finite() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Candidate {} in '{}' has non-finite metric '{}'",
                            candidate_idx, input_set, objective.metric
                        ),
                    });
                }
                row.push(value);
            }
            objective_values.push(row);
        }

        let mut dominated = vec![false; input.candidates.len()];
        for i in 0..input.candidates.len() {
            if dominated[i] {
                continue;
            }
            for j in 0..input.candidates.len() {
                if i == j {
                    continue;
                }
                if Self::candidate_dominates(&objective_values[j], &objective_values[i], &compiled)
                {
                    dominated[i] = true;
                    break;
                }
            }
        }

        let mut frontier = input
            .candidates
            .into_iter()
            .zip(dominated.into_iter())
            .filter_map(
                |(candidate, is_dominated)| if is_dominated { None } else { Some(candidate) },
            )
            .collect::<Vec<_>>();
        let raw_frontier_count = frontier.len();
        if let Some(limit) = max_candidates {
            if frontier.len() > limit {
                frontier.sort_by(|a, b| Self::compare_candidates_by_tie_break(a, b, tie_break));
                frontier.truncate(limit);
                result.warnings.push(format!(
                    "Pareto frontier for '{}' had {} candidates and was truncated to {} by tie-break policy '{}'",
                    input_set,
                    raw_frontier_count,
                    limit,
                    tie_break.as_str()
                ));
            }
        }

        let output_count = frontier.len();
        let replaced_existing = store
            .sets
            .insert(
                output_set.clone(),
                CandidateSet {
                    name: output_set.clone(),
                    created_at_unix_ms: Self::now_unix_ms(),
                    source_seq_ids: input.source_seq_ids,
                    candidates: frontier,
                },
            )
            .is_some();
        self.write_candidate_store(store)?;
        result.messages.push(format!(
            "Computed Pareto frontier from '{}' into '{}' ({} candidate(s), objectives={})",
            input_set,
            output_set,
            output_count,
            compiled
                .iter()
                .map(|objective| format!("{}:{}", objective.metric, objective.direction.as_str()))
                .collect::<Vec<_>>()
                .join(", ")
        ));
        if replaced_existing {
            result.warnings.push(format!(
                "ParetoFrontierCandidateSet output '{}' replaced existing candidate set",
                output_set
            ));
        }
        Ok(())
    }

    fn op_upsert_workflow_macro_template(
        &mut self,
        name: String,
        description: Option<String>,
        details_url: Option<String>,
        parameters: Vec<WorkflowMacroTemplateParam>,
        script: String,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let name = Self::normalize_workflow_macro_template_name(&name)?;
        let script = script.trim();
        if script.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Workflow macro template script cannot be empty".to_string(),
            });
        }
        let details_url = Self::normalize_workflow_macro_template_details_url(details_url)?;

        let mut normalized_parameters = Vec::with_capacity(parameters.len());
        let mut seen = HashSet::new();
        for parameter in parameters {
            let param_name = Self::normalize_workflow_macro_param_name(&parameter.name)?;
            if !seen.insert(param_name.clone()) {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Workflow macro template '{}' contains duplicate parameter '{}'",
                        name, param_name
                    ),
                });
            }
            let default_value = parameter
                .default_value
                .map(|value| value.trim().to_string())
                .filter(|value| !value.is_empty());
            let required = if default_value.is_some() {
                false
            } else {
                parameter.required
            };
            normalized_parameters.push(WorkflowMacroTemplateParam {
                name: param_name,
                default_value,
                required,
            });
        }

        let declared = normalized_parameters
            .iter()
            .map(|parameter| parameter.name.clone())
            .collect::<HashSet<_>>();
        let placeholder_regex =
            Regex::new(r"\$\{([A-Za-z_][A-Za-z0-9_]*)\}").map_err(|e| EngineError {
                code: ErrorCode::Internal,
                message: format!("Could not compile workflow macro placeholder regex: {e}"),
            })?;
        for captures in placeholder_regex.captures_iter(script) {
            if let Some(param_name) = captures.get(1).map(|m| m.as_str()) {
                if !declared.contains(param_name) {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Workflow macro template '{}' references undeclared parameter '{}' in script",
                            name, param_name
                        ),
                    });
                }
            }
        }

        let now = Self::now_unix_ms();
        let mut store = self.read_workflow_macro_template_store();
        let created_at_unix_ms = store
            .templates
            .get(&name)
            .map(|template| template.created_at_unix_ms)
            .unwrap_or(now);
        let replaced = store
            .templates
            .insert(
                name.clone(),
                WorkflowMacroTemplate {
                    name: name.clone(),
                    description: description
                        .map(|text| text.trim().to_string())
                        .filter(|text| !text.is_empty()),
                    details_url,
                    parameters: normalized_parameters,
                    template_schema: CLONING_MACRO_TEMPLATE_SCHEMA.to_string(),
                    script: script.to_string(),
                    created_at_unix_ms,
                    updated_at_unix_ms: now,
                },
            )
            .is_some();
        self.write_workflow_macro_template_store(store)?;
        if replaced {
            result
                .messages
                .push(format!("Updated workflow macro template '{}'", name));
        } else {
            result
                .messages
                .push(format!("Added workflow macro template '{}'", name));
        }
        Ok(())
    }

    fn op_delete_workflow_macro_template(
        &mut self,
        name: String,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let name = Self::normalize_workflow_macro_template_name(&name)?;
        let mut store = self.read_workflow_macro_template_store();
        let removed = store.templates.remove(&name).is_some();
        self.write_workflow_macro_template_store(store)?;
        if removed {
            result
                .messages
                .push(format!("Deleted workflow macro template '{}'", name));
        } else {
            result.warnings.push(format!(
                "Workflow macro template '{}' was not present",
                name
            ));
        }
        Ok(())
    }

    fn op_upsert_candidate_macro_template(
        &mut self,
        name: String,
        description: Option<String>,
        details_url: Option<String>,
        parameters: Vec<CandidateMacroTemplateParam>,
        script: String,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let name = Self::normalize_candidate_macro_template_name(&name)?;
        let script = script.trim();
        if script.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Candidate macro template script cannot be empty".to_string(),
            });
        }
        let details_url = Self::normalize_candidate_macro_template_details_url(details_url)?;

        let mut normalized_parameters = Vec::with_capacity(parameters.len());
        let mut seen = HashSet::new();
        for parameter in parameters {
            let param_name = Self::normalize_candidate_macro_param_name(&parameter.name)?;
            if !seen.insert(param_name.clone()) {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Candidate macro template '{}' contains duplicate parameter '{}'",
                        name, param_name
                    ),
                });
            }
            let default_value = parameter
                .default_value
                .map(|value| value.trim().to_string())
                .filter(|value| !value.is_empty());
            let required = if default_value.is_some() {
                false
            } else {
                parameter.required
            };
            normalized_parameters.push(CandidateMacroTemplateParam {
                name: param_name,
                default_value,
                required,
            });
        }

        let declared = normalized_parameters
            .iter()
            .map(|parameter| parameter.name.clone())
            .collect::<HashSet<_>>();
        let placeholder_regex =
            Regex::new(r"\$\{([A-Za-z_][A-Za-z0-9_]*)\}").map_err(|e| EngineError {
                code: ErrorCode::Internal,
                message: format!("Could not compile candidate macro placeholder regex: {e}"),
            })?;
        for captures in placeholder_regex.captures_iter(script) {
            if let Some(param_name) = captures.get(1).map(|m| m.as_str()) {
                if !declared.contains(param_name) {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Candidate macro template '{}' references undeclared parameter '{}' in script",
                            name, param_name
                        ),
                    });
                }
            }
        }

        let now = Self::now_unix_ms();
        let mut store = self.read_candidate_macro_template_store();
        let created_at_unix_ms = store
            .templates
            .get(&name)
            .map(|template| template.created_at_unix_ms)
            .unwrap_or(now);
        let replaced = store
            .templates
            .insert(
                name.clone(),
                CandidateMacroTemplate {
                    name: name.clone(),
                    description: description
                        .map(|text| text.trim().to_string())
                        .filter(|text| !text.is_empty()),
                    details_url,
                    parameters: normalized_parameters,
                    script: script.to_string(),
                    created_at_unix_ms,
                    updated_at_unix_ms: now,
                },
            )
            .is_some();
        self.write_candidate_macro_template_store(store)?;
        if replaced {
            result
                .messages
                .push(format!("Updated candidate macro template '{}'", name));
        } else {
            result
                .messages
                .push(format!("Added candidate macro template '{}'", name));
        }
        Ok(())
    }

    fn op_delete_candidate_macro_template(
        &mut self,
        name: String,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let name = Self::normalize_candidate_macro_template_name(&name)?;
        let mut store = self.read_candidate_macro_template_store();
        let removed = store.templates.remove(&name).is_some();
        self.write_candidate_macro_template_store(store)?;
        if removed {
            result
                .messages
                .push(format!("Deleted candidate macro template '{}'", name));
        } else {
            result.warnings.push(format!(
                "Candidate macro template '{}' was not present",
                name
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
            let catalog_path = entry
                .get("catalog_path")
                .and_then(|v| v.as_str())
                .map(str::trim)
                .filter(|v| !v.is_empty())
                .map(|v| v.to_string());
            let cache_dir = entry
                .get("cache_dir")
                .and_then(|v| v.as_str())
                .map(str::trim)
                .filter(|v| !v.is_empty())
                .map(|v| v.to_string());
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
                    catalog_path,
                    cache_dir,
                });
            }
        }

        latest.ok_or_else(|| EngineError {
            code: ErrorCode::NotFound,
            message: format!(
                "Sequence '{seq_id}' is not anchored to a genome interval; run ExtractGenomeRegion/ExtractGenomeGene/ExtendGenomeAnchor first"
            ),
        })
    }

    fn normalize_chromosome_alias(raw: &str) -> String {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return String::new();
        }
        let lower = trimmed.to_ascii_lowercase();
        let core = lower.strip_prefix("chr").unwrap_or(&lower);
        if core.is_empty() {
            return String::new();
        }
        if matches!(core, "m" | "mt" | "mitochondria" | "mitochondrion") {
            return "mt".to_string();
        }
        if core.chars().all(|ch| ch.is_ascii_digit()) {
            let normalized = core.trim_start_matches('0');
            if normalized.is_empty() {
                return "0".to_string();
            }
            return normalized.to_string();
        }
        if core.len() == 1 && matches!(core, "x" | "y" | "w" | "z") {
            return core.to_ascii_uppercase();
        }
        core.to_string()
    }

    fn chromosomes_match(left: &str, right: &str) -> bool {
        Self::normalize_chromosome_alias(left) == Self::normalize_chromosome_alias(right)
    }

    fn append_chromosome_mismatch_warning(
        report: &mut GenomeBedTrackImportReport,
        anchor_chromosome: &str,
        source_label: &str,
        mismatch_counts: &HashMap<String, usize>,
    ) {
        if mismatch_counts.is_empty() {
            return;
        }
        let mut sorted = mismatch_counts
            .iter()
            .map(|(chrom, count)| (chrom.clone(), *count))
            .collect::<Vec<_>>();
        sorted.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
        report.skipped_wrong_chromosome_examples = sorted
            .iter()
            .take(6)
            .map(|(chrom, count)| format!("{chrom} ({count})"))
            .collect();
        let seen = sorted
            .iter()
            .take(3)
            .map(|(chrom, count)| format!("{chrom} ({count})"))
            .collect::<Vec<_>>()
            .join(", ");
        report.warnings.push(format!(
            "{} record(s) in {} input did not match anchor chromosome '{}' (examples: {})",
            report.skipped_wrong_chromosome, source_label, anchor_chromosome, seen
        ));
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
        let format = fields
            .get(8)
            .map(|v| v.trim())
            .filter(|v| !v.is_empty() && *v != ".")
            .map(|v| v.to_string());
        let sample_columns = if fields.len() > 9 {
            fields[9..]
                .iter()
                .map(|value| value.trim().to_string())
                .collect::<Vec<_>>()
        } else {
            vec![]
        };

        Ok(VcfRecord {
            chromosome: chromosome.to_string(),
            pos_1based,
            id,
            reference,
            alternates,
            qual,
            filter,
            info,
            format,
            sample_columns,
        })
    }

    fn classify_vcf_alt(reference: &str, alt: &str) -> VcfVariantClass {
        let ref_u = reference.trim().to_ascii_uppercase();
        let alt_u = alt.trim().to_ascii_uppercase();
        if alt_u.is_empty() || alt_u == "." {
            return VcfVariantClass::Other;
        }
        if alt_u.starts_with('<')
            || alt_u.ends_with('>')
            || alt_u.contains('[')
            || alt_u.contains(']')
            || alt_u == "*"
        {
            return VcfVariantClass::Sv;
        }
        let ref_len = ref_u.len().max(1);
        let alt_len = alt_u.len().max(1);
        if ref_len == 1 && alt_len == 1 {
            return VcfVariantClass::Snp;
        }
        if alt_len > ref_len {
            return VcfVariantClass::Ins;
        }
        if alt_len < ref_len {
            return VcfVariantClass::Del;
        }
        VcfVariantClass::Other
    }

    fn summarize_vcf_alt_genotype(
        record: &VcfRecord,
        alt_allele_index_1based: usize,
        sample_names: &[String],
    ) -> Option<VcfAltGenotypeSummary> {
        let format = record.format.as_deref()?;
        let format_fields = format.split(':').collect::<Vec<_>>();
        let gt_idx = format_fields
            .iter()
            .position(|field| field.eq_ignore_ascii_case("GT"))?;
        let mut summary = VcfAltGenotypeSummary::default();
        for (sample_idx, sample_column) in record.sample_columns.iter().enumerate() {
            let value_fields = sample_column.split(':').collect::<Vec<_>>();
            let Some(raw_gt) = value_fields.get(gt_idx).map(|v| v.trim()) else {
                continue;
            };
            if raw_gt.is_empty() || raw_gt == "." {
                continue;
            }
            let (tokens, phased) = if raw_gt.contains('|') {
                (raw_gt.split('|').collect::<Vec<_>>(), true)
            } else if raw_gt.contains('/') {
                (raw_gt.split('/').collect::<Vec<_>>(), false)
            } else {
                (vec![raw_gt], false)
            };
            let mut parsed = vec![];
            for token in tokens {
                let trimmed = token.trim();
                if trimmed.is_empty() || trimmed == "." {
                    continue;
                }
                if let Ok(index) = trimmed.parse::<usize>() {
                    parsed.push(index);
                }
            }
            if parsed.is_empty() {
                continue;
            }
            if !parsed.iter().any(|index| *index == alt_allele_index_1based) {
                continue;
            }
            summary.carriers += 1;
            if phased {
                summary.phased_carriers += 1;
            } else {
                summary.unphased_carriers += 1;
            }
            let has_ref = parsed.contains(&0);
            let has_other_alt = parsed
                .iter()
                .any(|idx| *idx > 0 && *idx != alt_allele_index_1based);
            let all_target_alt = parsed.iter().all(|idx| *idx == alt_allele_index_1based);
            if parsed.len() == 1 && all_target_alt {
                summary.haploid_alt += 1;
            } else if all_target_alt {
                summary.hom_alt += 1;
            } else if has_ref {
                summary.het += 1;
            } else if has_other_alt {
                summary.mixed_alt += 1;
            } else {
                summary.het += 1;
            }
            let sample_label = sample_names
                .get(sample_idx)
                .filter(|v| !v.trim().is_empty())
                .cloned()
                .unwrap_or_else(|| format!("sample_{}", sample_idx + 1));
            summary.carrier_samples.push(sample_label);
        }
        Some(summary)
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
        alt_allele_index_1based: usize,
        sample_names: &[String],
        track_name: &str,
        path: &str,
        local_start_0based: usize,
        local_end_0based_exclusive: usize,
    ) -> gb_io::seq::Feature {
        let variant_class = Self::classify_vcf_alt(&record.reference, alt);
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
            (
                "vcf_alt_allele_index".into(),
                Some(alt_allele_index_1based.to_string()),
            ),
            (
                "vcf_variant_class".into(),
                Some(variant_class.as_str().to_string()),
            ),
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
        if let Some(format) = &record.format {
            qualifiers.push(("vcf_format".into(), Some(format.clone())));
        }
        if !record.sample_columns.is_empty() {
            qualifiers.push((
                "vcf_sample_count".into(),
                Some(record.sample_columns.len().to_string()),
            ));
        }
        if let Some(genotype) =
            Self::summarize_vcf_alt_genotype(record, alt_allele_index_1based, sample_names)
        {
            qualifiers.push((
                "vcf_alt_carriers".into(),
                Some(genotype.carriers.to_string()),
            ));
            qualifiers.push((
                "vcf_alt_carrier_phased".into(),
                Some(genotype.phased_carriers.to_string()),
            ));
            qualifiers.push((
                "vcf_alt_carrier_unphased".into(),
                Some(genotype.unphased_carriers.to_string()),
            ));
            qualifiers.push(("vcf_gt_het".into(), Some(genotype.het.to_string())));
            qualifiers.push(("vcf_gt_hom_alt".into(), Some(genotype.hom_alt.to_string())));
            qualifiers.push((
                "vcf_gt_mixed_alt".into(),
                Some(genotype.mixed_alt.to_string()),
            ));
            qualifiers.push((
                "vcf_gt_haploid_alt".into(),
                Some(genotype.haploid_alt.to_string()),
            ));
            let zygosity = if genotype.hom_alt > 0
                && genotype.het == 0
                && genotype.mixed_alt == 0
                && genotype.haploid_alt == 0
            {
                "hom_alt"
            } else if genotype.het > 0
                && genotype.hom_alt == 0
                && genotype.mixed_alt == 0
                && genotype.haploid_alt == 0
            {
                "het"
            } else if genotype.haploid_alt > 0
                && genotype.hom_alt == 0
                && genotype.het == 0
                && genotype.mixed_alt == 0
            {
                "haploid_alt"
            } else if genotype.carriers > 0 {
                "mixed"
            } else {
                "none"
            };
            qualifiers.push(("vcf_zygosity".into(), Some(zygosity.to_string())));
            let phase = if genotype.phased_carriers > 0 && genotype.unphased_carriers == 0 {
                "phased"
            } else if genotype.unphased_carriers > 0 && genotype.phased_carriers == 0 {
                "unphased"
            } else if genotype.phased_carriers > 0 && genotype.unphased_carriers > 0 {
                "mixed"
            } else {
                "unknown"
            };
            qualifiers.push(("vcf_phase".into(), Some(phase.to_string())));
            if !genotype.carrier_samples.is_empty() {
                qualifiers.push((
                    "vcf_alt_carrier_samples".into(),
                    Some(
                        genotype
                            .carrier_samples
                            .iter()
                            .take(20)
                            .cloned()
                            .collect::<Vec<_>>()
                            .join(","),
                    ),
                ));
            }
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
        mut progress_cb: Option<&mut dyn FnMut(usize, usize, usize, bool) -> bool>,
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
        let progress_stride = 250usize;
        let mut reader = Self::open_text_reader(path)?;
        let mut line = String::new();
        let mut line_no = 0usize;
        let mut mismatch_counts: HashMap<String, usize> = HashMap::new();

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
            if report.parsed_records % progress_stride == 0 {
                if let Some(cb) = progress_cb.as_mut() {
                    let should_continue = (**cb)(
                        report.parsed_records,
                        report.imported_features,
                        report.skipped_records,
                        false,
                    );
                    if !should_continue {
                        report.cancelled = true;
                        break;
                    }
                }
            }

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
                *mismatch_counts
                    .entry(record.chromosome.clone())
                    .or_insert(0) += 1;
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
        Self::append_chromosome_mismatch_warning(
            &mut report,
            &anchor.chromosome,
            "BED",
            &mismatch_counts,
        );
        if report.cancelled {
            report.warnings.push(format!(
                "BED import cancelled after parsed={}, imported={}, skipped={}",
                report.parsed_records, report.imported_features, report.skipped_records
            ));
        }
        if let Some(cb) = progress_cb.as_mut() {
            let _ = (**cb)(
                report.parsed_records,
                report.imported_features,
                report.skipped_records,
                true,
            );
        }

        Ok(report)
    }

    fn resolve_bigwig_to_bedgraph_executable() -> String {
        crate::tool_overrides::resolve_tool_executable(
            BIGWIG_TO_BEDGRAPH_ENV_BIN,
            DEFAULT_BIGWIG_TO_BEDGRAPH_BIN,
        )
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
        mut progress_cb: Option<&mut dyn FnMut(usize, usize, usize, bool) -> bool>,
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
        let progress_stride = 250usize;
        let bedgraph_file = Self::convert_bigwig_to_bedgraph(path)?;
        let bedgraph_path = bedgraph_file.path().to_string_lossy().to_string();
        let mut reader = Self::open_text_reader(&bedgraph_path)?;
        let mut line = String::new();
        let mut line_no = 0usize;
        let mut mismatch_counts: HashMap<String, usize> = HashMap::new();

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
            if report.parsed_records % progress_stride == 0 {
                if let Some(cb) = progress_cb.as_mut() {
                    let should_continue = (**cb)(
                        report.parsed_records,
                        report.imported_features,
                        report.skipped_records,
                        false,
                    );
                    if !should_continue {
                        report.cancelled = true;
                        break;
                    }
                }
            }

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
                *mismatch_counts
                    .entry(record.chromosome.clone())
                    .or_insert(0) += 1;
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
        Self::append_chromosome_mismatch_warning(
            &mut report,
            &anchor.chromosome,
            "BigWig",
            &mismatch_counts,
        );
        if report.cancelled {
            report.warnings.push(format!(
                "BigWig import cancelled after parsed={}, imported={}, skipped={}",
                report.parsed_records, report.imported_features, report.skipped_records
            ));
        }
        if let Some(cb) = progress_cb.as_mut() {
            let _ = (**cb)(
                report.parsed_records,
                report.imported_features,
                report.skipped_records,
                true,
            );
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
        mut progress_cb: Option<&mut dyn FnMut(usize, usize, usize, bool) -> bool>,
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
        let progress_stride = 250usize;
        let mut reader = Self::open_text_reader(path)?;
        let mut line = String::new();
        let mut line_no = 0usize;
        let mut sample_names: Vec<String> = vec![];
        let mut mismatch_counts: HashMap<String, usize> = HashMap::new();

        while {
            line.clear();
            reader.read_line(&mut line).map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not read VCF file '{path}': {e}"),
            })? > 0
        } {
            line_no += 1;
            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }
            if trimmed.starts_with("##") {
                continue;
            }
            if trimmed.starts_with("#CHROM") {
                let fields = trimmed.split('\t').collect::<Vec<_>>();
                sample_names = if fields.len() > 9 {
                    fields[9..]
                        .iter()
                        .map(|value| value.trim().to_string())
                        .collect::<Vec<_>>()
                } else {
                    vec![]
                };
                continue;
            }
            if trimmed.starts_with('#') {
                continue;
            }
            report.parsed_records += 1;
            if report.parsed_records % progress_stride == 0 {
                if let Some(cb) = progress_cb.as_mut() {
                    let should_continue = (**cb)(
                        report.parsed_records,
                        report.imported_features,
                        report.skipped_records,
                        false,
                    );
                    if !should_continue {
                        report.cancelled = true;
                        break;
                    }
                }
            }

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
                *mismatch_counts
                    .entry(record.chromosome.clone())
                    .or_insert(0) += 1;
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
            for (alt_idx, alt) in record.alternates.iter().enumerate() {
                if report.imported_features >= MAX_IMPORTED_SIGNAL_FEATURES {
                    report.truncated_at_limit = true;
                    stop = true;
                    break;
                }
                let feature = Self::build_genome_vcf_feature(
                    &record,
                    alt,
                    alt_idx + 1,
                    &sample_names,
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
        Self::append_chromosome_mismatch_warning(
            &mut report,
            &anchor.chromosome,
            "VCF",
            &mismatch_counts,
        );
        if report.cancelled {
            report.warnings.push(format!(
                "VCF import cancelled after parsed={}, imported={}, skipped={}",
                report.parsed_records, report.imported_features, report.skipped_records
            ));
        }
        if let Some(cb) = progress_cb.as_mut() {
            let _ = (**cb)(
                report.parsed_records,
                report.imported_features,
                report.skipped_records,
                true,
            );
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
        if block.is_empty() { None } else { Some(block) }
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
        if token.is_empty() { None } else { Some(token) }
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
            catalog_path: None,
            cache_dir: None,
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

    fn next_arrangement_id(&mut self) -> String {
        loop {
            self.state.container_state.next_arrangement_counter += 1;
            let id = format!(
                "arrangement-{}",
                self.state.container_state.next_arrangement_counter
            );
            if !self.state.container_state.arrangements.contains_key(&id) {
                return id;
            }
        }
    }

    fn add_serial_arrangement(
        &mut self,
        container_ids: &[ContainerId],
        arrangement_id: Option<String>,
        name: Option<String>,
        ladders: Option<Vec<String>>,
        created_by_op: Option<&str>,
    ) -> Result<String, EngineError> {
        if container_ids.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "CreateArrangementSerial requires at least one container id".to_string(),
            });
        }
        let mut lane_container_ids: Vec<ContainerId> = vec![];
        let mut seen: HashSet<String> = HashSet::new();
        for container_id in container_ids {
            let trimmed = container_id.trim();
            if trimmed.is_empty() {
                continue;
            }
            if !self.state.container_state.containers.contains_key(trimmed) {
                return Err(EngineError {
                    code: ErrorCode::NotFound,
                    message: format!("Container '{trimmed}' not found"),
                });
            }
            if seen.insert(trimmed.to_string()) {
                lane_container_ids.push(trimmed.to_string());
            }
        }
        if lane_container_ids.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "CreateArrangementSerial requires at least one non-empty container id"
                    .to_string(),
            });
        }
        let arrangement_id = arrangement_id
            .map(|v| v.trim().to_string())
            .filter(|v| !v.is_empty())
            .unwrap_or_else(|| self.next_arrangement_id());
        if self
            .state
            .container_state
            .arrangements
            .contains_key(&arrangement_id)
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Arrangement '{arrangement_id}' already exists"),
            });
        }
        let arrangement = Arrangement {
            arrangement_id: arrangement_id.clone(),
            mode: ArrangementMode::Serial,
            name: name.map(|v| v.trim().to_string()).filter(|v| !v.is_empty()),
            lane_container_ids,
            ladders: ladders
                .unwrap_or_default()
                .into_iter()
                .map(|v| v.trim().to_string())
                .filter(|v| !v.is_empty())
                .collect::<Vec<_>>(),
            created_by_op: created_by_op.map(ToString::to_string),
            created_at_unix_ms: Self::now_unix_ms(),
        };
        self.state
            .container_state
            .arrangements
            .insert(arrangement_id.clone(), arrangement);
        Ok(arrangement_id)
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
            Operation::ExtendGenomeAnchor { .. } => Some("Extended genome anchor".to_string()),
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

    fn gel_samples_from_container_ids(
        &self,
        container_ids: &[ContainerId],
    ) -> Result<Vec<GelSampleInput>, EngineError> {
        if container_ids.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "At least one container id is required for gel rendering".to_string(),
            });
        }
        let mut samples: Vec<GelSampleInput> = vec![];
        for container_id in container_ids {
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
            let mut members: Vec<(String, usize)> = Vec::with_capacity(container.members.len());
            for seq_id in &container.members {
                let dna = self
                    .state
                    .sequences
                    .get(seq_id)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!(
                            "Container '{container_id}' references unknown sequence '{seq_id}'"
                        ),
                    })?;
                members.push((seq_id.clone(), dna.len()));
            }
            let lane_name = container
                .name
                .as_deref()
                .map(str::trim)
                .filter(|v| !v.is_empty())
                .map(ToString::to_string)
                .unwrap_or_else(|| container.container_id.clone());
            samples.push(GelSampleInput {
                name: lane_name,
                members,
            });
        }
        Ok(samples)
    }

    fn gel_samples_from_arrangement(
        &self,
        arrangement_id: &str,
    ) -> Result<(Vec<GelSampleInput>, Vec<String>), EngineError> {
        let arrangement = self
            .state
            .container_state
            .arrangements
            .get(arrangement_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Arrangement '{arrangement_id}' not found"),
            })?;
        if arrangement.mode != ArrangementMode::Serial {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Arrangement '{}' is mode '{:?}', only serial arrangements can render gels",
                    arrangement_id, arrangement.mode
                ),
            });
        }
        let samples = self.gel_samples_from_container_ids(&arrangement.lane_container_ids)?;
        Ok((samples, arrangement.ladders.clone()))
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

        let mut arrangements: Vec<EngineArrangementSummary> = self
            .state
            .container_state
            .arrangements
            .iter()
            .map(|(id, arrangement)| EngineArrangementSummary {
                id: id.to_string(),
                mode: format!("{:?}", arrangement.mode),
                lane_count: arrangement.lane_container_ids.len(),
                lane_container_ids: arrangement.lane_container_ids.clone(),
                ladders: arrangement.ladders.clone(),
            })
            .collect();
        arrangements.sort_by(|a, b| a.id.cmp(&b.id));

        EngineStateSummary {
            sequence_count: sequences.len(),
            sequences,
            container_count: containers.len(),
            containers,
            arrangement_count: arrangements.len(),
            arrangements,
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

    fn resolve_sequence_anchor_position(
        dna: &DNAsequence,
        anchor: &SequenceAnchor,
        anchor_name: &str,
    ) -> Result<usize, EngineError> {
        match anchor {
            SequenceAnchor::Position { zero_based } => {
                if *zero_based > dna.len() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Sequence anchor '{anchor_name}' position {} is out of bounds for sequence length {}",
                            zero_based,
                            dna.len()
                        ),
                    });
                }
                Ok(*zero_based)
            }
            SequenceAnchor::FeatureBoundary {
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
                    let start = from as usize;
                    let end = to as usize;
                    let pos = match boundary {
                        AnchorBoundary::Start => start,
                        AnchorBoundary::End => end,
                        AnchorBoundary::Middle => start + (end.saturating_sub(start) / 2),
                    };
                    if pos <= dna.len() {
                        matches.push(pos);
                    }
                }
                if matches.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("No feature matched sequence anchor '{}'", anchor_name),
                    });
                }
                matches.sort_unstable();
                let idx = occurrence.unwrap_or(0);
                matches.get(idx).cloned().ok_or_else(|| EngineError {
                    code: ErrorCode::NotFound,
                    message: format!(
                        "Sequence anchor '{}' occurrence {} was requested, but only {} match(es) found",
                        anchor_name, idx, matches.len()
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
        on_progress: &mut dyn FnMut(OperationProgress) -> bool,
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
            Operation::RenderFeatureExpertSvg {
                seq_id,
                target,
                path,
            } => {
                self.render_feature_expert_svg_to_path(&seq_id, &target, &path)?;
                result.messages.push(format!(
                    "Wrote feature expert SVG for '{}' target={} to '{}'",
                    seq_id,
                    target.describe(),
                    path
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
            Operation::CreateArrangementSerial {
                container_ids,
                arrangement_id,
                name,
                ladders,
            } => {
                let created_id = self.add_serial_arrangement(
                    &container_ids,
                    arrangement_id,
                    name,
                    ladders,
                    Some(&result.op_id),
                )?;
                result.messages.push(format!(
                    "Created serial arrangement '{}' with {} lane container(s)",
                    created_id,
                    self.state
                        .container_state
                        .arrangements
                        .get(&created_id)
                        .map(|a| a.lane_container_ids.len())
                        .unwrap_or(0)
                ));
            }
            Operation::RenderPoolGelSvg {
                inputs,
                path,
                ladders,
                container_ids,
                arrangement_id,
            } => {
                let mut ladder_names = ladders
                    .unwrap_or_default()
                    .into_iter()
                    .map(|v| v.trim().to_string())
                    .filter(|v| !v.is_empty())
                    .collect::<Vec<_>>();
                let samples: Vec<GelSampleInput> = if let Some(arrangement_id) =
                    arrangement_id.as_deref().map(str::trim)
                {
                    if arrangement_id.is_empty() {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: "arrangement_id cannot be empty".to_string(),
                        });
                    }
                    let (arrangement_samples, arrangement_ladders) =
                        self.gel_samples_from_arrangement(arrangement_id)?;
                    if ladder_names.is_empty() {
                        ladder_names = arrangement_ladders;
                    }
                    if !inputs.is_empty() {
                        result.warnings.push(
                            "RenderPoolGelSvg ignored 'inputs' because arrangement_id was provided"
                                .to_string(),
                        );
                    }
                    if container_ids.as_ref().is_some_and(|ids| !ids.is_empty()) {
                        result.warnings.push(
                                "RenderPoolGelSvg ignored 'container_ids' because arrangement_id was provided"
                                    .to_string(),
                            );
                    }
                    arrangement_samples
                } else if let Some(container_ids) = container_ids {
                    if container_ids.is_empty() {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: "container_ids was provided but empty".to_string(),
                        });
                    }
                    if !inputs.is_empty() {
                        result.warnings.push(
                            "RenderPoolGelSvg ignored 'inputs' because container_ids were provided"
                                .to_string(),
                        );
                    }
                    self.gel_samples_from_container_ids(&container_ids)?
                } else {
                    if inputs.is_empty() {
                        return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "RenderPoolGelSvg requires either inputs, container_ids, or arrangement_id"
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
                    vec![GelSampleInput {
                        name: format!("Input tube (n={})", members.len()),
                        members,
                    }]
                };
                let layout =
                    build_serial_gel_layout(&samples, &ladder_names).map_err(|e| EngineError {
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
                    "Wrote serial gel SVG for {} sample lane(s), {} sequence(s) to '{}' (ladders: {})",
                    layout.sample_count,
                    layout.pool_member_count,
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
                timeout_seconds,
            } => {
                let report = Self::prepare_reference_genome_once(
                    &genome_id,
                    catalog_path.as_deref(),
                    cache_dir.as_deref(),
                    timeout_seconds,
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
            Operation::ExtendGenomeAnchor {
                seq_id,
                side,
                length_bp,
                output_id,
                catalog_path,
                cache_dir,
            } => {
                if length_bp == 0 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "ExtendGenomeAnchor requires length_bp >= 1".to_string(),
                    });
                }
                if !self.state.sequences.contains_key(&seq_id) {
                    return Err(EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{seq_id}' not found"),
                    });
                }
                let anchor = self.latest_genome_anchor_for_seq(&seq_id)?;
                let resolved_catalog_path = catalog_path
                    .or(anchor.catalog_path.clone())
                    .unwrap_or_else(|| DEFAULT_GENOME_CATALOG_PATH.to_string());
                let resolved_cache_dir = cache_dir.or(anchor.cache_dir.clone());
                let catalog =
                    GenomeCatalog::from_json_file(&resolved_catalog_path).map_err(|e| {
                        EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!(
                                "Could not open genome catalog '{}': {}",
                                resolved_catalog_path, e
                            ),
                        }
                    })?;

                let anchor_is_reverse = anchor.strand == Some('-');
                let (new_start_1based, new_end_1based) = match (anchor_is_reverse, side) {
                    (false, GenomeAnchorSide::FivePrime) => (
                        anchor.start_1based.saturating_sub(length_bp).max(1),
                        anchor.end_1based,
                    ),
                    (false, GenomeAnchorSide::ThreePrime) => (
                        anchor.start_1based,
                        anchor.end_1based.saturating_add(length_bp),
                    ),
                    (true, GenomeAnchorSide::FivePrime) => (
                        anchor.start_1based,
                        anchor.end_1based.saturating_add(length_bp),
                    ),
                    (true, GenomeAnchorSide::ThreePrime) => (
                        anchor.start_1based.saturating_sub(length_bp).max(1),
                        anchor.end_1based,
                    ),
                };
                let mut sequence = catalog
                    .get_sequence_region_with_cache(
                        &anchor.genome_id,
                        &anchor.chromosome,
                        new_start_1based,
                        new_end_1based,
                        resolved_cache_dir.as_deref(),
                    )
                    .map_err(|e| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!(
                            "Could not load extended genome region {}:{}-{} from '{}': {}",
                            anchor.chromosome,
                            new_start_1based,
                            new_end_1based,
                            anchor.genome_id,
                            e
                        ),
                    })?;
                if anchor_is_reverse {
                    sequence = Self::reverse_complement(&sequence);
                }

                let side_token = side.as_str();
                let default_id = format!("{seq_id}_ext_{side_token}_{length_bp}");
                let base = output_id.unwrap_or(default_id);
                let mut dna = DNAsequence::from_sequence(&sequence).map_err(|e| EngineError {
                    code: ErrorCode::Internal,
                    message: format!("Could not construct DNA sequence from extended anchor: {e}"),
                })?;
                Self::prepare_sequence(&mut dna);
                let extended_seq_id = self.unique_seq_id(&base);
                self.state.sequences.insert(extended_seq_id.clone(), dna);
                self.add_lineage_node(
                    &extended_seq_id,
                    SequenceOrigin::Derived,
                    Some(&result.op_id),
                );
                result.created_seq_ids.push(extended_seq_id.clone());
                parent_seq_ids.push(seq_id.clone());

                let source_plan = catalog
                    .source_plan(&anchor.genome_id, resolved_cache_dir.as_deref())
                    .ok();
                let inspection = catalog
                    .inspect_prepared_genome(&anchor.genome_id, resolved_cache_dir.as_deref())
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
                    seq_id: extended_seq_id.clone(),
                    recorded_at_unix_ms: Self::now_unix_ms(),
                    operation: "ExtendGenomeAnchor".to_string(),
                    genome_id: anchor.genome_id.clone(),
                    catalog_path: resolved_catalog_path.clone(),
                    cache_dir: resolved_cache_dir.clone(),
                    chromosome: Some(anchor.chromosome.clone()),
                    start_1based: Some(new_start_1based),
                    end_1based: Some(new_end_1based),
                    gene_query: None,
                    occurrence: None,
                    gene_id: None,
                    gene_name: None,
                    strand: anchor.strand,
                    anchor_strand: Some(anchor.strand.unwrap_or('+')),
                    sequence_source_type,
                    annotation_source_type,
                    sequence_source,
                    annotation_source,
                    sequence_sha1,
                    annotation_sha1,
                });
                let side_label = match side {
                    GenomeAnchorSide::FivePrime => "5'",
                    GenomeAnchorSide::ThreePrime => "3'",
                };
                let anchor_strand = anchor.strand.unwrap_or('+');
                result.messages.push(format!(
                    "Extended genome anchor '{}' on {} by {} bp (anchor strand {}) => {}:{}-{} as '{}'",
                    seq_id,
                    side_label,
                    length_bp,
                    anchor_strand,
                    anchor.chromosome,
                    new_start_1based,
                    new_end_1based,
                    extended_seq_id
                ));
                let lower_bound_clipped = matches!(
                    (anchor_is_reverse, side),
                    (false, GenomeAnchorSide::FivePrime) | (true, GenomeAnchorSide::ThreePrime)
                ) && new_start_1based == 1
                    && anchor.start_1based <= length_bp;
                if lower_bound_clipped {
                    result.warnings.push(format!(
                        "Requested {} bp {} extension for '{}' clipped at chromosome start position 1",
                        length_bp, side_label, seq_id
                    ));
                }
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
                let seq_id_for_progress = seq_id.clone();
                let path_for_progress = path.clone();
                let mut progress_cb = |parsed_records: usize,
                                       imported_features: usize,
                                       skipped_records: usize,
                                       done: bool| {
                    on_progress(OperationProgress::GenomeTrackImport(
                        GenomeTrackImportProgress {
                            seq_id: seq_id_for_progress.clone(),
                            source: "BED".to_string(),
                            path: path_for_progress.clone(),
                            parsed_records,
                            imported_features,
                            skipped_records,
                            done,
                        },
                    ))
                };

                let report = Self::import_genome_bed_track(
                    dna,
                    &anchor,
                    &path,
                    track_name.as_deref(),
                    min_score,
                    max_score,
                    clear_existing.unwrap_or(false),
                    Some(&mut progress_cb),
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
                let seq_id_for_progress = seq_id.clone();
                let path_for_progress = path.clone();
                let mut progress_cb = |parsed_records: usize,
                                       imported_features: usize,
                                       skipped_records: usize,
                                       done: bool| {
                    on_progress(OperationProgress::GenomeTrackImport(
                        GenomeTrackImportProgress {
                            seq_id: seq_id_for_progress.clone(),
                            source: "BigWig".to_string(),
                            path: path_for_progress.clone(),
                            parsed_records,
                            imported_features,
                            skipped_records,
                            done,
                        },
                    ))
                };

                let report = Self::import_genome_bigwig_track(
                    dna,
                    &anchor,
                    &path,
                    track_name.as_deref(),
                    min_score,
                    max_score,
                    clear_existing.unwrap_or(false),
                    Some(&mut progress_cb),
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
                let seq_id_for_progress = seq_id.clone();
                let path_for_progress = path.clone();
                let mut progress_cb = |parsed_records: usize,
                                       imported_features: usize,
                                       skipped_records: usize,
                                       done: bool| {
                    on_progress(OperationProgress::GenomeTrackImport(
                        GenomeTrackImportProgress {
                            seq_id: seq_id_for_progress.clone(),
                            source: "VCF".to_string(),
                            path: path_for_progress.clone(),
                            parsed_records,
                            imported_features,
                            skipped_records,
                            done,
                        },
                    ))
                };

                let report = Self::import_genome_vcf_track(
                    dna,
                    &anchor,
                    &path,
                    track_name.as_deref(),
                    min_score,
                    max_score,
                    clear_existing.unwrap_or(false),
                    Some(&mut progress_cb),
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
                let mut out = dna
                    .extract_region_preserving_features(from, to)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Could not extract region {}..{} from sequence '{}'",
                            from, to, input
                        ),
                    })?;
                if out.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Could not extract region {}..{} from sequence '{}'",
                            from, to, input
                        ),
                    });
                }
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

                let anchor_pos = Self::resolve_sequence_anchor_position(&dna, &anchor, "anchor")?;
                let min_len = target_length_bp.saturating_sub(length_tolerance_bp).max(1);
                let max_len = target_length_bp
                    .saturating_add(length_tolerance_bp)
                    .max(min_len);

                let forward_primer = match forward_primer {
                    Some(p) => {
                        let v = Self::normalize_iupac_text(&p)?;
                        if v.is_empty() { None } else { Some(v) }
                    }
                    None => None,
                };
                let reverse_primer = match reverse_primer {
                    Some(p) => {
                        let v = Self::normalize_iupac_text(&p)?;
                        if v.is_empty() { None } else { Some(v) }
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
                feature_geometry_mode,
                feature_boundary_mode,
                feature_strand_relation,
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
                    feature_geometry_mode,
                    feature_boundary_mode,
                    feature_strand_relation,
                    limit,
                    &mut result,
                )?;
            }
            Operation::GenerateCandidateSetBetweenAnchors {
                set_name,
                seq_id,
                anchor_a,
                anchor_b,
                length_bp,
                step_bp,
                limit,
            } => {
                self.op_generate_candidate_set_between_anchors(
                    set_name,
                    seq_id,
                    anchor_a,
                    anchor_b,
                    length_bp,
                    step_bp,
                    limit,
                    &mut result,
                )?;
            }
            Operation::DeleteCandidateSet { set_name } => {
                self.op_delete_candidate_set(set_name, &mut result)?;
            }
            Operation::UpsertGuideSet {
                guide_set_id,
                guides,
            } => {
                self.op_upsert_guide_set(guide_set_id, guides, &mut result)?;
            }
            Operation::DeleteGuideSet { guide_set_id } => {
                self.op_delete_guide_set(guide_set_id, &mut result)?;
            }
            Operation::FilterGuidesPractical {
                guide_set_id,
                config,
                output_guide_set_id,
            } => {
                self.op_filter_guides_practical(
                    guide_set_id,
                    config,
                    output_guide_set_id,
                    &mut result,
                )?;
            }
            Operation::GenerateGuideOligos {
                guide_set_id,
                template_id,
                apply_5prime_g_extension,
                output_oligo_set_id,
                passed_only,
            } => {
                self.op_generate_guide_oligos(
                    guide_set_id,
                    template_id,
                    apply_5prime_g_extension,
                    output_oligo_set_id,
                    passed_only,
                    &mut result,
                )?;
            }
            Operation::ExportGuideOligos {
                guide_set_id,
                oligo_set_id,
                format,
                path,
                plate_format,
            } => {
                self.op_export_guide_oligos(
                    guide_set_id,
                    oligo_set_id,
                    format,
                    path,
                    plate_format,
                    &mut result,
                )?;
            }
            Operation::ExportGuideProtocolText {
                guide_set_id,
                oligo_set_id,
                path,
                include_qc_checklist,
            } => {
                self.op_export_guide_protocol_text(
                    guide_set_id,
                    oligo_set_id,
                    path,
                    include_qc_checklist,
                    &mut result,
                )?;
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
                feature_geometry_mode,
                feature_boundary_mode,
                feature_strand_relation,
            } => {
                self.op_score_candidate_set_distance(
                    set_name,
                    metric,
                    feature_kinds,
                    feature_label_regex,
                    feature_geometry_mode,
                    feature_boundary_mode,
                    feature_strand_relation,
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
            Operation::ScoreCandidateSetWeightedObjective {
                set_name,
                metric,
                objectives,
                normalize_metrics,
            } => {
                self.op_score_candidate_set_weighted_objective(
                    set_name,
                    metric,
                    objectives,
                    normalize_metrics,
                    &mut result,
                )?;
            }
            Operation::TopKCandidateSet {
                input_set,
                output_set,
                metric,
                k,
                direction,
                tie_break,
            } => {
                self.op_top_k_candidate_set(
                    input_set,
                    output_set,
                    metric,
                    k,
                    direction,
                    tie_break,
                    &mut result,
                )?;
            }
            Operation::ParetoFrontierCandidateSet {
                input_set,
                output_set,
                objectives,
                max_candidates,
                tie_break,
            } => {
                self.op_pareto_frontier_candidate_set(
                    input_set,
                    output_set,
                    objectives,
                    max_candidates,
                    tie_break,
                    &mut result,
                )?;
            }
            Operation::UpsertWorkflowMacroTemplate {
                name,
                description,
                details_url,
                parameters,
                script,
            } => {
                self.op_upsert_workflow_macro_template(
                    name,
                    description,
                    details_url,
                    parameters,
                    script,
                    &mut result,
                )?;
            }
            Operation::DeleteWorkflowMacroTemplate { name } => {
                self.op_delete_workflow_macro_template(name, &mut result)?;
            }
            Operation::UpsertCandidateMacroTemplate {
                name,
                description,
                details_url,
                parameters,
                script,
            } => {
                self.op_upsert_candidate_macro_template(
                    name,
                    description,
                    details_url,
                    parameters,
                    script,
                    &mut result,
                )?;
            }
            Operation::DeleteCandidateMacroTemplate { name } => {
                self.op_delete_candidate_macro_template(name, &mut result)?;
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
            Operation::SetParameter { name, value } => match name.as_str() {
                "max_fragments_per_container" => {
                    let raw = value.as_u64().ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message:
                            "SetParameter max_fragments_per_container requires a positive integer"
                                .to_string(),
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
                "linear_external_feature_label_font_size"
                | "linear_feature_label_font_size"
                | "feature_label_font_size" => {
                    let raw = value.as_f64().ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("SetParameter {name} requires a number"),
                    })?;
                    if !raw.is_finite() {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message:
                                "linear_external_feature_label_font_size must be a finite number"
                                    .to_string(),
                        });
                    }
                    if !(8.0..=24.0).contains(&raw) {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: "linear_external_feature_label_font_size must be between 8.0 and 24.0".to_string(),
                        });
                    }
                    self.state.display.linear_external_feature_label_font_size = raw as f32;
                    result.messages.push(format!(
                        "Set parameter 'linear_external_feature_label_font_size' to {:.2}",
                        self.state.display.linear_external_feature_label_font_size
                    ));
                }
                "linear_external_feature_label_background_opacity"
                | "linear_feature_label_background_opacity"
                | "feature_label_background_opacity" => {
                    let raw = value.as_f64().ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("SetParameter {name} requires a number"),
                    })?;
                    if !raw.is_finite() {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: "linear_external_feature_label_background_opacity must be a finite number".to_string(),
                        });
                    }
                    if !(0.0..=1.0).contains(&raw) {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: "linear_external_feature_label_background_opacity must be between 0.0 and 1.0".to_string(),
                        });
                    }
                    self.state
                        .display
                        .linear_external_feature_label_background_opacity = raw as f32;
                    result.messages.push(format!(
                        "Set parameter 'linear_external_feature_label_background_opacity' to {:.3}",
                        self.state
                            .display
                            .linear_external_feature_label_background_opacity
                    ));
                }
                "regulatory_feature_max_view_span_bp"
                | "regulatory_max_view_span_bp"
                | "regulatory_max_span_bp" => {
                    let raw = value.as_u64().ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("SetParameter {name} requires a non-negative integer"),
                    })?;
                    self.state.display.regulatory_feature_max_view_span_bp = raw as usize;
                    result.messages.push(format!(
                        "Set parameter 'regulatory_feature_max_view_span_bp' to {}",
                        self.state.display.regulatory_feature_max_view_span_bp
                    ));
                }
                "linear_sequence_base_text_max_view_span_bp"
                | "linear_base_text_max_view_span_bp"
                | "sequence_base_text_max_view_span_bp" => {
                    let raw = value.as_u64().ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("SetParameter {name} requires a non-negative integer"),
                    })?;
                    self.state
                        .display
                        .linear_sequence_base_text_max_view_span_bp = raw as usize;
                    result.messages.push(format!(
                        "Set parameter 'linear_sequence_base_text_max_view_span_bp' to {}",
                        self.state
                            .display
                            .linear_sequence_base_text_max_view_span_bp
                    ));
                }
                "gc_content_bin_size_bp" | "gc_bin_size_bp" => {
                    let raw = value.as_u64().ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("SetParameter {name} requires a positive integer"),
                    })?;
                    if raw == 0 {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: "gc_content_bin_size_bp must be >= 1".to_string(),
                        });
                    }
                    self.state.display.gc_content_bin_size_bp = raw as usize;
                    result.messages.push(format!(
                        "Set parameter 'gc_content_bin_size_bp' to {}",
                        self.state.display.gc_content_bin_size_bp
                    ));
                }
                "linear_sequence_helical_max_view_span_bp" | "linear_helical_max_view_span_bp" => {
                    let raw = value.as_u64().ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("SetParameter {name} requires a non-negative integer"),
                    })?;
                    self.state.display.linear_sequence_helical_max_view_span_bp = raw as usize;
                    result.messages.push(format!(
                        "Set parameter 'linear_sequence_helical_max_view_span_bp' to {}",
                        self.state.display.linear_sequence_helical_max_view_span_bp
                    ));
                }
                "linear_sequence_condensed_max_view_span_bp"
                | "linear_condensed_max_view_span_bp" => {
                    let raw = value.as_u64().ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("SetParameter {name} requires a non-negative integer"),
                    })?;
                    self.state
                        .display
                        .linear_sequence_condensed_max_view_span_bp = raw as usize;
                    result.messages.push(format!(
                        "Set parameter 'linear_sequence_condensed_max_view_span_bp' to {}",
                        self.state
                            .display
                            .linear_sequence_condensed_max_view_span_bp
                    ));
                }
                "linear_sequence_letter_layout_mode" | "linear_helical_letter_layout_mode" => {
                    let raw = value.as_str().ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("SetParameter {name} requires a string value"),
                    })?;
                    let normalized = raw.trim().to_ascii_lowercase();
                    let layout_mode = match normalized.as_str() {
                        "continuous_helical" | "continuous-helical" | "continuous" => {
                            LinearSequenceLetterLayoutMode::ContinuousHelical
                        }
                        "condensed_10_row" | "condensed-10-row" | "condensed10row"
                        | "condensed" | "10_row" | "10-row" => {
                            LinearSequenceLetterLayoutMode::Condensed10Row
                        }
                        _ => {
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: format!(
                                    "Unsupported linear sequence letter layout mode '{raw}' (expected continuous_helical|condensed_10_row)"
                                ),
                            });
                        }
                    };
                    self.state.display.linear_sequence_letter_layout_mode = layout_mode;
                    result.messages.push(format!(
                        "Set parameter 'linear_sequence_letter_layout_mode' to {:?}",
                        self.state.display.linear_sequence_letter_layout_mode
                    ));
                }
                "linear_sequence_helical_phase_offset_bp" | "linear_helical_phase_offset_bp" => {
                    let raw = value.as_u64().ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("SetParameter {name} requires a non-negative integer"),
                    })?;
                    if raw > 9 {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message:
                                "linear_sequence_helical_phase_offset_bp must be between 0 and 9"
                                    .to_string(),
                        });
                    }
                    self.state.display.linear_sequence_helical_phase_offset_bp = raw as usize;
                    result.messages.push(format!(
                        "Set parameter 'linear_sequence_helical_phase_offset_bp' to {}",
                        self.state.display.linear_sequence_helical_phase_offset_bp
                    ));
                }
                "vcf_display_show_snp"
                | "vcf_display_show_ins"
                | "vcf_display_show_del"
                | "vcf_display_show_sv"
                | "vcf_display_show_other"
                | "vcf_display_pass_only"
                | "vcf_display_use_min_qual"
                | "vcf_display_use_max_qual"
                | "show_tfbs"
                | "tfbs_display_use_llr_bits"
                | "tfbs_display_use_llr_quantile"
                | "tfbs_display_use_true_log_odds_bits"
                | "tfbs_display_use_true_log_odds_quantile"
                | "auto_hide_sequence_panel_when_linear_bases_visible"
                | "linear_show_double_strand_bases"
                | "show_linear_double_strand_bases"
                | "linear_hide_backbone_when_sequence_bases_visible"
                | "hide_linear_backbone_when_bases_visible"
                | "linear_sequence_helical_letters_enabled"
                | "linear_reverse_strand_use_upside_down_letters"
                | "linear_reverse_strand_upside_down_letters" => {
                    let raw = value.as_bool().ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("SetParameter {} requires a boolean", name),
                    })?;
                    match name.as_str() {
                        "vcf_display_show_snp" => self.state.display.vcf_display_show_snp = raw,
                        "vcf_display_show_ins" => self.state.display.vcf_display_show_ins = raw,
                        "vcf_display_show_del" => self.state.display.vcf_display_show_del = raw,
                        "vcf_display_show_sv" => self.state.display.vcf_display_show_sv = raw,
                        "vcf_display_show_other" => self.state.display.vcf_display_show_other = raw,
                        "vcf_display_pass_only" => self.state.display.vcf_display_pass_only = raw,
                        "vcf_display_use_min_qual" => {
                            self.state.display.vcf_display_use_min_qual = raw
                        }
                        "vcf_display_use_max_qual" => {
                            self.state.display.vcf_display_use_max_qual = raw
                        }
                        "show_tfbs" => self.state.display.show_tfbs = raw,
                        "tfbs_display_use_llr_bits" => {
                            self.state.display.tfbs_display_use_llr_bits = raw
                        }
                        "tfbs_display_use_llr_quantile" => {
                            self.state.display.tfbs_display_use_llr_quantile = raw
                        }
                        "tfbs_display_use_true_log_odds_bits" => {
                            self.state.display.tfbs_display_use_true_log_odds_bits = raw
                        }
                        "tfbs_display_use_true_log_odds_quantile" => {
                            self.state.display.tfbs_display_use_true_log_odds_quantile = raw
                        }
                        "auto_hide_sequence_panel_when_linear_bases_visible" => {
                            self.state
                                .display
                                .auto_hide_sequence_panel_when_linear_bases_visible = raw
                        }
                        "linear_show_double_strand_bases" | "show_linear_double_strand_bases" => {
                            self.state.display.linear_show_double_strand_bases = raw
                        }
                        "linear_hide_backbone_when_sequence_bases_visible"
                        | "hide_linear_backbone_when_bases_visible" => {
                            self.state
                                .display
                                .linear_hide_backbone_when_sequence_bases_visible = raw
                        }
                        "linear_sequence_helical_letters_enabled" => {
                            self.state.display.linear_sequence_helical_letters_enabled = raw
                        }
                        "linear_reverse_strand_use_upside_down_letters"
                        | "linear_reverse_strand_upside_down_letters" => {
                            self.state
                                .display
                                .linear_reverse_strand_use_upside_down_letters = raw
                        }
                        _ => unreachable!(),
                    }
                    result
                        .messages
                        .push(format!("Set parameter '{}' to {}", name, raw));
                }
                "tfbs_display_min_llr_bits" | "tfbs_display_min_true_log_odds_bits" => {
                    let raw = value.as_f64().ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("SetParameter {} requires a number", name),
                    })?;
                    if !raw.is_finite() {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!("{} must be a finite number", name),
                        });
                    }
                    if name == "tfbs_display_min_llr_bits" {
                        self.state.display.tfbs_display_min_llr_bits = raw;
                    } else {
                        self.state.display.tfbs_display_min_true_log_odds_bits = raw;
                    }
                    result
                        .messages
                        .push(format!("Set parameter '{}' to {:.6}", name, raw));
                }
                "tfbs_display_min_llr_quantile" | "tfbs_display_min_true_log_odds_quantile" => {
                    let raw = value.as_f64().ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("SetParameter {} requires a number", name),
                    })?;
                    if !raw.is_finite() {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!("{} must be a finite number", name),
                        });
                    }
                    Self::validate_tf_thresholds(raw)?;
                    if name == "tfbs_display_min_llr_quantile" {
                        self.state.display.tfbs_display_min_llr_quantile = raw;
                    } else {
                        self.state.display.tfbs_display_min_true_log_odds_quantile = raw;
                    }
                    result
                        .messages
                        .push(format!("Set parameter '{}' to {:.6}", name, raw));
                }
                "vcf_display_min_qual" | "vcf_display_max_qual" => {
                    let raw = value.as_f64().ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("SetParameter {} requires a number", name),
                    })?;
                    if !raw.is_finite() {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!("{} must be a finite number", name),
                        });
                    }
                    if name == "vcf_display_min_qual" {
                        self.state.display.vcf_display_min_qual = raw;
                    } else {
                        self.state.display.vcf_display_max_qual = raw;
                    }
                    result
                        .messages
                        .push(format!("Set parameter '{}' to {:.6}", name, raw));
                }
                "vcf_display_required_info_keys"
                | "vcf_display_required_info_keys_csv"
                | "vcf_display_required_info" => {
                    let mut keys = if let Some(array) = value.as_array() {
                        let mut values = Vec::with_capacity(array.len());
                        for entry in array {
                            let Some(raw) = entry.as_str() else {
                                return Err(EngineError {
                                        code: ErrorCode::InvalidInput,
                                        message: "vcf_display_required_info_keys array entries must be strings".to_string(),
                                    });
                            };
                            values.push(raw.to_string());
                        }
                        values
                    } else if let Some(raw) = value.as_str() {
                        raw.split(',')
                            .map(str::trim)
                            .filter(|v| !v.is_empty())
                            .map(|v| v.to_string())
                            .collect::<Vec<_>>()
                    } else {
                        return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: "SetParameter vcf_display_required_info_keys requires a string (CSV) or string array".to_string(),
                            });
                    };
                    for key in &mut keys {
                        *key = key.trim().to_ascii_uppercase();
                    }
                    keys.retain(|key| !key.is_empty());
                    keys.sort();
                    keys.dedup();
                    self.state.display.vcf_display_required_info_keys = keys.clone();
                    result.messages.push(format!(
                        "Set parameter 'vcf_display_required_info_keys' to [{}]",
                        keys.join(",")
                    ));
                }
                _ => {
                    return Err(EngineError {
                        code: ErrorCode::Unsupported,
                        message: format!("Unknown parameter '{}'", name),
                    });
                }
            },
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
        let mut noop = |_p: OperationProgress| true;
        let checkpoint = self.maybe_capture_checkpoint(&op);
        let result = self.apply_internal(op.clone(), &run_id, &mut noop)?;
        self.journal.push(OperationRecord {
            run_id,
            op,
            result: result.clone(),
        });
        if let Some(checkpoint) = checkpoint {
            self.push_undo_checkpoint(checkpoint);
        }
        Ok(result)
    }

    fn apply_workflow(&mut self, wf: Workflow) -> Result<Vec<OpResult>, EngineError> {
        let mut results = Vec::new();
        for op in &wf.ops {
            let mut noop = |_p: OperationProgress| true;
            let checkpoint = self.maybe_capture_checkpoint(op);
            let result = self.apply_internal(op.clone(), &wf.run_id, &mut noop)?;
            self.journal.push(OperationRecord {
                run_id: wf.run_id.clone(),
                op: op.clone(),
                result: result.clone(),
            });
            if let Some(checkpoint) = checkpoint {
                self.push_undo_checkpoint(checkpoint);
            }
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
    use flate2::{Compression, write::GzEncoder};
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

    fn splicing_test_sequence() -> DNAsequence {
        let mut bases = vec![b'A'; 40];
        for (idx, base) in [
            (8usize, b'G'),
            (9, b'T'),
            (10, b'A'),
            (11, b'G'),
            (14, b'A'),
            (15, b'G'),
            (20, b'G'),
            (21, b'T'),
            (24, b'A'),
            (25, b'G'),
        ] {
            bases[idx] = base;
        }
        let sequence = String::from_utf8(bases).expect("valid ASCII DNA");
        let mut dna = DNAsequence::from_sequence(&sequence).expect("valid DNA sequence");
        dna.features_mut().push(gb_io::seq::Feature {
            kind: gb_io::seq::FeatureKind::from("mRNA"),
            location: gb_io::seq::Location::Join(vec![
                gb_io::seq::Location::simple_range(2, 8),
                gb_io::seq::Location::simple_range(12, 20),
                gb_io::seq::Location::simple_range(26, 34),
            ]),
            qualifiers: vec![
                ("gene".into(), Some("GENE1".to_string())),
                ("transcript_id".into(), Some("NM_TEST_1".to_string())),
                ("label".into(), Some("NM_TEST_1".to_string())),
            ],
        });
        dna.features_mut().push(gb_io::seq::Feature {
            kind: gb_io::seq::FeatureKind::from("mRNA"),
            location: gb_io::seq::Location::Join(vec![
                gb_io::seq::Location::simple_range(2, 8),
                gb_io::seq::Location::simple_range(16, 20),
                gb_io::seq::Location::simple_range(26, 34),
            ]),
            qualifiers: vec![
                ("gene".into(), Some("GENE1".to_string())),
                ("transcript_id".into(), Some("NM_TEST_2".to_string())),
                ("label".into(), Some("NM_TEST_2".to_string())),
            ],
        });
        dna
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

    fn candidate_store_env_lock() -> &'static std::sync::Mutex<()> {
        use std::sync::{Mutex, OnceLock};
        static LOCK: OnceLock<Mutex<()>> = OnceLock::new();
        LOCK.get_or_init(|| Mutex::new(()))
    }

    impl EnvVarGuard {
        fn set(key: &'static str, value: &str) -> Self {
            let previous = env::var(key).ok();
            unsafe {
                env::set_var(key, value);
            }
            Self { key, previous }
        }
    }

    impl Drop for EnvVarGuard {
        fn drop(&mut self) {
            match &self.previous {
                Some(value) => unsafe {
                    env::set_var(self.key, value);
                },
                None => unsafe {
                    env::remove_var(self.key);
                },
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
    fn test_extract_region_preserves_overlapping_features() {
        let mut state = ProjectState::default();
        let mut dna = seq("ATGCATGCATGC");
        dna.features_mut().push(gb_io::seq::Feature {
            kind: gb_io::seq::FeatureKind::from("gene"),
            location: gb_io::seq::Location::simple_range(2, 8),
            qualifiers: vec![("label".into(), Some("gene_a".to_string()))],
        });
        state.sequences.insert("x".to_string(), dna);
        let mut engine = GentleEngine::from_state(state);
        engine
            .apply(Operation::ExtractRegion {
                input: "x".to_string(),
                from: 2,
                to: 10,
                output_id: Some("part".to_string()),
            })
            .unwrap();
        let part = engine.state().sequences.get("part").expect("part sequence");
        let gene = part
            .features()
            .iter()
            .find(|feature| feature.kind.to_string().eq_ignore_ascii_case("gene"))
            .expect("gene should be preserved in extracted region");
        assert_eq!(gene.location.find_bounds().expect("gene bounds"), (0, 6));
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
        assert!(
            lineage
                .edges
                .iter()
                .any(|e| e.from_node_id == *x_node && e.to_node_id == *part_node)
        );
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
        assert!(
            lineage
                .edges
                .iter()
                .any(|e| e.from_node_id == *a_node && e.to_node_id == *ab_node)
        );
        assert!(
            lineage
                .edges
                .iter()
                .any(|e| e.from_node_id == *b_node && e.to_node_id == *ab_node)
        );
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
            assert!(
                lineage
                    .edges
                    .iter()
                    .any(|e| e.from_node_id == *s_node && e.to_node_id == *dnode)
            );
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
        assert!(
            res.messages
                .iter()
                .any(|m| m.contains("max_fragments_per_container"))
        );
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
        assert!(
            res.messages
                .iter()
                .any(|m| m.contains("feature_details_font_size"))
        );
        assert!((engine.state().display.feature_details_font_size - 9.5).abs() < f32::EPSILON);
    }

    #[test]
    fn test_set_parameter_linear_external_feature_label_style() {
        let mut engine = GentleEngine::new();
        engine
            .apply(Operation::SetParameter {
                name: "linear_external_feature_label_font_size".to_string(),
                value: serde_json::json!(12.5),
            })
            .unwrap();
        engine
            .apply(Operation::SetParameter {
                name: "linear_external_feature_label_background_opacity".to_string(),
                value: serde_json::json!(0.6),
            })
            .unwrap();
        assert!(
            (engine
                .state()
                .display
                .linear_external_feature_label_font_size
                - 12.5)
                .abs()
                < f32::EPSILON
        );
        assert!(
            (engine
                .state()
                .display
                .linear_external_feature_label_background_opacity
                - 0.6)
                .abs()
                < f32::EPSILON
        );
    }

    #[test]
    fn test_set_parameter_regulatory_feature_max_view_span_bp() {
        let mut engine = GentleEngine::new();
        let res = engine
            .apply(Operation::SetParameter {
                name: "regulatory_feature_max_view_span_bp".to_string(),
                value: serde_json::json!(50000),
            })
            .unwrap();
        assert!(
            res.messages
                .iter()
                .any(|m| m.contains("regulatory_feature_max_view_span_bp"))
        );
        assert_eq!(
            engine.state().display.regulatory_feature_max_view_span_bp,
            50_000
        );
    }

    #[test]
    fn test_set_parameter_gc_content_bin_size_bp() {
        let mut engine = GentleEngine::new();
        let res = engine
            .apply(Operation::SetParameter {
                name: "gc_content_bin_size_bp".to_string(),
                value: serde_json::json!(250),
            })
            .unwrap();
        assert!(
            res.messages
                .iter()
                .any(|m| m.contains("gc_content_bin_size_bp"))
        );
        assert_eq!(engine.state().display.gc_content_bin_size_bp, 250);
    }

    #[test]
    fn test_set_parameter_linear_sequence_base_text_max_view_span_bp() {
        let mut engine = GentleEngine::new();
        let res = engine
            .apply(Operation::SetParameter {
                name: "linear_sequence_base_text_max_view_span_bp".to_string(),
                value: serde_json::json!(1200),
            })
            .unwrap();
        assert!(
            res.messages
                .iter()
                .any(|m| m.contains("linear_sequence_base_text_max_view_span_bp"))
        );
        assert_eq!(
            engine
                .state()
                .display
                .linear_sequence_base_text_max_view_span_bp,
            1200
        );
    }

    #[test]
    fn test_set_parameter_linear_helical_display_controls() {
        let mut engine = GentleEngine::new();
        engine
            .apply(Operation::SetParameter {
                name: "linear_sequence_helical_letters_enabled".to_string(),
                value: serde_json::json!(true),
            })
            .unwrap();
        engine
            .apply(Operation::SetParameter {
                name: "linear_sequence_helical_max_view_span_bp".to_string(),
                value: serde_json::json!(2500),
            })
            .unwrap();
        engine
            .apply(Operation::SetParameter {
                name: "linear_sequence_condensed_max_view_span_bp".to_string(),
                value: serde_json::json!(1500),
            })
            .unwrap();
        engine
            .apply(Operation::SetParameter {
                name: "linear_sequence_helical_phase_offset_bp".to_string(),
                value: serde_json::json!(4),
            })
            .unwrap();
        engine
            .apply(Operation::SetParameter {
                name: "linear_sequence_letter_layout_mode".to_string(),
                value: serde_json::json!("condensed_10_row"),
            })
            .unwrap();
        engine
            .apply(Operation::SetParameter {
                name: "linear_hide_backbone_when_sequence_bases_visible".to_string(),
                value: serde_json::json!(true),
            })
            .unwrap();
        assert!(
            engine
                .state()
                .display
                .linear_sequence_helical_letters_enabled
        );
        assert_eq!(
            engine
                .state()
                .display
                .linear_sequence_helical_max_view_span_bp,
            2500
        );
        assert_eq!(
            engine
                .state()
                .display
                .linear_sequence_condensed_max_view_span_bp,
            1500
        );
        assert_eq!(
            engine
                .state()
                .display
                .linear_sequence_helical_phase_offset_bp,
            4
        );
        assert_eq!(
            engine.state().display.linear_sequence_letter_layout_mode,
            LinearSequenceLetterLayoutMode::Condensed10Row
        );
        assert!(
            engine
                .state()
                .display
                .linear_hide_backbone_when_sequence_bases_visible
        );
    }

    #[test]
    fn test_set_parameter_vcf_display_controls() {
        let mut engine = GentleEngine::new();
        engine
            .apply(Operation::SetParameter {
                name: "vcf_display_show_snp".to_string(),
                value: serde_json::json!(false),
            })
            .unwrap();
        engine
            .apply(Operation::SetParameter {
                name: "vcf_display_pass_only".to_string(),
                value: serde_json::json!(true),
            })
            .unwrap();
        engine
            .apply(Operation::SetParameter {
                name: "vcf_display_use_min_qual".to_string(),
                value: serde_json::json!(true),
            })
            .unwrap();
        engine
            .apply(Operation::SetParameter {
                name: "vcf_display_min_qual".to_string(),
                value: serde_json::json!(42.5),
            })
            .unwrap();
        engine
            .apply(Operation::SetParameter {
                name: "vcf_display_required_info_keys".to_string(),
                value: serde_json::json!("ac,ann"),
            })
            .unwrap();

        let display = &engine.state().display;
        assert!(!display.vcf_display_show_snp);
        assert!(display.vcf_display_pass_only);
        assert!(display.vcf_display_use_min_qual);
        assert!((display.vcf_display_min_qual - 42.5).abs() < f64::EPSILON);
        assert_eq!(
            display.vcf_display_required_info_keys,
            vec!["AC".to_string(), "ANN".to_string()]
        );
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
        assert!(
            err.message
                .contains("No amplicon introduced all requested mutations")
        );
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
        assert!(
            engine
                .list_sequences_with_genome_anchor()
                .iter()
                .any(|seq_id| seq_id == "tp73")
        );

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
        assert!(
            import_res
                .changed_seq_ids
                .iter()
                .any(|seq_id| seq_id == "tp73")
        );

        let tp73 = engine
            .state()
            .sequences
            .get("tp73")
            .expect("tp73 should exist");
        assert!(
            tp73.features()
                .iter()
                .any(GentleEngine::is_generated_genome_bed_feature)
        );
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
                    true
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
    fn test_inspect_tfbs_feature_expert_view() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("s".to_string(), seq("TTTACGTAAACGTGGG"));
        let mut engine = GentleEngine::from_state(state);
        engine
            .apply(Operation::AnnotateTfbs {
                seq_id: "s".to_string(),
                motifs: vec!["ACGT".to_string()],
                min_llr_bits: Some(0.0),
                min_llr_quantile: Some(0.0),
                per_tf_thresholds: vec![],
                clear_existing: Some(true),
                max_hits: Some(1),
            })
            .unwrap();
        let feature_id = engine
            .state()
            .sequences
            .get("s")
            .unwrap()
            .features()
            .iter()
            .position(|feature| {
                feature
                    .qualifier_values("gentle_generated".into())
                    .any(|v| v.eq_ignore_ascii_case("tfbs"))
            })
            .expect("a tfbs feature should exist");
        let view = engine
            .inspect_feature_expert("s", &FeatureExpertTarget::TfbsFeature { feature_id })
            .unwrap();
        match view {
            FeatureExpertView::Tfbs(tfbs) => {
                assert_eq!(tfbs.seq_id, "s");
                assert_eq!(tfbs.feature_id, feature_id);
                assert_eq!(tfbs.motif_length, 4);
                assert_eq!(tfbs.columns.len(), 4);
                assert_eq!(tfbs.matched_sequence.len(), 4);
                assert_eq!(tfbs.instruction, TFBS_EXPERT_INSTRUCTION);
                assert!(
                    tfbs.columns
                        .iter()
                        .all(|column| column.information_content_bits.is_finite())
                );
            }
            other => panic!("expected tfbs expert view, got {other:?}"),
        }
    }

    #[test]
    fn test_inspect_restriction_site_expert_view() {
        let mut dna = seq("AAGAATTCTT");
        *dna.restriction_enzymes_mut() = active_restriction_enzymes();
        dna.update_computed_features();
        let (key, names) = dna
            .restriction_enzyme_groups()
            .iter()
            .find(|(_, names)| names.iter().any(|name| name.eq_ignore_ascii_case("EcoRI")))
            .map(|(key, names)| (key.clone(), names.clone()))
            .expect("EcoRI site should exist in test sequence");
        let mut state = ProjectState::default();
        state.sequences.insert("s".to_string(), dna);
        let engine = GentleEngine::from_state(state);
        let view = engine
            .inspect_feature_expert(
                "s",
                &FeatureExpertTarget::RestrictionSite {
                    cut_pos_1based: key.pos() as usize + 1,
                    enzyme: Some("EcoRI".to_string()),
                    recognition_start_1based: Some(key.from() as usize + 1),
                    recognition_end_1based: Some(key.to() as usize),
                },
            )
            .unwrap();
        match view {
            FeatureExpertView::RestrictionSite(re) => {
                assert_eq!(re.seq_id, "s");
                assert_eq!(re.cut_pos_1based, key.pos() as usize + 1);
                assert_eq!(re.recognition_start_1based, key.from() as usize + 1);
                assert_eq!(re.recognition_end_1based, key.to() as usize);
                assert!(
                    re.enzyme_names
                        .iter()
                        .any(|name| name.eq_ignore_ascii_case("EcoRI"))
                );
                for name in names {
                    assert!(
                        re.enzyme_names
                            .iter()
                            .any(|entry| entry.eq_ignore_ascii_case(&name))
                    );
                }
                assert_eq!(re.instruction, RESTRICTION_EXPERT_INSTRUCTION);
            }
            other => panic!("expected restriction expert view, got {other:?}"),
        }
    }

    #[test]
    fn test_inspect_splicing_feature_expert_view() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("s".to_string(), splicing_test_sequence());
        let engine = GentleEngine::from_state(state);
        let view = engine
            .inspect_feature_expert("s", &FeatureExpertTarget::SplicingFeature { feature_id: 0 })
            .unwrap();
        match view {
            FeatureExpertView::Splicing(splicing) => {
                assert_eq!(splicing.seq_id, "s");
                assert_eq!(splicing.group_label, "GENE1");
                assert_eq!(splicing.transcript_count, 2);
                assert_eq!(splicing.transcripts.len(), 2);
                assert!(splicing.unique_exon_count >= 3);
                assert_eq!(splicing.matrix_rows.len(), 2);
                assert!(!splicing.unique_exons.is_empty());
                assert!(!splicing.boundaries.is_empty());
                assert!(!splicing.junctions.is_empty());
                assert!(
                    splicing
                        .boundaries
                        .iter()
                        .any(|m| m.side.eq_ignore_ascii_case("donor") && m.canonical)
                );
                assert!(
                    splicing
                        .boundaries
                        .iter()
                        .any(|m| m.side.eq_ignore_ascii_case("acceptor") && m.canonical)
                );
                assert_eq!(splicing.instruction, SPLICING_EXPERT_INSTRUCTION);
            }
            other => panic!("expected splicing expert view, got {other:?}"),
        }
    }

    #[test]
    fn test_render_feature_expert_svg_operation() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("s".to_string(), seq("TTTACGTAAACGTGGG"));
        let mut engine = GentleEngine::from_state(state);
        engine
            .apply(Operation::AnnotateTfbs {
                seq_id: "s".to_string(),
                motifs: vec!["ACGT".to_string()],
                min_llr_bits: Some(0.0),
                min_llr_quantile: Some(0.0),
                per_tf_thresholds: vec![],
                clear_existing: Some(true),
                max_hits: Some(1),
            })
            .unwrap();
        let feature_id = engine
            .state()
            .sequences
            .get("s")
            .unwrap()
            .features()
            .iter()
            .position(|feature| {
                feature
                    .qualifier_values("gentle_generated".into())
                    .any(|v| v.eq_ignore_ascii_case("tfbs"))
            })
            .expect("tfbs feature");
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = tmp.path().with_extension("feature.expert.svg");
        let path_text = path.display().to_string();
        let result = engine
            .apply(Operation::RenderFeatureExpertSvg {
                seq_id: "s".to_string(),
                target: FeatureExpertTarget::TfbsFeature { feature_id },
                path: path_text.clone(),
            })
            .unwrap();
        assert!(result.messages.iter().any(|m| m.contains("feature expert")));
        let text = std::fs::read_to_string(path_text).unwrap();
        assert!(text.contains("<svg"));
        assert!(text.contains("TFBS expert"));
    }

    #[test]
    fn test_render_splicing_feature_expert_svg_operation() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("s".to_string(), splicing_test_sequence());
        let mut engine = GentleEngine::from_state(state);
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = tmp.path().with_extension("splicing.feature.expert.svg");
        let path_text = path.display().to_string();
        let result = engine
            .apply(Operation::RenderFeatureExpertSvg {
                seq_id: "s".to_string(),
                target: FeatureExpertTarget::SplicingFeature { feature_id: 0 },
                path: path_text.clone(),
            })
            .unwrap();
        assert!(result.messages.iter().any(|m| m.contains("feature expert")));
        let text = std::fs::read_to_string(path_text).unwrap();
        assert!(text.contains("<svg"));
        assert!(text.contains("Splicing expert"));
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
    fn test_render_lineage_svg_includes_arrangement_nodes() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("a".to_string(), seq(&"ATGC".repeat(40)));
        state
            .sequences
            .insert("b".to_string(), seq(&"ATGC".repeat(55)));
        let mut engine = GentleEngine::from_state(state);

        let mut container_ids: Vec<String> = engine
            .state()
            .container_state
            .containers
            .keys()
            .cloned()
            .collect();
        container_ids.sort();
        engine
            .apply(Operation::CreateArrangementSerial {
                container_ids: container_ids.into_iter().take(2).collect(),
                arrangement_id: Some("arr-viz".to_string()),
                name: Some("Digest run".to_string()),
                ladders: Some(vec!["NEB 100bp DNA Ladder".to_string()]),
            })
            .unwrap();

        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = tmp.path().with_extension("lineage.arr.svg");
        let path_text = path.display().to_string();
        engine
            .apply(Operation::RenderLineageSvg {
                path: path_text.clone(),
            })
            .unwrap();
        let text = std::fs::read_to_string(path_text).unwrap();
        assert!(text.contains("arr-viz"));
        assert!(text.contains("lanes=2"));
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
                container_ids: None,
                arrangement_id: None,
            })
            .unwrap();
        assert!(res.messages.iter().any(|m| m.contains("serial gel SVG")));
        let text = std::fs::read_to_string(path_text).unwrap();
        assert!(text.contains("<svg"));
        assert!(text.contains("Serial Gel Preview"));
    }

    #[test]
    fn test_create_arrangement_serial_operation() {
        let mut state = ProjectState::default();
        state.sequences.insert("a".to_string(), seq("ATGCATGC"));
        state.sequences.insert("b".to_string(), seq("ATGCATGCATGC"));
        state.container_state.containers.insert(
            "container-1".to_string(),
            Container {
                container_id: "container-1".to_string(),
                kind: ContainerKind::Singleton,
                name: Some("Tube A".to_string()),
                members: vec!["a".to_string()],
                created_by_op: None,
                created_at_unix_ms: 0,
            },
        );
        state.container_state.containers.insert(
            "container-2".to_string(),
            Container {
                container_id: "container-2".to_string(),
                kind: ContainerKind::Singleton,
                name: Some("Tube B".to_string()),
                members: vec!["b".to_string()],
                created_by_op: None,
                created_at_unix_ms: 0,
            },
        );
        let mut engine = GentleEngine::from_state(state);
        let result = engine
            .apply(Operation::CreateArrangementSerial {
                container_ids: vec!["container-1".to_string(), "container-2".to_string()],
                arrangement_id: Some("arr-test".to_string()),
                name: Some("Digest run".to_string()),
                ladders: Some(vec!["NEB 1kb DNA Ladder".to_string()]),
            })
            .unwrap();
        assert!(
            result
                .messages
                .iter()
                .any(|m| m.contains("Created serial arrangement 'arr-test'"))
        );
        let arrangement = engine
            .state()
            .container_state
            .arrangements
            .get("arr-test")
            .expect("arrangement was created");
        assert_eq!(arrangement.mode, ArrangementMode::Serial);
        assert_eq!(
            arrangement.lane_container_ids,
            vec!["container-1".to_string(), "container-2".to_string()]
        );
        assert_eq!(arrangement.ladders, vec!["NEB 1kb DNA Ladder".to_string()]);
    }

    #[test]
    fn test_render_pool_gel_svg_operation_from_containers_and_arrangement() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("a".to_string(), seq(&"ATGC".repeat(70)));
        state
            .sequences
            .insert("b".to_string(), seq(&"ATGC".repeat(110)));
        state
            .sequences
            .insert("c".to_string(), seq(&"ATGC".repeat(150)));
        state.container_state.containers.insert(
            "container-1".to_string(),
            Container {
                container_id: "container-1".to_string(),
                kind: ContainerKind::Singleton,
                name: Some("Tube A".to_string()),
                members: vec!["a".to_string()],
                created_by_op: None,
                created_at_unix_ms: 0,
            },
        );
        state.container_state.containers.insert(
            "container-2".to_string(),
            Container {
                container_id: "container-2".to_string(),
                kind: ContainerKind::Pool,
                name: Some("Tube B".to_string()),
                members: vec!["b".to_string(), "c".to_string()],
                created_by_op: None,
                created_at_unix_ms: 0,
            },
        );
        state.container_state.arrangements.insert(
            "arr-1".to_string(),
            Arrangement {
                arrangement_id: "arr-1".to_string(),
                mode: ArrangementMode::Serial,
                name: Some("Run 1".to_string()),
                lane_container_ids: vec!["container-1".to_string(), "container-2".to_string()],
                ladders: vec!["NEB 100bp DNA Ladder".to_string()],
                created_by_op: None,
                created_at_unix_ms: 0,
            },
        );
        let mut engine = GentleEngine::from_state(state);

        let tmp_container = tempfile::NamedTempFile::new().unwrap();
        let path_container = tmp_container.path().with_extension("container.gel.svg");
        let path_container_text = path_container.display().to_string();
        let res_container = engine
            .apply(Operation::RenderPoolGelSvg {
                inputs: vec![],
                path: path_container_text.clone(),
                ladders: None,
                container_ids: Some(vec!["container-2".to_string()]),
                arrangement_id: None,
            })
            .unwrap();
        assert!(
            res_container
                .messages
                .iter()
                .any(|m| m.contains("serial gel SVG"))
        );
        let svg_container = std::fs::read_to_string(path_container_text).unwrap();
        assert!(svg_container.contains("Serial Gel Preview"));
        assert!(svg_container.contains("Tube B"));

        let tmp_arrangement = tempfile::NamedTempFile::new().unwrap();
        let path_arrangement = tmp_arrangement.path().with_extension("arrangement.gel.svg");
        let path_arrangement_text = path_arrangement.display().to_string();
        let res_arrangement = engine
            .apply(Operation::RenderPoolGelSvg {
                inputs: vec![],
                path: path_arrangement_text.clone(),
                ladders: None,
                container_ids: None,
                arrangement_id: Some("arr-1".to_string()),
            })
            .unwrap();
        assert!(
            res_arrangement
                .messages
                .iter()
                .any(|m| m.contains("2 sample lane(s)"))
        );
        let svg_arrangement = std::fs::read_to_string(path_arrangement_text).unwrap();
        assert!(svg_arrangement.contains("Serial Gel Preview"));
        assert!(svg_arrangement.contains("Tube A"));
        assert!(svg_arrangement.contains("Tube B"));
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
                container_ids: None,
                arrangement_id: None,
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
        assert!(
            catalog
                .ladders
                .iter()
                .any(|ladder| ladder.name == "NEB 100bp DNA Ladder")
        );
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
        assert!(
            res.messages
                .iter()
                .any(|m| m.contains("DNA ladders catalog"))
        );
        let text = std::fs::read_to_string(path_text).unwrap();
        let catalog: DnaLadderCatalog = serde_json::from_str(&text).unwrap();
        assert_eq!(catalog.schema, "gentle.dna_ladders.v1");
        assert!(catalog.ladder_count > 0);
        assert!(
            catalog
                .ladders
                .iter()
                .all(|ladder| ladder.name.to_ascii_lowercase().contains("neb"))
        );
    }

    #[test]
    fn test_inspect_rna_ladders() {
        let catalog = GentleEngine::inspect_rna_ladders(None);
        assert_eq!(catalog.schema, "gentle.rna_ladders.v1");
        assert!(catalog.ladder_count > 0);
        assert_eq!(catalog.ladder_count, catalog.ladders.len());
        assert!(
            catalog
                .ladders
                .iter()
                .any(|ladder| ladder.name.contains("RNA"))
        );
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
        assert!(
            res.messages
                .iter()
                .any(|m| m.contains("RNA ladders catalog"))
        );
        let text = std::fs::read_to_string(path_text).unwrap();
        let catalog: RnaLadderCatalog = serde_json::from_str(&text).unwrap();
        assert_eq!(catalog.schema, "gentle.rna_ladders.v1");
        assert!(catalog.ladder_count > 0);
        assert!(
            catalog
                .ladders
                .iter()
                .all(|ladder| ladder.name.to_ascii_lowercase().contains("neb"))
        );
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
    fn test_set_parameter_tfbs_display_filter_fields() {
        let mut engine = GentleEngine::new();
        engine
            .apply(Operation::SetParameter {
                name: "show_tfbs".to_string(),
                value: serde_json::json!(true),
            })
            .unwrap();
        engine
            .apply(Operation::SetParameter {
                name: "tfbs_display_use_llr_bits".to_string(),
                value: serde_json::json!(false),
            })
            .unwrap();
        engine
            .apply(Operation::SetParameter {
                name: "tfbs_display_min_llr_bits".to_string(),
                value: serde_json::json!(-2.5),
            })
            .unwrap();
        engine
            .apply(Operation::SetParameter {
                name: "tfbs_display_use_true_log_odds_quantile".to_string(),
                value: serde_json::json!(true),
            })
            .unwrap();
        engine
            .apply(Operation::SetParameter {
                name: "tfbs_display_min_true_log_odds_quantile".to_string(),
                value: serde_json::json!(0.8),
            })
            .unwrap();
        assert!(engine.state().display.show_tfbs);
        assert!(!engine.state().display.tfbs_display_use_llr_bits);
        assert!((engine.state().display.tfbs_display_min_llr_bits + 2.5).abs() < f64::EPSILON);
        assert!(
            engine
                .state()
                .display
                .tfbs_display_use_true_log_odds_quantile
        );
        assert!(
            (engine
                .state()
                .display
                .tfbs_display_min_true_log_odds_quantile
                - 0.8)
                .abs()
                < f64::EPSILON
        );
    }

    #[test]
    fn test_set_parameter_tfbs_display_quantile_out_of_range_fails() {
        let mut engine = GentleEngine::new();
        let err = engine
            .apply(Operation::SetParameter {
                name: "tfbs_display_min_llr_quantile".to_string(),
                value: serde_json::json!(1.1),
            })
            .unwrap_err();
        assert!(err.message.contains("between 0.0 and 1.0"));
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
        assert!(
            err.message
                .contains("feature_details_font_size requires a number")
        );
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
        assert!(
            err.message
                .contains("feature_details_font_size must be between 8.0 and 24.0")
        );
    }

    #[test]
    fn test_set_parameter_linear_external_feature_label_font_size_out_of_range_fails() {
        let mut engine = GentleEngine::new();
        let err = engine
            .apply(Operation::SetParameter {
                name: "linear_external_feature_label_font_size".to_string(),
                value: serde_json::json!(30.0),
            })
            .unwrap_err();
        assert!(
            err.message
                .contains("linear_external_feature_label_font_size must be between 8.0 and 24.0")
        );
    }

    #[test]
    fn test_set_parameter_linear_external_feature_label_background_opacity_out_of_range_fails() {
        let mut engine = GentleEngine::new();
        let err = engine
            .apply(Operation::SetParameter {
                name: "linear_external_feature_label_background_opacity".to_string(),
                value: serde_json::json!(1.5),
            })
            .unwrap_err();
        assert!(err.message.contains(
            "linear_external_feature_label_background_opacity must be between 0.0 and 1.0"
        ));
    }

    #[test]
    fn test_set_parameter_linear_sequence_helical_phase_offset_out_of_range_fails() {
        let mut engine = GentleEngine::new();
        let err = engine
            .apply(Operation::SetParameter {
                name: "linear_sequence_helical_phase_offset_bp".to_string(),
                value: serde_json::json!(10),
            })
            .unwrap_err();
        assert!(
            err.message
                .contains("linear_sequence_helical_phase_offset_bp must be between 0 and 9")
        );
    }

    #[test]
    fn test_set_parameter_regulatory_feature_max_view_span_bp_invalid_type_fails() {
        let mut engine = GentleEngine::new();
        let err = engine
            .apply(Operation::SetParameter {
                name: "regulatory_feature_max_view_span_bp".to_string(),
                value: serde_json::json!("wide"),
            })
            .unwrap_err();
        assert!(
            err.message
                .contains("regulatory_feature_max_view_span_bp requires a non-negative integer")
        );
    }

    #[test]
    fn test_set_parameter_gc_content_bin_size_bp_invalid_type_fails() {
        let mut engine = GentleEngine::new();
        let err = engine
            .apply(Operation::SetParameter {
                name: "gc_content_bin_size_bp".to_string(),
                value: serde_json::json!("fine"),
            })
            .unwrap_err();
        assert!(
            err.message
                .contains("gc_content_bin_size_bp requires a positive integer")
        );
    }

    #[test]
    fn test_set_parameter_gc_content_bin_size_bp_zero_fails() {
        let mut engine = GentleEngine::new();
        let err = engine
            .apply(Operation::SetParameter {
                name: "gc_content_bin_size_bp".to_string(),
                value: serde_json::json!(0),
            })
            .unwrap_err();
        assert!(err.message.contains("gc_content_bin_size_bp must be >= 1"));
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
                timeout_seconds: None,
            })
            .unwrap();
        assert!(
            prep.messages
                .iter()
                .any(|m| m.contains("Prepared genome 'ToyGenome'"))
        );
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
    fn test_prepare_genome_operation_supports_timeout_seconds() {
        let td = tempdir().unwrap();
        let root = td.path();
        let fasta = root.join("toy.fa");
        let ann = root.join("toy.gtf");
        fs::write(&fasta, ">chr1\nACGT\nACGT\n").unwrap();
        fs::write(
            &ann,
            "chr1\tsrc\tgene\t1\t8\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\";\n",
        )
        .unwrap();
        let cache_dir = root.join("cache");
        let catalog_path = root.join("catalog.json");
        let catalog_json = format!(
            r#"{{
  "ToyGenome": {{
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
            fasta.display(),
            ann.display(),
            cache_dir.display()
        );
        fs::write(&catalog_path, catalog_json).unwrap();

        let mut engine = GentleEngine::new();
        let err = engine
            .apply(Operation::PrepareGenome {
                genome_id: "ToyGenome".to_string(),
                catalog_path: Some(catalog_path.to_string_lossy().to_string()),
                cache_dir: None,
                timeout_seconds: Some(0),
            })
            .unwrap_err();
        assert!(err.message.to_ascii_lowercase().contains("timed out"));
    }

    #[test]
    fn test_extend_genome_anchor_plus_strand_adds_lineage_and_provenance() {
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
                timeout_seconds: None,
            })
            .unwrap();
        engine
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

        let extended = engine
            .apply(Operation::ExtendGenomeAnchor {
                seq_id: "toy_slice".to_string(),
                side: GenomeAnchorSide::FivePrime,
                length_bp: 2,
                output_id: Some("toy_slice_ext5".to_string()),
                catalog_path: Some(catalog_path_str),
                cache_dir: None,
            })
            .unwrap();
        assert_eq!(extended.created_seq_ids, vec!["toy_slice_ext5".to_string()]);
        let extended_seq = engine
            .state()
            .sequences
            .get("toy_slice_ext5")
            .expect("extended sequence should exist");
        assert_eq!(extended_seq.get_forward_string(), "ACGTACGTAC");

        let lineage = &engine.state().lineage;
        let parent = lineage
            .seq_to_node
            .get("toy_slice")
            .expect("parent lineage node should exist");
        let child = lineage
            .seq_to_node
            .get("toy_slice_ext5")
            .expect("child lineage node should exist");
        assert!(
            lineage
                .edges
                .iter()
                .any(|edge| edge.from_node_id == *parent && edge.to_node_id == *child)
        );

        let provenance = engine
            .state()
            .metadata
            .get(PROVENANCE_METADATA_KEY)
            .and_then(|v| v.as_object())
            .expect("provenance metadata object");
        let extractions = provenance
            .get(GENOME_EXTRACTIONS_METADATA_KEY)
            .and_then(|v| v.as_array())
            .expect("genome_extractions array");
        let entry = extractions
            .iter()
            .find(|entry| {
                entry
                    .get("seq_id")
                    .and_then(|v| v.as_str())
                    .map(|v| v == "toy_slice_ext5")
                    .unwrap_or(false)
            })
            .expect("extended provenance entry");
        assert_eq!(
            entry.get("operation").and_then(|v| v.as_str()),
            Some("ExtendGenomeAnchor")
        );
        assert_eq!(entry.get("start_1based").and_then(|v| v.as_u64()), Some(1));
        assert_eq!(entry.get("end_1based").and_then(|v| v.as_u64()), Some(10));
        assert_eq!(
            entry.get("anchor_strand").and_then(|v| v.as_str()),
            Some("+")
        );
    }

    #[test]
    fn test_extend_genome_anchor_reverse_strand_respects_5prime_and_3prime_physical_direction() {
        let td = tempdir().unwrap();
        let root = td.path();
        let fasta_gz = root.join("toy.fa.gz");
        let ann_gz = root.join("toy.gtf.gz");
        write_gzip(&fasta_gz, ">chr1\nACGTTGCAATGCCGTA\n");
        write_gzip(
            &ann_gz,
            "chr1\tsrc\tgene\t1\t16\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\";\n",
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

        let mut state = ProjectState::default();
        state
            .sequences
            .insert("rev_anchor".to_string(), seq("ATTGCA"));
        state.metadata.insert(
            PROVENANCE_METADATA_KEY.to_string(),
            serde_json::json!({
                GENOME_EXTRACTIONS_METADATA_KEY: [
                    {
                        "seq_id": "rev_anchor",
                        "recorded_at_unix_ms": 1,
                        "operation": "LoadFileGenBankRegion",
                        "genome_id": "ToyGenome",
                        "catalog_path": catalog_path_str,
                        "cache_dir": null,
                        "chromosome": "chr1",
                        "start_1based": 5,
                        "end_1based": 10,
                        "gene_query": null,
                        "occurrence": null,
                        "gene_id": null,
                        "gene_name": null,
                        "strand": null,
                        "anchor_strand": "-",
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
        engine
            .apply(Operation::PrepareGenome {
                genome_id: "ToyGenome".to_string(),
                catalog_path: Some(catalog_path.to_string_lossy().to_string()),
                cache_dir: None,
                timeout_seconds: None,
            })
            .unwrap();

        let ext5 = engine
            .apply(Operation::ExtendGenomeAnchor {
                seq_id: "rev_anchor".to_string(),
                side: GenomeAnchorSide::FivePrime,
                length_bp: 3,
                output_id: Some("rev_ext5".to_string()),
                catalog_path: None,
                cache_dir: None,
            })
            .unwrap();
        assert_eq!(ext5.created_seq_ids, vec!["rev_ext5".to_string()]);
        assert_eq!(
            engine
                .state()
                .sequences
                .get("rev_ext5")
                .expect("rev_ext5 sequence")
                .get_forward_string(),
            "GGCATTGCA"
        );

        let ext3 = engine
            .apply(Operation::ExtendGenomeAnchor {
                seq_id: "rev_anchor".to_string(),
                side: GenomeAnchorSide::ThreePrime,
                length_bp: 2,
                output_id: Some("rev_ext3".to_string()),
                catalog_path: None,
                cache_dir: None,
            })
            .unwrap();
        assert_eq!(ext3.created_seq_ids, vec!["rev_ext3".to_string()]);
        assert_eq!(
            engine
                .state()
                .sequences
                .get("rev_ext3")
                .expect("rev_ext3 sequence")
                .get_forward_string(),
            "ATTGCAAC"
        );

        let provenance = engine
            .state()
            .metadata
            .get(PROVENANCE_METADATA_KEY)
            .and_then(|v| v.as_object())
            .expect("provenance metadata object");
        let extractions = provenance
            .get(GENOME_EXTRACTIONS_METADATA_KEY)
            .and_then(|v| v.as_array())
            .expect("genome_extractions array");

        let ext5_entry = extractions
            .iter()
            .find(|entry| {
                entry
                    .get("seq_id")
                    .and_then(|v| v.as_str())
                    .map(|v| v == "rev_ext5")
                    .unwrap_or(false)
            })
            .expect("rev_ext5 provenance");
        assert_eq!(
            ext5_entry.get("start_1based").and_then(|v| v.as_u64()),
            Some(5)
        );
        assert_eq!(
            ext5_entry.get("end_1based").and_then(|v| v.as_u64()),
            Some(13)
        );
        assert_eq!(
            ext5_entry.get("anchor_strand").and_then(|v| v.as_str()),
            Some("-")
        );

        let ext3_entry = extractions
            .iter()
            .find(|entry| {
                entry
                    .get("seq_id")
                    .and_then(|v| v.as_str())
                    .map(|v| v == "rev_ext3")
                    .unwrap_or(false)
            })
            .expect("rev_ext3 provenance");
        assert_eq!(
            ext3_entry.get("start_1based").and_then(|v| v.as_u64()),
            Some(3)
        );
        assert_eq!(
            ext3_entry.get("end_1based").and_then(|v| v.as_u64()),
            Some(10)
        );
        assert_eq!(
            ext3_entry.get("anchor_strand").and_then(|v| v.as_str()),
            Some("-")
        );
    }

    #[test]
    fn test_extend_genome_anchor_rejects_zero_length() {
        let mut engine = GentleEngine::new();
        let err = engine
            .apply(Operation::ExtendGenomeAnchor {
                seq_id: "missing".to_string(),
                side: GenomeAnchorSide::FivePrime,
                length_bp: 0,
                output_id: None,
                catalog_path: None,
                cache_dir: None,
            })
            .unwrap_err();
        assert!(matches!(err.code, ErrorCode::InvalidInput));
        assert!(err.message.contains("length_bp >= 1"));
    }

    #[test]
    fn test_extend_genome_anchor_requires_genome_anchor_provenance() {
        let mut state = ProjectState::default();
        state.sequences.insert("plain".to_string(), seq("ACGTACGT"));
        let mut engine = GentleEngine::from_state(state);
        let err = engine
            .apply(Operation::ExtendGenomeAnchor {
                seq_id: "plain".to_string(),
                side: GenomeAnchorSide::ThreePrime,
                length_bp: 10,
                output_id: None,
                catalog_path: None,
                cache_dir: None,
            })
            .unwrap_err();
        assert!(matches!(err.code, ErrorCode::NotFound));
        assert!(err.message.contains("no genome anchor provenance"));
    }

    #[test]
    fn test_extend_genome_anchor_warns_when_clipped_at_chromosome_start() {
        let td = tempdir().unwrap();
        let root = td.path();
        let fasta = root.join("toy.fa");
        let gtf = root.join("toy.gtf");
        fs::write(&fasta, ">chr1\nACGTACGTACGT\n").unwrap();
        fs::write(
            &gtf,
            "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n",
        )
        .unwrap();
        let catalog_path = root.join("catalog.json");
        let cache_dir = root.join("cache");
        let catalog_json = format!(
            r#"{{
  "ToyGenome": {{
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
            fasta.display(),
            gtf.display(),
            cache_dir.display()
        );
        fs::write(&catalog_path, catalog_json).unwrap();
        let catalog_path_str = catalog_path.to_string_lossy().to_string();

        let mut state = ProjectState::default();
        state.sequences.insert("anch".to_string(), seq("CGTAC"));
        state.metadata.insert(
            PROVENANCE_METADATA_KEY.to_string(),
            serde_json::json!({
                GENOME_EXTRACTIONS_METADATA_KEY: [
                    {
                        "seq_id": "anch",
                        "recorded_at_unix_ms": 1,
                        "operation": "ExtractGenomeRegion",
                        "genome_id": "ToyGenome",
                        "catalog_path": catalog_path_str,
                        "cache_dir": null,
                        "chromosome": "chr1",
                        "start_1based": 2,
                        "end_1based": 6,
                        "gene_query": null,
                        "occurrence": null,
                        "gene_id": null,
                        "gene_name": null,
                        "strand": null,
                        "anchor_strand": "+",
                        "sequence_source_type": "local",
                        "annotation_source_type": "local",
                        "sequence_source": "local",
                        "annotation_source": "local",
                        "sequence_sha1": null,
                        "annotation_sha1": null
                    }
                ]
            }),
        );
        let mut engine = GentleEngine::from_state(state);
        engine
            .apply(Operation::PrepareGenome {
                genome_id: "ToyGenome".to_string(),
                catalog_path: Some(catalog_path.to_string_lossy().to_string()),
                cache_dir: None,
                timeout_seconds: None,
            })
            .unwrap();

        let result = engine
            .apply(Operation::ExtendGenomeAnchor {
                seq_id: "anch".to_string(),
                side: GenomeAnchorSide::FivePrime,
                length_bp: 10,
                output_id: Some("anch_ext".to_string()),
                catalog_path: None,
                cache_dir: None,
            })
            .unwrap();
        assert_eq!(result.created_seq_ids, vec!["anch_ext".to_string()]);
        assert!(
            result
                .warnings
                .iter()
                .any(|w| w.contains("clipped at chromosome start position 1"))
        );
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
                timeout_seconds: None,
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
        assert!(
            generated_gz[0]
                .qualifier_values("label".into())
                .next()
                .map(|v| v.contains("peak_a"))
                .unwrap_or(false)
        );
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
                timeout_seconds: None,
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
        assert!(
            plain
                .warnings
                .iter()
                .any(|w| w.contains("did not match anchor chromosome"))
        );
        assert!(plain.warnings.iter().any(|w| w.contains("chr2")));

        let dna_plain = engine.state().sequences.get("toy_slice").unwrap();
        let plain_features: Vec<_> = dna_plain
            .features()
            .iter()
            .filter(|f| GentleEngine::is_generated_genome_vcf_feature(f))
            .collect();
        assert_eq!(plain_features.len(), 3);
        assert!(
            plain_features
                .iter()
                .any(|f| { f.qualifier_values("vcf_alt".into()).any(|v| v == "<DEL>") })
        );

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
        assert!(
            filtered
                .warnings
                .iter()
                .any(|w| w.contains("QUAL-based score filters"))
        );

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
    fn test_import_genome_vcf_track_chrom_alias_and_genotype_summary() {
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
                        "chromosome": "1",
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
        let vcf_path = td.path().join("genotype_variants.vcf");
        fs::write(
            &vcf_path,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_a\tsample_b\nchr01\t101\trsGT\tA\tG,TT\t51\tPASS\tAC=2;AN=4\tGT:DP\t0/1:12\t1|1:20\n",
        )
        .unwrap();
        let result = engine
            .apply(Operation::ImportGenomeVcfTrack {
                seq_id: "toy_slice".to_string(),
                path: vcf_path.to_string_lossy().to_string(),
                track_name: Some("gt_track".to_string()),
                min_score: None,
                max_score: None,
                clear_existing: Some(true),
            })
            .unwrap();
        assert!(result.changed_seq_ids.contains(&"toy_slice".to_string()));
        let dna = engine.state().sequences.get("toy_slice").unwrap();
        let features: Vec<_> = dna
            .features()
            .iter()
            .filter(|f| GentleEngine::is_generated_genome_vcf_feature(f))
            .collect();
        assert_eq!(features.len(), 2);
        let alt1 = features
            .iter()
            .find(|f| f.qualifier_values("vcf_alt".into()).any(|v| v == "G"))
            .expect("ALT=G feature");
        assert_eq!(
            alt1.qualifier_values("vcf_variant_class".into())
                .next()
                .unwrap_or_default(),
            "SNP"
        );
        assert_eq!(
            alt1.qualifier_values("vcf_alt_allele_index".into())
                .next()
                .unwrap_or_default(),
            "1"
        );
        assert_eq!(
            alt1.qualifier_values("vcf_alt_carriers".into())
                .next()
                .unwrap_or_default(),
            "2"
        );
        assert_eq!(
            alt1.qualifier_values("vcf_gt_het".into())
                .next()
                .unwrap_or_default(),
            "1"
        );
        assert_eq!(
            alt1.qualifier_values("vcf_gt_hom_alt".into())
                .next()
                .unwrap_or_default(),
            "1"
        );
        assert_eq!(
            alt1.qualifier_values("vcf_phase".into())
                .next()
                .unwrap_or_default(),
            "mixed"
        );
        assert_eq!(
            alt1.qualifier_values("vcf_alt_carrier_samples".into())
                .next()
                .unwrap_or_default(),
            "sample_a,sample_b"
        );
    }

    #[test]
    fn test_import_genome_vcf_track_can_cancel_via_progress_callback() {
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
        let vcf_path = td.path().join("many_variants.vcf");
        let mut payload =
            String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        for i in 0..600usize {
            payload.push_str(&format!(
                "chr1\t{}\trs{}\tA\tG\t50\tPASS\tAC=1\n",
                100 + (i % 20),
                i + 1
            ));
        }
        fs::write(&vcf_path, payload).unwrap();

        let mut cancelled = false;
        let result = engine
            .apply_with_progress(
                Operation::ImportGenomeVcfTrack {
                    seq_id: "toy_slice".to_string(),
                    path: vcf_path.to_string_lossy().to_string(),
                    track_name: Some("many".to_string()),
                    min_score: None,
                    max_score: None,
                    clear_existing: Some(true),
                },
                |progress| match progress {
                    OperationProgress::GenomeTrackImport(p) => {
                        if !p.done && p.parsed_records >= 250 {
                            cancelled = true;
                            false
                        } else {
                            true
                        }
                    }
                    _ => true,
                },
            )
            .unwrap();
        assert!(cancelled);
        assert!(
            result
                .warnings
                .iter()
                .any(|w| w.contains("import cancelled"))
        );
        let dna = engine.state().sequences.get("toy_slice").unwrap();
        let features: Vec<_> = dna
            .features()
            .iter()
            .filter(|f| GentleEngine::is_generated_genome_vcf_feature(f))
            .collect();
        assert!(features.len() < 600);
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
                timeout_seconds: None,
            })
            .unwrap();
        assert!(
            prep.messages
                .iter()
                .any(|m| m.contains("[genbank_accession]"))
        );

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
        assert!(
            !engine
                .state()
                .metadata
                .contains_key(GENOME_TRACK_KNOWN_ANCHORS_METADATA_KEY)
        );
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

    #[test]
    fn test_candidate_store_save_externalizes_and_load_hydrates() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("seqA".to_string(), seq("ACGTACGTACGT"));
        let mut engine = GentleEngine::from_state(state);
        engine
            .apply(Operation::GenerateCandidateSet {
                set_name: "windows".to_string(),
                seq_id: "seqA".to_string(),
                length_bp: 4,
                step_bp: 4,
                feature_kinds: vec![],
                feature_label_regex: None,
                max_distance_bp: None,
                feature_geometry_mode: None,
                feature_boundary_mode: None,
                feature_strand_relation: None,
                limit: Some(32),
            })
            .expect("generate candidates");

        let td = tempdir().expect("tempdir");
        let project_path = td.path().join("demo.gentle.json");
        engine
            .state()
            .save_to_path(project_path.to_string_lossy().as_ref())
            .expect("save project");

        let project_text = std::fs::read_to_string(&project_path).expect("read project");
        let project_json: serde_json::Value =
            serde_json::from_str(&project_text).expect("parse project");
        let candidate_meta = project_json
            .get("metadata")
            .and_then(|m| m.get(CANDIDATE_SETS_METADATA_KEY))
            .expect("candidate metadata present");
        assert_eq!(
            candidate_meta
                .get("schema")
                .and_then(|v| v.as_str())
                .unwrap_or_default(),
            CANDIDATE_SETS_REF_SCHEMA
        );

        let index_rel = candidate_meta
            .get("index_path")
            .and_then(|v| v.as_str())
            .expect("index path present");
        let index_abs = project_path
            .parent()
            .expect("project parent")
            .join(index_rel);
        assert!(index_abs.exists(), "candidate index file should exist");

        let loaded =
            ProjectState::load_from_path(project_path.to_string_lossy().as_ref()).expect("load");
        let loaded_meta = loaded
            .metadata
            .get(CANDIDATE_SETS_METADATA_KEY)
            .expect("hydrated candidate metadata");
        assert_eq!(
            loaded_meta
                .get("schema")
                .and_then(|v| v.as_str())
                .unwrap_or_default(),
            CANDIDATE_SETS_SCHEMA
        );
        let loaded_engine = GentleEngine::from_state(loaded);
        let summaries = loaded_engine.list_candidate_sets();
        assert_eq!(summaries.len(), 1);
        assert_eq!(summaries[0].name, "windows");
        assert_eq!(summaries[0].candidate_count, 3);
    }

    #[test]
    fn test_project_state_load_degrades_when_candidate_sidecar_is_missing() {
        let _lock = candidate_store_env_lock().lock().unwrap();
        let _guard = EnvVarGuard::set(CANDIDATE_STORE_STRICT_LOAD_ENV, "0");
        let td = tempdir().expect("tempdir");
        let project_path = td.path().join("broken.gentle.json");
        let project_json = serde_json::json!({
            "sequences": {},
            "metadata": {
                CANDIDATE_SETS_METADATA_KEY: {
                    "schema": CANDIDATE_SETS_REF_SCHEMA,
                    "storage": "jsonl_indexed",
                    "index_path": "missing_sidecar/index.json",
                    "set_count": 1,
                    "updated_at_unix_ms": 0
                }
            }
        });
        std::fs::write(
            &project_path,
            serde_json::to_string_pretty(&project_json).expect("serialize project"),
        )
        .expect("write project");
        let loaded = ProjectState::load_from_path(project_path.to_string_lossy().as_ref())
            .expect("load in degraded mode");
        assert!(!loaded.metadata.contains_key(CANDIDATE_SETS_METADATA_KEY));
        let warning = loaded
            .metadata
            .get(CANDIDATE_SETS_LOAD_WARNING_METADATA_KEY)
            .expect("degraded-load warning metadata");
        let warning_message = warning
            .get("message")
            .and_then(|v| v.as_str())
            .unwrap_or_default();
        assert!(warning_message.contains("candidate-store index"));
    }

    #[test]
    fn test_project_state_load_strict_mode_errors_when_candidate_sidecar_is_missing() {
        let _lock = candidate_store_env_lock().lock().unwrap();
        let _guard = EnvVarGuard::set(CANDIDATE_STORE_STRICT_LOAD_ENV, "1");
        let td = tempdir().expect("tempdir");
        let project_path = td.path().join("broken_strict.gentle.json");
        let project_json = serde_json::json!({
            "sequences": {},
            "metadata": {
                CANDIDATE_SETS_METADATA_KEY: {
                    "schema": CANDIDATE_SETS_REF_SCHEMA,
                    "storage": "jsonl_indexed",
                    "index_path": "missing_sidecar/index.json",
                    "set_count": 1,
                    "updated_at_unix_ms": 0
                }
            }
        });
        std::fs::write(
            &project_path,
            serde_json::to_string_pretty(&project_json).expect("serialize project"),
        )
        .expect("write project");
        let err =
            ProjectState::load_from_path(project_path.to_string_lossy().as_ref()).unwrap_err();
        assert!(err.message.contains("candidate-store index"));
    }

    #[test]
    fn test_candidate_store_save_replaces_sidecar_and_removes_stale_files() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("seqA".to_string(), seq("ACGTACGTACGT"));
        let mut engine = GentleEngine::from_state(state);
        engine
            .apply(Operation::GenerateCandidateSet {
                set_name: "windows".to_string(),
                seq_id: "seqA".to_string(),
                length_bp: 4,
                step_bp: 2,
                feature_kinds: vec![],
                feature_label_regex: None,
                max_distance_bp: None,
                feature_geometry_mode: None,
                feature_boundary_mode: None,
                feature_strand_relation: None,
                limit: Some(16),
            })
            .expect("generate candidates");

        let td = tempdir().expect("tempdir");
        let project_path = td.path().join("replace_sidecar.gentle.json");
        engine
            .state()
            .save_to_path(project_path.to_string_lossy().as_ref())
            .expect("initial save");

        let sidecar_dir = ProjectState::candidate_store_sidecar_dir(&project_path);
        let stale_path = sidecar_dir.join("stale.jsonl");
        std::fs::write(&stale_path, "{\"stale\":true}\n").expect("write stale sidecar file");
        assert!(
            stale_path.exists(),
            "stale file should exist before re-save"
        );

        engine
            .state()
            .save_to_path(project_path.to_string_lossy().as_ref())
            .expect("second save");

        assert!(
            !stale_path.exists(),
            "stale sidecar file should be removed by directory replacement"
        );
    }

    #[test]
    fn test_candidate_generation_regex_anchor_and_filter_quantile_edges() {
        let mut state = ProjectState::default();
        let mut dna = DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence");
        dna.features_mut().push(gb_io::seq::Feature {
            kind: gb_io::seq::FeatureKind::from("gene"),
            location: gb_io::seq::Location::simple_range(0, 0),
            qualifiers: vec![("label".into(), Some("TP53".to_string()))],
        });
        dna.features_mut().push(gb_io::seq::Feature {
            kind: gb_io::seq::FeatureKind::from("gene"),
            location: gb_io::seq::Location::simple_range(11, 11),
            qualifiers: vec![("label".into(), Some("TP53-AS1".to_string()))],
        });
        state.sequences.insert("seqA".to_string(), dna);
        let mut engine = GentleEngine::from_state(state);

        engine
            .apply(Operation::GenerateCandidateSet {
                set_name: "gene_all".to_string(),
                seq_id: "seqA".to_string(),
                length_bp: 1,
                step_bp: 1,
                feature_kinds: vec!["gene".to_string()],
                feature_label_regex: None,
                max_distance_bp: Some(0),
                feature_geometry_mode: None,
                feature_boundary_mode: None,
                feature_strand_relation: None,
                limit: Some(64),
            })
            .expect("generate all gene-anchored candidates");
        let gene_all_count = engine
            .list_candidate_sets()
            .into_iter()
            .find(|s| s.name == "gene_all")
            .map(|s| s.candidate_count)
            .expect("gene_all set exists");

        engine
            .apply(Operation::GenerateCandidateSet {
                set_name: "tp53_only".to_string(),
                seq_id: "seqA".to_string(),
                length_bp: 1,
                step_bp: 1,
                feature_kinds: vec!["gene".to_string()],
                feature_label_regex: Some("^TP53$".to_string()),
                max_distance_bp: Some(0),
                feature_geometry_mode: None,
                feature_boundary_mode: None,
                feature_strand_relation: None,
                limit: Some(64),
            })
            .expect("generate regex-anchored candidates");
        let regex_count = engine
            .list_candidate_sets()
            .into_iter()
            .find(|s| s.name == "tp53_only")
            .map(|s| s.candidate_count)
            .expect("regex set exists");
        assert!(regex_count >= 1);
        assert!(regex_count < gene_all_count);

        engine
            .apply(Operation::GenerateCandidateSet {
                set_name: "windows".to_string(),
                seq_id: "seqA".to_string(),
                length_bp: 4,
                step_bp: 2,
                feature_kinds: vec![],
                feature_label_regex: None,
                max_distance_bp: None,
                feature_geometry_mode: None,
                feature_boundary_mode: None,
                feature_strand_relation: None,
                limit: Some(64),
            })
            .expect("generate windows");
        let windows_count = engine
            .list_candidate_sets()
            .into_iter()
            .find(|s| s.name == "windows")
            .map(|s| s.candidate_count)
            .expect("windows set exists");

        engine
            .apply(Operation::FilterCandidateSet {
                input_set: "windows".to_string(),
                output_set: "windows_all_q".to_string(),
                metric: "gc_fraction".to_string(),
                min: None,
                max: None,
                min_quantile: Some(0.0),
                max_quantile: Some(1.0),
            })
            .expect("quantile full-range filter");
        let windows_all_q_count = engine
            .list_candidate_sets()
            .into_iter()
            .find(|s| s.name == "windows_all_q")
            .map(|s| s.candidate_count)
            .expect("windows_all_q exists");
        assert_eq!(windows_all_q_count, windows_count);

        engine
            .apply(Operation::FilterCandidateSet {
                input_set: "windows".to_string(),
                output_set: "windows_top_q".to_string(),
                metric: "gc_fraction".to_string(),
                min: None,
                max: None,
                min_quantile: Some(1.0),
                max_quantile: Some(1.0),
            })
            .expect("top quantile filter");
        let top_count = engine
            .list_candidate_sets()
            .into_iter()
            .find(|s| s.name == "windows_top_q")
            .map(|s| s.candidate_count)
            .expect("windows_top_q exists");
        assert!(top_count >= 1);
        assert!(top_count <= windows_count);

        let err = engine
            .apply(Operation::FilterCandidateSet {
                input_set: "windows".to_string(),
                output_set: "bad_q".to_string(),
                metric: "gc_fraction".to_string(),
                min: None,
                max: None,
                min_quantile: Some(1.1),
                max_quantile: None,
            })
            .unwrap_err();
        assert!(err.message.contains("between 0 and 1"));
    }

    #[test]
    fn test_extract_anchored_region_supports_middle_feature_boundary() {
        let mut state = ProjectState::default();
        let mut dna = DNAsequence::from_sequence("ACGTACGTACGTACGTACGT").expect("sequence");
        dna.features_mut().push(gb_io::seq::Feature {
            kind: gb_io::seq::FeatureKind::from("gene"),
            location: gb_io::seq::Location::simple_range(6, 14),
            qualifiers: vec![("label".into(), Some("MID_GENE".to_string()))],
        });
        state.sequences.insert("seqA".to_string(), dna);
        let mut engine = GentleEngine::from_state(state);

        let res = engine
            .apply(Operation::ExtractAnchoredRegion {
                input: "seqA".to_string(),
                anchor: SequenceAnchor::FeatureBoundary {
                    feature_kind: Some("gene".to_string()),
                    feature_label: Some("MID_GENE".to_string()),
                    boundary: AnchorBoundary::Middle,
                    occurrence: Some(0),
                },
                direction: AnchorDirection::Upstream,
                target_length_bp: 4,
                length_tolerance_bp: 0,
                required_re_sites: vec![],
                required_tf_motifs: vec![],
                forward_primer: None,
                reverse_primer: None,
                output_prefix: Some("mid".to_string()),
                unique: Some(true),
                max_candidates: Some(1),
            })
            .expect("extract using middle boundary");
        assert_eq!(res.created_seq_ids, vec!["mid_1".to_string()]);
        let seq = engine
            .state()
            .sequences
            .get("mid_1")
            .expect("created middle-boundary extract");
        assert_eq!(seq.len(), 4);
    }

    #[test]
    fn test_generate_candidate_set_between_two_sequence_anchors() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "seqA".to_string(),
            DNAsequence::from_sequence("ACGTACGTACGTACGT").unwrap(),
        );
        let mut engine = GentleEngine::from_state(state);

        let res = engine
            .apply(Operation::GenerateCandidateSetBetweenAnchors {
                set_name: "between".to_string(),
                seq_id: "seqA".to_string(),
                anchor_a: SequenceAnchor::Position { zero_based: 2 },
                anchor_b: SequenceAnchor::Position { zero_based: 10 },
                length_bp: 4,
                step_bp: 2,
                limit: Some(32),
            })
            .expect("generate between anchors");
        assert!(res.messages.iter().any(|m| m.contains("between anchors")));

        let summary = engine
            .list_candidate_sets()
            .into_iter()
            .find(|set| set.name == "between")
            .expect("between set summary");
        assert_eq!(summary.candidate_count, 3);

        let (page, total, _) = engine
            .inspect_candidate_set_page("between", 64, 0)
            .expect("inspect between candidates");
        assert_eq!(total, 3);
        assert_eq!(page.candidates[0].start_0based, 2);
        assert_eq!(page.candidates[0].end_0based, 6);
        assert_eq!(
            page.candidates[0]
                .metrics
                .get("distance_to_anchor_a_bp")
                .copied(),
            Some(0.0)
        );
        assert_eq!(
            page.candidates[0]
                .metrics
                .get("anchor_interval_span_bp")
                .copied(),
            Some(8.0)
        );
    }

    #[test]
    fn test_candidate_generation_feature_parts_ignores_multipart_gaps() {
        let mut state = ProjectState::default();
        let mut dna = DNAsequence::from_sequence("ACGTACGTACGTACGTACGT").expect("sequence");
        dna.features_mut().push(gb_io::seq::Feature {
            kind: gb_io::seq::FeatureKind::from("exon"),
            location: gb_io::seq::Location::Join(vec![
                gb_io::seq::Location::simple_range(2, 4),
                gb_io::seq::Location::simple_range(10, 12),
            ]),
            qualifiers: vec![("label".into(), Some("EXON_JOIN".to_string()))],
        });
        state.sequences.insert("seqA".to_string(), dna);
        let mut engine = GentleEngine::from_state(state);

        engine
            .apply(Operation::GenerateCandidateSet {
                set_name: "span_mode".to_string(),
                seq_id: "seqA".to_string(),
                length_bp: 1,
                step_bp: 1,
                feature_kinds: vec!["exon".to_string()],
                feature_label_regex: Some("^EXON_JOIN$".to_string()),
                max_distance_bp: Some(0),
                feature_geometry_mode: Some(CandidateFeatureGeometryMode::FeatureSpan),
                feature_boundary_mode: None,
                feature_strand_relation: None,
                limit: Some(256),
            })
            .expect("generate span-mode candidates");

        engine
            .apply(Operation::GenerateCandidateSet {
                set_name: "parts_mode".to_string(),
                seq_id: "seqA".to_string(),
                length_bp: 1,
                step_bp: 1,
                feature_kinds: vec!["exon".to_string()],
                feature_label_regex: Some("^EXON_JOIN$".to_string()),
                max_distance_bp: Some(0),
                feature_geometry_mode: Some(CandidateFeatureGeometryMode::FeatureParts),
                feature_boundary_mode: None,
                feature_strand_relation: None,
                limit: Some(256),
            })
            .expect("generate parts-mode candidates");

        let span_count = engine
            .list_candidate_sets()
            .into_iter()
            .find(|s| s.name == "span_mode")
            .map(|s| s.candidate_count)
            .expect("span_mode exists");
        let parts_count = engine
            .list_candidate_sets()
            .into_iter()
            .find(|s| s.name == "parts_mode")
            .map(|s| s.candidate_count)
            .expect("parts_mode exists");
        assert!(
            parts_count < span_count,
            "feature_parts should not fill multipart feature gaps"
        );
    }

    #[test]
    fn test_candidate_distance_feature_boundaries_respects_five_prime_and_three_prime() {
        let mut state = ProjectState::default();
        let mut dna = DNAsequence::from_sequence("ACGTACGTACGTACGTACGT").expect("sequence");
        dna.features_mut().push(gb_io::seq::Feature {
            kind: gb_io::seq::FeatureKind::from("gene"),
            location: gb_io::seq::Location::simple_range(5, 8),
            qualifiers: vec![("label".into(), Some("PLUS_GENE".to_string()))],
        });
        dna.features_mut().push(gb_io::seq::Feature {
            kind: gb_io::seq::FeatureKind::from("gene"),
            location: gb_io::seq::Location::Complement(Box::new(
                gb_io::seq::Location::simple_range(12, 15),
            )),
            qualifiers: vec![("label".into(), Some("MINUS_GENE".to_string()))],
        });
        state.sequences.insert("seqA".to_string(), dna);
        let mut engine = GentleEngine::from_state(state);
        engine
            .apply(Operation::GenerateCandidateSet {
                set_name: "windows".to_string(),
                seq_id: "seqA".to_string(),
                length_bp: 1,
                step_bp: 1,
                feature_kinds: vec![],
                feature_label_regex: None,
                max_distance_bp: None,
                feature_geometry_mode: None,
                feature_boundary_mode: None,
                feature_strand_relation: None,
                limit: Some(512),
            })
            .expect("generate windows");

        engine
            .apply(Operation::ScoreCandidateSetDistance {
                set_name: "windows".to_string(),
                metric: "dist_5p".to_string(),
                feature_kinds: vec!["gene".to_string()],
                feature_label_regex: None,
                feature_geometry_mode: Some(CandidateFeatureGeometryMode::FeatureBoundaries),
                feature_boundary_mode: Some(CandidateFeatureBoundaryMode::FivePrime),
                feature_strand_relation: None,
            })
            .expect("score five-prime distance");
        engine
            .apply(Operation::ScoreCandidateSetDistance {
                set_name: "windows".to_string(),
                metric: "dist_3p".to_string(),
                feature_kinds: vec!["gene".to_string()],
                feature_label_regex: None,
                feature_geometry_mode: Some(CandidateFeatureGeometryMode::FeatureBoundaries),
                feature_boundary_mode: Some(CandidateFeatureBoundaryMode::ThreePrime),
                feature_strand_relation: None,
            })
            .expect("score three-prime distance");

        let (page, _, _) = engine
            .inspect_candidate_set_page("windows", 2048, 0)
            .expect("inspect windows");
        let at_pos = |pos: usize, metric: &str| -> f64 {
            page.candidates
                .iter()
                .find(|candidate| candidate.start_0based == pos)
                .and_then(|candidate| candidate.metrics.get(metric).copied())
                .unwrap_or(f64::NAN)
        };
        let plus_start_dist_5p = at_pos(5, "dist_5p");
        let plus_end_dist_5p = at_pos(7, "dist_5p");
        let plus_start_dist_3p = at_pos(5, "dist_3p");
        let plus_end_dist_3p = at_pos(7, "dist_3p");
        assert!(plus_start_dist_5p < plus_end_dist_5p);
        assert!(plus_end_dist_3p < plus_start_dist_3p);

        let minus_end_dist_5p = at_pos(14, "dist_5p");
        let minus_start_dist_5p = at_pos(12, "dist_5p");
        let minus_end_dist_3p = at_pos(14, "dist_3p");
        let minus_start_dist_3p = at_pos(12, "dist_3p");
        assert!(minus_end_dist_5p < minus_start_dist_5p);
        assert!(minus_start_dist_3p < minus_end_dist_3p);
    }

    #[test]
    fn test_candidate_generation_feature_strand_relation_filters_plus_and_minus() {
        let mut state = ProjectState::default();
        let mut dna = DNAsequence::from_sequence("ACGTACGTACGTACGTACGT").expect("sequence");
        dna.features_mut().push(gb_io::seq::Feature {
            kind: gb_io::seq::FeatureKind::from("gene"),
            location: gb_io::seq::Location::simple_range(2, 5),
            qualifiers: vec![("label".into(), Some("PLUS_GENE".to_string()))],
        });
        dna.features_mut().push(gb_io::seq::Feature {
            kind: gb_io::seq::FeatureKind::from("gene"),
            location: gb_io::seq::Location::Complement(Box::new(
                gb_io::seq::Location::simple_range(12, 15),
            )),
            qualifiers: vec![("label".into(), Some("MINUS_GENE".to_string()))],
        });
        state.sequences.insert("seqA".to_string(), dna);
        let mut engine = GentleEngine::from_state(state);

        for (name, strand_relation) in [
            ("any", CandidateFeatureStrandRelation::Any),
            ("same", CandidateFeatureStrandRelation::Same),
            ("opposite", CandidateFeatureStrandRelation::Opposite),
        ] {
            engine
                .apply(Operation::GenerateCandidateSet {
                    set_name: format!("gene_{name}"),
                    seq_id: "seqA".to_string(),
                    length_bp: 1,
                    step_bp: 1,
                    feature_kinds: vec!["gene".to_string()],
                    feature_label_regex: None,
                    max_distance_bp: Some(0),
                    feature_geometry_mode: None,
                    feature_boundary_mode: None,
                    feature_strand_relation: Some(strand_relation),
                    limit: Some(256),
                })
                .expect("generate gene candidates with strand relation");
        }

        let count_for = |set_name: &str| -> usize {
            engine
                .list_candidate_sets()
                .into_iter()
                .find(|set| set.name == set_name)
                .map(|set| set.candidate_count)
                .expect("candidate set exists")
        };
        let any_count = count_for("gene_any");
        let same_count = count_for("gene_same");
        let opposite_count = count_for("gene_opposite");
        assert!(same_count > 0);
        assert!(opposite_count > 0);
        assert_eq!(same_count + opposite_count, any_count);
    }

    #[test]
    fn test_candidate_distance_feature_strand_relation_prefers_matching_strand_features() {
        let mut state = ProjectState::default();
        let mut dna = DNAsequence::from_sequence("ACGTACGTACGTACGTACGT").expect("sequence");
        dna.features_mut().push(gb_io::seq::Feature {
            kind: gb_io::seq::FeatureKind::from("gene"),
            location: gb_io::seq::Location::simple_range(2, 5),
            qualifiers: vec![("label".into(), Some("PLUS_GENE".to_string()))],
        });
        dna.features_mut().push(gb_io::seq::Feature {
            kind: gb_io::seq::FeatureKind::from("gene"),
            location: gb_io::seq::Location::Complement(Box::new(
                gb_io::seq::Location::simple_range(12, 15),
            )),
            qualifiers: vec![("label".into(), Some("MINUS_GENE".to_string()))],
        });
        state.sequences.insert("seqA".to_string(), dna);
        let mut engine = GentleEngine::from_state(state);

        engine
            .apply(Operation::GenerateCandidateSet {
                set_name: "windows".to_string(),
                seq_id: "seqA".to_string(),
                length_bp: 1,
                step_bp: 1,
                feature_kinds: vec!["gene".to_string()],
                feature_label_regex: None,
                max_distance_bp: None,
                feature_geometry_mode: None,
                feature_boundary_mode: None,
                feature_strand_relation: None,
                limit: Some(256),
            })
            .expect("generate candidate windows");

        engine
            .apply(Operation::ScoreCandidateSetDistance {
                set_name: "windows".to_string(),
                metric: "dist_any".to_string(),
                feature_kinds: vec!["gene".to_string()],
                feature_label_regex: None,
                feature_geometry_mode: None,
                feature_boundary_mode: None,
                feature_strand_relation: Some(CandidateFeatureStrandRelation::Any),
            })
            .expect("score distance any");
        engine
            .apply(Operation::ScoreCandidateSetDistance {
                set_name: "windows".to_string(),
                metric: "dist_same".to_string(),
                feature_kinds: vec!["gene".to_string()],
                feature_label_regex: None,
                feature_geometry_mode: None,
                feature_boundary_mode: None,
                feature_strand_relation: Some(CandidateFeatureStrandRelation::Same),
            })
            .expect("score distance same");
        engine
            .apply(Operation::ScoreCandidateSetDistance {
                set_name: "windows".to_string(),
                metric: "dist_opposite".to_string(),
                feature_kinds: vec!["gene".to_string()],
                feature_label_regex: None,
                feature_geometry_mode: None,
                feature_boundary_mode: None,
                feature_strand_relation: Some(CandidateFeatureStrandRelation::Opposite),
            })
            .expect("score distance opposite");

        let (page, _, _) = engine
            .inspect_candidate_set_page("windows", 4096, 0)
            .expect("inspect scored windows");
        let at_pos = |pos: usize, metric: &str| -> f64 {
            page.candidates
                .iter()
                .find(|candidate| candidate.start_0based == pos)
                .and_then(|candidate| candidate.metrics.get(metric).copied())
                .unwrap_or(f64::NAN)
        };

        let plus_any = at_pos(2, "dist_any");
        let plus_same = at_pos(2, "dist_same");
        let plus_opposite = at_pos(2, "dist_opposite");
        assert_eq!(plus_any, 0.0);
        assert_eq!(plus_same, 0.0);
        assert!(plus_opposite > 0.0);

        let minus_any = at_pos(14, "dist_any");
        let minus_same = at_pos(14, "dist_same");
        let minus_opposite = at_pos(14, "dist_opposite");
        assert_eq!(minus_any, 0.0);
        assert!(minus_same > 0.0);
        assert_eq!(minus_opposite, 0.0);
    }

    #[test]
    fn test_candidate_optimizer_weighted_topk_and_pareto() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "seqA".to_string(),
            DNAsequence::from_sequence("GCATGAAA").expect("sequence"),
        );
        let mut engine = GentleEngine::from_state(state);

        engine
            .apply(Operation::GenerateCandidateSet {
                set_name: "cand".to_string(),
                seq_id: "seqA".to_string(),
                length_bp: 2,
                step_bp: 2,
                feature_kinds: vec![],
                feature_label_regex: None,
                max_distance_bp: None,
                feature_geometry_mode: None,
                feature_boundary_mode: None,
                feature_strand_relation: None,
                limit: Some(64),
            })
            .expect("generate candidates");

        engine
            .apply(Operation::ScoreCandidateSetWeightedObjective {
                set_name: "cand".to_string(),
                metric: "objective".to_string(),
                objectives: vec![
                    CandidateWeightedObjectiveTerm {
                        metric: "gc_fraction".to_string(),
                        weight: 0.7,
                        direction: CandidateObjectiveDirection::Maximize,
                    },
                    CandidateWeightedObjectiveTerm {
                        metric: "distance_to_seq_start_bp".to_string(),
                        weight: 0.3,
                        direction: CandidateObjectiveDirection::Minimize,
                    },
                ],
                normalize_metrics: Some(true),
            })
            .expect("score weighted objective");

        engine
            .apply(Operation::TopKCandidateSet {
                input_set: "cand".to_string(),
                output_set: "cand_top1".to_string(),
                metric: "objective".to_string(),
                k: 1,
                direction: Some(CandidateObjectiveDirection::Maximize),
                tie_break: Some(CandidateTieBreakPolicy::SeqStartEnd),
            })
            .expect("top-k selection");
        let (top_page, _, _) = engine
            .inspect_candidate_set_page("cand_top1", 10, 0)
            .expect("inspect top-k result");
        assert_eq!(top_page.candidates.len(), 1);
        assert_eq!(top_page.candidates[0].start_0based, 0);

        engine
            .apply(Operation::ParetoFrontierCandidateSet {
                input_set: "cand".to_string(),
                output_set: "cand_pareto".to_string(),
                objectives: vec![
                    CandidateObjectiveSpec {
                        metric: "gc_fraction".to_string(),
                        direction: CandidateObjectiveDirection::Maximize,
                    },
                    CandidateObjectiveSpec {
                        metric: "distance_to_seq_start_bp".to_string(),
                        direction: CandidateObjectiveDirection::Minimize,
                    },
                ],
                max_candidates: None,
                tie_break: Some(CandidateTieBreakPolicy::SeqStartEnd),
            })
            .expect("pareto frontier selection");
        let pareto_summary = engine
            .list_candidate_sets()
            .into_iter()
            .find(|set| set.name == "cand_pareto")
            .expect("pareto set summary");
        assert_eq!(pareto_summary.candidate_count, 1);
    }

    #[test]
    fn test_topk_tie_break_policy_is_deterministic() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "seqA".to_string(),
            DNAsequence::from_sequence("ACGT").expect("sequence"),
        );
        let mut engine = GentleEngine::from_state(state);

        engine
            .apply(Operation::GenerateCandidateSet {
                set_name: "cand".to_string(),
                seq_id: "seqA".to_string(),
                length_bp: 2,
                step_bp: 1,
                feature_kinds: vec![],
                feature_label_regex: None,
                max_distance_bp: None,
                feature_geometry_mode: None,
                feature_boundary_mode: None,
                feature_strand_relation: None,
                limit: Some(64),
            })
            .expect("generate candidates");

        engine
            .apply(Operation::TopKCandidateSet {
                input_set: "cand".to_string(),
                output_set: "cand_top2".to_string(),
                metric: "length_bp".to_string(),
                k: 2,
                direction: Some(CandidateObjectiveDirection::Maximize),
                tie_break: Some(CandidateTieBreakPolicy::SeqStartEnd),
            })
            .expect("top-k tie-break");
        let (page, _, _) = engine
            .inspect_candidate_set_page("cand_top2", 10, 0)
            .expect("inspect top2");
        let starts = page
            .candidates
            .iter()
            .map(|candidate| candidate.start_0based)
            .collect::<Vec<_>>();
        assert_eq!(starts, vec![0, 1]);
    }

    #[test]
    fn test_candidate_macro_template_store_and_render() {
        let mut engine = GentleEngine::from_state(ProjectState::default());
        engine
            .apply(Operation::UpsertCandidateMacroTemplate {
                name: "scan_tp53".to_string(),
                description: Some("demo template".to_string()),
                details_url: Some("https://example.org/candidates/scan-tp53".to_string()),
                parameters: vec![
                    CandidateMacroTemplateParam {
                        name: "set_name".to_string(),
                        default_value: None,
                        required: true,
                    },
                    CandidateMacroTemplateParam {
                        name: "seq_id".to_string(),
                        default_value: Some("seqA".to_string()),
                        required: false,
                    },
                    CandidateMacroTemplateParam {
                        name: "len".to_string(),
                        default_value: Some("20".to_string()),
                        required: false,
                    },
                ],
                script: "generate ${set_name} ${seq_id} --length ${len} --step 1".to_string(),
            })
            .expect("upsert template");

        let templates = engine.list_candidate_macro_templates();
        assert_eq!(templates.len(), 1);
        assert_eq!(templates[0].name, "scan_tp53");
        assert_eq!(
            templates[0].details_url.as_deref(),
            Some("https://example.org/candidates/scan-tp53")
        );
        let template = engine
            .get_candidate_macro_template("scan_tp53")
            .expect("get candidate template");
        assert_eq!(
            template.details_url.as_deref(),
            Some("https://example.org/candidates/scan-tp53")
        );

        let rendered = engine
            .render_candidate_macro_template_script(
                "scan_tp53",
                &HashMap::from([("set_name".to_string(), "tp53_hits".to_string())]),
            )
            .expect("render template");
        assert_eq!(rendered, "generate tp53_hits seqA --length 20 --step 1");

        engine
            .apply(Operation::DeleteCandidateMacroTemplate {
                name: "scan_tp53".to_string(),
            })
            .expect("delete template");
        assert!(engine.list_candidate_macro_templates().is_empty());
    }

    #[test]
    fn test_candidate_macro_template_rejects_non_http_details_url() {
        let mut engine = GentleEngine::from_state(ProjectState::default());
        let err = engine
            .apply(Operation::UpsertCandidateMacroTemplate {
                name: "bad_url".to_string(),
                description: Some("template".to_string()),
                details_url: Some("ftp://example.org/template".to_string()),
                parameters: vec![CandidateMacroTemplateParam {
                    name: "set_name".to_string(),
                    default_value: None,
                    required: true,
                }],
                script: "generate ${set_name} seqA --length 20".to_string(),
            })
            .expect_err("expected invalid details_url error");
        assert!(
            err.message.contains("must start with http:// or https://"),
            "unexpected error: {}",
            err.message
        );
    }

    #[test]
    fn test_workflow_macro_template_store_and_render() {
        let mut engine = GentleEngine::from_state(ProjectState::default());
        engine
            .apply(Operation::UpsertWorkflowMacroTemplate {
                name: "clone_step".to_string(),
                description: Some("demo workflow template".to_string()),
                details_url: Some("https://example.org/cloning/clone-step".to_string()),
                parameters: vec![
                    WorkflowMacroTemplateParam {
                        name: "seq_id".to_string(),
                        default_value: None,
                        required: true,
                    },
                    WorkflowMacroTemplateParam {
                        name: "out_id".to_string(),
                        default_value: Some("seqA_rev".to_string()),
                        required: false,
                    },
                ],
                script: "op {\"Reverse\":{\"input\":\"${seq_id}\",\"output_id\":\"${out_id}\"}}"
                    .to_string(),
            })
            .expect("upsert workflow template");

        let templates = engine.list_workflow_macro_templates();
        assert_eq!(templates.len(), 1);
        assert_eq!(templates[0].name, "clone_step");
        let template = engine
            .get_workflow_macro_template("clone_step")
            .expect("get workflow template");
        assert_eq!(template.template_schema, CLONING_MACRO_TEMPLATE_SCHEMA);
        assert_eq!(
            template.details_url.as_deref(),
            Some("https://example.org/cloning/clone-step")
        );

        let rendered = engine
            .render_workflow_macro_template_script(
                "clone_step",
                &HashMap::from([("seq_id".to_string(), "seqX".to_string())]),
            )
            .expect("render workflow template");
        assert_eq!(
            rendered,
            "op {\"Reverse\":{\"input\":\"seqX\",\"output_id\":\"seqA_rev\"}}"
        );

        engine
            .apply(Operation::DeleteWorkflowMacroTemplate {
                name: "clone_step".to_string(),
            })
            .expect("delete workflow template");
        assert!(engine.list_workflow_macro_templates().is_empty());
    }

    #[test]
    fn test_workflow_macro_template_rejects_non_http_details_url() {
        let mut engine = GentleEngine::from_state(ProjectState::default());
        let err = engine
            .apply(Operation::UpsertWorkflowMacroTemplate {
                name: "bad_url".to_string(),
                description: Some("template".to_string()),
                details_url: Some("ftp://example.org/template".to_string()),
                parameters: vec![WorkflowMacroTemplateParam {
                    name: "seq_id".to_string(),
                    default_value: None,
                    required: true,
                }],
                script: "op {\"Reverse\":{\"input\":\"${seq_id}\"}}".to_string(),
            })
            .expect_err("expected invalid details_url error");
        assert!(
            err.message.contains("must start with http:// or https://"),
            "unexpected error: {}",
            err.message
        );
    }

    #[test]
    fn test_guide_design_filter_generate_and_export() {
        let mut engine = GentleEngine::from_state(ProjectState::default());
        engine
            .apply(Operation::UpsertGuideSet {
                guide_set_id: "tp73_guides".to_string(),
                guides: vec![
                    GuideCandidate {
                        guide_id: "g1".to_string(),
                        seq_id: "tp73".to_string(),
                        start_0based: 100,
                        end_0based_exclusive: 120,
                        strand: "+".to_string(),
                        protospacer: "GACCTGTTGACGATGTTCCA".to_string(),
                        pam: "AGG".to_string(),
                        nuclease: "SpCas9".to_string(),
                        cut_offset_from_protospacer_start: 17,
                        rank: Some(1),
                    },
                    GuideCandidate {
                        guide_id: "g2".to_string(),
                        seq_id: "tp73".to_string(),
                        start_0based: 200,
                        end_0based_exclusive: 220,
                        strand: "+".to_string(),
                        protospacer: "ATATATATATATATATATAT".to_string(),
                        pam: "AGG".to_string(),
                        nuclease: "SpCas9".to_string(),
                        cut_offset_from_protospacer_start: 17,
                        rank: Some(2),
                    },
                    GuideCandidate {
                        guide_id: "g3".to_string(),
                        seq_id: "tp73".to_string(),
                        start_0based: 300,
                        end_0based_exclusive: 320,
                        strand: "+".to_string(),
                        protospacer: "GACTTTTGACTGACTGACT".to_string(),
                        pam: "AGG".to_string(),
                        nuclease: "SpCas9".to_string(),
                        cut_offset_from_protospacer_start: 17,
                        rank: Some(3),
                    },
                ],
            })
            .expect("upsert guide set");

        engine
            .apply(Operation::FilterGuidesPractical {
                guide_set_id: "tp73_guides".to_string(),
                config: GuidePracticalFilterConfig {
                    gc_min: Some(0.30),
                    gc_max: Some(0.70),
                    max_homopolymer_run: Some(4),
                    max_homopolymer_run_per_base: HashMap::new(),
                    reject_ambiguous_bases: true,
                    avoid_u6_terminator_tttt: true,
                    u6_terminator_window: GuideU6TerminatorWindow::SpacerPlusTail,
                    max_dinucleotide_repeat_units: None,
                    forbidden_motifs: vec![],
                    required_5prime_base: None,
                    allow_5prime_g_extension: true,
                },
                output_guide_set_id: Some("tp73_guides_pass".to_string()),
            })
            .expect("filter guides");
        let report = engine
            .get_guide_practical_filter_report("tp73_guides")
            .expect("guide report");
        assert_eq!(report.passed_count, 1);
        assert_eq!(report.rejected_count, 2);

        engine
            .apply(Operation::GenerateGuideOligos {
                guide_set_id: "tp73_guides".to_string(),
                template_id: "lenti_bsmbi_u6_default".to_string(),
                apply_5prime_g_extension: Some(true),
                output_oligo_set_id: Some("tp73_oligos".to_string()),
                passed_only: Some(true),
            })
            .expect("generate guide oligos");

        let dir = tempdir().expect("tempdir");
        let csv_path = dir.path().join("guides.csv").display().to_string();
        let fasta_path = dir.path().join("guides.fa").display().to_string();
        let protocol_path = dir.path().join("guides.protocol.txt").display().to_string();

        engine
            .apply(Operation::ExportGuideOligos {
                guide_set_id: "tp73_guides".to_string(),
                oligo_set_id: Some("tp73_oligos".to_string()),
                format: GuideOligoExportFormat::CsvTable,
                path: csv_path.clone(),
                plate_format: None,
            })
            .expect("export guide oligos csv");
        engine
            .apply(Operation::ExportGuideOligos {
                guide_set_id: "tp73_guides".to_string(),
                oligo_set_id: Some("tp73_oligos".to_string()),
                format: GuideOligoExportFormat::Fasta,
                path: fasta_path.clone(),
                plate_format: None,
            })
            .expect("export guide oligos fasta");
        engine
            .apply(Operation::ExportGuideProtocolText {
                guide_set_id: "tp73_guides".to_string(),
                oligo_set_id: Some("tp73_oligos".to_string()),
                path: protocol_path.clone(),
                include_qc_checklist: Some(true),
            })
            .expect("export guide protocol");

        let csv = fs::read_to_string(csv_path).expect("read csv export");
        assert!(csv.contains("guide_id,rank,forward_oligo,reverse_oligo,notes"));
        assert!(csv.contains("g1"));
        let fasta = fs::read_to_string(fasta_path).expect("read fasta export");
        assert!(fasta.contains(">g1|forward"));
        let protocol = fs::read_to_string(protocol_path).expect("read protocol export");
        assert!(protocol.contains("GENtle Guide Oligo Protocol"));
    }

    #[test]
    fn test_guide_set_duplicate_ids_rejected() {
        let mut engine = GentleEngine::from_state(ProjectState::default());
        let err = engine
            .apply(Operation::UpsertGuideSet {
                guide_set_id: "dup".to_string(),
                guides: vec![
                    GuideCandidate {
                        guide_id: "dup_1".to_string(),
                        seq_id: "s1".to_string(),
                        start_0based: 0,
                        end_0based_exclusive: 20,
                        strand: "+".to_string(),
                        protospacer: "GACCTGTTGACGATGTTCCA".to_string(),
                        pam: "AGG".to_string(),
                        nuclease: "SpCas9".to_string(),
                        cut_offset_from_protospacer_start: 17,
                        rank: Some(1),
                    },
                    GuideCandidate {
                        guide_id: "dup_1".to_string(),
                        seq_id: "s1".to_string(),
                        start_0based: 20,
                        end_0based_exclusive: 40,
                        strand: "+".to_string(),
                        protospacer: "GACCTGTTGACGATGTTCCC".to_string(),
                        pam: "AGG".to_string(),
                        nuclease: "SpCas9".to_string(),
                        cut_offset_from_protospacer_start: 17,
                        rank: Some(2),
                    },
                ],
            })
            .expect_err("duplicate guide ids should be rejected");
        assert!(err.message.contains("Duplicate guide_id"));
    }

    #[test]
    fn test_candidate_metric_expression_fuzz_smoke_does_not_panic() {
        let metrics = HashMap::from([
            ("gc_fraction".to_string(), 0.5),
            ("at_fraction".to_string(), 0.5),
            ("length_bp".to_string(), 20.0),
            ("candidate_index".to_string(), 3.0),
            ("distance_to_feature_start_bp".to_string(), 7.0),
        ]);
        let atoms = [
            "gc_fraction",
            "at_fraction",
            "length_bp",
            "candidate_index",
            "distance_to_feature_start_bp",
            "1",
            "2.5",
        ];
        let ops = ["+", "-", "*", "/"];

        for i in 0..512usize {
            let left = atoms[i % atoms.len()];
            let mid = atoms[(i / 7) % atoms.len()];
            let right = atoms[(i / 29) % atoms.len()];
            let op_a = ops[i % ops.len()];
            let op_b = ops[(i / 11) % ops.len()];
            let expr = format!("({left} {op_a} ({mid} {op_b} {right}))");
            let run = std::panic::catch_unwind(|| {
                if let Ok(parsed) = GentleEngine::parse_metric_expression(&expr) {
                    let _ = GentleEngine::evaluate_metric_expression(&parsed, &metrics);
                }
            });
            assert!(run.is_ok(), "expression caused panic: {expr}");
        }
    }
}
