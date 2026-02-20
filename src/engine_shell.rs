use crate::{
    dna_ladder::LadderMolecule,
    dna_sequence::DNAsequence,
    engine::{
        CandidateFeatureBoundaryMode, CandidateFeatureGeometryMode, CandidateFeatureStrandRelation,
        Engine, GenomeAnchorSide, GenomeTrackSource, GenomeTrackSubscription, GentleEngine,
        Operation, ProjectState, RenderSvgMode, SequenceAnchor, Workflow,
    },
    genomes::{GenomeGeneRecord, DEFAULT_GENOME_CATALOG_PATH, DEFAULT_HELPER_GENOME_CATALOG_PATH},
    resource_sync,
    shell_docs::{
        shell_help_json as render_shell_help_json,
        shell_help_markdown as render_shell_help_markdown,
        shell_help_text as render_shell_help_text,
        shell_topic_help_json as render_shell_topic_help_json,
        shell_topic_help_markdown as render_shell_topic_help_markdown,
        shell_topic_help_text as render_shell_topic_help_text, HelpOutputFormat,
    },
    tf_motifs,
};
#[cfg(all(target_os = "macos", feature = "screenshot-capture"))]
use objc2_app_kit::NSApplication;
#[cfg(all(target_os = "macos", feature = "screenshot-capture"))]
use objc2_foundation::MainThreadMarker;
use regex::{Regex, RegexBuilder};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
#[cfg(all(target_os = "macos", feature = "screenshot-capture"))]
use std::path::Path;
#[cfg(all(target_os = "macos", feature = "screenshot-capture"))]
use std::process::Command;
use std::{collections::BTreeSet, fs};

#[derive(Debug, Clone)]
pub enum ShellCommand {
    Help {
        topic: Vec<String>,
        format: HelpOutputFormat,
        interface_filter: Option<String>,
    },
    Capabilities,
    StateSummary,
    LoadProject {
        path: String,
    },
    SaveProject {
        path: String,
    },
    ScreenshotWindow {
        output: String,
    },
    RenderSvg {
        seq_id: String,
        mode: RenderSvgMode,
        output: String,
    },
    RenderRnaSvg {
        seq_id: String,
        output: String,
    },
    RnaInfo {
        seq_id: String,
    },
    RenderLineageSvg {
        output: String,
    },
    RenderPoolGelSvg {
        inputs: Vec<String>,
        output: String,
        ladders: Option<Vec<String>>,
    },
    LaddersList {
        molecule: LadderMolecule,
        name_filter: Option<String>,
    },
    LaddersExport {
        molecule: LadderMolecule,
        output: String,
        name_filter: Option<String>,
    },
    ExportPool {
        inputs: Vec<String>,
        output: String,
        human_id: Option<String>,
    },
    ImportPool {
        input: String,
        prefix: String,
    },
    ResourcesSyncRebase {
        input: String,
        output: Option<String>,
        commercial_only: bool,
    },
    ResourcesSyncJaspar {
        input: String,
        output: Option<String>,
    },
    ReferenceList {
        helper_mode: bool,
        catalog_path: Option<String>,
    },
    ReferenceValidateCatalog {
        helper_mode: bool,
        catalog_path: Option<String>,
    },
    ReferenceStatus {
        helper_mode: bool,
        genome_id: String,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
    },
    ReferenceGenes {
        helper_mode: bool,
        genome_id: String,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
        filter: String,
        biotypes: Vec<String>,
        limit: usize,
        offset: usize,
    },
    ReferencePrepare {
        helper_mode: bool,
        genome_id: String,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
        timeout_seconds: Option<u64>,
    },
    ReferenceExtractRegion {
        helper_mode: bool,
        genome_id: String,
        chromosome: String,
        start_1based: usize,
        end_1based: usize,
        output_id: Option<String>,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
    },
    ReferenceExtractGene {
        helper_mode: bool,
        genome_id: String,
        gene_query: String,
        occurrence: Option<usize>,
        output_id: Option<String>,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
    },
    ReferenceExtendAnchor {
        helper_mode: bool,
        seq_id: String,
        side: GenomeAnchorSide,
        length_bp: usize,
        output_id: Option<String>,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
    },
    ReferenceBlast {
        helper_mode: bool,
        genome_id: String,
        query_sequence: String,
        max_hits: usize,
        task: Option<String>,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
    },
    ReferenceBlastTrack {
        helper_mode: bool,
        genome_id: String,
        query_sequence: String,
        target_seq_id: String,
        max_hits: usize,
        task: Option<String>,
        track_name: Option<String>,
        clear_existing: bool,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
    },
    TracksImportBed {
        seq_id: String,
        path: String,
        track_name: Option<String>,
        min_score: Option<f64>,
        max_score: Option<f64>,
        clear_existing: bool,
    },
    TracksImportBigWig {
        seq_id: String,
        path: String,
        track_name: Option<String>,
        min_score: Option<f64>,
        max_score: Option<f64>,
        clear_existing: bool,
    },
    TracksImportVcf {
        seq_id: String,
        path: String,
        track_name: Option<String>,
        min_score: Option<f64>,
        max_score: Option<f64>,
        clear_existing: bool,
    },
    TracksTrackedList,
    TracksTrackedAdd {
        subscription: GenomeTrackSubscription,
    },
    TracksTrackedRemove {
        index: usize,
    },
    TracksTrackedClear,
    TracksTrackedApply {
        index: Option<usize>,
        only_new_anchors: bool,
    },
    CandidatesList,
    CandidatesDelete {
        set_name: String,
    },
    CandidatesGenerate {
        set_name: String,
        seq_id: String,
        length_bp: usize,
        step_bp: usize,
        feature_kinds: Vec<String>,
        feature_label_regex: Option<String>,
        max_distance_bp: Option<usize>,
        feature_geometry_mode: Option<CandidateFeatureGeometryMode>,
        feature_boundary_mode: Option<CandidateFeatureBoundaryMode>,
        feature_strand_relation: Option<CandidateFeatureStrandRelation>,
        limit: usize,
    },
    CandidatesGenerateBetweenAnchors {
        set_name: String,
        seq_id: String,
        anchor_a: SequenceAnchor,
        anchor_b: SequenceAnchor,
        length_bp: usize,
        step_bp: usize,
        limit: usize,
    },
    CandidatesShow {
        set_name: String,
        limit: usize,
        offset: usize,
    },
    CandidatesMetrics {
        set_name: String,
    },
    CandidatesScoreExpression {
        set_name: String,
        metric: String,
        expression: String,
    },
    CandidatesScoreDistance {
        set_name: String,
        metric: String,
        feature_kinds: Vec<String>,
        feature_label_regex: Option<String>,
        feature_geometry_mode: Option<CandidateFeatureGeometryMode>,
        feature_boundary_mode: Option<CandidateFeatureBoundaryMode>,
        feature_strand_relation: Option<CandidateFeatureStrandRelation>,
    },
    CandidatesFilter {
        input_set: String,
        output_set: String,
        metric: String,
        min: Option<f64>,
        max: Option<f64>,
        min_quantile: Option<f64>,
        max_quantile: Option<f64>,
    },
    CandidatesSetOp {
        op: CandidateSetOperator,
        left_set: String,
        right_set: String,
        output_set: String,
    },
    CandidatesMacro {
        script: String,
        transactional: bool,
    },
    SetParameter {
        name: String,
        value_json: String,
    },
    Op {
        payload: String,
    },
    Workflow {
        payload: String,
    },
}

const SCREENSHOT_DISABLED_MESSAGE: &str = "screenshot-window is disabled by security policy";

#[derive(Debug, Clone)]
pub struct ShellRunResult {
    pub state_changed: bool,
    pub output: Value,
}

#[derive(Debug, Clone, Copy, Default)]
pub struct ShellExecutionOptions {
    pub allow_screenshots: bool,
}

impl ShellExecutionOptions {
    pub fn from_env() -> Self {
        let raw = std::env::var("GENTLE_ALLOW_SCREENSHOTS").unwrap_or_default();
        let allow_screenshots = cfg!(feature = "screenshot-capture")
            && matches!(
                raw.trim().to_ascii_lowercase().as_str(),
                "1" | "true" | "yes" | "on"
            );
        Self { allow_screenshots }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[allow(dead_code)]
struct ScreenshotReport {
    schema: String,
    path: String,
    window_title: String,
    captured_at_unix_ms: u128,
    pixel_width: u32,
    pixel_height: u32,
    backend: String,
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
    seq_id: String,
    #[serde(default)]
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
    #[serde(default)]
    human_id: String,
    member_count: usize,
    members: Vec<PoolMember>,
}

const CANDIDATE_SETS_METADATA_KEY: &str = "candidate_sets";
const DEFAULT_CANDIDATE_PAGE_SIZE: usize = 100;
const DEFAULT_CANDIDATE_SET_LIMIT: usize = 50_000;

#[derive(Debug, Clone, Copy)]
pub enum CandidateSetOperator {
    Union,
    Intersect,
    Subtract,
}

impl CandidateSetOperator {
    fn parse(raw: &str) -> Option<Self> {
        match raw.trim().to_ascii_lowercase().as_str() {
            "union" => Some(Self::Union),
            "intersect" | "intersection" => Some(Self::Intersect),
            "subtract" | "difference" => Some(Self::Subtract),
            _ => None,
        }
    }

    fn as_str(self) -> &'static str {
        match self {
            Self::Union => "union",
            Self::Intersect => "intersect",
            Self::Subtract => "subtract",
        }
    }
}

fn parse_candidate_feature_geometry_mode(
    raw: &str,
) -> Result<CandidateFeatureGeometryMode, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "feature_span" | "span" => Ok(CandidateFeatureGeometryMode::FeatureSpan),
        "feature_parts" | "parts" | "segments" => Ok(CandidateFeatureGeometryMode::FeatureParts),
        "feature_boundaries" | "boundaries" | "boundary" => {
            Ok(CandidateFeatureGeometryMode::FeatureBoundaries)
        }
        other => Err(format!(
            "Unsupported feature geometry mode '{other}' (expected feature_span|feature_parts|feature_boundaries)"
        )),
    }
}

fn parse_candidate_feature_boundary_mode(
    raw: &str,
) -> Result<CandidateFeatureBoundaryMode, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "any" => Ok(CandidateFeatureBoundaryMode::Any),
        "five_prime" | "five-prime" | "5p" | "5'" => Ok(CandidateFeatureBoundaryMode::FivePrime),
        "three_prime" | "three-prime" | "3p" | "3'" => {
            Ok(CandidateFeatureBoundaryMode::ThreePrime)
        }
        "start" => Ok(CandidateFeatureBoundaryMode::Start),
        "end" => Ok(CandidateFeatureBoundaryMode::End),
        other => Err(format!(
            "Unsupported feature boundary mode '{other}' (expected any|five_prime|three_prime|start|end)"
        )),
    }
}

fn parse_candidate_feature_strand_relation(
    raw: &str,
) -> Result<CandidateFeatureStrandRelation, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "any" => Ok(CandidateFeatureStrandRelation::Any),
        "same" => Ok(CandidateFeatureStrandRelation::Same),
        "opposite" => Ok(CandidateFeatureStrandRelation::Opposite),
        other => Err(format!(
            "Unsupported strand relation '{other}' (expected any|same|opposite)"
        )),
    }
}

fn parse_sequence_anchor_json(raw: &str, option_name: &str) -> Result<SequenceAnchor, String> {
    serde_json::from_str::<SequenceAnchor>(raw)
        .map_err(|e| format!("Invalid {} JSON '{}': {}", option_name, raw, e))
}

impl ShellCommand {
    pub fn preview(&self) -> String {
        match self {
            Self::Help {
                topic,
                format,
                interface_filter,
            } => {
                let topic_label = if topic.is_empty() {
                    "all".to_string()
                } else {
                    topic.join(" ")
                };
                let interface = interface_filter
                    .as_deref()
                    .filter(|v| !v.trim().is_empty())
                    .unwrap_or("all");
                format!(
                    "show shell command help (topic='{topic_label}', format={}, interface={interface})",
                    format.as_str()
                )
            }
            Self::Capabilities => "inspect engine capabilities".to_string(),
            Self::StateSummary => "show sequence/container state summary".to_string(),
            Self::LoadProject { path } => format!("load project state from '{path}'"),
            Self::SaveProject { path } => format!("save current project state to '{path}'"),
            Self::ScreenshotWindow { output } => {
                format!("capture active GENtle window screenshot to '{output}'")
            }
            Self::RenderSvg {
                seq_id,
                mode,
                output,
            } => format!("render {mode:?} SVG for '{seq_id}' to '{output}'"),
            Self::RenderRnaSvg { seq_id, output } => {
                format!("render RNA structure SVG for '{seq_id}' to '{output}'")
            }
            Self::RnaInfo { seq_id } => {
                format!("inspect rnapkin textual RNA report for '{seq_id}'")
            }
            Self::RenderLineageSvg { output } => format!("render lineage SVG to '{output}'"),
            Self::RenderPoolGelSvg {
                inputs,
                output,
                ladders,
            } => {
                let ladders = ladders
                    .as_ref()
                    .map(|v| v.join(","))
                    .unwrap_or_else(|| "auto".to_string());
                format!(
                    "render pool gel SVG from {} input(s) to '{output}' with ladders {ladders}",
                    inputs.len()
                )
            }
            Self::LaddersList {
                molecule,
                name_filter,
            } => {
                let filter = name_filter
                    .as_deref()
                    .map(str::trim)
                    .filter(|v| !v.is_empty())
                    .unwrap_or("-");
                format!(
                    "inspect {} ladders (filter={filter})",
                    molecule.display_name()
                )
            }
            Self::LaddersExport {
                molecule,
                output,
                name_filter,
            } => {
                let filter = name_filter
                    .as_deref()
                    .map(str::trim)
                    .filter(|v| !v.is_empty())
                    .unwrap_or("-");
                format!(
                    "export {} ladders to '{output}' (filter={filter})",
                    molecule.display_name()
                )
            }
            Self::ExportPool {
                inputs,
                output,
                human_id,
            } => {
                let human = human_id.clone().unwrap_or_else(|| "-".to_string());
                format!(
                    "export pool with {} input(s) to '{output}' (human_id={human})",
                    inputs.len()
                )
            }
            Self::ImportPool { input, prefix } => {
                format!("import pool from '{input}' with prefix '{prefix}'")
            }
            Self::ResourcesSyncRebase {
                input,
                output,
                commercial_only,
            } => {
                let output = output
                    .clone()
                    .unwrap_or_else(|| resource_sync::DEFAULT_REBASE_RESOURCE_PATH.to_string());
                format!(
                    "sync REBASE from '{input}' to '{output}'{}",
                    if *commercial_only {
                        " (commercial-only)"
                    } else {
                        ""
                    }
                )
            }
            Self::ResourcesSyncJaspar { input, output } => {
                let output = output
                    .clone()
                    .unwrap_or_else(|| resource_sync::DEFAULT_JASPAR_RESOURCE_PATH.to_string());
                format!("sync JASPAR from '{input}' to '{output}'")
            }
            Self::ReferenceList {
                helper_mode,
                catalog_path,
            } => {
                let label = if *helper_mode { "helpers" } else { "genomes" };
                let catalog = catalog_path
                    .clone()
                    .unwrap_or_else(|| default_catalog_path(*helper_mode).to_string());
                format!("list {label} from catalog '{catalog}'")
            }
            Self::ReferenceValidateCatalog {
                helper_mode,
                catalog_path,
            } => {
                let label = if *helper_mode { "helpers" } else { "genomes" };
                let catalog = catalog_path
                    .clone()
                    .unwrap_or_else(|| default_catalog_path(*helper_mode).to_string());
                format!("validate {label} catalog '{catalog}'")
            }
            Self::ReferenceStatus {
                helper_mode,
                genome_id,
                catalog_path,
                cache_dir,
            } => {
                let label = if *helper_mode { "helper" } else { "genome" };
                let catalog = catalog_path
                    .clone()
                    .unwrap_or_else(|| default_catalog_path(*helper_mode).to_string());
                let cache = cache_dir.clone().unwrap_or_else(|| "-".to_string());
                format!("check {label} '{genome_id}' status (catalog='{catalog}', cache='{cache}')")
            }
            Self::ReferenceGenes {
                helper_mode,
                genome_id,
                catalog_path,
                cache_dir,
                filter,
                biotypes,
                limit,
                offset,
            } => {
                let label = if *helper_mode { "helper" } else { "genome" };
                let catalog = catalog_path
                    .clone()
                    .unwrap_or_else(|| default_catalog_path(*helper_mode).to_string());
                let cache = cache_dir.clone().unwrap_or_else(|| "-".to_string());
                let biotypes = if biotypes.is_empty() {
                    "-".to_string()
                } else {
                    biotypes.join(",")
                };
                format!(
                    "list {label} genes for '{genome_id}' (catalog='{catalog}', cache='{cache}', filter='{filter}', biotypes='{biotypes}', limit={limit}, offset={offset})"
                )
            }
            Self::ReferencePrepare {
                helper_mode,
                genome_id,
                catalog_path,
                cache_dir,
                timeout_seconds,
            } => {
                let label = if *helper_mode { "helper" } else { "genome" };
                let catalog = catalog_path
                    .clone()
                    .unwrap_or_else(|| default_catalog_path(*helper_mode).to_string());
                let cache = cache_dir.clone().unwrap_or_else(|| "-".to_string());
                let timeout = timeout_seconds
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string());
                format!(
                    "prepare {label} '{genome_id}' (catalog='{catalog}', cache='{cache}', timeout='{timeout}')"
                )
            }
            Self::ReferenceExtractRegion {
                helper_mode,
                genome_id,
                chromosome,
                start_1based,
                end_1based,
                output_id,
                catalog_path,
                cache_dir,
            } => {
                let label = if *helper_mode { "helper" } else { "genome" };
                let catalog = catalog_path
                    .clone()
                    .unwrap_or_else(|| default_catalog_path(*helper_mode).to_string());
                let cache = cache_dir.clone().unwrap_or_else(|| "-".to_string());
                let output = output_id.clone().unwrap_or_else(|| "-".to_string());
                format!(
                    "extract {label} region {genome_id}:{chromosome}:{start_1based}-{end_1based} (output='{output}', catalog='{catalog}', cache='{cache}')"
                )
            }
            Self::ReferenceExtractGene {
                helper_mode,
                genome_id,
                gene_query,
                occurrence,
                output_id,
                catalog_path,
                cache_dir,
            } => {
                let label = if *helper_mode { "helper" } else { "genome" };
                let catalog = catalog_path
                    .clone()
                    .unwrap_or_else(|| default_catalog_path(*helper_mode).to_string());
                let cache = cache_dir.clone().unwrap_or_else(|| "-".to_string());
                let occ = occurrence
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string());
                let output = output_id.clone().unwrap_or_else(|| "-".to_string());
                format!(
                    "extract {label} gene '{gene_query}' from '{genome_id}' (occurrence={occ}, output='{output}', catalog='{catalog}', cache='{cache}')"
                )
            }
            Self::ReferenceExtendAnchor {
                helper_mode,
                seq_id,
                side,
                length_bp,
                output_id,
                catalog_path,
                cache_dir,
            } => {
                let label = if *helper_mode { "helper" } else { "genome" };
                let catalog = catalog_path
                    .clone()
                    .unwrap_or_else(|| default_catalog_path(*helper_mode).to_string());
                let cache = cache_dir.clone().unwrap_or_else(|| "-".to_string());
                let output = output_id.clone().unwrap_or_else(|| "-".to_string());
                let side_label = match side {
                    GenomeAnchorSide::FivePrime => "5'",
                    GenomeAnchorSide::ThreePrime => "3'",
                };
                format!(
                    "extend {label}-anchored sequence '{seq_id}' on {side_label} by {length_bp} bp (output='{output}', catalog='{catalog}', cache='{cache}')"
                )
            }
            Self::ReferenceBlast {
                helper_mode,
                genome_id,
                query_sequence,
                max_hits,
                task,
                catalog_path,
                cache_dir,
            } => {
                let label = if *helper_mode { "helper" } else { "genome" };
                let catalog = catalog_path
                    .clone()
                    .unwrap_or_else(|| default_catalog_path(*helper_mode).to_string());
                let cache = cache_dir.clone().unwrap_or_else(|| "-".to_string());
                let task = task.clone().unwrap_or_else(|| "blastn-short".to_string());
                format!(
                    "blast query (len={}) against {label} '{genome_id}' (max_hits={max_hits}, task='{task}', catalog='{catalog}', cache='{cache}')",
                    query_sequence.len()
                )
            }
            Self::ReferenceBlastTrack {
                helper_mode,
                genome_id,
                query_sequence,
                target_seq_id,
                max_hits,
                task,
                track_name,
                clear_existing,
                catalog_path,
                cache_dir,
            } => {
                let label = if *helper_mode { "helper" } else { "genome" };
                let catalog = catalog_path
                    .clone()
                    .unwrap_or_else(|| default_catalog_path(*helper_mode).to_string());
                let cache = cache_dir.clone().unwrap_or_else(|| "-".to_string());
                let task = task.clone().unwrap_or_else(|| "blastn-short".to_string());
                format!(
                    "blast query (len={}) against {label} '{genome_id}' and import hits to '{}' (max_hits={max_hits}, task='{}', track_name='{}', clear_existing={}, catalog='{}', cache='{}')",
                    query_sequence.len(),
                    target_seq_id,
                    task,
                    track_name.clone().unwrap_or_else(|| "blast_hits".to_string()),
                    clear_existing,
                    catalog,
                    cache
                )
            }
            Self::TracksImportBed {
                seq_id,
                path,
                track_name,
                min_score,
                max_score,
                clear_existing,
            } => format!(
                "import BED track for '{seq_id}' from '{}' (track_name='{}', min_score={}, max_score={}, clear_existing={})",
                path,
                track_name
                    .clone()
                    .unwrap_or_else(|| "-".to_string()),
                min_score
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                max_score
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                clear_existing
            ),
            Self::TracksImportBigWig {
                seq_id,
                path,
                track_name,
                min_score,
                max_score,
                clear_existing,
            } => format!(
                "import BigWig track for '{seq_id}' from '{}' (track_name='{}', min_score={}, max_score={}, clear_existing={})",
                path,
                track_name
                    .clone()
                    .unwrap_or_else(|| "-".to_string()),
                min_score
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                max_score
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                clear_existing
            ),
            Self::TracksImportVcf {
                seq_id,
                path,
                track_name,
                min_score,
                max_score,
                clear_existing,
            } => format!(
                "import VCF track for '{seq_id}' from '{}' (track_name='{}', min_score={}, max_score={}, clear_existing={})",
                path,
                track_name
                    .clone()
                    .unwrap_or_else(|| "-".to_string()),
                min_score
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                max_score
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                clear_existing
            ),
            Self::TracksTrackedList => "list tracked genome signal files".to_string(),
            Self::TracksTrackedAdd { subscription } => format!(
                "add tracked {} file '{}' (track_name='{}', min_score={}, max_score={}, clear_existing={})",
                subscription.source.label(),
                subscription.path,
                subscription
                    .track_name
                    .clone()
                    .unwrap_or_else(|| "-".to_string()),
                subscription
                    .min_score
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                subscription
                    .max_score
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                subscription.clear_existing
            ),
            Self::TracksTrackedRemove { index } => {
                format!("remove tracked genome signal file at index {}", index)
            }
            Self::TracksTrackedClear => "clear all tracked genome signal files".to_string(),
            Self::TracksTrackedApply {
                index,
                only_new_anchors,
            } => match index {
                Some(idx) => format!(
                    "apply tracked genome signal file index {} to all anchored sequences (only_new_anchors={})",
                    idx, only_new_anchors
                ),
                None => format!(
                    "apply tracked genome signal files (only_new_anchors={})",
                    only_new_anchors
                ),
            },
            Self::CandidatesList => "list persisted candidate sets".to_string(),
            Self::CandidatesDelete { set_name } => {
                format!("delete candidate set '{set_name}'")
            }
            Self::CandidatesGenerate {
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
                let kinds = if feature_kinds.is_empty() {
                    "-".to_string()
                } else {
                    feature_kinds.join(",")
                };
                let label = feature_label_regex
                    .as_deref()
                    .filter(|v| !v.trim().is_empty())
                    .unwrap_or("-");
                let distance = max_distance_bp
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string());
                let geometry = feature_geometry_mode
                    .map(|mode| mode.as_str())
                    .unwrap_or("-");
                let boundary = feature_boundary_mode
                    .map(|mode| mode.as_str())
                    .unwrap_or("-");
                let strand_relation = feature_strand_relation
                    .map(|mode| mode.as_str())
                    .unwrap_or("-");
                format!(
                    "generate candidate set '{set_name}' from '{seq_id}' (length={length_bp}, step={step_bp}, feature_kinds={kinds}, feature_label_regex={label}, max_distance={distance}, feature_geometry_mode={geometry}, feature_boundary_mode={boundary}, feature_strand_relation={strand_relation}, limit={limit})"
                )
            }
            Self::CandidatesGenerateBetweenAnchors {
                set_name,
                seq_id,
                anchor_a,
                anchor_b,
                length_bp,
                step_bp,
                limit,
            } => format!(
                "generate candidate set '{set_name}' from '{seq_id}' between anchors A={anchor_a:?} and B={anchor_b:?} (length={length_bp}, step={step_bp}, limit={limit})"
            ),
            Self::CandidatesShow {
                set_name,
                limit,
                offset,
            } => format!("show candidates in '{set_name}' (limit={limit}, offset={offset})"),
            Self::CandidatesMetrics { set_name } => {
                format!("list available metrics for candidate set '{set_name}'")
            }
            Self::CandidatesScoreExpression {
                set_name,
                metric,
                expression,
            } => format!(
                "compute derived metric '{metric}' for candidate set '{set_name}' using expression '{expression}'"
            ),
            Self::CandidatesScoreDistance {
                set_name,
                metric,
                feature_kinds,
                feature_label_regex,
                feature_geometry_mode,
                feature_boundary_mode,
                feature_strand_relation,
            } => {
                let kinds = if feature_kinds.is_empty() {
                    "-".to_string()
                } else {
                    feature_kinds.join(",")
                };
                let label = feature_label_regex
                    .as_deref()
                    .filter(|v| !v.trim().is_empty())
                    .unwrap_or("-");
                let geometry = feature_geometry_mode
                    .map(|mode| mode.as_str())
                    .unwrap_or("-");
                let boundary = feature_boundary_mode
                    .map(|mode| mode.as_str())
                    .unwrap_or("-");
                let strand_relation = feature_strand_relation
                    .map(|mode| mode.as_str())
                    .unwrap_or("-");
                format!(
                    "compute distance metric '{metric}' for candidate set '{set_name}' (feature_kinds={kinds}, feature_label_regex={label}, feature_geometry_mode={geometry}, feature_boundary_mode={boundary}, feature_strand_relation={strand_relation})"
                )
            }
            Self::CandidatesFilter {
                input_set,
                output_set,
                metric,
                min,
                max,
                min_quantile,
                max_quantile,
            } => format!(
                "filter candidate set '{input_set}' by metric '{metric}' into '{output_set}' (min={}, max={}, min_quantile={}, max_quantile={})",
                min.map(|v| v.to_string()).unwrap_or_else(|| "-".to_string()),
                max.map(|v| v.to_string()).unwrap_or_else(|| "-".to_string()),
                min_quantile
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                max_quantile
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string())
            ),
            Self::CandidatesSetOp {
                op,
                left_set,
                right_set,
                output_set,
            } => format!(
                "candidate set {} '{}' and '{}' into '{}'",
                op.as_str(),
                left_set,
                right_set,
                output_set
            ),
            Self::CandidatesMacro {
                script,
                transactional,
            } => {
                let trimmed = script.trim();
                let preview = if trimmed.starts_with('@') {
                    trimmed.to_string()
                } else {
                    let single_line = trimmed.replace('\n', " ");
                    if single_line.len() > 80 {
                        format!("{}...", &single_line[..80])
                    } else {
                        single_line
                    }
                };
                let mode = if *transactional {
                    "transactional"
                } else {
                    "best-effort"
                };
                format!("run {mode} candidates macro '{preview}'")
            }
            Self::SetParameter { name, value_json } => {
                format!("set parameter '{}' to {}", name, value_json)
            }
            Self::Op { .. } => "apply one engine operation from JSON".to_string(),
            Self::Workflow { .. } => "apply engine workflow from JSON".to_string(),
        }
    }

    pub fn is_state_mutating(&self) -> bool {
        matches!(
            self,
            Self::LoadProject { .. }
                | Self::ImportPool { .. }
                | Self::ReferencePrepare { .. }
                | Self::ReferenceExtractRegion { .. }
                | Self::ReferenceExtractGene { .. }
                | Self::ReferenceExtendAnchor { .. }
                | Self::ReferenceBlastTrack { .. }
                | Self::TracksImportBed { .. }
                | Self::TracksImportBigWig { .. }
                | Self::TracksImportVcf { .. }
                | Self::TracksTrackedAdd { .. }
                | Self::TracksTrackedRemove { .. }
                | Self::TracksTrackedClear
                | Self::TracksTrackedApply { .. }
                | Self::CandidatesDelete { .. }
                | Self::CandidatesGenerate { .. }
                | Self::CandidatesGenerateBetweenAnchors { .. }
                | Self::CandidatesScoreExpression { .. }
                | Self::CandidatesScoreDistance { .. }
                | Self::CandidatesFilter { .. }
                | Self::CandidatesSetOp { .. }
                | Self::CandidatesMacro { .. }
                | Self::SetParameter { .. }
                | Self::Op { .. }
                | Self::Workflow { .. }
        )
    }
}

pub fn shell_help_text() -> String {
    render_shell_help_text(None).unwrap_or_else(|e| format!("Could not render shell help: {e}"))
}

fn split_ids(input: &str) -> Vec<String> {
    input
        .split(',')
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .map(|s| s.to_string())
        .collect()
}

fn unique_id(existing: &std::collections::HashMap<String, DNAsequence>, base: &str) -> String {
    if !existing.contains_key(base) {
        return base.to_string();
    }
    let mut i = 2usize;
    loop {
        let candidate = format!("{base}_{i}");
        if !existing.contains_key(&candidate) {
            return candidate;
        }
        i += 1;
    }
}

fn apply_member_overhang(member: &PoolMember, dna: &mut DNAsequence) -> Result<(), String> {
    let mut value =
        serde_json::to_value(&*dna).map_err(|e| format!("Could not serialize sequence: {e}"))?;
    let obj = value
        .as_object_mut()
        .ok_or_else(|| "Internal serialization shape error".to_string())?;
    obj.insert(
        "overhang".to_string(),
        json!({
            "forward_3": member.ends.forward_3.as_bytes(),
            "forward_5": member.ends.forward_5.as_bytes(),
            "reverse_3": member.ends.reverse_3.as_bytes(),
            "reverse_5": member.ends.reverse_5.as_bytes(),
        }),
    );
    let patched: DNAsequence =
        serde_json::from_value(value).map_err(|e| format!("Could not restore overhang: {e}"))?;
    *dna = patched;
    Ok(())
}

fn default_catalog_path(helper_mode: bool) -> &'static str {
    if helper_mode {
        DEFAULT_HELPER_GENOME_CATALOG_PATH
    } else {
        DEFAULT_GENOME_CATALOG_PATH
    }
}

fn operation_catalog_path(catalog_path: &Option<String>, helper_mode: bool) -> Option<String> {
    catalog_path
        .as_deref()
        .map(str::trim)
        .filter(|v| !v.is_empty())
        .map(|v| v.to_string())
        .or_else(|| helper_mode.then(|| default_catalog_path(helper_mode).to_string()))
}

fn resolved_catalog_path<'a>(
    catalog_path: &'a Option<String>,
    helper_mode: bool,
) -> Option<&'a str> {
    if let Some(value) = catalog_path
        .as_deref()
        .map(str::trim)
        .filter(|v| !v.is_empty())
    {
        Some(value)
    } else if helper_mode {
        Some(default_catalog_path(helper_mode))
    } else {
        None
    }
}

fn effective_catalog_path(catalog_path: &Option<String>, helper_mode: bool) -> String {
    catalog_path
        .as_deref()
        .map(str::trim)
        .filter(|v| !v.is_empty())
        .map(|v| v.to_string())
        .unwrap_or_else(|| default_catalog_path(helper_mode).to_string())
}

fn compile_gene_filter_regex(filter: &str) -> Result<Option<Regex>, String> {
    let pattern = filter.trim();
    if pattern.is_empty() {
        return Ok(None);
    }
    RegexBuilder::new(pattern)
        .case_insensitive(true)
        .build()
        .map(Some)
        .map_err(|e| format!("Invalid --filter regex '{}': {}", pattern, e))
}

fn genome_gene_matches_regex_filter(gene: &GenomeGeneRecord, regex: &Regex) -> bool {
    gene.gene_name
        .as_ref()
        .map(|name| regex.is_match(name))
        .unwrap_or(false)
        || gene
            .gene_id
            .as_ref()
            .map(|id| regex.is_match(id))
            .unwrap_or(false)
        || regex.is_match(&gene.chromosome)
}

fn collect_biotypes(genes: &[GenomeGeneRecord]) -> Vec<String> {
    let mut biotypes: BTreeSet<String> = BTreeSet::new();
    for gene in genes {
        let Some(biotype) = gene.biotype.as_ref() else {
            continue;
        };
        let trimmed = biotype.trim();
        if !trimmed.is_empty() {
            biotypes.insert(trimmed.to_string());
        }
    }
    biotypes.into_iter().collect()
}

fn genome_gene_matches_filter(
    gene: &GenomeGeneRecord,
    filter_regex: Option<&Regex>,
    allowed_biotypes_lower: &[String],
) -> bool {
    let regex_ok = filter_regex
        .map(|re| genome_gene_matches_regex_filter(gene, re))
        .unwrap_or(true);
    if !regex_ok {
        return false;
    }
    if allowed_biotypes_lower.is_empty() {
        return true;
    }
    gene.biotype
        .as_ref()
        .map(|b| b.trim().to_ascii_lowercase())
        .map(|probe| allowed_biotypes_lower.iter().any(|b| b == &probe))
        .unwrap_or(false)
}

fn parse_option_path(
    tokens: &[String],
    idx: &mut usize,
    option_name: &str,
    context: &str,
) -> Result<String, String> {
    if *idx + 1 >= tokens.len() {
        return Err(format!("Missing value after {option_name} for {context}"));
    }
    let value = tokens[*idx + 1].clone();
    *idx += 2;
    Ok(value)
}

fn parse_mode(mode: &str) -> Result<RenderSvgMode, String> {
    match mode {
        "linear" => Ok(RenderSvgMode::Linear),
        "circular" => Ok(RenderSvgMode::Circular),
        other => Err(format!(
            "Unknown render mode '{other}', expected 'linear' or 'circular'"
        )),
    }
}

fn parse_ladder_molecule(value: &str) -> Result<LadderMolecule, String> {
    LadderMolecule::parse(value)
        .ok_or_else(|| format!("Unknown ladder molecule '{value}', expected 'dna' or 'rna'"))
}

fn parse_anchor_side(value: &str) -> Result<GenomeAnchorSide, String> {
    match value.trim().to_ascii_lowercase().as_str() {
        "5" | "5p" | "5prime" | "5'" | "five_prime" | "five-prime" => {
            Ok(GenomeAnchorSide::FivePrime)
        }
        "3" | "3p" | "3prime" | "3'" | "three_prime" | "three-prime" => {
            Ok(GenomeAnchorSide::ThreePrime)
        }
        _ => Err(format!(
            "Unknown anchor side '{}'; expected 5p or 3p",
            value
        )),
    }
}

#[cfg(all(target_os = "macos", feature = "screenshot-capture"))]
#[allow(dead_code)]
fn now_unix_ms() -> u128 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_millis())
        .unwrap_or(0)
}

#[cfg(all(target_os = "macos", feature = "screenshot-capture"))]
#[allow(dead_code)]
fn ensure_parent_dir(path: &str) -> Result<(), String> {
    let parent = Path::new(path)
        .parent()
        .map(|p| p.to_path_buf())
        .unwrap_or_else(|| Path::new(".").to_path_buf());
    fs::create_dir_all(parent)
        .map_err(|e| format!("Could not create output directory for '{path}': {e}"))
}

#[cfg(all(target_os = "macos", feature = "screenshot-capture"))]
#[allow(dead_code)]
fn active_window_info_from_appkit() -> Result<(u64, String), String> {
    let Some(mtm) = MainThreadMarker::new() else {
        return Err("Screenshot capture requires the macOS main thread".to_string());
    };
    let app = NSApplication::sharedApplication(mtm);
    let Some(window) = app.keyWindow().or_else(|| unsafe { app.mainWindow() }) else {
        return Err(
            "No active GENtle window in this process; run from GUI shell with the target window focused"
                .to_string(),
        );
    };
    let raw_window_id = unsafe { window.windowNumber() };
    if raw_window_id <= 0 {
        return Err(format!(
            "Active GENtle window reported invalid window number {raw_window_id}"
        ));
    }
    let window_id = u64::try_from(raw_window_id)
        .map_err(|_| format!("Could not convert window number {raw_window_id} to u64"))?;
    let title = window.title().to_string();
    let window_title = if title.trim().is_empty() {
        "<untitled>".to_string()
    } else {
        title.trim().to_string()
    };
    Ok((window_id, window_title))
}

#[cfg(all(target_os = "macos", feature = "screenshot-capture"))]
#[allow(dead_code)]
fn capture_active_window_screenshot(output: &str) -> Result<ScreenshotReport, String> {
    ensure_parent_dir(output)?;
    let (window_id, window_title) = active_window_info_from_appkit()?;

    let status = Command::new("screencapture")
        .arg("-x")
        .arg("-l")
        .arg(window_id.to_string())
        .arg(output)
        .status()
        .map_err(|e| format!("Could not run macOS screencapture command: {e}"))?;
    if !status.success() {
        return Err(format!(
            "macOS screencapture failed with status {:?}",
            status.code()
        ));
    }
    let (pixel_width, pixel_height) = image::image_dimensions(output)
        .map_err(|e| format!("Screenshot written but dimensions could not be read: {e}"))?;

    Ok(ScreenshotReport {
        schema: "gentle.screenshot.v1".to_string(),
        path: output.to_string(),
        window_title,
        captured_at_unix_ms: now_unix_ms(),
        pixel_width,
        pixel_height,
        backend: "macos.screencapture".to_string(),
    })
}

#[cfg(any(not(target_os = "macos"), not(feature = "screenshot-capture")))]
#[allow(dead_code)]
fn capture_active_window_screenshot(_output: &str) -> Result<ScreenshotReport, String> {
    if !cfg!(feature = "screenshot-capture") {
        Err(
            "screenshot-window is unavailable in this build; enable feature 'screenshot-capture'"
                .to_string(),
        )
    } else {
        Err("screenshot-window is currently supported only on macOS".to_string())
    }
}

fn parse_json_payload(raw: &str) -> Result<String, String> {
    if let Some(path) = raw.strip_prefix('@') {
        fs::read_to_string(path).map_err(|e| format!("Could not read JSON file '{path}': {e}"))
    } else {
        Ok(raw.to_string())
    }
}

fn token_error(command: &str) -> String {
    format!("Invalid '{command}' usage. Try: help")
}

fn parse_help_command(tokens: &[String]) -> Result<ShellCommand, String> {
    let mut topic = vec![];
    let mut format = HelpOutputFormat::Text;
    let mut interface_filter: Option<String> = None;
    let mut idx = 1usize;
    while idx < tokens.len() {
        match tokens[idx].as_str() {
            "--format" => {
                let raw = parse_option_path(tokens, &mut idx, "--format", "help")?;
                format = HelpOutputFormat::parse(&raw)?;
            }
            "--interface" => {
                let raw = parse_option_path(tokens, &mut idx, "--interface", "help")?;
                interface_filter = Some(raw);
            }
            other if other.starts_with("--") => {
                return Err(format!("Unknown option '{other}' for help"));
            }
            value => {
                topic.push(value.to_string());
                idx += 1;
            }
        }
    }
    Ok(ShellCommand::Help {
        topic,
        format,
        interface_filter,
    })
}

fn parse_reference_command(tokens: &[String], helper_mode: bool) -> Result<ShellCommand, String> {
    let label = if helper_mode { "helpers" } else { "genomes" };
    if tokens.len() < 2 {
        return Err(format!("{label} requires a subcommand"));
    }
    match tokens[1].as_str() {
        "list" => {
            let mut catalog_path: Option<String> = None;
            let mut idx = 2usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--catalog" => {
                        catalog_path = Some(parse_option_path(tokens, &mut idx, "--catalog", label)?)
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for {label} list"));
                    }
                }
            }
            Ok(ShellCommand::ReferenceList {
                helper_mode,
                catalog_path,
            })
        }
        "validate-catalog" => {
            let mut catalog_path: Option<String> = None;
            let mut idx = 2usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--catalog" => {
                        catalog_path = Some(parse_option_path(tokens, &mut idx, "--catalog", label)?)
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for {label} validate-catalog"
                        ));
                    }
                }
            }
            Ok(ShellCommand::ReferenceValidateCatalog {
                helper_mode,
                catalog_path,
            })
        }
        "status" => {
            if tokens.len() < 3 {
                return Err(format!(
                    "{label} status requires GENOME_ID [--catalog PATH] [--cache-dir PATH]"
                ));
            }
            let genome_id = tokens[2].clone();
            let mut catalog_path: Option<String> = None;
            let mut cache_dir: Option<String> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--catalog" => {
                        catalog_path =
                            Some(parse_option_path(tokens, &mut idx, "--catalog", label)?)
                    }
                    "--cache-dir" => {
                        cache_dir = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--cache-dir",
                            label,
                        )?)
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for {label} status"));
                    }
                }
            }
            Ok(ShellCommand::ReferenceStatus {
                helper_mode,
                genome_id,
                catalog_path,
                cache_dir,
            })
        }
        "genes" => {
            if tokens.len() < 3 {
                return Err(format!(
                    "{label} genes requires GENOME_ID [--catalog PATH] [--cache-dir PATH] [--filter REGEX] [--biotype NAME] [--limit N] [--offset N]"
                ));
            }
            let genome_id = tokens[2].clone();
            let mut catalog_path: Option<String> = None;
            let mut cache_dir: Option<String> = None;
            let mut filter = String::new();
            let mut biotypes: Vec<String> = vec![];
            let mut limit: usize = 200;
            let mut offset: usize = 0;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--catalog" => {
                        catalog_path =
                            Some(parse_option_path(tokens, &mut idx, "--catalog", label)?)
                    }
                    "--cache-dir" => {
                        cache_dir = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--cache-dir",
                            label,
                        )?)
                    }
                    "--filter" => {
                        filter = parse_option_path(tokens, &mut idx, "--filter", label)?
                    }
                    "--biotype" => {
                        let biotype = parse_option_path(tokens, &mut idx, "--biotype", label)?;
                        let trimmed = biotype.trim();
                        if !trimmed.is_empty() {
                            biotypes.push(trimmed.to_string());
                        }
                    }
                    "--limit" => {
                        let raw = parse_option_path(tokens, &mut idx, "--limit", label)?;
                        limit = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --limit value '{raw}': {e}"))?;
                        if limit == 0 {
                            return Err("--limit must be >= 1".to_string());
                        }
                    }
                    "--offset" => {
                        let raw = parse_option_path(tokens, &mut idx, "--offset", label)?;
                        offset = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --offset value '{raw}': {e}"))?;
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for {label} genes"));
                    }
                }
            }
            Ok(ShellCommand::ReferenceGenes {
                helper_mode,
                genome_id,
                catalog_path,
                cache_dir,
                filter,
                biotypes,
                limit,
                offset,
            })
        }
        "prepare" => {
            if tokens.len() < 3 {
                return Err(format!(
                    "{label} prepare requires GENOME_ID [--catalog PATH] [--cache-dir PATH] [--timeout-secs N]"
                ));
            }
            let genome_id = tokens[2].clone();
            let mut catalog_path: Option<String> = None;
            let mut cache_dir: Option<String> = None;
            let mut timeout_seconds: Option<u64> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--catalog" => {
                        catalog_path =
                            Some(parse_option_path(tokens, &mut idx, "--catalog", label)?)
                    }
                    "--cache-dir" => {
                        cache_dir = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--cache-dir",
                            label,
                        )?)
                    }
                    "--timeout-secs" => {
                        let raw = parse_option_path(tokens, &mut idx, "--timeout-secs", label)?;
                        let parsed = raw
                            .parse::<u64>()
                            .map_err(|e| format!("Invalid --timeout-secs value '{raw}': {e}"))?;
                        if parsed == 0 {
                            timeout_seconds = None;
                        } else {
                            timeout_seconds = Some(parsed);
                        }
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for {label} prepare"));
                    }
                }
            }
            Ok(ShellCommand::ReferencePrepare {
                helper_mode,
                genome_id,
                catalog_path,
                cache_dir,
                timeout_seconds,
            })
        }
        "blast" => {
            if tokens.len() < 4 {
                return Err(format!(
                    "{label} blast requires GENOME_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--catalog PATH] [--cache-dir PATH]"
                ));
            }
            let genome_id = tokens[2].clone();
            let query_sequence = tokens[3].clone();
            let mut max_hits: usize = 25;
            let mut task: Option<String> = None;
            let mut catalog_path: Option<String> = None;
            let mut cache_dir: Option<String> = None;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--max-hits" => {
                        let raw = parse_option_path(tokens, &mut idx, "--max-hits", label)?;
                        max_hits = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --max-hits value '{raw}': {e}"))?;
                        if max_hits == 0 {
                            return Err("--max-hits must be >= 1".to_string());
                        }
                    }
                    "--task" => {
                        let raw = parse_option_path(tokens, &mut idx, "--task", label)?;
                        let normalized = raw.trim().to_ascii_lowercase();
                        match normalized.as_str() {
                            "blastn-short" | "blastn" => task = Some(normalized),
                            _ => {
                                return Err(format!(
                                    "Unsupported --task value '{}'; expected blastn-short or blastn",
                                    raw
                                ))
                            }
                        }
                    }
                    "--catalog" => {
                        catalog_path =
                            Some(parse_option_path(tokens, &mut idx, "--catalog", label)?)
                    }
                    "--cache-dir" => {
                        cache_dir = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--cache-dir",
                            label,
                        )?)
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for {label} blast"));
                    }
                }
            }
            Ok(ShellCommand::ReferenceBlast {
                helper_mode,
                genome_id,
                query_sequence,
                max_hits,
                task,
                catalog_path,
                cache_dir,
            })
        }
        "blast-track" => {
            if tokens.len() < 5 {
                return Err(format!(
                    "{label} blast-track requires GENOME_ID QUERY_SEQUENCE TARGET_SEQ_ID [--max-hits N] [--task blastn-short|blastn] [--track-name NAME] [--clear-existing] [--catalog PATH] [--cache-dir PATH]"
                ));
            }
            let genome_id = tokens[2].clone();
            let query_sequence = tokens[3].clone();
            let target_seq_id = tokens[4].clone();
            let mut max_hits: usize = 25;
            let mut task: Option<String> = None;
            let mut track_name: Option<String> = None;
            let mut clear_existing = false;
            let mut catalog_path: Option<String> = None;
            let mut cache_dir: Option<String> = None;
            let mut idx = 5usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--max-hits" => {
                        let raw = parse_option_path(tokens, &mut idx, "--max-hits", label)?;
                        max_hits = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --max-hits value '{raw}': {e}"))?;
                        if max_hits == 0 {
                            return Err("--max-hits must be >= 1".to_string());
                        }
                    }
                    "--task" => {
                        let raw = parse_option_path(tokens, &mut idx, "--task", label)?;
                        let normalized = raw.trim().to_ascii_lowercase();
                        match normalized.as_str() {
                            "blastn-short" | "blastn" => task = Some(normalized),
                            _ => {
                                return Err(format!(
                                    "Unsupported --task value '{}'; expected blastn-short or blastn",
                                    raw
                                ))
                            }
                        }
                    }
                    "--track-name" => {
                        track_name =
                            Some(parse_option_path(tokens, &mut idx, "--track-name", label)?)
                    }
                    "--clear-existing" => {
                        clear_existing = true;
                        idx += 1;
                    }
                    "--catalog" => {
                        catalog_path =
                            Some(parse_option_path(tokens, &mut idx, "--catalog", label)?)
                    }
                    "--cache-dir" => {
                        cache_dir = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--cache-dir",
                            label,
                        )?)
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for {label} blast-track"));
                    }
                }
            }
            Ok(ShellCommand::ReferenceBlastTrack {
                helper_mode,
                genome_id,
                query_sequence,
                target_seq_id,
                max_hits,
                task,
                track_name,
                clear_existing,
                catalog_path,
                cache_dir,
            })
        }
        "extract-region" => {
            if tokens.len() < 6 {
                return Err(format!(
                    "{label} extract-region requires GENOME_ID CHR START END [--output-id ID] [--catalog PATH] [--cache-dir PATH]"
                ));
            }
            let genome_id = tokens[2].clone();
            let chromosome = tokens[3].clone();
            let start_1based = tokens[4]
                .parse::<usize>()
                .map_err(|e| format!("Invalid START coordinate '{}': {e}", tokens[4]))?;
            let end_1based = tokens[5]
                .parse::<usize>()
                .map_err(|e| format!("Invalid END coordinate '{}': {e}", tokens[5]))?;
            let mut output_id: Option<String> = None;
            let mut catalog_path: Option<String> = None;
            let mut cache_dir: Option<String> = None;
            let mut idx = 6usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--output-id" => {
                        output_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--output-id",
                            label,
                        )?)
                    }
                    "--catalog" => {
                        catalog_path =
                            Some(parse_option_path(tokens, &mut idx, "--catalog", label)?)
                    }
                    "--cache-dir" => {
                        cache_dir = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--cache-dir",
                            label,
                        )?)
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for {label} extract-region"
                        ));
                    }
                }
            }
            Ok(ShellCommand::ReferenceExtractRegion {
                helper_mode,
                genome_id,
                chromosome,
                start_1based,
                end_1based,
                output_id,
                catalog_path,
                cache_dir,
            })
        }
        "extract-gene" => {
            if tokens.len() < 4 {
                return Err(format!(
                    "{label} extract-gene requires GENOME_ID QUERY [--occurrence N] [--output-id ID] [--catalog PATH] [--cache-dir PATH]"
                ));
            }
            let genome_id = tokens[2].clone();
            let gene_query = tokens[3].clone();
            let mut occurrence: Option<usize> = None;
            let mut output_id: Option<String> = None;
            let mut catalog_path: Option<String> = None;
            let mut cache_dir: Option<String> = None;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--occurrence" => {
                        let raw = parse_option_path(tokens, &mut idx, "--occurrence", label)?;
                        let value = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --occurrence value '{raw}': {e}"))?;
                        if value == 0 {
                            return Err("--occurrence must be >= 1".to_string());
                        }
                        occurrence = Some(value);
                    }
                    "--output-id" => {
                        output_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--output-id",
                            label,
                        )?)
                    }
                    "--catalog" => {
                        catalog_path =
                            Some(parse_option_path(tokens, &mut idx, "--catalog", label)?)
                    }
                    "--cache-dir" => {
                        cache_dir = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--cache-dir",
                            label,
                        )?)
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for {label} extract-gene"));
                    }
                }
            }
            Ok(ShellCommand::ReferenceExtractGene {
                helper_mode,
                genome_id,
                gene_query,
                occurrence,
                output_id,
                catalog_path,
                cache_dir,
            })
        }
        "extend-anchor" => {
            if tokens.len() < 5 {
                return Err(format!(
                    "{label} extend-anchor requires SEQ_ID 5p|3p LENGTH_BP [--output-id ID] [--catalog PATH] [--cache-dir PATH]"
                ));
            }
            let seq_id = tokens[2].clone();
            let side = parse_anchor_side(&tokens[3])?;
            let length_bp = tokens[4]
                .parse::<usize>()
                .map_err(|e| format!("Invalid LENGTH_BP '{}': {e}", tokens[4]))?;
            if length_bp == 0 {
                return Err("LENGTH_BP must be >= 1".to_string());
            }
            let mut output_id: Option<String> = None;
            let mut catalog_path: Option<String> = None;
            let mut cache_dir: Option<String> = None;
            let mut idx = 5usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--output-id" => {
                        output_id = Some(parse_option_path(tokens, &mut idx, "--output-id", label)?)
                    }
                    "--catalog" => {
                        catalog_path =
                            Some(parse_option_path(tokens, &mut idx, "--catalog", label)?)
                    }
                    "--cache-dir" => {
                        cache_dir = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--cache-dir",
                            label,
                        )?)
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for {label} extend-anchor"
                        ));
                    }
                }
            }
            Ok(ShellCommand::ReferenceExtendAnchor {
                helper_mode,
                seq_id,
                side,
                length_bp,
                output_id,
                catalog_path,
                cache_dir,
            })
        }
        other => Err(format!(
            "Unknown {label} subcommand '{other}' (expected list, validate-catalog, status, genes, prepare, blast, blast-track, extract-region, extract-gene, extend-anchor)"
        )),
    }
}

fn parse_candidates_command(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.len() < 2 {
        return Err(
            "candidates requires a subcommand: list, delete, generate, generate-between-anchors, show, metrics, score, score-distance, filter, set-op, macro"
                .to_string(),
        );
    }
    match tokens[1].as_str() {
        "list" => {
            if tokens.len() > 2 {
                return Err("candidates list takes no options".to_string());
            }
            Ok(ShellCommand::CandidatesList)
        }
        "delete" => {
            if tokens.len() != 3 {
                return Err("candidates delete requires SET_NAME".to_string());
            }
            Ok(ShellCommand::CandidatesDelete {
                set_name: tokens[2].clone(),
            })
        }
        "generate" => {
            if tokens.len() < 4 {
                return Err(
                    "candidates generate requires SET_NAME SEQ_ID --length N [--step N] [--feature-kind KIND] [--feature-label-regex REGEX] [--max-distance N] [--feature-geometry MODE] [--feature-boundary MODE] [--strand-relation MODE] [--limit N]"
                        .to_string(),
                );
            }
            let set_name = tokens[2].clone();
            let seq_id = tokens[3].clone();
            let mut length_bp: Option<usize> = None;
            let mut step_bp: usize = 1;
            let mut feature_kinds: Vec<String> = vec![];
            let mut feature_label_regex: Option<String> = None;
            let mut max_distance_bp: Option<usize> = None;
            let mut feature_geometry_mode: Option<CandidateFeatureGeometryMode> = None;
            let mut feature_boundary_mode: Option<CandidateFeatureBoundaryMode> = None;
            let mut feature_strand_relation: Option<CandidateFeatureStrandRelation> = None;
            let mut limit: usize = DEFAULT_CANDIDATE_SET_LIMIT;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--length" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--length", "candidates generate")?;
                        let parsed = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --length value '{raw}': {e}"))?;
                        if parsed == 0 {
                            return Err("--length must be >= 1".to_string());
                        }
                        length_bp = Some(parsed);
                    }
                    "--step" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--step", "candidates generate")?;
                        step_bp = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --step value '{raw}': {e}"))?;
                        if step_bp == 0 {
                            return Err("--step must be >= 1".to_string());
                        }
                    }
                    "--feature-kind" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--feature-kind",
                            "candidates generate",
                        )?;
                        let trimmed = raw.trim();
                        if !trimmed.is_empty() {
                            feature_kinds.push(trimmed.to_string());
                        }
                    }
                    "--feature-label-regex" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--feature-label-regex",
                            "candidates generate",
                        )?;
                        if raw.trim().is_empty() {
                            feature_label_regex = None;
                        } else {
                            feature_label_regex = Some(raw);
                        }
                    }
                    "--max-distance" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--max-distance",
                            "candidates generate",
                        )?;
                        max_distance_bp = Some(
                            raw.parse::<usize>()
                                .map_err(|e| format!("Invalid --max-distance value '{raw}': {e}"))?,
                        );
                    }
                    "--feature-geometry" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--feature-geometry",
                            "candidates generate",
                        )?;
                        feature_geometry_mode = Some(parse_candidate_feature_geometry_mode(&raw)?);
                    }
                    "--feature-boundary" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--feature-boundary",
                            "candidates generate",
                        )?;
                        feature_boundary_mode = Some(parse_candidate_feature_boundary_mode(&raw)?);
                    }
                    "--strand-relation" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--strand-relation",
                            "candidates generate",
                        )?;
                        feature_strand_relation =
                            Some(parse_candidate_feature_strand_relation(&raw)?);
                    }
                    "--limit" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--limit", "candidates generate")?;
                        limit = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --limit value '{raw}': {e}"))?;
                        if limit == 0 {
                            return Err("--limit must be >= 1".to_string());
                        }
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for candidates generate"
                        ));
                    }
                }
            }
            let Some(length_bp) = length_bp else {
                return Err("candidates generate requires --length N".to_string());
            };
            Ok(ShellCommand::CandidatesGenerate {
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
            })
        }
        "generate-between-anchors" | "generate-between" => {
            if tokens.len() < 4 {
                return Err(
                    "candidates generate-between-anchors requires SET_NAME SEQ_ID --length N (--anchor-a-pos N | --anchor-a-json JSON) (--anchor-b-pos N | --anchor-b-json JSON) [--step N] [--limit N]"
                        .to_string(),
                );
            }
            let set_name = tokens[2].clone();
            let seq_id = tokens[3].clone();
            let mut length_bp: Option<usize> = None;
            let mut step_bp: usize = 1;
            let mut limit: usize = DEFAULT_CANDIDATE_SET_LIMIT;
            let mut anchor_a: Option<SequenceAnchor> = None;
            let mut anchor_b: Option<SequenceAnchor> = None;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--length" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--length",
                            "candidates generate-between-anchors",
                        )?;
                        let parsed = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --length value '{raw}': {e}"))?;
                        if parsed == 0 {
                            return Err("--length must be >= 1".to_string());
                        }
                        length_bp = Some(parsed);
                    }
                    "--step" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--step",
                            "candidates generate-between-anchors",
                        )?;
                        step_bp = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --step value '{raw}': {e}"))?;
                        if step_bp == 0 {
                            return Err("--step must be >= 1".to_string());
                        }
                    }
                    "--limit" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--limit",
                            "candidates generate-between-anchors",
                        )?;
                        limit = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --limit value '{raw}': {e}"))?;
                        if limit == 0 {
                            return Err("--limit must be >= 1".to_string());
                        }
                    }
                    "--anchor-a-pos" => {
                        if anchor_a.is_some() {
                            return Err("Anchor A was already specified".to_string());
                        }
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--anchor-a-pos",
                            "candidates generate-between-anchors",
                        )?;
                        let parsed = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --anchor-a-pos value '{raw}': {e}"))?;
                        anchor_a = Some(SequenceAnchor::Position { zero_based: parsed });
                    }
                    "--anchor-b-pos" => {
                        if anchor_b.is_some() {
                            return Err("Anchor B was already specified".to_string());
                        }
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--anchor-b-pos",
                            "candidates generate-between-anchors",
                        )?;
                        let parsed = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --anchor-b-pos value '{raw}': {e}"))?;
                        anchor_b = Some(SequenceAnchor::Position { zero_based: parsed });
                    }
                    "--anchor-a-json" => {
                        if anchor_a.is_some() {
                            return Err("Anchor A was already specified".to_string());
                        }
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--anchor-a-json",
                            "candidates generate-between-anchors",
                        )?;
                        anchor_a = Some(parse_sequence_anchor_json(&raw, "--anchor-a-json")?);
                    }
                    "--anchor-b-json" => {
                        if anchor_b.is_some() {
                            return Err("Anchor B was already specified".to_string());
                        }
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--anchor-b-json",
                            "candidates generate-between-anchors",
                        )?;
                        anchor_b = Some(parse_sequence_anchor_json(&raw, "--anchor-b-json")?);
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for candidates generate-between-anchors"
                        ));
                    }
                }
            }
            let Some(length_bp) = length_bp else {
                return Err("candidates generate-between-anchors requires --length N".to_string());
            };
            let anchor_a = anchor_a.ok_or_else(|| {
                "candidates generate-between-anchors requires --anchor-a-pos N or --anchor-a-json JSON"
                    .to_string()
            })?;
            let anchor_b = anchor_b.ok_or_else(|| {
                "candidates generate-between-anchors requires --anchor-b-pos N or --anchor-b-json JSON"
                    .to_string()
            })?;
            Ok(ShellCommand::CandidatesGenerateBetweenAnchors {
                set_name,
                seq_id,
                anchor_a,
                anchor_b,
                length_bp,
                step_bp,
                limit,
            })
        }
        "show" => {
            if tokens.len() < 3 {
                return Err("candidates show requires SET_NAME [--limit N] [--offset N]".to_string());
            }
            let set_name = tokens[2].clone();
            let mut limit = DEFAULT_CANDIDATE_PAGE_SIZE;
            let mut offset = 0usize;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--limit" => {
                        let raw = parse_option_path(tokens, &mut idx, "--limit", "candidates show")?;
                        limit = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --limit value '{raw}': {e}"))?;
                        if limit == 0 {
                            return Err("--limit must be >= 1".to_string());
                        }
                    }
                    "--offset" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--offset", "candidates show")?;
                        offset = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --offset value '{raw}': {e}"))?;
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for candidates show"));
                    }
                }
            }
            Ok(ShellCommand::CandidatesShow {
                set_name,
                limit,
                offset,
            })
        }
        "metrics" => {
            if tokens.len() != 3 {
                return Err("candidates metrics requires SET_NAME".to_string());
            }
            Ok(ShellCommand::CandidatesMetrics {
                set_name: tokens[2].clone(),
            })
        }
        "score" => {
            if tokens.len() < 5 {
                return Err("candidates score requires SET_NAME METRIC_NAME EXPRESSION".to_string());
            }
            let set_name = tokens[2].clone();
            let metric = tokens[3].clone();
            let expression = tokens[4..].join(" ");
            if expression.trim().is_empty() {
                return Err("candidates score requires non-empty EXPRESSION".to_string());
            }
            Ok(ShellCommand::CandidatesScoreExpression {
                set_name,
                metric,
                expression,
            })
        }
        "score-distance" => {
            if tokens.len() < 4 {
                return Err(
                    "candidates score-distance requires SET_NAME METRIC_NAME [--feature-kind KIND] [--feature-label-regex REGEX] [--feature-geometry MODE] [--feature-boundary MODE] [--strand-relation MODE]"
                        .to_string(),
                );
            }
            let set_name = tokens[2].clone();
            let metric = tokens[3].clone();
            let mut feature_kinds: Vec<String> = vec![];
            let mut feature_label_regex: Option<String> = None;
            let mut feature_geometry_mode: Option<CandidateFeatureGeometryMode> = None;
            let mut feature_boundary_mode: Option<CandidateFeatureBoundaryMode> = None;
            let mut feature_strand_relation: Option<CandidateFeatureStrandRelation> = None;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--feature-kind" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--feature-kind",
                            "candidates score-distance",
                        )?;
                        let trimmed = raw.trim();
                        if !trimmed.is_empty() {
                            feature_kinds.push(trimmed.to_string());
                        }
                    }
                    "--feature-label-regex" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--feature-label-regex",
                            "candidates score-distance",
                        )?;
                        feature_label_regex = Some(raw);
                    }
                    "--feature-geometry" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--feature-geometry",
                            "candidates score-distance",
                        )?;
                        feature_geometry_mode = Some(parse_candidate_feature_geometry_mode(&raw)?);
                    }
                    "--feature-boundary" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--feature-boundary",
                            "candidates score-distance",
                        )?;
                        feature_boundary_mode = Some(parse_candidate_feature_boundary_mode(&raw)?);
                    }
                    "--strand-relation" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--strand-relation",
                            "candidates score-distance",
                        )?;
                        feature_strand_relation =
                            Some(parse_candidate_feature_strand_relation(&raw)?);
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for candidates score-distance"
                        ));
                    }
                }
            }
            Ok(ShellCommand::CandidatesScoreDistance {
                set_name,
                metric,
                feature_kinds,
                feature_label_regex,
                feature_geometry_mode,
                feature_boundary_mode,
                feature_strand_relation,
            })
        }
        "filter" => {
            if tokens.len() < 5 {
                return Err("candidates filter requires INPUT_SET OUTPUT_SET --metric METRIC_NAME [--min N] [--max N] [--min-quantile Q] [--max-quantile Q]".to_string());
            }
            let input_set = tokens[2].clone();
            let output_set = tokens[3].clone();
            let mut metric: Option<String> = None;
            let mut min: Option<f64> = None;
            let mut max: Option<f64> = None;
            let mut min_quantile: Option<f64> = None;
            let mut max_quantile: Option<f64> = None;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--metric" => {
                        metric = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--metric",
                            "candidates filter",
                        )?);
                    }
                    "--min" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--min", "candidates filter")?;
                        min = Some(
                            raw.parse::<f64>()
                                .map_err(|e| format!("Invalid --min value '{raw}': {e}"))?,
                        );
                    }
                    "--max" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--max", "candidates filter")?;
                        max = Some(
                            raw.parse::<f64>()
                                .map_err(|e| format!("Invalid --max value '{raw}': {e}"))?,
                        );
                    }
                    "--min-quantile" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--min-quantile",
                            "candidates filter",
                        )?;
                        min_quantile = Some(raw.parse::<f64>().map_err(|e| {
                            format!("Invalid --min-quantile value '{raw}': {e}")
                        })?);
                    }
                    "--max-quantile" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--max-quantile",
                            "candidates filter",
                        )?;
                        max_quantile = Some(raw.parse::<f64>().map_err(|e| {
                            format!("Invalid --max-quantile value '{raw}': {e}")
                        })?);
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for candidates filter"));
                    }
                }
            }
            let metric = metric.ok_or_else(|| "candidates filter requires --metric".to_string())?;
            if min
                .zip(max)
                .map(|(lo, hi)| lo > hi)
                .unwrap_or(false)
            {
                return Err("--min must be <= --max".to_string());
            }
            if min_quantile
                .zip(max_quantile)
                .map(|(lo, hi)| lo > hi)
                .unwrap_or(false)
            {
                return Err("--min-quantile must be <= --max-quantile".to_string());
            }
            for (name, value) in [("min-quantile", min_quantile), ("max-quantile", max_quantile)] {
                if let Some(q) = value {
                    if !(0.0..=1.0).contains(&q) {
                        return Err(format!("--{name} must be between 0 and 1"));
                    }
                }
            }
            if min.is_none() && max.is_none() && min_quantile.is_none() && max_quantile.is_none() {
                return Err(
                    "candidates filter requires at least one of --min/--max/--min-quantile/--max-quantile"
                        .to_string(),
                );
            }
            Ok(ShellCommand::CandidatesFilter {
                input_set,
                output_set,
                metric,
                min,
                max,
                min_quantile,
                max_quantile,
            })
        }
        "set-op" => {
            if tokens.len() != 6 {
                return Err(
                    "candidates set-op requires union|intersect|subtract LEFT_SET RIGHT_SET OUTPUT_SET"
                        .to_string(),
                );
            }
            let op = CandidateSetOperator::parse(&tokens[2]).ok_or_else(|| {
                format!(
                    "Unsupported candidates set-op '{}'; expected union|intersect|subtract",
                    tokens[2]
                )
            })?;
            Ok(ShellCommand::CandidatesSetOp {
                op,
                left_set: tokens[3].clone(),
                right_set: tokens[4].clone(),
                output_set: tokens[5].clone(),
            })
        }
        "macro" => {
            if tokens.len() < 3 {
                return Err(
                    "candidates macro requires SCRIPT_OR_@FILE (or --file PATH), optionally with --transactional".to_string(),
                );
            }
            let mut idx = 2usize;
            let mut transactional = false;
            let mut script_file: Option<String> = None;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--transactional" | "--atomic" => {
                        transactional = true;
                        idx += 1;
                    }
                    "--file" => {
                        if script_file.is_some() {
                            return Err(
                                "candidates macro --file may only be specified once".to_string(),
                            );
                        }
                        idx += 1;
                        if idx >= tokens.len() {
                            return Err("candidates macro --file requires PATH".to_string());
                        }
                        script_file = Some(tokens[idx].trim().to_string());
                        idx += 1;
                    }
                    _ => break,
                }
            }
            let script = if let Some(path) = script_file {
                if idx != tokens.len() {
                    return Err(
                        "candidates macro does not accept inline script after --file PATH"
                            .to_string(),
                    );
                }
                format!("@{path}")
            } else {
                if idx >= tokens.len() {
                    return Err("candidates macro requires SCRIPT_OR_@FILE".to_string());
                }
                tokens[idx..].join(" ")
            };
            if script.trim().is_empty() {
                return Err("candidates macro requires non-empty script".to_string());
            }
            Ok(ShellCommand::CandidatesMacro {
                script,
                transactional,
            })
        }
        other => Err(format!(
            "Unknown candidates subcommand '{other}' (expected list, delete, generate, generate-between-anchors, show, metrics, score, score-distance, filter, set-op, macro)"
        )),
    }
}

pub fn parse_shell_tokens(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.is_empty() {
        return Err("Missing shell command".to_string());
    }
    let cmd = tokens[0].as_str();
    match cmd {
        "help" | "-h" | "--help" => parse_help_command(tokens),
        "capabilities" => {
            if tokens.len() == 1 {
                Ok(ShellCommand::Capabilities)
            } else {
                Err(token_error(cmd))
            }
        }
        "state-summary" => {
            if tokens.len() == 1 {
                Ok(ShellCommand::StateSummary)
            } else {
                Err(token_error(cmd))
            }
        }
        "load-project" | "import-state" => {
            if tokens.len() == 2 {
                Ok(ShellCommand::LoadProject {
                    path: tokens[1].clone(),
                })
            } else {
                Err(token_error(cmd))
            }
        }
        "save-project" | "export-state" => {
            if tokens.len() == 2 {
                Ok(ShellCommand::SaveProject {
                    path: tokens[1].clone(),
                })
            } else {
                Err(token_error(cmd))
            }
        }
        "screenshot-window" => {
            let _ = tokens;
            Err(SCREENSHOT_DISABLED_MESSAGE.to_string())
        }
        "render-svg" => {
            if tokens.len() != 4 {
                return Err(token_error(cmd));
            }
            Ok(ShellCommand::RenderSvg {
                seq_id: tokens[1].clone(),
                mode: parse_mode(&tokens[2])?,
                output: tokens[3].clone(),
            })
        }
        "render-rna-svg" => {
            if tokens.len() != 3 {
                return Err(token_error(cmd));
            }
            Ok(ShellCommand::RenderRnaSvg {
                seq_id: tokens[1].clone(),
                output: tokens[2].clone(),
            })
        }
        "rna-info" => {
            if tokens.len() != 2 {
                return Err(token_error(cmd));
            }
            Ok(ShellCommand::RnaInfo {
                seq_id: tokens[1].clone(),
            })
        }
        "render-lineage-svg" => {
            if tokens.len() == 2 {
                Ok(ShellCommand::RenderLineageSvg {
                    output: tokens[1].clone(),
                })
            } else {
                Err(token_error(cmd))
            }
        }
        "render-pool-gel-svg" => {
            if tokens.len() < 3 {
                return Err(token_error(cmd));
            }
            let inputs = split_ids(&tokens[1]);
            if inputs.is_empty() {
                return Err("render-pool-gel-svg requires at least one sequence id".to_string());
            }
            let output = tokens[2].clone();
            let mut ladders: Option<Vec<String>> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--ladders" => {
                        if idx + 1 >= tokens.len() {
                            return Err("Missing value after --ladders".to_string());
                        }
                        let parsed = split_ids(&tokens[idx + 1]);
                        if !parsed.is_empty() {
                            ladders = Some(parsed);
                        }
                        idx += 2;
                    }
                    other => {
                        return Err(format!(
                            "Unknown argument '{other}' for render-pool-gel-svg"
                        ));
                    }
                }
            }
            Ok(ShellCommand::RenderPoolGelSvg {
                inputs,
                output,
                ladders,
            })
        }
        "ladders" => {
            if tokens.len() < 2 {
                return Err("ladders requires a subcommand: list or export".to_string());
            }
            match tokens[1].as_str() {
                "list" => {
                    let mut molecule = LadderMolecule::Dna;
                    let mut name_filter: Option<String> = None;
                    let mut idx = 2usize;
                    while idx < tokens.len() {
                        match tokens[idx].as_str() {
                            "--molecule" => {
                                if idx + 1 >= tokens.len() {
                                    return Err("Missing value after --molecule".to_string());
                                }
                                molecule = parse_ladder_molecule(&tokens[idx + 1])?;
                                idx += 2;
                            }
                            "--filter" => {
                                if idx + 1 >= tokens.len() {
                                    return Err("Missing value after --filter".to_string());
                                }
                                name_filter = Some(tokens[idx + 1].clone());
                                idx += 2;
                            }
                            other => {
                                return Err(format!("Unknown argument '{other}' for ladders list"));
                            }
                        }
                    }
                    Ok(ShellCommand::LaddersList {
                        molecule,
                        name_filter,
                    })
                }
                "export" => {
                    if tokens.len() < 3 {
                        return Err(
                            "ladders export requires: OUTPUT.json [--molecule dna|rna] [--filter TEXT]".to_string()
                        );
                    }
                    let output = tokens[2].clone();
                    let mut molecule = LadderMolecule::Dna;
                    let mut name_filter: Option<String> = None;
                    let mut idx = 3usize;
                    while idx < tokens.len() {
                        match tokens[idx].as_str() {
                            "--molecule" => {
                                if idx + 1 >= tokens.len() {
                                    return Err("Missing value after --molecule".to_string());
                                }
                                molecule = parse_ladder_molecule(&tokens[idx + 1])?;
                                idx += 2;
                            }
                            "--filter" => {
                                if idx + 1 >= tokens.len() {
                                    return Err("Missing value after --filter".to_string());
                                }
                                name_filter = Some(tokens[idx + 1].clone());
                                idx += 2;
                            }
                            other => {
                                return Err(format!(
                                    "Unknown argument '{other}' for ladders export"
                                ));
                            }
                        }
                    }
                    Ok(ShellCommand::LaddersExport {
                        molecule,
                        output,
                        name_filter,
                    })
                }
                other => Err(format!(
                    "Unknown ladders subcommand '{other}' (expected list or export)"
                )),
            }
        }
        "export-pool" => {
            if tokens.len() < 3 {
                return Err(token_error(cmd));
            }
            let inputs = split_ids(&tokens[1]);
            if inputs.is_empty() {
                return Err("export-pool requires at least one sequence id".to_string());
            }
            Ok(ShellCommand::ExportPool {
                inputs,
                output: tokens[2].clone(),
                human_id: tokens.get(3).cloned(),
            })
        }
        "import-pool" => {
            if tokens.len() > 3 {
                return Err(token_error(cmd));
            }
            if tokens.len() < 2 {
                return Err("import-pool requires: INPUT.pool.gentle.json [PREFIX]".to_string());
            }
            Ok(ShellCommand::ImportPool {
                input: tokens[1].clone(),
                prefix: tokens.get(2).cloned().unwrap_or_else(|| "pool".to_string()),
            })
        }
        "resources" => {
            if tokens.len() < 2 {
                return Err(
                    "resources requires a subcommand: sync-rebase or sync-jaspar".to_string(),
                );
            }
            match tokens[1].as_str() {
                "sync-rebase" => {
                    if tokens.len() < 3 {
                        return Err(
                            "resources sync-rebase requires INPUT.withrefm_or_URL".to_string()
                        );
                    }
                    let input = tokens[2].clone();
                    let mut output: Option<String> = None;
                    let mut commercial_only = false;
                    let mut idx = 3usize;
                    while idx < tokens.len() {
                        match tokens[idx].as_str() {
                            "--commercial-only" => {
                                commercial_only = true;
                                idx += 1;
                            }
                            value if value.starts_with("--") => {
                                return Err(format!(
                                    "Unknown option '{value}' for resources sync-rebase"
                                ));
                            }
                            value => {
                                if output.is_some() {
                                    return Err(format!(
                                        "Unexpected extra positional argument '{value}' for resources sync-rebase"
                                    ));
                                }
                                output = Some(value.to_string());
                                idx += 1;
                            }
                        }
                    }
                    Ok(ShellCommand::ResourcesSyncRebase {
                        input,
                        output,
                        commercial_only,
                    })
                }
                "sync-jaspar" => {
                    if tokens.len() < 3 {
                        return Err(
                            "resources sync-jaspar requires INPUT.jaspar_or_URL".to_string()
                        );
                    }
                    let input = tokens[2].clone();
                    let mut output: Option<String> = None;
                    let mut idx = 3usize;
                    while idx < tokens.len() {
                        let value = tokens[idx].as_str();
                        if value.starts_with("--") {
                            return Err(format!(
                                "Unknown option '{value}' for resources sync-jaspar"
                            ));
                        }
                        if output.is_some() {
                            return Err(format!(
                                "Unexpected extra positional argument '{value}' for resources sync-jaspar"
                            ));
                        }
                        output = Some(value.to_string());
                        idx += 1;
                    }
                    Ok(ShellCommand::ResourcesSyncJaspar { input, output })
                }
                other => Err(format!(
                    "Unknown resources subcommand '{other}' (expected sync-rebase or sync-jaspar)"
                )),
            }
        }
        "tracks" => {
            if tokens.len() < 2 {
                return Err(
                    "tracks requires a subcommand: import-bed, import-bigwig, import-vcf, or tracked"
                        .to_string(),
                );
            }
            match tokens[1].as_str() {
                "import-bed" => {
                    if tokens.len() < 4 {
                        return Err(
                            "tracks import-bed requires SEQ_ID PATH [--name NAME] [--min-score N] [--max-score N] [--clear-existing]".to_string()
                        );
                    }
                    let seq_id = tokens[2].clone();
                    let path = tokens[3].clone();
                    let mut track_name: Option<String> = None;
                    let mut min_score: Option<f64> = None;
                    let mut max_score: Option<f64> = None;
                    let mut clear_existing = false;
                    let mut idx = 4usize;
                    while idx < tokens.len() {
                        match tokens[idx].as_str() {
                            "--name" => {
                                track_name = Some(parse_option_path(
                                    tokens,
                                    &mut idx,
                                    "--name",
                                    "tracks import-bed",
                                )?);
                            }
                            "--min-score" => {
                                let raw = parse_option_path(
                                    tokens,
                                    &mut idx,
                                    "--min-score",
                                    "tracks import-bed",
                                )?;
                                min_score = Some(raw.parse::<f64>().map_err(|e| {
                                    format!("Invalid --min-score value '{raw}': {e}")
                                })?);
                            }
                            "--max-score" => {
                                let raw = parse_option_path(
                                    tokens,
                                    &mut idx,
                                    "--max-score",
                                    "tracks import-bed",
                                )?;
                                max_score = Some(raw.parse::<f64>().map_err(|e| {
                                    format!("Invalid --max-score value '{raw}': {e}")
                                })?);
                            }
                            "--clear-existing" => {
                                clear_existing = true;
                                idx += 1;
                            }
                            other => {
                                return Err(format!(
                                    "Unknown option '{other}' for tracks import-bed"
                                ));
                            }
                        }
                    }
                    if min_score
                        .zip(max_score)
                        .map(|(min, max)| min > max)
                        .unwrap_or(false)
                    {
                        return Err("--min-score must be <= --max-score".to_string());
                    }
                    Ok(ShellCommand::TracksImportBed {
                        seq_id,
                        path,
                        track_name,
                        min_score,
                        max_score,
                        clear_existing,
                    })
                }
                "import-bigwig" => {
                    if tokens.len() < 4 {
                        return Err(
                            "tracks import-bigwig requires SEQ_ID PATH [--name NAME] [--min-score N] [--max-score N] [--clear-existing]".to_string()
                        );
                    }
                    let seq_id = tokens[2].clone();
                    let path = tokens[3].clone();
                    let mut track_name: Option<String> = None;
                    let mut min_score: Option<f64> = None;
                    let mut max_score: Option<f64> = None;
                    let mut clear_existing = false;
                    let mut idx = 4usize;
                    while idx < tokens.len() {
                        match tokens[idx].as_str() {
                            "--name" => {
                                track_name = Some(parse_option_path(
                                    tokens,
                                    &mut idx,
                                    "--name",
                                    "tracks import-bigwig",
                                )?);
                            }
                            "--min-score" => {
                                let raw = parse_option_path(
                                    tokens,
                                    &mut idx,
                                    "--min-score",
                                    "tracks import-bigwig",
                                )?;
                                min_score = Some(raw.parse::<f64>().map_err(|e| {
                                    format!("Invalid --min-score value '{raw}': {e}")
                                })?);
                            }
                            "--max-score" => {
                                let raw = parse_option_path(
                                    tokens,
                                    &mut idx,
                                    "--max-score",
                                    "tracks import-bigwig",
                                )?;
                                max_score = Some(raw.parse::<f64>().map_err(|e| {
                                    format!("Invalid --max-score value '{raw}': {e}")
                                })?);
                            }
                            "--clear-existing" => {
                                clear_existing = true;
                                idx += 1;
                            }
                            other => {
                                return Err(format!(
                                    "Unknown option '{other}' for tracks import-bigwig"
                                ));
                            }
                        }
                    }
                    if min_score
                        .zip(max_score)
                        .map(|(min, max)| min > max)
                        .unwrap_or(false)
                    {
                        return Err("--min-score must be <= --max-score".to_string());
                    }
                    Ok(ShellCommand::TracksImportBigWig {
                        seq_id,
                        path,
                        track_name,
                        min_score,
                        max_score,
                        clear_existing,
                    })
                }
                "import-vcf" => {
                    if tokens.len() < 4 {
                        return Err(
                            "tracks import-vcf requires SEQ_ID PATH [--name NAME] [--min-score N] [--max-score N] [--clear-existing]".to_string()
                        );
                    }
                    let seq_id = tokens[2].clone();
                    let path = tokens[3].clone();
                    let mut track_name: Option<String> = None;
                    let mut min_score: Option<f64> = None;
                    let mut max_score: Option<f64> = None;
                    let mut clear_existing = false;
                    let mut idx = 4usize;
                    while idx < tokens.len() {
                        match tokens[idx].as_str() {
                            "--name" => {
                                track_name = Some(parse_option_path(
                                    tokens,
                                    &mut idx,
                                    "--name",
                                    "tracks import-vcf",
                                )?);
                            }
                            "--min-score" => {
                                let raw = parse_option_path(
                                    tokens,
                                    &mut idx,
                                    "--min-score",
                                    "tracks import-vcf",
                                )?;
                                min_score = Some(raw.parse::<f64>().map_err(|e| {
                                    format!("Invalid --min-score value '{raw}': {e}")
                                })?);
                            }
                            "--max-score" => {
                                let raw = parse_option_path(
                                    tokens,
                                    &mut idx,
                                    "--max-score",
                                    "tracks import-vcf",
                                )?;
                                max_score = Some(raw.parse::<f64>().map_err(|e| {
                                    format!("Invalid --max-score value '{raw}': {e}")
                                })?);
                            }
                            "--clear-existing" => {
                                clear_existing = true;
                                idx += 1;
                            }
                            other => {
                                return Err(format!(
                                    "Unknown option '{other}' for tracks import-vcf"
                                ));
                            }
                        }
                    }
                    if min_score
                        .zip(max_score)
                        .map(|(min, max)| min > max)
                        .unwrap_or(false)
                    {
                        return Err("--min-score must be <= --max-score".to_string());
                    }
                    Ok(ShellCommand::TracksImportVcf {
                        seq_id,
                        path,
                        track_name,
                        min_score,
                        max_score,
                        clear_existing,
                    })
                }
                "tracked" => {
                    if tokens.len() < 3 {
                        return Err(
                            "tracks tracked requires a subcommand: list, add, remove, clear, apply, sync".to_string(),
                        );
                    }
                    match tokens[2].as_str() {
                        "list" => {
                            if tokens.len() > 3 {
                                return Err("tracks tracked list takes no options".to_string());
                            }
                            Ok(ShellCommand::TracksTrackedList)
                        }
                        "add" => {
                            if tokens.len() < 4 {
                                return Err(
                                    "tracks tracked add requires PATH [--source auto|bed|bigwig|vcf] [--name NAME] [--min-score N] [--max-score N] [--clear-existing]".to_string(),
                                );
                            }
                            let path = tokens[3].clone();
                            let mut source: Option<GenomeTrackSource> = None;
                            let mut track_name: Option<String> = None;
                            let mut min_score: Option<f64> = None;
                            let mut max_score: Option<f64> = None;
                            let mut clear_existing = false;
                            let mut idx = 4usize;
                            while idx < tokens.len() {
                                match tokens[idx].as_str() {
                                    "--source" => {
                                        let raw = parse_option_path(
                                            tokens,
                                            &mut idx,
                                            "--source",
                                            "tracks tracked add",
                                        )?;
                                        let normalized = raw.trim().to_ascii_lowercase();
                                        source = Some(match normalized.as_str() {
                                            "auto" => GenomeTrackSource::from_path(&path),
                                            "bed" => GenomeTrackSource::Bed,
                                            "bigwig" | "bw" => GenomeTrackSource::BigWig,
                                            "vcf" => GenomeTrackSource::Vcf,
                                            _ => {
                                                return Err(format!(
                                                    "Unsupported --source value '{}'; expected auto|bed|bigwig|vcf",
                                                    raw
                                                ))
                                            }
                                        });
                                    }
                                    "--name" => {
                                        track_name = Some(parse_option_path(
                                            tokens,
                                            &mut idx,
                                            "--name",
                                            "tracks tracked add",
                                        )?);
                                    }
                                    "--min-score" => {
                                        let raw = parse_option_path(
                                            tokens,
                                            &mut idx,
                                            "--min-score",
                                            "tracks tracked add",
                                        )?;
                                        min_score = Some(raw.parse::<f64>().map_err(|e| {
                                            format!("Invalid --min-score value '{raw}': {e}")
                                        })?);
                                    }
                                    "--max-score" => {
                                        let raw = parse_option_path(
                                            tokens,
                                            &mut idx,
                                            "--max-score",
                                            "tracks tracked add",
                                        )?;
                                        max_score = Some(raw.parse::<f64>().map_err(|e| {
                                            format!("Invalid --max-score value '{raw}': {e}")
                                        })?);
                                    }
                                    "--clear-existing" => {
                                        clear_existing = true;
                                        idx += 1;
                                    }
                                    other => {
                                        return Err(format!(
                                            "Unknown option '{other}' for tracks tracked add"
                                        ));
                                    }
                                }
                            }
                            if min_score
                                .zip(max_score)
                                .map(|(min, max)| min > max)
                                .unwrap_or(false)
                            {
                                return Err("--min-score must be <= --max-score".to_string());
                            }
                            let subscription = GenomeTrackSubscription {
                                source: source.unwrap_or_else(|| GenomeTrackSource::from_path(&path)),
                                path,
                                track_name,
                                min_score,
                                max_score,
                                clear_existing,
                            };
                            Ok(ShellCommand::TracksTrackedAdd { subscription })
                        }
                        "remove" => {
                            if tokens.len() != 4 {
                                return Err(
                                    "tracks tracked remove requires INDEX".to_string(),
                                );
                            }
                            let index = tokens[3]
                                .parse::<usize>()
                                .map_err(|e| format!("Invalid INDEX '{}': {e}", tokens[3]))?;
                            Ok(ShellCommand::TracksTrackedRemove { index })
                        }
                        "clear" => {
                            if tokens.len() > 3 {
                                return Err("tracks tracked clear takes no options".to_string());
                            }
                            Ok(ShellCommand::TracksTrackedClear)
                        }
                        "apply" => {
                            let mut index: Option<usize> = None;
                            let mut only_new_anchors = false;
                            let mut idx = 3usize;
                            while idx < tokens.len() {
                                match tokens[idx].as_str() {
                                    "--index" => {
                                        let raw = parse_option_path(
                                            tokens,
                                            &mut idx,
                                            "--index",
                                            "tracks tracked apply",
                                        )?;
                                        index = Some(raw.parse::<usize>().map_err(|e| {
                                            format!("Invalid --index value '{raw}': {e}")
                                        })?);
                                    }
                                    "--only-new-anchors" => {
                                        only_new_anchors = true;
                                        idx += 1;
                                    }
                                    other => {
                                        return Err(format!(
                                            "Unknown option '{other}' for tracks tracked apply"
                                        ));
                                    }
                                }
                            }
                            Ok(ShellCommand::TracksTrackedApply {
                                index,
                                only_new_anchors,
                            })
                        }
                        "sync" => {
                            if tokens.len() > 3 {
                                return Err("tracks tracked sync takes no options".to_string());
                            }
                            Ok(ShellCommand::TracksTrackedApply {
                                index: None,
                                only_new_anchors: true,
                            })
                        }
                        other => Err(format!(
                            "Unknown tracks tracked subcommand '{other}' (expected list, add, remove, clear, apply, sync)"
                        )),
                    }
                }
                other => Err(format!(
                    "Unknown tracks subcommand '{other}' (expected import-bed, import-bigwig, import-vcf, or tracked)"
                )),
            }
        }
        "genomes" => parse_reference_command(tokens, false),
        "helpers" => parse_reference_command(tokens, true),
        "candidates" => parse_candidates_command(tokens),
        "set-param" => {
            if tokens.len() < 3 {
                return Err("set-param requires NAME JSON_VALUE".to_string());
            }
            let name = tokens[1].trim().to_string();
            if name.is_empty() {
                return Err("set-param NAME must not be empty".to_string());
            }
            let value_json = tokens[2..].join(" ");
            if value_json.trim().is_empty() {
                return Err("set-param JSON_VALUE must not be empty".to_string());
            }
            Ok(ShellCommand::SetParameter { name, value_json })
        }
        "op" => {
            let payload = tokens[1..].join(" ");
            if payload.trim().is_empty() {
                return Err("Missing operation JSON".to_string());
            }
            Ok(ShellCommand::Op { payload })
        }
        "workflow" => {
            let payload = tokens[1..].join(" ");
            if payload.trim().is_empty() {
                return Err("Missing workflow JSON".to_string());
            }
            Ok(ShellCommand::Workflow { payload })
        }
        other => Err(format!("Unknown shell command '{other}'. Try: help")),
    }
}

pub fn parse_shell_line(line: &str) -> Result<ShellCommand, String> {
    let tokens = split_shell_words(line)?;
    parse_shell_tokens(&tokens)
}

pub fn split_shell_words(line: &str) -> Result<Vec<String>, String> {
    #[derive(Clone, Copy, PartialEq, Eq)]
    enum Mode {
        Normal,
        SingleQuoted,
        DoubleQuoted,
    }

    let mut out = Vec::new();
    let mut current = String::new();
    let mut mode = Mode::Normal;
    let mut chars = line.chars().peekable();

    while let Some(ch) = chars.next() {
        match mode {
            Mode::Normal => match ch {
                '\'' => mode = Mode::SingleQuoted,
                '"' => mode = Mode::DoubleQuoted,
                '\\' => {
                    if let Some(next) = chars.next() {
                        current.push(next);
                    }
                }
                c if c.is_whitespace() => {
                    if !current.is_empty() {
                        out.push(current.clone());
                        current.clear();
                    }
                }
                _ => current.push(ch),
            },
            Mode::SingleQuoted => {
                if ch == '\'' {
                    mode = Mode::Normal;
                } else {
                    current.push(ch);
                }
            }
            Mode::DoubleQuoted => {
                if ch == '"' {
                    mode = Mode::Normal;
                } else if ch == '\\' {
                    if let Some(next) = chars.next() {
                        current.push(next);
                    }
                } else {
                    current.push(ch);
                }
            }
        }
    }

    if mode != Mode::Normal {
        return Err("Unterminated quoted string in shell command".to_string());
    }
    if !current.is_empty() {
        out.push(current);
    }
    if out.is_empty() {
        return Err("Empty shell command".to_string());
    }
    Ok(out)
}

fn load_candidates_macro_script(script_or_ref: &str) -> Result<String, String> {
    let trimmed = script_or_ref.trim();
    if let Some(path) = trimmed.strip_prefix('@') {
        let path = path.trim();
        if path.is_empty() {
            return Err("candidates macro @FILE requires a non-empty file path".to_string());
        }
        fs::read_to_string(path)
            .map_err(|e| format!("Could not read candidates macro file '{path}': {e}"))
    } else {
        Ok(trimmed.to_string())
    }
}

fn split_candidates_macro_statements(script: &str) -> Result<Vec<String>, String> {
    #[derive(Clone, Copy, PartialEq, Eq)]
    enum Mode {
        Normal,
        SingleQuoted,
        DoubleQuoted,
    }
    let mut out: Vec<String> = vec![];
    let mut current = String::new();
    let mut mode = Mode::Normal;
    let mut chars = script.chars().peekable();
    while let Some(ch) = chars.next() {
        match mode {
            Mode::Normal => match ch {
                '\'' => {
                    mode = Mode::SingleQuoted;
                    current.push(ch);
                }
                '"' => {
                    mode = Mode::DoubleQuoted;
                    current.push(ch);
                }
                ';' | '\n' | '\r' => {
                    let stmt = current.trim();
                    if !stmt.is_empty() {
                        out.push(stmt.to_string());
                    }
                    current.clear();
                }
                '#' if current.trim().is_empty() => {
                    for next in chars.by_ref() {
                        if next == '\n' {
                            break;
                        }
                    }
                    let stmt = current.trim();
                    if !stmt.is_empty() {
                        out.push(stmt.to_string());
                    }
                    current.clear();
                }
                _ => current.push(ch),
            },
            Mode::SingleQuoted => {
                current.push(ch);
                if ch == '\'' {
                    mode = Mode::Normal;
                }
            }
            Mode::DoubleQuoted => {
                current.push(ch);
                if ch == '\\' {
                    if let Some(next) = chars.next() {
                        current.push(next);
                    }
                } else if ch == '"' {
                    mode = Mode::Normal;
                }
            }
        }
    }
    if mode != Mode::Normal {
        return Err("Unterminated quoted string in candidates macro script".to_string());
    }
    let tail = current.trim();
    if !tail.is_empty() {
        out.push(tail.to_string());
    }
    Ok(out)
}

fn run_candidates_macro(
    engine: &mut GentleEngine,
    script_or_ref: &str,
    transactional: bool,
    options: &ShellExecutionOptions,
) -> Result<ShellRunResult, String> {
    let script = load_candidates_macro_script(script_or_ref)?;
    let statements = split_candidates_macro_statements(&script)?;
    if statements.is_empty() {
        return Err("candidates macro script is empty".to_string());
    }
    let rollback_state = if transactional {
        Some(engine.state().clone())
    } else {
        None
    };
    let mut executed = 0usize;
    let mut changed = false;
    let mut rows: Vec<Value> = vec![];
    for statement in statements {
        let statement = statement.trim();
        if statement.is_empty() {
            continue;
        }
        let prefixed = if statement.starts_with("candidates ") {
            statement.to_string()
        } else {
            format!("candidates {statement}")
        };
        let tokens = split_shell_words(&prefixed)?;
        let cmd = parse_candidates_command(&tokens)?;
        if matches!(cmd, ShellCommand::CandidatesMacro { .. }) {
            return Err("Nested candidates macro calls are not allowed".to_string());
        }
        let run = match execute_shell_command_with_options(engine, &cmd, options) {
            Ok(run) => run,
            Err(err) => {
                if transactional {
                    if let Some(state) = rollback_state {
                        *engine = GentleEngine::from_state(state);
                    }
                    return Err(format!(
                        "candidates macro failed at statement {} ('{}'): {err}; all macro changes were rolled back",
                        executed + 1,
                        statement
                    ));
                }
                return Err(format!(
                    "candidates macro failed at statement {} ('{}'): {err}",
                    executed + 1,
                    statement
                ));
            }
        };
        executed = executed.saturating_add(1);
        changed |= run.state_changed;
        rows.push(json!({
            "statement": statement,
            "state_changed": run.state_changed,
            "output": run.output
        }));
    }
    if executed == 0 {
        return Err("candidates macro script has no executable statements".to_string());
    }
    Ok(ShellRunResult {
        state_changed: changed,
        output: json!({
            "executed": executed,
            "transactional": transactional,
            "state_changed": changed,
            "results": rows
        }),
    })
}

pub fn execute_shell_command(
    engine: &mut GentleEngine,
    command: &ShellCommand,
) -> Result<ShellRunResult, String> {
    execute_shell_command_with_options(engine, command, &ShellExecutionOptions::default())
}

pub fn execute_shell_command_with_options(
    engine: &mut GentleEngine,
    command: &ShellCommand,
    options: &ShellExecutionOptions,
) -> Result<ShellRunResult, String> {
    let result = match command {
        ShellCommand::Help {
            topic,
            format,
            interface_filter,
        } => {
            let help_output = if topic.is_empty() {
                match format {
                    HelpOutputFormat::Text => {
                        json!({ "help": render_shell_help_text(interface_filter.as_deref())? })
                    }
                    HelpOutputFormat::Json => render_shell_help_json(interface_filter.as_deref())?,
                    HelpOutputFormat::Markdown => json!({
                        "help_markdown": render_shell_help_markdown(interface_filter.as_deref())?
                    }),
                }
            } else {
                match format {
                    HelpOutputFormat::Text => json!({
                        "topic": topic.join(" "),
                        "help": render_shell_topic_help_text(topic, interface_filter.as_deref())?
                    }),
                    HelpOutputFormat::Json => {
                        render_shell_topic_help_json(topic, interface_filter.as_deref())?
                    }
                    HelpOutputFormat::Markdown => json!({
                        "topic": topic.join(" "),
                        "help_markdown": render_shell_topic_help_markdown(topic, interface_filter.as_deref())?
                    }),
                }
            };
            ShellRunResult {
                state_changed: false,
                output: help_output,
            }
        }
        ShellCommand::Capabilities => ShellRunResult {
            state_changed: false,
            output: serde_json::to_value(GentleEngine::capabilities())
                .map_err(|e| format!("Could not serialize capabilities: {e}"))?,
        },
        ShellCommand::StateSummary => ShellRunResult {
            state_changed: false,
            output: serde_json::to_value(engine.summarize_state())
                .map_err(|e| format!("Could not serialize state summary: {e}"))?,
        },
        ShellCommand::LoadProject { path } => {
            let state = ProjectState::load_from_path(path).map_err(|e| e.to_string())?;
            *engine = GentleEngine::from_state(state);
            ShellRunResult {
                state_changed: true,
                output: json!({
                    "message": format!("Loaded project from '{path}'"),
                    "summary": engine.summarize_state()
                }),
            }
        }
        ShellCommand::SaveProject { path } => {
            engine
                .state()
                .save_to_path(path)
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({ "message": format!("Saved project to '{path}'") }),
            }
        }
        ShellCommand::ScreenshotWindow { output: _ } => {
            let _ = options;
            return Err(SCREENSHOT_DISABLED_MESSAGE.to_string());
        }
        ShellCommand::RenderSvg {
            seq_id,
            mode,
            output,
        } => {
            let op_result = engine
                .apply(Operation::RenderSequenceSvg {
                    seq_id: seq_id.clone(),
                    mode: mode.clone(),
                    path: output.clone(),
                })
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::RenderRnaSvg { seq_id, output } => {
            let op_result = engine
                .apply(Operation::RenderRnaStructureSvg {
                    seq_id: seq_id.clone(),
                    path: output.clone(),
                })
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::RnaInfo { seq_id } => {
            let report = engine
                .inspect_rna_structure(seq_id)
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: serde_json::to_value(report)
                    .map_err(|e| format!("Could not serialize RNA report: {e}"))?,
            }
        }
        ShellCommand::RenderLineageSvg { output } => {
            let op_result = engine
                .apply(Operation::RenderLineageSvg {
                    path: output.clone(),
                })
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::RenderPoolGelSvg {
            inputs,
            output,
            ladders,
        } => {
            let op_result = engine
                .apply(Operation::RenderPoolGelSvg {
                    inputs: inputs.clone(),
                    path: output.clone(),
                    ladders: ladders.clone(),
                })
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::LaddersList {
            molecule,
            name_filter,
        } => ShellRunResult {
            state_changed: false,
            output: match molecule {
                LadderMolecule::Dna => {
                    serde_json::to_value(GentleEngine::inspect_dna_ladders(name_filter.as_deref()))
                        .map_err(|e| format!("Could not serialize DNA ladders catalog: {e}"))?
                }
                LadderMolecule::Rna => {
                    serde_json::to_value(GentleEngine::inspect_rna_ladders(name_filter.as_deref()))
                        .map_err(|e| format!("Could not serialize RNA ladders catalog: {e}"))?
                }
            },
        },
        ShellCommand::LaddersExport {
            molecule,
            output,
            name_filter,
        } => {
            let op_result = match molecule {
                LadderMolecule::Dna => engine
                    .apply(Operation::ExportDnaLadders {
                        path: output.clone(),
                        name_filter: name_filter.clone(),
                    })
                    .map_err(|e| e.to_string())?,
                LadderMolecule::Rna => engine
                    .apply(Operation::ExportRnaLadders {
                        path: output.clone(),
                        name_filter: name_filter.clone(),
                    })
                    .map_err(|e| e.to_string())?,
            };
            ShellRunResult {
                state_changed: false,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::ExportPool {
            inputs,
            output,
            human_id,
        } => {
            let op_result = engine
                .apply(Operation::ExportPool {
                    inputs: inputs.clone(),
                    path: output.clone(),
                    pool_id: Some("pool_export".to_string()),
                    human_id: human_id.clone(),
                })
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::ImportPool { input, prefix } => {
            let text = fs::read_to_string(input)
                .map_err(|e| format!("Could not read pool file '{input}': {e}"))?;
            let pool: PoolExport = serde_json::from_str(&text)
                .map_err(|e| format!("Invalid pool JSON '{input}': {e}"))?;
            if pool.schema != "gentle.pool.v1" {
                return Err(format!(
                    "Unsupported pool schema '{}', expected 'gentle.pool.v1'",
                    pool.schema
                ));
            }

            let mut state = engine.state().clone();
            let mut imported_ids = Vec::new();
            for (idx, member) in pool.members.iter().enumerate() {
                let mut dna = DNAsequence::from_sequence(&member.sequence)
                    .map_err(|e| format!("Invalid DNA in pool member '{}': {e}", member.seq_id))?;
                if let Some(name) = &member.name {
                    let mut value = serde_json::to_value(&dna)
                        .map_err(|e| format!("Could not serialize sequence: {e}"))?;
                    if let Some(obj) = value.as_object_mut() {
                        if let Some(seq_obj) = obj.get_mut("seq").and_then(|v| v.as_object_mut()) {
                            seq_obj.insert("name".to_string(), json!(name));
                        }
                    }
                    dna = serde_json::from_value(value)
                        .map_err(|e| format!("Could not set sequence name: {e}"))?;
                }
                if member.topology.eq_ignore_ascii_case("circular") {
                    dna.set_circular(true);
                }
                apply_member_overhang(member, &mut dna)?;
                dna.update_computed_features();
                let base = format!("{prefix}_{}", idx + 1);
                let id = unique_id(&state.sequences, &base);
                state.sequences.insert(id.clone(), dna);
                imported_ids.push(id);
            }
            *engine = GentleEngine::from_state(state);
            ShellRunResult {
                state_changed: true,
                output: json!({
                    "message": format!("Imported pool '{}' ({} members)", pool.pool_id, pool.member_count),
                    "input": input,
                    "pool_id": pool.pool_id,
                    "member_count": pool.member_count,
                    "imported_ids": imported_ids,
                }),
            }
        }
        ShellCommand::ResourcesSyncRebase {
            input,
            output,
            commercial_only,
        } => {
            let report = resource_sync::sync_rebase(input, output.as_deref(), *commercial_only)?;
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "message": format!("Synced {} {} entries to '{}'", report.item_count, report.resource, report.output),
                    "report": report,
                }),
            }
        }
        ShellCommand::ResourcesSyncJaspar { input, output } => {
            let report = resource_sync::sync_jaspar(input, output.as_deref())?;
            tf_motifs::reload();
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "message": format!("Synced {} {} entries to '{}'", report.item_count, report.resource, report.output),
                    "report": report,
                }),
            }
        }
        ShellCommand::ReferenceList {
            helper_mode,
            catalog_path,
        } => {
            let resolved_catalog = resolved_catalog_path(catalog_path, *helper_mode);
            let genomes = GentleEngine::list_reference_genomes(resolved_catalog)
                .map_err(|e| e.to_string())?;
            let effective_catalog = effective_catalog_path(catalog_path, *helper_mode);
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "catalog_path": effective_catalog,
                    "genome_count": genomes.len(),
                    "genomes": genomes,
                }),
            }
        }
        ShellCommand::ReferenceValidateCatalog {
            helper_mode,
            catalog_path,
        } => {
            let resolved_catalog = resolved_catalog_path(catalog_path, *helper_mode);
            let effective_catalog = effective_catalog_path(catalog_path, *helper_mode);
            let genomes = GentleEngine::list_reference_genomes(resolved_catalog)
                .map_err(|e| e.to_string())?;
            for genome_id in &genomes {
                GentleEngine::describe_reference_genome_sources(resolved_catalog, genome_id, None)
                    .map_err(|e| e.to_string())?;
            }
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "catalog_path": effective_catalog,
                    "valid": true,
                    "genome_count": genomes.len(),
                    "validated_sources": genomes.len(),
                    "genomes": genomes,
                }),
            }
        }
        ShellCommand::ReferenceStatus {
            helper_mode,
            genome_id,
            catalog_path,
            cache_dir,
        } => {
            let resolved_catalog = resolved_catalog_path(catalog_path, *helper_mode);
            let prepared = GentleEngine::is_reference_genome_prepared(
                resolved_catalog,
                genome_id,
                cache_dir.as_deref(),
            )
            .map_err(|e| e.to_string())?;
            let source_plan = GentleEngine::describe_reference_genome_sources(
                resolved_catalog,
                genome_id,
                cache_dir.as_deref(),
            )
            .map_err(|e| e.to_string())?;
            let effective_catalog = effective_catalog_path(catalog_path, *helper_mode);
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "genome_id": genome_id,
                    "catalog_path": effective_catalog,
                    "cache_dir": cache_dir,
                    "prepared": prepared,
                    "sequence_source_type": source_plan.sequence_source_type,
                    "annotation_source_type": source_plan.annotation_source_type,
                    "sequence_source": source_plan.sequence_source,
                    "annotation_source": source_plan.annotation_source,
                }),
            }
        }
        ShellCommand::ReferenceGenes {
            helper_mode,
            genome_id,
            catalog_path,
            cache_dir,
            filter,
            biotypes,
            limit,
            offset,
        } => {
            let resolved_catalog = resolved_catalog_path(catalog_path, *helper_mode);
            let genes = GentleEngine::list_reference_genome_genes(
                resolved_catalog,
                genome_id,
                cache_dir.as_deref(),
            )
            .map_err(|e| e.to_string())?;
            let filter_regex = compile_gene_filter_regex(filter)?;
            let biotype_filter: Vec<String> = biotypes
                .iter()
                .map(|v| v.trim().to_ascii_lowercase())
                .filter(|v| !v.is_empty())
                .collect();
            let available_biotypes = collect_biotypes(&genes);
            let filtered: Vec<GenomeGeneRecord> = genes
                .into_iter()
                .filter(|g| genome_gene_matches_filter(g, filter_regex.as_ref(), &biotype_filter))
                .collect();
            let total = filtered.len();
            let clamped_offset = (*offset).min(total);
            let returned: Vec<GenomeGeneRecord> = filtered
                .into_iter()
                .skip(clamped_offset)
                .take(*limit)
                .collect();
            let effective_catalog = effective_catalog_path(catalog_path, *helper_mode);
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "genome_id": genome_id,
                    "catalog_path": effective_catalog,
                    "cache_dir": cache_dir,
                    "filter": filter,
                    "biotype_filter": biotypes,
                    "available_biotypes": available_biotypes,
                    "offset": clamped_offset,
                    "limit": limit,
                    "total_matches": total,
                    "returned": returned.len(),
                    "genes": returned,
                }),
            }
        }
        ShellCommand::ReferencePrepare {
            helper_mode,
            genome_id,
            catalog_path,
            cache_dir,
            timeout_seconds,
        } => {
            let op = Operation::PrepareGenome {
                genome_id: genome_id.clone(),
                catalog_path: operation_catalog_path(catalog_path, *helper_mode),
                cache_dir: cache_dir.clone(),
                timeout_seconds: *timeout_seconds,
            };
            let op_result = engine.apply(op).map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: true,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::ReferenceBlast {
            helper_mode,
            genome_id,
            query_sequence,
            max_hits,
            task,
            catalog_path,
            cache_dir,
        } => {
            let resolved_catalog = resolved_catalog_path(catalog_path, *helper_mode);
            let report = if *helper_mode {
                GentleEngine::blast_helper_genome(
                    genome_id,
                    query_sequence,
                    *max_hits,
                    task.as_deref(),
                    resolved_catalog,
                    cache_dir.as_deref(),
                )
            } else {
                GentleEngine::blast_reference_genome(
                    resolved_catalog,
                    genome_id,
                    query_sequence,
                    *max_hits,
                    task.as_deref(),
                    cache_dir.as_deref(),
                )
            }
            .map_err(|e| e.to_string())?;
            let effective_catalog = effective_catalog_path(catalog_path, *helper_mode);
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "catalog_path": effective_catalog,
                    "cache_dir": cache_dir,
                    "report": report,
                }),
            }
        }
        ShellCommand::ReferenceBlastTrack {
            helper_mode,
            genome_id,
            query_sequence,
            target_seq_id,
            max_hits,
            task,
            track_name,
            clear_existing,
            catalog_path,
            cache_dir,
        } => {
            let resolved_catalog = resolved_catalog_path(catalog_path, *helper_mode);
            let report = if *helper_mode {
                GentleEngine::blast_helper_genome(
                    genome_id,
                    query_sequence,
                    *max_hits,
                    task.as_deref(),
                    resolved_catalog,
                    cache_dir.as_deref(),
                )
            } else {
                GentleEngine::blast_reference_genome(
                    resolved_catalog,
                    genome_id,
                    query_sequence,
                    *max_hits,
                    task.as_deref(),
                    cache_dir.as_deref(),
                )
            }
            .map_err(|e| e.to_string())?;
            let hit_inputs = report
                .hits
                .iter()
                .map(|hit| crate::engine::BlastHitFeatureInput {
                    subject_id: hit.subject_id.clone(),
                    query_start_1based: hit.query_start,
                    query_end_1based: hit.query_end,
                    subject_start_1based: hit.subject_start,
                    subject_end_1based: hit.subject_end,
                    identity_percent: hit.identity_percent,
                    bit_score: hit.bit_score,
                    evalue: hit.evalue,
                    query_coverage_percent: hit.query_coverage_percent,
                })
                .collect::<Vec<_>>();
            let op_result = engine
                .apply(Operation::ImportBlastHitsTrack {
                    seq_id: target_seq_id.clone(),
                    hits: hit_inputs,
                    track_name: track_name.clone(),
                    clear_existing: Some(*clear_existing),
                })
                .map_err(|e| e.to_string())?;
            let state_changed =
                !op_result.created_seq_ids.is_empty() || !op_result.changed_seq_ids.is_empty();
            let effective_catalog = effective_catalog_path(catalog_path, *helper_mode);
            ShellRunResult {
                state_changed,
                output: json!({
                    "catalog_path": effective_catalog,
                    "cache_dir": cache_dir,
                    "report": report,
                    "result": op_result
                }),
            }
        }
        ShellCommand::ReferenceExtractRegion {
            helper_mode,
            genome_id,
            chromosome,
            start_1based,
            end_1based,
            output_id,
            catalog_path,
            cache_dir,
        } => {
            let op_result = engine
                .apply(Operation::ExtractGenomeRegion {
                    genome_id: genome_id.clone(),
                    chromosome: chromosome.clone(),
                    start_1based: *start_1based,
                    end_1based: *end_1based,
                    output_id: output_id.clone(),
                    catalog_path: operation_catalog_path(catalog_path, *helper_mode),
                    cache_dir: cache_dir.clone(),
                })
                .map_err(|e| e.to_string())?;
            let state_changed =
                !op_result.created_seq_ids.is_empty() || !op_result.changed_seq_ids.is_empty();
            ShellRunResult {
                state_changed,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::ReferenceExtractGene {
            helper_mode,
            genome_id,
            gene_query,
            occurrence,
            output_id,
            catalog_path,
            cache_dir,
        } => {
            let op_result = engine
                .apply(Operation::ExtractGenomeGene {
                    genome_id: genome_id.clone(),
                    gene_query: gene_query.clone(),
                    occurrence: *occurrence,
                    output_id: output_id.clone(),
                    catalog_path: operation_catalog_path(catalog_path, *helper_mode),
                    cache_dir: cache_dir.clone(),
                })
                .map_err(|e| e.to_string())?;
            let state_changed =
                !op_result.created_seq_ids.is_empty() || !op_result.changed_seq_ids.is_empty();
            ShellRunResult {
                state_changed,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::ReferenceExtendAnchor {
            helper_mode,
            seq_id,
            side,
            length_bp,
            output_id,
            catalog_path,
            cache_dir,
        } => {
            let op_result = engine
                .apply(Operation::ExtendGenomeAnchor {
                    seq_id: seq_id.clone(),
                    side: *side,
                    length_bp: *length_bp,
                    output_id: output_id.clone(),
                    catalog_path: operation_catalog_path(catalog_path, *helper_mode),
                    cache_dir: cache_dir.clone(),
                })
                .map_err(|e| e.to_string())?;
            let state_changed =
                !op_result.created_seq_ids.is_empty() || !op_result.changed_seq_ids.is_empty();
            ShellRunResult {
                state_changed,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::TracksImportBed {
            seq_id,
            path,
            track_name,
            min_score,
            max_score,
            clear_existing,
        } => {
            let op_result = engine
                .apply(Operation::ImportGenomeBedTrack {
                    seq_id: seq_id.clone(),
                    path: path.clone(),
                    track_name: track_name.clone(),
                    min_score: *min_score,
                    max_score: *max_score,
                    clear_existing: Some(*clear_existing),
                })
                .map_err(|e| e.to_string())?;
            let state_changed =
                !op_result.created_seq_ids.is_empty() || !op_result.changed_seq_ids.is_empty();
            ShellRunResult {
                state_changed,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::TracksImportBigWig {
            seq_id,
            path,
            track_name,
            min_score,
            max_score,
            clear_existing,
        } => {
            let op_result = engine
                .apply(Operation::ImportGenomeBigWigTrack {
                    seq_id: seq_id.clone(),
                    path: path.clone(),
                    track_name: track_name.clone(),
                    min_score: *min_score,
                    max_score: *max_score,
                    clear_existing: Some(*clear_existing),
                })
                .map_err(|e| e.to_string())?;
            let state_changed =
                !op_result.created_seq_ids.is_empty() || !op_result.changed_seq_ids.is_empty();
            ShellRunResult {
                state_changed,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::TracksImportVcf {
            seq_id,
            path,
            track_name,
            min_score,
            max_score,
            clear_existing,
        } => {
            let op_result = engine
                .apply(Operation::ImportGenomeVcfTrack {
                    seq_id: seq_id.clone(),
                    path: path.clone(),
                    track_name: track_name.clone(),
                    min_score: *min_score,
                    max_score: *max_score,
                    clear_existing: Some(*clear_existing),
                })
                .map_err(|e| e.to_string())?;
            let state_changed =
                !op_result.created_seq_ids.is_empty() || !op_result.changed_seq_ids.is_empty();
            ShellRunResult {
                state_changed,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::TracksTrackedList => {
            let subscriptions = engine.list_genome_track_subscriptions();
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "subscriptions": subscriptions,
                    "count": subscriptions.len()
                }),
            }
        }
        ShellCommand::TracksTrackedAdd { subscription } => {
            let inserted = engine
                .add_genome_track_subscription(subscription.clone())
                .map_err(|e| e.to_string())?;
            let subscriptions = engine.list_genome_track_subscriptions();
            ShellRunResult {
                state_changed: inserted,
                output: json!({
                    "inserted": inserted,
                    "subscription": subscription,
                    "subscriptions": subscriptions,
                    "count": subscriptions.len()
                }),
            }
        }
        ShellCommand::TracksTrackedRemove { index } => {
            let removed = engine
                .remove_genome_track_subscription(*index)
                .map_err(|e| e.to_string())?;
            let subscriptions = engine.list_genome_track_subscriptions();
            ShellRunResult {
                state_changed: true,
                output: json!({
                    "removed": removed,
                    "subscriptions": subscriptions,
                    "count": subscriptions.len()
                }),
            }
        }
        ShellCommand::TracksTrackedClear => {
            engine.clear_genome_track_subscriptions();
            ShellRunResult {
                state_changed: true,
                output: json!({
                    "cleared": true,
                    "subscriptions": [],
                    "count": 0
                }),
            }
        }
        ShellCommand::TracksTrackedApply {
            index,
            only_new_anchors,
        } => {
            let report = match index {
                Some(idx) => engine
                    .apply_tracked_genome_track_subscription(*idx)
                    .map_err(|e| e.to_string())?,
                None => engine
                    .sync_tracked_genome_track_subscriptions(*only_new_anchors)
                    .map_err(|e| e.to_string())?,
            };
            let state_changed = report.applied_imports > 0;
            ShellRunResult {
                state_changed,
                output: json!({
                    "report": report
                }),
            }
        }
        ShellCommand::CandidatesList => {
            let sets = engine.list_candidate_sets();
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": "gentle.candidate_sets.v1",
                    "set_count": sets.len(),
                    "sets": sets
                }),
            }
        }
        ShellCommand::CandidatesDelete { set_name } => {
            let before = engine
                .state()
                .metadata
                .get(CANDIDATE_SETS_METADATA_KEY)
                .cloned();
            let op_result = engine
                .apply(Operation::DeleteCandidateSet {
                    set_name: set_name.clone(),
                })
                .map_err(|e| e.to_string())?;
            let remaining = engine.list_candidate_sets().len();
            let after = engine
                .state()
                .metadata
                .get(CANDIDATE_SETS_METADATA_KEY)
                .cloned();
            let removed = before != after;
            ShellRunResult {
                state_changed: removed,
                output: json!({
                    "removed": removed,
                    "set_name": set_name,
                    "remaining_set_count": remaining,
                    "result": op_result
                }),
            }
        }
        ShellCommand::CandidatesGenerate {
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
            let before = engine
                .state()
                .metadata
                .get(CANDIDATE_SETS_METADATA_KEY)
                .cloned();
            let op_result = engine
                .apply(Operation::GenerateCandidateSet {
                    set_name: set_name.clone(),
                    seq_id: seq_id.clone(),
                    length_bp: *length_bp,
                    step_bp: *step_bp,
                    feature_kinds: feature_kinds.clone(),
                    feature_label_regex: feature_label_regex.clone(),
                    max_distance_bp: *max_distance_bp,
                    feature_geometry_mode: *feature_geometry_mode,
                    feature_boundary_mode: *feature_boundary_mode,
                    feature_strand_relation: *feature_strand_relation,
                    limit: Some(*limit),
                })
                .map_err(|e| e.to_string())?;
            let after = engine
                .state()
                .metadata
                .get(CANDIDATE_SETS_METADATA_KEY)
                .cloned();
            ShellRunResult {
                state_changed: before != after,
                output: json!({
                    "set_name": set_name,
                    "seq_id": seq_id,
                    "feature_geometry_mode": feature_geometry_mode.map(|mode| mode.as_str()),
                    "feature_boundary_mode": feature_boundary_mode.map(|mode| mode.as_str()),
                    "feature_strand_relation": feature_strand_relation.map(|mode| mode.as_str()),
                    "result": op_result
                }),
            }
        }
        ShellCommand::CandidatesGenerateBetweenAnchors {
            set_name,
            seq_id,
            anchor_a,
            anchor_b,
            length_bp,
            step_bp,
            limit,
        } => {
            let before = engine
                .state()
                .metadata
                .get(CANDIDATE_SETS_METADATA_KEY)
                .cloned();
            let op_result = engine
                .apply(Operation::GenerateCandidateSetBetweenAnchors {
                    set_name: set_name.clone(),
                    seq_id: seq_id.clone(),
                    anchor_a: anchor_a.clone(),
                    anchor_b: anchor_b.clone(),
                    length_bp: *length_bp,
                    step_bp: *step_bp,
                    limit: Some(*limit),
                })
                .map_err(|e| e.to_string())?;
            let after = engine
                .state()
                .metadata
                .get(CANDIDATE_SETS_METADATA_KEY)
                .cloned();
            ShellRunResult {
                state_changed: before != after,
                output: json!({
                    "set_name": set_name,
                    "seq_id": seq_id,
                    "anchor_a": anchor_a,
                    "anchor_b": anchor_b,
                    "result": op_result
                }),
            }
        }
        ShellCommand::CandidatesShow {
            set_name,
            limit,
            offset,
        } => {
            let (page, total, clamped_offset) = engine
                .inspect_candidate_set_page(set_name, *limit, *offset)
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "set_name": page.name,
                    "candidate_count": total,
                    "offset": clamped_offset,
                    "limit": limit,
                    "returned": page.candidates.len(),
                    "source_seq_ids": page.source_seq_ids,
                    "created_at_unix_ms": page.created_at_unix_ms,
                    "rows": page.candidates
                }),
            }
        }
        ShellCommand::CandidatesMetrics { set_name } => {
            let metrics = engine
                .list_candidate_set_metrics(set_name)
                .map_err(|e| e.to_string())?;
            let candidate_count = engine
                .list_candidate_sets()
                .into_iter()
                .find(|summary| summary.name == *set_name)
                .map(|summary| summary.candidate_count)
                .unwrap_or(0);
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "set_name": set_name,
                    "candidate_count": candidate_count,
                    "metric_count": metrics.len(),
                    "metrics": metrics
                }),
            }
        }
        ShellCommand::CandidatesScoreExpression {
            set_name,
            metric,
            expression,
        } => {
            let before = engine
                .state()
                .metadata
                .get(CANDIDATE_SETS_METADATA_KEY)
                .cloned();
            let op_result = engine
                .apply(Operation::ScoreCandidateSetExpression {
                    set_name: set_name.clone(),
                    metric: metric.clone(),
                    expression: expression.clone(),
                })
                .map_err(|e| e.to_string())?;
            let after = engine
                .state()
                .metadata
                .get(CANDIDATE_SETS_METADATA_KEY)
                .cloned();
            ShellRunResult {
                state_changed: before != after,
                output: json!({
                    "set_name": set_name,
                    "metric": metric,
                    "expression": expression,
                    "result": op_result
                }),
            }
        }
        ShellCommand::CandidatesScoreDistance {
            set_name,
            metric,
            feature_kinds,
            feature_label_regex,
            feature_geometry_mode,
            feature_boundary_mode,
            feature_strand_relation,
        } => {
            let before = engine
                .state()
                .metadata
                .get(CANDIDATE_SETS_METADATA_KEY)
                .cloned();
            let op_result = engine
                .apply(Operation::ScoreCandidateSetDistance {
                    set_name: set_name.clone(),
                    metric: metric.clone(),
                    feature_kinds: feature_kinds.clone(),
                    feature_label_regex: feature_label_regex.clone(),
                    feature_geometry_mode: *feature_geometry_mode,
                    feature_boundary_mode: *feature_boundary_mode,
                    feature_strand_relation: *feature_strand_relation,
                })
                .map_err(|e| e.to_string())?;
            let after = engine
                .state()
                .metadata
                .get(CANDIDATE_SETS_METADATA_KEY)
                .cloned();
            ShellRunResult {
                state_changed: before != after,
                output: json!({
                    "set_name": set_name,
                    "metric": metric,
                    "feature_kind_filter": feature_kinds,
                    "feature_label_regex": feature_label_regex,
                    "feature_geometry_mode": feature_geometry_mode.map(|mode| mode.as_str()),
                    "feature_boundary_mode": feature_boundary_mode.map(|mode| mode.as_str()),
                    "feature_strand_relation": feature_strand_relation.map(|mode| mode.as_str()),
                    "result": op_result
                }),
            }
        }
        ShellCommand::CandidatesFilter {
            input_set,
            output_set,
            metric,
            min,
            max,
            min_quantile,
            max_quantile,
        } => {
            let before = engine
                .state()
                .metadata
                .get(CANDIDATE_SETS_METADATA_KEY)
                .cloned();
            let op_result = engine
                .apply(Operation::FilterCandidateSet {
                    input_set: input_set.clone(),
                    output_set: output_set.clone(),
                    metric: metric.clone(),
                    min: *min,
                    max: *max,
                    min_quantile: *min_quantile,
                    max_quantile: *max_quantile,
                })
                .map_err(|e| e.to_string())?;
            let after = engine
                .state()
                .metadata
                .get(CANDIDATE_SETS_METADATA_KEY)
                .cloned();
            ShellRunResult {
                state_changed: before != after,
                output: json!({
                    "input_set": input_set,
                    "output_set": output_set,
                    "metric": metric,
                    "min_quantile": min_quantile,
                    "max_quantile": max_quantile,
                    "result": op_result
                }),
            }
        }
        ShellCommand::CandidatesSetOp {
            op,
            left_set,
            right_set,
            output_set,
        } => {
            let before = engine
                .state()
                .metadata
                .get(CANDIDATE_SETS_METADATA_KEY)
                .cloned();
            let op_result = engine
                .apply(Operation::CandidateSetOp {
                    op: match op {
                        CandidateSetOperator::Union => crate::engine::CandidateSetOperator::Union,
                        CandidateSetOperator::Intersect => {
                            crate::engine::CandidateSetOperator::Intersect
                        }
                        CandidateSetOperator::Subtract => {
                            crate::engine::CandidateSetOperator::Subtract
                        }
                    },
                    left_set: left_set.clone(),
                    right_set: right_set.clone(),
                    output_set: output_set.clone(),
                })
                .map_err(|e| e.to_string())?;
            let after = engine
                .state()
                .metadata
                .get(CANDIDATE_SETS_METADATA_KEY)
                .cloned();
            ShellRunResult {
                state_changed: before != after,
                output: json!({
                    "operator": op.as_str(),
                    "left_set": left_set,
                    "right_set": right_set,
                    "output_set": output_set,
                    "result": op_result
                }),
            }
        }
        ShellCommand::CandidatesMacro {
            script,
            transactional,
        } => run_candidates_macro(engine, script, *transactional, options)?,
        ShellCommand::SetParameter { name, value_json } => {
            let raw = parse_json_payload(value_json)?;
            let value: serde_json::Value = serde_json::from_str(&raw)
                .map_err(|e| format!("Invalid JSON value for set-param '{}': {e}", name))?;
            let op_result = engine
                .apply(Operation::SetParameter {
                    name: name.clone(),
                    value,
                })
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: true,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::Op { payload } => {
            let json_text = parse_json_payload(payload)?;
            let op: Operation = serde_json::from_str(&json_text)
                .map_err(|e| format!("Invalid operation JSON: {e}"))?;
            let op_result = engine.apply(op).map_err(|e| e.to_string())?;
            let state_changed =
                !op_result.created_seq_ids.is_empty() || !op_result.changed_seq_ids.is_empty();
            ShellRunResult {
                state_changed,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::Workflow { payload } => {
            let json_text = parse_json_payload(payload)?;
            let workflow: Workflow = serde_json::from_str(&json_text)
                .map_err(|e| format!("Invalid workflow JSON: {e}"))?;
            let results = engine.apply_workflow(workflow).map_err(|e| e.to_string())?;
            let state_changed = results
                .iter()
                .any(|r| !r.created_seq_ids.is_empty() || !r.changed_seq_ids.is_empty());
            ShellRunResult {
                state_changed,
                output: json!({ "results": results }),
            }
        }
    };
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dna_sequence::DNAsequence;
    use gb_io::seq::{Feature, FeatureKind, Location};
    use std::fs;
    use tempfile::tempdir;

    #[test]
    fn parse_help_with_topic_and_options() {
        let cmd = parse_shell_line("help candidates generate --format json --interface cli-shell")
            .expect("parse help with topic/options");
        match cmd {
            ShellCommand::Help {
                topic,
                format,
                interface_filter,
            } => {
                assert_eq!(
                    topic,
                    vec!["candidates".to_string(), "generate".to_string()]
                );
                assert_eq!(format, HelpOutputFormat::Json);
                assert_eq!(interface_filter.as_deref(), Some("cli-shell"));
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn execute_help_topic_json() {
        let mut engine = GentleEngine::from_state(ProjectState::default());
        let out = execute_shell_command(
            &mut engine,
            &ShellCommand::Help {
                topic: vec!["candidates".to_string(), "score-distance".to_string()],
                format: HelpOutputFormat::Json,
                interface_filter: None,
            },
        )
        .expect("execute help topic json");
        assert!(!out.state_changed);
        assert_eq!(
            out.output["doc"]["path"].as_str(),
            Some("candidates score-distance")
        );
    }

    #[test]
    fn parse_workflow_payload_keeps_whitespace() {
        let cmd = parse_shell_line("workflow '{ \"run_id\": \"x\", \"ops\": [] }'")
            .expect("workflow command parse");
        match cmd {
            ShellCommand::Workflow { payload } => {
                assert!(payload.contains("\"run_id\""));
                assert!(payload.contains("\"ops\""));
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_render_pool_gel_with_ladders() {
        let cmd = parse_shell_line("render-pool-gel-svg a,b out.svg --ladders 1kb,100bp")
            .expect("parse command");
        match cmd {
            ShellCommand::RenderPoolGelSvg {
                inputs,
                output,
                ladders,
            } => {
                assert_eq!(inputs, vec!["a".to_string(), "b".to_string()]);
                assert_eq!(output, "out.svg".to_string());
                assert_eq!(ladders, Some(vec!["1kb".to_string(), "100bp".to_string()]));
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_render_rna_svg() {
        let cmd = parse_shell_line("render-rna-svg rna_seq rna.svg").expect("parse command");
        match cmd {
            ShellCommand::RenderRnaSvg { seq_id, output } => {
                assert_eq!(seq_id, "rna_seq".to_string());
                assert_eq!(output, "rna.svg".to_string());
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_rna_info() {
        let cmd = parse_shell_line("rna-info rna_seq").expect("parse command");
        match cmd {
            ShellCommand::RnaInfo { seq_id } => {
                assert_eq!(seq_id, "rna_seq".to_string());
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_ladders_list_with_filter() {
        let cmd = parse_shell_line("ladders list --filter NEB").expect("parse command");
        match cmd {
            ShellCommand::LaddersList {
                molecule,
                name_filter,
            } => {
                assert_eq!(molecule, LadderMolecule::Dna);
                assert_eq!(name_filter, Some("NEB".to_string()));
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_ladders_export_with_filter() {
        let cmd =
            parse_shell_line("ladders export ladders.json --filter ruler").expect("parse command");
        match cmd {
            ShellCommand::LaddersExport {
                molecule,
                output,
                name_filter,
            } => {
                assert_eq!(molecule, LadderMolecule::Dna);
                assert_eq!(output, "ladders.json".to_string());
                assert_eq!(name_filter, Some("ruler".to_string()));
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_ladders_list_rna() {
        let cmd = parse_shell_line("ladders list --molecule rna --filter ss").expect("parse");
        match cmd {
            ShellCommand::LaddersList {
                molecule,
                name_filter,
            } => {
                assert_eq!(molecule, LadderMolecule::Rna);
                assert_eq!(name_filter, Some("ss".to_string()));
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_reference_genes_with_regex_and_biotypes() {
        let cmd = parse_shell_line(
            "helpers genes Helper --filter '^bla$' --biotype promoter --biotype cds --limit 10 --offset 3",
        )
        .expect("parse command");
        match cmd {
            ShellCommand::ReferenceGenes {
                helper_mode,
                genome_id,
                filter,
                biotypes,
                limit,
                offset,
                ..
            } => {
                assert!(helper_mode);
                assert_eq!(genome_id, "Helper");
                assert_eq!(filter, "^bla$");
                assert_eq!(biotypes, vec!["promoter".to_string(), "cds".to_string()]);
                assert_eq!(limit, 10);
                assert_eq!(offset, 3);
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_genomes_blast_with_options() {
        let cmd = parse_shell_line(
            "genomes blast ToyGenome ACGTACGT --max-hits 12 --task blastn --catalog c.json --cache-dir cache",
        )
        .expect("parse command");
        match cmd {
            ShellCommand::ReferenceBlast {
                helper_mode,
                genome_id,
                query_sequence,
                max_hits,
                task,
                catalog_path,
                cache_dir,
            } => {
                assert!(!helper_mode);
                assert_eq!(genome_id, "ToyGenome");
                assert_eq!(query_sequence, "ACGTACGT");
                assert_eq!(max_hits, 12);
                assert_eq!(task.as_deref(), Some("blastn"));
                assert_eq!(catalog_path.as_deref(), Some("c.json"));
                assert_eq!(cache_dir.as_deref(), Some("cache"));
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_helpers_blast_defaults() {
        let cmd = parse_shell_line("helpers blast pUC19 ACGTAG").expect("parse command");
        match cmd {
            ShellCommand::ReferenceBlast {
                helper_mode,
                genome_id,
                query_sequence,
                max_hits,
                task,
                ..
            } => {
                assert!(helper_mode);
                assert_eq!(genome_id, "pUC19");
                assert_eq!(query_sequence, "ACGTAG");
                assert_eq!(max_hits, 25);
                assert!(task.is_none());
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_genomes_blast_track_with_options() {
        let cmd = parse_shell_line(
            "genomes blast-track ToyGenome ACGTACGT query_seq --max-hits 8 --task blastn --track-name Hits --clear-existing --catalog c.json --cache-dir cache",
        )
        .expect("parse command");
        match cmd {
            ShellCommand::ReferenceBlastTrack {
                helper_mode,
                genome_id,
                query_sequence,
                target_seq_id,
                max_hits,
                task,
                track_name,
                clear_existing,
                catalog_path,
                cache_dir,
            } => {
                assert!(!helper_mode);
                assert_eq!(genome_id, "ToyGenome");
                assert_eq!(query_sequence, "ACGTACGT");
                assert_eq!(target_seq_id, "query_seq");
                assert_eq!(max_hits, 8);
                assert_eq!(task.as_deref(), Some("blastn"));
                assert_eq!(track_name.as_deref(), Some("Hits"));
                assert!(clear_existing);
                assert_eq!(catalog_path.as_deref(), Some("c.json"));
                assert_eq!(cache_dir.as_deref(), Some("cache"));
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_helpers_blast_track_defaults() {
        let cmd =
            parse_shell_line("helpers blast-track pUC19 ACGTAG target_seq").expect("parse command");
        match cmd {
            ShellCommand::ReferenceBlastTrack {
                helper_mode,
                genome_id,
                query_sequence,
                target_seq_id,
                max_hits,
                task,
                track_name,
                clear_existing,
                ..
            } => {
                assert!(helper_mode);
                assert_eq!(genome_id, "pUC19");
                assert_eq!(query_sequence, "ACGTAG");
                assert_eq!(target_seq_id, "target_seq");
                assert_eq!(max_hits, 25);
                assert!(task.is_none());
                assert!(track_name.is_none());
                assert!(!clear_existing);
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_tracks_import_bed() {
        let cmd = parse_shell_line(
            "tracks import-bed toy_slice test_files/data/peaks.bed.gz --name ChIP --min-score 5 --max-score 50 --clear-existing",
        )
        .expect("parse command");
        match cmd {
            ShellCommand::TracksImportBed {
                seq_id,
                path,
                track_name,
                min_score,
                max_score,
                clear_existing,
            } => {
                assert_eq!(seq_id, "toy_slice".to_string());
                assert_eq!(path, "test_files/data/peaks.bed.gz".to_string());
                assert_eq!(track_name, Some("ChIP".to_string()));
                assert_eq!(min_score, Some(5.0));
                assert_eq!(max_score, Some(50.0));
                assert!(clear_existing);
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_tracks_import_bigwig() {
        let cmd = parse_shell_line(
            "tracks import-bigwig toy_slice test_files/data/signal.bw --name RNA --min-score 0.5 --max-score 2.5 --clear-existing",
        )
        .expect("parse command");
        match cmd {
            ShellCommand::TracksImportBigWig {
                seq_id,
                path,
                track_name,
                min_score,
                max_score,
                clear_existing,
            } => {
                assert_eq!(seq_id, "toy_slice".to_string());
                assert_eq!(path, "test_files/data/signal.bw".to_string());
                assert_eq!(track_name, Some("RNA".to_string()));
                assert_eq!(min_score, Some(0.5));
                assert_eq!(max_score, Some(2.5));
                assert!(clear_existing);
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_tracks_import_vcf() {
        let cmd = parse_shell_line(
            "tracks import-vcf toy_slice test_files/data/variants.vcf.gz --name SNPs --min-score 10 --max-score 60 --clear-existing",
        )
        .expect("parse command");
        match cmd {
            ShellCommand::TracksImportVcf {
                seq_id,
                path,
                track_name,
                min_score,
                max_score,
                clear_existing,
            } => {
                assert_eq!(seq_id, "toy_slice".to_string());
                assert_eq!(path, "test_files/data/variants.vcf.gz".to_string());
                assert_eq!(track_name, Some("SNPs".to_string()));
                assert_eq!(min_score, Some(10.0));
                assert_eq!(max_score, Some(60.0));
                assert!(clear_existing);
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_tracks_tracked_add_with_options() {
        let cmd = parse_shell_line(
            "tracks tracked add test_files/data/signal.bw --source bed --name ChIP --min-score 0.5 --max-score 2.5 --clear-existing",
        )
        .expect("parse command");
        match cmd {
            ShellCommand::TracksTrackedAdd { subscription } => {
                assert_eq!(subscription.source, GenomeTrackSource::Bed);
                assert_eq!(subscription.path, "test_files/data/signal.bw");
                assert_eq!(subscription.track_name.as_deref(), Some("ChIP"));
                assert_eq!(subscription.min_score, Some(0.5));
                assert_eq!(subscription.max_score, Some(2.5));
                assert!(subscription.clear_existing);
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_tracks_tracked_add_with_vcf_source() {
        let cmd = parse_shell_line(
            "tracks tracked add test_files/data/variants.vcf.gz --source vcf --name Variants",
        )
        .expect("parse command");
        match cmd {
            ShellCommand::TracksTrackedAdd { subscription } => {
                assert_eq!(subscription.source, GenomeTrackSource::Vcf);
                assert_eq!(subscription.path, "test_files/data/variants.vcf.gz");
                assert_eq!(subscription.track_name.as_deref(), Some("Variants"));
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_tracks_tracked_apply_with_index() {
        let cmd = parse_shell_line("tracks tracked apply --index 3 --only-new-anchors")
            .expect("parse command");
        match cmd {
            ShellCommand::TracksTrackedApply {
                index,
                only_new_anchors,
            } => {
                assert_eq!(index, Some(3));
                assert!(only_new_anchors);
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_tracks_tracked_sync_alias() {
        let cmd = parse_shell_line("tracks tracked sync").expect("parse command");
        match cmd {
            ShellCommand::TracksTrackedApply {
                index,
                only_new_anchors,
            } => {
                assert_eq!(index, None);
                assert!(only_new_anchors);
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_candidates_generate_with_feature_filters() {
        let cmd = parse_shell_line(
            "candidates generate set1 seqA --length 20 --step 5 --feature-kind gene --feature-kind CDS --feature-label-regex '^TP53$' --max-distance 100 --limit 50",
        )
        .expect("parse candidates generate");
        match cmd {
            ShellCommand::CandidatesGenerate {
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
                assert_eq!(set_name, "set1");
                assert_eq!(seq_id, "seqA");
                assert_eq!(length_bp, 20);
                assert_eq!(step_bp, 5);
                assert_eq!(feature_kinds, vec!["gene".to_string(), "CDS".to_string()]);
                assert_eq!(feature_label_regex.as_deref(), Some("^TP53$"));
                assert_eq!(max_distance_bp, Some(100));
                assert_eq!(feature_geometry_mode, None);
                assert_eq!(feature_boundary_mode, None);
                assert_eq!(feature_strand_relation, None);
                assert_eq!(limit, 50);
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_candidates_generate_with_strand_relation() {
        let cmd =
            parse_shell_line("candidates generate set1 seqA --length 20 --strand-relation same")
                .expect("parse candidates generate with strand relation");
        match cmd {
            ShellCommand::CandidatesGenerate {
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
                assert_eq!(set_name, "set1");
                assert_eq!(seq_id, "seqA");
                assert_eq!(length_bp, 20);
                assert_eq!(step_bp, 1);
                assert!(feature_kinds.is_empty());
                assert_eq!(feature_label_regex, None);
                assert_eq!(max_distance_bp, None);
                assert_eq!(feature_geometry_mode, None);
                assert_eq!(feature_boundary_mode, None);
                assert_eq!(
                    feature_strand_relation,
                    Some(CandidateFeatureStrandRelation::Same)
                );
                assert_eq!(limit, DEFAULT_CANDIDATE_SET_LIMIT);
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_candidates_generate_between_anchors_with_positions() {
        let cmd = parse_shell_line(
            "candidates generate-between-anchors set1 seqA --length 20 --anchor-a-pos 5 --anchor-b-pos 105 --step 5 --limit 50",
        )
        .expect("parse candidates generate-between-anchors");
        match cmd {
            ShellCommand::CandidatesGenerateBetweenAnchors {
                set_name,
                seq_id,
                anchor_a,
                anchor_b,
                length_bp,
                step_bp,
                limit,
            } => {
                assert_eq!(set_name, "set1");
                assert_eq!(seq_id, "seqA");
                assert_eq!(length_bp, 20);
                assert_eq!(step_bp, 5);
                assert_eq!(limit, 50);
                assert_eq!(anchor_a, SequenceAnchor::Position { zero_based: 5 });
                assert_eq!(anchor_b, SequenceAnchor::Position { zero_based: 105 });
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_candidates_score_distance_with_geometry_and_boundary() {
        let cmd = parse_shell_line(
            "candidates score-distance set1 dist --feature-kind exon --feature-geometry feature_boundaries --feature-boundary five_prime",
        )
        .expect("parse candidates score-distance");
        match cmd {
            ShellCommand::CandidatesScoreDistance {
                set_name,
                metric,
                feature_kinds,
                feature_label_regex,
                feature_geometry_mode,
                feature_boundary_mode,
                feature_strand_relation,
            } => {
                assert_eq!(set_name, "set1");
                assert_eq!(metric, "dist");
                assert_eq!(feature_kinds, vec!["exon".to_string()]);
                assert_eq!(feature_label_regex, None);
                assert_eq!(
                    feature_geometry_mode,
                    Some(CandidateFeatureGeometryMode::FeatureBoundaries)
                );
                assert_eq!(
                    feature_boundary_mode,
                    Some(CandidateFeatureBoundaryMode::FivePrime)
                );
                assert_eq!(feature_strand_relation, None);
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_candidates_score_distance_with_strand_relation() {
        let cmd = parse_shell_line(
            "candidates score-distance set1 dist --feature-kind gene --strand-relation opposite",
        )
        .expect("parse candidates score-distance with strand relation");
        match cmd {
            ShellCommand::CandidatesScoreDistance {
                set_name,
                metric,
                feature_kinds,
                feature_label_regex,
                feature_geometry_mode,
                feature_boundary_mode,
                feature_strand_relation,
            } => {
                assert_eq!(set_name, "set1");
                assert_eq!(metric, "dist");
                assert_eq!(feature_kinds, vec!["gene".to_string()]);
                assert_eq!(feature_label_regex, None);
                assert_eq!(feature_geometry_mode, None);
                assert_eq!(feature_boundary_mode, None);
                assert_eq!(
                    feature_strand_relation,
                    Some(CandidateFeatureStrandRelation::Opposite)
                );
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_candidates_score_expression_preserves_formula() {
        let cmd = parse_shell_line(
            "candidates score set1 custom_score '(gc_fraction*100)+(length_bp/2)'",
        )
        .expect("parse candidates score");
        match cmd {
            ShellCommand::CandidatesScoreExpression {
                set_name,
                metric,
                expression,
            } => {
                assert_eq!(set_name, "set1");
                assert_eq!(metric, "custom_score");
                assert_eq!(expression, "(gc_fraction*100)+(length_bp/2)");
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_candidates_macro_file_reference() {
        let cmd = parse_shell_line("candidates macro --file test_files/candidates_plan.gsh")
            .expect("parse candidates macro --file");
        match cmd {
            ShellCommand::CandidatesMacro {
                script,
                transactional,
            } => {
                assert_eq!(script, "@test_files/candidates_plan.gsh");
                assert!(!transactional);
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_candidates_macro_transactional_flag() {
        let cmd = parse_shell_line(
            "candidates macro --transactional --file test_files/candidates_plan.gsh",
        )
        .expect("parse transactional candidates macro");
        match cmd {
            ShellCommand::CandidatesMacro {
                script,
                transactional,
            } => {
                assert_eq!(script, "@test_files/candidates_plan.gsh");
                assert!(transactional);
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn split_candidates_macro_statements_handles_quotes_and_comments() {
        let script = r#"
# comment line
generate set1 seqA --length 4 --step 2;
score set1 score "(gc_fraction*100)+1";
filter set1 set2 --metric score --min 10
"#;
        let statements =
            split_candidates_macro_statements(script).expect("split candidates macro statements");
        assert_eq!(statements.len(), 3);
        assert_eq!(statements[0], "generate set1 seqA --length 4 --step 2");
        assert_eq!(statements[1], "score set1 score \"(gc_fraction*100)+1\"");
        assert_eq!(statements[2], "filter set1 set2 --metric score --min 10");
    }

    #[test]
    fn execute_candidates_generate_score_distance_and_filter() {
        let mut state = ProjectState::default();
        let mut dna = DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence");
        dna.features_mut().push(Feature {
            kind: FeatureKind::from("gene"),
            location: Location::simple_range(5, 6),
            qualifiers: vec![("label".into(), Some("target".to_string()))],
        });
        state.sequences.insert("seqA".to_string(), dna);
        let mut engine = GentleEngine::from_state(state);

        let generated = execute_shell_command(
            &mut engine,
            &ShellCommand::CandidatesGenerate {
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
                limit: 10,
            },
        )
        .expect("generate candidates");
        assert!(generated.state_changed);
        assert_eq!(generated.output["set_name"].as_str(), Some("windows"));
        assert!(generated.output["result"]["messages"]
            .as_array()
            .map(|messages| !messages.is_empty())
            .unwrap_or(false));

        let score_distance = execute_shell_command(
            &mut engine,
            &ShellCommand::CandidatesScoreDistance {
                set_name: "windows".to_string(),
                metric: "dist_gene".to_string(),
                feature_kinds: vec!["gene".to_string()],
                feature_label_regex: None,
                feature_geometry_mode: None,
                feature_boundary_mode: None,
                feature_strand_relation: None,
            },
        )
        .expect("score distance");
        assert!(score_distance.state_changed);

        let filtered = execute_shell_command(
            &mut engine,
            &ShellCommand::CandidatesFilter {
                input_set: "windows".to_string(),
                output_set: "near_gene".to_string(),
                metric: "dist_gene".to_string(),
                min: None,
                max: Some(0.0),
                min_quantile: None,
                max_quantile: None,
            },
        )
        .expect("filter by distance");
        assert!(filtered.state_changed);
        assert_eq!(filtered.output["output_set"].as_str(), Some("near_gene"));

        let shown = execute_shell_command(
            &mut engine,
            &ShellCommand::CandidatesShow {
                set_name: "near_gene".to_string(),
                limit: 10,
                offset: 0,
            },
        )
        .expect("show near_gene");
        assert!(!shown.state_changed);
        assert_eq!(shown.output["candidate_count"].as_u64(), Some(1));
    }

    #[test]
    fn execute_candidates_generate_between_anchors_creates_set() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "seqA".to_string(),
            DNAsequence::from_sequence("ACGTACGTACGTACGT").expect("sequence"),
        );
        let mut engine = GentleEngine::from_state(state);

        let out = execute_shell_command(
            &mut engine,
            &ShellCommand::CandidatesGenerateBetweenAnchors {
                set_name: "between".to_string(),
                seq_id: "seqA".to_string(),
                anchor_a: SequenceAnchor::Position { zero_based: 2 },
                anchor_b: SequenceAnchor::Position { zero_based: 10 },
                length_bp: 4,
                step_bp: 2,
                limit: 64,
            },
        )
        .expect("generate candidates between anchors");
        assert!(out.state_changed);
        assert_eq!(out.output["set_name"].as_str(), Some("between"));

        let shown = execute_shell_command(
            &mut engine,
            &ShellCommand::CandidatesShow {
                set_name: "between".to_string(),
                limit: 64,
                offset: 0,
            },
        )
        .expect("show between");
        assert_eq!(shown.output["candidate_count"].as_u64(), Some(3));
    }

    #[test]
    fn execute_candidates_score_distance_honors_strand_relation() {
        let mut state = ProjectState::default();
        let mut dna = DNAsequence::from_sequence("ACGTACGTACGTACGTACGT").expect("sequence");
        dna.features_mut().push(Feature {
            kind: FeatureKind::from("gene"),
            location: Location::simple_range(2, 5),
            qualifiers: vec![("label".into(), Some("PLUS_GENE".to_string()))],
        });
        dna.features_mut().push(Feature {
            kind: FeatureKind::from("gene"),
            location: Location::Complement(Box::new(Location::simple_range(12, 15))),
            qualifiers: vec![("label".into(), Some("MINUS_GENE".to_string()))],
        });
        state.sequences.insert("seqA".to_string(), dna);
        let mut engine = GentleEngine::from_state(state);

        execute_shell_command(
            &mut engine,
            &ShellCommand::CandidatesGenerate {
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
                limit: 128,
            },
        )
        .expect("generate candidate windows");

        execute_shell_command(
            &mut engine,
            &ShellCommand::CandidatesScoreDistance {
                set_name: "windows".to_string(),
                metric: "dist_same".to_string(),
                feature_kinds: vec!["gene".to_string()],
                feature_label_regex: None,
                feature_geometry_mode: None,
                feature_boundary_mode: None,
                feature_strand_relation: Some(CandidateFeatureStrandRelation::Same),
            },
        )
        .expect("score distance same");
        execute_shell_command(
            &mut engine,
            &ShellCommand::CandidatesScoreDistance {
                set_name: "windows".to_string(),
                metric: "dist_opposite".to_string(),
                feature_kinds: vec!["gene".to_string()],
                feature_label_regex: None,
                feature_geometry_mode: None,
                feature_boundary_mode: None,
                feature_strand_relation: Some(CandidateFeatureStrandRelation::Opposite),
            },
        )
        .expect("score distance opposite");

        let (page, _, _) = engine
            .inspect_candidate_set_page("windows", 1024, 0)
            .expect("inspect scored windows");
        let at_pos = |pos: usize, metric: &str| -> f64 {
            page.candidates
                .iter()
                .find(|candidate| candidate.start_0based == pos)
                .and_then(|candidate| candidate.metrics.get(metric).copied())
                .unwrap_or(f64::NAN)
        };

        let plus_same = at_pos(2, "dist_same");
        let plus_opposite = at_pos(2, "dist_opposite");
        assert_eq!(plus_same, 0.0);
        assert!(plus_opposite > 0.0);

        let minus_same = at_pos(14, "dist_same");
        let minus_opposite = at_pos(14, "dist_opposite");
        assert!(minus_same > 0.0);
        assert_eq!(minus_opposite, 0.0);
    }

    #[test]
    fn execute_candidates_set_operations() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "seqA".to_string(),
            DNAsequence::from_sequence("ACGTACGT").expect("sequence"),
        );
        let mut engine = GentleEngine::from_state(state);

        execute_shell_command(
            &mut engine,
            &ShellCommand::CandidatesGenerate {
                set_name: "left".to_string(),
                seq_id: "seqA".to_string(),
                length_bp: 4,
                step_bp: 4,
                feature_kinds: vec![],
                feature_label_regex: None,
                max_distance_bp: None,
                feature_geometry_mode: None,
                feature_boundary_mode: None,
                feature_strand_relation: None,
                limit: 10,
            },
        )
        .expect("generate left");
        execute_shell_command(
            &mut engine,
            &ShellCommand::CandidatesGenerate {
                set_name: "right".to_string(),
                seq_id: "seqA".to_string(),
                length_bp: 4,
                step_bp: 2,
                feature_kinds: vec![],
                feature_label_regex: None,
                max_distance_bp: None,
                feature_geometry_mode: None,
                feature_boundary_mode: None,
                feature_strand_relation: None,
                limit: 10,
            },
        )
        .expect("generate right");

        let intersect = execute_shell_command(
            &mut engine,
            &ShellCommand::CandidatesSetOp {
                op: CandidateSetOperator::Intersect,
                left_set: "left".to_string(),
                right_set: "right".to_string(),
                output_set: "inter".to_string(),
            },
        )
        .expect("set intersect");
        assert!(intersect.state_changed);
        assert_eq!(intersect.output["operator"].as_str(), Some("intersect"));
        assert_eq!(intersect.output["output_set"].as_str(), Some("inter"));
    }

    #[test]
    fn execute_candidates_macro_runs_multiple_statements() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "seqA".to_string(),
            DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
        );
        let mut engine = GentleEngine::from_state(state);
        let out = execute_shell_command(
            &mut engine,
            &ShellCommand::CandidatesMacro {
                script: "generate win seqA --length 4 --step 2; score win x '(gc_fraction+1)'; filter win win_top --metric x --min-quantile 1.0".to_string(),
                transactional: false,
            },
        )
        .expect("execute candidates macro");
        assert!(out.state_changed);
        assert_eq!(out.output["executed"].as_u64(), Some(3));
        assert_eq!(out.output["transactional"].as_bool(), Some(false));
        let sets = engine.list_candidate_sets();
        assert!(sets.iter().any(|s| s.name == "win"));
        assert!(sets.iter().any(|s| s.name == "win_top"));
    }

    #[test]
    fn execute_candidates_macro_transactional_rolls_back_on_error() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "seqA".to_string(),
            DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
        );
        let mut engine = GentleEngine::from_state(state);
        let err = execute_shell_command(
            &mut engine,
            &ShellCommand::CandidatesMacro {
                script: "generate win seqA --length 4 --step 2; filter missing out --metric gc_fraction --min 0.1".to_string(),
                transactional: true,
            },
        )
        .unwrap_err();
        assert!(err.contains("rolled back"));
        let sets = engine.list_candidate_sets();
        assert!(
            !sets.iter().any(|s| s.name == "win"),
            "transactional macro should roll back earlier statements"
        );
    }

    #[test]
    fn execute_candidates_macro_rejects_nested_macro() {
        let mut engine = GentleEngine::new();
        let err = execute_shell_command(
            &mut engine,
            &ShellCommand::CandidatesMacro {
                script: "macro x".to_string(),
                transactional: false,
            },
        )
        .unwrap_err();
        assert!(err.contains("Nested candidates macro"));
    }

    #[test]
    fn execute_tracks_tracked_add_and_list() {
        let mut engine = GentleEngine::from_state(ProjectState::default());
        let subscription = GenomeTrackSubscription {
            source: GenomeTrackSource::Bed,
            path: "test_files/data/peaks.bed.gz".to_string(),
            track_name: Some("ChIP".to_string()),
            min_score: Some(1.0),
            max_score: Some(10.0),
            clear_existing: true,
        };
        let add = execute_shell_command(
            &mut engine,
            &ShellCommand::TracksTrackedAdd {
                subscription: subscription.clone(),
            },
        )
        .expect("execute tracked add");
        assert!(add.state_changed);
        assert_eq!(add.output["inserted"].as_bool(), Some(true));
        assert_eq!(add.output["count"].as_u64(), Some(1));

        let list = execute_shell_command(&mut engine, &ShellCommand::TracksTrackedList)
            .expect("execute tracked list");
        assert!(!list.state_changed);
        assert_eq!(list.output["count"].as_u64(), Some(1));
    }

    #[test]
    fn execute_state_summary_returns_json() {
        let mut engine = GentleEngine::from_state(ProjectState::default());
        let out = execute_shell_command(&mut engine, &ShellCommand::StateSummary)
            .expect("execute state summary");
        assert!(!out.state_changed);
        assert!(out.output.get("sequence_count").is_some());
    }

    #[test]
    fn parse_resources_sync_rebase_with_flag() {
        let cmd = parse_shell_line(
            "resources sync-rebase https://example.org/rebase.withrefm out.json --commercial-only",
        )
        .expect("parse resources sync-rebase");
        match cmd {
            ShellCommand::ResourcesSyncRebase {
                input,
                output,
                commercial_only,
            } => {
                assert_eq!(input, "https://example.org/rebase.withrefm".to_string());
                assert_eq!(output, Some("out.json".to_string()));
                assert!(commercial_only);
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_resources_sync_jaspar_with_output() {
        let cmd = parse_shell_line("resources sync-jaspar motifs.pfm out.motifs.json")
            .expect("parse resources sync-jaspar");
        match cmd {
            ShellCommand::ResourcesSyncJaspar { input, output } => {
                assert_eq!(input, "motifs.pfm".to_string());
                assert_eq!(output, Some("out.motifs.json".to_string()));
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_import_pool_with_prefix() {
        let cmd = parse_shell_line("import-pool test_files/demo.pool.gentle.json imported")
            .expect("parse import-pool");
        match cmd {
            ShellCommand::ImportPool { input, prefix } => {
                assert_eq!(input, "test_files/demo.pool.gentle.json".to_string());
                assert_eq!(prefix, "imported".to_string());
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn execute_import_pool_loads_members_into_engine_state() {
        let td = tempdir().expect("tempdir");
        let pool_path = td.path().join("demo.pool.gentle.json");
        let pool_json = json!({
            "schema": "gentle.pool.v1",
            "pool_id": "demo_pool",
            "human_id": "demo",
            "member_count": 1,
            "members": [
                {
                    "seq_id": "member_1",
                    "human_id": "member_1",
                    "name": "Member One",
                    "sequence": "ATGCATGC",
                    "length_bp": 8,
                    "topology": "linear",
                    "ends": {
                        "end_type": "blunt",
                        "forward_5": "",
                        "forward_3": "",
                        "reverse_5": "",
                        "reverse_3": ""
                    }
                }
            ]
        });
        fs::write(
            &pool_path,
            serde_json::to_string_pretty(&pool_json).expect("serialize pool json"),
        )
        .expect("write pool json");

        let mut state = ProjectState::default();
        state.sequences.insert(
            "imported_1".to_string(),
            DNAsequence::from_sequence("AAAA").expect("seed sequence"),
        );
        let mut engine = GentleEngine::from_state(state);
        let out = execute_shell_command(
            &mut engine,
            &ShellCommand::ImportPool {
                input: pool_path.to_string_lossy().to_string(),
                prefix: "imported".to_string(),
            },
        )
        .expect("execute import-pool");

        assert!(out.state_changed);
        assert_eq!(out.output["pool_id"].as_str(), Some("demo_pool"));
        assert_eq!(out.output["member_count"].as_u64(), Some(1));
        let imported_ids = out
            .output
            .get("imported_ids")
            .and_then(|v| v.as_array())
            .expect("imported_ids array");
        assert_eq!(imported_ids.len(), 1);
        let imported_id = imported_ids[0].as_str().expect("imported id string");
        assert!(
            imported_id.starts_with("imported_1"),
            "import id should be a collision-resolved variant of imported_1, got {imported_id}"
        );
        assert_ne!(imported_id, "imported_1");
        assert!(engine.state().sequences.contains_key(imported_id));
    }

    #[test]
    fn parse_helpers_extract_gene() {
        let cmd = parse_shell_line(
            "helpers extract-gene pUC19 bla --occurrence 2 --output-id out --cache-dir cache",
        )
        .expect("parse helpers extract-gene");
        match cmd {
            ShellCommand::ReferenceExtractGene {
                helper_mode,
                genome_id,
                gene_query,
                occurrence,
                output_id,
                cache_dir,
                ..
            } => {
                assert!(helper_mode);
                assert_eq!(genome_id, "pUC19".to_string());
                assert_eq!(gene_query, "bla".to_string());
                assert_eq!(occurrence, Some(2));
                assert_eq!(output_id, Some("out".to_string()));
                assert_eq!(cache_dir, Some("cache".to_string()));
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_genomes_prepare_with_timeout() {
        let cmd = parse_shell_line(
            "genomes prepare ToyGenome --catalog c.json --cache-dir cache --timeout-secs 90",
        )
        .expect("parse genomes prepare");
        match cmd {
            ShellCommand::ReferencePrepare {
                helper_mode,
                genome_id,
                catalog_path,
                cache_dir,
                timeout_seconds,
            } => {
                assert!(!helper_mode);
                assert_eq!(genome_id, "ToyGenome".to_string());
                assert_eq!(catalog_path, Some("c.json".to_string()));
                assert_eq!(cache_dir, Some("cache".to_string()));
                assert_eq!(timeout_seconds, Some(90));
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_genomes_extend_anchor_five_prime() {
        let cmd = parse_shell_line(
            "genomes extend-anchor tp73 5p 150 --output-id tp73_ext --catalog c.json --cache-dir cache",
        )
        .expect("parse genomes extend-anchor");
        match cmd {
            ShellCommand::ReferenceExtendAnchor {
                helper_mode,
                seq_id,
                side,
                length_bp,
                output_id,
                catalog_path,
                cache_dir,
            } => {
                assert!(!helper_mode);
                assert_eq!(seq_id, "tp73".to_string());
                assert_eq!(side, GenomeAnchorSide::FivePrime);
                assert_eq!(length_bp, 150);
                assert_eq!(output_id, Some("tp73_ext".to_string()));
                assert_eq!(catalog_path, Some("c.json".to_string()));
                assert_eq!(cache_dir, Some("cache".to_string()));
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn execute_genomes_extend_anchor_creates_sequence() {
        let td = tempdir().expect("tempdir");
        let fasta = td.path().join("toy.fa");
        let gtf = td.path().join("toy.gtf");
        fs::write(&fasta, ">chr1\nACGTACGTACGT\n").expect("write fasta");
        fs::write(
            &gtf,
            "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n",
        )
        .expect("write gtf");
        let catalog = td.path().join("catalog.json");
        let cache_dir = td.path().join("cache");
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
        fs::write(&catalog, catalog_json).expect("write catalog");
        let catalog_path = catalog.to_string_lossy().to_string();

        let mut engine = GentleEngine::new();
        execute_shell_command(
            &mut engine,
            &ShellCommand::ReferencePrepare {
                helper_mode: false,
                genome_id: "ToyGenome".to_string(),
                catalog_path: Some(catalog_path.clone()),
                cache_dir: None,
                timeout_seconds: None,
            },
        )
        .expect("prepare genome");

        execute_shell_command(
            &mut engine,
            &ShellCommand::ReferenceExtractRegion {
                helper_mode: false,
                genome_id: "ToyGenome".to_string(),
                chromosome: "chr1".to_string(),
                start_1based: 3,
                end_1based: 10,
                output_id: Some("slice".to_string()),
                catalog_path: Some(catalog_path.clone()),
                cache_dir: None,
            },
        )
        .expect("extract region");

        let out = execute_shell_command(
            &mut engine,
            &ShellCommand::ReferenceExtendAnchor {
                helper_mode: false,
                seq_id: "slice".to_string(),
                side: GenomeAnchorSide::FivePrime,
                length_bp: 2,
                output_id: Some("slice_ext5".to_string()),
                catalog_path: Some(catalog_path),
                cache_dir: None,
            },
        )
        .expect("execute extend-anchor");

        assert!(out.state_changed);
        let created = out
            .output
            .get("result")
            .and_then(|v| v.get("created_seq_ids"))
            .and_then(|v| v.as_array())
            .expect("created_seq_ids");
        assert!(created
            .iter()
            .any(|v| v.as_str().map(|id| id == "slice_ext5").unwrap_or(false)));
        let seq = engine
            .state()
            .sequences
            .get("slice_ext5")
            .expect("extended sequence in state");
        assert_eq!(seq.get_forward_string(), "ACGTACGTAC");
    }

    #[test]
    fn parse_genomes_validate_catalog() {
        let cmd = parse_shell_line("genomes validate-catalog --catalog assets/genomes.json")
            .expect("parse command");
        match cmd {
            ShellCommand::ReferenceValidateCatalog {
                helper_mode,
                catalog_path,
            } => {
                assert!(!helper_mode);
                assert_eq!(catalog_path, Some("assets/genomes.json".to_string()));
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn execute_genomes_validate_catalog_reports_valid() {
        let td = tempdir().expect("tempdir");
        let fasta = td.path().join("toy.fa");
        let gtf = td.path().join("toy.gtf");
        let cache = td.path().join("cache");
        fs::write(&fasta, ">chr1\nACGT\n").expect("write fasta");
        fs::write(
            &gtf,
            "chr1\tsrc\tgene\t1\t4\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n",
        )
        .expect("write gtf");
        let catalog = td.path().join("catalog.json");
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
            cache.display()
        );
        fs::write(&catalog, catalog_json).expect("write catalog");
        let mut engine = GentleEngine::new();
        let out = execute_shell_command(
            &mut engine,
            &ShellCommand::ReferenceValidateCatalog {
                helper_mode: false,
                catalog_path: Some(catalog.to_string_lossy().to_string()),
            },
        )
        .expect("execute validate-catalog");
        assert!(!out.state_changed);
        assert_eq!(out.output["valid"].as_bool(), Some(true));
        assert_eq!(out.output["genome_count"].as_u64(), Some(1));
    }

    #[test]
    fn parse_set_param_command() {
        let cmd =
            parse_shell_line("set-param vcf_display_pass_only true").expect("parse set-param");
        match cmd {
            ShellCommand::SetParameter { name, value_json } => {
                assert_eq!(name, "vcf_display_pass_only");
                assert_eq!(value_json, "true");
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn execute_set_param_updates_display_state() {
        let mut engine = GentleEngine::new();
        let out = execute_shell_command(
            &mut engine,
            &ShellCommand::SetParameter {
                name: "vcf_display_min_qual".to_string(),
                value_json: "33.5".to_string(),
            },
        )
        .expect("execute set-param");
        assert!(out.state_changed);
        assert!(
            (engine.state().display.vcf_display_min_qual - 33.5).abs() < f64::EPSILON,
            "vcf_display_min_qual should be updated by set-param"
        );
    }
}
