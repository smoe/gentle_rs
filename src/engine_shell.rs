//! Shared shell command grammar and executor for GENtle.
//!
//! This module parses textual shell input into typed commands and executes them
//! against the shared engine. It is intentionally reused by both the GUI Shell
//! panel and `gentle_cli shell` so command behavior stays adapter-equivalent.
//!
//! Core responsibilities:
//! - Tokenize/parse shell syntax into structured `ShellCommand` values.
//! - Execute commands through `GentleEngine` and shared operation paths.
//! - Return deterministic, machine-readable command results and diagnostics.
//! - Provide help/introspection renderers backed by glossary-driven docs.
//! - Bridge adapter-facing command families (for example resources, agents,
//!   and UI intents) without duplicating core biology logic.
//!
//! Safety and consistency constraints:
//! - Suggested agent commands execute through the same parser/executor path.
//! - Recursive/nested `agents ask` execution from suggested commands is blocked.
//! - Screenshot-capture routes remain explicitly policy-gated.

use crate::{
    agent_bridge::{
        AGENT_BASE_URL_ENV, AGENT_CONNECT_TIMEOUT_SECS_ENV, AGENT_MAX_RESPONSE_BYTES_ENV,
        AGENT_MAX_RETRIES_ENV, AGENT_MODEL_ENV, AGENT_READ_TIMEOUT_SECS_ENV,
        AGENT_TIMEOUT_SECS_ENV, AgentExecutionIntent, agent_system_availability,
        invoke_agent_support_with_env_overrides, load_agent_system_catalog,
    },
    dna_ladder::LadderMolecule,
    dna_sequence::DNAsequence,
    engine::{
        CANDIDATE_MACRO_TEMPLATES_METADATA_KEY, CandidateFeatureBoundaryMode,
        CandidateFeatureGeometryMode, CandidateFeatureStrandRelation, CandidateMacroTemplateParam,
        CandidateObjectiveDirection, CandidateObjectiveSpec, CandidateTieBreakPolicy,
        CandidateWeightedObjectiveTerm, Engine, FeatureExpertTarget, GUIDE_DESIGN_METADATA_KEY,
        GenomeAnchorSide, GenomeTrackSource, GenomeTrackSubscription, GentleEngine, GuideCandidate,
        GuideOligoExportFormat, GuideOligoPlateFormat, GuidePracticalFilterConfig, Operation,
        ProjectState, RenderSvgMode, SequenceAnchor, WORKFLOW_MACRO_TEMPLATES_METADATA_KEY,
        Workflow, WorkflowMacroTemplateParam,
    },
    genomes::{
        DEFAULT_GENOME_CATALOG_PATH, DEFAULT_HELPER_GENOME_CATALOG_PATH, GenomeCatalog,
        GenomeGeneRecord,
    },
    resource_sync,
    shell_docs::{
        HelpOutputFormat, shell_help_json as render_shell_help_json,
        shell_help_markdown as render_shell_help_markdown,
        shell_help_text as render_shell_help_text,
        shell_topic_help_json as render_shell_topic_help_json,
        shell_topic_help_markdown as render_shell_topic_help_markdown,
        shell_topic_help_text as render_shell_topic_help_text,
    },
    tf_motifs,
};
#[cfg(all(target_os = "macos", feature = "screenshot-capture"))]
use objc2_app_kit::NSApplication;
#[cfg(all(target_os = "macos", feature = "screenshot-capture"))]
use objc2_foundation::MainThreadMarker;
use regex::{Regex, RegexBuilder};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
#[cfg(all(target_os = "macos", feature = "screenshot-capture"))]
use std::path::Path;
#[cfg(all(target_os = "macos", feature = "screenshot-capture"))]
use std::process::Command;
use std::{
    collections::{BTreeSet, HashMap},
    fs,
    path::{Path, PathBuf},
};

const CLONING_PATTERN_FILE_SCHEMA: &str = "gentle.cloning_patterns.v1";
const CLONING_PATTERN_TEMPLATE_FILE_SCHEMA: &str = "gentle.cloning_pattern_template.v1";
const CLONING_ROUTINE_CATALOG_SCHEMA: &str = "gentle.cloning_routines.v1";
const CLONING_ROUTINE_LIST_SCHEMA: &str = "gentle.cloning_routines_list.v1";
const DEFAULT_CLONING_ROUTINE_CATALOG_PATH: &str = "assets/cloning_routines.json";

#[derive(Debug, Clone, Deserialize)]
#[serde(default)]
struct CloningPatternFile {
    schema: String,
    templates: Vec<CloningPatternTemplate>,
}

impl Default for CloningPatternFile {
    fn default() -> Self {
        Self {
            schema: CLONING_PATTERN_FILE_SCHEMA.to_string(),
            templates: vec![],
        }
    }
}

#[derive(Debug, Clone, Deserialize, Default)]
#[serde(default)]
struct CloningPatternTemplate {
    name: String,
    description: Option<String>,
    details_url: Option<String>,
    parameters: Vec<WorkflowMacroTemplateParam>,
    script: String,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(default)]
struct CloningPatternTemplateFile {
    schema: String,
    name: String,
    description: Option<String>,
    details_url: Option<String>,
    parameters: Vec<WorkflowMacroTemplateParam>,
    script: String,
}

impl Default for CloningPatternTemplateFile {
    fn default() -> Self {
        Self {
            schema: CLONING_PATTERN_TEMPLATE_FILE_SCHEMA.to_string(),
            name: String::new(),
            description: None,
            details_url: None,
            parameters: vec![],
            script: String::new(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
struct CloningRoutineCatalog {
    schema: String,
    routines: Vec<CloningRoutineDefinition>,
}

impl Default for CloningRoutineCatalog {
    fn default() -> Self {
        Self {
            schema: CLONING_ROUTINE_CATALOG_SCHEMA.to_string(),
            routines: vec![],
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
struct CloningRoutinePort {
    port_id: String,
    kind: String,
    required: bool,
    cardinality: String,
    description: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
struct CloningRoutineDefinition {
    routine_id: String,
    title: String,
    family: String,
    status: String,
    vocabulary_tags: Vec<String>,
    summary: Option<String>,
    details_url: Option<String>,
    template_name: String,
    template_path: Option<String>,
    input_ports: Vec<CloningRoutinePort>,
    output_ports: Vec<CloningRoutinePort>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum UiIntentAction {
    Open,
    Focus,
}

impl UiIntentAction {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Open => "open",
            Self::Focus => "focus",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum UiIntentTarget {
    PreparedReferences,
    PrepareReferenceGenome,
    RetrieveGenomeSequence,
    BlastGenomeSequence,
    ImportGenomeTrack,
    AgentAssistant,
    PrepareHelperGenome,
    RetrieveHelperSequence,
    BlastHelperSequence,
}

impl UiIntentTarget {
    fn parse(raw: &str) -> Option<Self> {
        match raw.trim().to_ascii_lowercase().as_str() {
            "prepared-references" | "prepared_references" | "prepared" => {
                Some(Self::PreparedReferences)
            }
            "prepare-reference-genome"
            | "prepare_reference_genome"
            | "prepare-genome"
            | "prepare_genome" => Some(Self::PrepareReferenceGenome),
            "retrieve-genome-sequence"
            | "retrieve_genome_sequence"
            | "retrieve-genome"
            | "retrieve_genome" => Some(Self::RetrieveGenomeSequence),
            "blast-genome-sequence" | "blast_genome_sequence" | "blast-genome" | "blast_genome" => {
                Some(Self::BlastGenomeSequence)
            }
            "import-genome-track" | "import_genome_track" | "genome-track" | "tracks-import" => {
                Some(Self::ImportGenomeTrack)
            }
            "agent-assistant" | "agent_assistant" | "agent" => Some(Self::AgentAssistant),
            "prepare-helper-genome"
            | "prepare_helper_genome"
            | "prepare-helper"
            | "prepare_helper" => Some(Self::PrepareHelperGenome),
            "retrieve-helper-sequence"
            | "retrieve_helper_sequence"
            | "retrieve-helper"
            | "retrieve_helper" => Some(Self::RetrieveHelperSequence),
            "blast-helper-sequence" | "blast_helper_sequence" | "blast-helper" | "blast_helper" => {
                Some(Self::BlastHelperSequence)
            }
            _ => None,
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::PreparedReferences => "prepared-references",
            Self::PrepareReferenceGenome => "prepare-reference-genome",
            Self::RetrieveGenomeSequence => "retrieve-genome-sequence",
            Self::BlastGenomeSequence => "blast-genome-sequence",
            Self::ImportGenomeTrack => "import-genome-track",
            Self::AgentAssistant => "agent-assistant",
            Self::PrepareHelperGenome => "prepare-helper-genome",
            Self::RetrieveHelperSequence => "retrieve-helper-sequence",
            Self::BlastHelperSequence => "blast-helper-sequence",
        }
    }
}

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
    InspectFeatureExpert {
        seq_id: String,
        target: FeatureExpertTarget,
    },
    RenderFeatureExpertSvg {
        seq_id: String,
        target: FeatureExpertTarget,
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
        container_ids: Option<Vec<String>>,
        arrangement_id: Option<String>,
    },
    CreateArrangementSerial {
        container_ids: Vec<String>,
        arrangement_id: Option<String>,
        name: Option<String>,
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
    RoutinesList {
        catalog_path: Option<String>,
        family: Option<String>,
        status: Option<String>,
        tag: Option<String>,
        query: Option<String>,
    },
    AgentsList {
        catalog_path: Option<String>,
    },
    AgentsAsk {
        system_id: String,
        prompt: String,
        catalog_path: Option<String>,
        base_url_override: Option<String>,
        model_override: Option<String>,
        timeout_seconds: Option<u64>,
        connect_timeout_seconds: Option<u64>,
        read_timeout_seconds: Option<u64>,
        max_retries: Option<usize>,
        max_response_bytes: Option<usize>,
        include_state_summary: bool,
        allow_auto_exec: bool,
        execute_all: bool,
        execute_indices: Vec<usize>,
    },
    UiListIntents,
    UiIntent {
        action: UiIntentAction,
        target: UiIntentTarget,
        genome_id: Option<String>,
        helper_mode: bool,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
        filter: Option<String>,
        species: Option<String>,
        latest: bool,
    },
    UiPreparedGenomes {
        helper_mode: bool,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
        filter: Option<String>,
        species: Option<String>,
        latest: bool,
    },
    UiLatestPrepared {
        helper_mode: bool,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
        species: String,
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
    MacrosRun {
        script: String,
        transactional: bool,
    },
    MacrosTemplateList,
    MacrosTemplateShow {
        name: String,
    },
    MacrosTemplateUpsert {
        name: String,
        description: Option<String>,
        details_url: Option<String>,
        parameters: Vec<WorkflowMacroTemplateParam>,
        script: String,
    },
    MacrosTemplateDelete {
        name: String,
    },
    MacrosTemplateImport {
        path: String,
    },
    MacrosTemplateRun {
        name: String,
        bindings: HashMap<String, String>,
        transactional: bool,
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
    CandidatesScoreWeightedObjective {
        set_name: String,
        metric: String,
        objectives: Vec<CandidateWeightedObjectiveTerm>,
        normalize_metrics: bool,
    },
    CandidatesTopK {
        input_set: String,
        output_set: String,
        metric: String,
        k: usize,
        direction: CandidateObjectiveDirection,
        tie_break: CandidateTieBreakPolicy,
    },
    CandidatesParetoFrontier {
        input_set: String,
        output_set: String,
        objectives: Vec<CandidateObjectiveSpec>,
        max_candidates: Option<usize>,
        tie_break: CandidateTieBreakPolicy,
    },
    CandidatesMacro {
        script: String,
        transactional: bool,
    },
    CandidatesTemplateList,
    CandidatesTemplateShow {
        name: String,
    },
    CandidatesTemplateUpsert {
        name: String,
        description: Option<String>,
        details_url: Option<String>,
        parameters: Vec<CandidateMacroTemplateParam>,
        script: String,
    },
    CandidatesTemplateDelete {
        name: String,
    },
    CandidatesTemplateRun {
        name: String,
        bindings: HashMap<String, String>,
        transactional: bool,
    },
    GuidesList,
    GuidesShow {
        guide_set_id: String,
        limit: usize,
        offset: usize,
    },
    GuidesPut {
        guide_set_id: String,
        guides_json: String,
    },
    GuidesDelete {
        guide_set_id: String,
    },
    GuidesFilter {
        guide_set_id: String,
        config_json: Option<String>,
        output_guide_set_id: Option<String>,
    },
    GuidesFilterShow {
        guide_set_id: String,
    },
    GuidesOligosGenerate {
        guide_set_id: String,
        template_id: String,
        apply_5prime_g_extension: bool,
        output_oligo_set_id: Option<String>,
        passed_only: bool,
    },
    GuidesOligosList {
        guide_set_id: Option<String>,
    },
    GuidesOligosShow {
        oligo_set_id: String,
    },
    GuidesOligosExport {
        guide_set_id: String,
        oligo_set_id: Option<String>,
        format: GuideOligoExportFormat,
        path: String,
        plate_format: Option<GuideOligoPlateFormat>,
    },
    GuidesProtocolExport {
        guide_set_id: String,
        oligo_set_id: Option<String>,
        path: String,
        include_qc_checklist: bool,
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

#[derive(Debug, Clone, Copy)]
pub struct ShellExecutionOptions {
    pub allow_screenshots: bool,
    pub allow_agent_commands: bool,
}

impl Default for ShellExecutionOptions {
    fn default() -> Self {
        Self {
            allow_screenshots: false,
            allow_agent_commands: true,
        }
    }
}

impl ShellExecutionOptions {
    pub fn from_env() -> Self {
        let raw = std::env::var("GENTLE_ALLOW_SCREENSHOTS").unwrap_or_default();
        let allow_screenshots = cfg!(feature = "screenshot-capture")
            && matches!(
                raw.trim().to_ascii_lowercase().as_str(),
                "1" | "true" | "yes" | "on"
            );
        Self {
            allow_screenshots,
            allow_agent_commands: true,
        }
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
const DEFAULT_GUIDE_PAGE_SIZE: usize = 100;

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
        "three_prime" | "three-prime" | "3p" | "3'" => Ok(CandidateFeatureBoundaryMode::ThreePrime),
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

fn parse_candidate_objective_direction(raw: &str) -> Result<CandidateObjectiveDirection, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "max" | "maximize" => Ok(CandidateObjectiveDirection::Maximize),
        "min" | "minimize" => Ok(CandidateObjectiveDirection::Minimize),
        other => Err(format!(
            "Unsupported objective direction '{other}' (expected max|min)"
        )),
    }
}

fn parse_candidate_tie_break_policy(raw: &str) -> Result<CandidateTieBreakPolicy, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "seq_start_end" | "default" | "seq" => Ok(CandidateTieBreakPolicy::SeqStartEnd),
        "seq_end_start" => Ok(CandidateTieBreakPolicy::SeqEndStart),
        "length_ascending" | "length_asc" | "shortest" => {
            Ok(CandidateTieBreakPolicy::LengthAscending)
        }
        "length_descending" | "length_desc" | "longest" => {
            Ok(CandidateTieBreakPolicy::LengthDescending)
        }
        "sequence_lexicographic" | "sequence" | "lexicographic" => {
            Ok(CandidateTieBreakPolicy::SequenceLexicographic)
        }
        other => Err(format!(
            "Unsupported tie-break policy '{other}' (expected seq_start_end|seq_end_start|length_ascending|length_descending|sequence_lexicographic)"
        )),
    }
}

fn parse_weighted_objective_term(raw: &str) -> Result<CandidateWeightedObjectiveTerm, String> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return Err("Weighted objective term cannot be empty".to_string());
    }
    let parts = trimmed.split(':').collect::<Vec<_>>();
    if parts.len() < 2 || parts.len() > 3 {
        return Err(format!(
            "Invalid weighted objective term '{}'; expected METRIC:WEIGHT[:max|min]",
            trimmed
        ));
    }
    let metric = parts[0].trim();
    if metric.is_empty() {
        return Err(format!(
            "Invalid weighted objective term '{}': missing metric",
            trimmed
        ));
    }
    let weight_raw = parts[1].trim();
    let weight = weight_raw
        .parse::<f64>()
        .map_err(|e| format!("Invalid weighted objective weight '{}': {e}", weight_raw))?;
    let direction = if parts.len() == 3 {
        parse_candidate_objective_direction(parts[2])?
    } else {
        CandidateObjectiveDirection::Maximize
    };
    Ok(CandidateWeightedObjectiveTerm {
        metric: metric.to_string(),
        weight,
        direction,
    })
}

fn parse_pareto_objective(raw: &str) -> Result<CandidateObjectiveSpec, String> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return Err("Pareto objective cannot be empty".to_string());
    }
    let parts = trimmed.split(':').collect::<Vec<_>>();
    if parts.is_empty() || parts.len() > 2 {
        return Err(format!(
            "Invalid Pareto objective '{}'; expected METRIC[:max|min]",
            trimmed
        ));
    }
    let metric = parts[0].trim();
    if metric.is_empty() {
        return Err(format!(
            "Invalid Pareto objective '{}': missing metric",
            trimmed
        ));
    }
    let direction = if parts.len() == 2 {
        parse_candidate_objective_direction(parts[1])?
    } else {
        CandidateObjectiveDirection::Maximize
    };
    Ok(CandidateObjectiveSpec {
        metric: metric.to_string(),
        direction,
    })
}

fn parse_template_param_spec_parts(raw: &str) -> Result<(String, Option<String>, bool), String> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return Err("Template parameter cannot be empty".to_string());
    }
    if let Some((name, default_value)) = trimmed.split_once('=') {
        let name = name.trim();
        if name.is_empty() {
            return Err(format!(
                "Invalid template parameter '{}': missing parameter name",
                trimmed
            ));
        }
        Ok((name.to_string(), Some(default_value.to_string()), false))
    } else {
        Ok((trimmed.to_string(), None, true))
    }
}

fn parse_candidate_template_param_spec(raw: &str) -> Result<CandidateMacroTemplateParam, String> {
    let (name, default_value, required) = parse_template_param_spec_parts(raw)?;
    Ok(CandidateMacroTemplateParam {
        name,
        default_value,
        required,
    })
}

fn parse_workflow_template_param_spec(raw: &str) -> Result<WorkflowMacroTemplateParam, String> {
    let (name, default_value, required) = parse_template_param_spec_parts(raw)?;
    Ok(WorkflowMacroTemplateParam {
        name,
        default_value,
        required,
    })
}

fn parse_template_binding(raw: &str) -> Result<(String, String), String> {
    let trimmed = raw.trim();
    let (key, value) = trimmed
        .split_once('=')
        .ok_or_else(|| format!("Invalid --bind '{}': expected KEY=VALUE", raw))?;
    let key = key.trim();
    if key.is_empty() {
        return Err(format!(
            "Invalid --bind '{}': binding key cannot be empty",
            raw
        ));
    }
    Ok((key.to_string(), value.to_string()))
}

fn load_cloning_pattern_file(path: &str) -> Result<CloningPatternFile, String> {
    let raw = fs::read_to_string(path)
        .map_err(|e| format!("Could not read pattern file '{path}': {e}"))?;
    let file: CloningPatternFile = serde_json::from_str(&raw)
        .map_err(|e| format!("Could not parse pattern file '{path}' as JSON: {e}"))?;
    let schema = file.schema.trim();
    if !schema.eq_ignore_ascii_case(CLONING_PATTERN_FILE_SCHEMA) {
        return Err(format!(
            "Unsupported pattern file schema '{}' in '{}' (expected '{}')",
            file.schema, path, CLONING_PATTERN_FILE_SCHEMA
        ));
    }
    if file.templates.is_empty() {
        return Err(format!(
            "Pattern file '{}' does not contain any templates",
            path
        ));
    }
    for template in &file.templates {
        if template.name.trim().is_empty() {
            return Err(format!(
                "Pattern file '{}' contains a template with empty name",
                path
            ));
        }
        if template.script.trim().is_empty() {
            return Err(format!(
                "Pattern file '{}' template '{}' has empty script",
                path, template.name
            ));
        }
    }
    Ok(file)
}

fn load_cloning_pattern_template_file(path: &str) -> Result<CloningPatternTemplateFile, String> {
    let raw = fs::read_to_string(path)
        .map_err(|e| format!("Could not read pattern template file '{path}': {e}"))?;
    let file: CloningPatternTemplateFile = serde_json::from_str(&raw)
        .map_err(|e| format!("Could not parse pattern template file '{path}' as JSON: {e}"))?;
    let schema = file.schema.trim();
    if !schema.eq_ignore_ascii_case(CLONING_PATTERN_TEMPLATE_FILE_SCHEMA) {
        return Err(format!(
            "Unsupported pattern template file schema '{}' in '{}' (expected '{}')",
            file.schema, path, CLONING_PATTERN_TEMPLATE_FILE_SCHEMA
        ));
    }
    if file.name.trim().is_empty() {
        return Err(format!(
            "Pattern template file '{}' contains empty template name",
            path
        ));
    }
    if file.script.trim().is_empty() {
        return Err(format!(
            "Pattern template file '{}' template '{}' has empty script",
            path, file.name
        ));
    }
    Ok(file)
}

fn collect_json_files_recursive(path: &Path) -> Result<Vec<PathBuf>, String> {
    let mut out = vec![];
    if path.is_file() {
        if path
            .extension()
            .and_then(|ext| ext.to_str())
            .is_some_and(|ext| ext.eq_ignore_ascii_case("json"))
        {
            out.push(path.to_path_buf());
        }
        return Ok(out);
    }
    if !path.is_dir() {
        return Err(format!(
            "Pattern import path '{}' is neither a JSON file nor a directory",
            path.display()
        ));
    }

    let mut entries = fs::read_dir(path)
        .map_err(|e| format!("Could not read pattern directory '{}': {e}", path.display()))?
        .filter_map(|entry| entry.ok())
        .collect::<Vec<_>>();
    entries.sort_by_key(|entry| entry.file_name().to_string_lossy().to_string());
    for entry in entries {
        let child = entry.path();
        if child.is_dir() {
            out.extend(collect_json_files_recursive(&child)?);
            continue;
        }
        if child
            .extension()
            .and_then(|ext| ext.to_str())
            .is_some_and(|ext| ext.eq_ignore_ascii_case("json"))
        {
            out.push(child);
        }
    }
    Ok(out)
}

fn load_cloning_pattern_templates_from_path(
    path: &str,
) -> Result<(Vec<CloningPatternTemplate>, Vec<String>), String> {
    let target_path = Path::new(path);
    let json_files = collect_json_files_recursive(target_path)?;
    if json_files.is_empty() {
        return Err(format!(
            "Pattern import path '{}' does not contain any JSON files",
            path
        ));
    }

    let mut templates: Vec<CloningPatternTemplate> = vec![];
    let mut source_files: Vec<String> = vec![];
    for file_path in json_files {
        let file_label = file_path.display().to_string();
        let file_str = file_path.to_string_lossy().to_string();
        let raw = fs::read_to_string(&file_path)
            .map_err(|e| format!("Could not read pattern file '{}': {e}", file_label))?;
        let value: Value = serde_json::from_str(&raw)
            .map_err(|e| format!("Could not parse pattern file '{}': {e}", file_label))?;
        let schema = value
            .get("schema")
            .and_then(|v| v.as_str())
            .unwrap_or_default()
            .trim()
            .to_string();
        if schema.eq_ignore_ascii_case(CLONING_PATTERN_FILE_SCHEMA) {
            let file = load_cloning_pattern_file(&file_str)?;
            templates.extend(file.templates);
            source_files.push(file_label);
        } else if schema.eq_ignore_ascii_case(CLONING_PATTERN_TEMPLATE_FILE_SCHEMA) {
            let file = load_cloning_pattern_template_file(&file_str)?;
            templates.push(CloningPatternTemplate {
                name: file.name,
                description: file.description,
                details_url: file.details_url,
                parameters: file.parameters,
                script: file.script,
            });
            source_files.push(file_label);
        } else {
            return Err(format!(
                "Unsupported pattern schema '{}' in '{}' (expected '{}' or '{}')",
                schema,
                file_label,
                CLONING_PATTERN_FILE_SCHEMA,
                CLONING_PATTERN_TEMPLATE_FILE_SCHEMA
            ));
        }
    }
    if templates.is_empty() {
        return Err(format!(
            "Pattern import path '{}' did not yield any templates",
            path
        ));
    }
    Ok((templates, source_files))
}

fn load_cloning_routine_catalog(path: &str) -> Result<CloningRoutineCatalog, String> {
    let raw = fs::read_to_string(path)
        .map_err(|e| format!("Could not read cloning routine catalog '{}': {e}", path))?;
    let mut catalog: CloningRoutineCatalog = serde_json::from_str(&raw)
        .map_err(|e| format!("Could not parse cloning routine catalog '{}': {e}", path))?;
    let schema = catalog.schema.trim();
    if !schema.eq_ignore_ascii_case(CLONING_ROUTINE_CATALOG_SCHEMA) {
        return Err(format!(
            "Unsupported cloning routine catalog schema '{}' in '{}' (expected '{}')",
            catalog.schema, path, CLONING_ROUTINE_CATALOG_SCHEMA
        ));
    }

    let mut seen_ids = BTreeSet::new();
    for routine in &mut catalog.routines {
        if routine.routine_id.trim().is_empty() {
            return Err(format!(
                "Cloning routine catalog '{}' contains a routine with empty routine_id",
                path
            ));
        }
        if !seen_ids.insert(routine.routine_id.trim().to_string()) {
            return Err(format!(
                "Cloning routine catalog '{}' contains duplicate routine_id '{}'",
                path, routine.routine_id
            ));
        }
        if routine.title.trim().is_empty() {
            return Err(format!(
                "Cloning routine '{}' in '{}' has empty title",
                routine.routine_id, path
            ));
        }
        if routine.family.trim().is_empty() {
            return Err(format!(
                "Cloning routine '{}' in '{}' has empty family",
                routine.routine_id, path
            ));
        }
        if routine.status.trim().is_empty() {
            return Err(format!(
                "Cloning routine '{}' in '{}' has empty status",
                routine.routine_id, path
            ));
        }
        if routine.template_name.trim().is_empty() {
            return Err(format!(
                "Cloning routine '{}' in '{}' has empty template_name",
                routine.routine_id, path
            ));
        }
        routine.vocabulary_tags = routine
            .vocabulary_tags
            .iter()
            .map(|tag| tag.trim().to_string())
            .filter(|tag| !tag.is_empty())
            .collect::<Vec<_>>();
        routine
            .vocabulary_tags
            .sort_by_key(|tag| tag.to_ascii_lowercase());
        routine
            .vocabulary_tags
            .dedup_by(|a, b| a.eq_ignore_ascii_case(b));

        for port in routine
            .input_ports
            .iter()
            .chain(routine.output_ports.iter())
        {
            if port.port_id.trim().is_empty() {
                return Err(format!(
                    "Cloning routine '{}' in '{}' has a port with empty port_id",
                    routine.routine_id, path
                ));
            }
            if port.kind.trim().is_empty() {
                return Err(format!(
                    "Cloning routine '{}' in '{}' has port '{}' with empty kind",
                    routine.routine_id, path, port.port_id
                ));
            }
        }
    }
    Ok(catalog)
}

fn routine_matches_filter(
    routine: &CloningRoutineDefinition,
    family: Option<&str>,
    status: Option<&str>,
    tag: Option<&str>,
    query: Option<&str>,
) -> bool {
    if let Some(family_filter) = family.map(str::trim).filter(|v| !v.is_empty()) {
        if !routine.family.eq_ignore_ascii_case(family_filter) {
            return false;
        }
    }
    if let Some(status_filter) = status.map(str::trim).filter(|v| !v.is_empty()) {
        if !routine.status.eq_ignore_ascii_case(status_filter) {
            return false;
        }
    }
    if let Some(tag_filter) = tag.map(str::trim).filter(|v| !v.is_empty()) {
        if !routine
            .vocabulary_tags
            .iter()
            .any(|entry| entry.eq_ignore_ascii_case(tag_filter))
        {
            return false;
        }
    }
    if let Some(query_filter) = query
        .map(str::trim)
        .filter(|v| !v.is_empty())
        .map(|v| v.to_ascii_lowercase())
    {
        let mut haystack = vec![
            routine.routine_id.to_ascii_lowercase(),
            routine.title.to_ascii_lowercase(),
            routine.family.to_ascii_lowercase(),
            routine.status.to_ascii_lowercase(),
            routine.template_name.to_ascii_lowercase(),
            routine
                .summary
                .as_deref()
                .unwrap_or_default()
                .to_ascii_lowercase(),
        ];
        haystack.extend(
            routine
                .vocabulary_tags
                .iter()
                .map(|tag| tag.to_ascii_lowercase()),
        );
        if !haystack.iter().any(|entry| entry.contains(&query_filter)) {
            return false;
        }
    }
    true
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
            Self::InspectFeatureExpert { seq_id, target } => {
                format!(
                    "inspect feature expert view for '{seq_id}' target={}",
                    target.describe()
                )
            }
            Self::RenderFeatureExpertSvg {
                seq_id,
                target,
                output,
            } => format!(
                "render feature expert SVG for '{seq_id}' target={} to '{output}'",
                target.describe()
            ),
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
                container_ids,
                arrangement_id,
            } => {
                let ladders = ladders
                    .as_ref()
                    .map(|v| v.join(","))
                    .unwrap_or_else(|| "auto".to_string());
                let containers = container_ids
                    .as_ref()
                    .map(|v| v.join(","))
                    .unwrap_or_else(|| "-".to_string());
                let arrangement = arrangement_id
                    .as_deref()
                    .map(str::trim)
                    .filter(|v| !v.is_empty())
                    .unwrap_or("-");
                format!(
                    "render serial gel SVG to '{output}' (inputs={}, containers={}, arrangement={}, ladders={ladders})",
                    inputs.len(),
                    containers,
                    arrangement
                )
            }
            Self::CreateArrangementSerial {
                container_ids,
                arrangement_id,
                name,
                ladders,
            } => {
                let arrangement_id = arrangement_id.as_deref().unwrap_or("auto");
                let name = name.as_deref().unwrap_or("-");
                let ladders = ladders
                    .as_ref()
                    .map(|v| v.join(","))
                    .unwrap_or_else(|| "auto".to_string());
                format!(
                    "create serial arrangement id={arrangement_id}, name={name}, lanes={} (ladders={ladders})",
                    container_ids.len()
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
            Self::RoutinesList {
                catalog_path,
                family,
                status,
                tag,
                query,
            } => {
                let catalog = catalog_path
                    .clone()
                    .unwrap_or_else(|| DEFAULT_CLONING_ROUTINE_CATALOG_PATH.to_string());
                let family = family
                    .as_deref()
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .unwrap_or("-");
                let status = status
                    .as_deref()
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .unwrap_or("-");
                let tag = tag
                    .as_deref()
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .unwrap_or("-");
                let query = query
                    .as_deref()
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .unwrap_or("-");
                format!(
                    "list cloning routines from '{catalog}' (family={family}, status={status}, tag={tag}, query={query})"
                )
            }
            Self::AgentsList { catalog_path } => {
                let catalog = catalog_path
                    .clone()
                    .unwrap_or_else(|| "assets/agent_systems.json".to_string());
                format!("list agent systems from catalog '{catalog}'")
            }
            Self::AgentsAsk {
                system_id,
                prompt,
                catalog_path,
                base_url_override,
                model_override,
                timeout_seconds,
                connect_timeout_seconds,
                read_timeout_seconds,
                max_retries,
                max_response_bytes,
                include_state_summary,
                allow_auto_exec,
                execute_all,
                execute_indices,
            } => {
                let catalog = catalog_path
                    .clone()
                    .unwrap_or_else(|| "assets/agent_systems.json".to_string());
                let execute_mode = if *execute_all {
                    "all".to_string()
                } else if !execute_indices.is_empty() {
                    format!("indices={}", execute_indices.len())
                } else if *allow_auto_exec {
                    "auto-only".to_string()
                } else {
                    "none".to_string()
                };
                let base_url = base_url_override
                    .as_deref()
                    .map(str::trim)
                    .filter(|v| !v.is_empty())
                    .unwrap_or("-");
                let model = model_override
                    .as_deref()
                    .map(str::trim)
                    .filter(|v| !v.is_empty())
                    .unwrap_or("-");
                let timeout = timeout_seconds
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "-".to_string());
                let connect_timeout = connect_timeout_seconds
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "-".to_string());
                let read_timeout = read_timeout_seconds
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "-".to_string());
                let retries = max_retries
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "-".to_string());
                let max_bytes = max_response_bytes
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "-".to_string());
                format!(
                    "ask agent '{system_id}' (catalog='{catalog}', prompt_len={}, include_state_summary={}, base_url_override={base_url}, model_override={model}, timeout_secs={timeout}, connect_timeout_secs={connect_timeout}, read_timeout_secs={read_timeout}, max_retries={retries}, max_response_bytes={max_bytes}, execute={execute_mode})",
                    prompt.len(),
                    include_state_summary
                )
            }
            Self::UiListIntents => "list supported GUI intent commands".to_string(),
            Self::UiIntent {
                action,
                target,
                genome_id,
                helper_mode,
                catalog_path,
                cache_dir,
                filter,
                species,
                latest,
            } => {
                let genome = genome_id.as_deref().unwrap_or("-");
                let scope = if *helper_mode { "helpers" } else { "genomes" };
                let catalog = catalog_path.as_deref().unwrap_or("-");
                let cache = cache_dir.as_deref().unwrap_or("-");
                let filter = filter.as_deref().unwrap_or("-");
                let species = species.as_deref().unwrap_or("-");
                format!(
                    "request GUI {} for '{}' (genome_id={genome}, scope={scope}, catalog='{catalog}', cache='{cache}', filter='{filter}', species='{species}', latest={latest})",
                    action.as_str(),
                    target.as_str()
                )
            }
            Self::UiPreparedGenomes {
                helper_mode,
                catalog_path,
                cache_dir,
                filter,
                species,
                latest,
            } => {
                let scope = if *helper_mode { "helpers" } else { "genomes" };
                let catalog = catalog_path
                    .clone()
                    .unwrap_or_else(|| default_catalog_path(*helper_mode).to_string());
                let cache = cache_dir.clone().unwrap_or_else(|| "-".to_string());
                let filter = filter.as_deref().unwrap_or("-");
                let species = species.as_deref().unwrap_or("-");
                format!(
                    "list prepared {scope} (catalog='{catalog}', cache='{cache}', filter='{filter}', species='{species}', latest={latest})"
                )
            }
            Self::UiLatestPrepared {
                helper_mode,
                catalog_path,
                cache_dir,
                species,
            } => {
                let scope = if *helper_mode { "helpers" } else { "genomes" };
                let catalog = catalog_path
                    .clone()
                    .unwrap_or_else(|| default_catalog_path(*helper_mode).to_string());
                let cache = cache_dir.clone().unwrap_or_else(|| "-".to_string());
                format!(
                    "select latest prepared {scope} for species '{species}' (catalog='{catalog}', cache='{cache}')"
                )
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
                    track_name
                        .clone()
                        .unwrap_or_else(|| "blast_hits".to_string()),
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
                track_name.clone().unwrap_or_else(|| "-".to_string()),
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
                track_name.clone().unwrap_or_else(|| "-".to_string()),
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
                track_name.clone().unwrap_or_else(|| "-".to_string()),
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
            Self::MacrosRun {
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
                format!("run {mode} workflow macro '{preview}'")
            }
            Self::MacrosTemplateList => "list workflow macro templates".to_string(),
            Self::MacrosTemplateShow { name } => {
                format!("show workflow macro template '{}'", name)
            }
            Self::MacrosTemplateUpsert {
                name,
                description,
                details_url,
                parameters,
                script,
            } => format!(
                "upsert workflow macro template '{}' (description='{}', details_url='{}', params={}, script_len={})",
                name,
                description
                    .as_deref()
                    .filter(|value| !value.trim().is_empty())
                    .unwrap_or("-"),
                details_url
                    .as_deref()
                    .filter(|value| !value.trim().is_empty())
                    .unwrap_or("-"),
                parameters
                    .iter()
                    .map(|p| {
                        if let Some(default_value) = &p.default_value {
                            format!("{}={}", p.name, default_value)
                        } else {
                            p.name.clone()
                        }
                    })
                    .collect::<Vec<_>>()
                    .join(","),
                script.len()
            ),
            Self::MacrosTemplateDelete { name } => {
                format!("delete workflow macro template '{}'", name)
            }
            Self::MacrosTemplateImport { path } => {
                format!("import workflow macro templates from '{}'", path)
            }
            Self::MacrosTemplateRun {
                name,
                bindings,
                transactional,
            } => format!(
                "run workflow macro template '{}' (bindings={}, transactional={})",
                name,
                bindings
                    .iter()
                    .map(|(key, value)| format!("{key}={value}"))
                    .collect::<Vec<_>>()
                    .join(","),
                transactional
            ),
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
                min.map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                max.map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string()),
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
            Self::CandidatesScoreWeightedObjective {
                set_name,
                metric,
                objectives,
                normalize_metrics,
            } => format!(
                "compute weighted objective metric '{}' for candidate set '{}' (terms={}, normalize_metrics={})",
                metric,
                set_name,
                objectives
                    .iter()
                    .map(|term| format!(
                        "{}:{}:{}",
                        term.metric,
                        term.weight,
                        term.direction.as_str()
                    ))
                    .collect::<Vec<_>>()
                    .join(","),
                normalize_metrics
            ),
            Self::CandidatesTopK {
                input_set,
                output_set,
                metric,
                k,
                direction,
                tie_break,
            } => format!(
                "select top-k candidate set from '{}' into '{}' by metric '{}' (k={}, direction={}, tie_break={})",
                input_set,
                output_set,
                metric,
                k,
                direction.as_str(),
                tie_break.as_str()
            ),
            Self::CandidatesParetoFrontier {
                input_set,
                output_set,
                objectives,
                max_candidates,
                tie_break,
            } => format!(
                "compute pareto frontier from '{}' into '{}' (objectives={}, max_candidates={}, tie_break={})",
                input_set,
                output_set,
                objectives
                    .iter()
                    .map(|objective| format!(
                        "{}:{}",
                        objective.metric,
                        objective.direction.as_str()
                    ))
                    .collect::<Vec<_>>()
                    .join(","),
                max_candidates
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                tie_break.as_str()
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
            Self::CandidatesTemplateList => "list candidate macro templates".to_string(),
            Self::CandidatesTemplateShow { name } => {
                format!("show candidate macro template '{}'", name)
            }
            Self::CandidatesTemplateUpsert {
                name,
                description,
                details_url,
                parameters,
                script,
            } => format!(
                "upsert candidate macro template '{}' (description='{}', details_url='{}', params={}, script_len={})",
                name,
                description
                    .as_deref()
                    .filter(|value| !value.trim().is_empty())
                    .unwrap_or("-"),
                details_url
                    .as_deref()
                    .filter(|value| !value.trim().is_empty())
                    .unwrap_or("-"),
                parameters
                    .iter()
                    .map(|p| {
                        if let Some(default_value) = &p.default_value {
                            format!("{}={}", p.name, default_value)
                        } else {
                            p.name.clone()
                        }
                    })
                    .collect::<Vec<_>>()
                    .join(","),
                script.len()
            ),
            Self::CandidatesTemplateDelete { name } => {
                format!("delete candidate macro template '{}'", name)
            }
            Self::CandidatesTemplateRun {
                name,
                bindings,
                transactional,
            } => format!(
                "run candidate macro template '{}' (bindings={}, transactional={})",
                name,
                bindings
                    .iter()
                    .map(|(key, value)| format!("{key}={value}"))
                    .collect::<Vec<_>>()
                    .join(","),
                transactional
            ),
            Self::GuidesList => "list persisted guide sets".to_string(),
            Self::GuidesShow {
                guide_set_id,
                limit,
                offset,
            } => format!(
                "show guide set '{}' (limit={}, offset={})",
                guide_set_id, limit, offset
            ),
            Self::GuidesPut {
                guide_set_id,
                guides_json,
            } => format!(
                "upsert guide set '{}' from JSON payload (len={})",
                guide_set_id,
                guides_json.len()
            ),
            Self::GuidesDelete { guide_set_id } => {
                format!("delete guide set '{}'", guide_set_id)
            }
            Self::GuidesFilter {
                guide_set_id,
                config_json,
                output_guide_set_id,
            } => format!(
                "filter practical guide constraints for '{}' (config_len={}, output_set='{}')",
                guide_set_id,
                config_json.as_ref().map(|v| v.len()).unwrap_or(2),
                output_guide_set_id
                    .as_deref()
                    .filter(|v| !v.trim().is_empty())
                    .unwrap_or("-")
            ),
            Self::GuidesFilterShow { guide_set_id } => {
                format!("show practical guide filter report for '{}'", guide_set_id)
            }
            Self::GuidesOligosGenerate {
                guide_set_id,
                template_id,
                apply_5prime_g_extension,
                output_oligo_set_id,
                passed_only,
            } => format!(
                "generate guide oligos for '{}' (template='{}', apply_5prime_g_extension={}, output_oligo_set='{}', passed_only={})",
                guide_set_id,
                template_id,
                apply_5prime_g_extension,
                output_oligo_set_id
                    .as_deref()
                    .filter(|v| !v.trim().is_empty())
                    .unwrap_or("-"),
                passed_only
            ),
            Self::GuidesOligosList { guide_set_id } => format!(
                "list guide oligo sets{}",
                guide_set_id
                    .as_deref()
                    .map(|id| format!(" for '{}'", id))
                    .unwrap_or_default()
            ),
            Self::GuidesOligosShow { oligo_set_id } => {
                format!("show guide oligo set '{}'", oligo_set_id)
            }
            Self::GuidesOligosExport {
                guide_set_id,
                oligo_set_id,
                format,
                path,
                plate_format,
            } => format!(
                "export guide oligos for '{}' to '{}' (oligo_set='{}', format={}, plate={})",
                guide_set_id,
                path,
                oligo_set_id
                    .as_deref()
                    .filter(|v| !v.trim().is_empty())
                    .unwrap_or("-"),
                format.as_str(),
                plate_format
                    .map(|v| v.as_str().to_string())
                    .unwrap_or_else(|| "-".to_string())
            ),
            Self::GuidesProtocolExport {
                guide_set_id,
                oligo_set_id,
                path,
                include_qc_checklist,
            } => format!(
                "export guide protocol text for '{}' to '{}' (oligo_set='{}', include_qc_checklist={})",
                guide_set_id,
                path,
                oligo_set_id
                    .as_deref()
                    .filter(|v| !v.trim().is_empty())
                    .unwrap_or("-"),
                include_qc_checklist
            ),
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
                | Self::MacrosRun { .. }
                | Self::MacrosTemplateUpsert { .. }
                | Self::MacrosTemplateDelete { .. }
                | Self::MacrosTemplateImport { .. }
                | Self::MacrosTemplateRun { .. }
                | Self::CandidatesDelete { .. }
                | Self::CandidatesGenerate { .. }
                | Self::CandidatesGenerateBetweenAnchors { .. }
                | Self::CandidatesScoreExpression { .. }
                | Self::CandidatesScoreDistance { .. }
                | Self::CandidatesFilter { .. }
                | Self::CandidatesSetOp { .. }
                | Self::CandidatesScoreWeightedObjective { .. }
                | Self::CandidatesTopK { .. }
                | Self::CandidatesParetoFrontier { .. }
                | Self::CandidatesMacro { .. }
                | Self::CandidatesTemplateUpsert { .. }
                | Self::CandidatesTemplateDelete { .. }
                | Self::CandidatesTemplateRun { .. }
                | Self::GuidesPut { .. }
                | Self::GuidesDelete { .. }
                | Self::GuidesFilter { .. }
                | Self::GuidesOligosGenerate { .. }
                | Self::GuidesOligosExport { .. }
                | Self::GuidesProtocolExport { .. }
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

fn parse_feature_expert_target_tokens(
    tokens: &[String],
    context: &str,
) -> Result<FeatureExpertTarget, String> {
    if tokens.len() < 2 {
        return Err(format!(
            "{context} requires target syntax: tfbs FEATURE_ID | restriction CUT_POS_1BASED [--enzyme NAME] [--start START_1BASED] [--end END_1BASED] | splicing FEATURE_ID"
        ));
    }
    match tokens[0].trim().to_ascii_lowercase().as_str() {
        "tfbs" => {
            if tokens.len() != 2 {
                return Err(format!(
                    "{context} tfbs target expects exactly: tfbs FEATURE_ID"
                ));
            }
            let feature_id = tokens[1]
                .parse::<usize>()
                .map_err(|e| format!("Invalid TFBS feature id '{}': {e}", tokens[1]))?;
            Ok(FeatureExpertTarget::TfbsFeature { feature_id })
        }
        "restriction" | "re" => {
            let cut_pos_1based = tokens[1]
                .parse::<usize>()
                .map_err(|e| format!("Invalid restriction cut position '{}': {e}", tokens[1]))?;
            if cut_pos_1based == 0 {
                return Err("Restriction cut position must be >= 1".to_string());
            }
            let mut enzyme: Option<String> = None;
            let mut recognition_start_1based: Option<usize> = None;
            let mut recognition_end_1based: Option<usize> = None;
            let mut idx = 2usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--enzyme" => {
                        let raw = parse_option_path(tokens, &mut idx, "--enzyme", context)?;
                        let trimmed = raw.trim();
                        if trimmed.is_empty() {
                            return Err("--enzyme value must not be empty".to_string());
                        }
                        enzyme = Some(trimmed.to_string());
                    }
                    "--start" => {
                        let raw = parse_option_path(tokens, &mut idx, "--start", context)?;
                        let parsed = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --start value '{raw}': {e}"))?;
                        if parsed == 0 {
                            return Err("--start must be >= 1".to_string());
                        }
                        recognition_start_1based = Some(parsed);
                    }
                    "--end" => {
                        let raw = parse_option_path(tokens, &mut idx, "--end", context)?;
                        let parsed = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --end value '{raw}': {e}"))?;
                        if parsed == 0 {
                            return Err("--end must be >= 1".to_string());
                        }
                        recognition_end_1based = Some(parsed);
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for {context} restriction"
                        ));
                    }
                }
            }
            if let (Some(start), Some(end)) = (recognition_start_1based, recognition_end_1based) {
                if start > end {
                    return Err("--start must be <= --end".to_string());
                }
            }
            Ok(FeatureExpertTarget::RestrictionSite {
                cut_pos_1based,
                enzyme,
                recognition_start_1based,
                recognition_end_1based,
            })
        }
        "splicing" => {
            if tokens.len() != 2 {
                return Err(format!(
                    "{context} splicing target expects exactly: splicing FEATURE_ID"
                ));
            }
            let feature_id = tokens[1]
                .parse::<usize>()
                .map_err(|e| format!("Invalid splicing feature id '{}': {e}", tokens[1]))?;
            Ok(FeatureExpertTarget::SplicingFeature { feature_id })
        }
        other => Err(format!(
            "Unknown feature target '{other}' (expected tfbs|restriction|splicing)"
        )),
    }
}

fn parse_ladder_molecule(value: &str) -> Result<LadderMolecule, String> {
    LadderMolecule::parse(value)
        .ok_or_else(|| format!("Unknown ladder molecule '{value}', expected 'dna' or 'rna'"))
}

fn parse_guide_export_format(value: &str) -> Result<GuideOligoExportFormat, String> {
    match value.trim().to_ascii_lowercase().as_str() {
        "csv_table" | "csv" | "table" => Ok(GuideOligoExportFormat::CsvTable),
        "plate_csv" | "plate" => Ok(GuideOligoExportFormat::PlateCsv),
        "fasta" | "fa" => Ok(GuideOligoExportFormat::Fasta),
        other => Err(format!(
            "Unsupported guide export format '{other}' (expected csv_table|plate_csv|fasta)"
        )),
    }
}

fn parse_guide_plate_format(value: &str) -> Result<GuideOligoPlateFormat, String> {
    match value.trim().to_ascii_lowercase().as_str() {
        "96" | "plate96" | "p96" => Ok(GuideOligoPlateFormat::Plate96),
        "384" | "plate384" | "p384" => Ok(GuideOligoPlateFormat::Plate384),
        other => Err(format!(
            "Unsupported plate format '{other}' (expected 96|384)"
        )),
    }
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
                        catalog_path =
                            Some(parse_option_path(tokens, &mut idx, "--catalog", label)?)
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
                        catalog_path =
                            Some(parse_option_path(tokens, &mut idx, "--catalog", label)?)
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
                        cache_dir = Some(parse_option_path(tokens, &mut idx, "--cache-dir", label)?)
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
                        cache_dir = Some(parse_option_path(tokens, &mut idx, "--cache-dir", label)?)
                    }
                    "--filter" => filter = parse_option_path(tokens, &mut idx, "--filter", label)?,
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
                        cache_dir = Some(parse_option_path(tokens, &mut idx, "--cache-dir", label)?)
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
                                ));
                            }
                        }
                    }
                    "--catalog" => {
                        catalog_path =
                            Some(parse_option_path(tokens, &mut idx, "--catalog", label)?)
                    }
                    "--cache-dir" => {
                        cache_dir = Some(parse_option_path(tokens, &mut idx, "--cache-dir", label)?)
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
                                ));
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
                        cache_dir = Some(parse_option_path(tokens, &mut idx, "--cache-dir", label)?)
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
                        output_id = Some(parse_option_path(tokens, &mut idx, "--output-id", label)?)
                    }
                    "--catalog" => {
                        catalog_path =
                            Some(parse_option_path(tokens, &mut idx, "--catalog", label)?)
                    }
                    "--cache-dir" => {
                        cache_dir = Some(parse_option_path(tokens, &mut idx, "--cache-dir", label)?)
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
                        output_id = Some(parse_option_path(tokens, &mut idx, "--output-id", label)?)
                    }
                    "--catalog" => {
                        catalog_path =
                            Some(parse_option_path(tokens, &mut idx, "--catalog", label)?)
                    }
                    "--cache-dir" => {
                        cache_dir = Some(parse_option_path(tokens, &mut idx, "--cache-dir", label)?)
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
                        cache_dir = Some(parse_option_path(tokens, &mut idx, "--cache-dir", label)?)
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
            "candidates requires a subcommand: list, delete, generate, generate-between-anchors, show, metrics, score, score-distance, score-weighted, top-k, pareto, filter, set-op, macro, template-list, template-show, template-put, template-delete, template-run"
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
                        max_distance_bp =
                            Some(raw.parse::<usize>().map_err(|e| {
                                format!("Invalid --max-distance value '{raw}': {e}")
                            })?);
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
                        return Err(format!("Unknown option '{other}' for candidates generate"));
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
                return Err(
                    "candidates show requires SET_NAME [--limit N] [--offset N]".to_string()
                );
            }
            let set_name = tokens[2].clone();
            let mut limit = DEFAULT_CANDIDATE_PAGE_SIZE;
            let mut offset = 0usize;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--limit" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--limit", "candidates show")?;
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
        "score-weighted" => {
            if tokens.len() < 5 {
                return Err(
                    "candidates score-weighted requires SET_NAME METRIC_NAME --term METRIC:WEIGHT[:max|min] [--term ...] [--normalize|--no-normalize]"
                        .to_string(),
                );
            }
            let set_name = tokens[2].clone();
            let metric = tokens[3].clone();
            let mut objectives: Vec<CandidateWeightedObjectiveTerm> = vec![];
            let mut normalize_metrics = true;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--term" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--term",
                            "candidates score-weighted",
                        )?;
                        objectives.push(parse_weighted_objective_term(&raw)?);
                    }
                    "--normalize" => {
                        normalize_metrics = true;
                        idx += 1;
                    }
                    "--no-normalize" => {
                        normalize_metrics = false;
                        idx += 1;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for candidates score-weighted"
                        ));
                    }
                }
            }
            if objectives.is_empty() {
                return Err(
                    "candidates score-weighted requires at least one --term METRIC:WEIGHT[:max|min]"
                        .to_string(),
                );
            }
            Ok(ShellCommand::CandidatesScoreWeightedObjective {
                set_name,
                metric,
                objectives,
                normalize_metrics,
            })
        }
        "top-k" => {
            if tokens.len() < 4 {
                return Err(
                    "candidates top-k requires INPUT_SET OUTPUT_SET --metric METRIC_NAME --k N [--direction max|min] [--tie-break POLICY]"
                        .to_string(),
                );
            }
            let input_set = tokens[2].clone();
            let output_set = tokens[3].clone();
            let mut metric: Option<String> = None;
            let mut k: Option<usize> = None;
            let mut direction = CandidateObjectiveDirection::Maximize;
            let mut tie_break = CandidateTieBreakPolicy::SeqStartEnd;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--metric" => {
                        metric = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--metric",
                            "candidates top-k",
                        )?);
                    }
                    "--k" => {
                        let raw = parse_option_path(tokens, &mut idx, "--k", "candidates top-k")?;
                        let parsed = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --k value '{raw}': {e}"))?;
                        if parsed == 0 {
                            return Err("--k must be >= 1".to_string());
                        }
                        k = Some(parsed);
                    }
                    "--direction" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--direction", "candidates top-k")?;
                        direction = parse_candidate_objective_direction(&raw)?;
                    }
                    "--tie-break" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--tie-break", "candidates top-k")?;
                        tie_break = parse_candidate_tie_break_policy(&raw)?;
                    }
                    other => return Err(format!("Unknown option '{other}' for candidates top-k")),
                }
            }
            let metric = metric.ok_or_else(|| "candidates top-k requires --metric".to_string())?;
            let k = k.ok_or_else(|| "candidates top-k requires --k N".to_string())?;
            Ok(ShellCommand::CandidatesTopK {
                input_set,
                output_set,
                metric,
                k,
                direction,
                tie_break,
            })
        }
        "pareto" => {
            if tokens.len() < 4 {
                return Err(
                    "candidates pareto requires INPUT_SET OUTPUT_SET --objective METRIC[:max|min] [--objective ...] [--max-candidates N] [--tie-break POLICY]"
                        .to_string(),
                );
            }
            let input_set = tokens[2].clone();
            let output_set = tokens[3].clone();
            let mut objectives: Vec<CandidateObjectiveSpec> = vec![];
            let mut max_candidates: Option<usize> = None;
            let mut tie_break = CandidateTieBreakPolicy::SeqStartEnd;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--objective" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--objective",
                            "candidates pareto",
                        )?;
                        objectives.push(parse_pareto_objective(&raw)?);
                    }
                    "--max-candidates" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--max-candidates",
                            "candidates pareto",
                        )?;
                        let parsed = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --max-candidates value '{raw}': {e}"))?;
                        if parsed == 0 {
                            return Err("--max-candidates must be >= 1".to_string());
                        }
                        max_candidates = Some(parsed);
                    }
                    "--tie-break" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--tie-break",
                            "candidates pareto",
                        )?;
                        tie_break = parse_candidate_tie_break_policy(&raw)?;
                    }
                    other => return Err(format!("Unknown option '{other}' for candidates pareto")),
                }
            }
            if objectives.is_empty() {
                return Err(
                    "candidates pareto requires at least one --objective METRIC[:max|min]"
                        .to_string(),
                );
            }
            Ok(ShellCommand::CandidatesParetoFrontier {
                input_set,
                output_set,
                objectives,
                max_candidates,
                tie_break,
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
                        min_quantile =
                            Some(raw.parse::<f64>().map_err(|e| {
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
                        max_quantile =
                            Some(raw.parse::<f64>().map_err(|e| {
                                format!("Invalid --max-quantile value '{raw}': {e}")
                            })?);
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for candidates filter"));
                    }
                }
            }
            let metric = metric.ok_or_else(|| "candidates filter requires --metric".to_string())?;
            if min.zip(max).map(|(lo, hi)| lo > hi).unwrap_or(false) {
                return Err("--min must be <= --max".to_string());
            }
            if min_quantile
                .zip(max_quantile)
                .map(|(lo, hi)| lo > hi)
                .unwrap_or(false)
            {
                return Err("--min-quantile must be <= --max-quantile".to_string());
            }
            for (name, value) in [
                ("min-quantile", min_quantile),
                ("max-quantile", max_quantile),
            ] {
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
                                "candidates macro --file may only be specified once".to_string()
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
        "template-list" => {
            if tokens.len() != 2 {
                return Err("candidates template-list takes no options".to_string());
            }
            Ok(ShellCommand::CandidatesTemplateList)
        }
        "template-show" => {
            if tokens.len() != 3 {
                return Err("candidates template-show requires TEMPLATE_NAME".to_string());
            }
            Ok(ShellCommand::CandidatesTemplateShow {
                name: tokens[2].clone(),
            })
        }
        "template-put" | "template-upsert" => {
            if tokens.len() < 4 {
                return Err(
                    "candidates template-put requires TEMPLATE_NAME (--script SCRIPT_OR_@FILE | --file PATH) [--description TEXT] [--details-url URL] [--param NAME|NAME=DEFAULT ...]"
                        .to_string(),
                );
            }
            let name = tokens[2].clone();
            let mut description: Option<String> = None;
            let mut details_url: Option<String> = None;
            let mut parameters: Vec<CandidateMacroTemplateParam> = vec![];
            let mut script: Option<String> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--description" => {
                        description = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--description",
                            "candidates template-put",
                        )?);
                    }
                    "--details-url" | "--url" => {
                        if details_url.is_some() {
                            return Err(
                                "candidates template-put details URL was already specified"
                                    .to_string(),
                            );
                        }
                        details_url = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--details-url",
                            "candidates template-put",
                        )?);
                    }
                    "--param" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--param",
                            "candidates template-put",
                        )?;
                        parameters.push(parse_candidate_template_param_spec(&raw)?);
                    }
                    "--script" => {
                        if script.is_some() {
                            return Err(
                                "candidates template-put script was already specified".to_string()
                            );
                        }
                        script = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--script",
                            "candidates template-put",
                        )?);
                    }
                    "--file" => {
                        if script.is_some() {
                            return Err(
                                "candidates template-put script was already specified".to_string()
                            );
                        }
                        let path = parse_option_path(
                            tokens,
                            &mut idx,
                            "--file",
                            "candidates template-put",
                        )?;
                        script = Some(format!("@{path}"));
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for candidates template-put"
                        ));
                    }
                }
            }
            let script = script.ok_or_else(|| {
                "candidates template-put requires --script SCRIPT_OR_@FILE or --file PATH"
                    .to_string()
            })?;
            Ok(ShellCommand::CandidatesTemplateUpsert {
                name,
                description,
                details_url,
                parameters,
                script,
            })
        }
        "template-delete" => {
            if tokens.len() != 3 {
                return Err("candidates template-delete requires TEMPLATE_NAME".to_string());
            }
            Ok(ShellCommand::CandidatesTemplateDelete {
                name: tokens[2].clone(),
            })
        }
        "template-run" => {
            if tokens.len() < 3 {
                return Err(
                    "candidates template-run requires TEMPLATE_NAME [--bind KEY=VALUE ...] [--transactional]"
                        .to_string(),
                );
            }
            let name = tokens[2].clone();
            let mut bindings: HashMap<String, String> = HashMap::new();
            let mut transactional = false;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--bind" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--bind",
                            "candidates template-run",
                        )?;
                        let (key, value) = parse_template_binding(&raw)?;
                        if bindings.insert(key.clone(), value).is_some() {
                            return Err(format!(
                                "Duplicate --bind key '{}' in candidates template-run",
                                key
                            ));
                        }
                    }
                    "--transactional" | "--atomic" => {
                        transactional = true;
                        idx += 1;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for candidates template-run"
                        ));
                    }
                }
            }
            Ok(ShellCommand::CandidatesTemplateRun {
                name,
                bindings,
                transactional,
            })
        }
        other => Err(format!(
            "Unknown candidates subcommand '{other}' (expected list, delete, generate, generate-between-anchors, show, metrics, score, score-distance, score-weighted, top-k, pareto, filter, set-op, macro, template-list, template-show, template-put, template-delete, template-run)"
        )),
    }
}

fn parse_guides_command(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.len() < 2 {
        return Err(
            "guides requires a subcommand: list, show, put, delete, filter, filter-show, oligos-generate, oligos-list, oligos-show, oligos-export, protocol-export"
                .to_string(),
        );
    }
    match tokens[1].as_str() {
        "list" => {
            if tokens.len() > 2 {
                return Err("guides list takes no options".to_string());
            }
            Ok(ShellCommand::GuidesList)
        }
        "show" => {
            if tokens.len() < 3 {
                return Err(
                    "guides show requires GUIDE_SET_ID [--limit N] [--offset N]".to_string()
                );
            }
            let guide_set_id = tokens[2].clone();
            let mut limit = DEFAULT_GUIDE_PAGE_SIZE;
            let mut offset = 0usize;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--limit" => {
                        let raw = parse_option_path(tokens, &mut idx, "--limit", "guides show")?;
                        limit = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --limit value '{raw}': {e}"))?;
                        if limit == 0 {
                            return Err("--limit must be >= 1".to_string());
                        }
                    }
                    "--offset" => {
                        let raw = parse_option_path(tokens, &mut idx, "--offset", "guides show")?;
                        offset = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --offset value '{raw}': {e}"))?;
                    }
                    other => return Err(format!("Unknown option '{other}' for guides show")),
                }
            }
            Ok(ShellCommand::GuidesShow {
                guide_set_id,
                limit,
                offset,
            })
        }
        "put" | "upsert" => {
            if tokens.len() < 4 {
                return Err(
                    "guides put requires GUIDE_SET_ID (--json JSON|@FILE | --file PATH)"
                        .to_string(),
                );
            }
            let guide_set_id = tokens[2].clone();
            let mut guides_json: Option<String> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--json" => {
                        if guides_json.is_some() {
                            return Err("guides put JSON payload was already specified".to_string());
                        }
                        guides_json =
                            Some(parse_option_path(tokens, &mut idx, "--json", "guides put")?);
                    }
                    "--file" => {
                        if guides_json.is_some() {
                            return Err("guides put JSON payload was already specified".to_string());
                        }
                        let path = parse_option_path(tokens, &mut idx, "--file", "guides put")?;
                        guides_json = Some(format!("@{path}"));
                    }
                    other => return Err(format!("Unknown option '{other}' for guides put")),
                }
            }
            let guides_json = guides_json.ok_or_else(|| {
                "guides put requires --json JSON|@FILE or --file PATH".to_string()
            })?;
            Ok(ShellCommand::GuidesPut {
                guide_set_id,
                guides_json,
            })
        }
        "delete" => {
            if tokens.len() != 3 {
                return Err("guides delete requires GUIDE_SET_ID".to_string());
            }
            Ok(ShellCommand::GuidesDelete {
                guide_set_id: tokens[2].clone(),
            })
        }
        "filter" => {
            if tokens.len() < 3 {
                return Err(
                    "guides filter requires GUIDE_SET_ID [--config JSON|@FILE] [--config-file PATH] [--output-set GUIDE_SET_ID]"
                        .to_string(),
                );
            }
            let guide_set_id = tokens[2].clone();
            let mut config_json: Option<String> = None;
            let mut output_guide_set_id: Option<String> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--config" => {
                        if config_json.is_some() {
                            return Err("guides filter config was already specified".to_string());
                        }
                        config_json = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--config",
                            "guides filter",
                        )?);
                    }
                    "--config-file" => {
                        if config_json.is_some() {
                            return Err("guides filter config was already specified".to_string());
                        }
                        let path =
                            parse_option_path(tokens, &mut idx, "--config-file", "guides filter")?;
                        config_json = Some(format!("@{path}"));
                    }
                    "--output-set" => {
                        output_guide_set_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--output-set",
                            "guides filter",
                        )?);
                    }
                    other => return Err(format!("Unknown option '{other}' for guides filter")),
                }
            }
            Ok(ShellCommand::GuidesFilter {
                guide_set_id,
                config_json,
                output_guide_set_id,
            })
        }
        "filter-show" => {
            if tokens.len() != 3 {
                return Err("guides filter-show requires GUIDE_SET_ID".to_string());
            }
            Ok(ShellCommand::GuidesFilterShow {
                guide_set_id: tokens[2].clone(),
            })
        }
        "oligos-generate" => {
            if tokens.len() < 4 {
                return Err(
                    "guides oligos-generate requires GUIDE_SET_ID TEMPLATE_ID [--apply-5prime-g-extension] [--output-oligo-set ID] [--passed-only]"
                        .to_string(),
                );
            }
            let guide_set_id = tokens[2].clone();
            let template_id = tokens[3].clone();
            let mut apply_5prime_g_extension = false;
            let mut output_oligo_set_id: Option<String> = None;
            let mut passed_only = false;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--apply-5prime-g-extension" => {
                        apply_5prime_g_extension = true;
                        idx += 1;
                    }
                    "--output-oligo-set" => {
                        output_oligo_set_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--output-oligo-set",
                            "guides oligos-generate",
                        )?);
                    }
                    "--passed-only" => {
                        passed_only = true;
                        idx += 1;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for guides oligos-generate"
                        ));
                    }
                }
            }
            Ok(ShellCommand::GuidesOligosGenerate {
                guide_set_id,
                template_id,
                apply_5prime_g_extension,
                output_oligo_set_id,
                passed_only,
            })
        }
        "oligos-list" => {
            let mut guide_set_id: Option<String> = None;
            let mut idx = 2usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--guide-set" => {
                        guide_set_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--guide-set",
                            "guides oligos-list",
                        )?);
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for guides oligos-list"));
                    }
                }
            }
            Ok(ShellCommand::GuidesOligosList { guide_set_id })
        }
        "oligos-show" => {
            if tokens.len() != 3 {
                return Err("guides oligos-show requires OLIGO_SET_ID".to_string());
            }
            Ok(ShellCommand::GuidesOligosShow {
                oligo_set_id: tokens[2].clone(),
            })
        }
        "oligos-export" => {
            if tokens.len() < 4 {
                return Err(
                    "guides oligos-export requires GUIDE_SET_ID OUTPUT_PATH [--format csv_table|plate_csv|fasta] [--plate 96|384] [--oligo-set ID]"
                        .to_string(),
                );
            }
            let guide_set_id = tokens[2].clone();
            let path = tokens[3].clone();
            let mut format = GuideOligoExportFormat::CsvTable;
            let mut plate_format: Option<GuideOligoPlateFormat> = None;
            let mut oligo_set_id: Option<String> = None;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--format" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--format",
                            "guides oligos-export",
                        )?;
                        format = parse_guide_export_format(&raw)?;
                    }
                    "--plate" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--plate", "guides oligos-export")?;
                        plate_format = Some(parse_guide_plate_format(&raw)?);
                    }
                    "--oligo-set" => {
                        oligo_set_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--oligo-set",
                            "guides oligos-export",
                        )?);
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for guides oligos-export"));
                    }
                }
            }
            Ok(ShellCommand::GuidesOligosExport {
                guide_set_id,
                oligo_set_id,
                format,
                path,
                plate_format,
            })
        }
        "protocol-export" => {
            if tokens.len() < 4 {
                return Err(
                    "guides protocol-export requires GUIDE_SET_ID OUTPUT_PATH [--oligo-set ID] [--no-qc]"
                        .to_string(),
                );
            }
            let guide_set_id = tokens[2].clone();
            let path = tokens[3].clone();
            let mut oligo_set_id: Option<String> = None;
            let mut include_qc_checklist = true;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--oligo-set" => {
                        oligo_set_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--oligo-set",
                            "guides protocol-export",
                        )?);
                    }
                    "--no-qc" => {
                        include_qc_checklist = false;
                        idx += 1;
                    }
                    "--with-qc" => {
                        include_qc_checklist = true;
                        idx += 1;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for guides protocol-export"
                        ));
                    }
                }
            }
            Ok(ShellCommand::GuidesProtocolExport {
                guide_set_id,
                oligo_set_id,
                path,
                include_qc_checklist,
            })
        }
        other => Err(format!(
            "Unknown guides subcommand '{other}' (expected list, show, put, delete, filter, filter-show, oligos-generate, oligos-list, oligos-show, oligos-export, protocol-export)"
        )),
    }
}

fn parse_macros_command(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.len() < 2 {
        return Err(
            "macros requires a subcommand: run, template-list, template-show, template-put, template-delete, template-import, template-run"
                .to_string(),
        );
    }
    match tokens[1].as_str() {
        "run" => {
            if tokens.len() < 3 {
                return Err(
                    "macros run requires SCRIPT_OR_@FILE (or --file PATH), optionally with --transactional".to_string(),
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
                            return Err("macros run --file may only be specified once".to_string());
                        }
                        idx += 1;
                        if idx >= tokens.len() {
                            return Err("macros run --file requires PATH".to_string());
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
                        "macros run does not accept inline script after --file PATH".to_string()
                    );
                }
                format!("@{path}")
            } else {
                if idx >= tokens.len() {
                    return Err("macros run requires SCRIPT_OR_@FILE".to_string());
                }
                tokens[idx..].join(" ")
            };
            if script.trim().is_empty() {
                return Err("macros run requires non-empty script".to_string());
            }
            Ok(ShellCommand::MacrosRun {
                script,
                transactional,
            })
        }
        "template-list" => {
            if tokens.len() != 2 {
                return Err("macros template-list takes no options".to_string());
            }
            Ok(ShellCommand::MacrosTemplateList)
        }
        "template-show" => {
            if tokens.len() != 3 {
                return Err("macros template-show requires TEMPLATE_NAME".to_string());
            }
            Ok(ShellCommand::MacrosTemplateShow {
                name: tokens[2].clone(),
            })
        }
        "template-put" | "template-upsert" => {
            if tokens.len() < 4 {
                return Err(
                    "macros template-put requires TEMPLATE_NAME (--script SCRIPT_OR_@FILE | --file PATH) [--description TEXT] [--details-url URL] [--param NAME|NAME=DEFAULT ...]"
                        .to_string(),
                );
            }
            let name = tokens[2].clone();
            let mut description: Option<String> = None;
            let mut details_url: Option<String> = None;
            let mut parameters: Vec<WorkflowMacroTemplateParam> = vec![];
            let mut script: Option<String> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--description" => {
                        description = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--description",
                            "macros template-put",
                        )?);
                    }
                    "--details-url" | "--url" => {
                        if details_url.is_some() {
                            return Err(
                                "macros template-put details URL was already specified".to_string()
                            );
                        }
                        details_url = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--details-url",
                            "macros template-put",
                        )?);
                    }
                    "--param" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--param", "macros template-put")?;
                        parameters.push(parse_workflow_template_param_spec(&raw)?);
                    }
                    "--script" => {
                        if script.is_some() {
                            return Err(
                                "macros template-put script was already specified".to_string()
                            );
                        }
                        script = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--script",
                            "macros template-put",
                        )?);
                    }
                    "--file" => {
                        if script.is_some() {
                            return Err(
                                "macros template-put script was already specified".to_string()
                            );
                        }
                        let path =
                            parse_option_path(tokens, &mut idx, "--file", "macros template-put")?;
                        script = Some(format!("@{path}"));
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for macros template-put"));
                    }
                }
            }
            let script = script.ok_or_else(|| {
                "macros template-put requires --script SCRIPT_OR_@FILE or --file PATH".to_string()
            })?;
            Ok(ShellCommand::MacrosTemplateUpsert {
                name,
                description,
                details_url,
                parameters,
                script,
            })
        }
        "template-delete" => {
            if tokens.len() != 3 {
                return Err("macros template-delete requires TEMPLATE_NAME".to_string());
            }
            Ok(ShellCommand::MacrosTemplateDelete {
                name: tokens[2].clone(),
            })
        }
        "template-import" => {
            if tokens.len() != 3 {
                return Err("macros template-import requires PATH".to_string());
            }
            Ok(ShellCommand::MacrosTemplateImport {
                path: tokens[2].clone(),
            })
        }
        "template-run" => {
            if tokens.len() < 3 {
                return Err(
                    "macros template-run requires TEMPLATE_NAME [--bind KEY=VALUE ...] [--transactional]"
                        .to_string(),
                );
            }
            let name = tokens[2].clone();
            let mut bindings: HashMap<String, String> = HashMap::new();
            let mut transactional = false;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--bind" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--bind", "macros template-run")?;
                        let (key, value) = parse_template_binding(&raw)?;
                        if bindings.insert(key.clone(), value).is_some() {
                            return Err(format!(
                                "Duplicate --bind key '{}' in macros template-run",
                                key
                            ));
                        }
                    }
                    "--transactional" | "--atomic" => {
                        transactional = true;
                        idx += 1;
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for macros template-run"));
                    }
                }
            }
            Ok(ShellCommand::MacrosTemplateRun {
                name,
                bindings,
                transactional,
            })
        }
        other => Err(format!(
            "Unknown macros subcommand '{other}' (expected run, template-list, template-show, template-put, template-delete, template-import, template-run)"
        )),
    }
}

fn parse_routines_command(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.len() < 2 {
        return Err("routines requires a subcommand: list".to_string());
    }
    match tokens[1].as_str() {
        "list" => {
            let mut catalog_path: Option<String> = None;
            let mut family: Option<String> = None;
            let mut status: Option<String> = None;
            let mut tag: Option<String> = None;
            let mut query: Option<String> = None;
            let mut idx = 2usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--catalog" => {
                        catalog_path = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--catalog",
                            "routines list",
                        )?);
                    }
                    "--family" => {
                        family = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--family",
                            "routines list",
                        )?);
                    }
                    "--status" => {
                        status = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--status",
                            "routines list",
                        )?);
                    }
                    "--tag" => {
                        tag = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--tag",
                            "routines list",
                        )?);
                    }
                    "--query" => {
                        query = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--query",
                            "routines list",
                        )?);
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for routines list"));
                    }
                }
            }
            Ok(ShellCommand::RoutinesList {
                catalog_path,
                family,
                status,
                tag,
                query,
            })
        }
        other => Err(format!(
            "Unknown routines subcommand '{other}' (expected list)"
        )),
    }
}

fn parse_agents_command(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.len() < 2 {
        return Err("agents requires a subcommand: list or ask".to_string());
    }
    match tokens[1].as_str() {
        "list" => {
            let mut catalog_path: Option<String> = None;
            let mut idx = 2usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--catalog" => {
                        catalog_path = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--catalog",
                            "agents list",
                        )?);
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for agents list"));
                    }
                }
            }
            Ok(ShellCommand::AgentsList { catalog_path })
        }
        "ask" => {
            if tokens.len() < 3 {
                return Err(
                    "agents ask requires SYSTEM_ID --prompt TEXT [--catalog PATH] [--base-url URL] [--model MODEL] [--timeout-secs N] [--connect-timeout-secs N] [--read-timeout-secs N] [--max-retries N] [--max-response-bytes N] [--allow-auto-exec] [--execute-all] [--execute-index N ...] [--no-state-summary]".to_string(),
                );
            }
            let system_id = tokens[2].trim().to_string();
            if system_id.is_empty() {
                return Err("agents ask SYSTEM_ID cannot be empty".to_string());
            }

            let mut prompt: Option<String> = None;
            let mut catalog_path: Option<String> = None;
            let mut base_url_override: Option<String> = None;
            let mut model_override: Option<String> = None;
            let mut timeout_seconds: Option<u64> = None;
            let mut connect_timeout_seconds: Option<u64> = None;
            let mut read_timeout_seconds: Option<u64> = None;
            let mut max_retries: Option<usize> = None;
            let mut max_response_bytes: Option<usize> = None;
            let mut include_state_summary = true;
            let mut allow_auto_exec = false;
            let mut execute_all = false;
            let mut execute_indices: Vec<usize> = vec![];
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--prompt" => {
                        prompt = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--prompt",
                            "agents ask",
                        )?);
                    }
                    "--catalog" => {
                        catalog_path = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--catalog",
                            "agents ask",
                        )?);
                    }
                    "--base-url" => {
                        base_url_override = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--base-url",
                            "agents ask",
                        )?);
                    }
                    "--model" => {
                        model_override = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--model",
                            "agents ask",
                        )?);
                    }
                    "--timeout-secs" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--timeout-secs", "agents ask")?;
                        let parsed = raw
                            .parse::<u64>()
                            .map_err(|e| format!("Invalid --timeout-secs value '{raw}': {e}"))?;
                        if parsed == 0 {
                            timeout_seconds = None;
                        } else {
                            timeout_seconds = Some(parsed);
                        }
                    }
                    "--connect-timeout-secs" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--connect-timeout-secs",
                            "agents ask",
                        )?;
                        let parsed = raw.parse::<u64>().map_err(|e| {
                            format!("Invalid --connect-timeout-secs value '{raw}': {e}")
                        })?;
                        if parsed == 0 {
                            connect_timeout_seconds = None;
                        } else {
                            connect_timeout_seconds = Some(parsed);
                        }
                    }
                    "--read-timeout-secs" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--read-timeout-secs",
                            "agents ask",
                        )?;
                        let parsed = raw.parse::<u64>().map_err(|e| {
                            format!("Invalid --read-timeout-secs value '{raw}': {e}")
                        })?;
                        if parsed == 0 {
                            read_timeout_seconds = None;
                        } else {
                            read_timeout_seconds = Some(parsed);
                        }
                    }
                    "--max-retries" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--max-retries", "agents ask")?;
                        let parsed = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --max-retries value '{raw}': {e}"))?;
                        max_retries = Some(parsed);
                    }
                    "--max-response-bytes" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--max-response-bytes",
                            "agents ask",
                        )?;
                        let parsed = raw.parse::<usize>().map_err(|e| {
                            format!("Invalid --max-response-bytes value '{raw}': {e}")
                        })?;
                        if parsed == 0 {
                            max_response_bytes = None;
                        } else {
                            max_response_bytes = Some(parsed);
                        }
                    }
                    "--allow-auto-exec" => {
                        allow_auto_exec = true;
                        idx += 1;
                    }
                    "--execute-all" => {
                        execute_all = true;
                        idx += 1;
                    }
                    "--execute-index" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--execute-index", "agents ask")?;
                        let index = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --execute-index value '{raw}': {e}"))?;
                        if index == 0 {
                            return Err("--execute-index must be >= 1".to_string());
                        }
                        execute_indices.push(index);
                    }
                    "--no-state-summary" => {
                        include_state_summary = false;
                        idx += 1;
                    }
                    "--with-state-summary" => {
                        include_state_summary = true;
                        idx += 1;
                    }
                    other => return Err(format!("Unknown option '{other}' for agents ask")),
                }
            }

            let prompt = prompt
                .map(|value| value.trim().to_string())
                .filter(|value| !value.is_empty())
                .ok_or_else(|| "agents ask requires non-empty --prompt TEXT".to_string())?;
            execute_indices.sort_unstable();
            execute_indices.dedup();

            Ok(ShellCommand::AgentsAsk {
                system_id,
                prompt,
                catalog_path,
                base_url_override,
                model_override,
                timeout_seconds,
                connect_timeout_seconds,
                read_timeout_seconds,
                max_retries,
                max_response_bytes,
                include_state_summary,
                allow_auto_exec,
                execute_all,
                execute_indices,
            })
        }
        other => Err(format!(
            "Unknown agents subcommand '{other}' (expected list or ask)"
        )),
    }
}

fn parse_ui_command(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.len() < 2 {
        return Err(
            "ui requires a subcommand: intents, open, focus, prepared-genomes, latest-prepared"
                .to_string(),
        );
    }
    match tokens[1].as_str() {
        "intents" | "list-intents" => {
            if tokens.len() > 2 {
                return Err("ui intents takes no options".to_string());
            }
            Ok(ShellCommand::UiListIntents)
        }
        "open" | "focus" => {
            if tokens.len() < 3 {
                return Err(
                    "ui open|focus requires TARGET [--genome-id GENOME_ID] [--helpers] [--catalog PATH] [--cache-dir PATH] [--filter TEXT] [--species TEXT] [--latest]"
                        .to_string(),
                );
            }
            let action = if tokens[1] == "open" {
                UiIntentAction::Open
            } else {
                UiIntentAction::Focus
            };
            let target = UiIntentTarget::parse(&tokens[2]).ok_or_else(|| {
                "Unknown ui target. Expected one of: prepared-references, prepare-reference-genome, retrieve-genome-sequence, blast-genome-sequence, import-genome-track, agent-assistant, prepare-helper-genome, retrieve-helper-sequence, blast-helper-sequence".to_string()
            })?;
            let mut genome_id: Option<String> = None;
            let mut helper_mode = false;
            let mut catalog_path: Option<String> = None;
            let mut cache_dir: Option<String> = None;
            let mut filter: Option<String> = None;
            let mut species: Option<String> = None;
            let mut latest = false;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--genome-id" => {
                        genome_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--genome-id",
                            "ui open|focus",
                        )?);
                    }
                    "--helpers" => {
                        helper_mode = true;
                        idx += 1;
                    }
                    "--catalog" => {
                        catalog_path = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--catalog",
                            "ui open|focus",
                        )?);
                    }
                    "--cache-dir" => {
                        cache_dir = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--cache-dir",
                            "ui open|focus",
                        )?);
                    }
                    "--filter" => {
                        filter = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--filter",
                            "ui open|focus",
                        )?);
                    }
                    "--species" => {
                        species = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--species",
                            "ui open|focus",
                        )?);
                    }
                    "--latest" => {
                        latest = true;
                        idx += 1;
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for ui {}", tokens[1]));
                    }
                }
            }
            let genome_id = genome_id
                .map(|value| value.trim().to_string())
                .filter(|value| !value.is_empty());
            let filter = filter
                .map(|value| value.trim().to_string())
                .filter(|value| !value.is_empty());
            let species = species
                .map(|value| value.trim().to_string())
                .filter(|value| !value.is_empty());
            if target != UiIntentTarget::PreparedReferences
                && (helper_mode
                    || catalog_path.is_some()
                    || cache_dir.is_some()
                    || filter.is_some()
                    || species.is_some()
                    || latest)
            {
                return Err(
                    "ui open|focus TARGET only supports --helpers/--catalog/--cache-dir/--filter/--species/--latest when TARGET is prepared-references"
                        .to_string(),
                );
            }
            Ok(ShellCommand::UiIntent {
                action,
                target,
                genome_id,
                helper_mode,
                catalog_path,
                cache_dir,
                filter,
                species,
                latest,
            })
        }
        "prepared-genomes" => {
            let mut helper_mode = false;
            let mut catalog_path: Option<String> = None;
            let mut cache_dir: Option<String> = None;
            let mut filter: Option<String> = None;
            let mut species: Option<String> = None;
            let mut latest = false;
            let mut idx = 2usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--helpers" => {
                        helper_mode = true;
                        idx += 1;
                    }
                    "--catalog" => {
                        catalog_path = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--catalog",
                            "ui prepared-genomes",
                        )?);
                    }
                    "--cache-dir" => {
                        cache_dir = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--cache-dir",
                            "ui prepared-genomes",
                        )?);
                    }
                    "--filter" => {
                        filter = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--filter",
                            "ui prepared-genomes",
                        )?);
                    }
                    "--species" => {
                        species = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--species",
                            "ui prepared-genomes",
                        )?);
                    }
                    "--latest" => {
                        latest = true;
                        idx += 1;
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for ui prepared-genomes"));
                    }
                }
            }
            let filter = filter
                .map(|value| value.trim().to_string())
                .filter(|value| !value.is_empty());
            let species = species
                .map(|value| value.trim().to_string())
                .filter(|value| !value.is_empty());
            Ok(ShellCommand::UiPreparedGenomes {
                helper_mode,
                catalog_path,
                cache_dir,
                filter,
                species,
                latest,
            })
        }
        "latest-prepared" => {
            if tokens.len() < 3 {
                return Err(
                    "ui latest-prepared requires SPECIES [--helpers] [--catalog PATH] [--cache-dir PATH]"
                        .to_string(),
                );
            }
            if tokens[2].starts_with("--") {
                return Err("ui latest-prepared requires SPECIES before option flags".to_string());
            }
            let species = tokens[2].trim().to_string();
            if species.is_empty() {
                return Err("ui latest-prepared SPECIES must not be empty".to_string());
            }
            let mut helper_mode = false;
            let mut catalog_path: Option<String> = None;
            let mut cache_dir: Option<String> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--helpers" => {
                        helper_mode = true;
                        idx += 1;
                    }
                    "--catalog" => {
                        catalog_path = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--catalog",
                            "ui latest-prepared",
                        )?);
                    }
                    "--cache-dir" => {
                        cache_dir = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--cache-dir",
                            "ui latest-prepared",
                        )?);
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for ui latest-prepared"));
                    }
                }
            }
            Ok(ShellCommand::UiLatestPrepared {
                helper_mode,
                catalog_path,
                cache_dir,
                species,
            })
        }
        other => Err(format!(
            "Unknown ui subcommand '{other}' (expected intents, open, focus, prepared-genomes, latest-prepared)"
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
        "inspect-feature-expert" => {
            if tokens.len() < 4 {
                return Err(token_error(cmd));
            }
            let seq_id = tokens[1].clone();
            let target =
                parse_feature_expert_target_tokens(&tokens[2..], "inspect-feature-expert")?;
            Ok(ShellCommand::InspectFeatureExpert { seq_id, target })
        }
        "render-feature-expert-svg" => {
            if tokens.len() < 5 {
                return Err(token_error(cmd));
            }
            let seq_id = tokens[1].clone();
            let output = tokens[tokens.len() - 1].clone();
            let target = parse_feature_expert_target_tokens(
                &tokens[2..tokens.len() - 1],
                "render-feature-expert-svg",
            )?;
            Ok(ShellCommand::RenderFeatureExpertSvg {
                seq_id,
                target,
                output,
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
        "render-pool-gel-svg" | "render-gel-svg" => {
            if tokens.len() < 3 {
                return Err(token_error(cmd));
            }
            let cmd_name = tokens[0].as_str();
            let inputs = match tokens[1].trim() {
                "-" | "_" => vec![],
                raw => split_ids(raw),
            };
            let output = tokens[2].clone();
            let mut ladders: Option<Vec<String>> = None;
            let mut container_ids: Option<Vec<String>> = None;
            let mut arrangement_id: Option<String> = None;
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
                    "--containers" => {
                        if idx + 1 >= tokens.len() {
                            return Err("Missing value after --containers".to_string());
                        }
                        let parsed = split_ids(&tokens[idx + 1]);
                        if !parsed.is_empty() {
                            container_ids = Some(parsed);
                        }
                        idx += 2;
                    }
                    "--arrangement" => {
                        if idx + 1 >= tokens.len() {
                            return Err("Missing value after --arrangement".to_string());
                        }
                        let value = tokens[idx + 1].trim();
                        if !value.is_empty() {
                            arrangement_id = Some(value.to_string());
                        }
                        idx += 2;
                    }
                    other => {
                        return Err(format!("Unknown argument '{other}' for {cmd_name}"));
                    }
                }
            }
            if inputs.is_empty()
                && container_ids.as_ref().map_or(true, |v| v.is_empty())
                && arrangement_id.is_none()
            {
                return Err(format!(
                    "{cmd_name} requires inputs, --containers, or --arrangement"
                ));
            }
            Ok(ShellCommand::RenderPoolGelSvg {
                inputs,
                output,
                ladders,
                container_ids,
                arrangement_id,
            })
        }
        "arrange-serial" => {
            if tokens.len() < 2 {
                return Err(
                    "arrange-serial requires: CONTAINER_IDS [--id ARR_ID] [--name TEXT] [--ladders NAME[,NAME]]"
                        .to_string(),
                );
            }
            let container_ids = split_ids(&tokens[1]);
            if container_ids.is_empty() {
                return Err("arrange-serial requires at least one container id".to_string());
            }
            let mut arrangement_id: Option<String> = None;
            let mut name: Option<String> = None;
            let mut ladders: Option<Vec<String>> = None;
            let mut idx = 2usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--id" => {
                        if idx + 1 >= tokens.len() {
                            return Err("Missing value after --id".to_string());
                        }
                        let value = tokens[idx + 1].trim();
                        if !value.is_empty() {
                            arrangement_id = Some(value.to_string());
                        }
                        idx += 2;
                    }
                    "--name" => {
                        if idx + 1 >= tokens.len() {
                            return Err("Missing value after --name".to_string());
                        }
                        let value = tokens[idx + 1].trim();
                        if !value.is_empty() {
                            name = Some(value.to_string());
                        }
                        idx += 2;
                    }
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
                        return Err(format!("Unknown argument '{other}' for arrange-serial"));
                    }
                }
            }
            Ok(ShellCommand::CreateArrangementSerial {
                container_ids,
                arrangement_id,
                name,
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
        "ui" => parse_ui_command(tokens),
        "agents" => parse_agents_command(tokens),
        "routines" => parse_routines_command(tokens),
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
                                                ));
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
                                source: source
                                    .unwrap_or_else(|| GenomeTrackSource::from_path(&path)),
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
                                return Err("tracks tracked remove requires INDEX".to_string());
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
        "macros" => parse_macros_command(tokens),
        "candidates" => parse_candidates_command(tokens),
        "guides" => parse_guides_command(tokens),
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

fn load_macro_script(script_or_ref: &str, label: &str) -> Result<String, String> {
    let trimmed = script_or_ref.trim();
    if let Some(path) = trimmed.strip_prefix('@') {
        let path = path.trim();
        if path.is_empty() {
            return Err(format!("{label} @FILE requires a non-empty file path"));
        }
        fs::read_to_string(path).map_err(|e| format!("Could not read {label} file '{path}': {e}"))
    } else {
        Ok(trimmed.to_string())
    }
}

fn load_candidates_macro_script(script_or_ref: &str) -> Result<String, String> {
    load_macro_script(script_or_ref, "candidates macro")
}

fn load_workflow_macro_script(script_or_ref: &str) -> Result<String, String> {
    load_macro_script(script_or_ref, "macros run")
}

fn split_macro_statements(script: &str) -> Result<Vec<String>, String> {
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

fn split_candidates_macro_statements(script: &str) -> Result<Vec<String>, String> {
    split_macro_statements(script)
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
        if matches!(
            cmd,
            ShellCommand::CandidatesMacro { .. }
                | ShellCommand::CandidatesTemplateRun { .. }
                | ShellCommand::MacrosRun { .. }
                | ShellCommand::MacrosTemplateRun { .. }
        ) {
            return Err("Nested candidates macro/template-run calls are not allowed".to_string());
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

fn run_workflow_macro(
    engine: &mut GentleEngine,
    script_or_ref: &str,
    transactional: bool,
    options: &ShellExecutionOptions,
) -> Result<ShellRunResult, String> {
    let script = load_workflow_macro_script(script_or_ref)?;
    let statements = split_macro_statements(&script)?;
    if statements.is_empty() {
        return Err("macros run script is empty".to_string());
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
        let cmd = if let Some(raw_payload) = statement.strip_prefix("op ") {
            let payload = raw_payload.trim();
            if payload.is_empty() {
                return Err("macros run statement 'op' requires a JSON payload".to_string());
            }
            ShellCommand::Op {
                payload: payload.to_string(),
            }
        } else if let Some(raw_payload) = statement.strip_prefix("workflow ") {
            let payload = raw_payload.trim();
            if payload.is_empty() {
                return Err("macros run statement 'workflow' requires a JSON payload".to_string());
            }
            ShellCommand::Workflow {
                payload: payload.to_string(),
            }
        } else {
            let tokens = split_shell_words(statement)?;
            parse_shell_tokens(&tokens)?
        };
        if matches!(
            cmd,
            ShellCommand::CandidatesMacro { .. }
                | ShellCommand::CandidatesTemplateRun { .. }
                | ShellCommand::MacrosRun { .. }
                | ShellCommand::MacrosTemplateRun { .. }
        ) {
            return Err("Nested macros/template-run calls are not allowed".to_string());
        }
        let run = match execute_shell_command_with_options(engine, &cmd, options) {
            Ok(run) => run,
            Err(err) => {
                if transactional {
                    if let Some(state) = rollback_state {
                        *engine = GentleEngine::from_state(state);
                    }
                    return Err(format!(
                        "macros run failed at statement {} ('{}'): {err}; all macro changes were rolled back",
                        executed + 1,
                        statement
                    ));
                }
                return Err(format!(
                    "macros run failed at statement {} ('{}'): {err}",
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
        return Err("macros run script has no executable statements".to_string());
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

#[derive(Debug, Clone, Serialize, Deserialize)]
struct AgentSuggestedExecutionReport {
    index: usize,
    command: String,
    execution_intent: String,
    trigger: String,
    executed: bool,
    ok: bool,
    state_changed: bool,
    error: Option<String>,
    output: Option<Value>,
}

fn should_execute_agent_suggestion(
    index_1based: usize,
    intent: AgentExecutionIntent,
    execute_all: bool,
    execute_indices: &BTreeSet<usize>,
    allow_auto_exec: bool,
) -> Option<&'static str> {
    if execute_all {
        return Some("execute_all");
    }
    if execute_indices.contains(&index_1based) {
        return Some("execute_index");
    }
    if allow_auto_exec && intent == AgentExecutionIntent::Auto {
        return Some("allow_auto_exec");
    }
    None
}

fn execute_agent_suggested_commands(
    engine: &mut GentleEngine,
    suggestions: &[crate::agent_bridge::AgentSuggestedCommand],
    execute_all: bool,
    execute_indices: &BTreeSet<usize>,
    allow_auto_exec: bool,
    options: &ShellExecutionOptions,
) -> (bool, Vec<AgentSuggestedExecutionReport>) {
    let mut changed = false;
    let mut rows: Vec<AgentSuggestedExecutionReport> = vec![];
    let nested_options = ShellExecutionOptions {
        allow_screenshots: options.allow_screenshots,
        allow_agent_commands: false,
    };
    for (idx, suggestion) in suggestions.iter().enumerate() {
        let index_1based = idx + 1;
        let trigger = should_execute_agent_suggestion(
            index_1based,
            suggestion.execution,
            execute_all,
            execute_indices,
            allow_auto_exec,
        );
        if trigger.is_none() {
            rows.push(AgentSuggestedExecutionReport {
                index: index_1based,
                command: suggestion.command.clone(),
                execution_intent: suggestion.execution.as_str().to_string(),
                trigger: "none".to_string(),
                executed: false,
                ok: true,
                state_changed: false,
                error: None,
                output: None,
            });
            continue;
        }
        let command_text = suggestion.command.trim().to_string();
        if command_text.is_empty() {
            rows.push(AgentSuggestedExecutionReport {
                index: index_1based,
                command: suggestion.command.clone(),
                execution_intent: suggestion.execution.as_str().to_string(),
                trigger: trigger.unwrap_or("unknown").to_string(),
                executed: false,
                ok: false,
                state_changed: false,
                error: Some("Suggested command is empty".to_string()),
                output: None,
            });
            continue;
        }
        let parsed = match parse_shell_line(&command_text) {
            Ok(command) => command,
            Err(err) => {
                rows.push(AgentSuggestedExecutionReport {
                    index: index_1based,
                    command: command_text,
                    execution_intent: suggestion.execution.as_str().to_string(),
                    trigger: trigger.unwrap_or("unknown").to_string(),
                    executed: true,
                    ok: false,
                    state_changed: false,
                    error: Some(format!("Could not parse suggested command: {err}")),
                    output: None,
                });
                continue;
            }
        };
        if matches!(parsed, ShellCommand::AgentsAsk { .. }) {
            rows.push(AgentSuggestedExecutionReport {
                index: index_1based,
                command: command_text,
                execution_intent: suggestion.execution.as_str().to_string(),
                trigger: trigger.unwrap_or("unknown").to_string(),
                executed: true,
                ok: false,
                state_changed: false,
                error: Some(
                    "Agent-to-agent 'agents ask' execution is blocked for suggested commands"
                        .to_string(),
                ),
                output: None,
            });
            continue;
        }
        match execute_shell_command_with_options(engine, &parsed, &nested_options) {
            Ok(run) => {
                changed |= run.state_changed;
                rows.push(AgentSuggestedExecutionReport {
                    index: index_1based,
                    command: command_text,
                    execution_intent: suggestion.execution.as_str().to_string(),
                    trigger: trigger.unwrap_or("unknown").to_string(),
                    executed: true,
                    ok: true,
                    state_changed: run.state_changed,
                    error: None,
                    output: Some(run.output),
                });
            }
            Err(err) => {
                rows.push(AgentSuggestedExecutionReport {
                    index: index_1based,
                    command: command_text,
                    execution_intent: suggestion.execution.as_str().to_string(),
                    trigger: trigger.unwrap_or("unknown").to_string(),
                    executed: true,
                    ok: false,
                    state_changed: false,
                    error: Some(err),
                    output: None,
                });
            }
        }
    }
    (changed, rows)
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
        ShellCommand::InspectFeatureExpert { seq_id, target } => {
            let view = engine
                .inspect_feature_expert(seq_id, target)
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: serde_json::to_value(view)
                    .map_err(|e| format!("Could not serialize feature expert view: {e}"))?,
            }
        }
        ShellCommand::RenderFeatureExpertSvg {
            seq_id,
            target,
            output,
        } => {
            let op_result = engine
                .apply(Operation::RenderFeatureExpertSvg {
                    seq_id: seq_id.clone(),
                    target: target.clone(),
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
            container_ids,
            arrangement_id,
        } => {
            let op_result = engine
                .apply(Operation::RenderPoolGelSvg {
                    inputs: inputs.clone(),
                    path: output.clone(),
                    ladders: ladders.clone(),
                    container_ids: container_ids.clone(),
                    arrangement_id: arrangement_id.clone(),
                })
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::CreateArrangementSerial {
            container_ids,
            arrangement_id,
            name,
            ladders,
        } => {
            let op_result = engine
                .apply(Operation::CreateArrangementSerial {
                    container_ids: container_ids.clone(),
                    arrangement_id: arrangement_id.clone(),
                    name: name.clone(),
                    ladders: ladders.clone(),
                })
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: true,
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
            tf_motifs::reload_from_path(Some(report.output.as_str()));
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "message": format!("Synced {} {} entries to '{}'", report.item_count, report.resource, report.output),
                    "report": report,
                }),
            }
        }
        ShellCommand::RoutinesList {
            catalog_path,
            family,
            status,
            tag,
            query,
        } => {
            let resolved_catalog = catalog_path
                .as_deref()
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .unwrap_or(DEFAULT_CLONING_ROUTINE_CATALOG_PATH)
                .to_string();
            let catalog = load_cloning_routine_catalog(&resolved_catalog)?;
            let catalog_schema = catalog.schema.clone();

            let mut available_families = catalog
                .routines
                .iter()
                .map(|routine| routine.family.trim().to_string())
                .filter(|family| !family.is_empty())
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect::<Vec<_>>();
            available_families.sort_by_key(|value| value.to_ascii_lowercase());

            let mut available_statuses = catalog
                .routines
                .iter()
                .map(|routine| routine.status.trim().to_string())
                .filter(|status| !status.is_empty())
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect::<Vec<_>>();
            available_statuses.sort_by_key(|value| value.to_ascii_lowercase());

            let mut routines = catalog
                .routines
                .into_iter()
                .filter(|routine| {
                    routine_matches_filter(
                        routine,
                        family.as_deref(),
                        status.as_deref(),
                        tag.as_deref(),
                        query.as_deref(),
                    )
                })
                .collect::<Vec<_>>();
            routines.sort_by(|left, right| {
                left.family
                    .to_ascii_lowercase()
                    .cmp(&right.family.to_ascii_lowercase())
                    .then(
                        left.title
                            .to_ascii_lowercase()
                            .cmp(&right.title.to_ascii_lowercase()),
                    )
                    .then(
                        left.routine_id
                            .to_ascii_lowercase()
                            .cmp(&right.routine_id.to_ascii_lowercase()),
                    )
            });
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": CLONING_ROUTINE_LIST_SCHEMA,
                    "catalog_path": resolved_catalog,
                    "catalog_schema": catalog_schema,
                    "filters": {
                        "family": family,
                        "status": status,
                        "tag": tag,
                        "query": query,
                    },
                    "available_families": available_families,
                    "available_statuses": available_statuses,
                    "routine_count": routines.len(),
                    "routines": routines,
                }),
            }
        }
        ShellCommand::AgentsList { catalog_path } => {
            let (resolved_catalog_path, catalog) =
                load_agent_system_catalog(catalog_path.as_deref())?;
            let systems = catalog
                .systems
                .iter()
                .map(|system| {
                    let availability = agent_system_availability(system);
                    json!({
                        "id": system.id,
                        "label": system.label,
                        "description": system.description,
                        "transport": system.transport.as_str(),
                        "command": system.command,
                        "model": system.model,
                        "base_url": system.base_url,
                        "working_dir": system.working_dir,
                        "available": availability.available,
                        "availability_reason": availability.reason,
                    })
                })
                .collect::<Vec<_>>();
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": "gentle.agent_systems_list.v1",
                    "catalog_path": resolved_catalog_path,
                    "catalog_schema": catalog.schema,
                    "system_count": systems.len(),
                    "systems": systems
                }),
            }
        }
        ShellCommand::AgentsAsk {
            system_id,
            prompt,
            catalog_path,
            base_url_override,
            model_override,
            timeout_seconds,
            connect_timeout_seconds,
            read_timeout_seconds,
            max_retries,
            max_response_bytes,
            include_state_summary,
            allow_auto_exec,
            execute_all,
            execute_indices,
        } => {
            if !options.allow_agent_commands {
                return Err(
                    "agents ask execution is blocked in this context (agent-to-agent recursion guardrail)"
                        .to_string(),
                );
            }
            let state_summary = if *include_state_summary {
                Some(engine.summarize_state())
            } else {
                None
            };
            let mut env_overrides = HashMap::new();
            if let Some(base_url) = base_url_override
                .as_deref()
                .map(str::trim)
                .filter(|value| !value.is_empty())
            {
                env_overrides.insert(AGENT_BASE_URL_ENV.to_string(), base_url.to_string());
            }
            if let Some(model) = model_override
                .as_deref()
                .map(str::trim)
                .filter(|value| !value.is_empty())
            {
                env_overrides.insert(AGENT_MODEL_ENV.to_string(), model.to_string());
            }
            if let Some(timeout_seconds) = timeout_seconds.filter(|value| *value > 0) {
                env_overrides.insert(
                    AGENT_TIMEOUT_SECS_ENV.to_string(),
                    timeout_seconds.to_string(),
                );
            }
            if let Some(connect_timeout_seconds) =
                connect_timeout_seconds.filter(|value| *value > 0)
            {
                env_overrides.insert(
                    AGENT_CONNECT_TIMEOUT_SECS_ENV.to_string(),
                    connect_timeout_seconds.to_string(),
                );
            }
            if let Some(read_timeout_seconds) = read_timeout_seconds.filter(|value| *value > 0) {
                env_overrides.insert(
                    AGENT_READ_TIMEOUT_SECS_ENV.to_string(),
                    read_timeout_seconds.to_string(),
                );
            }
            if let Some(max_retries) = max_retries {
                env_overrides.insert(AGENT_MAX_RETRIES_ENV.to_string(), max_retries.to_string());
            }
            if let Some(max_response_bytes) = max_response_bytes.filter(|value| *value > 0) {
                env_overrides.insert(
                    AGENT_MAX_RESPONSE_BYTES_ENV.to_string(),
                    max_response_bytes.to_string(),
                );
            }
            let invocation = invoke_agent_support_with_env_overrides(
                catalog_path.as_deref(),
                system_id,
                prompt,
                state_summary.as_ref(),
                if env_overrides.is_empty() {
                    None
                } else {
                    Some(&env_overrides)
                },
            )?;
            let suggested = invocation.response.suggested_commands.clone();
            let suggested_count = suggested.len();
            let requested_indices = execute_indices.iter().copied().collect::<BTreeSet<usize>>();
            let mut invalid_execute_indices = requested_indices
                .iter()
                .copied()
                .filter(|idx| *idx == 0 || *idx > suggested_count)
                .collect::<Vec<_>>();
            invalid_execute_indices.sort_unstable();
            let execute_index_set = requested_indices
                .iter()
                .copied()
                .filter(|idx| *idx >= 1 && *idx <= suggested_count)
                .collect::<BTreeSet<_>>();
            let (state_changed, execution_reports) = execute_agent_suggested_commands(
                engine,
                &suggested,
                *execute_all,
                &execute_index_set,
                *allow_auto_exec,
                options,
            );
            let executed_count = execution_reports.iter().filter(|row| row.executed).count();
            let executed_error_count = execution_reports
                .iter()
                .filter(|row| row.executed && !row.ok)
                .count();
            let auto_suggestion_count = suggested
                .iter()
                .filter(|item| item.execution == AgentExecutionIntent::Auto)
                .count();
            ShellRunResult {
                state_changed,
                output: json!({
                    "schema": "gentle.agent_ask_result.v1",
                    "invocation": invocation,
                    "request": {
                        "include_state_summary": include_state_summary,
                        "allow_auto_exec": allow_auto_exec,
                        "execute_all": execute_all,
                        "execute_indices": execute_indices,
                        "base_url_override": base_url_override,
                        "model_override": model_override,
                        "timeout_seconds": timeout_seconds,
                        "connect_timeout_seconds": connect_timeout_seconds,
                        "read_timeout_seconds": read_timeout_seconds,
                        "max_retries": max_retries,
                        "max_response_bytes": max_response_bytes,
                        "invalid_execute_indices": invalid_execute_indices,
                    },
                    "summary": {
                        "suggested_command_count": suggested_count,
                        "auto_suggestion_count": auto_suggestion_count,
                        "executed_count": executed_count,
                        "executed_error_count": executed_error_count,
                    },
                    "executions": execution_reports
                }),
            }
        }
        ShellCommand::UiListIntents => ShellRunResult {
            state_changed: false,
            output: json!({
                "schema": "gentle.ui_intents.v1",
                "targets": [
                    "prepared-references",
                    "prepare-reference-genome",
                    "retrieve-genome-sequence",
                    "blast-genome-sequence",
                    "import-genome-track",
                    "agent-assistant",
                    "prepare-helper-genome",
                    "retrieve-helper-sequence",
                    "blast-helper-sequence"
                ],
                "commands": [
                    "ui intents",
                    "ui open TARGET [--genome-id GENOME_ID] [--helpers] [--catalog PATH] [--cache-dir PATH] [--filter TEXT] [--species TEXT] [--latest]",
                    "ui focus TARGET [--genome-id GENOME_ID] [--helpers] [--catalog PATH] [--cache-dir PATH] [--filter TEXT] [--species TEXT] [--latest]",
                    "ui prepared-genomes [--helpers] [--catalog PATH] [--cache-dir PATH] [--filter TEXT] [--species TEXT] [--latest]",
                    "ui latest-prepared SPECIES [--helpers] [--catalog PATH] [--cache-dir PATH]"
                ],
                "notes": [
                    "UI intent commands are host-application intents and require GUI host integration to apply.",
                    "CLI execution returns deterministic intent/query payloads for agents/automation.",
                    "prepared-references target accepts query flags to resolve selected_genome_id; explicit --genome-id overrides query selection."
                ]
            }),
        },
        ShellCommand::UiIntent {
            action,
            target,
            genome_id,
            helper_mode,
            catalog_path,
            cache_dir,
            filter,
            species,
            latest,
        } => {
            let mut selected_genome_id = genome_id
                .as_ref()
                .map(|value| value.trim().to_string())
                .filter(|value| !value.is_empty());
            let mut prepared_query: Option<Value> = None;
            if matches!(target, UiIntentTarget::PreparedReferences) && selected_genome_id.is_none()
            {
                let prepared = execute_shell_command_with_options(
                    engine,
                    &ShellCommand::UiPreparedGenomes {
                        helper_mode: *helper_mode,
                        catalog_path: catalog_path.clone(),
                        cache_dir: cache_dir.clone(),
                        filter: filter.clone(),
                        species: species.clone(),
                        latest: *latest,
                    },
                    options,
                )?;
                selected_genome_id = prepared
                    .output
                    .get("selected_genome_id")
                    .and_then(|v| v.as_str())
                    .map(|v| v.to_string());
                prepared_query = Some(prepared.output);
            }
            let message = if matches!(target, UiIntentTarget::PreparedReferences) {
                if selected_genome_id.is_some() {
                    "UI intent recorded; prepared-references selection resolved deterministically."
                } else {
                    "UI intent recorded; prepared-references selection found no prepared genome."
                }
            } else {
                "UI intent recorded; requires GUI host integration to apply."
            };
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": "gentle.ui_intent.v1",
                    "ui_intent": {
                        "action": action.as_str(),
                        "target": target.as_str(),
                        "genome_id": genome_id,
                        "helper_mode": helper_mode,
                        "catalog_path": catalog_path,
                        "cache_dir": cache_dir,
                        "filter": filter,
                        "species": species,
                        "latest": latest
                    },
                    "selected_genome_id": selected_genome_id,
                    "prepared_query": prepared_query,
                    "applied": false,
                    "message": message
                }),
            }
        }
        ShellCommand::UiPreparedGenomes {
            helper_mode,
            catalog_path,
            cache_dir,
            filter,
            species,
            latest,
        } => {
            let catalog_path_effective = effective_catalog_path(catalog_path, *helper_mode);
            let catalog = GenomeCatalog::from_json_file(&catalog_path_effective).map_err(|e| {
                format!("Could not open genome catalog '{catalog_path_effective}': {e}")
            })?;
            let cache_dir = cache_dir
                .as_deref()
                .map(str::trim)
                .filter(|v| !v.is_empty())
                .map(|v| v.to_string());
            let filter_lc = filter
                .as_ref()
                .map(|v| v.trim().to_ascii_lowercase())
                .filter(|v| !v.is_empty());
            let species_lc = species
                .as_ref()
                .map(|v| v.trim().to_ascii_lowercase())
                .filter(|v| !v.is_empty());

            let mut prepared_rows: Vec<Value> = vec![];
            for genome_id in catalog.list_genomes() {
                let prepared = catalog
                    .is_prepared(&genome_id, cache_dir.as_deref())
                    .map_err(|e| {
                        format!("Could not inspect preparation state for '{genome_id}': {e}")
                    })?;
                if !prepared {
                    continue;
                }
                let id_lc = genome_id.to_ascii_lowercase();
                if let Some(filter_lc) = filter_lc.as_ref() {
                    if !id_lc.contains(filter_lc) {
                        continue;
                    }
                }
                if let Some(species_lc) = species_lc.as_ref() {
                    if !id_lc.contains(species_lc) {
                        continue;
                    }
                }
                let installed_at_unix_ms = catalog
                    .inspect_prepared_genome(&genome_id, cache_dir.as_deref())
                    .ok()
                    .flatten()
                    .map(|inspection| inspection.installed_at_unix_ms);
                prepared_rows.push(json!({
                    "genome_id": genome_id,
                    "installed_at_unix_ms": installed_at_unix_ms
                }));
            }

            prepared_rows.sort_by(|left, right| {
                let left_id = left
                    .get("genome_id")
                    .and_then(|v| v.as_str())
                    .unwrap_or_default();
                let right_id = right
                    .get("genome_id")
                    .and_then(|v| v.as_str())
                    .unwrap_or_default();
                let left_installed = left
                    .get("installed_at_unix_ms")
                    .and_then(|v| v.as_u64())
                    .unwrap_or(0);
                let right_installed = right
                    .get("installed_at_unix_ms")
                    .and_then(|v| v.as_u64())
                    .unwrap_or(0);
                right_installed
                    .cmp(&left_installed)
                    .then(left_id.cmp(right_id))
            });

            let latest_row = prepared_rows.first().cloned();
            let rows = if *latest {
                latest_row.clone().into_iter().collect::<Vec<_>>()
            } else {
                prepared_rows.clone()
            };
            let selected_genome_id = latest_row
                .as_ref()
                .and_then(|v| v.get("genome_id"))
                .and_then(|v| v.as_str())
                .map(|v| v.to_string());
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": "gentle.ui_prepared_genomes.v1",
                    "helper_mode": helper_mode,
                    "catalog_path": catalog_path_effective,
                    "cache_dir": cache_dir,
                    "filter": filter,
                    "species": species,
                    "latest_only": latest,
                    "prepared_count": prepared_rows.len(),
                    "selected_genome_id": selected_genome_id,
                    "genomes": rows
                }),
            }
        }
        ShellCommand::UiLatestPrepared {
            helper_mode,
            catalog_path,
            cache_dir,
            species,
        } => {
            let prepared = execute_shell_command_with_options(
                engine,
                &ShellCommand::UiPreparedGenomes {
                    helper_mode: *helper_mode,
                    catalog_path: catalog_path.clone(),
                    cache_dir: cache_dir.clone(),
                    filter: None,
                    species: Some(species.clone()),
                    latest: true,
                },
                options,
            )?;
            let selected_genome_id = prepared
                .output
                .get("selected_genome_id")
                .and_then(|v| v.as_str())
                .map(|v| v.to_string());
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": "gentle.ui_latest_prepared.v1",
                    "helper_mode": helper_mode,
                    "species": species,
                    "catalog_path": prepared.output.get("catalog_path").cloned().unwrap_or(json!(null)),
                    "cache_dir": prepared.output.get("cache_dir").cloned().unwrap_or(json!(null)),
                    "selected_genome_id": selected_genome_id,
                    "prepared_query": prepared.output
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
                    "nucleotide_length_bp": source_plan.nucleotide_length_bp,
                    "molecular_mass_da": source_plan.molecular_mass_da,
                    "molecular_mass_source": source_plan.molecular_mass_source,
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
        ShellCommand::MacrosRun {
            script,
            transactional,
        } => run_workflow_macro(engine, script, *transactional, options)?,
        ShellCommand::MacrosTemplateList => {
            let templates = engine.list_workflow_macro_templates();
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": "gentle.workflow_macro_templates.v1",
                    "template_count": templates.len(),
                    "templates": templates
                }),
            }
        }
        ShellCommand::MacrosTemplateShow { name } => {
            let template = engine
                .get_workflow_macro_template(name)
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "template": template
                }),
            }
        }
        ShellCommand::MacrosTemplateUpsert {
            name,
            description,
            details_url,
            parameters,
            script,
        } => {
            let before = engine
                .state()
                .metadata
                .get(WORKFLOW_MACRO_TEMPLATES_METADATA_KEY)
                .cloned();
            let loaded_script = load_workflow_macro_script(script)?;
            let op_result = engine
                .apply(Operation::UpsertWorkflowMacroTemplate {
                    name: name.clone(),
                    description: description.clone(),
                    details_url: details_url.clone(),
                    parameters: parameters.clone(),
                    script: loaded_script.clone(),
                })
                .map_err(|e| e.to_string())?;
            let after = engine
                .state()
                .metadata
                .get(WORKFLOW_MACRO_TEMPLATES_METADATA_KEY)
                .cloned();
            ShellRunResult {
                state_changed: before != after,
                output: json!({
                    "name": name,
                    "description": description,
                    "details_url": details_url,
                    "parameter_count": parameters.len(),
                    "script_length": loaded_script.len(),
                    "result": op_result
                }),
            }
        }
        ShellCommand::MacrosTemplateDelete { name } => {
            let before = engine
                .state()
                .metadata
                .get(WORKFLOW_MACRO_TEMPLATES_METADATA_KEY)
                .cloned();
            let op_result = engine
                .apply(Operation::DeleteWorkflowMacroTemplate { name: name.clone() })
                .map_err(|e| e.to_string())?;
            let after = engine
                .state()
                .metadata
                .get(WORKFLOW_MACRO_TEMPLATES_METADATA_KEY)
                .cloned();
            ShellRunResult {
                state_changed: before != after,
                output: json!({
                    "name": name,
                    "result": op_result
                }),
            }
        }
        ShellCommand::MacrosTemplateImport { path } => {
            let (templates, source_files) = load_cloning_pattern_templates_from_path(path)?;
            let before = engine.state().clone();
            let mut imported_names: Vec<String> = vec![];
            for template in &templates {
                let script = load_workflow_macro_script(&template.script)?;
                let op = Operation::UpsertWorkflowMacroTemplate {
                    name: template.name.clone(),
                    description: template.description.clone(),
                    details_url: template.details_url.clone(),
                    parameters: template.parameters.clone(),
                    script,
                };
                if let Err(err) = engine.apply(op) {
                    *engine.state_mut() = before;
                    return Err(format!(
                        "Failed to import template '{}' from '{}': {}",
                        template.name, path, err.message
                    ));
                }
                imported_names.push(template.name.clone());
            }
            imported_names.sort();
            imported_names.dedup();
            ShellRunResult {
                state_changed: !imported_names.is_empty(),
                output: json!({
                    "schema": CLONING_PATTERN_FILE_SCHEMA,
                    "path": path,
                    "source_files": source_files,
                    "imported_count": imported_names.len(),
                    "templates": imported_names,
                }),
            }
        }
        ShellCommand::MacrosTemplateRun {
            name,
            bindings,
            transactional,
        } => {
            let script = engine
                .render_workflow_macro_template_script(name, bindings)
                .map_err(|e| e.to_string())?;
            let mut run = run_workflow_macro(engine, &script, *transactional, options)?;
            run.output = json!({
                "template_name": name,
                "bindings": bindings,
                "expanded_script": script,
                "run": run.output
            });
            run
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
        ShellCommand::CandidatesScoreWeightedObjective {
            set_name,
            metric,
            objectives,
            normalize_metrics,
        } => {
            let before = engine
                .state()
                .metadata
                .get(CANDIDATE_SETS_METADATA_KEY)
                .cloned();
            let op_result = engine
                .apply(Operation::ScoreCandidateSetWeightedObjective {
                    set_name: set_name.clone(),
                    metric: metric.clone(),
                    objectives: objectives.clone(),
                    normalize_metrics: Some(*normalize_metrics),
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
                    "normalize_metrics": normalize_metrics,
                    "objective_count": objectives.len(),
                    "result": op_result
                }),
            }
        }
        ShellCommand::CandidatesTopK {
            input_set,
            output_set,
            metric,
            k,
            direction,
            tie_break,
        } => {
            let before = engine
                .state()
                .metadata
                .get(CANDIDATE_SETS_METADATA_KEY)
                .cloned();
            let op_result = engine
                .apply(Operation::TopKCandidateSet {
                    input_set: input_set.clone(),
                    output_set: output_set.clone(),
                    metric: metric.clone(),
                    k: *k,
                    direction: Some(*direction),
                    tie_break: Some(*tie_break),
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
                    "k": k,
                    "direction": direction.as_str(),
                    "tie_break": tie_break.as_str(),
                    "result": op_result
                }),
            }
        }
        ShellCommand::CandidatesParetoFrontier {
            input_set,
            output_set,
            objectives,
            max_candidates,
            tie_break,
        } => {
            let before = engine
                .state()
                .metadata
                .get(CANDIDATE_SETS_METADATA_KEY)
                .cloned();
            let op_result = engine
                .apply(Operation::ParetoFrontierCandidateSet {
                    input_set: input_set.clone(),
                    output_set: output_set.clone(),
                    objectives: objectives.clone(),
                    max_candidates: *max_candidates,
                    tie_break: Some(*tie_break),
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
                    "objective_count": objectives.len(),
                    "max_candidates": max_candidates,
                    "tie_break": tie_break.as_str(),
                    "result": op_result
                }),
            }
        }
        ShellCommand::CandidatesMacro {
            script,
            transactional,
        } => run_candidates_macro(engine, script, *transactional, options)?,
        ShellCommand::CandidatesTemplateList => {
            let templates = engine.list_candidate_macro_templates();
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": "gentle.candidate_macro_templates.v1",
                    "template_count": templates.len(),
                    "templates": templates
                }),
            }
        }
        ShellCommand::CandidatesTemplateShow { name } => {
            let template = engine
                .get_candidate_macro_template(name)
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "template": template
                }),
            }
        }
        ShellCommand::CandidatesTemplateUpsert {
            name,
            description,
            details_url,
            parameters,
            script,
        } => {
            let before = engine
                .state()
                .metadata
                .get(CANDIDATE_MACRO_TEMPLATES_METADATA_KEY)
                .cloned();
            let loaded_script = load_candidates_macro_script(script)?;
            let op_result = engine
                .apply(Operation::UpsertCandidateMacroTemplate {
                    name: name.clone(),
                    description: description.clone(),
                    details_url: details_url.clone(),
                    parameters: parameters.clone(),
                    script: loaded_script.clone(),
                })
                .map_err(|e| e.to_string())?;
            let after = engine
                .state()
                .metadata
                .get(CANDIDATE_MACRO_TEMPLATES_METADATA_KEY)
                .cloned();
            ShellRunResult {
                state_changed: before != after,
                output: json!({
                    "name": name,
                    "description": description,
                    "details_url": details_url,
                    "parameter_count": parameters.len(),
                    "script_length": loaded_script.len(),
                    "result": op_result
                }),
            }
        }
        ShellCommand::CandidatesTemplateDelete { name } => {
            let before = engine
                .state()
                .metadata
                .get(CANDIDATE_MACRO_TEMPLATES_METADATA_KEY)
                .cloned();
            let op_result = engine
                .apply(Operation::DeleteCandidateMacroTemplate { name: name.clone() })
                .map_err(|e| e.to_string())?;
            let after = engine
                .state()
                .metadata
                .get(CANDIDATE_MACRO_TEMPLATES_METADATA_KEY)
                .cloned();
            ShellRunResult {
                state_changed: before != after,
                output: json!({
                    "name": name,
                    "result": op_result
                }),
            }
        }
        ShellCommand::CandidatesTemplateRun {
            name,
            bindings,
            transactional,
        } => {
            let script = engine
                .render_candidate_macro_template_script(name, bindings)
                .map_err(|e| e.to_string())?;
            let mut run = run_candidates_macro(engine, &script, *transactional, options)?;
            run.output = json!({
                "template_name": name,
                "bindings": bindings,
                "expanded_script": script,
                "run": run.output
            });
            run
        }
        ShellCommand::GuidesList => {
            let sets = engine.list_guide_sets();
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": "gentle.guide_design.v1",
                    "set_count": sets.len(),
                    "sets": sets
                }),
            }
        }
        ShellCommand::GuidesShow {
            guide_set_id,
            limit,
            offset,
        } => {
            let (set, total, clamped_offset) = engine
                .inspect_guide_set_page(guide_set_id, *limit, *offset)
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "guide_set_id": set.guide_set_id,
                    "guide_count": total,
                    "limit": limit,
                    "offset": clamped_offset,
                    "returned": set.guides.len(),
                    "created_at_unix_ms": set.created_at_unix_ms,
                    "updated_at_unix_ms": set.updated_at_unix_ms,
                    "guides": set.guides,
                }),
            }
        }
        ShellCommand::GuidesPut {
            guide_set_id,
            guides_json,
        } => {
            let before = engine
                .state()
                .metadata
                .get(GUIDE_DESIGN_METADATA_KEY)
                .cloned();
            let raw = parse_json_payload(guides_json)?;
            let guides: Vec<GuideCandidate> = serde_json::from_str(&raw)
                .map_err(|e| format!("Invalid guides JSON payload: {e}"))?;
            let op_result = engine
                .apply(Operation::UpsertGuideSet {
                    guide_set_id: guide_set_id.clone(),
                    guides,
                })
                .map_err(|e| e.to_string())?;
            let after = engine
                .state()
                .metadata
                .get(GUIDE_DESIGN_METADATA_KEY)
                .cloned();
            ShellRunResult {
                state_changed: before != after,
                output: json!({
                    "guide_set_id": guide_set_id,
                    "result": op_result
                }),
            }
        }
        ShellCommand::GuidesDelete { guide_set_id } => {
            let before = engine
                .state()
                .metadata
                .get(GUIDE_DESIGN_METADATA_KEY)
                .cloned();
            let op_result = engine
                .apply(Operation::DeleteGuideSet {
                    guide_set_id: guide_set_id.clone(),
                })
                .map_err(|e| e.to_string())?;
            let after = engine
                .state()
                .metadata
                .get(GUIDE_DESIGN_METADATA_KEY)
                .cloned();
            ShellRunResult {
                state_changed: before != after,
                output: json!({
                    "guide_set_id": guide_set_id,
                    "result": op_result
                }),
            }
        }
        ShellCommand::GuidesFilter {
            guide_set_id,
            config_json,
            output_guide_set_id,
        } => {
            let before = engine
                .state()
                .metadata
                .get(GUIDE_DESIGN_METADATA_KEY)
                .cloned();
            let config = if let Some(raw) = config_json {
                let loaded = parse_json_payload(raw)?;
                serde_json::from_str::<GuidePracticalFilterConfig>(&loaded)
                    .map_err(|e| format!("Invalid practical guide filter config JSON: {e}"))?
            } else {
                GuidePracticalFilterConfig::default()
            };
            let op_result = engine
                .apply(Operation::FilterGuidesPractical {
                    guide_set_id: guide_set_id.clone(),
                    config,
                    output_guide_set_id: output_guide_set_id.clone(),
                })
                .map_err(|e| e.to_string())?;
            let after = engine
                .state()
                .metadata
                .get(GUIDE_DESIGN_METADATA_KEY)
                .cloned();
            ShellRunResult {
                state_changed: before != after,
                output: json!({
                    "guide_set_id": guide_set_id,
                    "output_guide_set_id": output_guide_set_id,
                    "result": op_result
                }),
            }
        }
        ShellCommand::GuidesFilterShow { guide_set_id } => {
            let report = engine
                .get_guide_practical_filter_report(guide_set_id)
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "report": report
                }),
            }
        }
        ShellCommand::GuidesOligosGenerate {
            guide_set_id,
            template_id,
            apply_5prime_g_extension,
            output_oligo_set_id,
            passed_only,
        } => {
            let before = engine
                .state()
                .metadata
                .get(GUIDE_DESIGN_METADATA_KEY)
                .cloned();
            let op_result = engine
                .apply(Operation::GenerateGuideOligos {
                    guide_set_id: guide_set_id.clone(),
                    template_id: template_id.clone(),
                    apply_5prime_g_extension: Some(*apply_5prime_g_extension),
                    output_oligo_set_id: output_oligo_set_id.clone(),
                    passed_only: Some(*passed_only),
                })
                .map_err(|e| e.to_string())?;
            let after = engine
                .state()
                .metadata
                .get(GUIDE_DESIGN_METADATA_KEY)
                .cloned();
            ShellRunResult {
                state_changed: before != after,
                output: json!({
                    "guide_set_id": guide_set_id,
                    "template_id": template_id,
                    "result": op_result
                }),
            }
        }
        ShellCommand::GuidesOligosList { guide_set_id } => {
            let sets = engine.list_guide_oligo_sets(guide_set_id.as_deref());
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": "gentle.guide_design.v1",
                    "oligo_set_count": sets.len(),
                    "oligo_sets": sets
                }),
            }
        }
        ShellCommand::GuidesOligosShow { oligo_set_id } => {
            let set = engine
                .get_guide_oligo_set(oligo_set_id)
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "oligo_set": set
                }),
            }
        }
        ShellCommand::GuidesOligosExport {
            guide_set_id,
            oligo_set_id,
            format,
            path,
            plate_format,
        } => {
            let before = engine
                .state()
                .metadata
                .get(GUIDE_DESIGN_METADATA_KEY)
                .cloned();
            let op_result = engine
                .apply(Operation::ExportGuideOligos {
                    guide_set_id: guide_set_id.clone(),
                    oligo_set_id: oligo_set_id.clone(),
                    format: *format,
                    path: path.clone(),
                    plate_format: *plate_format,
                })
                .map_err(|e| e.to_string())?;
            let after = engine
                .state()
                .metadata
                .get(GUIDE_DESIGN_METADATA_KEY)
                .cloned();
            ShellRunResult {
                state_changed: before != after,
                output: json!({
                    "guide_set_id": guide_set_id,
                    "oligo_set_id": oligo_set_id,
                    "format": format.as_str(),
                    "path": path,
                    "result": op_result
                }),
            }
        }
        ShellCommand::GuidesProtocolExport {
            guide_set_id,
            oligo_set_id,
            path,
            include_qc_checklist,
        } => {
            let before = engine
                .state()
                .metadata
                .get(GUIDE_DESIGN_METADATA_KEY)
                .cloned();
            let op_result = engine
                .apply(Operation::ExportGuideProtocolText {
                    guide_set_id: guide_set_id.clone(),
                    oligo_set_id: oligo_set_id.clone(),
                    path: path.clone(),
                    include_qc_checklist: Some(*include_qc_checklist),
                })
                .map_err(|e| e.to_string())?;
            let after = engine
                .state()
                .metadata
                .get(GUIDE_DESIGN_METADATA_KEY)
                .cloned();
            ShellRunResult {
                state_changed: before != after,
                output: json!({
                    "guide_set_id": guide_set_id,
                    "oligo_set_id": oligo_set_id,
                    "path": path,
                    "result": op_result
                }),
            }
        }
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
            let before_state = serde_json::to_value(engine.snapshot()).ok();
            let op_result = engine.apply(op).map_err(|e| e.to_string())?;
            let state_changed = if let Some(before) = before_state {
                serde_json::to_value(engine.snapshot())
                    .map(|after| after != before)
                    .unwrap_or_else(|_| {
                        !op_result.created_seq_ids.is_empty()
                            || !op_result.changed_seq_ids.is_empty()
                    })
            } else {
                !op_result.created_seq_ids.is_empty() || !op_result.changed_seq_ids.is_empty()
            };
            ShellRunResult {
                state_changed,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::Workflow { payload } => {
            let json_text = parse_json_payload(payload)?;
            let workflow: Workflow = serde_json::from_str(&json_text)
                .map_err(|e| format!("Invalid workflow JSON: {e}"))?;
            let before_state = serde_json::to_value(engine.snapshot()).ok();
            let results = engine.apply_workflow(workflow).map_err(|e| e.to_string())?;
            let state_changed = if let Some(before) = before_state {
                serde_json::to_value(engine.snapshot())
                    .map(|after| after != before)
                    .unwrap_or_else(|_| {
                        results
                            .iter()
                            .any(|r| !r.created_seq_ids.is_empty() || !r.changed_seq_ids.is_empty())
                    })
            } else {
                results
                    .iter()
                    .any(|r| !r.created_seq_ids.is_empty() || !r.changed_seq_ids.is_empty())
            };
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
    use std::sync::atomic::{AtomicUsize, Ordering};
    use tempfile::tempdir;

    static JASPAR_RELOAD_TEST_COUNTER: AtomicUsize = AtomicUsize::new(1);

    fn resource_fixture_path(name: &str) -> String {
        format!(
            "{}/test_files/fixtures/resources/{name}",
            env!("CARGO_MANIFEST_DIR")
        )
    }

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
                container_ids,
                arrangement_id,
            } => {
                assert_eq!(inputs, vec!["a".to_string(), "b".to_string()]);
                assert_eq!(output, "out.svg".to_string());
                assert_eq!(ladders, Some(vec!["1kb".to_string(), "100bp".to_string()]));
                assert_eq!(container_ids, None);
                assert_eq!(arrangement_id, None);
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_render_pool_gel_from_arrangement() {
        let cmd = parse_shell_line("render-pool-gel-svg - out.svg --arrangement arrangement-2")
            .expect("parse command");
        match cmd {
            ShellCommand::RenderPoolGelSvg {
                inputs,
                output,
                ladders,
                container_ids,
                arrangement_id,
            } => {
                assert!(inputs.is_empty());
                assert_eq!(output, "out.svg".to_string());
                assert_eq!(ladders, None);
                assert_eq!(container_ids, None);
                assert_eq!(arrangement_id, Some("arrangement-2".to_string()));
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_render_gel_svg_alias() {
        let cmd = parse_shell_line("render-gel-svg - out.svg --arrangement arrangement-2")
            .expect("parse command");
        match cmd {
            ShellCommand::RenderPoolGelSvg {
                inputs,
                output,
                ladders,
                container_ids,
                arrangement_id,
            } => {
                assert!(inputs.is_empty());
                assert_eq!(output, "out.svg".to_string());
                assert_eq!(ladders, None);
                assert_eq!(container_ids, None);
                assert_eq!(arrangement_id, Some("arrangement-2".to_string()));
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_arrange_serial_command() {
        let cmd = parse_shell_line(
            "arrange-serial container-1,container-2 --id arr-x --name test --ladders 100bp,1kb",
        )
        .expect("parse command");
        match cmd {
            ShellCommand::CreateArrangementSerial {
                container_ids,
                arrangement_id,
                name,
                ladders,
            } => {
                assert_eq!(
                    container_ids,
                    vec!["container-1".to_string(), "container-2".to_string()]
                );
                assert_eq!(arrangement_id, Some("arr-x".to_string()));
                assert_eq!(name, Some("test".to_string()));
                assert_eq!(ladders, Some(vec!["100bp".to_string(), "1kb".to_string()]));
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
    fn parse_routines_list_with_filters() {
        let cmd = parse_shell_line(
            "routines list --catalog assets/cloning_routines.json --family crispr --status implemented --tag guide --query scan",
        )
        .expect("parse routines list");
        match cmd {
            ShellCommand::RoutinesList {
                catalog_path,
                family,
                status,
                tag,
                query,
            } => {
                assert_eq!(
                    catalog_path.as_deref(),
                    Some("assets/cloning_routines.json")
                );
                assert_eq!(family.as_deref(), Some("crispr"));
                assert_eq!(status.as_deref(), Some("implemented"));
                assert_eq!(tag.as_deref(), Some("guide"));
                assert_eq!(query.as_deref(), Some("scan"));
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
    fn parse_candidates_score_weighted() {
        let cmd = parse_shell_line(
            "candidates score-weighted set1 objective --term gc_fraction:0.7:max --term distance_to_seq_start_bp:0.3:min --normalize",
        )
        .expect("parse candidates score-weighted");
        match cmd {
            ShellCommand::CandidatesScoreWeightedObjective {
                set_name,
                metric,
                objectives,
                normalize_metrics,
            } => {
                assert_eq!(set_name, "set1");
                assert_eq!(metric, "objective");
                assert_eq!(objectives.len(), 2);
                assert!(normalize_metrics);
                assert_eq!(objectives[0].metric, "gc_fraction");
                assert_eq!(
                    objectives[1].direction,
                    CandidateObjectiveDirection::Minimize
                );
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_candidates_top_k() {
        let cmd = parse_shell_line(
            "candidates top-k in_set out_set --metric objective --k 5 --direction max --tie-break length_descending",
        )
        .expect("parse candidates top-k");
        match cmd {
            ShellCommand::CandidatesTopK {
                input_set,
                output_set,
                metric,
                k,
                direction,
                tie_break,
            } => {
                assert_eq!(input_set, "in_set");
                assert_eq!(output_set, "out_set");
                assert_eq!(metric, "objective");
                assert_eq!(k, 5);
                assert_eq!(direction, CandidateObjectiveDirection::Maximize);
                assert_eq!(tie_break, CandidateTieBreakPolicy::LengthDescending);
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_candidates_pareto() {
        let cmd = parse_shell_line(
            "candidates pareto in_set out_set --objective gc_fraction:max --objective distance_to_seq_start_bp:min --max-candidates 10 --tie-break seq_start_end",
        )
        .expect("parse candidates pareto");
        match cmd {
            ShellCommand::CandidatesParetoFrontier {
                input_set,
                output_set,
                objectives,
                max_candidates,
                tie_break,
            } => {
                assert_eq!(input_set, "in_set");
                assert_eq!(output_set, "out_set");
                assert_eq!(objectives.len(), 2);
                assert_eq!(max_candidates, Some(10));
                assert_eq!(tie_break, CandidateTieBreakPolicy::SeqStartEnd);
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_candidates_template_put_and_run() {
        let put = parse_shell_line(
            "candidates template-put scan --script 'generate ${set_name} ${seq_id} --length ${len}' --details-url https://example.org/candidates/scan --param set_name --param seq_id=seqA --param len=20",
        )
        .expect("parse template-put");
        match put {
            ShellCommand::CandidatesTemplateUpsert {
                name,
                details_url,
                parameters,
                script,
                ..
            } => {
                assert_eq!(name, "scan");
                assert_eq!(
                    details_url.as_deref(),
                    Some("https://example.org/candidates/scan")
                );
                assert_eq!(parameters.len(), 3);
                assert_eq!(script, "generate ${set_name} ${seq_id} --length ${len}");
            }
            other => panic!("unexpected command: {other:?}"),
        }

        let run = parse_shell_line(
            "candidates template-run scan --bind set_name=hits --bind seq_id=seqB --transactional",
        )
        .expect("parse template-run");
        match run {
            ShellCommand::CandidatesTemplateRun {
                name,
                bindings,
                transactional,
            } => {
                assert_eq!(name, "scan");
                assert_eq!(bindings.get("set_name"), Some(&"hits".to_string()));
                assert_eq!(bindings.get("seq_id"), Some(&"seqB".to_string()));
                assert!(transactional);
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_candidates_template_put_accepts_url_alias() {
        let put = parse_shell_line(
            "candidates template-put scan --script 'generate ${set_name} ${seq_id} --length 20' --url https://example.org/candidates/scan --param set_name --param seq_id",
        )
        .expect("parse candidates template-put with --url");
        match put {
            ShellCommand::CandidatesTemplateUpsert { details_url, .. } => {
                assert_eq!(
                    details_url.as_deref(),
                    Some("https://example.org/candidates/scan")
                );
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_guides_filter_and_oligos_export_with_options() {
        let filter = parse_shell_line(
            "guides filter tp73 --config '{\"gc_min\":0.35,\"gc_max\":0.7}' --output-set tp73_pass",
        )
        .expect("parse guides filter");
        match filter {
            ShellCommand::GuidesFilter {
                guide_set_id,
                config_json,
                output_guide_set_id,
            } => {
                assert_eq!(guide_set_id, "tp73".to_string());
                assert!(config_json.is_some());
                assert_eq!(output_guide_set_id, Some("tp73_pass".to_string()));
            }
            other => panic!("unexpected command: {other:?}"),
        }

        let export = parse_shell_line(
            "guides oligos-export tp73 out.csv --format plate_csv --plate 96 --oligo-set tp73_set",
        )
        .expect("parse guides oligos-export");
        match export {
            ShellCommand::GuidesOligosExport {
                guide_set_id,
                oligo_set_id,
                format,
                path,
                plate_format,
            } => {
                assert_eq!(guide_set_id, "tp73");
                assert_eq!(oligo_set_id.as_deref(), Some("tp73_set"));
                assert_eq!(format, GuideOligoExportFormat::PlateCsv);
                assert_eq!(path, "out.csv");
                assert_eq!(plate_format, Some(GuideOligoPlateFormat::Plate96));
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_macros_run_and_template_commands() {
        let run =
            parse_shell_line("macros run --transactional --file test_files/workflow_plan.gsh")
                .expect("parse macros run");
        match run {
            ShellCommand::MacrosRun {
                script,
                transactional,
            } => {
                assert_eq!(script, "@test_files/workflow_plan.gsh");
                assert!(transactional);
            }
            other => panic!("unexpected command: {other:?}"),
        }

        let put = parse_shell_line(
            r#"macros template-put clone --script 'op {"Reverse":{"input":"${seq_id}","output_id":"${out_id}"}}' --details-url https://example.org/clone --param seq_id --param out_id=seqA_rev"#,
        )
        .expect("parse macros template-put");
        match put {
            ShellCommand::MacrosTemplateUpsert {
                name,
                details_url,
                parameters,
                script,
                ..
            } => {
                assert_eq!(name, "clone");
                assert_eq!(details_url.as_deref(), Some("https://example.org/clone"));
                assert_eq!(parameters.len(), 2);
                assert_eq!(
                    script,
                    r#"op {"Reverse":{"input":"${seq_id}","output_id":"${out_id}"}}"#
                );
            }
            other => panic!("unexpected command: {other:?}"),
        }

        let run_template = parse_shell_line(
            "macros template-run clone --bind seq_id=seqB --bind out_id=seqB_rev --transactional",
        )
        .expect("parse macros template-run");
        match run_template {
            ShellCommand::MacrosTemplateRun {
                name,
                bindings,
                transactional,
            } => {
                assert_eq!(name, "clone");
                assert_eq!(bindings.get("seq_id"), Some(&"seqB".to_string()));
                assert_eq!(bindings.get("out_id"), Some(&"seqB_rev".to_string()));
                assert!(transactional);
            }
            other => panic!("unexpected command: {other:?}"),
        }

        let import = parse_shell_line("macros template-import assets/cloning_patterns.json")
            .expect("parse macros template-import");
        match import {
            ShellCommand::MacrosTemplateImport { path } => {
                assert_eq!(path, "assets/cloning_patterns.json");
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_macros_template_put_accepts_url_alias() {
        let put = parse_shell_line(
            r#"macros template-put clone --script 'op {"Reverse":{"input":"${seq_id}"}}' --url https://example.org/clone --param seq_id"#,
        )
        .expect("parse macros template-put with --url");
        match put {
            ShellCommand::MacrosTemplateUpsert { details_url, .. } => {
                assert_eq!(details_url.as_deref(), Some("https://example.org/clone"));
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
        assert!(
            generated.output["result"]["messages"]
                .as_array()
                .map(|messages| !messages.is_empty())
                .unwrap_or(false)
        );

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
    fn execute_candidates_optimizer_primitives() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "seqA".to_string(),
            DNAsequence::from_sequence("GCATGAAA").expect("sequence"),
        );
        let mut engine = GentleEngine::from_state(state);

        execute_shell_command(
            &mut engine,
            &ShellCommand::CandidatesGenerate {
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
                limit: 64,
            },
        )
        .expect("generate candidate set");

        let weighted = execute_shell_command(
            &mut engine,
            &ShellCommand::CandidatesScoreWeightedObjective {
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
                normalize_metrics: true,
            },
        )
        .expect("weighted objective");
        assert!(weighted.state_changed);

        let topk = execute_shell_command(
            &mut engine,
            &ShellCommand::CandidatesTopK {
                input_set: "cand".to_string(),
                output_set: "top".to_string(),
                metric: "objective".to_string(),
                k: 1,
                direction: CandidateObjectiveDirection::Maximize,
                tie_break: CandidateTieBreakPolicy::SeqStartEnd,
            },
        )
        .expect("top-k");
        assert!(topk.state_changed);

        let pareto = execute_shell_command(
            &mut engine,
            &ShellCommand::CandidatesParetoFrontier {
                input_set: "cand".to_string(),
                output_set: "front".to_string(),
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
                tie_break: CandidateTieBreakPolicy::SeqStartEnd,
            },
        )
        .expect("pareto");
        assert!(pareto.state_changed);

        let top_set = engine
            .inspect_candidate_set_page("top", 10, 0)
            .expect("inspect top")
            .0;
        assert_eq!(top_set.candidates.len(), 1);
        assert_eq!(top_set.candidates[0].start_0based, 0);
    }

    #[test]
    fn execute_candidates_template_registry_and_run() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "seqA".to_string(),
            DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
        );
        let mut engine = GentleEngine::from_state(state);

        execute_shell_command(
            &mut engine,
            &ShellCommand::CandidatesTemplateUpsert {
                name: "scan".to_string(),
                description: Some("scan template".to_string()),
                details_url: Some("https://example.org/candidates/scan".to_string()),
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
                ],
                script: "generate ${set_name} ${seq_id} --length 4 --step 2".to_string(),
            },
        )
        .expect("upsert template");

        let listed = execute_shell_command(&mut engine, &ShellCommand::CandidatesTemplateList)
            .expect("list templates");
        assert_eq!(listed.output["template_count"].as_u64(), Some(1));
        let details_url = listed
            .output
            .get("templates")
            .and_then(|v| v.as_array())
            .and_then(|arr| arr.first())
            .and_then(|template| template.get("details_url"))
            .and_then(|v| v.as_str());
        assert_eq!(details_url, Some("https://example.org/candidates/scan"));

        let run = execute_shell_command(
            &mut engine,
            &ShellCommand::CandidatesTemplateRun {
                name: "scan".to_string(),
                bindings: HashMap::from([("set_name".to_string(), "hits".to_string())]),
                transactional: false,
            },
        )
        .expect("run template");
        assert!(run.state_changed);
        assert_eq!(run.output["template_name"].as_str(), Some("scan"));
        assert!(
            engine
                .list_candidate_sets()
                .iter()
                .any(|set| set.name == "hits")
        );
    }

    #[test]
    fn execute_routines_list_filters_and_searches() {
        let mut engine = GentleEngine::default();
        let tmp = tempdir().expect("tempdir");
        let catalog_path = tmp.path().join("cloning_routines.json");
        fs::write(
            &catalog_path,
            r#"{
  "schema": "gentle.cloning_routines.v1",
  "routines": [
    {
      "routine_id": "restriction.digest_basic",
      "title": "Restriction Digest Basic",
      "family": "restriction",
      "status": "implemented",
      "vocabulary_tags": ["restriction", "digest", "sticky"],
      "template_name": "digest_ligate_extract_sticky",
      "input_ports": [
        { "port_id": "seq_id", "kind": "sequence", "required": true, "cardinality": "one" }
      ],
      "output_ports": [
        { "port_id": "output_id", "kind": "sequence", "required": false, "cardinality": "one" }
      ]
    },
    {
      "routine_id": "crispr.guides.scan_basic",
      "title": "Guide Candidate Scan",
      "family": "crispr",
      "status": "planned",
      "vocabulary_tags": ["crispr", "guide", "scan"],
      "template_name": "grna_candidate_priority_scan",
      "input_ports": [
        { "port_id": "seq_id", "kind": "sequence", "required": true, "cardinality": "one" }
      ],
      "output_ports": [
        { "port_id": "set_name", "kind": "candidate_set", "required": true, "cardinality": "one" }
      ]
    }
  ]
}"#,
        )
        .expect("write routines catalog");

        let run = execute_shell_command(
            &mut engine,
            &ShellCommand::RoutinesList {
                catalog_path: Some(catalog_path.to_string_lossy().to_string()),
                family: Some("crispr".to_string()),
                status: Some("planned".to_string()),
                tag: Some("guide".to_string()),
                query: Some("scan".to_string()),
            },
        )
        .expect("list routines");
        assert!(!run.state_changed);
        assert_eq!(
            run.output["schema"].as_str(),
            Some(CLONING_ROUTINE_LIST_SCHEMA)
        );
        assert_eq!(run.output["routine_count"].as_u64(), Some(1));
        let rows = run
            .output
            .get("routines")
            .and_then(|value| value.as_array())
            .cloned()
            .unwrap_or_default();
        assert_eq!(rows.len(), 1);
        assert_eq!(
            rows[0].get("routine_id").and_then(|value| value.as_str()),
            Some("crispr.guides.scan_basic")
        );
    }

    #[test]
    fn execute_macros_template_registry_and_run() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "seqA".to_string(),
            DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
        );
        let mut engine = GentleEngine::from_state(state);

        execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateUpsert {
                name: "clone".to_string(),
                description: Some("reverse helper".to_string()),
                details_url: Some("https://example.org/macros/clone".to_string()),
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
                script: r#"op {"Reverse":{"input":"${seq_id}","output_id":"${out_id}"}}"#
                    .to_string(),
            },
        )
        .expect("upsert macros template");

        let listed = execute_shell_command(&mut engine, &ShellCommand::MacrosTemplateList)
            .expect("list macros templates");
        assert_eq!(listed.output["template_count"].as_u64(), Some(1));
        let details_url = listed
            .output
            .get("templates")
            .and_then(|v| v.as_array())
            .and_then(|arr| arr.first())
            .and_then(|template| template.get("details_url"))
            .and_then(|v| v.as_str());
        assert_eq!(details_url, Some("https://example.org/macros/clone"));

        let run = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateRun {
                name: "clone".to_string(),
                bindings: HashMap::from([("seq_id".to_string(), "seqA".to_string())]),
                transactional: false,
            },
        )
        .expect("run macros template");
        assert!(run.state_changed);
        assert_eq!(run.output["template_name"].as_str(), Some("clone"));
        assert!(engine.state().sequences.contains_key("seqA_rev"));
    }

    #[test]
    fn execute_macros_template_cloning_digest_ligation_extract_fixture() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "x".to_string(),
            DNAsequence::from_sequence("ATGGATCCGCATGGATCCGCATGGATCCGC").expect("sequence"),
        );
        let mut engine = GentleEngine::from_state(state);

        let put = parse_shell_line(
            "macros template-put clone_slice --file test_files/cloning_digest_ligation_extract.gsh --param seq_id=x --param digest_prefix=d --param ligation_prefix=lig --param extract_from=0 --param extract_to=1 --param output_id=slice",
        )
        .expect("parse macros template-put clone fixture");
        let upsert =
            execute_shell_command(&mut engine, &put).expect("upsert clone fixture template");
        assert!(upsert.state_changed);

        let run_cmd = parse_shell_line("macros template-run clone_slice --transactional")
            .expect("parse macros template-run clone fixture");
        let run = execute_shell_command(&mut engine, &run_cmd).expect("run clone fixture template");
        assert!(run.state_changed);
        assert_eq!(run.output["template_name"].as_str(), Some("clone_slice"));
        assert!(engine.state().sequences.contains_key("d_1"));
        assert!(engine.state().sequences.contains_key("lig_1"));
        let slice = engine
            .state()
            .sequences
            .get("slice")
            .expect("expected extracted slice output");
        assert_eq!(slice.len(), 1);
    }

    #[test]
    fn execute_macros_template_import_file_and_run() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "seqA".to_string(),
            DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
        );
        let mut engine = GentleEngine::from_state(state);

        let tmp = tempdir().expect("tempdir");
        let pattern_path = tmp.path().join("patterns.json");
        fs::write(
            &pattern_path,
            r#"{
  "schema": "gentle.cloning_patterns.v1",
  "templates": [
    {
      "name": "reverse_default",
      "description": "reverse one sequence",
      "details_url": "https://example.org/reverse-default",
      "parameters": [
        { "name": "seq_id", "default_value": "seqA", "required": false },
        { "name": "output_id", "default_value": "seqA_rev", "required": false }
      ],
      "script": "op {\"Reverse\":{\"input\":\"${seq_id}\",\"output_id\":\"${output_id}\"}}"
    }
  ]
}"#,
        )
        .expect("write patterns file");

        let import = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateImport {
                path: pattern_path.to_string_lossy().to_string(),
            },
        )
        .expect("import patterns");
        assert!(import.state_changed);
        assert_eq!(import.output["imported_count"].as_u64(), Some(1));
        let imported_template = engine
            .get_workflow_macro_template("reverse_default")
            .expect("get imported template");
        assert_eq!(
            imported_template.details_url.as_deref(),
            Some("https://example.org/reverse-default")
        );

        let run = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateRun {
                name: "reverse_default".to_string(),
                bindings: HashMap::new(),
                transactional: false,
            },
        )
        .expect("run imported template");
        assert!(run.state_changed);
        assert!(engine.state().sequences.contains_key("seqA_rev"));
    }

    #[test]
    fn execute_macros_template_import_single_template_schema_file_and_run() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "seqA".to_string(),
            DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
        );
        let mut engine = GentleEngine::from_state(state);

        let tmp = tempdir().expect("tempdir");
        let pattern_path = tmp.path().join("reverse_one.json");
        fs::write(
            &pattern_path,
            r#"{
  "schema": "gentle.cloning_pattern_template.v1",
  "name": "reverse_one",
  "description": "reverse one sequence",
  "details_url": "https://example.org/reverse-one",
  "parameters": [
    { "name": "seq_id", "default_value": "seqA", "required": false },
    { "name": "output_id", "default_value": "seqA_rev", "required": false }
  ],
  "script": "op {\"Reverse\":{\"input\":\"${seq_id}\",\"output_id\":\"${output_id}\"}}"
}"#,
        )
        .expect("write single-template file");

        let import = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateImport {
                path: pattern_path.to_string_lossy().to_string(),
            },
        )
        .expect("import single-template file");
        assert!(import.state_changed);
        assert_eq!(import.output["imported_count"].as_u64(), Some(1));
        assert_eq!(
            import
                .output
                .get("source_files")
                .and_then(|v| v.as_array())
                .map(|rows| rows.len()),
            Some(1)
        );

        let run = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateRun {
                name: "reverse_one".to_string(),
                bindings: HashMap::new(),
                transactional: false,
            },
        )
        .expect("run imported single-template");
        assert!(run.state_changed);
        assert!(engine.state().sequences.contains_key("seqA_rev"));
    }

    #[test]
    fn execute_macros_template_import_directory_recursive() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "seqA".to_string(),
            DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
        );
        let mut engine = GentleEngine::from_state(state);

        let tmp = tempdir().expect("tempdir");
        let root = tmp.path().join("catalog");
        fs::create_dir_all(root.join("a")).expect("mkdir a");
        fs::create_dir_all(root.join("b").join("nested")).expect("mkdir nested");
        fs::write(
            root.join("a").join("reverse_one.json"),
            r#"{
  "schema": "gentle.cloning_pattern_template.v1",
  "name": "reverse_one",
  "parameters": [
    { "name": "seq_id", "default_value": "seqA", "required": false },
    { "name": "output_id", "default_value": "seqA_rev", "required": false }
  ],
  "script": "op {\"Reverse\":{\"input\":\"${seq_id}\",\"output_id\":\"${output_id}\"}}"
}"#,
        )
        .expect("write reverse");
        fs::write(
            root.join("b").join("nested").join("complement_one.json"),
            r#"{
  "schema": "gentle.cloning_pattern_template.v1",
  "name": "complement_one",
  "parameters": [
    { "name": "seq_id", "default_value": "seqA", "required": false },
    { "name": "output_id", "default_value": "seqA_comp", "required": false }
  ],
  "script": "op {\"Complement\":{\"input\":\"${seq_id}\",\"output_id\":\"${output_id}\"}}"
}"#,
        )
        .expect("write complement");

        let import = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateImport {
                path: root.to_string_lossy().to_string(),
            },
        )
        .expect("import catalog directory");
        assert!(import.state_changed);
        assert_eq!(import.output["imported_count"].as_u64(), Some(2));
        assert_eq!(
            import
                .output
                .get("source_files")
                .and_then(|v| v.as_array())
                .map(|rows| rows.len()),
            Some(2)
        );
        let templates = import
            .output
            .get("templates")
            .and_then(|v| v.as_array())
            .cloned()
            .unwrap_or_default();
        assert!(templates.iter().any(|v| v.as_str() == Some("reverse_one")));
        assert!(
            templates
                .iter()
                .any(|v| v.as_str() == Some("complement_one"))
        );
    }

    #[test]
    fn execute_macros_template_import_builtin_patterns_and_run_template() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "seqA".to_string(),
            DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
        );
        let mut engine = GentleEngine::from_state(state);

        let path = format!(
            "{}/assets/cloning_patterns.json",
            env!("CARGO_MANIFEST_DIR")
        );
        let import =
            execute_shell_command(&mut engine, &ShellCommand::MacrosTemplateImport { path })
                .expect("import built-in patterns");
        assert!(import.state_changed);
        assert!(import.output["imported_count"].as_u64().unwrap_or(0) >= 6);
        let templates = import
            .output
            .get("templates")
            .and_then(|v| v.as_array())
            .cloned()
            .unwrap_or_default();
        assert!(templates.iter().any(|v| {
            v.as_str()
                .map(|s| s == "grna_candidate_priority_scan")
                .unwrap_or(false)
        }));
        let branch_template = engine
            .get_workflow_macro_template("branch_reverse_complement")
            .expect("imported branch template should exist");
        assert_eq!(
            branch_template.details_url.as_deref(),
            Some("https://www.bioinformatics.org/sms/rev_comp.html")
        );

        let mut bindings = HashMap::new();
        bindings.insert("seq_id".to_string(), "seqA".to_string());
        bindings.insert("branch_copy_id".to_string(), "seqA_branch".to_string());
        bindings.insert(
            "reverse_complement_id".to_string(),
            "seqA_branch_rc".to_string(),
        );
        let run = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateRun {
                name: "branch_reverse_complement".to_string(),
                bindings,
                transactional: false,
            },
        )
        .expect("run imported branch template");
        assert!(run.state_changed);
        assert!(engine.state().sequences.contains_key("seqA_branch"));
        assert!(engine.state().sequences.contains_key("seqA_branch_rc"));
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
    fn execute_guides_commands_end_to_end() {
        let mut engine = GentleEngine::from_state(ProjectState::default());
        let tmp = tempdir().expect("tempdir");
        let csv_path = tmp.path().join("guides.csv");
        let protocol_path = tmp.path().join("guides.protocol.txt");

        let guides_json = serde_json::to_string(&vec![
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
                start_0based: 220,
                end_0based_exclusive: 240,
                strand: "+".to_string(),
                protospacer: "TTTTGCCATGTTGACCTGAA".to_string(),
                pam: "TGG".to_string(),
                nuclease: "SpCas9".to_string(),
                cut_offset_from_protospacer_start: 17,
                rank: Some(2),
            },
        ])
        .expect("serialize guides");

        let put = execute_shell_command(
            &mut engine,
            &ShellCommand::GuidesPut {
                guide_set_id: "tp73_guides".to_string(),
                guides_json,
            },
        )
        .expect("guides put");
        assert!(put.state_changed);

        let filter = execute_shell_command(
            &mut engine,
            &ShellCommand::GuidesFilter {
                guide_set_id: "tp73_guides".to_string(),
                config_json: Some(
                    "{\"gc_min\":0.3,\"gc_max\":0.7,\"avoid_u6_terminator_tttt\":true}".to_string(),
                ),
                output_guide_set_id: Some("tp73_pass".to_string()),
            },
        )
        .expect("guides filter");
        assert!(filter.state_changed);

        let report = execute_shell_command(
            &mut engine,
            &ShellCommand::GuidesFilterShow {
                guide_set_id: "tp73_guides".to_string(),
            },
        )
        .expect("guides filter-show");
        assert!(!report.state_changed);
        let report_rows = report.output["report"]["results"]
            .as_array()
            .map(|rows| rows.len())
            .unwrap_or_default();
        assert_eq!(report_rows, 2, "expected filter report rows for all guides");

        let generated = execute_shell_command(
            &mut engine,
            &ShellCommand::GuidesOligosGenerate {
                guide_set_id: "tp73_guides".to_string(),
                template_id: "lenti_bsmbi_u6_default".to_string(),
                apply_5prime_g_extension: true,
                output_oligo_set_id: Some("tp73_oligos".to_string()),
                passed_only: true,
            },
        )
        .expect("guides oligos-generate");
        assert!(generated.state_changed);

        let listed = execute_shell_command(
            &mut engine,
            &ShellCommand::GuidesOligosList {
                guide_set_id: Some("tp73_guides".to_string()),
            },
        )
        .expect("guides oligos-list");
        assert!(!listed.state_changed);
        assert_eq!(listed.output["oligo_set_count"].as_u64(), Some(1));

        let exported = execute_shell_command(
            &mut engine,
            &ShellCommand::GuidesOligosExport {
                guide_set_id: "tp73_guides".to_string(),
                oligo_set_id: Some("tp73_oligos".to_string()),
                format: GuideOligoExportFormat::CsvTable,
                path: csv_path.to_string_lossy().to_string(),
                plate_format: None,
            },
        )
        .expect("guides oligos-export");
        assert!(exported.state_changed);
        let csv = fs::read_to_string(&csv_path).expect("read csv export");
        assert!(csv.contains("guide_id,rank,forward_oligo,reverse_oligo,notes"));

        let protocol = execute_shell_command(
            &mut engine,
            &ShellCommand::GuidesProtocolExport {
                guide_set_id: "tp73_guides".to_string(),
                oligo_set_id: Some("tp73_oligos".to_string()),
                path: protocol_path.to_string_lossy().to_string(),
                include_qc_checklist: true,
            },
        )
        .expect("guides protocol-export");
        assert!(protocol.state_changed);
        let protocol_text = fs::read_to_string(&protocol_path).expect("read protocol export");
        assert!(protocol_text.contains("GENtle Guide Oligo Protocol"));
    }

    #[test]
    fn execute_workflow_macro_transactional_rolls_back_on_error() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "seqA".to_string(),
            DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
        );
        let mut engine = GentleEngine::from_state(state);
        let err = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosRun {
                script: r#"op {"Reverse":{"input":"seqA","output_id":"tmp_rev"}}
op {"Reverse":{"input":"missing","output_id":"bad"}}"#
                    .to_string(),
                transactional: true,
            },
        )
        .unwrap_err();
        assert!(err.contains("rolled back"));
        assert!(!engine.state().sequences.contains_key("tmp_rev"));
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
    fn parse_agents_ask_command() {
        let cmd = parse_shell_line(
            "agents ask builtin_echo --prompt 'auto: state-summary' --catalog catalog.json --base-url http://localhost:11964 --model deepseek-r1:8b --timeout-secs 600 --connect-timeout-secs 20 --read-timeout-secs 900 --max-retries 4 --max-response-bytes 2097152 --allow-auto-exec --execute-index 2 --no-state-summary",
        )
        .expect("parse agents ask");
        match cmd {
            ShellCommand::AgentsAsk {
                system_id,
                prompt,
                catalog_path,
                base_url_override,
                model_override,
                timeout_seconds,
                connect_timeout_seconds,
                read_timeout_seconds,
                max_retries,
                max_response_bytes,
                include_state_summary,
                allow_auto_exec,
                execute_all,
                execute_indices,
            } => {
                assert_eq!(system_id, "builtin_echo");
                assert_eq!(prompt, "auto: state-summary");
                assert_eq!(catalog_path.as_deref(), Some("catalog.json"));
                assert_eq!(base_url_override.as_deref(), Some("http://localhost:11964"));
                assert_eq!(model_override.as_deref(), Some("deepseek-r1:8b"));
                assert_eq!(timeout_seconds, Some(600));
                assert_eq!(connect_timeout_seconds, Some(20));
                assert_eq!(read_timeout_seconds, Some(900));
                assert_eq!(max_retries, Some(4));
                assert_eq!(max_response_bytes, Some(2_097_152));
                assert!(!include_state_summary);
                assert!(allow_auto_exec);
                assert!(!execute_all);
                assert_eq!(execute_indices, vec![2]);
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn execute_agents_ask_runs_auto_suggestion_when_enabled() {
        let tmp = tempdir().expect("tempdir");
        let catalog_path = tmp.path().join("agents.json");
        let catalog_json = r#"{
  "schema": "gentle.agent_systems.v1",
  "systems": [
    {
      "id": "builtin_echo",
      "label": "Builtin Echo",
      "transport": "builtin_echo"
    }
  ]
}"#;
        fs::write(&catalog_path, catalog_json).expect("write catalog");

        let mut engine = GentleEngine::from_state(ProjectState::default());
        let out = execute_shell_command(
            &mut engine,
            &ShellCommand::AgentsAsk {
                system_id: "builtin_echo".to_string(),
                prompt: "auto: state-summary\nask: capabilities".to_string(),
                catalog_path: Some(catalog_path.display().to_string()),
                base_url_override: None,
                model_override: None,
                timeout_seconds: None,
                connect_timeout_seconds: None,
                read_timeout_seconds: None,
                max_retries: None,
                max_response_bytes: None,
                include_state_summary: true,
                allow_auto_exec: true,
                execute_all: false,
                execute_indices: vec![],
            },
        )
        .expect("execute agents ask");
        assert!(!out.state_changed);
        assert_eq!(
            out.output["summary"]["suggested_command_count"].as_u64(),
            Some(2)
        );
        assert_eq!(out.output["summary"]["executed_count"].as_u64(), Some(1));
    }

    #[test]
    fn execute_agents_ask_rejected_when_context_disallows_agent_commands() {
        let tmp = tempdir().expect("tempdir");
        let catalog_path = tmp.path().join("agents.json");
        let catalog_json = r#"{
  "schema": "gentle.agent_systems.v1",
  "systems": [
    {
      "id": "builtin_echo",
      "label": "Builtin Echo",
      "transport": "builtin_echo"
    }
  ]
}"#;
        fs::write(&catalog_path, catalog_json).expect("write catalog");
        let mut engine = GentleEngine::from_state(ProjectState::default());
        let err = execute_shell_command_with_options(
            &mut engine,
            &ShellCommand::AgentsAsk {
                system_id: "builtin_echo".to_string(),
                prompt: "ask: state-summary".to_string(),
                catalog_path: Some(catalog_path.display().to_string()),
                base_url_override: None,
                model_override: None,
                timeout_seconds: None,
                connect_timeout_seconds: None,
                read_timeout_seconds: None,
                max_retries: None,
                max_response_bytes: None,
                include_state_summary: true,
                allow_auto_exec: false,
                execute_all: false,
                execute_indices: vec![],
            },
            &ShellExecutionOptions {
                allow_screenshots: false,
                allow_agent_commands: false,
            },
        )
        .expect_err("agents ask should be blocked in this execution context");
        assert!(
            err.contains("agent-to-agent recursion guardrail"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn execute_agents_ask_blocks_nested_agent_call_inside_macro_suggestion() {
        let tmp = tempdir().expect("tempdir");
        let catalog_path = tmp.path().join("agents.json");
        let catalog_json = r#"{
  "schema": "gentle.agent_systems.v1",
  "systems": [
    {
      "id": "builtin_echo",
      "label": "Builtin Echo",
      "transport": "builtin_echo"
    }
  ]
}"#;
        fs::write(&catalog_path, catalog_json).expect("write catalog");

        let mut engine = GentleEngine::from_state(ProjectState::default());
        let out = execute_shell_command(
            &mut engine,
            &ShellCommand::AgentsAsk {
                system_id: "builtin_echo".to_string(),
                prompt: "auto: macros run 'agents ask builtin_echo --prompt \"ask: capabilities\"'"
                    .to_string(),
                catalog_path: Some(catalog_path.display().to_string()),
                base_url_override: None,
                model_override: None,
                timeout_seconds: None,
                connect_timeout_seconds: None,
                read_timeout_seconds: None,
                max_retries: None,
                max_response_bytes: None,
                include_state_summary: true,
                allow_auto_exec: true,
                execute_all: false,
                execute_indices: vec![],
            },
        )
        .expect("execute agents ask");
        assert!(!out.state_changed);
        assert_eq!(out.output["summary"]["executed_count"].as_u64(), Some(1));
        assert_eq!(
            out.output["summary"]["executed_error_count"].as_u64(),
            Some(1)
        );
        assert_eq!(out.output["executions"][0]["ok"].as_bool(), Some(false));
        let error = out.output["executions"][0]["error"]
            .as_str()
            .expect("error string");
        assert!(
            error.contains("agent-to-agent recursion guardrail"),
            "unexpected error: {error}"
        );
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
    fn execute_resources_sync_rebase_with_local_fixture() {
        let td = tempdir().expect("tempdir");
        let output_path = td.path().join("rebase.sync.json");
        let mut engine = GentleEngine::from_state(ProjectState::default());

        let out = execute_shell_command(
            &mut engine,
            &ShellCommand::ResourcesSyncRebase {
                input: resource_fixture_path("rebase.edge.withrefm"),
                output: Some(output_path.to_string_lossy().to_string()),
                commercial_only: true,
            },
        )
        .expect("execute resources sync-rebase");

        assert!(!out.state_changed);
        assert_eq!(
            out.output["report"]["resource"].as_str(),
            Some("rebase-commercial")
        );
        assert_eq!(out.output["report"]["item_count"].as_u64(), Some(2));
        assert_eq!(
            out.output["report"]["output"].as_str(),
            Some(output_path.to_string_lossy().as_ref())
        );

        let written = fs::read_to_string(&output_path).expect("read synced REBASE output");
        let enzymes = serde_json::from_str::<serde_json::Value>(&written).expect("parse JSON");
        let names = enzymes
            .as_array()
            .expect("enzyme array")
            .iter()
            .filter_map(|entry| entry.get("name").and_then(|v| v.as_str()))
            .collect::<Vec<_>>();
        assert_eq!(names, vec!["BsaI", "EcoRI"]);
    }

    #[test]
    fn execute_resources_sync_jaspar_with_local_fixture() {
        let td = tempdir().expect("tempdir");
        let output_path = td.path().join("jaspar.sync.json");
        let mut engine = GentleEngine::from_state(ProjectState::default());

        let out = execute_shell_command(
            &mut engine,
            &ShellCommand::ResourcesSyncJaspar {
                input: resource_fixture_path("jaspar.edge.pfm"),
                output: Some(output_path.to_string_lossy().to_string()),
            },
        )
        .expect("execute resources sync-jaspar");

        assert!(!out.state_changed);
        assert_eq!(out.output["report"]["resource"].as_str(), Some("jaspar"));
        assert_eq!(out.output["report"]["item_count"].as_u64(), Some(2));
        assert_eq!(
            out.output["report"]["output"].as_str(),
            Some(output_path.to_string_lossy().as_ref())
        );

        let written = fs::read_to_string(&output_path).expect("read synced JASPAR output");
        let snapshot = serde_json::from_str::<serde_json::Value>(&written).expect("parse JSON");
        assert_eq!(
            snapshot.get("schema").and_then(|v| v.as_str()),
            Some("gentle.tf_motifs.v2")
        );
        assert_eq!(
            snapshot.get("motif_count").and_then(|v| v.as_u64()),
            Some(2)
        );
    }

    #[test]
    fn execute_resources_sync_jaspar_reloads_motif_registry_from_synced_output() {
        struct ReloadResetGuard;
        impl Drop for ReloadResetGuard {
            fn drop(&mut self) {
                crate::tf_motifs::reload();
            }
        }
        let _guard = ReloadResetGuard;

        let td = tempdir().expect("tempdir");
        let unique = JASPAR_RELOAD_TEST_COUNTER.fetch_add(1, Ordering::Relaxed);
        let motif_id = format!("MTEST{unique}.1");
        let motif_name = format!("CODEX_TEST_MOTIF_{unique}");
        let input_path = td.path().join("custom_reload.pfm");
        let output_path = td.path().join("custom_reload.motifs.json");
        let motif_text = format!(
            ">{motif_id} {motif_name}\nA [1 0 0 0]\nC [0 1 0 0]\nG [0 0 1 0]\nT [0 0 0 1]\n"
        );
        fs::write(&input_path, motif_text).expect("write custom JASPAR input");

        assert_eq!(crate::tf_motifs::resolve_motif(&motif_name), None);

        let mut engine = GentleEngine::from_state(ProjectState::default());
        let out = execute_shell_command(
            &mut engine,
            &ShellCommand::ResourcesSyncJaspar {
                input: input_path.to_string_lossy().to_string(),
                output: Some(output_path.to_string_lossy().to_string()),
            },
        )
        .expect("execute resources sync-jaspar");

        assert!(!out.state_changed);
        assert_eq!(out.output["report"]["item_count"].as_u64(), Some(1));
        assert_eq!(
            crate::tf_motifs::resolve_motif(&motif_name).as_deref(),
            Some("ACGT")
        );
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
        assert!(
            created
                .iter()
                .any(|v| v.as_str().map(|id| id == "slice_ext5").unwrap_or(false))
        );
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
    fn execute_genomes_status_reports_length_and_mass_metadata() {
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
    "cache_dir": "{}",
    "nucleotide_length_bp": 4
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
            &ShellCommand::ReferenceStatus {
                helper_mode: false,
                genome_id: "ToyGenome".to_string(),
                catalog_path: Some(catalog.to_string_lossy().to_string()),
                cache_dir: None,
            },
        )
        .expect("execute status");
        assert!(!out.state_changed);
        assert_eq!(out.output["nucleotide_length_bp"].as_u64(), Some(4));
        let expected_mass = 4.0 * 617.96 + 36.04;
        let observed_mass = out.output["molecular_mass_da"]
            .as_f64()
            .expect("molecular_mass_da should be present");
        assert!((observed_mass - expected_mass).abs() < 1e-6);
        assert_eq!(
            out.output["molecular_mass_source"].as_str(),
            Some("estimated_from_nucleotide_length")
        );
    }

    #[test]
    fn parse_ui_open_and_prepared_commands() {
        let open = parse_shell_line("ui open prepared-references --genome-id \"Human GRCh38\"")
            .expect("parse ui open");
        match open {
            ShellCommand::UiIntent {
                action,
                target,
                genome_id,
                helper_mode,
                catalog_path,
                cache_dir,
                filter,
                species,
                latest,
            } => {
                assert_eq!(action, UiIntentAction::Open);
                assert_eq!(target, UiIntentTarget::PreparedReferences);
                assert_eq!(genome_id.as_deref(), Some("Human GRCh38"));
                assert!(!helper_mode);
                assert!(catalog_path.is_none());
                assert!(cache_dir.is_none());
                assert!(filter.is_none());
                assert!(species.is_none());
                assert!(!latest);
            }
            other => panic!("unexpected command: {other:?}"),
        }

        let open_latest = parse_shell_line("ui open prepared-references --species human --latest")
            .expect("parse ui open latest");
        match open_latest {
            ShellCommand::UiIntent {
                action,
                target,
                genome_id,
                helper_mode,
                species,
                latest,
                ..
            } => {
                assert_eq!(action, UiIntentAction::Open);
                assert_eq!(target, UiIntentTarget::PreparedReferences);
                assert!(genome_id.is_none());
                assert!(!helper_mode);
                assert_eq!(species.as_deref(), Some("human"));
                assert!(latest);
            }
            other => panic!("unexpected command: {other:?}"),
        }

        let prepared = parse_shell_line("ui prepared-genomes --species human --latest")
            .expect("parse ui prepared-genomes");
        match prepared {
            ShellCommand::UiPreparedGenomes {
                helper_mode,
                species,
                latest,
                ..
            } => {
                assert!(!helper_mode);
                assert_eq!(species.as_deref(), Some("human"));
                assert!(latest);
            }
            other => panic!("unexpected command: {other:?}"),
        }

        let err = parse_shell_line("ui open agent-assistant --latest")
            .expect_err("ui open agent-assistant --latest should fail");
        assert!(
            err.contains("only supports --helpers/--catalog/--cache-dir/--filter/--species/--latest when TARGET is prepared-references"),
            "unexpected parse error: {err}"
        );
    }

    #[test]
    fn execute_ui_prepared_and_latest_prepared_queries() {
        let td = tempdir().expect("tempdir");
        let fasta_113 = td.path().join("h113.fa");
        let ann_113 = td.path().join("h113.gtf");
        let fasta_116 = td.path().join("h116.fa");
        let ann_116 = td.path().join("h116.gtf");
        let cache_dir = td.path().join("cache");
        fs::write(&fasta_113, ">chr1\nACGT\n").expect("write fasta 113");
        fs::write(
            &ann_113,
            "chr1\tsrc\tgene\t1\t4\t.\t+\t.\tgene_id \"G113\"; gene_name \"G113\";\n",
        )
        .expect("write ann 113");
        fs::write(&fasta_116, ">chr1\nACGTACGT\n").expect("write fasta 116");
        fs::write(
            &ann_116,
            "chr1\tsrc\tgene\t1\t8\t.\t+\t.\tgene_id \"G116\"; gene_name \"G116\";\n",
        )
        .expect("write ann 116");
        let catalog_path = td.path().join("catalog.json");
        let catalog_json = format!(
            r#"{{
  "Human GRCh38 Ensembl 113": {{
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }},
  "Human GRCh38 Ensembl 116": {{
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
            fasta_113.display(),
            ann_113.display(),
            cache_dir.display(),
            fasta_116.display(),
            ann_116.display(),
            cache_dir.display()
        );
        fs::write(&catalog_path, catalog_json).expect("write catalog");

        let catalog =
            GenomeCatalog::from_json_file(catalog_path.to_string_lossy().as_ref()).unwrap();
        catalog
            .prepare_genome_once("Human GRCh38 Ensembl 113")
            .expect("prepare 113");
        std::thread::sleep(std::time::Duration::from_millis(2));
        catalog
            .prepare_genome_once("Human GRCh38 Ensembl 116")
            .expect("prepare 116");

        let mut engine = GentleEngine::new();
        let prepared = execute_shell_command(
            &mut engine,
            &ShellCommand::UiPreparedGenomes {
                helper_mode: false,
                catalog_path: Some(catalog_path.to_string_lossy().to_string()),
                cache_dir: None,
                filter: None,
                species: Some("human".to_string()),
                latest: false,
            },
        )
        .expect("execute ui prepared-genomes");
        assert!(!prepared.state_changed);
        assert_eq!(prepared.output["prepared_count"].as_u64(), Some(2));
        assert_eq!(
            prepared.output["genomes"][0]["genome_id"].as_str(),
            Some("Human GRCh38 Ensembl 116")
        );

        let latest = execute_shell_command(
            &mut engine,
            &ShellCommand::UiLatestPrepared {
                helper_mode: false,
                catalog_path: Some(catalog_path.to_string_lossy().to_string()),
                cache_dir: None,
                species: "human".to_string(),
            },
        )
        .expect("execute ui latest-prepared");
        assert!(!latest.state_changed);
        assert_eq!(
            latest.output["selected_genome_id"].as_str(),
            Some("Human GRCh38 Ensembl 116")
        );

        let intent_latest = execute_shell_command(
            &mut engine,
            &ShellCommand::UiIntent {
                action: UiIntentAction::Open,
                target: UiIntentTarget::PreparedReferences,
                genome_id: None,
                helper_mode: false,
                catalog_path: Some(catalog_path.to_string_lossy().to_string()),
                cache_dir: None,
                filter: None,
                species: Some("human".to_string()),
                latest: true,
            },
        )
        .expect("execute ui intent latest");
        assert!(!intent_latest.state_changed);
        assert_eq!(
            intent_latest.output["selected_genome_id"].as_str(),
            Some("Human GRCh38 Ensembl 116")
        );

        let intent_explicit = execute_shell_command(
            &mut engine,
            &ShellCommand::UiIntent {
                action: UiIntentAction::Open,
                target: UiIntentTarget::PreparedReferences,
                genome_id: Some("Human GRCh38 Ensembl 113".to_string()),
                helper_mode: false,
                catalog_path: Some(catalog_path.to_string_lossy().to_string()),
                cache_dir: None,
                filter: None,
                species: Some("human".to_string()),
                latest: true,
            },
        )
        .expect("execute ui intent explicit genome");
        assert!(!intent_explicit.state_changed);
        assert_eq!(
            intent_explicit.output["selected_genome_id"].as_str(),
            Some("Human GRCh38 Ensembl 113")
        );
        assert!(
            intent_explicit.output["prepared_query"].is_null(),
            "explicit genome_id should bypass prepared query resolution"
        );
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
    fn parse_feature_expert_commands() {
        let inspect = parse_shell_line("inspect-feature-expert s tfbs 7")
            .expect("parse inspect-feature-expert");
        match inspect {
            ShellCommand::InspectFeatureExpert { seq_id, target } => {
                assert_eq!(seq_id, "s");
                assert_eq!(target, FeatureExpertTarget::TfbsFeature { feature_id: 7 });
            }
            other => panic!("unexpected command: {other:?}"),
        }

        let splicing =
            parse_shell_line("inspect-feature-expert s splicing 11").expect("parse splicing");
        match splicing {
            ShellCommand::InspectFeatureExpert { seq_id, target } => {
                assert_eq!(seq_id, "s");
                assert_eq!(
                    target,
                    FeatureExpertTarget::SplicingFeature { feature_id: 11 }
                );
            }
            other => panic!("unexpected command: {other:?}"),
        }

        let render = parse_shell_line(
            "render-feature-expert-svg s restriction 123 --enzyme EcoRI --start 100 --end 106 out.svg",
        )
        .expect("parse render-feature-expert-svg");
        match render {
            ShellCommand::RenderFeatureExpertSvg {
                seq_id,
                target,
                output,
            } => {
                assert_eq!(seq_id, "s");
                assert_eq!(output, "out.svg");
                assert_eq!(
                    target,
                    FeatureExpertTarget::RestrictionSite {
                        cut_pos_1based: 123,
                        enzyme: Some("EcoRI".to_string()),
                        recognition_start_1based: Some(100),
                        recognition_end_1based: Some(106),
                    }
                );
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn execute_inspect_feature_expert_returns_view_json() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "s".to_string(),
            DNAsequence::from_sequence("TTTACGTAAACGTGGG").expect("valid dna"),
        );
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
            .expect("annotate tfbs");
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

        let output = execute_shell_command(
            &mut engine,
            &ShellCommand::InspectFeatureExpert {
                seq_id: "s".to_string(),
                target: FeatureExpertTarget::TfbsFeature { feature_id },
            },
        )
        .expect("execute inspect-feature-expert");
        assert!(!output.state_changed);
        assert_eq!(output.output["kind"].as_str(), Some("tfbs"));
        assert_eq!(
            output.output["data"]["feature_id"].as_u64(),
            Some(feature_id as u64)
        );
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

    #[test]
    fn execute_set_param_updates_tfbs_display_state() {
        let mut engine = GentleEngine::new();
        let out = execute_shell_command(
            &mut engine,
            &ShellCommand::SetParameter {
                name: "tfbs_display_min_llr_quantile".to_string(),
                value_json: "0.85".to_string(),
            },
        )
        .expect("execute set-param");
        assert!(out.state_changed);
        assert!(
            (engine.state().display.tfbs_display_min_llr_quantile - 0.85).abs() < f64::EPSILON,
            "tfbs_display_min_llr_quantile should be updated by set-param"
        );
    }

    #[test]
    fn execute_op_set_display_visibility_marks_state_changed() {
        let mut engine = GentleEngine::new();
        let out = execute_shell_command(
            &mut engine,
            &ShellCommand::Op {
                payload: r#"{"SetDisplayVisibility":{"target":"Tfbs","visible":true}}"#.to_string(),
            },
        )
        .expect("execute op");
        assert!(out.state_changed);
        assert!(engine.state().display.show_tfbs);
    }
}
