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
        CandidateWeightedObjectiveTerm, DOTPLOT_ANALYSIS_METADATA_KEY, DotplotMode, Engine,
        FeatureExpertTarget, FeatureExpertView, FlexibilityModel, GUIDE_DESIGN_METADATA_KEY,
        GenomeAnchorSide, GenomeAnnotationScope, GenomeTrackSource, GenomeTrackSubscription,
        GentleEngine, GuideCandidate, GuideOligoExportFormat, GuideOligoPlateFormat,
        GuidePracticalFilterConfig, LineageMacroInstance, LineageMacroPortBinding,
        MacroInstanceStatus, Operation, PLANNING_ESTIMATE_SCHEMA, PLANNING_OBJECTIVE_SCHEMA,
        PLANNING_PROFILE_SCHEMA, PLANNING_SUGGESTION_SCHEMA, PLANNING_SYNC_STATUS_SCHEMA,
        PRIMER_DESIGN_REPORTS_METADATA_KEY, PlanningEstimate, PlanningObjective, PlanningProfile,
        PlanningProfileScope, PlanningSuggestionStatus, PrimerDesignBackend,
        PrimerDesignPairConstraint, PrimerDesignSideConstraint, ProjectState, RenderSvgMode,
        RnaReadAlignConfig, RnaReadHitSelection, RnaReadInputFormat, RnaReadInterpretationProfile,
        RnaReadOriginMode, RnaReadReportMode, RnaReadScoreDensityScale,
        RnaReadSeedFilterConfig, SequenceAnchor, SplicingScopePreset,
        WORKFLOW_MACRO_TEMPLATES_METADATA_KEY, Workflow, WorkflowMacroTemplate,
        WorkflowMacroTemplateParam, WorkflowMacroTemplatePort,
    },
    enzymes::active_restriction_enzymes,
    feature_location::collect_location_ranges_usize,
    genomes::{
        DEFAULT_GENOME_CATALOG_PATH, DEFAULT_HELPER_GENOME_CATALOG_PATH, GenomeBlastReport,
        GenomeCatalog, GenomeGeneRecord,
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
use lazy_static::lazy_static;
#[cfg(all(target_os = "macos", feature = "screenshot-capture"))]
use objc2_app_kit::NSApplication;
#[cfg(all(target_os = "macos", feature = "screenshot-capture"))]
use objc2_foundation::MainThreadMarker;
use regex::{Regex, RegexBuilder};
use serde::{Deserialize, Serialize, de::DeserializeOwned};
use serde_json::{Value, json};
#[cfg(all(target_os = "macos", feature = "screenshot-capture"))]
use std::process::Command;
#[cfg(test)]
use std::sync::atomic::AtomicUsize;
#[cfg(test)]
use std::time::Duration;
use std::{
    collections::{BTreeSet, HashMap, HashSet},
    fs,
    path::{Path, PathBuf},
    sync::{
        Arc, Mutex,
        atomic::{AtomicBool, AtomicU64, Ordering},
        mpsc,
    },
    thread,
    time::{SystemTime, UNIX_EPOCH},
};

const CLONING_PATTERN_FILE_SCHEMA: &str = "gentle.cloning_patterns.v1";
const CLONING_PATTERN_TEMPLATE_FILE_SCHEMA: &str = "gentle.cloning_pattern_template.v1";
const CLONING_ROUTINE_CATALOG_SCHEMA: &str = "gentle.cloning_routines.v1";
const CLONING_ROUTINE_LIST_SCHEMA: &str = "gentle.cloning_routines_list.v1";
const CLONING_ROUTINE_EXPLAIN_SCHEMA: &str = "gentle.cloning_routine_explain.v1";
const CLONING_ROUTINE_COMPARE_SCHEMA: &str = "gentle.cloning_routine_compare.v1";
const DEFAULT_CLONING_ROUTINE_CATALOG_PATH: &str = "assets/cloning_routines.json";
const BLAST_ASYNC_JOB_SCHEMA: &str = "gentle.blast_async_job_status.v1";
const BLAST_ASYNC_JOB_HISTORY_LIMIT: usize = 200;
const BLAST_ASYNC_MAX_CONCURRENT_ENV: &str = "GENTLE_BLAST_ASYNC_MAX_CONCURRENT";
const BLAST_ASYNC_MAX_CONCURRENT_HARD_LIMIT: usize = 256;
static BLAST_ASYNC_JOB_COUNTER: AtomicU64 = AtomicU64::new(1);

#[derive(Debug)]
enum BlastAsyncWorkerMessage {
    Done(Result<GenomeBlastReport, String>),
}

#[derive(Debug)]
struct BlastAsyncJobRecord {
    status: BlastAsyncJobStatus,
    cancel_requested: Arc<AtomicBool>,
    receiver: Option<mpsc::Receiver<BlastAsyncWorkerMessage>>,
    launch_spec: Option<BlastAsyncLaunchSpec>,
    report: Option<GenomeBlastReport>,
}

#[derive(Debug, Clone)]
struct BlastAsyncLaunchSpec {
    engine_snapshot: GentleEngine,
    helper_mode: bool,
    genome_id: String,
    query_sequence: String,
    max_hits: usize,
    max_hits_explicit: bool,
    task: Option<String>,
    request_options_json: Option<Value>,
    resolved_catalog: Option<String>,
    cache_dir: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
struct BlastAsyncJobStatus {
    schema: String,
    job_id: String,
    state: String,
    helper_mode: bool,
    genome_id: String,
    query_length: usize,
    task: String,
    max_hits: usize,
    created_at_unix_ms: u128,
    started_at_unix_ms: Option<u128>,
    finished_at_unix_ms: Option<u128>,
    cancel_requested: bool,
    done_queries: usize,
    total_queries: usize,
    result_available: bool,
    error: Option<String>,
    max_concurrent_jobs: usize,
    running_jobs: usize,
    queued_jobs: usize,
    queue_position: Option<usize>,
}

lazy_static! {
    static ref BLAST_ASYNC_JOBS: Mutex<HashMap<String, BlastAsyncJobRecord>> =
        Mutex::new(HashMap::new());
}

#[cfg(test)]
static BLAST_ASYNC_MAX_CONCURRENT_TEST_OVERRIDE: AtomicUsize = AtomicUsize::new(0);
#[cfg(test)]
static BLAST_ASYNC_WORKER_DELAY_MS_TEST_OVERRIDE: AtomicU64 = AtomicU64::new(0);

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
    input_ports: Vec<WorkflowMacroTemplatePort>,
    output_ports: Vec<WorkflowMacroTemplatePort>,
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
    input_ports: Vec<WorkflowMacroTemplatePort>,
    output_ports: Vec<WorkflowMacroTemplatePort>,
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
            input_ports: vec![],
            output_ports: vec![],
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
struct CloningRoutineDifferenceAxis {
    axis: String,
    value: String,
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
    purpose: Option<String>,
    mechanism: Option<String>,
    requires: Vec<String>,
    contraindications: Vec<String>,
    confusing_alternatives: Vec<String>,
    difference_matrix: Vec<CloningRoutineDifferenceAxis>,
    disambiguation_questions: Vec<String>,
    failure_modes: Vec<String>,
    details_url: Option<String>,
    template_name: String,
    template_path: Option<String>,
    input_ports: Vec<CloningRoutinePort>,
    output_ports: Vec<CloningRoutinePort>,
    base_time_hours: Option<f64>,
    base_cost: Option<f64>,
    required_material_classes: Vec<String>,
    required_machine_classes: Vec<String>,
    required_capabilities: Vec<String>,
}

#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "snake_case")]
/// Validation status for one typed preflight port check row.
enum RoutinePortValidationStatus {
    Ok,
    Missing,
    Invalid,
    Skipped,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
/// Port direction in typed routine/template contracts.
enum RoutinePortDirection {
    Input,
    Output,
}

#[derive(Debug, Clone, Serialize)]
/// Deterministic validation record for one routine/template port binding.
///
/// This row is emitted in `gentle.macro_template_preflight.v1` reports so GUI,
/// CLI, and automation adapters can render equivalent diagnostics.
struct RoutinePortValidationRow {
    direction: RoutinePortDirection,
    source: String,
    port_id: String,
    kind: String,
    required: bool,
    cardinality: String,
    values: Vec<String>,
    status: RoutinePortValidationStatus,
    message: Option<String>,
}

#[derive(Debug, Clone, Serialize)]
/// Preflight report emitted by `macros template-run` (including `--validate-only`).
///
/// This report is intentionally adapter-neutral and machine-readable. It carries
/// typed port diagnostics, warnings, and hard errors that determine whether a
/// macro template can execute safely.
struct MacroTemplatePreflightReport {
    schema: String,
    template_name: String,
    routine_id: Option<String>,
    routine_title: Option<String>,
    routine_family: Option<String>,
    routine_status: Option<String>,
    catalog_path: Option<String>,
    contract_source: String,
    checked_ports: Vec<RoutinePortValidationRow>,
    warnings: Vec<String>,
    errors: Vec<String>,
}

impl MacroTemplatePreflightReport {
    /// Returns true when no hard preflight errors are present.
    fn can_execute(&self) -> bool {
        self.errors.is_empty()
    }
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
    PanelsImportIsoform {
        seq_id: String,
        panel_path: String,
        panel_id: Option<String>,
        strict: bool,
    },
    PanelsInspectIsoform {
        seq_id: String,
        panel_id: String,
    },
    PanelsRenderIsoformSvg {
        seq_id: String,
        panel_id: String,
        output: String,
    },
    PanelsValidateIsoform {
        panel_path: String,
        panel_id: Option<String>,
    },
    UniprotFetch {
        query: String,
        entry_id: Option<String>,
    },
    GenbankFetch {
        accession: String,
        as_id: Option<String>,
    },
    UniprotImportSwissProt {
        path: String,
        entry_id: Option<String>,
    },
    UniprotList,
    UniprotShow {
        entry_id: String,
    },
    UniprotMap {
        entry_id: String,
        seq_id: String,
        projection_id: Option<String>,
        transcript_id: Option<String>,
    },
    UniprotProjectionList {
        seq_id: Option<String>,
    },
    UniprotProjectionShow {
        projection_id: String,
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
    ExportRunBundle {
        output: String,
        run_id: Option<String>,
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
    RoutinesExplain {
        catalog_path: Option<String>,
        routine_id: String,
    },
    RoutinesCompare {
        catalog_path: Option<String>,
        left_routine_id: String,
        right_routine_id: String,
    },
    PlanningProfileShow {
        scope: PlanningProfileScope,
    },
    PlanningProfileSet {
        scope: PlanningProfileScope,
        payload_json: String,
    },
    PlanningObjectiveShow,
    PlanningObjectiveSet {
        payload_json: String,
    },
    PlanningSuggestionsList {
        status: Option<PlanningSuggestionStatus>,
    },
    PlanningSuggestionAccept {
        suggestion_id: String,
    },
    PlanningSuggestionReject {
        suggestion_id: String,
        reason: Option<String>,
    },
    PlanningSyncStatus,
    PlanningSyncPull {
        payload_json: String,
        source: Option<String>,
        confidence: Option<f64>,
        snapshot_id: Option<String>,
    },
    PlanningSyncPush {
        payload_json: String,
        source: Option<String>,
        confidence: Option<f64>,
        snapshot_id: Option<String>,
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
        annotation_scope: Option<GenomeAnnotationScope>,
        max_annotation_features: Option<usize>,
        include_genomic_annotation: Option<bool>,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
    },
    ReferenceExtractGene {
        helper_mode: bool,
        genome_id: String,
        gene_query: String,
        occurrence: Option<usize>,
        output_id: Option<String>,
        annotation_scope: Option<GenomeAnnotationScope>,
        max_annotation_features: Option<usize>,
        include_genomic_annotation: Option<bool>,
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
        prepared_genome_id: Option<String>,
    },
    ReferenceVerifyAnchor {
        helper_mode: bool,
        seq_id: String,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
        prepared_genome_id: Option<String>,
    },
    ReferenceBlast {
        helper_mode: bool,
        genome_id: String,
        query_sequence: String,
        max_hits: usize,
        max_hits_explicit: bool,
        task: Option<String>,
        request_options_json: Option<Value>,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
    },
    ReferenceBlastAsyncStart {
        helper_mode: bool,
        genome_id: String,
        query_sequence: String,
        max_hits: usize,
        max_hits_explicit: bool,
        task: Option<String>,
        request_options_json: Option<Value>,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
    },
    ReferenceBlastAsyncStatus {
        helper_mode: bool,
        job_id: String,
        include_report: bool,
    },
    ReferenceBlastAsyncCancel {
        helper_mode: bool,
        job_id: String,
    },
    ReferenceBlastAsyncList {
        helper_mode: bool,
    },
    ReferenceBlastTrack {
        helper_mode: bool,
        genome_id: String,
        query_sequence: String,
        target_seq_id: String,
        max_hits: usize,
        max_hits_explicit: bool,
        task: Option<String>,
        request_options_json: Option<Value>,
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
    MacrosInstanceList,
    MacrosInstanceShow {
        macro_instance_id: String,
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
        input_ports: Vec<WorkflowMacroTemplatePort>,
        output_ports: Vec<WorkflowMacroTemplatePort>,
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
        validate_only: bool,
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
    PrimersDesign {
        request_json: String,
        backend: Option<PrimerDesignBackend>,
        primer3_executable: Option<String>,
    },
    PrimersSeedFromFeature {
        seq_id: String,
        feature_id: usize,
    },
    PrimersSeedFromSplicing {
        seq_id: String,
        feature_id: usize,
    },
    PrimersDesignQpcr {
        request_json: String,
        backend: Option<PrimerDesignBackend>,
        primer3_executable: Option<String>,
    },
    PrimersPreflight {
        backend: Option<PrimerDesignBackend>,
        primer3_executable: Option<String>,
    },
    PrimersListReports,
    PrimersShowReport {
        report_id: String,
    },
    PrimersExportReport {
        report_id: String,
        path: String,
    },
    PrimersListQpcrReports,
    PrimersShowQpcrReport {
        report_id: String,
    },
    PrimersExportQpcrReport {
        report_id: String,
        path: String,
    },
    DotplotCompute {
        seq_id: String,
        span_start_0based: Option<usize>,
        span_end_0based: Option<usize>,
        mode: DotplotMode,
        word_size: usize,
        step_bp: usize,
        max_mismatches: usize,
        tile_bp: Option<usize>,
        dotplot_id: Option<String>,
    },
    DotplotList {
        seq_id: Option<String>,
    },
    DotplotShow {
        dotplot_id: String,
    },
    FlexCompute {
        seq_id: String,
        span_start_0based: Option<usize>,
        span_end_0based: Option<usize>,
        model: FlexibilityModel,
        bin_bp: usize,
        smoothing_bp: Option<usize>,
        track_id: Option<String>,
    },
    FlexList {
        seq_id: Option<String>,
    },
    FlexShow {
        track_id: String,
    },
    RnaReadsInterpret {
        seq_id: String,
        seed_feature_id: usize,
        input_path: String,
        profile: RnaReadInterpretationProfile,
        input_format: RnaReadInputFormat,
        scope: SplicingScopePreset,
        origin_mode: RnaReadOriginMode,
        target_gene_ids: Vec<String>,
        roi_seed_capture_enabled: bool,
        seed_filter: RnaReadSeedFilterConfig,
        align_config: RnaReadAlignConfig,
        report_id: Option<String>,
        report_mode: RnaReadReportMode,
        checkpoint_path: Option<String>,
        checkpoint_every_reads: usize,
        resume_from_checkpoint: bool,
    },
    RnaReadsListReports {
        seq_id: Option<String>,
    },
    RnaReadsShowReport {
        report_id: String,
    },
    RnaReadsExportReport {
        report_id: String,
        path: String,
    },
    RnaReadsExportHitsFasta {
        report_id: String,
        path: String,
        selection: RnaReadHitSelection,
    },
    RnaReadsExportSampleSheet {
        path: String,
        seq_id: Option<String>,
        report_ids: Vec<String>,
        append: bool,
    },
    RnaReadsExportExonPathsTsv {
        report_id: String,
        path: String,
        selection: RnaReadHitSelection,
    },
    RnaReadsExportExonAbundanceTsv {
        report_id: String,
        path: String,
        selection: RnaReadHitSelection,
    },
    RnaReadsExportScoreDensitySvg {
        report_id: String,
        path: String,
        scale: RnaReadScoreDensityScale,
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

fn parse_workflow_template_port_spec(raw: &str) -> Result<WorkflowMacroTemplatePort, String> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return Err("Invalid template port: value is empty".to_string());
    }
    let parts = trimmed.split(':').map(str::trim).collect::<Vec<_>>();
    if parts.len() < 2 {
        return Err(format!(
            "Invalid template port '{}': expected PORT_ID:KIND[:one|many][:required|optional][:description]",
            raw
        ));
    }
    let port_id = parts[0];
    if port_id.is_empty() {
        return Err(format!("Invalid template port '{}': missing PORT_ID", raw));
    }
    let kind = parts[1].to_ascii_lowercase();
    if kind.is_empty() {
        return Err(format!("Invalid template port '{}': missing KIND", raw));
    }
    let mut cardinality = "one".to_string();
    let mut required = true;
    let mut description_parts: Vec<String> = vec![];
    for token in parts.into_iter().skip(2) {
        if token.is_empty() {
            continue;
        }
        let lower = token.to_ascii_lowercase();
        if lower == "one" || lower == "many" {
            cardinality = lower;
            continue;
        }
        if lower == "required" || lower == "req" {
            required = true;
            continue;
        }
        if lower == "optional" || lower == "opt" {
            required = false;
            continue;
        }
        if let Some(parsed) = parse_bool_binding(&lower) {
            required = parsed;
            continue;
        }
        description_parts.push(token.to_string());
    }
    let description = if description_parts.is_empty() {
        None
    } else {
        Some(description_parts.join(":"))
    };
    Ok(WorkflowMacroTemplatePort {
        port_id: port_id.to_string(),
        kind,
        required,
        cardinality,
        description,
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
                input_ports: file.input_ports,
                output_ports: file.output_ports,
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
        if !seen_ids.insert(routine.routine_id.trim().to_ascii_lowercase()) {
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
        normalize_optional_string(&mut routine.summary);
        normalize_optional_string(&mut routine.purpose);
        normalize_optional_string(&mut routine.mechanism);
        normalize_optional_string(&mut routine.details_url);
        normalize_optional_string(&mut routine.template_path);

        normalize_string_list_ordered(&mut routine.requires);
        normalize_string_list_ordered(&mut routine.contraindications);
        normalize_string_list_ordered(&mut routine.confusing_alternatives);
        normalize_string_list_ordered(&mut routine.disambiguation_questions);
        normalize_string_list_ordered(&mut routine.failure_modes);
        normalize_routine_difference_axes(&mut routine.difference_matrix);

        routine
            .confusing_alternatives
            .retain(|entry| !entry.eq_ignore_ascii_case(routine.routine_id.as_str()));

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
    for routine in &catalog.routines {
        for alt in &routine.confusing_alternatives {
            if !seen_ids.contains(&alt.to_ascii_lowercase()) {
                return Err(format!(
                    "Cloning routine '{}' in '{}' references unknown confusing_alternative '{}'",
                    routine.routine_id, path, alt
                ));
            }
        }
    }
    Ok(catalog)
}

fn normalize_optional_string(value: &mut Option<String>) {
    *value = value
        .as_deref()
        .map(str::trim)
        .filter(|entry| !entry.is_empty())
        .map(ToOwned::to_owned);
}

fn normalize_string_list_ordered(values: &mut Vec<String>) {
    let mut out = Vec::with_capacity(values.len());
    let mut seen = HashSet::new();
    for entry in values
        .iter()
        .map(|value| value.trim())
        .filter(|v| !v.is_empty())
    {
        let key = entry.to_ascii_lowercase();
        if seen.insert(key) {
            out.push(entry.to_string());
        }
    }
    *values = out;
}

fn normalize_routine_difference_axes(values: &mut Vec<CloningRoutineDifferenceAxis>) {
    let mut out = Vec::with_capacity(values.len());
    let mut seen = HashSet::new();
    for row in values.iter() {
        let axis = row.axis.trim();
        let value = row.value.trim();
        if axis.is_empty() || value.is_empty() {
            continue;
        }
        let key = axis.to_ascii_lowercase();
        if seen.insert(key) {
            out.push(CloningRoutineDifferenceAxis {
                axis: axis.to_string(),
                value: value.to_string(),
            });
        }
    }
    *values = out;
}

fn resolve_catalog_routine<'a>(
    catalog: &'a CloningRoutineCatalog,
    routine_id: &str,
) -> Option<&'a CloningRoutineDefinition> {
    let needle = routine_id.trim();
    if needle.is_empty() {
        return None;
    }
    catalog
        .routines
        .iter()
        .find(|routine| routine.routine_id.eq_ignore_ascii_case(needle))
}

fn build_default_routine_requirements(routine: &CloningRoutineDefinition) -> Vec<String> {
    routine
        .input_ports
        .iter()
        .filter(|port| port.required)
        .map(|port| {
            let desc = port
                .description
                .as_deref()
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .unwrap_or("required input");
            format!("{} ({}, {})", port.port_id.trim(), port.kind.trim(), desc)
        })
        .collect::<Vec<_>>()
}

fn routine_summary_row(routine: &CloningRoutineDefinition) -> Value {
    json!({
        "routine_id": routine.routine_id,
        "title": routine.title,
        "family": routine.family,
        "status": routine.status,
        "summary": routine.summary,
        "template_name": routine.template_name,
    })
}

fn build_routine_axis_map(routine: &CloningRoutineDefinition) -> HashMap<String, (String, String)> {
    let mut map = HashMap::new();
    for row in &routine.difference_matrix {
        let key = row.axis.to_ascii_lowercase();
        map.insert(key, (row.axis.clone(), row.value.clone()));
    }
    map
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
            routine
                .purpose
                .as_deref()
                .unwrap_or_default()
                .to_ascii_lowercase(),
            routine
                .mechanism
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
        haystack.extend(
            routine
                .requires
                .iter()
                .map(|value| value.to_ascii_lowercase()),
        );
        haystack.extend(
            routine
                .contraindications
                .iter()
                .map(|value| value.to_ascii_lowercase()),
        );
        haystack.extend(
            routine
                .disambiguation_questions
                .iter()
                .map(|value| value.to_ascii_lowercase()),
        );
        haystack.extend(
            routine
                .failure_modes
                .iter()
                .map(|value| value.to_ascii_lowercase()),
        );
        haystack.extend(routine.difference_matrix.iter().flat_map(|row| {
            [
                row.axis.to_ascii_lowercase(),
                row.value.to_ascii_lowercase(),
            ]
        }));
        if !haystack.iter().any(|entry| entry.contains(&query_filter)) {
            return false;
        }
    }
    true
}

fn normalize_planning_class_key(raw: &str) -> String {
    let mut out = String::with_capacity(raw.len());
    for ch in raw.trim().chars() {
        if ch.is_ascii_alphanumeric() {
            out.push(ch.to_ascii_lowercase());
        } else if matches!(ch, '_' | '-' | '.' | ':') {
            out.push('_');
        }
    }
    let trimmed = out.trim_matches('_');
    if trimmed.is_empty() {
        "item".to_string()
    } else {
        trimmed.to_string()
    }
}

fn normalize_planning_class_set(values: &[String]) -> Vec<String> {
    let mut out = values
        .iter()
        .map(|value| normalize_planning_class_key(value))
        .filter(|value| !value.is_empty())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect::<Vec<_>>();
    out.sort_by_key(|value| value.to_ascii_lowercase());
    out
}

fn routine_requirement_haystack(routine: &CloningRoutineDefinition) -> String {
    let mut fields = vec![
        routine.routine_id.as_str(),
        routine.title.as_str(),
        routine.family.as_str(),
        routine.summary.as_deref().unwrap_or_default(),
        routine.purpose.as_deref().unwrap_or_default(),
        routine.mechanism.as_deref().unwrap_or_default(),
    ];
    for req in &routine.requires {
        fields.push(req.as_str());
    }
    for tag in &routine.vocabulary_tags {
        fields.push(tag.as_str());
    }
    fields.join(" ").to_ascii_lowercase()
}

fn infer_routine_material_classes(routine: &CloningRoutineDefinition) -> Vec<String> {
    if !routine.required_material_classes.is_empty() {
        return normalize_planning_class_set(&routine.required_material_classes);
    }
    let mut classes: BTreeSet<String> = BTreeSet::new();
    let family = routine.family.trim().to_ascii_lowercase();
    match family.as_str() {
        "restriction" => {
            classes.insert("restriction_enzymes".to_string());
            classes.insert("ligase".to_string());
        }
        "gibson" => {
            classes.insert("gibson_master_mix".to_string());
        }
        "infusion" => {
            classes.insert("infusion_kit".to_string());
        }
        "nebuilder_hifi" | "nebuilder" => {
            classes.insert("nebuilder_hifi_mix".to_string());
        }
        "gateway" => {
            classes.insert("gateway_clonase".to_string());
        }
        _ => {}
    }
    let text = routine_requirement_haystack(routine);
    if text.contains("restriction") || text.contains("enzyme") {
        classes.insert("restriction_enzymes".to_string());
    }
    if text.contains("gibson") {
        classes.insert("gibson_master_mix".to_string());
    }
    if text.contains("ligase") || text.contains("ligation") {
        classes.insert("ligase".to_string());
    }
    if text.contains("primer") {
        classes.insert("primers".to_string());
    }
    if text.contains("polymerase") || text.contains("pcr") {
        classes.insert("pcr_master_mix".to_string());
    }
    classes
        .into_iter()
        .map(|value| normalize_planning_class_key(&value))
        .collect::<Vec<_>>()
}

fn infer_routine_machine_classes(routine: &CloningRoutineDefinition) -> Vec<String> {
    if !routine.required_machine_classes.is_empty() {
        return normalize_planning_class_set(&routine.required_machine_classes);
    }
    let mut classes: BTreeSet<String> = BTreeSet::new();
    let family = routine.family.trim().to_ascii_lowercase();
    if matches!(
        family.as_str(),
        "restriction" | "gibson" | "infusion" | "nebuilder_hifi" | "nebuilder" | "gateway"
    ) {
        classes.insert("thermocycler".to_string());
    }
    let text = routine_requirement_haystack(routine);
    if text.contains("gel") {
        classes.insert("gel_imager".to_string());
    }
    if text.contains("sequenc") {
        classes.insert("sequencer".to_string());
    }
    classes
        .into_iter()
        .map(|value| normalize_planning_class_key(&value))
        .collect::<Vec<_>>()
}

fn infer_routine_capabilities(
    routine: &CloningRoutineDefinition,
    machine_classes: &[String],
) -> Vec<String> {
    let mut caps = normalize_planning_class_set(&routine.required_capabilities);
    for machine in machine_classes {
        if !caps.iter().any(|existing| existing == machine) {
            caps.push(machine.clone());
        }
    }
    caps.sort_by_key(|value| value.to_ascii_lowercase());
    caps.dedup();
    caps
}

fn default_routine_base_time_hours(routine: &CloningRoutineDefinition) -> f64 {
    if let Some(value) = routine
        .base_time_hours
        .filter(|value| value.is_finite() && *value > 0.0)
    {
        return value;
    }
    match routine.family.trim().to_ascii_lowercase().as_str() {
        "restriction" => 6.0,
        "gibson" => 5.0,
        "infusion" => 4.0,
        "nebuilder_hifi" | "nebuilder" => 4.0,
        "gateway" => 8.0,
        "sequence" => 1.0,
        _ => 4.0,
    }
}

fn default_routine_base_cost(routine: &CloningRoutineDefinition) -> f64 {
    routine
        .base_cost
        .filter(|value| value.is_finite() && *value >= 0.0)
        .unwrap_or(0.0)
}

fn planning_business_days_to_elapsed_hours(days: f64) -> f64 {
    if !days.is_finite() || days <= 0.0 {
        return 0.0;
    }
    // v1 approximation without a calendar anchor:
    // treat business days as Monday-Friday, then convert to elapsed hours.
    days * 24.0 * (7.0 / 5.0)
}

fn estimate_routine_planning(
    engine: &GentleEngine,
    routine: &CloningRoutineDefinition,
) -> PlanningEstimate {
    let profile = engine.planning_effective_profile();
    let objective = engine.planning_objective();
    let material_classes = infer_routine_material_classes(routine);
    let machine_classes = infer_routine_machine_classes(routine);
    let required_capabilities = infer_routine_capabilities(routine, &machine_classes);

    let base_time_hours = default_routine_base_time_hours(routine);
    let mut estimated_time_hours = base_time_hours;
    let mut estimated_cost = default_routine_base_cost(routine);
    let procurement_default_days = if profile.procurement_business_days_default.is_finite()
        && profile.procurement_business_days_default > 0.0
    {
        profile.procurement_business_days_default
    } else {
        10.0
    };

    let mut missing_materials: Vec<String> = vec![];
    let mut unknown_materials: Vec<String> = vec![];
    let mut procurement_delay_business_days = 0.0;
    for class in &material_classes {
        if let Some(item) = profile.inventory.get(class) {
            if let Some(unit_cost) = item
                .unit_cost
                .filter(|value| value.is_finite() && *value >= 0.0)
            {
                estimated_cost += unit_cost;
            }
            if !item.available {
                missing_materials.push(class.clone());
                procurement_delay_business_days += item
                    .procurement_business_days
                    .filter(|value| value.is_finite() && *value > 0.0)
                    .unwrap_or(procurement_default_days);
            }
        } else {
            unknown_materials.push(class.clone());
        }
    }

    let mut missing_machines: Vec<String> = vec![];
    for class in &machine_classes {
        if let Some(machine) = profile.machine_availability.get(class) {
            if !machine.available {
                missing_machines.push(class.clone());
            }
            if machine.queue_business_days.is_finite() && machine.queue_business_days > 0.0 {
                estimated_time_hours +=
                    planning_business_days_to_elapsed_hours(machine.queue_business_days);
            }
            if let Some(run_cost_per_hour) = machine
                .run_cost_per_hour
                .filter(|value| value.is_finite() && *value >= 0.0)
            {
                estimated_cost += run_cost_per_hour * base_time_hours;
            }
        }
    }
    estimated_time_hours += planning_business_days_to_elapsed_hours(procurement_delay_business_days);

    let capability_set = profile
        .capabilities
        .iter()
        .map(|value| normalize_planning_class_key(value))
        .collect::<HashSet<_>>();
    let mut missing_capabilities: Vec<String> = vec![];
    for capability in &required_capabilities {
        if !capability_set.contains(capability) {
            missing_capabilities.push(capability.clone());
        }
    }

    let mut local_fit_score = 1.0
        - (missing_materials.len() as f64 * 0.20)
        - (missing_machines.len() as f64 * 0.25)
        - (missing_capabilities.len() as f64 * 0.15);
    if !local_fit_score.is_finite() {
        local_fit_score = 0.0;
    }
    local_fit_score = local_fit_score.clamp(0.0, 1.0);

    let mut guardrail_failures = vec![];
    if objective.enforce_guardrails {
        if let Some(max_cost) = objective.max_cost {
            if estimated_cost > max_cost {
                guardrail_failures.push(format!(
                    "max_cost_exceeded ({estimated_cost:.2} > {max_cost:.2})"
                ));
            }
        }
        if let Some(max_time_hours) = objective.max_time_hours {
            if estimated_time_hours > max_time_hours {
                guardrail_failures.push(format!(
                    "max_time_hours_exceeded ({estimated_time_hours:.1} > {max_time_hours:.1})"
                ));
            }
        }
        for capability in &objective.required_capabilities {
            if !capability_set.contains(capability) {
                guardrail_failures.push(format!("missing_required_capability ({capability})"));
            }
        }
    }
    let passes_guardrails = guardrail_failures.is_empty();
    let composite_meta_score = (objective.weight_local_fit * local_fit_score)
        - (objective.weight_time * estimated_time_hours)
        - (objective.weight_cost * estimated_cost);

    PlanningEstimate {
        schema: PLANNING_ESTIMATE_SCHEMA.to_string(),
        estimated_time_hours,
        estimated_cost,
        local_fit_score,
        composite_meta_score,
        passes_guardrails,
        guardrail_failures,
        explanation: json!({
            "objective_schema": PLANNING_OBJECTIVE_SCHEMA,
            "profile_schema": PLANNING_PROFILE_SCHEMA,
            "base_time_hours": base_time_hours,
            "base_cost": default_routine_base_cost(routine),
            "procurement_business_days_default": procurement_default_days,
            "procurement_delay_business_days": procurement_delay_business_days,
            "required_material_classes": material_classes,
            "required_machine_classes": machine_classes,
            "required_capabilities": required_capabilities,
            "missing_material_classes": missing_materials,
            "unknown_material_classes": unknown_materials,
            "missing_machine_classes": missing_machines,
            "missing_capabilities": missing_capabilities,
            "weights": {
                "time": objective.weight_time,
                "cost": objective.weight_cost,
                "local_fit": objective.weight_local_fit,
            },
            "guardrails_enforced": objective.enforce_guardrails,
        }),
    }
}

fn resolve_workflow_template_bindings_for_preflight(
    template: &WorkflowMacroTemplate,
    bindings: &HashMap<String, String>,
) -> Result<HashMap<String, String>, String> {
    let declared = template
        .parameters
        .iter()
        .map(|p| p.name.clone())
        .collect::<HashSet<_>>();
    for key in bindings.keys() {
        if !declared.contains(key) {
            return Err(format!(
                "Workflow macro template '{}' does not define parameter '{}'",
                template.name, key
            ));
        }
    }
    let mut resolved: HashMap<String, String> = HashMap::new();
    for param in &template.parameters {
        if let Some(value) = bindings.get(&param.name) {
            resolved.insert(param.name.clone(), value.clone());
        } else if let Some(default_value) = &param.default_value {
            resolved.insert(param.name.clone(), default_value.clone());
        } else if param.required {
            return Err(format!(
                "Workflow macro template '{}' is missing required parameter '{}'",
                template.name, param.name
            ));
        }
    }
    Ok(resolved)
}

fn parse_bool_binding(raw: &str) -> Option<bool> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "true" | "1" | "yes" | "y" => Some(true),
        "false" | "0" | "no" | "n" => Some(false),
        _ => None,
    }
}

fn validate_sequence_anchor_binding(raw: &str) -> Result<SequenceAnchor, String> {
    if let Ok(anchor) = serde_json::from_str::<SequenceAnchor>(raw) {
        return Ok(anchor);
    }
    if let Ok(zero_based) = raw.trim().parse::<usize>() {
        return Ok(SequenceAnchor::Position { zero_based });
    }
    Err(format!(
        "Invalid sequence_anchor value '{}': expected SequenceAnchor JSON or zero-based integer",
        raw
    ))
}

fn split_port_values_for_validation(raw: &str, cardinality: &str) -> Result<Vec<String>, String> {
    let normalized = cardinality.trim().to_ascii_lowercase();
    match normalized.as_str() {
        "one" | "" => {
            let value = raw.trim();
            if value.is_empty() {
                return Ok(vec![]);
            }
            if value.contains(',') {
                return Err(format!(
                    "Cardinality 'one' requires exactly one value, got comma-separated '{}'",
                    raw
                ));
            }
            Ok(vec![value.to_string()])
        }
        "many" => Ok(raw
            .split(',')
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_string())
            .collect::<Vec<_>>()),
        other => Err(format!(
            "Unsupported port cardinality '{}' (expected one|many)",
            other
        )),
    }
}

fn validate_port_values(
    engine: &GentleEngine,
    kind: &str,
    values: &[String],
    direction: RoutinePortDirection,
) -> Result<(), String> {
    let normalized = kind.trim().to_ascii_lowercase();
    match normalized.as_str() {
        "sequence" => {
            if direction == RoutinePortDirection::Input {
                for value in values {
                    if !engine.state().sequences.contains_key(value) {
                        return Err(format!("Sequence '{}' was not found", value));
                    }
                }
            }
            Ok(())
        }
        "container" => {
            if direction == RoutinePortDirection::Input {
                for value in values {
                    if !engine
                        .state()
                        .container_state
                        .containers
                        .contains_key(value)
                    {
                        return Err(format!("Container '{}' was not found", value));
                    }
                }
            }
            Ok(())
        }
        "candidate_set" => {
            if direction == RoutinePortDirection::Input {
                let names = engine
                    .list_candidate_sets()
                    .into_iter()
                    .map(|set| set.name)
                    .collect::<HashSet<_>>();
                for value in values {
                    if !names.contains(value) {
                        return Err(format!("Candidate set '{}' was not found", value));
                    }
                }
            }
            Ok(())
        }
        "guide_set" => {
            if direction == RoutinePortDirection::Input {
                let names = engine
                    .list_guide_sets()
                    .into_iter()
                    .map(|set| set.guide_set_id)
                    .collect::<HashSet<_>>();
                for value in values {
                    if !names.contains(value) {
                        return Err(format!("Guide set '{}' was not found", value));
                    }
                }
            }
            Ok(())
        }
        "guide_oligo_set" => {
            if direction == RoutinePortDirection::Input {
                let names = engine
                    .list_guide_oligo_sets(None)
                    .into_iter()
                    .map(|set| set.oligo_set_id)
                    .collect::<HashSet<_>>();
                for value in values {
                    if !names.contains(value) {
                        return Err(format!("Guide oligo set '{}' was not found", value));
                    }
                }
            }
            Ok(())
        }
        "string" | "path" => Ok(()),
        "number" => {
            for value in values {
                value
                    .parse::<f64>()
                    .map_err(|e| format!("Number value '{}' is invalid: {}", value, e))?;
            }
            Ok(())
        }
        "bool" | "boolean" => {
            for value in values {
                if parse_bool_binding(value).is_none() {
                    return Err(format!(
                        "Boolean value '{}' is invalid (expected true|false|1|0|yes|no)",
                        value
                    ));
                }
            }
            Ok(())
        }
        "sequence_anchor" => {
            for value in values {
                let _ = validate_sequence_anchor_binding(value)?;
            }
            Ok(())
        }
        other => Err(format!("Unsupported routine port kind '{}'", other)),
    }
}

fn normalize_port_kind(kind: &str) -> String {
    kind.trim().to_ascii_lowercase()
}

fn is_entity_port_kind(kind: &str) -> bool {
    matches!(
        normalize_port_kind(kind).as_str(),
        "sequence" | "container" | "candidate_set" | "guide_set" | "guide_oligo_set"
    )
}

fn existing_values_for_port_kind(engine: &GentleEngine, kind: &str) -> HashSet<String> {
    match normalize_port_kind(kind).as_str() {
        "sequence" => engine
            .state()
            .sequences
            .keys()
            .cloned()
            .collect::<HashSet<_>>(),
        "container" => engine
            .state()
            .container_state
            .containers
            .keys()
            .cloned()
            .collect::<HashSet<_>>(),
        "candidate_set" => engine
            .list_candidate_sets()
            .into_iter()
            .map(|set| set.name)
            .collect::<HashSet<_>>(),
        "guide_set" => engine
            .list_guide_sets()
            .into_iter()
            .map(|set| set.guide_set_id)
            .collect::<HashSet<_>>(),
        "guide_oligo_set" => engine
            .list_guide_oligo_sets(None)
            .into_iter()
            .map(|set| set.oligo_set_id)
            .collect::<HashSet<_>>(),
        _ => HashSet::new(),
    }
}

fn sequence_anchor_matches_filter(
    anchor_kind: Option<&String>,
    anchor_label: Option<&String>,
    feature: &gb_io::seq::Feature,
) -> bool {
    if feature.kind.to_string().eq_ignore_ascii_case("SOURCE") {
        return false;
    }
    if let Some(expected_kind) = anchor_kind {
        if !feature.kind.to_string().eq_ignore_ascii_case(expected_kind) {
            return false;
        }
    }
    if let Some(expected_label) = anchor_label {
        let expected = expected_label.to_ascii_uppercase();
        let mut found = false;
        for key in [
            "label",
            "gene",
            "locus_tag",
            "product",
            "standard_name",
            "note",
        ] {
            for value in feature.qualifier_values(key.into()) {
                let upper = value.to_ascii_uppercase();
                if upper == expected || upper.contains(&expected) {
                    found = true;
                    break;
                }
            }
            if found {
                break;
            }
        }
        if !found {
            return false;
        }
    }
    true
}

fn validate_sequence_anchor_against_sequence(
    dna: &DNAsequence,
    anchor: &SequenceAnchor,
    anchor_port_id: &str,
) -> Result<(), String> {
    match anchor {
        SequenceAnchor::Position { zero_based } => {
            if *zero_based > dna.len() {
                return Err(format!(
                    "position {} is out of bounds for sequence length {}",
                    zero_based,
                    dna.len()
                ));
            }
            Ok(())
        }
        SequenceAnchor::FeatureBoundary {
            feature_kind,
            feature_label,
            occurrence,
            ..
        } => {
            let mut match_count = 0usize;
            for feature in dna.features() {
                if !sequence_anchor_matches_filter(
                    feature_kind.as_ref(),
                    feature_label.as_ref(),
                    feature,
                ) {
                    continue;
                }
                let Ok((from, to)) = feature.location.find_bounds() else {
                    continue;
                };
                if from < 0 || to < 0 {
                    continue;
                }
                let from = from as usize;
                let to = to as usize;
                if from > dna.len() || to > dna.len() {
                    continue;
                }
                match_count += 1;
            }
            if match_count == 0 {
                return Err(format!(
                    "anchor '{}' did not match any feature on this sequence",
                    anchor_port_id
                ));
            }
            let idx = occurrence.unwrap_or(0);
            if idx >= match_count {
                return Err(format!(
                    "anchor '{}' requested occurrence {} but only {} match(es) exist",
                    anchor_port_id, idx, match_count
                ));
            }
            Ok(())
        }
    }
}

fn apply_cross_port_semantics(engine: &GentleEngine, report: &mut MacroTemplatePreflightReport) {
    // Cross-port checks stay family-agnostic; they run for any typed contract.
    let checked_rows = report
        .checked_ports
        .iter()
        .filter(|row| matches!(row.status, RoutinePortValidationStatus::Ok))
        .collect::<Vec<_>>();
    if checked_rows.is_empty() {
        return;
    }

    let input_rows = checked_rows
        .iter()
        .copied()
        .filter(|row| row.direction == RoutinePortDirection::Input)
        .collect::<Vec<_>>();
    let output_rows = checked_rows
        .iter()
        .copied()
        .filter(|row| row.direction == RoutinePortDirection::Output)
        .collect::<Vec<_>>();

    // 1) Output entity alias collisions (same concrete identifier used by multiple output ports).
    let mut output_aliases: HashMap<(String, String), Vec<String>> = HashMap::new();
    for row in &output_rows {
        if !is_entity_port_kind(&row.kind) {
            continue;
        }
        let kind = normalize_port_kind(&row.kind);
        for value in &row.values {
            output_aliases
                .entry((kind.clone(), value.clone()))
                .or_default()
                .push(row.port_id.clone());
        }
    }
    for ((kind, value), ports) in output_aliases {
        if ports.len() <= 1 {
            continue;
        }
        report.errors.push(format!(
            "Output alias conflict: {} value '{}' is bound by multiple output ports ({})",
            kind,
            value,
            ports.join(", ")
        ));
    }

    // 2) Input/output alias relationship diagnostics (same entity id reused across directions).
    let mut input_values_by_kind: HashMap<String, HashSet<String>> = HashMap::new();
    for row in &input_rows {
        if !is_entity_port_kind(&row.kind) {
            continue;
        }
        let kind = normalize_port_kind(&row.kind);
        for value in &row.values {
            input_values_by_kind
                .entry(kind.clone())
                .or_default()
                .insert(value.clone());
        }
    }
    for row in &output_rows {
        if !is_entity_port_kind(&row.kind) {
            continue;
        }
        let kind = normalize_port_kind(&row.kind);
        let Some(input_values) = input_values_by_kind.get(&kind) else {
            continue;
        };
        for value in &row.values {
            if input_values.contains(value) {
                report.warnings.push(format!(
                    "Output alias: output port '{}' reuses existing {} id '{}'",
                    row.port_id, kind, value
                ));
            }
        }
    }

    // 3) Output collision diagnostics against current state.
    let mut existing_cache: HashMap<String, HashSet<String>> = HashMap::new();
    for row in &output_rows {
        if !is_entity_port_kind(&row.kind) {
            continue;
        }
        let kind = normalize_port_kind(&row.kind);
        let existing = existing_cache
            .entry(kind.clone())
            .or_insert_with(|| existing_values_for_port_kind(engine, &kind));
        for value in &row.values {
            if existing.contains(value) {
                report.warnings.push(format!(
                    "Output {} id '{}' already exists in current state",
                    kind, value
                ));
            }
        }
    }

    // 4) Sequence/container cross-port compatibility checks.
    let input_sequences: HashSet<String> = input_rows
        .iter()
        .filter(|row| normalize_port_kind(&row.kind) == "sequence")
        .flat_map(|row| row.values.iter().cloned())
        .collect();
    let input_containers: HashSet<String> = input_rows
        .iter()
        .filter(|row| normalize_port_kind(&row.kind) == "container")
        .flat_map(|row| row.values.iter().cloned())
        .collect();
    if !input_sequences.is_empty() && !input_containers.is_empty() {
        let mut members: HashSet<String> = HashSet::new();
        for container_id in &input_containers {
            if let Some(container) = engine.state().container_state.containers.get(container_id) {
                members.extend(container.members.iter().cloned());
            }
        }
        for seq_id in &input_sequences {
            if !members.contains(seq_id) {
                report.warnings.push(format!(
                    "Cross-port mismatch: input sequence '{}' is not present in any bound input container",
                    seq_id
                ));
            }
        }
    }

    // 5) Sequence-anchor semantic validation when one concrete sequence input is bound.
    let anchor_rows = input_rows
        .iter()
        .filter(|row| normalize_port_kind(&row.kind) == "sequence_anchor")
        .copied()
        .collect::<Vec<_>>();
    if !anchor_rows.is_empty() {
        let bound_sequences = input_sequences.into_iter().collect::<Vec<_>>();
        if bound_sequences.len() == 1 {
            let seq_id = &bound_sequences[0];
            if let Some(dna) = engine.state().sequences.get(seq_id) {
                for row in anchor_rows {
                    for raw in &row.values {
                        match validate_sequence_anchor_binding(raw) {
                            Ok(anchor) => {
                                if let Err(err) = validate_sequence_anchor_against_sequence(
                                    dna,
                                    &anchor,
                                    &row.port_id,
                                ) {
                                    report.errors.push(format!(
                                        "Port '{}' anchor check failed for sequence '{}': {}",
                                        row.port_id, seq_id, err
                                    ));
                                }
                            }
                            Err(err) => {
                                report.errors.push(format!(
                                    "Port '{}' anchor parsing failed during semantic check: {}",
                                    row.port_id, err
                                ));
                            }
                        }
                    }
                }
            }
        } else if bound_sequences.is_empty() {
            report.warnings.push(
                "Anchor semantic checks skipped: no input sequence binding was available"
                    .to_string(),
            );
        } else {
            report.warnings.push(format!(
                "Anchor semantic checks skipped: {} input sequences were bound (expected one)",
                bound_sequences.len()
            ));
        }
    }
}

fn summarize_overlap_segment(raw: &str, max_len: usize) -> String {
    if raw.chars().count() <= max_len {
        return raw.to_string();
    }
    let prefix = raw.chars().take(max_len).collect::<String>();
    format!("{prefix}...")
}

fn parse_optional_usize_binding_with_context(
    bindings: &HashMap<String, String>,
    key: &str,
    context_label: &str,
    errors: &mut Vec<String>,
) -> Option<usize> {
    let raw = bindings.get(key)?;
    let value = raw.trim();
    if value.is_empty() {
        return None;
    }
    match value.parse::<usize>() {
        Ok(parsed) => Some(parsed),
        Err(err) => {
            errors.push(format!(
                "{context_label} expects '{}' to be a non-negative integer, got '{}': {}",
                key, raw, err
            ));
            None
        }
    }
}

fn parse_optional_usize_binding(
    bindings: &HashMap<String, String>,
    key: &str,
    errors: &mut Vec<String>,
) -> Option<usize> {
    parse_optional_usize_binding_with_context(bindings, key, "Restriction preflight", errors)
}

fn collect_restriction_enzyme_bindings(bindings: &HashMap<String, String>) -> Vec<String> {
    let mut ordered_values: Vec<String> = vec![];
    let preferred_keys = [
        "enzyme_a",
        "enzyme_b",
        "enzymes",
        "digest_enzymes",
        "restriction_enzymes",
        "enzyme",
    ];
    for key in preferred_keys {
        if let Some(raw) = bindings.get(key) {
            ordered_values.push(raw.clone());
        }
    }
    let mut extra_keys = bindings.keys().cloned().collect::<Vec<_>>();
    extra_keys.sort_by_key(|k| k.to_ascii_lowercase());
    for key in extra_keys {
        if preferred_keys.contains(&key.as_str()) {
            continue;
        }
        if key.to_ascii_lowercase().contains("enzyme") {
            if let Some(raw) = bindings.get(&key) {
                ordered_values.push(raw.clone());
            }
        }
    }

    let mut names: Vec<String> = vec![];
    for raw in ordered_values {
        for token in raw.split(',') {
            let name = token.trim();
            if name.is_empty() {
                continue;
            }
            names.push(name.to_string());
        }
    }
    names
}

fn collect_sequence_input_values(report: &MacroTemplatePreflightReport) -> Vec<String> {
    report
        .checked_ports
        .iter()
        .filter(|row| {
            row.direction == RoutinePortDirection::Input
                && matches!(row.status, RoutinePortValidationStatus::Ok)
                && normalize_port_kind(&row.kind) == "sequence"
        })
        .flat_map(|row| row.values.iter().cloned())
        .collect::<Vec<_>>()
}

fn split_csv_tokens_with_empty_error(raw: &str) -> Result<Vec<String>, String> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return Ok(vec![]);
    }
    let mut tokens = Vec::new();
    for token in trimmed.split(',') {
        let value = token.trim();
        if value.is_empty() {
            return Err(format!(
                "Token list '{}' contains an empty token (double comma or trailing comma)",
                raw
            ));
        }
        tokens.push(value.to_string());
    }
    Ok(tokens)
}

fn normalize_compact_token(raw: &str) -> String {
    raw.chars()
        .filter(|c| c.is_ascii_alphanumeric())
        .map(|c| c.to_ascii_lowercase())
        .collect::<String>()
}

fn is_type_iis_capable_enzyme_name(name: &str) -> bool {
    matches!(
        normalize_compact_token(name).as_str(),
        "eco31"
            | "eco31i"
            | "bsai"
            | "bsmbi"
            | "esp3i"
            | "bbsi"
            | "aari"
            | "sapi"
            | "btgzi"
            | "bsmai"
            | "bfuai"
    )
}

fn apply_overlap_family_preflight_semantics(
    engine: &GentleEngine,
    report: &mut MacroTemplatePreflightReport,
    family_label: &str,
    recommended_min_overlap: usize,
    recommended_max_overlap: usize,
) {
    let checked_rows = report
        .checked_ports
        .iter()
        .filter(|row| {
            row.direction == RoutinePortDirection::Input
                && matches!(row.status, RoutinePortValidationStatus::Ok)
        })
        .collect::<Vec<_>>();
    if checked_rows.is_empty() {
        return;
    }

    let sequence_values = collect_sequence_input_values(report);
    if sequence_values.len() < 2 {
        report.errors.push(format!(
            "{family_label} preflight requires at least two sequence inputs, found {}",
            sequence_values.len()
        ));
        return;
    }

    let overlap_bp = checked_rows
        .iter()
        .find(|row| normalize_port_kind(&row.kind) == "number")
        .and_then(|row| row.values.first())
        .and_then(|raw| raw.trim().parse::<usize>().ok());
    let Some(overlap_bp) = overlap_bp else {
        report.errors.push(format!(
            "{family_label} preflight requires one numeric overlap input port (for example overlap_bp)"
        ));
        return;
    };
    if overlap_bp == 0 {
        report
            .errors
            .push(format!("{family_label} preflight requires overlap_bp >= 1"));
        return;
    }
    if overlap_bp < recommended_min_overlap {
        report.warnings.push(format!(
            "{family_label} overlap_bp={} is short for typical design (often >={} bp)",
            overlap_bp, recommended_min_overlap
        ));
    }
    if overlap_bp > recommended_max_overlap {
        report.warnings.push(format!(
            "{family_label} overlap_bp={} is unusually long; verify primer and synthesis practicality",
            overlap_bp
        ));
    }

    for pair in sequence_values.windows(2) {
        let left_id = &pair[0];
        let right_id = &pair[1];
        let Some(left) = engine.state().sequences.get(left_id) else {
            continue;
        };
        let Some(right) = engine.state().sequences.get(right_id) else {
            continue;
        };
        let left_seq = left.get_forward_string().to_ascii_uppercase();
        let right_seq = right.get_forward_string().to_ascii_uppercase();
        if left_seq.len() < overlap_bp {
            report.errors.push(format!(
                "{family_label} overlap check failed: sequence '{}' length {} is shorter than overlap_bp={}",
                left_id,
                left_seq.len(),
                overlap_bp
            ));
            continue;
        }
        if right_seq.len() < overlap_bp {
            report.errors.push(format!(
                "{family_label} overlap check failed: sequence '{}' length {} is shorter than overlap_bp={}",
                right_id,
                right_seq.len(),
                overlap_bp
            ));
            continue;
        }
        let left_suffix = &left_seq[left_seq.len() - overlap_bp..];
        let right_prefix = &right_seq[..overlap_bp];
        if left_suffix != right_prefix {
            report.errors.push(format!(
                "{family_label} overlap mismatch between '{}' and '{}': left suffix '{}' != right prefix '{}' (overlap_bp={})",
                left_id,
                right_id,
                summarize_overlap_segment(left_suffix, 24),
                summarize_overlap_segment(right_prefix, 24),
                overlap_bp
            ));
        }
    }
}

fn apply_gibson_family_preflight_semantics(
    engine: &GentleEngine,
    report: &mut MacroTemplatePreflightReport,
) {
    apply_overlap_family_preflight_semantics(engine, report, "Gibson", 15, 120);
}

fn apply_restriction_family_preflight_semantics(
    engine: &GentleEngine,
    report: &mut MacroTemplatePreflightReport,
    bindings: &HashMap<String, String>,
) {
    // Restriction-family checks focus on enzyme resolution/site presence and
    // digest parameter sanity before any mutating macro execution.
    let checked_rows = report
        .checked_ports
        .iter()
        .filter(|row| {
            row.direction == RoutinePortDirection::Input
                && matches!(row.status, RoutinePortValidationStatus::Ok)
        })
        .collect::<Vec<_>>();
    if checked_rows.is_empty() {
        return;
    }

    let sequence_values = checked_rows
        .iter()
        .filter(|row| normalize_port_kind(&row.kind) == "sequence")
        .flat_map(|row| row.values.iter().cloned())
        .collect::<Vec<_>>();
    if sequence_values.is_empty() {
        report
            .errors
            .push("Restriction preflight requires at least one sequence input".to_string());
        return;
    }

    let enzyme_names = collect_restriction_enzyme_bindings(bindings);
    if enzyme_names.is_empty() {
        report.warnings.push(
            "Restriction preflight could not infer enzyme bindings; enzyme/site checks were skipped"
                .to_string(),
        );
    } else {
        let catalog = active_restriction_enzymes();
        if catalog.is_empty() {
            report.warnings.push(
                "Restriction preflight could not load active restriction-enzyme catalog; enzyme/site checks were skipped"
                    .to_string(),
            );
        } else {
            let mut by_name_lower = HashMap::new();
            for enzyme in catalog {
                by_name_lower.insert(enzyme.name.to_ascii_lowercase(), enzyme);
            }

            let mut valid_enzymes = vec![];
            let mut unknown_names = vec![];
            for name in &enzyme_names {
                let lookup = name.to_ascii_lowercase();
                if let Some(enzyme) = by_name_lower.get(&lookup) {
                    valid_enzymes.push((name.clone(), enzyme.clone()));
                } else {
                    unknown_names.push(name.clone());
                }
            }
            if !unknown_names.is_empty() {
                report.errors.push(format!(
                    "Restriction preflight unknown enzyme name(s): {}",
                    unknown_names.join(", ")
                ));
            }

            if valid_enzymes.len() >= 2 {
                let first = valid_enzymes[0].0.to_ascii_lowercase();
                let second = valid_enzymes[1].0.to_ascii_lowercase();
                if first == second {
                    report.errors.push(format!(
                        "Restriction preflight expects distinct enzyme inputs, got duplicate enzyme '{}'",
                        valid_enzymes[0].0
                    ));
                }
            }

            if !valid_enzymes.is_empty() {
                let mut per_sequence = Vec::with_capacity(sequence_values.len());
                for seq_id in &sequence_values {
                    if let Some(dna) = engine.state().sequences.get(seq_id) {
                        per_sequence.push((seq_id.as_str(), dna));
                    }
                }
                for (raw_name, enzyme) in &valid_enzymes {
                    let total_sites = per_sequence
                        .iter()
                        .map(|(_, dna)| enzyme.get_sites(dna, None).len())
                        .sum::<usize>();
                    if total_sites == 0 {
                        report.errors.push(format!(
                            "Restriction preflight: enzyme '{}' has no recognition site across input sequence(s): {}",
                            raw_name,
                            sequence_values.join(", ")
                        ));
                    } else if total_sites == 1 {
                        report.warnings.push(format!(
                            "Restriction preflight: enzyme '{}' has only one recognition site across input sequence(s)",
                            raw_name
                        ));
                    }
                }

                if valid_enzymes.len() >= 2 {
                    let has_shared_cuttable_sequence = per_sequence.iter().any(|(_, dna)| {
                        valid_enzymes
                            .iter()
                            .all(|(_, enzyme)| !enzyme.get_sites(dna, None).is_empty())
                    });
                    if !has_shared_cuttable_sequence {
                        report.errors.push(
                            "Restriction preflight: no single input sequence has recognition sites for all selected enzymes"
                                .to_string(),
                        );
                    }
                }
            }
        }
    }

    let mut local_errors: Vec<String> = vec![];
    let left_fragment = parse_optional_usize_binding(bindings, "left_fragment", &mut local_errors);
    let right_fragment =
        parse_optional_usize_binding(bindings, "right_fragment", &mut local_errors);
    let extract_from = parse_optional_usize_binding(bindings, "extract_from", &mut local_errors);
    let extract_to = parse_optional_usize_binding(bindings, "extract_to", &mut local_errors);
    report.errors.extend(local_errors);

    if left_fragment == Some(0) {
        report.errors.push(
            "Restriction preflight expects left_fragment >= 1 (digest products are 1-indexed)"
                .to_string(),
        );
    }
    if right_fragment == Some(0) {
        report.errors.push(
            "Restriction preflight expects right_fragment >= 1 (digest products are 1-indexed)"
                .to_string(),
        );
    }
    if let (Some(left), Some(right)) = (left_fragment, right_fragment) {
        if left == right {
            report.warnings.push(format!(
                "Restriction preflight: left_fragment and right_fragment are both {}",
                left
            ));
        }
    }
    if let (Some(from), Some(to)) = (extract_from, extract_to) {
        if from == to {
            report.warnings.push(
                "Restriction preflight: extract_from equals extract_to (full-length circular-style extract)"
                    .to_string(),
            );
        } else if from > to {
            report.warnings.push(
                "Restriction preflight: extract_from > extract_to (wrap-around extract semantics)"
                    .to_string(),
            );
        }
        let max_input_len = sequence_values
            .iter()
            .filter_map(|seq_id| engine.state().sequences.get(seq_id).map(|dna| dna.len()))
            .max()
            .unwrap_or(0);
        if max_input_len > 0 && to > max_input_len {
            report.warnings.push(format!(
                "Restriction preflight: extract_to={} exceeds longest input sequence length {}",
                to, max_input_len
            ));
        }
    }
}

fn apply_golden_gate_family_preflight_semantics(
    engine: &GentleEngine,
    report: &mut MacroTemplatePreflightReport,
    bindings: &HashMap<String, String>,
) {
    let sequence_values = collect_sequence_input_values(report);
    if sequence_values.is_empty() {
        report
            .errors
            .push("Golden Gate preflight requires at least one sequence input".to_string());
        return;
    }

    let enzyme_names = collect_restriction_enzyme_bindings(bindings);
    if enzyme_names.is_empty() {
        report
            .errors
            .push("Golden Gate preflight requires at least one enzyme binding".to_string());
        return;
    }

    let catalog = active_restriction_enzymes();
    if catalog.is_empty() {
        report.errors.push(
            "Golden Gate preflight could not load active restriction-enzyme catalog".to_string(),
        );
        return;
    }
    let mut by_name_lower = HashMap::new();
    for enzyme in catalog {
        by_name_lower.insert(enzyme.name.to_ascii_lowercase(), enzyme);
    }

    let mut valid_enzymes = vec![];
    let mut unknown_names: Vec<String> = vec![];
    let mut non_type_iis: Vec<String> = vec![];
    for name in &enzyme_names {
        let lookup = name.to_ascii_lowercase();
        let Some(enzyme) = by_name_lower.get(&lookup) else {
            unknown_names.push(name.clone());
            continue;
        };
        if !is_type_iis_capable_enzyme_name(&enzyme.name) {
            non_type_iis.push(enzyme.name.clone());
            continue;
        }
        valid_enzymes.push((name.clone(), enzyme.clone()));
    }
    if !unknown_names.is_empty() {
        report.errors.push(format!(
            "Golden Gate preflight unknown enzyme name(s): {}",
            unknown_names.join(", ")
        ));
    }
    if !non_type_iis.is_empty() {
        report.errors.push(format!(
            "Golden Gate preflight requires Type IIS-capable enzymes; got non-Type-IIS: {}",
            non_type_iis.join(", ")
        ));
    }
    if valid_enzymes.is_empty() {
        return;
    }

    let mut per_sequence = Vec::with_capacity(sequence_values.len());
    for seq_id in &sequence_values {
        if let Some(dna) = engine.state().sequences.get(seq_id) {
            per_sequence.push((seq_id.as_str(), dna));
        }
    }
    for (raw_name, enzyme) in &valid_enzymes {
        let mut any_site = false;
        for (seq_id, dna) in &per_sequence {
            let sites = enzyme.get_sites(dna, None);
            if sites.is_empty() {
                report.errors.push(format!(
                    "Golden Gate preflight: enzyme '{}' has no recognition site on input sequence '{}'",
                    raw_name, seq_id
                ));
            } else {
                any_site = true;
            }
        }
        if !any_site {
            report.errors.push(format!(
                "Golden Gate preflight: enzyme '{}' has no recognition site across input sequence(s): {}",
                raw_name,
                sequence_values.join(", ")
            ));
        }
    }

    let expected_junctions = sequence_values.len().saturating_sub(1);
    let junction_binding = bindings
        .get("junction_overhangs")
        .or_else(|| bindings.get("junctions"))
        .or_else(|| bindings.get("overhangs"));
    if expected_junctions >= 2 && junction_binding.is_none() {
        report.errors.push(format!(
            "Golden Gate preflight requires explicit junction_overhangs for multipart assembly (expected {} token(s))",
            expected_junctions
        ));
    }
    if let Some(raw) = junction_binding {
        match split_csv_tokens_with_empty_error(raw) {
            Ok(tokens) => {
                if expected_junctions > 0 && tokens.len() != expected_junctions {
                    report.errors.push(format!(
                        "Golden Gate preflight expected {} junction_overhang token(s), got {}",
                        expected_junctions,
                        tokens.len()
                    ));
                }
                if tokens.len() > 2 {
                    let mut seen = HashSet::new();
                    for token in tokens.iter().skip(1).take(tokens.len().saturating_sub(2)) {
                        let normalized = token.to_ascii_uppercase();
                        if !seen.insert(normalized.clone()) {
                            report.errors.push(format!(
                                "Golden Gate preflight duplicate non-terminal junction overhang '{}'",
                                normalized
                            ));
                        }
                    }
                }
            }
            Err(err) => report.errors.push(format!("Golden Gate preflight: {err}")),
        }
    }

    let mut local_errors: Vec<String> = vec![];
    let mut fragment_keys = bindings
        .keys()
        .filter(|key| key.to_ascii_lowercase().contains("fragment"))
        .cloned()
        .collect::<Vec<_>>();
    fragment_keys.sort_by_key(|value| value.to_ascii_lowercase());
    for key in fragment_keys {
        let parsed = parse_optional_usize_binding_with_context(
            bindings,
            &key,
            "Golden Gate preflight",
            &mut local_errors,
        );
        if parsed == Some(0) {
            local_errors.push(format!(
                "Golden Gate preflight expects {} >= 1 (digest products are 1-indexed)",
                key
            ));
        }
    }
    let extract_from = parse_optional_usize_binding_with_context(
        bindings,
        "extract_from",
        "Golden Gate preflight",
        &mut local_errors,
    );
    let extract_to = parse_optional_usize_binding_with_context(
        bindings,
        "extract_to",
        "Golden Gate preflight",
        &mut local_errors,
    );
    if let (Some(from), Some(to)) = (extract_from, extract_to) {
        if from == to {
            report
                .warnings
                .push("Golden Gate preflight: extract_from equals extract_to".to_string());
        } else if from > to {
            report
                .warnings
                .push("Golden Gate preflight: extract_from > extract_to".to_string());
        }
    }
    report.errors.extend(local_errors);
}

fn parse_att_token_set(raw: &str) -> Result<HashSet<String>, String> {
    let tokens = split_csv_tokens_with_empty_error(raw)?;
    Ok(tokens
        .into_iter()
        .map(|token| normalize_compact_token(&token))
        .filter(|token| !token.is_empty())
        .collect::<HashSet<_>>())
}

fn apply_gateway_family_preflight_semantics(
    report: &mut MacroTemplatePreflightReport,
    bindings: &HashMap<String, String>,
) {
    let sequence_values = collect_sequence_input_values(report);
    if sequence_values.len() < 2 {
        report.errors.push(format!(
            "Gateway preflight requires at least two sequence inputs, found {}",
            sequence_values.len()
        ));
    }

    let raw_phase = bindings
        .get("gateway_phase")
        .map(|value| value.trim().to_ascii_lowercase())
        .unwrap_or_else(|| {
            let template_name = report.template_name.to_ascii_lowercase();
            if template_name.contains("bp_lr") {
                "bp_lr".to_string()
            } else if template_name.contains("_lr_") {
                "lr".to_string()
            } else {
                "bp".to_string()
            }
        });
    let phase = match raw_phase.as_str() {
        "bp" | "bp_only" => "bp",
        "lr" | "lr_only" => "lr",
        "bp_lr" | "bp+lr" => "bp_lr",
        other => {
            report.errors.push(format!(
                "Gateway preflight unsupported gateway_phase '{}'; expected bp|lr|bp_lr",
                other
            ));
            return;
        }
    };

    let raw_tokens = bindings
        .get("att_tokens")
        .or_else(|| bindings.get("att_site_tokens"))
        .or_else(|| bindings.get("att_sites"))
        .map(|value| value.as_str())
        .unwrap_or("");
    let att_tokens = match parse_att_token_set(raw_tokens) {
        Ok(tokens) if !tokens.is_empty() => tokens,
        Ok(_) => {
            report.errors.push(
                "Gateway preflight requires non-empty att_tokens (for example 'attB,attP')"
                    .to_string(),
            );
            return;
        }
        Err(err) => {
            report
                .errors
                .push(format!("Gateway preflight att_tokens parse error: {}", err));
            return;
        }
    };

    let has_attb = att_tokens.contains("attb");
    let has_attp = att_tokens.contains("attp");
    let has_attl = att_tokens.contains("attl");
    let has_attr = att_tokens.contains("attr");
    match phase {
        "bp" => {
            if !(has_attb && has_attp) {
                report.errors.push(
                    "Gateway BP preflight requires att_tokens to include attB and attP".to_string(),
                );
            }
            if has_attl || has_attr {
                report
                    .errors
                    .push("Gateway BP preflight rejects mixed-phase tokens attL/attR".to_string());
            }
        }
        "lr" => {
            if !(has_attl && has_attr) {
                report.errors.push(
                    "Gateway LR preflight requires att_tokens to include attL and attR".to_string(),
                );
            }
            if has_attb || has_attp {
                report
                    .errors
                    .push("Gateway LR preflight rejects mixed-phase tokens attB/attP".to_string());
            }
        }
        "bp_lr" => {
            if !(has_attb && has_attp && has_attl && has_attr) {
                report.errors.push(
                    "Gateway BP+LR preflight requires att_tokens to include attB,attP,attL,attR"
                        .to_string(),
                );
            }
        }
        _ => {}
    }
}

fn apply_topo_family_preflight_semantics(
    engine: &GentleEngine,
    report: &mut MacroTemplatePreflightReport,
    bindings: &HashMap<String, String>,
) {
    let insert_id = bindings
        .get("insert_seq_id")
        .or_else(|| bindings.get("donor_seq_id"))
        .map(|value| value.trim().to_string());
    let vector_id = bindings
        .get("vector_seq_id")
        .or_else(|| bindings.get("entry_vector_seq_id"))
        .map(|value| value.trim().to_string());
    let Some(insert_id) = insert_id.filter(|value| !value.is_empty()) else {
        report
            .errors
            .push("TOPO preflight requires insert_seq_id binding".to_string());
        return;
    };
    let Some(vector_id) = vector_id.filter(|value| !value.is_empty()) else {
        report
            .errors
            .push("TOPO preflight requires vector_seq_id binding".to_string());
        return;
    };
    let Some(insert_seq) = engine.state().sequences.get(&insert_id) else {
        report.errors.push(format!(
            "TOPO preflight insert sequence '{}' was not found",
            insert_id
        ));
        return;
    };
    let Some(vector_seq) = engine.state().sequences.get(&vector_id) else {
        report.errors.push(format!(
            "TOPO preflight vector sequence '{}' was not found",
            vector_id
        ));
        return;
    };
    let insert_text = insert_seq.get_forward_string().to_ascii_uppercase();
    let vector_text = vector_seq.get_forward_string().to_ascii_uppercase();
    let mode = bindings
        .get("topo_mode")
        .map(|value| normalize_compact_token(value))
        .unwrap_or_else(|| {
            let template = report.template_name.to_ascii_lowercase();
            if template.contains("directional") {
                "directionalcacc".to_string()
            } else if template.contains("blunt") {
                "blunt".to_string()
            } else {
                "ta".to_string()
            }
        });
    match mode.as_str() {
        "ta" => {
            if !insert_text.ends_with('A') {
                report.errors.push(format!(
                    "TOPO TA preflight expects insert '{}' to end with 'A'",
                    insert_id
                ));
            }
            if !vector_text.ends_with('T') {
                report.errors.push(format!(
                    "TOPO TA preflight expects vector '{}' to end with 'T'",
                    vector_id
                ));
            }
        }
        "blunt" => {}
        "directionalcacc" | "directional" => {
            if !insert_text.starts_with("CACC") {
                report.errors.push(format!(
                    "TOPO directional preflight expects insert '{}' to start with 'CACC'",
                    insert_id
                ));
            }
        }
        other => report.errors.push(format!(
            "TOPO preflight unsupported topo_mode '{}'; expected ta|blunt|directional_cacc",
            other
        )),
    }
}

fn apply_ta_gc_family_preflight_semantics(
    engine: &GentleEngine,
    report: &mut MacroTemplatePreflightReport,
    bindings: &HashMap<String, String>,
) {
    let insert_id = bindings
        .get("insert_seq_id")
        .map(|value| value.trim().to_string());
    let vector_id = bindings
        .get("vector_seq_id")
        .map(|value| value.trim().to_string());
    let Some(insert_id) = insert_id.filter(|value| !value.is_empty()) else {
        report
            .errors
            .push("TA/GC preflight requires insert_seq_id binding".to_string());
        return;
    };
    let Some(vector_id) = vector_id.filter(|value| !value.is_empty()) else {
        report
            .errors
            .push("TA/GC preflight requires vector_seq_id binding".to_string());
        return;
    };
    let Some(insert_seq) = engine.state().sequences.get(&insert_id) else {
        report.errors.push(format!(
            "TA/GC preflight insert sequence '{}' was not found",
            insert_id
        ));
        return;
    };
    let Some(vector_seq) = engine.state().sequences.get(&vector_id) else {
        report.errors.push(format!(
            "TA/GC preflight vector sequence '{}' was not found",
            vector_id
        ));
        return;
    };
    let mode = bindings
        .get("tail_mode")
        .map(|value| normalize_compact_token(value))
        .unwrap_or_else(|| {
            if report.template_name.to_ascii_lowercase().contains("gc_") {
                "gc".to_string()
            } else {
                "ta".to_string()
            }
        });
    let insert_text = insert_seq.get_forward_string().to_ascii_uppercase();
    let vector_text = vector_seq.get_forward_string().to_ascii_uppercase();
    match mode.as_str() {
        "ta" => {
            if !insert_text.ends_with('A') {
                report.errors.push(format!(
                    "TA/GC preflight expects TA insert '{}' to end with 'A'",
                    insert_id
                ));
            }
            if !vector_text.ends_with('T') {
                report.errors.push(format!(
                    "TA/GC preflight expects TA vector '{}' to end with 'T'",
                    vector_id
                ));
            }
        }
        "gc" => {
            if !insert_text.ends_with('G') {
                report.errors.push(format!(
                    "TA/GC preflight expects GC insert '{}' to end with 'G'",
                    insert_id
                ));
            }
            if !vector_text.ends_with('C') {
                report.errors.push(format!(
                    "TA/GC preflight expects GC vector '{}' to end with 'C'",
                    vector_id
                ));
            }
        }
        other => report.errors.push(format!(
            "TA/GC preflight unsupported tail_mode '{}'; expected ta|gc",
            other
        )),
    }
}

fn apply_family_specific_preflight_semantics(
    engine: &GentleEngine,
    report: &mut MacroTemplatePreflightReport,
    bindings: &HashMap<String, String>,
) {
    // Family hooks remain explicit to avoid hidden adapter-side semantics.
    let Some(family) = report.routine_family.as_deref() else {
        return;
    };
    if family.eq_ignore_ascii_case("gibson") {
        apply_gibson_family_preflight_semantics(engine, report);
    } else if family.eq_ignore_ascii_case("infusion")
        || family.eq_ignore_ascii_case("in_fusion")
        || family.eq_ignore_ascii_case("in-fusion")
    {
        apply_overlap_family_preflight_semantics(engine, report, "In-Fusion", 12, 80);
    } else if family.eq_ignore_ascii_case("nebuilder_hifi")
        || family.eq_ignore_ascii_case("nebuilder")
        || family.eq_ignore_ascii_case("nebuilder-hifi")
    {
        apply_overlap_family_preflight_semantics(engine, report, "NEBuilder HiFi", 15, 120);
    } else if family.eq_ignore_ascii_case("restriction") {
        apply_restriction_family_preflight_semantics(engine, report, bindings);
    } else if family.eq_ignore_ascii_case("golden_gate")
        || family.eq_ignore_ascii_case("goldengate")
    {
        apply_golden_gate_family_preflight_semantics(engine, report, bindings);
    } else if family.eq_ignore_ascii_case("gateway") {
        apply_gateway_family_preflight_semantics(report, bindings);
    } else if family.eq_ignore_ascii_case("topo") {
        apply_topo_family_preflight_semantics(engine, report, bindings);
    } else if family.eq_ignore_ascii_case("ta_gc")
        || family.eq_ignore_ascii_case("tagc")
        || family.eq_ignore_ascii_case("ta-gc")
    {
        apply_ta_gc_family_preflight_semantics(engine, report, bindings);
    }
}

#[derive(Debug, Clone)]
struct PreflightPortContract {
    direction: RoutinePortDirection,
    source: String,
    port: WorkflowMacroTemplatePort,
}

fn routine_port_to_template_port(port: &CloningRoutinePort) -> WorkflowMacroTemplatePort {
    WorkflowMacroTemplatePort {
        port_id: port.port_id.clone(),
        kind: port.kind.clone(),
        required: port.required,
        cardinality: if port.cardinality.trim().is_empty() {
            "one".to_string()
        } else {
            port.cardinality.clone()
        },
        description: port.description.clone(),
    }
}

fn preflight_workflow_macro_template_run(
    engine: &GentleEngine,
    template_name: &str,
    bindings: &HashMap<String, String>,
) -> Result<MacroTemplatePreflightReport, String> {
    let template = engine
        .get_workflow_macro_template(template_name)
        .map_err(|e| e.to_string())?;
    let mut report = MacroTemplatePreflightReport {
        schema: "gentle.macro_template_preflight.v1".to_string(),
        template_name: template.name.clone(),
        routine_id: None,
        routine_title: None,
        routine_family: None,
        routine_status: None,
        catalog_path: None,
        contract_source: "none".to_string(),
        checked_ports: vec![],
        warnings: vec![],
        errors: vec![],
    };

    let resolved_bindings =
        match resolve_workflow_template_bindings_for_preflight(&template, bindings) {
            Ok(resolved) => resolved,
            Err(err) => {
                report.errors.push(err);
                return Ok(report);
            }
        };

    let catalog_path = DEFAULT_CLONING_ROUTINE_CATALOG_PATH.to_string();
    report.catalog_path = Some(catalog_path.clone());
    let mut matched_routine: Option<CloningRoutineDefinition> = None;
    match load_cloning_routine_catalog(&catalog_path) {
        Ok(catalog) => {
            let matching = catalog
                .routines
                .into_iter()
                .filter(|routine| routine.template_name == template.name)
                .collect::<Vec<_>>();
            if matching.is_empty() {
                report.warnings.push(format!(
                    "No routine catalog entry references template '{}'",
                    template.name
                ));
            } else {
                if matching.len() > 1 {
                    let candidate_labels = matching
                        .iter()
                        .map(|routine| format!("{} ({})", routine.routine_id, routine.family))
                        .collect::<Vec<_>>();
                    report.warnings.push(format!(
                        "Multiple routine catalog entries reference template '{}': {}",
                        template.name,
                        candidate_labels.join(", ")
                    ));
                    if template.input_ports.is_empty() && template.output_ports.is_empty() {
                        report.errors.push(format!(
                            "Ambiguous routine binding for template '{}' ({} matches) without explicit template port contract",
                            template.name,
                            matching.len()
                        ));
                    }
                }
                let routine = matching[0].clone();
                report.routine_id = Some(routine.routine_id.clone());
                report.routine_title = Some(routine.title.clone());
                report.routine_family = Some(routine.family.clone());
                report.routine_status = Some(routine.status.clone());
                matched_routine = Some(routine);
            }
        }
        Err(err) => {
            report.warnings.push(format!(
                "Routine catalog '{}' unavailable: {}",
                catalog_path, err
            ));
        }
    }

    let mut contracts: Vec<PreflightPortContract> = vec![];
    if !template.input_ports.is_empty() || !template.output_ports.is_empty() {
        report.contract_source = "template_ports".to_string();
        contracts.extend(
            template
                .input_ports
                .iter()
                .cloned()
                .map(|port| PreflightPortContract {
                    direction: RoutinePortDirection::Input,
                    source: "template".to_string(),
                    port,
                }),
        );
        contracts.extend(
            template
                .output_ports
                .iter()
                .cloned()
                .map(|port| PreflightPortContract {
                    direction: RoutinePortDirection::Output,
                    source: "template".to_string(),
                    port,
                }),
        );

        if let Some(routine) = matched_routine.as_ref() {
            let routine_port_ids = routine
                .input_ports
                .iter()
                .map(|p| format!("input:{}", p.port_id))
                .chain(
                    routine
                        .output_ports
                        .iter()
                        .map(|p| format!("output:{}", p.port_id)),
                )
                .collect::<HashSet<_>>();
            let template_port_ids = template
                .input_ports
                .iter()
                .map(|p| format!("input:{}", p.port_id))
                .chain(
                    template
                        .output_ports
                        .iter()
                        .map(|p| format!("output:{}", p.port_id)),
                )
                .collect::<HashSet<_>>();
            for missing in routine_port_ids.difference(&template_port_ids) {
                report.warnings.push(format!(
                    "Template contract is missing routine-declared port '{}'",
                    missing
                ));
            }
            for extra in template_port_ids.difference(&routine_port_ids) {
                report.warnings.push(format!(
                    "Template contract contains extra port '{}' not present in routine catalog",
                    extra
                ));
            }
        }
    } else if let Some(routine) = matched_routine {
        report.contract_source = "routine_catalog".to_string();
        contracts.extend(
            routine
                .input_ports
                .into_iter()
                .map(|port| PreflightPortContract {
                    direction: RoutinePortDirection::Input,
                    source: "routine_catalog".to_string(),
                    port: routine_port_to_template_port(&port),
                }),
        );
        contracts.extend(
            routine
                .output_ports
                .into_iter()
                .map(|port| PreflightPortContract {
                    direction: RoutinePortDirection::Output,
                    source: "routine_catalog".to_string(),
                    port: routine_port_to_template_port(&port),
                }),
        );
    } else {
        report.contract_source = "none".to_string();
        report.warnings.push(
            "No typed contract available (template has no ports and no matching routine catalog entry)"
                .to_string(),
        );
    }

    for contract in contracts {
        let port = contract.port;
        let value = resolved_bindings
            .get(&port.port_id)
            .or_else(|| bindings.get(&port.port_id))
            .cloned();
        match value {
            None => {
                if port.required {
                    let direction_label = match contract.direction {
                        RoutinePortDirection::Input => "input",
                        RoutinePortDirection::Output => "output",
                    };
                    let message = format!(
                        "Missing required {} port '{}' (kind={}, cardinality={})",
                        direction_label, port.port_id, port.kind, port.cardinality
                    );
                    report.errors.push(message.clone());
                    report.checked_ports.push(RoutinePortValidationRow {
                        direction: contract.direction,
                        source: contract.source,
                        port_id: port.port_id.clone(),
                        kind: port.kind.clone(),
                        required: port.required,
                        cardinality: port.cardinality.clone(),
                        values: vec![],
                        status: RoutinePortValidationStatus::Missing,
                        message: Some(message),
                    });
                } else {
                    report.checked_ports.push(RoutinePortValidationRow {
                        direction: contract.direction,
                        source: contract.source,
                        port_id: port.port_id.clone(),
                        kind: port.kind.clone(),
                        required: port.required,
                        cardinality: port.cardinality.clone(),
                        values: vec![],
                        status: RoutinePortValidationStatus::Skipped,
                        message: Some("Optional port not bound".to_string()),
                    });
                }
            }
            Some(raw_value) => {
                let values = match split_port_values_for_validation(&raw_value, &port.cardinality) {
                    Ok(values) => values,
                    Err(err) => {
                        let message = format!(
                            "Port '{}' cardinality validation failed: {}",
                            port.port_id, err
                        );
                        report.errors.push(message.clone());
                        report.checked_ports.push(RoutinePortValidationRow {
                            direction: contract.direction,
                            source: contract.source,
                            port_id: port.port_id.clone(),
                            kind: port.kind.clone(),
                            required: port.required,
                            cardinality: port.cardinality.clone(),
                            values: vec![],
                            status: RoutinePortValidationStatus::Invalid,
                            message: Some(message),
                        });
                        continue;
                    }
                };
                if values.is_empty() && port.required {
                    let message =
                        format!("Port '{}' is required and must not be empty", port.port_id);
                    report.errors.push(message.clone());
                    report.checked_ports.push(RoutinePortValidationRow {
                        direction: contract.direction,
                        source: contract.source,
                        port_id: port.port_id.clone(),
                        kind: port.kind.clone(),
                        required: port.required,
                        cardinality: port.cardinality.clone(),
                        values: vec![],
                        status: RoutinePortValidationStatus::Missing,
                        message: Some(message),
                    });
                    continue;
                }
                match validate_port_values(engine, &port.kind, &values, contract.direction) {
                    Ok(_) => {
                        report.checked_ports.push(RoutinePortValidationRow {
                            direction: contract.direction,
                            source: contract.source,
                            port_id: port.port_id.clone(),
                            kind: port.kind.clone(),
                            required: port.required,
                            cardinality: port.cardinality.clone(),
                            values,
                            status: RoutinePortValidationStatus::Ok,
                            message: None,
                        });
                    }
                    Err(err) => {
                        let message = format!("Port '{}' validation failed: {}", port.port_id, err);
                        report.errors.push(message.clone());
                        report.checked_ports.push(RoutinePortValidationRow {
                            direction: contract.direction,
                            source: contract.source,
                            port_id: port.port_id.clone(),
                            kind: port.kind.clone(),
                            required: port.required,
                            cardinality: port.cardinality.clone(),
                            values,
                            status: RoutinePortValidationStatus::Invalid,
                            message: Some(message),
                        });
                    }
                }
            }
        }
    }

    apply_cross_port_semantics(engine, &mut report);
    apply_family_specific_preflight_semantics(engine, &mut report, &resolved_bindings);

    Ok(report)
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
            Self::PanelsImportIsoform {
                seq_id,
                panel_path,
                panel_id,
                strict,
            } => {
                let panel = panel_id
                    .as_deref()
                    .filter(|v| !v.trim().is_empty())
                    .unwrap_or("resource default");
                format!(
                    "import isoform panel from '{panel_path}' for '{seq_id}' (panel_id={panel}, strict={strict})"
                )
            }
            Self::PanelsInspectIsoform { seq_id, panel_id } => {
                format!("inspect isoform architecture for '{seq_id}' panel='{panel_id}'")
            }
            Self::PanelsRenderIsoformSvg {
                seq_id,
                panel_id,
                output,
            } => format!(
                "render isoform architecture SVG for '{seq_id}' panel='{panel_id}' to '{output}'"
            ),
            Self::PanelsValidateIsoform {
                panel_path,
                panel_id,
            } => {
                let panel = panel_id
                    .as_deref()
                    .filter(|v| !v.trim().is_empty())
                    .unwrap_or("resource default");
                format!("validate isoform panel '{panel_path}' (panel_id={panel})")
            }
            Self::UniprotFetch { query, entry_id } => format!(
                "fetch UniProt SWISS-PROT text for '{}' (entry_id={})",
                query,
                entry_id
                    .as_deref()
                    .filter(|v| !v.trim().is_empty())
                    .unwrap_or("auto")
            ),
            Self::GenbankFetch { accession, as_id } => format!(
                "fetch GenBank accession '{}' (as_id={})",
                accession,
                as_id
                    .as_deref()
                    .filter(|v| !v.trim().is_empty())
                    .unwrap_or("auto")
            ),
            Self::UniprotImportSwissProt { path, entry_id } => format!(
                "import UniProt SWISS-PROT text from '{}' (entry_id={})",
                path,
                entry_id
                    .as_deref()
                    .filter(|v| !v.trim().is_empty())
                    .unwrap_or("auto")
            ),
            Self::UniprotList => "list imported UniProt entries".to_string(),
            Self::UniprotShow { entry_id } => format!("show imported UniProt entry '{}'", entry_id),
            Self::UniprotMap {
                entry_id,
                seq_id,
                projection_id,
                transcript_id,
            } => format!(
                "project UniProt entry '{}' onto '{}' (projection_id={}, transcript={})",
                entry_id,
                seq_id,
                projection_id
                    .as_deref()
                    .filter(|v| !v.trim().is_empty())
                    .unwrap_or("auto"),
                transcript_id
                    .as_deref()
                    .filter(|v| !v.trim().is_empty())
                    .unwrap_or("auto"),
            ),
            Self::UniprotProjectionList { seq_id } => format!(
                "list stored UniProt genome projections{}",
                seq_id
                    .as_deref()
                    .filter(|v| !v.trim().is_empty())
                    .map(|v| format!(" for '{}'", v))
                    .unwrap_or_default()
            ),
            Self::UniprotProjectionShow { projection_id } => {
                format!("show UniProt genome projection '{}'", projection_id)
            }
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
            Self::ExportRunBundle { output, run_id } => {
                let run_id = run_id
                    .as_deref()
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .unwrap_or("all");
                format!("export process run bundle to '{output}' (run_id={run_id})")
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
            Self::RoutinesExplain {
                catalog_path,
                routine_id,
            } => {
                let catalog = catalog_path
                    .clone()
                    .unwrap_or_else(|| DEFAULT_CLONING_ROUTINE_CATALOG_PATH.to_string());
                format!("explain cloning routine '{routine_id}' from '{catalog}'")
            }
            Self::RoutinesCompare {
                catalog_path,
                left_routine_id,
                right_routine_id,
            } => {
                let catalog = catalog_path
                    .clone()
                    .unwrap_or_else(|| DEFAULT_CLONING_ROUTINE_CATALOG_PATH.to_string());
                format!(
                    "compare cloning routines '{}' vs '{}' from '{}'",
                    left_routine_id, right_routine_id, catalog
                )
            }
            Self::PlanningProfileShow { scope } => {
                format!("show planning profile (scope={})", scope.as_str())
            }
            Self::PlanningProfileSet {
                scope,
                payload_json,
            } => format!(
                "set planning profile (scope={}, payload_len={})",
                scope.as_str(),
                payload_json.len()
            ),
            Self::PlanningObjectiveShow => "show planning objective".to_string(),
            Self::PlanningObjectiveSet { payload_json } => {
                format!(
                    "set planning objective (payload_len={})",
                    payload_json.len()
                )
            }
            Self::PlanningSuggestionsList { status } => {
                let status = status.map(|value| value.as_str()).unwrap_or("all");
                format!("list planning suggestions (status={status})")
            }
            Self::PlanningSuggestionAccept { suggestion_id } => {
                format!("accept planning suggestion '{suggestion_id}'")
            }
            Self::PlanningSuggestionReject {
                suggestion_id,
                reason,
            } => {
                let reason = reason.as_deref().unwrap_or("-");
                format!("reject planning suggestion '{suggestion_id}' (reason={reason})")
            }
            Self::PlanningSyncStatus => "show planning sync status".to_string(),
            Self::PlanningSyncPull {
                payload_json,
                source,
                confidence,
                snapshot_id,
            } => {
                let source = source.as_deref().unwrap_or("-");
                let confidence = confidence
                    .map(|value| format!("{value:.3}"))
                    .unwrap_or_else(|| "-".to_string());
                let snapshot = snapshot_id.as_deref().unwrap_or("-");
                format!(
                    "register planning sync pull suggestion (source={source}, confidence={confidence}, snapshot={snapshot}, payload_len={})",
                    payload_json.len()
                )
            }
            Self::PlanningSyncPush {
                payload_json,
                source,
                confidence,
                snapshot_id,
            } => {
                let source = source.as_deref().unwrap_or("-");
                let confidence = confidence
                    .map(|value| format!("{value:.3}"))
                    .unwrap_or_else(|| "-".to_string());
                let snapshot = snapshot_id.as_deref().unwrap_or("-");
                format!(
                    "register planning sync push suggestion (source={source}, confidence={confidence}, snapshot={snapshot}, payload_len={})",
                    payload_json.len()
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
                annotation_scope,
                max_annotation_features,
                include_genomic_annotation,
                catalog_path,
                cache_dir,
            } => {
                let label = if *helper_mode { "helper" } else { "genome" };
                let catalog = catalog_path
                    .clone()
                    .unwrap_or_else(|| default_catalog_path(*helper_mode).to_string());
                let cache = cache_dir.clone().unwrap_or_else(|| "-".to_string());
                let output = output_id.clone().unwrap_or_else(|| "-".to_string());
                let scope = annotation_scope
                    .map(|value| value.as_str().to_string())
                    .or_else(|| {
                        include_genomic_annotation.map(|v| {
                            if v {
                                "core(legacy-flag)".to_string()
                            } else {
                                "none(legacy-flag)".to_string()
                            }
                        })
                    })
                    .unwrap_or_else(|| "core(default)".to_string());
                let max_features = max_annotation_features
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "-".to_string());
                let include_annotation = include_genomic_annotation
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "-".to_string());
                format!(
                    "extract {label} region {genome_id}:{chromosome}:{start_1based}-{end_1based} (output='{output}', annotation_scope={scope}, max_annotation_features={max_features}, include_genomic_annotation={include_annotation}, catalog='{catalog}', cache='{cache}')"
                )
            }
            Self::ReferenceExtractGene {
                helper_mode,
                genome_id,
                gene_query,
                occurrence,
                output_id,
                annotation_scope,
                max_annotation_features,
                include_genomic_annotation,
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
                let scope = annotation_scope
                    .map(|value| value.as_str().to_string())
                    .or_else(|| {
                        include_genomic_annotation.map(|v| {
                            if v {
                                "core(legacy-flag)".to_string()
                            } else {
                                "none(legacy-flag)".to_string()
                            }
                        })
                    })
                    .unwrap_or_else(|| "core(default)".to_string());
                let max_features = max_annotation_features
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "-".to_string());
                let include_annotation = include_genomic_annotation
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "-".to_string());
                format!(
                    "extract {label} gene '{gene_query}' from '{genome_id}' (occurrence={occ}, output='{output}', annotation_scope={scope}, max_annotation_features={max_features}, include_genomic_annotation={include_annotation}, catalog='{catalog}', cache='{cache}')"
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
                prepared_genome_id,
            } => {
                let label = if *helper_mode { "helper" } else { "genome" };
                let catalog = catalog_path
                    .clone()
                    .unwrap_or_else(|| default_catalog_path(*helper_mode).to_string());
                let cache = cache_dir.clone().unwrap_or_else(|| "-".to_string());
                let output = output_id.clone().unwrap_or_else(|| "-".to_string());
                let prepared = prepared_genome_id
                    .clone()
                    .unwrap_or_else(|| "-".to_string());
                let side_label = match side {
                    GenomeAnchorSide::FivePrime => "5'",
                    GenomeAnchorSide::ThreePrime => "3'",
                };
                format!(
                    "extend {label}-anchored sequence '{seq_id}' on {side_label} by {length_bp} bp (output='{output}', catalog='{catalog}', cache='{cache}', prepared_genome_id='{prepared}')"
                )
            }
            Self::ReferenceVerifyAnchor {
                helper_mode,
                seq_id,
                catalog_path,
                cache_dir,
                prepared_genome_id,
            } => {
                let label = if *helper_mode { "helper" } else { "genome" };
                let catalog = catalog_path
                    .clone()
                    .unwrap_or_else(|| default_catalog_path(*helper_mode).to_string());
                let cache = cache_dir.clone().unwrap_or_else(|| "-".to_string());
                let prepared = prepared_genome_id
                    .clone()
                    .unwrap_or_else(|| "-".to_string());
                format!(
                    "verify {label}-anchored sequence '{seq_id}' (catalog='{catalog}', cache='{cache}', prepared_genome_id='{prepared}')"
                )
            }
            Self::ReferenceBlast {
                helper_mode,
                genome_id,
                query_sequence,
                max_hits,
                max_hits_explicit,
                task,
                request_options_json,
                catalog_path,
                cache_dir,
            } => {
                let label = if *helper_mode { "helper" } else { "genome" };
                let catalog = catalog_path
                    .clone()
                    .unwrap_or_else(|| default_catalog_path(*helper_mode).to_string());
                let cache = cache_dir.clone().unwrap_or_else(|| "-".to_string());
                let task = task.clone().unwrap_or_else(|| "blastn-short".to_string());
                let max_hits_label = if *max_hits_explicit {
                    max_hits.to_string()
                } else {
                    "layered-default".to_string()
                };
                format!(
                    "blast query (len={}) against {label} '{genome_id}' (max_hits={}, task='{}', request_options={}, catalog='{}', cache='{}')",
                    query_sequence.len(),
                    max_hits_label,
                    task,
                    if request_options_json.is_some() {
                        "yes"
                    } else {
                        "no"
                    },
                    catalog,
                    cache
                )
            }
            Self::ReferenceBlastAsyncStart {
                helper_mode,
                genome_id,
                query_sequence,
                max_hits,
                max_hits_explicit,
                task,
                request_options_json,
                catalog_path,
                cache_dir,
            } => {
                let label = if *helper_mode { "helper" } else { "genome" };
                let catalog = catalog_path
                    .clone()
                    .unwrap_or_else(|| default_catalog_path(*helper_mode).to_string());
                let cache = cache_dir.clone().unwrap_or_else(|| "-".to_string());
                let task = task.clone().unwrap_or_else(|| "blastn-short".to_string());
                let max_hits_label = if *max_hits_explicit {
                    max_hits.to_string()
                } else {
                    "layered-default".to_string()
                };
                format!(
                    "start async blast query (len={}) against {label} '{genome_id}' (max_hits={}, task='{}', request_options={}, catalog='{}', cache='{}')",
                    query_sequence.len(),
                    max_hits_label,
                    task,
                    if request_options_json.is_some() {
                        "yes"
                    } else {
                        "no"
                    },
                    catalog,
                    cache
                )
            }
            Self::ReferenceBlastAsyncStatus {
                helper_mode,
                job_id,
                include_report,
            } => {
                let label = if *helper_mode { "helper" } else { "genome" };
                format!(
                    "inspect async {label} blast job '{}' (include_report={})",
                    job_id, include_report
                )
            }
            Self::ReferenceBlastAsyncCancel {
                helper_mode,
                job_id,
            } => {
                let label = if *helper_mode { "helper" } else { "genome" };
                format!("cancel async {label} blast job '{}'", job_id)
            }
            Self::ReferenceBlastAsyncList { helper_mode } => {
                let label = if *helper_mode { "helper" } else { "genome" };
                format!("list async {label} blast jobs")
            }
            Self::ReferenceBlastTrack {
                helper_mode,
                genome_id,
                query_sequence,
                target_seq_id,
                max_hits,
                max_hits_explicit,
                task,
                request_options_json,
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
                let max_hits_label = if *max_hits_explicit {
                    max_hits.to_string()
                } else {
                    "layered-default".to_string()
                };
                format!(
                    "blast query (len={}) against {label} '{genome_id}' and import hits to '{}' (max_hits={}, task='{}', request_options={}, track_name='{}', clear_existing={}, catalog='{}', cache='{}')",
                    query_sequence.len(),
                    target_seq_id,
                    max_hits_label,
                    task,
                    if request_options_json.is_some() {
                        "yes"
                    } else {
                        "no"
                    },
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
            Self::MacrosInstanceList => "list recorded workflow macro instances".to_string(),
            Self::MacrosInstanceShow { macro_instance_id } => {
                format!(
                    "show recorded workflow macro instance '{}'",
                    macro_instance_id
                )
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
                input_ports,
                output_ports,
                script,
            } => format!(
                "upsert workflow macro template '{}' (description='{}', details_url='{}', params={}, input_ports={}, output_ports={}, script_len={})",
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
                input_ports.len(),
                output_ports.len(),
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
                validate_only,
            } => format!(
                "run workflow macro template '{}' (bindings={}, transactional={}, validate_only={})",
                name,
                bindings
                    .iter()
                    .map(|(key, value)| format!("{key}={value}"))
                    .collect::<Vec<_>>()
                    .join(","),
                transactional,
                validate_only
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
            Self::PrimersDesign {
                request_json,
                backend,
                primer3_executable,
            } => format!(
                "design primer pairs from JSON request payload (len={}, backend='{}', primer3_executable='{}')",
                request_json.len(),
                backend
                    .map(PrimerDesignBackend::as_str)
                    .unwrap_or("default"),
                primer3_executable
                    .as_deref()
                    .map(str::trim)
                    .filter(|v| !v.is_empty())
                    .unwrap_or("default"),
            ),
            Self::PrimersSeedFromFeature { seq_id, feature_id } => format!(
                "seed primer/qPCR design ROI payloads from feature n-{} on '{}'",
                feature_id, seq_id
            ),
            Self::PrimersSeedFromSplicing { seq_id, feature_id } => format!(
                "seed primer/qPCR design ROI payloads from splicing group for feature n-{} on '{}'",
                feature_id, seq_id
            ),
            Self::PrimersDesignQpcr {
                request_json,
                backend,
                primer3_executable,
            } => format!(
                "design qPCR assays (forward/reverse/probe) from JSON request payload (len={}, backend='{}', primer3_executable='{}')",
                request_json.len(),
                backend
                    .map(PrimerDesignBackend::as_str)
                    .unwrap_or("default"),
                primer3_executable
                    .as_deref()
                    .map(str::trim)
                    .filter(|v| !v.is_empty())
                    .unwrap_or("default"),
            ),
            Self::PrimersPreflight {
                backend,
                primer3_executable,
            } => format!(
                "probe Primer3 backend availability/version (backend='{}', primer3_executable='{}')",
                backend
                    .map(PrimerDesignBackend::as_str)
                    .unwrap_or("default"),
                primer3_executable
                    .as_deref()
                    .map(str::trim)
                    .filter(|v| !v.is_empty())
                    .unwrap_or("default"),
            ),
            Self::PrimersListReports => "list stored primer-design reports".to_string(),
            Self::PrimersShowReport { report_id } => {
                format!("show stored primer-design report '{}'", report_id)
            }
            Self::PrimersExportReport { report_id, path } => format!(
                "export stored primer-design report '{}' to '{}'",
                report_id, path
            ),
            Self::PrimersListQpcrReports => "list stored qPCR-design reports".to_string(),
            Self::PrimersShowQpcrReport { report_id } => {
                format!("show stored qPCR-design report '{}'", report_id)
            }
            Self::PrimersExportQpcrReport { report_id, path } => format!(
                "export stored qPCR-design report '{}' to '{}'",
                report_id, path
            ),
            Self::DotplotCompute {
                seq_id,
                span_start_0based,
                span_end_0based,
                mode,
                word_size,
                step_bp,
                max_mismatches,
                tile_bp,
                dotplot_id,
            } => format!(
                "compute dotplot for '{}' (span={}..{}, mode={}, word_size={}, step_bp={}, max_mismatches={}, tile_bp={}, id='{}')",
                seq_id,
                span_start_0based
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "0".to_string()),
                span_end_0based
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "seq_end".to_string()),
                mode.as_str(),
                word_size,
                step_bp,
                max_mismatches,
                tile_bp
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                dotplot_id
                    .as_deref()
                    .filter(|v| !v.trim().is_empty())
                    .unwrap_or("auto")
            ),
            Self::DotplotList { seq_id } => format!(
                "list stored dotplot views{}",
                seq_id
                    .as_deref()
                    .map(|id| format!(" for '{}'", id))
                    .unwrap_or_default()
            ),
            Self::DotplotShow { dotplot_id } => {
                format!("show stored dotplot view '{}'", dotplot_id)
            }
            Self::FlexCompute {
                seq_id,
                span_start_0based,
                span_end_0based,
                model,
                bin_bp,
                smoothing_bp,
                track_id,
            } => format!(
                "compute flexibility track for '{}' (span={}..{}, model={}, bin_bp={}, smoothing_bp={}, id='{}')",
                seq_id,
                span_start_0based
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "0".to_string()),
                span_end_0based
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "seq_end".to_string()),
                model.as_str(),
                bin_bp,
                smoothing_bp
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                track_id
                    .as_deref()
                    .filter(|v| !v.trim().is_empty())
                    .unwrap_or("auto")
            ),
            Self::FlexList { seq_id } => format!(
                "list stored flexibility tracks{}",
                seq_id
                    .as_deref()
                    .map(|id| format!(" for '{}'", id))
                    .unwrap_or_default()
            ),
            Self::FlexShow { track_id } => {
                format!("show stored flexibility track '{}'", track_id)
            }
            Self::RnaReadsInterpret {
                seq_id,
                seed_feature_id,
                input_path,
                profile,
                input_format,
                scope,
                origin_mode,
                target_gene_ids,
                roi_seed_capture_enabled,
                seed_filter,
                align_config,
                report_id,
                report_mode,
                checkpoint_path,
                checkpoint_every_reads,
                resume_from_checkpoint,
            } => format!(
                "interpret RNA reads from '{}' for '{}' feature={} (profile={}, format={}, scope={}, origin_mode={}, report_mode={}, target_genes={}, roi_seed_capture={}, k={}, short_max={}, long_window={}x{}, min_seed_hit_fraction={:.2}, min_weighted_seed_hit_fraction={:.2}, min_unique_matched_kmers={}, min_chain_consistency_fraction={:.2}, max_median_transcript_gap={:.2}, min_confirmed_exon_transitions={}, min_transition_support_fraction={:.2}, cdna_poly_t_flip={}, poly_t_prefix_min_bp={}, align_band={}, align_min_identity={:.2}, max_secondary={}, report_id='{}', checkpoint_path='{}', checkpoint_every_reads={}, resume_from_checkpoint={})",
                input_path,
                seq_id,
                seed_feature_id,
                profile.as_str(),
                input_format.as_str(),
                scope.as_str(),
                origin_mode.as_str(),
                report_mode.as_str(),
                target_gene_ids.len(),
                roi_seed_capture_enabled,
                seed_filter.kmer_len,
                seed_filter.short_full_hash_max_bp,
                seed_filter.long_window_bp,
                seed_filter.long_window_count,
                seed_filter.min_seed_hit_fraction,
                seed_filter.min_weighted_seed_hit_fraction,
                seed_filter.min_unique_matched_kmers,
                seed_filter.min_chain_consistency_fraction,
                seed_filter.max_median_transcript_gap,
                seed_filter.min_confirmed_exon_transitions,
                seed_filter.min_transition_support_fraction,
                seed_filter.cdna_poly_t_flip_enabled,
                seed_filter.poly_t_prefix_min_bp,
                align_config.band_width_bp,
                align_config.min_identity_fraction,
                align_config.max_secondary_mappings,
                report_id
                    .as_deref()
                    .filter(|v| !v.trim().is_empty())
                    .unwrap_or("auto"),
                checkpoint_path
                    .as_deref()
                    .filter(|v| !v.trim().is_empty())
                    .unwrap_or("-"),
                checkpoint_every_reads,
                resume_from_checkpoint
            ),
            Self::RnaReadsListReports { seq_id } => format!(
                "list stored RNA-read reports{}",
                seq_id
                    .as_deref()
                    .map(|id| format!(" for '{}'", id))
                    .unwrap_or_default()
            ),
            Self::RnaReadsShowReport { report_id } => {
                format!("show stored RNA-read report '{}'", report_id)
            }
            Self::RnaReadsExportReport { report_id, path } => format!(
                "export stored RNA-read report '{}' to '{}'",
                report_id, path
            ),
            Self::RnaReadsExportHitsFasta {
                report_id,
                path,
                selection,
            } => format!(
                "export RNA-read hits from '{}' to '{}' (selection={})",
                report_id,
                path,
                selection.as_str()
            ),
            Self::RnaReadsExportSampleSheet {
                path,
                seq_id,
                report_ids,
                append,
            } => format!(
                "export RNA-read sample sheet to '{}' (seq_id='{}', report_ids={}, append={})",
                path,
                seq_id.as_deref().unwrap_or("*"),
                if report_ids.is_empty() {
                    "*".to_string()
                } else {
                    report_ids.join(",")
                },
                append
            ),
            Self::RnaReadsExportExonPathsTsv {
                report_id,
                path,
                selection,
            } => format!(
                "export RNA-read exon paths from '{}' to '{}' (selection={})",
                report_id,
                path,
                selection.as_str()
            ),
            Self::RnaReadsExportExonAbundanceTsv {
                report_id,
                path,
                selection,
            } => format!(
                "export RNA-read exon abundance from '{}' to '{}' (selection={})",
                report_id,
                path,
                selection.as_str()
            ),
            Self::RnaReadsExportScoreDensitySvg {
                report_id,
                path,
                scale,
            } => format!(
                "export RNA-read score-density SVG from '{}' to '{}' (scale={})",
                report_id,
                path,
                scale.as_str()
            ),
            Self::SetParameter { name, value_json } => match name.as_str() {
                "genome_anchor_prepared_fallback_policy"
                | "genome_anchor_fallback_mode"
                | "genome_anchor_prepared_mode" => format!(
                    "set genome-anchor prepared-fallback policy to {} (off|single_compatible|always_explicit)",
                    value_json
                ),
                "require_verified_genome_anchor_for_extension"
                | "strict_genome_anchor_verification"
                | "strict_anchor_verification" => format!(
                    "set strict genome-anchor verification for extension to {}",
                    value_json
                ),
                "primer_design_backend" | "primers_design_backend" => format!(
                    "set primer design backend to {} (auto|internal|primer3)",
                    value_json
                ),
                "primer3_executable" | "primer3_backend_executable" | "primer3_path" => {
                    format!("set primer3 executable path to {}", value_json)
                }
                "linear_sequence_letter_layout_mode" | "linear_helical_letter_layout_mode" => {
                    format!(
                        "set adaptive linear DNA letter mode '{}' (auto|standard|helical|condensed_10_row)",
                        value_json
                    )
                }
                "linear_sequence_helical_letters_enabled"
                | "linear_helical_letters_enabled"
                | "helical_letters_enabled" => format!(
                    "set compressed linear DNA letters '{}' (applies in auto mode)",
                    value_json
                ),
                "linear_show_double_strand_bases"
                | "show_linear_double_strand_bases"
                | "linear_show_reverse_strand_bases"
                | "show_linear_reverse_strand_bases" => format!(
                    "set reverse-strand linear DNA letter visibility to {}",
                    value_json
                ),
                "linear_helical_parallel_strands" | "linear_sequence_helical_parallel_strands" => {
                    format!(
                        "set parallel-vs-mirrored helical strand geometry to {}",
                        value_json
                    )
                }
                "reverse_strand_visual_opacity"
                | "linear_reverse_strand_visual_opacity"
                | "linear_reverse_strand_letter_opacity"
                | "reverse_strand_letter_opacity" => format!(
                    "set reverse-strand letter opacity to {} (0.2..1.0)",
                    value_json
                ),
                "sequence_panel_max_text_length_bp"
                | "sequence_text_panel_max_length_bp"
                | "sequence_panel_max_length_bp" => format!(
                    "set sequence text-panel max length to {} bp (0=unlimited)",
                    value_json
                ),
                "linear_sequence_base_text_max_view_span_bp"
                | "linear_base_text_max_view_span_bp"
                | "sequence_base_text_max_view_span_bp"
                | "linear_sequence_helical_max_view_span_bp"
                | "linear_helical_max_view_span_bp"
                | "linear_sequence_condensed_max_view_span_bp"
                | "linear_condensed_max_view_span_bp" => format!(
                    "set parameter '{}' to {} (deprecated no-op under adaptive routing)",
                    name, value_json
                ),
                _ => format!("set parameter '{}' to {}", name, value_json),
            },
            Self::Op { .. } => "apply one engine operation from JSON".to_string(),
            Self::Workflow { .. } => "apply engine workflow from JSON".to_string(),
        }
    }

    pub fn is_state_mutating(&self) -> bool {
        if let Self::MacrosTemplateRun { validate_only, .. } = self {
            if *validate_only {
                return false;
            }
        }
        matches!(
            self,
            Self::LoadProject { .. }
                | Self::ImportPool { .. }
                | Self::ReferencePrepare { .. }
                | Self::ReferenceExtractRegion { .. }
                | Self::ReferenceExtractGene { .. }
                | Self::ReferenceExtendAnchor { .. }
                | Self::ReferenceVerifyAnchor { .. }
                | Self::ReferenceBlastTrack { .. }
                | Self::TracksImportBed { .. }
                | Self::TracksImportBigWig { .. }
                | Self::TracksImportVcf { .. }
                | Self::TracksTrackedAdd { .. }
                | Self::TracksTrackedRemove { .. }
                | Self::TracksTrackedClear
                | Self::TracksTrackedApply { .. }
                | Self::UniprotFetch { .. }
                | Self::GenbankFetch { .. }
                | Self::UniprotImportSwissProt { .. }
                | Self::UniprotMap { .. }
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
                | Self::PrimersDesign { .. }
                | Self::PrimersDesignQpcr { .. }
                | Self::DotplotCompute { .. }
                | Self::FlexCompute { .. }
                | Self::RnaReadsInterpret { .. }
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

fn shell_now_unix_ms() -> u128 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_millis())
        .unwrap_or(0)
}

fn blast_async_default_max_concurrent_jobs() -> usize {
    thread::available_parallelism()
        .map(|value| value.get())
        .unwrap_or(1)
        .clamp(1, BLAST_ASYNC_MAX_CONCURRENT_HARD_LIMIT)
}

fn blast_async_max_concurrent_jobs() -> usize {
    #[cfg(test)]
    {
        let override_value = BLAST_ASYNC_MAX_CONCURRENT_TEST_OVERRIDE.load(Ordering::Relaxed);
        if override_value > 0 {
            return override_value.clamp(1, BLAST_ASYNC_MAX_CONCURRENT_HARD_LIMIT);
        }
    }
    let default = blast_async_default_max_concurrent_jobs();
    let Some(raw) = std::env::var(BLAST_ASYNC_MAX_CONCURRENT_ENV).ok() else {
        return default;
    };
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return default;
    }
    trimmed
        .parse::<usize>()
        .ok()
        .map(|value| value.clamp(1, BLAST_ASYNC_MAX_CONCURRENT_HARD_LIMIT))
        .unwrap_or(default)
}

fn next_blast_async_job_id() -> String {
    let next = BLAST_ASYNC_JOB_COUNTER.fetch_add(1, Ordering::Relaxed);
    format!("blast-job-{}", next.max(1))
}

fn refresh_blast_async_job_record(record: &mut BlastAsyncJobRecord) {
    let Some(receiver) = &record.receiver else {
        return;
    };
    loop {
        match receiver.try_recv() {
            Ok(BlastAsyncWorkerMessage::Done(result)) => {
                record.status.done_queries = 1;
                record.status.total_queries = 1;
                record.status.finished_at_unix_ms = Some(shell_now_unix_ms());
                match result {
                    Ok(report) => {
                        record.status.state = "completed".to_string();
                        record.status.error = None;
                        record.status.result_available = true;
                        record.report = Some(report);
                    }
                    Err(err) => {
                        record.status.state = if record.cancel_requested.load(Ordering::Relaxed) {
                            "cancelled".to_string()
                        } else {
                            "failed".to_string()
                        };
                        record.status.error = Some(err);
                        record.status.result_available = false;
                        record.report = None;
                    }
                }
                record.receiver = None;
                break;
            }
            Err(mpsc::TryRecvError::Empty) => break,
            Err(mpsc::TryRecvError::Disconnected) => {
                record.status.state = "failed".to_string();
                record.status.finished_at_unix_ms = Some(shell_now_unix_ms());
                record.status.error = Some("BLAST async worker disconnected".to_string());
                record.status.result_available = false;
                record.receiver = None;
                break;
            }
        }
    }
}

fn spawn_blast_async_worker(record: &mut BlastAsyncJobRecord) {
    if record.receiver.is_some() || record.status.state != "queued" {
        return;
    }
    if record.cancel_requested.load(Ordering::Relaxed) {
        record.status.state = "cancelled".to_string();
        record.status.finished_at_unix_ms = Some(shell_now_unix_ms());
        record.status.error = Some("BLAST search cancelled by caller".to_string());
        record.status.result_available = false;
        record.launch_spec = None;
        return;
    }

    let Some(spec) = record.launch_spec.take() else {
        record.status.state = "failed".to_string();
        record.status.finished_at_unix_ms = Some(shell_now_unix_ms());
        record.status.error = Some("BLAST async worker launch payload missing".to_string());
        record.status.result_available = false;
        return;
    };

    let cancel_for_thread = Arc::clone(&record.cancel_requested);
    let (tx, rx) = mpsc::channel::<BlastAsyncWorkerMessage>();
    let started_at = shell_now_unix_ms();
    record.status.state = "running".to_string();
    record.status.started_at_unix_ms = Some(started_at);
    record.receiver = Some(rx);

    thread::spawn(move || {
        #[cfg(test)]
        {
            let delay_ms = BLAST_ASYNC_WORKER_DELAY_MS_TEST_OVERRIDE.load(Ordering::Relaxed);
            if delay_ms > 0 {
                thread::sleep(Duration::from_millis(delay_ms));
            }
        }
        let mut should_cancel = || cancel_for_thread.load(Ordering::Relaxed);
        let result = if spec.helper_mode {
            spec.engine_snapshot
                .blast_helper_genome_with_project_and_request_options_and_cancel(
                    &spec.genome_id,
                    &spec.query_sequence,
                    spec.request_options_json.as_ref(),
                    spec.task.as_deref(),
                    if spec.max_hits_explicit {
                        Some(spec.max_hits)
                    } else {
                        None
                    },
                    spec.resolved_catalog.as_deref(),
                    spec.cache_dir.as_deref(),
                    &mut should_cancel,
                )
        } else {
            spec.engine_snapshot
                .blast_reference_genome_with_project_and_request_options_and_cancel(
                    spec.resolved_catalog.as_deref(),
                    &spec.genome_id,
                    &spec.query_sequence,
                    spec.request_options_json.as_ref(),
                    spec.task.as_deref(),
                    if spec.max_hits_explicit {
                        Some(spec.max_hits)
                    } else {
                        None
                    },
                    spec.cache_dir.as_deref(),
                    &mut should_cancel,
                )
        };
        let _ = tx.send(BlastAsyncWorkerMessage::Done(
            result.map_err(|e| e.to_string()),
        ));
        kick_blast_async_scheduler();
    });
}

fn refresh_and_dispatch_blast_async_jobs_locked(jobs: &mut HashMap<String, BlastAsyncJobRecord>) {
    for record in jobs.values_mut() {
        refresh_blast_async_job_record(record);
    }

    let max_concurrent = blast_async_max_concurrent_jobs();
    let mut running_jobs = jobs
        .values()
        .filter(|record| record.status.state == "running")
        .count();
    if running_jobs < max_concurrent {
        let mut queued_job_ids: Vec<(u128, String)> = jobs
            .iter()
            .filter_map(|(job_id, record)| {
                if record.status.state == "queued" {
                    Some((record.status.created_at_unix_ms, job_id.clone()))
                } else {
                    None
                }
            })
            .collect();
        queued_job_ids.sort_by(|left, right| left.0.cmp(&right.0).then(left.1.cmp(&right.1)));
        for (_, job_id) in queued_job_ids {
            if running_jobs >= max_concurrent {
                break;
            }
            let Some(record) = jobs.get_mut(&job_id) else {
                continue;
            };
            if record.status.state != "queued" {
                continue;
            }
            spawn_blast_async_worker(record);
            if record.status.state == "running" {
                running_jobs += 1;
            }
        }
    }

    let running_jobs = jobs
        .values()
        .filter(|record| record.status.state == "running")
        .count();
    let mut queued_job_ids: Vec<(u128, String)> = jobs
        .iter()
        .filter_map(|(job_id, record)| {
            if record.status.state == "queued" {
                Some((record.status.created_at_unix_ms, job_id.clone()))
            } else {
                None
            }
        })
        .collect();
    queued_job_ids.sort_by(|left, right| left.0.cmp(&right.0).then(left.1.cmp(&right.1)));
    let queued_jobs = queued_job_ids.len();
    let mut queue_positions: HashMap<String, usize> = HashMap::new();
    for (index, (_, job_id)) in queued_job_ids.into_iter().enumerate() {
        queue_positions.insert(job_id, index + 1);
    }
    for (job_id, record) in jobs.iter_mut() {
        record.status.max_concurrent_jobs = max_concurrent;
        record.status.running_jobs = running_jobs;
        record.status.queued_jobs = queued_jobs;
        record.status.queue_position = queue_positions.get(job_id).copied();
    }
}

fn kick_blast_async_scheduler() {
    let Ok(mut jobs) = BLAST_ASYNC_JOBS.lock() else {
        return;
    };
    refresh_and_dispatch_blast_async_jobs_locked(&mut jobs);
    if jobs.len() > BLAST_ASYNC_JOB_HISTORY_LIMIT {
        prune_blast_async_jobs_locked(&mut jobs);
    }
}

fn prune_blast_async_jobs_locked(jobs: &mut HashMap<String, BlastAsyncJobRecord>) {
    refresh_and_dispatch_blast_async_jobs_locked(jobs);
    if jobs.len() <= BLAST_ASYNC_JOB_HISTORY_LIMIT {
        return;
    }
    let mut removable: Vec<(u128, String)> = jobs
        .iter()
        .filter_map(|(job_id, record)| match record.status.state.as_str() {
            "completed" | "failed" | "cancelled" => {
                Some((record.status.created_at_unix_ms, job_id.clone()))
            }
            _ => None,
        })
        .collect();
    removable.sort_by_key(|(created_at, _)| *created_at);
    let excess = jobs.len().saturating_sub(BLAST_ASYNC_JOB_HISTORY_LIMIT);
    for (_, job_id) in removable.into_iter().take(excess) {
        jobs.remove(&job_id);
    }
}

fn collect_blast_async_job_snapshots(
    helper_mode_filter: Option<bool>,
) -> Result<Vec<BlastAsyncJobStatus>, String> {
    let mut jobs = BLAST_ASYNC_JOBS
        .lock()
        .map_err(|_| "Could not lock BLAST async job registry".to_string())?;
    refresh_and_dispatch_blast_async_jobs_locked(&mut jobs);
    let mut statuses: Vec<BlastAsyncJobStatus> = vec![];
    for record in jobs.values_mut() {
        if let Some(helper_mode) = helper_mode_filter {
            if record.status.helper_mode != helper_mode {
                continue;
            }
        }
        statuses.push(record.status.clone());
    }
    statuses.sort_by(|a, b| {
        a.created_at_unix_ms
            .cmp(&b.created_at_unix_ms)
            .then(a.job_id.cmp(&b.job_id))
    });
    Ok(statuses)
}

fn get_blast_async_job_snapshot(
    job_id: &str,
    include_report: bool,
) -> Result<(BlastAsyncJobStatus, Option<GenomeBlastReport>), String> {
    let mut jobs = BLAST_ASYNC_JOBS
        .lock()
        .map_err(|_| "Could not lock BLAST async job registry".to_string())?;
    refresh_and_dispatch_blast_async_jobs_locked(&mut jobs);
    let record = jobs
        .get_mut(job_id)
        .ok_or_else(|| format!("BLAST async job '{}' not found", job_id))?;
    Ok((
        record.status.clone(),
        if include_report {
            record.report.clone()
        } else {
            None
        },
    ))
}

fn cancel_blast_async_job(job_id: &str) -> Result<BlastAsyncJobStatus, String> {
    let mut jobs = BLAST_ASYNC_JOBS
        .lock()
        .map_err(|_| "Could not lock BLAST async job registry".to_string())?;
    let record = jobs
        .get_mut(job_id)
        .ok_or_else(|| format!("BLAST async job '{}' not found", job_id))?;
    record.cancel_requested.store(true, Ordering::Relaxed);
    record.status.cancel_requested = true;
    if record.receiver.is_none() {
        if record.status.state == "queued" || record.status.state == "running" {
            record.status.state = "cancelled".to_string();
            record.status.finished_at_unix_ms = Some(shell_now_unix_ms());
            record.status.error = Some("BLAST search cancelled by caller".to_string());
            record.status.result_available = false;
            record.launch_spec = None;
        }
    } else if record.status.started_at_unix_ms.is_none() {
        record.status.started_at_unix_ms = Some(shell_now_unix_ms());
    }
    refresh_and_dispatch_blast_async_jobs_locked(&mut jobs);
    let record = jobs
        .get(job_id)
        .ok_or_else(|| format!("BLAST async job '{}' not found", job_id))?;
    Ok(record.status.clone())
}

#[allow(clippy::too_many_arguments)]
fn start_blast_async_job(
    engine: &GentleEngine,
    helper_mode: bool,
    genome_id: &str,
    query_sequence: &str,
    max_hits: usize,
    max_hits_explicit: bool,
    task: Option<String>,
    request_options_json: Option<Value>,
    catalog_path: Option<String>,
    cache_dir: Option<String>,
) -> Result<BlastAsyncJobStatus, String> {
    if genome_id.trim().is_empty() {
        return Err("BLAST async start requires a non-empty genome_id".to_string());
    }
    if query_sequence.trim().is_empty() {
        return Err("BLAST async start requires a non-empty query_sequence".to_string());
    }
    if max_hits == 0 {
        return Err("BLAST async start requires max_hits >= 1".to_string());
    }

    let resolved_catalog = resolved_catalog_path(&catalog_path, helper_mode).map(str::to_string);
    let job_id = next_blast_async_job_id();
    let created_at = shell_now_unix_ms();
    let cancel_requested = Arc::new(AtomicBool::new(false));
    let task_name = task.clone().unwrap_or_else(|| "blastn-short".to_string());
    let launch_spec = BlastAsyncLaunchSpec {
        engine_snapshot: engine.clone(),
        helper_mode,
        genome_id: genome_id.to_string(),
        query_sequence: query_sequence.to_string(),
        max_hits,
        max_hits_explicit,
        task: task.clone(),
        request_options_json,
        resolved_catalog,
        cache_dir,
    };

    let status = BlastAsyncJobStatus {
        schema: BLAST_ASYNC_JOB_SCHEMA.to_string(),
        job_id: job_id.clone(),
        state: "queued".to_string(),
        helper_mode,
        genome_id: genome_id.to_string(),
        query_length: query_sequence.len(),
        task: task_name,
        max_hits,
        created_at_unix_ms: created_at,
        started_at_unix_ms: None,
        finished_at_unix_ms: None,
        cancel_requested: false,
        done_queries: 0,
        total_queries: 1,
        result_available: false,
        error: None,
        max_concurrent_jobs: 0,
        running_jobs: 0,
        queued_jobs: 0,
        queue_position: None,
    };
    let mut jobs = BLAST_ASYNC_JOBS
        .lock()
        .map_err(|_| "Could not lock BLAST async job registry".to_string())?;
    jobs.insert(
        job_id,
        BlastAsyncJobRecord {
            status: status.clone(),
            cancel_requested,
            receiver: None,
            launch_spec: Some(launch_spec),
            report: None,
        },
    );
    refresh_and_dispatch_blast_async_jobs_locked(&mut jobs);
    let returned_status = jobs
        .get(&status.job_id)
        .map(|record| record.status.clone())
        .unwrap_or(status);
    prune_blast_async_jobs_locked(&mut jobs);
    Ok(returned_status)
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
            "{context} requires target syntax: tfbs FEATURE_ID | restriction CUT_POS_1BASED [--enzyme NAME] [--start START_1BASED] [--end END_1BASED] | splicing FEATURE_ID | isoform PANEL_ID"
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
            Ok(FeatureExpertTarget::SplicingFeature {
                feature_id,
                scope: SplicingScopePreset::AllOverlappingBothStrands,
            })
        }
        "isoform" => {
            if tokens.len() != 2 {
                return Err(format!(
                    "{context} isoform target expects exactly: isoform PANEL_ID"
                ));
            }
            let panel_id = tokens[1].trim().to_string();
            if panel_id.is_empty() {
                return Err("isoform PANEL_ID must not be empty".to_string());
            }
            Ok(FeatureExpertTarget::IsoformArchitecture { panel_id })
        }
        other => Err(format!(
            "Unknown feature target '{other}' (expected tfbs|restriction|splicing|isoform)"
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

fn parse_primer_design_backend(value: &str) -> Result<PrimerDesignBackend, String> {
    match value.trim().to_ascii_lowercase().as_str() {
        "auto" => Ok(PrimerDesignBackend::Auto),
        "internal" | "baseline" => Ok(PrimerDesignBackend::Internal),
        "primer3" | "primer3_core" | "external-primer3" | "external_primer3" => {
            Ok(PrimerDesignBackend::Primer3)
        }
        other => Err(format!(
            "Unsupported primer backend '{}' (expected auto|internal|primer3)",
            other
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

fn parse_optional_json_payload<T>(raw: &str, context: &str) -> Result<Option<T>, String>
where
    T: DeserializeOwned,
{
    let loaded = parse_json_payload(raw)?;
    let trimmed = loaded.trim();
    if trimmed.is_empty() || trimmed.eq_ignore_ascii_case("null") {
        return Ok(None);
    }
    serde_json::from_str::<T>(trimmed)
        .map(Some)
        .map_err(|e| format!("Invalid {context} JSON payload: {e}"))
}

fn parse_required_json_payload<T>(raw: &str, context: &str) -> Result<T, String>
where
    T: DeserializeOwned,
{
    let loaded = parse_json_payload(raw)?;
    serde_json::from_str::<T>(loaded.trim())
        .map_err(|e| format!("Invalid {context} JSON payload: {e}"))
}

#[derive(Debug, Clone, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
struct PlanningSyncSuggestionPayload {
    profile_patch: Option<PlanningProfile>,
    objective_patch: Option<PlanningObjective>,
    message: Option<String>,
}

#[derive(Debug, Clone, Deserialize)]
struct WorkflowPayloadWrapper {
    workflow: Workflow,
}

/// Parse workflow payload JSON accepted by direct and shell adapters.
///
/// Accepted shapes:
/// - Raw workflow object: `{ "run_id": "...", "ops": [...] }`
/// - Wrapped example object: `{ ..., "workflow": { "run_id": "...", "ops": [...] } }`
pub fn parse_workflow_json_payload(raw_json: &str) -> Result<Workflow, String> {
    match serde_json::from_str::<Workflow>(raw_json) {
        Ok(workflow) => Ok(workflow),
        Err(primary_error) => {
            let wrapped: WorkflowPayloadWrapper = serde_json::from_str(raw_json)
                .map_err(|_| format!("Invalid workflow JSON: {primary_error}"))?;
            Ok(wrapped.workflow)
        }
    }
}

fn parse_blast_options_override(raw: &str, label: &str, flag: &str) -> Result<Value, String> {
    let loaded = parse_json_payload(raw)?;
    let parsed: Value = serde_json::from_str(&loaded)
        .map_err(|e| format!("Invalid {flag} payload for {label} blast command: {e}"))?;
    if !parsed.is_object() {
        return Err(format!(
            "{flag} payload for {label} blast command must decode to a JSON object"
        ));
    }
    Ok(parsed)
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
                    "{label} blast requires GENOME_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--options-json JSON_OR_@FILE | --options-file PATH] [--catalog PATH] [--cache-dir PATH]"
                ));
            }
            let genome_id = tokens[2].clone();
            let query_sequence = tokens[3].clone();
            let mut max_hits: usize = 25;
            let mut max_hits_explicit = false;
            let mut task: Option<String> = None;
            let mut request_options_json: Option<Value> = None;
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
                        max_hits_explicit = true;
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
                    "--options-json" => {
                        if request_options_json.is_some() {
                            return Err(format!(
                                "Only one of --options-json/--options-file may be provided for {label} blast"
                            ));
                        }
                        let raw = parse_option_path(tokens, &mut idx, "--options-json", label)?;
                        request_options_json =
                            Some(parse_blast_options_override(&raw, label, "--options-json")?);
                    }
                    "--options-file" => {
                        if request_options_json.is_some() {
                            return Err(format!(
                                "Only one of --options-json/--options-file may be provided for {label} blast"
                            ));
                        }
                        let raw = parse_option_path(tokens, &mut idx, "--options-file", label)?;
                        let file_ref = format!("@{raw}");
                        request_options_json = Some(parse_blast_options_override(
                            &file_ref,
                            label,
                            "--options-file",
                        )?);
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
                max_hits_explicit,
                task,
                request_options_json,
                catalog_path,
                cache_dir,
            })
        }
        "blast-start" => {
            if tokens.len() < 4 {
                return Err(format!(
                    "{label} blast-start requires GENOME_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--options-json JSON_OR_@FILE | --options-file PATH] [--catalog PATH] [--cache-dir PATH]"
                ));
            }
            let genome_id = tokens[2].clone();
            let query_sequence = tokens[3].clone();
            let mut max_hits: usize = 25;
            let mut max_hits_explicit = false;
            let mut task: Option<String> = None;
            let mut request_options_json: Option<Value> = None;
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
                        max_hits_explicit = true;
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
                    "--options-json" => {
                        if request_options_json.is_some() {
                            return Err(format!(
                                "Only one of --options-json/--options-file may be provided for {label} blast-start"
                            ));
                        }
                        let raw = parse_option_path(tokens, &mut idx, "--options-json", label)?;
                        request_options_json =
                            Some(parse_blast_options_override(&raw, label, "--options-json")?);
                    }
                    "--options-file" => {
                        if request_options_json.is_some() {
                            return Err(format!(
                                "Only one of --options-json/--options-file may be provided for {label} blast-start"
                            ));
                        }
                        let raw = parse_option_path(tokens, &mut idx, "--options-file", label)?;
                        let file_ref = format!("@{raw}");
                        request_options_json = Some(parse_blast_options_override(
                            &file_ref,
                            label,
                            "--options-file",
                        )?);
                    }
                    "--catalog" => {
                        catalog_path =
                            Some(parse_option_path(tokens, &mut idx, "--catalog", label)?)
                    }
                    "--cache-dir" => {
                        cache_dir = Some(parse_option_path(tokens, &mut idx, "--cache-dir", label)?)
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for {label} blast-start"));
                    }
                }
            }
            Ok(ShellCommand::ReferenceBlastAsyncStart {
                helper_mode,
                genome_id,
                query_sequence,
                max_hits,
                max_hits_explicit,
                task,
                request_options_json,
                catalog_path,
                cache_dir,
            })
        }
        "blast-status" => {
            if tokens.len() < 3 {
                return Err(format!(
                    "{label} blast-status requires JOB_ID [--with-report]"
                ));
            }
            let job_id = tokens[2].trim().to_string();
            if job_id.is_empty() {
                return Err(format!("{label} blast-status requires non-empty JOB_ID"));
            }
            let mut include_report = false;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--with-report" => {
                        include_report = true;
                        idx += 1;
                    }
                    "--no-report" => {
                        include_report = false;
                        idx += 1;
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for {label} blast-status"));
                    }
                }
            }
            Ok(ShellCommand::ReferenceBlastAsyncStatus {
                helper_mode,
                job_id,
                include_report,
            })
        }
        "blast-cancel" => {
            if tokens.len() != 3 {
                return Err(format!("{label} blast-cancel requires JOB_ID"));
            }
            let job_id = tokens[2].trim().to_string();
            if job_id.is_empty() {
                return Err(format!("{label} blast-cancel requires non-empty JOB_ID"));
            }
            Ok(ShellCommand::ReferenceBlastAsyncCancel {
                helper_mode,
                job_id,
            })
        }
        "blast-list" => {
            if tokens.len() != 2 {
                return Err(format!("{label} blast-list takes no options"));
            }
            Ok(ShellCommand::ReferenceBlastAsyncList { helper_mode })
        }
        "blast-track" => {
            if tokens.len() < 5 {
                return Err(format!(
                    "{label} blast-track requires GENOME_ID QUERY_SEQUENCE TARGET_SEQ_ID [--max-hits N] [--task blastn-short|blastn] [--options-json JSON_OR_@FILE | --options-file PATH] [--track-name NAME] [--clear-existing] [--catalog PATH] [--cache-dir PATH]"
                ));
            }
            let genome_id = tokens[2].clone();
            let query_sequence = tokens[3].clone();
            let target_seq_id = tokens[4].clone();
            let mut max_hits: usize = 25;
            let mut max_hits_explicit = false;
            let mut task: Option<String> = None;
            let mut request_options_json: Option<Value> = None;
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
                        max_hits_explicit = true;
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
                    "--options-json" => {
                        if request_options_json.is_some() {
                            return Err(format!(
                                "Only one of --options-json/--options-file may be provided for {label} blast-track"
                            ));
                        }
                        let raw = parse_option_path(tokens, &mut idx, "--options-json", label)?;
                        request_options_json =
                            Some(parse_blast_options_override(&raw, label, "--options-json")?);
                    }
                    "--options-file" => {
                        if request_options_json.is_some() {
                            return Err(format!(
                                "Only one of --options-json/--options-file may be provided for {label} blast-track"
                            ));
                        }
                        let raw = parse_option_path(tokens, &mut idx, "--options-file", label)?;
                        let file_ref = format!("@{raw}");
                        request_options_json = Some(parse_blast_options_override(
                            &file_ref,
                            label,
                            "--options-file",
                        )?);
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
                max_hits_explicit,
                task,
                request_options_json,
                track_name,
                clear_existing,
                catalog_path,
                cache_dir,
            })
        }
        "extract-region" => {
            if tokens.len() < 6 {
                return Err(format!(
                    "{label} extract-region requires GENOME_ID CHR START END [--output-id ID] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--catalog PATH] [--cache-dir PATH]"
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
            let mut annotation_scope: Option<GenomeAnnotationScope> = None;
            let mut max_annotation_features: Option<usize> = None;
            let mut include_genomic_annotation: Option<bool> = None;
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
                    "--annotation-scope" => {
                        let raw = parse_option_path(tokens, &mut idx, "--annotation-scope", label)?;
                        let scope = match raw.trim().to_ascii_lowercase().as_str() {
                            "none" => GenomeAnnotationScope::None,
                            "core" => GenomeAnnotationScope::Core,
                            "full" => GenomeAnnotationScope::Full,
                            other => {
                                return Err(format!(
                                    "Invalid --annotation-scope value '{other}' for {label} extract-region (expected none|core|full)"
                                ));
                            }
                        };
                        annotation_scope = Some(scope);
                    }
                    "--max-annotation-features" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--max-annotation-features",
                            label,
                        )?;
                        let parsed = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --max-annotation-features value '{raw}' for {label} extract-region: {e}"
                            )
                        })?;
                        max_annotation_features = Some(parsed);
                    }
                    "--include-genomic-annotation" => {
                        include_genomic_annotation = Some(true);
                        idx += 1;
                    }
                    "--no-include-genomic-annotation" => {
                        include_genomic_annotation = Some(false);
                        idx += 1;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for {label} extract-region"
                        ));
                    }
                }
            }
            if let Some(include) = include_genomic_annotation {
                let mapped_scope = if include {
                    GenomeAnnotationScope::Core
                } else {
                    GenomeAnnotationScope::None
                };
                if let Some(explicit_scope) = annotation_scope {
                    if explicit_scope != mapped_scope {
                        return Err(format!(
                            "Conflicting annotation options for {label} extract-region: --annotation-scope={} with legacy include/no-include flag",
                            explicit_scope.as_str()
                        ));
                    }
                } else {
                    annotation_scope = Some(mapped_scope);
                }
            }
            Ok(ShellCommand::ReferenceExtractRegion {
                helper_mode,
                genome_id,
                chromosome,
                start_1based,
                end_1based,
                output_id,
                annotation_scope,
                max_annotation_features,
                include_genomic_annotation,
                catalog_path,
                cache_dir,
            })
        }
        "extract-gene" => {
            if tokens.len() < 4 {
                return Err(format!(
                    "{label} extract-gene requires GENOME_ID QUERY [--occurrence N] [--output-id ID] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--catalog PATH] [--cache-dir PATH]"
                ));
            }
            let genome_id = tokens[2].clone();
            let gene_query = tokens[3].clone();
            let mut occurrence: Option<usize> = None;
            let mut output_id: Option<String> = None;
            let mut annotation_scope: Option<GenomeAnnotationScope> = None;
            let mut max_annotation_features: Option<usize> = None;
            let mut include_genomic_annotation: Option<bool> = None;
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
                    "--annotation-scope" => {
                        let raw = parse_option_path(tokens, &mut idx, "--annotation-scope", label)?;
                        let scope = match raw.trim().to_ascii_lowercase().as_str() {
                            "none" => GenomeAnnotationScope::None,
                            "core" => GenomeAnnotationScope::Core,
                            "full" => GenomeAnnotationScope::Full,
                            other => {
                                return Err(format!(
                                    "Invalid --annotation-scope value '{other}' for {label} extract-gene (expected none|core|full)"
                                ));
                            }
                        };
                        annotation_scope = Some(scope);
                    }
                    "--max-annotation-features" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--max-annotation-features",
                            label,
                        )?;
                        let parsed = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --max-annotation-features value '{raw}' for {label} extract-gene: {e}"
                            )
                        })?;
                        max_annotation_features = Some(parsed);
                    }
                    "--include-genomic-annotation" => {
                        include_genomic_annotation = Some(true);
                        idx += 1;
                    }
                    "--no-include-genomic-annotation" => {
                        include_genomic_annotation = Some(false);
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
                        return Err(format!("Unknown option '{other}' for {label} extract-gene"));
                    }
                }
            }
            if let Some(include) = include_genomic_annotation {
                let mapped_scope = if include {
                    GenomeAnnotationScope::Core
                } else {
                    GenomeAnnotationScope::None
                };
                if let Some(explicit_scope) = annotation_scope {
                    if explicit_scope != mapped_scope {
                        return Err(format!(
                            "Conflicting annotation options for {label} extract-gene: --annotation-scope={} with legacy include/no-include flag",
                            explicit_scope.as_str()
                        ));
                    }
                } else {
                    annotation_scope = Some(mapped_scope);
                }
            }
            Ok(ShellCommand::ReferenceExtractGene {
                helper_mode,
                genome_id,
                gene_query,
                occurrence,
                output_id,
                annotation_scope,
                max_annotation_features,
                include_genomic_annotation,
                catalog_path,
                cache_dir,
            })
        }
        "extend-anchor" => {
            if tokens.len() < 5 {
                return Err(format!(
                    "{label} extend-anchor requires SEQ_ID 5p|3p LENGTH_BP [--output-id ID] [--catalog PATH] [--cache-dir PATH] [--prepared-genome GENOME_ID]"
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
            let mut prepared_genome_id: Option<String> = None;
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
                    "--prepared-genome" => {
                        prepared_genome_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--prepared-genome",
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
                prepared_genome_id,
            })
        }
        "verify-anchor" => {
            if tokens.len() < 3 {
                return Err(format!(
                    "{label} verify-anchor requires SEQ_ID [--catalog PATH] [--cache-dir PATH] [--prepared-genome GENOME_ID]"
                ));
            }
            let seq_id = tokens[2].clone();
            let mut catalog_path: Option<String> = None;
            let mut cache_dir: Option<String> = None;
            let mut prepared_genome_id: Option<String> = None;
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
                    "--prepared-genome" => {
                        prepared_genome_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--prepared-genome",
                            label,
                        )?)
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for {label} verify-anchor"
                        ));
                    }
                }
            }
            Ok(ShellCommand::ReferenceVerifyAnchor {
                helper_mode,
                seq_id,
                catalog_path,
                cache_dir,
                prepared_genome_id,
            })
        }
        other => Err(format!(
            "Unknown {label} subcommand '{other}' (expected list, validate-catalog, status, genes, prepare, blast, blast-start, blast-status, blast-cancel, blast-list, blast-track, extract-region, extract-gene, extend-anchor, verify-anchor)"
        )),
    }
}

fn parse_panels_command(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.len() < 2 {
        return Err(
            "panels requires a subcommand: import-isoform, inspect-isoform, render-isoform-svg, validate-isoform"
                .to_string(),
        );
    }
    match tokens[1].as_str() {
        "import-isoform" | "import" => {
            if tokens.len() < 4 {
                return Err(
                    "panels import-isoform requires SEQ_ID PANEL_PATH [--panel-id ID] [--strict]"
                        .to_string(),
                );
            }
            let seq_id = tokens[2].clone();
            let panel_path = tokens[3].clone();
            let mut panel_id: Option<String> = None;
            let mut strict = false;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--panel-id" => {
                        panel_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--panel-id",
                            "panels import-isoform",
                        )?);
                    }
                    "--strict" => {
                        strict = true;
                        idx += 1;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for panels import-isoform"
                        ));
                    }
                }
            }
            Ok(ShellCommand::PanelsImportIsoform {
                seq_id,
                panel_path,
                panel_id,
                strict,
            })
        }
        "inspect-isoform" | "inspect" => {
            if tokens.len() != 4 {
                return Err("panels inspect-isoform requires SEQ_ID PANEL_ID".to_string());
            }
            let seq_id = tokens[2].clone();
            let panel_id = tokens[3].trim().to_string();
            if panel_id.is_empty() {
                return Err("panels inspect-isoform PANEL_ID must not be empty".to_string());
            }
            Ok(ShellCommand::PanelsInspectIsoform { seq_id, panel_id })
        }
        "render-isoform-svg" | "render-svg" => {
            if tokens.len() != 5 {
                return Err(
                    "panels render-isoform-svg requires SEQ_ID PANEL_ID OUTPUT.svg".to_string(),
                );
            }
            let seq_id = tokens[2].clone();
            let panel_id = tokens[3].trim().to_string();
            if panel_id.is_empty() {
                return Err("panels render-isoform-svg PANEL_ID must not be empty".to_string());
            }
            let output = tokens[4].clone();
            Ok(ShellCommand::PanelsRenderIsoformSvg {
                seq_id,
                panel_id,
                output,
            })
        }
        "validate-isoform" | "validate" => {
            if tokens.len() < 3 {
                return Err(
                    "panels validate-isoform requires PANEL_PATH [--panel-id ID]".to_string(),
                );
            }
            let panel_path = tokens[2].clone();
            let mut panel_id: Option<String> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--panel-id" => {
                        panel_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--panel-id",
                            "panels validate-isoform",
                        )?);
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for panels validate-isoform"
                        ));
                    }
                }
            }
            Ok(ShellCommand::PanelsValidateIsoform {
                panel_path,
                panel_id,
            })
        }
        other => Err(format!(
            "Unknown panels subcommand '{other}' (expected import-isoform, inspect-isoform, render-isoform-svg, validate-isoform)"
        )),
    }
}

fn parse_genbank_command(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.len() < 2 {
        return Err("genbank requires a subcommand: fetch".to_string());
    }
    match tokens[1].as_str() {
        "fetch" => {
            if tokens.len() < 3 {
                return Err("genbank fetch requires ACCESSION [--as-id ID]".to_string());
            }
            let accession = tokens[2].trim().to_string();
            if accession.is_empty() {
                return Err("genbank fetch ACCESSION must not be empty".to_string());
            }
            let mut as_id: Option<String> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--as-id" => {
                        as_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--as-id",
                            "genbank fetch",
                        )?);
                    }
                    other => return Err(format!("Unknown option '{other}' for genbank fetch")),
                }
            }
            Ok(ShellCommand::GenbankFetch { accession, as_id })
        }
        other => Err(format!(
            "Unknown genbank subcommand '{other}' (expected fetch)"
        )),
    }
}

fn parse_uniprot_command(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.len() < 2 {
        return Err(
            "uniprot requires a subcommand: fetch, import-swissprot, list, show, map, projection-list, projection-show"
                .to_string(),
        );
    }
    match tokens[1].as_str() {
        "fetch" => {
            if tokens.len() < 3 {
                return Err("uniprot fetch requires QUERY [--entry-id ID]".to_string());
            }
            let query = tokens[2].trim().to_string();
            if query.is_empty() {
                return Err("uniprot fetch QUERY must not be empty".to_string());
            }
            let mut entry_id: Option<String> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--entry-id" => {
                        entry_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--entry-id",
                            "uniprot fetch",
                        )?);
                    }
                    other => return Err(format!("Unknown option '{other}' for uniprot fetch")),
                }
            }
            Ok(ShellCommand::UniprotFetch { query, entry_id })
        }
        "import-swissprot" | "import" => {
            if tokens.len() < 3 {
                return Err("uniprot import-swissprot requires PATH [--entry-id ID]".to_string());
            }
            let path = tokens[2].trim().to_string();
            if path.is_empty() {
                return Err("uniprot import-swissprot PATH must not be empty".to_string());
            }
            let mut entry_id: Option<String> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--entry-id" => {
                        entry_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--entry-id",
                            "uniprot import-swissprot",
                        )?);
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for uniprot import-swissprot"
                        ));
                    }
                }
            }
            Ok(ShellCommand::UniprotImportSwissProt { path, entry_id })
        }
        "list" => {
            if tokens.len() > 2 {
                return Err("uniprot list takes no options".to_string());
            }
            Ok(ShellCommand::UniprotList)
        }
        "show" => {
            if tokens.len() != 3 {
                return Err("uniprot show requires ENTRY_ID".to_string());
            }
            let entry_id = tokens[2].trim().to_string();
            if entry_id.is_empty() {
                return Err("uniprot show ENTRY_ID must not be empty".to_string());
            }
            Ok(ShellCommand::UniprotShow { entry_id })
        }
        "map" => {
            if tokens.len() < 4 {
                return Err(
                    "uniprot map requires ENTRY_ID SEQ_ID [--projection-id ID] [--transcript ID]"
                        .to_string(),
                );
            }
            let entry_id = tokens[2].trim().to_string();
            let seq_id = tokens[3].trim().to_string();
            if entry_id.is_empty() {
                return Err("uniprot map ENTRY_ID must not be empty".to_string());
            }
            if seq_id.is_empty() {
                return Err("uniprot map SEQ_ID must not be empty".to_string());
            }
            let mut projection_id: Option<String> = None;
            let mut transcript_id: Option<String> = None;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--projection-id" => {
                        projection_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--projection-id",
                            "uniprot map",
                        )?);
                    }
                    "--transcript" => {
                        transcript_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--transcript",
                            "uniprot map",
                        )?);
                    }
                    other => return Err(format!("Unknown option '{other}' for uniprot map")),
                }
            }
            Ok(ShellCommand::UniprotMap {
                entry_id,
                seq_id,
                projection_id,
                transcript_id,
            })
        }
        "projection-list" | "projections" => {
            let mut seq_id: Option<String> = None;
            let mut idx = 2usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--seq" => {
                        seq_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--seq",
                            "uniprot projection-list",
                        )?);
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for uniprot projection-list"
                        ));
                    }
                }
            }
            Ok(ShellCommand::UniprotProjectionList { seq_id })
        }
        "projection-show" | "show-projection" => {
            if tokens.len() != 3 {
                return Err("uniprot projection-show requires PROJECTION_ID".to_string());
            }
            let projection_id = tokens[2].trim().to_string();
            if projection_id.is_empty() {
                return Err("uniprot projection-show PROJECTION_ID must not be empty".to_string());
            }
            Ok(ShellCommand::UniprotProjectionShow { projection_id })
        }
        other => Err(format!(
            "Unknown uniprot subcommand '{other}' (expected fetch, import-swissprot, list, show, map, projection-list, projection-show)"
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

fn sequence_feature_roi_range_0based(
    dna: &DNAsequence,
    feature_id: usize,
) -> Result<(usize, usize), String> {
    let feature = dna
        .features()
        .get(feature_id)
        .ok_or_else(|| format!("Feature id {feature_id} is out of range"))?;
    let mut ranges: Vec<(usize, usize)> = Vec::new();
    collect_location_ranges_usize(&feature.location, &mut ranges);
    if ranges.is_empty() {
        let (from, to) = feature
            .location
            .find_bounds()
            .map_err(|e| format!("Could not read bounds for feature id {feature_id}: {e}"))?;
        if from < 0 || to < 0 {
            return Err(format!("Feature id {feature_id} has negative bounds"));
        }
        ranges.push((from as usize, to as usize));
    }
    let seq_len = dna.len();
    if seq_len == 0 {
        return Err("Template sequence is empty".to_string());
    }
    let start = ranges
        .iter()
        .map(|(start, _)| *start)
        .min()
        .ok_or_else(|| format!("Feature id {feature_id} has no usable ranges"))?;
    if start >= seq_len {
        return Err(format!(
            "Feature id {feature_id} starts outside sequence length {seq_len}"
        ));
    }
    let end_exclusive = ranges
        .iter()
        .map(|(_, end)| *end)
        .max()
        .ok_or_else(|| format!("Feature id {feature_id} has no usable ranges"))?
        .min(seq_len);
    if end_exclusive <= start {
        return Err(format!(
            "Feature id {feature_id} has invalid range {}..{} for sequence length {seq_len}",
            start, end_exclusive
        ));
    }
    Ok((start, end_exclusive))
}

fn build_seeded_primer_pair_operation(
    template: &str,
    roi_start_0based: usize,
    roi_end_0based_exclusive: usize,
) -> Operation {
    Operation::DesignPrimerPairs {
        template: template.to_string(),
        roi_start_0based,
        roi_end_0based: roi_end_0based_exclusive,
        forward: PrimerDesignSideConstraint::default(),
        reverse: PrimerDesignSideConstraint::default(),
        min_amplicon_bp: 120,
        max_amplicon_bp: 1200,
        pair_constraints: PrimerDesignPairConstraint::default(),
        max_tm_delta_c: Some(2.0),
        max_pairs: Some(200),
        report_id: None,
    }
}

fn build_seeded_qpcr_operation(
    template: &str,
    roi_start_0based: usize,
    roi_end_0based_exclusive: usize,
) -> Operation {
    Operation::DesignQpcrAssays {
        template: template.to_string(),
        roi_start_0based,
        roi_end_0based: roi_end_0based_exclusive,
        forward: PrimerDesignSideConstraint::default(),
        reverse: PrimerDesignSideConstraint::default(),
        probe: PrimerDesignSideConstraint::default(),
        min_amplicon_bp: 120,
        max_amplicon_bp: 1200,
        pair_constraints: PrimerDesignPairConstraint::default(),
        max_tm_delta_c: Some(2.0),
        max_probe_tm_delta_c: Some(10.0),
        max_assays: Some(200),
        report_id: None,
    }
}

fn parse_primers_command(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.len() < 2 {
        return Err(
            "primers requires a subcommand: design, design-qpcr, preflight, seed-from-feature, seed-from-splicing, list-reports, show-report, export-report, list-qpcr-reports, show-qpcr-report, export-qpcr-report"
                .to_string(),
        );
    }
    match tokens[1].as_str() {
        "design" => {
            if tokens.len() < 3 {
                return Err(
                    "primers design requires REQUEST_JSON_OR_@FILE [--backend auto|internal|primer3] [--primer3-exec PATH]"
                        .to_string(),
                );
            }
            let request_json = tokens[2].clone();
            let mut backend: Option<PrimerDesignBackend> = None;
            let mut primer3_executable: Option<String> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--backend" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--backend", "primers design")?;
                        backend = Some(parse_primer_design_backend(&raw)?);
                    }
                    "--primer3-exec" | "--primer3-executable" => {
                        let flag = tokens[idx].clone();
                        primer3_executable = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            &flag,
                            "primers design",
                        )?);
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for primers design"));
                    }
                }
            }
            Ok(ShellCommand::PrimersDesign {
                request_json,
                backend,
                primer3_executable,
            })
        }
        "design-qpcr" => {
            if tokens.len() < 3 {
                return Err(
                    "primers design-qpcr requires REQUEST_JSON_OR_@FILE [--backend auto|internal|primer3] [--primer3-exec PATH]"
                        .to_string(),
                );
            }
            let request_json = tokens[2].clone();
            let mut backend: Option<PrimerDesignBackend> = None;
            let mut primer3_executable: Option<String> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--backend" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--backend",
                            "primers design-qpcr",
                        )?;
                        backend = Some(parse_primer_design_backend(&raw)?);
                    }
                    "--primer3-exec" | "--primer3-executable" => {
                        let flag = tokens[idx].clone();
                        primer3_executable = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            &flag,
                            "primers design-qpcr",
                        )?);
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for primers design-qpcr"));
                    }
                }
            }
            Ok(ShellCommand::PrimersDesignQpcr {
                request_json,
                backend,
                primer3_executable,
            })
        }
        "preflight" => {
            let mut backend: Option<PrimerDesignBackend> = None;
            let mut primer3_executable: Option<String> = None;
            let mut idx = 2usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--backend" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--backend", "primers preflight")?;
                        backend = Some(parse_primer_design_backend(&raw)?);
                    }
                    "--primer3-exec" | "--primer3-executable" => {
                        let flag = tokens[idx].clone();
                        primer3_executable = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            &flag,
                            "primers preflight",
                        )?);
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for primers preflight"));
                    }
                }
            }
            Ok(ShellCommand::PrimersPreflight {
                backend,
                primer3_executable,
            })
        }
        "seed-from-feature" => {
            if tokens.len() != 4 {
                return Err("primers seed-from-feature requires SEQ_ID FEATURE_ID".to_string());
            }
            let seq_id = tokens[2].clone();
            let feature_id = tokens[3].parse::<usize>().map_err(|e| {
                format!(
                    "Invalid feature id '{}' for primers seed-from-feature: {e}",
                    tokens[3]
                )
            })?;
            Ok(ShellCommand::PrimersSeedFromFeature { seq_id, feature_id })
        }
        "seed-from-splicing" => {
            if tokens.len() != 4 {
                return Err("primers seed-from-splicing requires SEQ_ID FEATURE_ID".to_string());
            }
            let seq_id = tokens[2].clone();
            let feature_id = tokens[3].parse::<usize>().map_err(|e| {
                format!(
                    "Invalid feature id '{}' for primers seed-from-splicing: {e}",
                    tokens[3]
                )
            })?;
            Ok(ShellCommand::PrimersSeedFromSplicing { seq_id, feature_id })
        }
        "list-reports" => {
            if tokens.len() != 2 {
                return Err("primers list-reports takes no options".to_string());
            }
            Ok(ShellCommand::PrimersListReports)
        }
        "show-report" => {
            if tokens.len() != 3 {
                return Err("primers show-report requires REPORT_ID".to_string());
            }
            Ok(ShellCommand::PrimersShowReport {
                report_id: tokens[2].clone(),
            })
        }
        "export-report" => {
            if tokens.len() != 4 {
                return Err("primers export-report requires REPORT_ID OUTPUT.json".to_string());
            }
            Ok(ShellCommand::PrimersExportReport {
                report_id: tokens[2].clone(),
                path: tokens[3].clone(),
            })
        }
        "list-qpcr-reports" => {
            if tokens.len() != 2 {
                return Err("primers list-qpcr-reports takes no options".to_string());
            }
            Ok(ShellCommand::PrimersListQpcrReports)
        }
        "show-qpcr-report" => {
            if tokens.len() != 3 {
                return Err("primers show-qpcr-report requires REPORT_ID".to_string());
            }
            Ok(ShellCommand::PrimersShowQpcrReport {
                report_id: tokens[2].clone(),
            })
        }
        "export-qpcr-report" => {
            if tokens.len() != 4 {
                return Err("primers export-qpcr-report requires REPORT_ID OUTPUT.json".to_string());
            }
            Ok(ShellCommand::PrimersExportQpcrReport {
                report_id: tokens[2].clone(),
                path: tokens[3].clone(),
            })
        }
        other => Err(format!(
            "Unknown primers subcommand '{other}' (expected design, design-qpcr, preflight, seed-from-feature, seed-from-splicing, list-reports, show-report, export-report, list-qpcr-reports, show-qpcr-report, export-qpcr-report)"
        )),
    }
}

fn parse_dotplot_mode(raw: &str) -> Result<DotplotMode, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "self_forward" | "self-forward" | "forward" | "self" => Ok(DotplotMode::SelfForward),
        "self_reverse_complement"
        | "self-reverse-complement"
        | "self_revcomp"
        | "self-revcomp"
        | "revcomp"
        | "reverse_complement"
        | "reverse-complement" => Ok(DotplotMode::SelfReverseComplement),
        other => Err(format!(
            "Unsupported dotplot mode '{other}', expected self_forward or self_reverse_complement"
        )),
    }
}

fn parse_flexibility_model(raw: &str) -> Result<FlexibilityModel, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "at_richness" | "at-richness" | "at" | "at_content" | "at-content" => {
            Ok(FlexibilityModel::AtRichness)
        }
        "at_skew" | "at-skew" | "skew" => Ok(FlexibilityModel::AtSkew),
        other => Err(format!(
            "Unsupported flexibility model '{other}', expected at_richness or at_skew"
        )),
    }
}

fn parse_rna_read_profile(raw: &str) -> Result<RnaReadInterpretationProfile, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "nanopore_cdna_v1" | "nanopore" | "nanopore_cdna" => {
            Ok(RnaReadInterpretationProfile::NanoporeCdnaV1)
        }
        "short_read_v1" | "shortread" | "short_read" => {
            Ok(RnaReadInterpretationProfile::ShortReadV1)
        }
        "transposon_v1" | "transposon" => Ok(RnaReadInterpretationProfile::TransposonV1),
        other => Err(format!(
            "Unsupported RNA-read profile '{other}', expected nanopore_cdna_v1|short_read_v1|transposon_v1"
        )),
    }
}

fn parse_rna_read_input_format(raw: &str) -> Result<RnaReadInputFormat, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "fasta" | "fa" => Ok(RnaReadInputFormat::Fasta),
        other => Err(format!(
            "Unsupported RNA-read input format '{other}', expected fasta"
        )),
    }
}

fn parse_splicing_scope_preset(raw: &str) -> Result<SplicingScopePreset, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "all_overlapping_both_strands" | "all" | "broad" => {
            Ok(SplicingScopePreset::AllOverlappingBothStrands)
        }
        "target_group_any_strand" | "target_group" | "group" => {
            Ok(SplicingScopePreset::TargetGroupAnyStrand)
        }
        "all_overlapping_target_strand" | "target_strand" | "strand" => {
            Ok(SplicingScopePreset::AllOverlappingTargetStrand)
        }
        "target_group_target_strand" | "group_strand" | "legacy" => {
            Ok(SplicingScopePreset::TargetGroupTargetStrand)
        }
        other => Err(format!(
            "Unsupported scope preset '{other}', expected all_overlapping_both_strands|target_group_any_strand|all_overlapping_target_strand|target_group_target_strand"
        )),
    }
}

fn parse_rna_read_origin_mode(raw: &str) -> Result<RnaReadOriginMode, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "single_gene" | "single-gene" | "single" | "legacy" => {
            Ok(RnaReadOriginMode::SingleGene)
        }
        "multi_gene_sparse" | "multi-gene-sparse" | "multi_sparse" | "multi" => {
            Ok(RnaReadOriginMode::MultiGeneSparse)
        }
        other => Err(format!(
            "Unsupported origin mode '{other}', expected single_gene|multi_gene_sparse"
        )),
    }
}

fn parse_rna_read_report_mode(raw: &str) -> Result<RnaReadReportMode, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "full" => Ok(RnaReadReportMode::Full),
        "seed_passed_only" | "seed-passed-only" | "seed_passed" | "seed-only" => {
            Ok(RnaReadReportMode::SeedPassedOnly)
        }
        other => Err(format!(
            "Unsupported report mode '{other}', expected full|seed_passed_only"
        )),
    }
}

fn parse_rna_read_hit_selection(raw: &str) -> Result<RnaReadHitSelection, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "all" => Ok(RnaReadHitSelection::All),
        "seed_passed" | "seed" => Ok(RnaReadHitSelection::SeedPassed),
        "aligned" => Ok(RnaReadHitSelection::Aligned),
        other => Err(format!(
            "Unsupported hit selection '{other}', expected all|seed_passed|aligned"
        )),
    }
}

fn parse_rna_read_score_density_scale(raw: &str) -> Result<RnaReadScoreDensityScale, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "linear" | "lin" => Ok(RnaReadScoreDensityScale::Linear),
        "log" | "log1p" => Ok(RnaReadScoreDensityScale::Log),
        other => Err(format!(
            "Unsupported score-density scale '{other}', expected linear|log"
        )),
    }
}

fn parse_dotplot_command(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.len() < 2 {
        return Err("dotplot requires a subcommand: compute, list, show".to_string());
    }
    match tokens[1].as_str() {
        "compute" => {
            if tokens.len() < 3 {
                return Err(
                    "dotplot compute requires SEQ_ID [--start N] [--end N] [--mode MODE] [--word-size N] [--step N] [--max-mismatches N] [--tile-bp N] [--id DOTPLOT_ID]"
                        .to_string(),
                );
            }
            let seq_id = tokens[2].trim().to_string();
            if seq_id.is_empty() {
                return Err("dotplot compute SEQ_ID must not be empty".to_string());
            }
            let mut span_start_0based: Option<usize> = None;
            let mut span_end_0based: Option<usize> = None;
            let mut mode = DotplotMode::SelfForward;
            let mut word_size = 12usize;
            let mut step_bp = 2usize;
            let mut max_mismatches = 0usize;
            let mut tile_bp: Option<usize> = None;
            let mut dotplot_id: Option<String> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--start" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--start", "dotplot compute")?;
                        span_start_0based = Some(raw.parse::<usize>().map_err(|e| {
                            format!("Invalid --start value '{raw}' for dotplot compute: {e}")
                        })?);
                    }
                    "--end" => {
                        let raw = parse_option_path(tokens, &mut idx, "--end", "dotplot compute")?;
                        span_end_0based = Some(raw.parse::<usize>().map_err(|e| {
                            format!("Invalid --end value '{raw}' for dotplot compute: {e}")
                        })?);
                    }
                    "--mode" => {
                        let raw = parse_option_path(tokens, &mut idx, "--mode", "dotplot compute")?;
                        mode = parse_dotplot_mode(&raw)?;
                    }
                    "--word-size" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--word-size", "dotplot compute")?;
                        word_size = raw.parse::<usize>().map_err(|e| {
                            format!("Invalid --word-size value '{raw}' for dotplot compute: {e}")
                        })?;
                    }
                    "--step" | "--step-bp" => {
                        let flag = tokens[idx].clone();
                        let raw = parse_option_path(tokens, &mut idx, &flag, "dotplot compute")?;
                        step_bp = raw.parse::<usize>().map_err(|e| {
                            format!("Invalid {flag} value '{raw}' for dotplot compute: {e}")
                        })?;
                    }
                    "--max-mismatches" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--max-mismatches",
                            "dotplot compute",
                        )?;
                        max_mismatches = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --max-mismatches value '{raw}' for dotplot compute: {e}"
                            )
                        })?;
                    }
                    "--tile-bp" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--tile-bp", "dotplot compute")?;
                        tile_bp = Some(raw.parse::<usize>().map_err(|e| {
                            format!("Invalid --tile-bp value '{raw}' for dotplot compute: {e}")
                        })?);
                    }
                    "--id" => {
                        let raw = parse_option_path(tokens, &mut idx, "--id", "dotplot compute")?;
                        dotplot_id = Some(raw);
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for dotplot compute"));
                    }
                }
            }
            Ok(ShellCommand::DotplotCompute {
                seq_id,
                span_start_0based,
                span_end_0based,
                mode,
                word_size,
                step_bp,
                max_mismatches,
                tile_bp,
                dotplot_id,
            })
        }
        "list" => {
            if tokens.len() > 3 {
                return Err("dotplot list expects at most one optional SEQ_ID".to_string());
            }
            let seq_id = if tokens.len() == 3 {
                let value = tokens[2].trim();
                if value.is_empty() {
                    None
                } else {
                    Some(value.to_string())
                }
            } else {
                None
            };
            Ok(ShellCommand::DotplotList { seq_id })
        }
        "show" => {
            if tokens.len() != 3 {
                return Err("dotplot show requires DOTPLOT_ID".to_string());
            }
            Ok(ShellCommand::DotplotShow {
                dotplot_id: tokens[2].clone(),
            })
        }
        other => Err(format!(
            "Unknown dotplot subcommand '{other}' (expected compute, list, show)"
        )),
    }
}

fn parse_flex_command(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.len() < 2 {
        return Err("flex requires a subcommand: compute, list, show".to_string());
    }
    match tokens[1].as_str() {
        "compute" => {
            if tokens.len() < 3 {
                return Err(
                    "flex compute requires SEQ_ID [--start N] [--end N] [--model MODEL] [--bin-bp N] [--smoothing-bp N] [--id TRACK_ID]"
                        .to_string(),
                );
            }
            let seq_id = tokens[2].trim().to_string();
            if seq_id.is_empty() {
                return Err("flex compute SEQ_ID must not be empty".to_string());
            }
            let mut span_start_0based: Option<usize> = None;
            let mut span_end_0based: Option<usize> = None;
            let mut model = FlexibilityModel::AtRichness;
            let mut bin_bp = 25usize;
            let mut smoothing_bp: Option<usize> = None;
            let mut track_id: Option<String> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--start" => {
                        let raw = parse_option_path(tokens, &mut idx, "--start", "flex compute")?;
                        span_start_0based = Some(raw.parse::<usize>().map_err(|e| {
                            format!("Invalid --start value '{raw}' for flex compute: {e}")
                        })?);
                    }
                    "--end" => {
                        let raw = parse_option_path(tokens, &mut idx, "--end", "flex compute")?;
                        span_end_0based = Some(raw.parse::<usize>().map_err(|e| {
                            format!("Invalid --end value '{raw}' for flex compute: {e}")
                        })?);
                    }
                    "--model" => {
                        let raw = parse_option_path(tokens, &mut idx, "--model", "flex compute")?;
                        model = parse_flexibility_model(&raw)?;
                    }
                    "--bin-bp" => {
                        let raw = parse_option_path(tokens, &mut idx, "--bin-bp", "flex compute")?;
                        bin_bp = raw.parse::<usize>().map_err(|e| {
                            format!("Invalid --bin-bp value '{raw}' for flex compute: {e}")
                        })?;
                    }
                    "--smoothing-bp" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--smoothing-bp", "flex compute")?;
                        smoothing_bp = Some(raw.parse::<usize>().map_err(|e| {
                            format!("Invalid --smoothing-bp value '{raw}' for flex compute: {e}")
                        })?);
                    }
                    "--id" => {
                        let raw = parse_option_path(tokens, &mut idx, "--id", "flex compute")?;
                        track_id = Some(raw);
                    }
                    other => return Err(format!("Unknown option '{other}' for flex compute")),
                }
            }
            Ok(ShellCommand::FlexCompute {
                seq_id,
                span_start_0based,
                span_end_0based,
                model,
                bin_bp,
                smoothing_bp,
                track_id,
            })
        }
        "list" => {
            if tokens.len() > 3 {
                return Err("flex list expects at most one optional SEQ_ID".to_string());
            }
            let seq_id = if tokens.len() == 3 {
                let value = tokens[2].trim();
                if value.is_empty() {
                    None
                } else {
                    Some(value.to_string())
                }
            } else {
                None
            };
            Ok(ShellCommand::FlexList { seq_id })
        }
        "show" => {
            if tokens.len() != 3 {
                return Err("flex show requires TRACK_ID".to_string());
            }
            Ok(ShellCommand::FlexShow {
                track_id: tokens[2].clone(),
            })
        }
        other => Err(format!(
            "Unknown flex subcommand '{other}' (expected compute, list, show)"
        )),
    }
}

fn parse_rna_reads_command(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.len() < 2 {
        return Err(
            "rna-reads requires a subcommand: interpret, list-reports, show-report, export-report, export-hits-fasta, export-sample-sheet, export-paths-tsv, export-abundance-tsv, export-score-density-svg"
                .to_string(),
        );
    }
    match tokens[1].as_str() {
        "interpret" => {
            if tokens.len() < 5 {
                return Err(
                    "rna-reads interpret requires SEQ_ID FEATURE_ID INPUT.fa[.gz] [--report-id ID] [--report-mode full|seed_passed_only] [--checkpoint-path PATH] [--checkpoint-every-reads N] [--resume-from-checkpoint|--no-resume-from-checkpoint] [--profile PROFILE] [--format fasta] [--scope SCOPE] [--origin-mode single_gene|multi_gene_sparse] [--target-gene GENE_ID]... [--roi-seed-capture|--no-roi-seed-capture] [--kmer-len N] [--short-max-bp N] [--long-window-bp N] [--long-window-count N] [--min-seed-hit-fraction F] [--min-weighted-seed-hit-fraction F] [--min-unique-matched-kmers N] [--min-chain-consistency-fraction F] [--max-median-transcript-gap F] [--min-confirmed-transitions N] [--min-transition-support-fraction F] [--cdna-poly-t-flip|--no-cdna-poly-t-flip] [--poly-t-prefix-min-bp N] [--align-band-bp N] [--align-min-identity F] [--max-secondary-mappings N]"
                        .to_string(),
                );
            }
            let seq_id = tokens[2].trim().to_string();
            if seq_id.is_empty() {
                return Err("rna-reads interpret SEQ_ID must not be empty".to_string());
            }
            let seed_feature_id = tokens[3].parse::<usize>().map_err(|e| {
                format!(
                    "Invalid FEATURE_ID '{}' for rna-reads interpret: {e}",
                    tokens[3]
                )
            })?;
            let input_path = tokens[4].trim().to_string();
            if input_path.is_empty() {
                return Err("rna-reads interpret INPUT.fa[.gz] must not be empty".to_string());
            }
            let mut profile = RnaReadInterpretationProfile::NanoporeCdnaV1;
            let mut input_format = RnaReadInputFormat::Fasta;
            let mut scope = SplicingScopePreset::AllOverlappingBothStrands;
            let mut origin_mode = RnaReadOriginMode::SingleGene;
            let mut target_gene_ids: Vec<String> = vec![];
            let mut roi_seed_capture_enabled = false;
            let mut seed_filter = RnaReadSeedFilterConfig::default();
            let mut align_config = RnaReadAlignConfig::default();
            let mut report_id: Option<String> = None;
            let mut report_mode = RnaReadReportMode::Full;
            let mut checkpoint_path: Option<String> = None;
            let mut checkpoint_every_reads = 10_000usize;
            let mut resume_from_checkpoint = false;
            let mut idx = 5usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--report-id" => {
                        report_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--report-id",
                            "rna-reads interpret",
                        )?);
                    }
                    "--profile" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--profile",
                            "rna-reads interpret",
                        )?;
                        profile = parse_rna_read_profile(&raw)?;
                    }
                    "--report-mode" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--report-mode",
                            "rna-reads interpret",
                        )?;
                        report_mode = parse_rna_read_report_mode(&raw)?;
                    }
                    "--checkpoint-path" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--checkpoint-path",
                            "rna-reads interpret",
                        )?;
                        checkpoint_path = Some(raw);
                    }
                    "--checkpoint-every-reads" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--checkpoint-every-reads",
                            "rna-reads interpret",
                        )?;
                        checkpoint_every_reads = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --checkpoint-every-reads value '{raw}' for rna-reads interpret: {e}"
                            )
                        })?;
                    }
                    "--resume-from-checkpoint" => {
                        resume_from_checkpoint = true;
                        idx += 1;
                    }
                    "--no-resume-from-checkpoint" => {
                        resume_from_checkpoint = false;
                        idx += 1;
                    }
                    "--format" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--format", "rna-reads interpret")?;
                        input_format = parse_rna_read_input_format(&raw)?;
                    }
                    "--scope" => {
                        let raw =
                            parse_option_path(tokens, &mut idx, "--scope", "rna-reads interpret")?;
                        scope = parse_splicing_scope_preset(&raw)?;
                    }
                    "--origin-mode" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--origin-mode",
                            "rna-reads interpret",
                        )?;
                        origin_mode = parse_rna_read_origin_mode(&raw)?;
                    }
                    "--target-gene" | "--target-gene-id" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--target-gene",
                            "rna-reads interpret",
                        )?;
                        let gene_id = raw.trim();
                        if gene_id.is_empty() {
                            return Err(
                                "--target-gene requires a non-empty gene identifier".to_string()
                            );
                        }
                        target_gene_ids.push(gene_id.to_string());
                    }
                    "--roi-seed-capture" => {
                        roi_seed_capture_enabled = true;
                        idx += 1;
                    }
                    "--no-roi-seed-capture" => {
                        roi_seed_capture_enabled = false;
                        idx += 1;
                    }
                    "--kmer-len" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--kmer-len",
                            "rna-reads interpret",
                        )?;
                        seed_filter.kmer_len = raw.parse::<usize>().map_err(|e| {
                            format!("Invalid --kmer-len value '{raw}' for rna-reads interpret: {e}")
                        })?;
                    }
                    "--short-max-bp" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--short-max-bp",
                            "rna-reads interpret",
                        )?;
                        seed_filter.short_full_hash_max_bp = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --short-max-bp value '{raw}' for rna-reads interpret: {e}"
                            )
                        })?;
                    }
                    "--long-window-bp" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--long-window-bp",
                            "rna-reads interpret",
                        )?;
                        seed_filter.long_window_bp = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --long-window-bp value '{raw}' for rna-reads interpret: {e}"
                            )
                        })?;
                    }
                    "--long-window-count" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--long-window-count",
                            "rna-reads interpret",
                        )?;
                        seed_filter.long_window_count = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --long-window-count value '{raw}' for rna-reads interpret: {e}"
                            )
                        })?;
                    }
                    "--min-seed-hit-fraction" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--min-seed-hit-fraction",
                            "rna-reads interpret",
                        )?;
                        seed_filter.min_seed_hit_fraction =
                            raw.parse::<f64>().map_err(|e| {
                                format!(
                                    "Invalid --min-seed-hit-fraction value '{raw}' for rna-reads interpret: {e}"
                                )
                            })?;
                    }
                    "--min-weighted-seed-hit-fraction" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--min-weighted-seed-hit-fraction",
                            "rna-reads interpret",
                        )?;
                        seed_filter.min_weighted_seed_hit_fraction =
                            raw.parse::<f64>().map_err(|e| {
                                format!(
                                    "Invalid --min-weighted-seed-hit-fraction value '{raw}' for rna-reads interpret: {e}"
                                )
                            })?;
                    }
                    "--min-unique-matched-kmers" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--min-unique-matched-kmers",
                            "rna-reads interpret",
                        )?;
                        seed_filter.min_unique_matched_kmers =
                            raw.parse::<usize>().map_err(|e| {
                                format!(
                                    "Invalid --min-unique-matched-kmers value '{raw}' for rna-reads interpret: {e}"
                                )
                            })?;
                    }
                    "--max-median-transcript-gap" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--max-median-transcript-gap",
                            "rna-reads interpret",
                        )?;
                        seed_filter.max_median_transcript_gap =
                            raw.parse::<f64>().map_err(|e| {
                                format!(
                                    "Invalid --max-median-transcript-gap value '{raw}' for rna-reads interpret: {e}"
                                )
                            })?;
                    }
                    "--min-chain-consistency-fraction" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--min-chain-consistency-fraction",
                            "rna-reads interpret",
                        )?;
                        seed_filter.min_chain_consistency_fraction =
                            raw.parse::<f64>().map_err(|e| {
                                format!(
                                    "Invalid --min-chain-consistency-fraction value '{raw}' for rna-reads interpret: {e}"
                                )
                            })?;
                    }
                    "--min-confirmed-transitions" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--min-confirmed-transitions",
                            "rna-reads interpret",
                        )?;
                        seed_filter.min_confirmed_exon_transitions =
                            raw.parse::<usize>().map_err(|e| {
                                format!(
                                    "Invalid --min-confirmed-transitions value '{raw}' for rna-reads interpret: {e}"
                                )
                            })?;
                    }
                    "--min-transition-support-fraction" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--min-transition-support-fraction",
                            "rna-reads interpret",
                        )?;
                        seed_filter.min_transition_support_fraction =
                            raw.parse::<f64>().map_err(|e| {
                                format!(
                                    "Invalid --min-transition-support-fraction value '{raw}' for rna-reads interpret: {e}"
                                )
                            })?;
                    }
                    "--cdna-poly-t-flip" => {
                        seed_filter.cdna_poly_t_flip_enabled = true;
                        idx += 1;
                    }
                    "--no-cdna-poly-t-flip" => {
                        seed_filter.cdna_poly_t_flip_enabled = false;
                        idx += 1;
                    }
                    "--poly-t-prefix-min-bp" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--poly-t-prefix-min-bp",
                            "rna-reads interpret",
                        )?;
                        seed_filter.poly_t_prefix_min_bp = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --poly-t-prefix-min-bp value '{raw}' for rna-reads interpret: {e}"
                            )
                        })?;
                    }
                    "--align-band-bp" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--align-band-bp",
                            "rna-reads interpret",
                        )?;
                        align_config.band_width_bp = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --align-band-bp value '{raw}' for rna-reads interpret: {e}"
                            )
                        })?;
                    }
                    "--align-min-identity" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--align-min-identity",
                            "rna-reads interpret",
                        )?;
                        align_config.min_identity_fraction = raw.parse::<f64>().map_err(|e| {
                            format!(
                                "Invalid --align-min-identity value '{raw}' for rna-reads interpret: {e}"
                            )
                        })?;
                    }
                    "--max-secondary-mappings" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--max-secondary-mappings",
                            "rna-reads interpret",
                        )?;
                        align_config.max_secondary_mappings =
                            raw.parse::<usize>().map_err(|e| {
                                format!(
                                    "Invalid --max-secondary-mappings value '{raw}' for rna-reads interpret: {e}"
                                )
                            })?;
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for rna-reads interpret"));
                    }
                }
            }
            Ok(ShellCommand::RnaReadsInterpret {
                seq_id,
                seed_feature_id,
                input_path,
                profile,
                input_format,
                scope,
                origin_mode,
                target_gene_ids,
                roi_seed_capture_enabled,
                seed_filter,
                align_config,
                report_id,
                report_mode,
                checkpoint_path,
                checkpoint_every_reads,
                resume_from_checkpoint,
            })
        }
        "list-reports" => {
            if tokens.len() > 3 {
                return Err(
                    "rna-reads list-reports expects at most one optional SEQ_ID".to_string()
                );
            }
            let seq_id = if tokens.len() == 3 {
                let value = tokens[2].trim();
                if value.is_empty() {
                    None
                } else {
                    Some(value.to_string())
                }
            } else {
                None
            };
            Ok(ShellCommand::RnaReadsListReports { seq_id })
        }
        "show-report" => {
            if tokens.len() != 3 {
                return Err("rna-reads show-report requires REPORT_ID".to_string());
            }
            Ok(ShellCommand::RnaReadsShowReport {
                report_id: tokens[2].clone(),
            })
        }
        "export-report" => {
            if tokens.len() != 4 {
                return Err("rna-reads export-report requires REPORT_ID OUTPUT.json".to_string());
            }
            Ok(ShellCommand::RnaReadsExportReport {
                report_id: tokens[2].clone(),
                path: tokens[3].clone(),
            })
        }
        "export-hits-fasta" => {
            if tokens.len() < 4 {
                return Err(
                    "rna-reads export-hits-fasta requires REPORT_ID OUTPUT.fa [--selection all|seed_passed|aligned]"
                        .to_string(),
                );
            }
            let report_id = tokens[2].clone();
            let path = tokens[3].clone();
            let mut selection = RnaReadHitSelection::Aligned;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--selection" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--selection",
                            "rna-reads export-hits-fasta",
                        )?;
                        selection = parse_rna_read_hit_selection(&raw)?;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for rna-reads export-hits-fasta"
                        ));
                    }
                }
            }
            Ok(ShellCommand::RnaReadsExportHitsFasta {
                report_id,
                path,
                selection,
            })
        }
        "export-sample-sheet" => {
            if tokens.len() < 3 {
                return Err(
                    "rna-reads export-sample-sheet requires OUTPUT.tsv [--seq-id ID] [--report-id ID]... [--append]"
                        .to_string(),
                );
            }
            let path = tokens[2].clone();
            let mut seq_id: Option<String> = None;
            let mut report_ids: Vec<String> = vec![];
            let mut append = false;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--seq-id" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--seq-id",
                            "rna-reads export-sample-sheet",
                        )?;
                        seq_id = Some(raw);
                    }
                    "--report-id" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--report-id",
                            "rna-reads export-sample-sheet",
                        )?;
                        report_ids.push(raw);
                    }
                    "--append" => {
                        append = true;
                        idx += 1;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for rna-reads export-sample-sheet"
                        ));
                    }
                }
            }
            Ok(ShellCommand::RnaReadsExportSampleSheet {
                path,
                seq_id,
                report_ids,
                append,
            })
        }
        "export-paths-tsv" => {
            if tokens.len() < 4 {
                return Err(
                    "rna-reads export-paths-tsv requires REPORT_ID OUTPUT.tsv [--selection all|seed_passed|aligned]"
                        .to_string(),
                );
            }
            let report_id = tokens[2].clone();
            let path = tokens[3].clone();
            let mut selection = RnaReadHitSelection::All;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--selection" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--selection",
                            "rna-reads export-paths-tsv",
                        )?;
                        selection = parse_rna_read_hit_selection(&raw)?;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for rna-reads export-paths-tsv"
                        ));
                    }
                }
            }
            Ok(ShellCommand::RnaReadsExportExonPathsTsv {
                report_id,
                path,
                selection,
            })
        }
        "export-abundance-tsv" => {
            if tokens.len() < 4 {
                return Err(
                    "rna-reads export-abundance-tsv requires REPORT_ID OUTPUT.tsv [--selection all|seed_passed|aligned]"
                        .to_string(),
                );
            }
            let report_id = tokens[2].clone();
            let path = tokens[3].clone();
            let mut selection = RnaReadHitSelection::All;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--selection" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--selection",
                            "rna-reads export-abundance-tsv",
                        )?;
                        selection = parse_rna_read_hit_selection(&raw)?;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for rna-reads export-abundance-tsv"
                        ));
                    }
                }
            }
            Ok(ShellCommand::RnaReadsExportExonAbundanceTsv {
                report_id,
                path,
                selection,
            })
        }
        "export-score-density-svg" => {
            if tokens.len() < 4 {
                return Err(
                    "rna-reads export-score-density-svg requires REPORT_ID OUTPUT.svg [--scale linear|log]"
                        .to_string(),
                );
            }
            let report_id = tokens[2].clone();
            let path = tokens[3].clone();
            let mut scale = RnaReadScoreDensityScale::Log;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--scale" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--scale",
                            "rna-reads export-score-density-svg",
                        )?;
                        scale = parse_rna_read_score_density_scale(&raw)?;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for rna-reads export-score-density-svg"
                        ));
                    }
                }
            }
            Ok(ShellCommand::RnaReadsExportScoreDensitySvg {
                report_id,
                path,
                scale,
            })
        }
        other => Err(format!(
            "Unknown rna-reads subcommand '{other}' (expected interpret, list-reports, show-report, export-report, export-hits-fasta, export-sample-sheet, export-paths-tsv, export-abundance-tsv, export-score-density-svg)"
        )),
    }
}

fn parse_macros_command(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.len() < 2 {
        return Err(
            "macros requires a subcommand: run, instance-list, instance-show, template-list, template-show, template-put, template-delete, template-import, template-run"
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
        "instance-list" => {
            if tokens.len() != 2 {
                return Err("macros instance-list takes no options".to_string());
            }
            Ok(ShellCommand::MacrosInstanceList)
        }
        "instance-show" => {
            if tokens.len() != 3 {
                return Err("macros instance-show requires MACRO_INSTANCE_ID".to_string());
            }
            Ok(ShellCommand::MacrosInstanceShow {
                macro_instance_id: tokens[2].clone(),
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
                    "macros template-put requires TEMPLATE_NAME (--script SCRIPT_OR_@FILE | --file PATH) [--description TEXT] [--details-url URL] [--param NAME|NAME=DEFAULT ...] [--input-port PORT_ID:KIND[:one|many][:required|optional][:description]] [--output-port ...]"
                        .to_string(),
                );
            }
            let name = tokens[2].clone();
            let mut description: Option<String> = None;
            let mut details_url: Option<String> = None;
            let mut parameters: Vec<WorkflowMacroTemplateParam> = vec![];
            let mut input_ports: Vec<WorkflowMacroTemplatePort> = vec![];
            let mut output_ports: Vec<WorkflowMacroTemplatePort> = vec![];
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
                    "--input-port" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--input-port",
                            "macros template-put",
                        )?;
                        input_ports.push(parse_workflow_template_port_spec(&raw)?);
                    }
                    "--output-port" => {
                        let raw = parse_option_path(
                            tokens,
                            &mut idx,
                            "--output-port",
                            "macros template-put",
                        )?;
                        output_ports.push(parse_workflow_template_port_spec(&raw)?);
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
                input_ports,
                output_ports,
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
                    "macros template-run requires TEMPLATE_NAME [--bind KEY=VALUE ...] [--transactional] [--validate-only]"
                        .to_string(),
                );
            }
            let name = tokens[2].clone();
            let mut bindings: HashMap<String, String> = HashMap::new();
            let mut transactional = false;
            let mut validate_only = false;
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
                    "--validate-only" | "--dry-run" => {
                        validate_only = true;
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
                validate_only,
            })
        }
        other => Err(format!(
            "Unknown macros subcommand '{other}' (expected run, instance-list, instance-show, template-list, template-show, template-put, template-delete, template-import, template-run)"
        )),
    }
}

fn parse_routines_command(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.len() < 2 {
        return Err("routines requires a subcommand: list, explain, compare".to_string());
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
        "explain" => {
            if tokens.len() < 3 {
                return Err("routines explain requires ROUTINE_ID [--catalog PATH]".to_string());
            }
            let routine_id = tokens[2].trim().to_string();
            if routine_id.is_empty() {
                return Err("routines explain ROUTINE_ID cannot be empty".to_string());
            }
            let mut catalog_path: Option<String> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--catalog" => {
                        catalog_path = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--catalog",
                            "routines explain",
                        )?);
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for routines explain"));
                    }
                }
            }
            Ok(ShellCommand::RoutinesExplain {
                catalog_path,
                routine_id,
            })
        }
        "compare" => {
            if tokens.len() < 4 {
                return Err(
                    "routines compare requires ROUTINE_A ROUTINE_B [--catalog PATH]".to_string(),
                );
            }
            let left_routine_id = tokens[2].trim().to_string();
            let right_routine_id = tokens[3].trim().to_string();
            if left_routine_id.is_empty() || right_routine_id.is_empty() {
                return Err("routines compare routine ids cannot be empty".to_string());
            }
            let mut catalog_path: Option<String> = None;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--catalog" => {
                        catalog_path = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--catalog",
                            "routines compare",
                        )?);
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for routines compare"));
                    }
                }
            }
            Ok(ShellCommand::RoutinesCompare {
                catalog_path,
                left_routine_id,
                right_routine_id,
            })
        }
        other => Err(format!(
            "Unknown routines subcommand '{other}' (expected list, explain, compare)"
        )),
    }
}

fn parse_planning_profile_scope(raw: &str) -> Result<PlanningProfileScope, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "global" => Ok(PlanningProfileScope::Global),
        "project" | "project_override" | "project-override" => {
            Ok(PlanningProfileScope::ProjectOverride)
        }
        "agent" | "agent_overlay" | "confirmed_agent_overlay" | "confirmed-agent-overlay" => {
            Ok(PlanningProfileScope::ConfirmedAgentOverlay)
        }
        "effective" => Ok(PlanningProfileScope::Effective),
        other => Err(format!(
            "Unsupported planning profile scope '{other}' (expected global|project_override|confirmed_agent_overlay|effective)"
        )),
    }
}

fn parse_planning_suggestion_status(raw: &str) -> Result<PlanningSuggestionStatus, String> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "pending" => Ok(PlanningSuggestionStatus::Pending),
        "accepted" => Ok(PlanningSuggestionStatus::Accepted),
        "rejected" => Ok(PlanningSuggestionStatus::Rejected),
        other => Err(format!(
            "Unsupported planning suggestion status '{other}' (expected pending|accepted|rejected)"
        )),
    }
}

fn parse_planning_command(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.len() < 2 {
        return Err(
            "planning requires a subcommand: profile, objective, suggestions, sync".to_string(),
        );
    }
    match tokens[1].as_str() {
        "profile" => {
            if tokens.len() < 3 {
                return Err("planning profile requires a subcommand: show, set, clear".to_string());
            }
            match tokens[2].as_str() {
                "show" => {
                    let mut scope = PlanningProfileScope::Effective;
                    let mut idx = 3usize;
                    while idx < tokens.len() {
                        match tokens[idx].as_str() {
                            "--scope" => {
                                let raw = parse_option_path(
                                    tokens,
                                    &mut idx,
                                    "--scope",
                                    "planning profile show",
                                )?;
                                scope = parse_planning_profile_scope(&raw)?;
                            }
                            other => {
                                return Err(format!(
                                    "Unknown option '{other}' for planning profile show"
                                ));
                            }
                        }
                    }
                    Ok(ShellCommand::PlanningProfileShow { scope })
                }
                "set" => {
                    if tokens.len() < 4 {
                        return Err(
                            "planning profile set requires JSON_OR_@FILE [--scope SCOPE]"
                                .to_string(),
                        );
                    }
                    let mut scope = PlanningProfileScope::ProjectOverride;
                    let payload_json = tokens[3].clone();
                    let mut idx = 4usize;
                    while idx < tokens.len() {
                        match tokens[idx].as_str() {
                            "--scope" => {
                                let raw = parse_option_path(
                                    tokens,
                                    &mut idx,
                                    "--scope",
                                    "planning profile set",
                                )?;
                                scope = parse_planning_profile_scope(&raw)?;
                            }
                            other => {
                                return Err(format!(
                                    "Unknown option '{other}' for planning profile set"
                                ));
                            }
                        }
                    }
                    if scope == PlanningProfileScope::Effective {
                        return Err(
                            "planning profile set does not support --scope effective".to_string()
                        );
                    }
                    Ok(ShellCommand::PlanningProfileSet {
                        scope,
                        payload_json,
                    })
                }
                "clear" => {
                    let mut scope = PlanningProfileScope::ProjectOverride;
                    let mut idx = 3usize;
                    while idx < tokens.len() {
                        match tokens[idx].as_str() {
                            "--scope" => {
                                let raw = parse_option_path(
                                    tokens,
                                    &mut idx,
                                    "--scope",
                                    "planning profile clear",
                                )?;
                                scope = parse_planning_profile_scope(&raw)?;
                            }
                            other => {
                                return Err(format!(
                                    "Unknown option '{other}' for planning profile clear"
                                ));
                            }
                        }
                    }
                    if scope == PlanningProfileScope::Effective {
                        return Err(
                            "planning profile clear does not support --scope effective".to_string()
                        );
                    }
                    Ok(ShellCommand::PlanningProfileSet {
                        scope,
                        payload_json: "null".to_string(),
                    })
                }
                other => Err(format!(
                    "Unknown planning profile subcommand '{other}' (expected show, set, clear)"
                )),
            }
        }
        "objective" => {
            if tokens.len() < 3 {
                return Err(
                    "planning objective requires a subcommand: show, set, clear".to_string()
                );
            }
            match tokens[2].as_str() {
                "show" => {
                    if tokens.len() > 3 {
                        return Err("planning objective show takes no options".to_string());
                    }
                    Ok(ShellCommand::PlanningObjectiveShow)
                }
                "set" => {
                    if tokens.len() != 4 {
                        return Err("planning objective set requires JSON_OR_@FILE".to_string());
                    }
                    Ok(ShellCommand::PlanningObjectiveSet {
                        payload_json: tokens[3].clone(),
                    })
                }
                "clear" => {
                    if tokens.len() > 3 {
                        return Err("planning objective clear takes no options".to_string());
                    }
                    Ok(ShellCommand::PlanningObjectiveSet {
                        payload_json: "null".to_string(),
                    })
                }
                other => Err(format!(
                    "Unknown planning objective subcommand '{other}' (expected show, set, clear)"
                )),
            }
        }
        "suggestions" => {
            if tokens.len() < 3 {
                return Err(
                    "planning suggestions requires a subcommand: list, accept, reject".to_string(),
                );
            }
            match tokens[2].as_str() {
                "list" => {
                    let mut status: Option<PlanningSuggestionStatus> = None;
                    let mut idx = 3usize;
                    while idx < tokens.len() {
                        match tokens[idx].as_str() {
                            "--status" => {
                                let raw = parse_option_path(
                                    tokens,
                                    &mut idx,
                                    "--status",
                                    "planning suggestions list",
                                )?;
                                status = Some(parse_planning_suggestion_status(&raw)?);
                            }
                            other => {
                                return Err(format!(
                                    "Unknown option '{other}' for planning suggestions list"
                                ));
                            }
                        }
                    }
                    Ok(ShellCommand::PlanningSuggestionsList { status })
                }
                "accept" => {
                    if tokens.len() != 4 {
                        return Err(
                            "planning suggestions accept requires SUGGESTION_ID".to_string()
                        );
                    }
                    let suggestion_id = tokens[3].trim().to_string();
                    if suggestion_id.is_empty() {
                        return Err(
                            "planning suggestions accept SUGGESTION_ID cannot be empty".to_string()
                        );
                    }
                    Ok(ShellCommand::PlanningSuggestionAccept { suggestion_id })
                }
                "reject" => {
                    if tokens.len() < 4 {
                        return Err(
                            "planning suggestions reject requires SUGGESTION_ID [--reason TEXT]"
                                .to_string(),
                        );
                    }
                    let suggestion_id = tokens[3].trim().to_string();
                    if suggestion_id.is_empty() {
                        return Err(
                            "planning suggestions reject SUGGESTION_ID cannot be empty".to_string()
                        );
                    }
                    let mut reason: Option<String> = None;
                    let mut idx = 4usize;
                    while idx < tokens.len() {
                        match tokens[idx].as_str() {
                            "--reason" => {
                                reason = Some(parse_option_path(
                                    tokens,
                                    &mut idx,
                                    "--reason",
                                    "planning suggestions reject",
                                )?);
                            }
                            other => {
                                return Err(format!(
                                    "Unknown option '{other}' for planning suggestions reject"
                                ));
                            }
                        }
                    }
                    Ok(ShellCommand::PlanningSuggestionReject {
                        suggestion_id,
                        reason,
                    })
                }
                other => Err(format!(
                    "Unknown planning suggestions subcommand '{other}' (expected list, accept, reject)"
                )),
            }
        }
        "sync" => {
            if tokens.len() < 3 {
                return Err("planning sync requires a subcommand: status, pull, push".to_string());
            }
            match tokens[2].as_str() {
                "status" => {
                    if tokens.len() > 3 {
                        return Err("planning sync status takes no options".to_string());
                    }
                    Ok(ShellCommand::PlanningSyncStatus)
                }
                "pull" | "push" => {
                    if tokens.len() < 4 {
                        return Err(format!(
                            "planning sync {} requires JSON_OR_@FILE [--source ID] [--confidence N] [--snapshot-id ID]",
                            tokens[2]
                        ));
                    }
                    let payload_json = tokens[3].clone();
                    let mut source: Option<String> = None;
                    let mut confidence: Option<f64> = None;
                    let mut snapshot_id: Option<String> = None;
                    let mut idx = 4usize;
                    while idx < tokens.len() {
                        match tokens[idx].as_str() {
                            "--source" => {
                                source = Some(parse_option_path(
                                    tokens,
                                    &mut idx,
                                    "--source",
                                    "planning sync",
                                )?);
                            }
                            "--confidence" => {
                                let raw = parse_option_path(
                                    tokens,
                                    &mut idx,
                                    "--confidence",
                                    "planning sync",
                                )?;
                                let parsed = raw.parse::<f64>().map_err(|e| {
                                    format!("Invalid --confidence value '{raw}': {e}")
                                })?;
                                confidence = Some(parsed);
                            }
                            "--snapshot-id" => {
                                snapshot_id = Some(parse_option_path(
                                    tokens,
                                    &mut idx,
                                    "--snapshot-id",
                                    "planning sync",
                                )?);
                            }
                            other => {
                                return Err(format!("Unknown option '{other}' for planning sync"));
                            }
                        }
                    }
                    if tokens[2].eq_ignore_ascii_case("pull") {
                        Ok(ShellCommand::PlanningSyncPull {
                            payload_json,
                            source,
                            confidence,
                            snapshot_id,
                        })
                    } else {
                        Ok(ShellCommand::PlanningSyncPush {
                            payload_json,
                            source,
                            confidence,
                            snapshot_id,
                        })
                    }
                }
                other => Err(format!(
                    "Unknown planning sync subcommand '{other}' (expected status, pull, push)"
                )),
            }
        }
        other => Err(format!(
            "Unknown planning subcommand '{other}' (expected profile, objective, suggestions, sync)"
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
        "export-run-bundle" => {
            if tokens.len() < 2 {
                return Err("export-run-bundle requires: OUTPUT.json [--run-id RUN_ID]".to_string());
            }
            let output = tokens[1].clone();
            let mut run_id: Option<String> = None;
            let mut idx = 2usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--run-id" => {
                        if idx + 1 >= tokens.len() {
                            return Err("Missing value after --run-id".to_string());
                        }
                        let value = tokens[idx + 1].trim();
                        if !value.is_empty() {
                            run_id = Some(value.to_string());
                        }
                        idx += 2;
                    }
                    other => {
                        return Err(format!("Unknown argument '{other}' for export-run-bundle"));
                    }
                }
            }
            Ok(ShellCommand::ExportRunBundle { output, run_id })
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
        "dotplot" => parse_dotplot_command(tokens),
        "flex" => parse_flex_command(tokens),
        "rna-reads" | "rna_reads" | "rnareads" => parse_rna_reads_command(tokens),
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
        "panels" => parse_panels_command(tokens),
        "genbank" => parse_genbank_command(tokens),
        "uniprot" => parse_uniprot_command(tokens),
        "macros" => parse_macros_command(tokens),
        "candidates" => parse_candidates_command(tokens),
        "planning" => parse_planning_command(tokens),
        "guides" => parse_guides_command(tokens),
        "primers" => parse_primers_command(tokens),
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

fn macro_bindings_from_preflight(
    preflight: &MacroTemplatePreflightReport,
    direction: RoutinePortDirection,
) -> Vec<LineageMacroPortBinding> {
    preflight
        .checked_ports
        .iter()
        .filter(|row| row.direction == direction)
        .filter(|row| matches!(row.status, RoutinePortValidationStatus::Ok))
        .map(|row| LineageMacroPortBinding {
            port_id: row.port_id.clone(),
            kind: row.kind.clone(),
            required: row.required,
            cardinality: row.cardinality.clone(),
            values: row.values.clone(),
            description: None,
        })
        .collect::<Vec<_>>()
}

fn macro_run_id_from_operation_records(records: &[crate::engine::OperationRecord]) -> String {
    let mut run_ids = records
        .iter()
        .map(|record| record.run_id.clone())
        .collect::<Vec<_>>();
    run_ids.sort();
    run_ids.dedup();
    if run_ids.is_empty() {
        "macro".to_string()
    } else if run_ids.len() == 1 {
        run_ids[0].clone()
    } else {
        format!("mixed:{}", run_ids.join("+"))
    }
}

fn record_macro_lineage_instance(
    engine: &mut GentleEngine,
    start_operation_index: usize,
    template_name: Option<String>,
    routine_id: Option<String>,
    routine_title: Option<String>,
    bound_inputs: Vec<LineageMacroPortBinding>,
    bound_outputs: Vec<LineageMacroPortBinding>,
    status: MacroInstanceStatus,
    status_message: Option<String>,
) -> String {
    let operation_records = engine
        .operation_log()
        .iter()
        .skip(start_operation_index)
        .cloned()
        .collect::<Vec<_>>();
    let expanded_op_ids = operation_records
        .iter()
        .map(|record| record.result.op_id.clone())
        .collect::<Vec<_>>();
    let run_id = macro_run_id_from_operation_records(&operation_records);
    let instance = LineageMacroInstance {
        macro_instance_id: String::new(),
        routine_id,
        routine_title,
        template_name,
        run_id,
        created_at_unix_ms: 0,
        bound_inputs,
        bound_outputs,
        expanded_op_ids,
        status,
        status_message,
    };
    engine.record_lineage_macro_instance(instance)
}

fn macro_status_from_error(err: &str) -> MacroInstanceStatus {
    let lower = err.to_ascii_lowercase();
    if lower.contains("cancelled")
        || lower.contains("canceled")
        || lower.contains("timeout")
        || lower.contains("interrupted")
    {
        MacroInstanceStatus::Cancelled
    } else {
        MacroInstanceStatus::Failed
    }
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
        ShellCommand::PanelsImportIsoform {
            seq_id,
            panel_path,
            panel_id,
            strict,
        } => {
            let op_result = engine
                .apply(Operation::ImportIsoformPanel {
                    seq_id: seq_id.clone(),
                    panel_path: panel_path.clone(),
                    panel_id: panel_id.clone(),
                    strict: *strict,
                })
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: true,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::PanelsInspectIsoform { seq_id, panel_id } => {
            let view = engine
                .inspect_feature_expert(
                    seq_id,
                    &FeatureExpertTarget::IsoformArchitecture {
                        panel_id: panel_id.clone(),
                    },
                )
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: serde_json::to_value(view)
                    .map_err(|e| format!("Could not serialize isoform architecture view: {e}"))?,
            }
        }
        ShellCommand::PanelsRenderIsoformSvg {
            seq_id,
            panel_id,
            output,
        } => {
            let op_result = engine
                .apply(Operation::RenderIsoformArchitectureSvg {
                    seq_id: seq_id.clone(),
                    panel_id: panel_id.clone(),
                    path: output.clone(),
                })
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::PanelsValidateIsoform {
            panel_path,
            panel_id,
        } => {
            let report =
                GentleEngine::validate_isoform_panel_resource(panel_path, panel_id.as_deref())
                    .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: serde_json::to_value(report)
                    .map_err(|e| format!("Could not serialize isoform validation report: {e}"))?,
            }
        }
        ShellCommand::GenbankFetch { accession, as_id } => {
            let op_result = engine
                .apply(Operation::FetchGenBankAccession {
                    accession: accession.clone(),
                    as_id: as_id.clone(),
                })
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: true,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::UniprotFetch { query, entry_id } => {
            let op_result = engine
                .apply(Operation::FetchUniprotSwissProt {
                    query: query.clone(),
                    entry_id: entry_id.clone(),
                })
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: true,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::UniprotImportSwissProt { path, entry_id } => {
            let op_result = engine
                .apply(Operation::ImportUniprotSwissProt {
                    path: path.clone(),
                    entry_id: entry_id.clone(),
                })
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: true,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::UniprotList => {
            let rows = engine.list_uniprot_entries();
            ShellRunResult {
                state_changed: false,
                output: serde_json::to_value(rows)
                    .map_err(|e| format!("Could not serialize UniProt entry list: {e}"))?,
            }
        }
        ShellCommand::UniprotShow { entry_id } => {
            let entry = engine
                .get_uniprot_entry(entry_id)
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: serde_json::to_value(entry)
                    .map_err(|e| format!("Could not serialize UniProt entry: {e}"))?,
            }
        }
        ShellCommand::UniprotMap {
            entry_id,
            seq_id,
            projection_id,
            transcript_id,
        } => {
            let op_result = engine
                .apply(Operation::ProjectUniprotToGenome {
                    seq_id: seq_id.clone(),
                    entry_id: entry_id.clone(),
                    projection_id: projection_id.clone(),
                    transcript_id: transcript_id.clone(),
                })
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: true,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::UniprotProjectionList { seq_id } => {
            let rows = engine.list_uniprot_genome_projections(seq_id.as_deref());
            ShellRunResult {
                state_changed: false,
                output: serde_json::to_value(rows)
                    .map_err(|e| format!("Could not serialize UniProt projection list: {e}"))?,
            }
        }
        ShellCommand::UniprotProjectionShow { projection_id } => {
            let projection = engine
                .get_uniprot_genome_projection(projection_id)
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: serde_json::to_value(projection)
                    .map_err(|e| format!("Could not serialize UniProt projection: {e}"))?,
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
        ShellCommand::ExportRunBundle { output, run_id } => {
            let op_result = engine
                .apply(Operation::ExportProcessRunBundle {
                    path: output.clone(),
                    run_id: run_id.clone(),
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
            let planning_enabled = engine.planning_meta_enabled();
            let planning_objective = engine.planning_objective();
            let planning_profile = engine.planning_effective_profile();
            let mut planning_rows = routines
                .drain(..)
                .map(|routine| {
                    let estimate = estimate_routine_planning(engine, &routine);
                    let mut payload = serde_json::to_value(&routine).unwrap_or_else(|_| json!({}));
                    if let Some(obj) = payload.as_object_mut() {
                        obj.insert(
                            "estimated_time_hours".to_string(),
                            json!(estimate.estimated_time_hours),
                        );
                        obj.insert("estimated_cost".to_string(), json!(estimate.estimated_cost));
                        obj.insert(
                            "local_fit_score".to_string(),
                            json!(estimate.local_fit_score),
                        );
                        obj.insert(
                            "composite_meta_score".to_string(),
                            json!(estimate.composite_meta_score),
                        );
                        obj.insert(
                            "planning_estimate".to_string(),
                            serde_json::to_value(&estimate).unwrap_or_else(|_| json!({})),
                        );
                    }
                    (routine, estimate, payload)
                })
                .collect::<Vec<_>>();
            if planning_enabled {
                planning_rows.sort_by(
                    |(left_routine, left_estimate, _), (right_routine, right_estimate, _)| {
                        right_estimate
                            .passes_guardrails
                            .cmp(&left_estimate.passes_guardrails)
                            .then_with(|| {
                                right_estimate
                                    .composite_meta_score
                                    .total_cmp(&left_estimate.composite_meta_score)
                            })
                            .then(
                                left_routine
                                    .family
                                    .to_ascii_lowercase()
                                    .cmp(&right_routine.family.to_ascii_lowercase()),
                            )
                            .then(
                                left_routine
                                    .title
                                    .to_ascii_lowercase()
                                    .cmp(&right_routine.title.to_ascii_lowercase()),
                            )
                            .then(
                                left_routine
                                    .routine_id
                                    .to_ascii_lowercase()
                                    .cmp(&right_routine.routine_id.to_ascii_lowercase()),
                            )
                    },
                );
            } else {
                planning_rows.sort_by(|(left, _, _), (right, _, _)| {
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
            }
            let guardrail_blocked_count = planning_rows
                .iter()
                .filter(|(_, estimate, _)| !estimate.passes_guardrails)
                .count();
            let routines_payload = planning_rows
                .into_iter()
                .map(|(_, _, payload)| payload)
                .collect::<Vec<_>>();
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
                    "routine_count": routines_payload.len(),
                    "routines": routines_payload,
                    "planning": {
                        "enabled": planning_enabled,
                        "profile_merge_order": [
                            "global_profile",
                            "confirmed_agent_overlay",
                            "project_override"
                        ],
                        "profile_procurement_business_days_default": planning_profile.procurement_business_days_default,
                        "objective": planning_objective,
                        "guardrail_blocked_count": guardrail_blocked_count,
                        "estimate_schema": PLANNING_ESTIMATE_SCHEMA,
                    },
                }),
            }
        }
        ShellCommand::RoutinesExplain {
            catalog_path,
            routine_id,
        } => {
            let resolved_catalog = catalog_path
                .as_deref()
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .unwrap_or(DEFAULT_CLONING_ROUTINE_CATALOG_PATH)
                .to_string();
            let catalog = load_cloning_routine_catalog(&resolved_catalog)?;
            let catalog_schema = catalog.schema.clone();
            let available_routine_ids = catalog
                .routines
                .iter()
                .map(|row| row.routine_id.clone())
                .collect::<Vec<_>>();
            let Some(routine) = resolve_catalog_routine(&catalog, routine_id) else {
                return Err(format!(
                    "Routine '{}' was not found in catalog '{}'; available routine_id values: {}",
                    routine_id.trim(),
                    resolved_catalog,
                    available_routine_ids.join(", ")
                ));
            };
            let routine = routine.clone();

            let mut alternatives = vec![];
            if routine.confusing_alternatives.is_empty() {
                for alt in catalog
                    .routines
                    .iter()
                    .filter(|candidate| {
                        !candidate
                            .routine_id
                            .eq_ignore_ascii_case(routine.routine_id.as_str())
                            && candidate
                                .family
                                .eq_ignore_ascii_case(routine.family.as_str())
                    })
                    .take(4)
                {
                    alternatives.push(routine_summary_row(alt));
                }
            } else {
                for alt_id in &routine.confusing_alternatives {
                    if let Some(alt) = resolve_catalog_routine(&catalog, alt_id) {
                        alternatives.push(routine_summary_row(alt));
                    }
                }
            }

            let requires = if routine.requires.is_empty() {
                build_default_routine_requirements(&routine)
            } else {
                routine.requires.clone()
            };

            let disambiguation_questions = if routine.disambiguation_questions.is_empty() {
                vec![
                    "What molecule types and termini are expected for insert and vector?"
                        .to_string(),
                    "Do you need directionality constraints or compatible overhang control?"
                        .to_string(),
                    "Is this intended as planning-only preflight or as mutating execution?"
                        .to_string(),
                ]
            } else {
                routine.disambiguation_questions.clone()
            };
            let purpose = routine.purpose.clone().or_else(|| routine.summary.clone());
            let mechanism = routine.mechanism.clone();
            let contraindications = routine.contraindications.clone();
            let failure_modes = routine.failure_modes.clone();
            let routine_payload = routine.clone();

            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": CLONING_ROUTINE_EXPLAIN_SCHEMA,
                    "catalog_path": resolved_catalog,
                    "catalog_schema": catalog_schema,
                    "routine": routine_payload,
                    "explanation": {
                        "purpose": purpose,
                        "mechanism": mechanism,
                        "requires": requires,
                        "contraindications": contraindications,
                        "disambiguation_questions": disambiguation_questions,
                        "failure_modes": failure_modes,
                    },
                    "alternatives": alternatives,
                }),
            }
        }
        ShellCommand::RoutinesCompare {
            catalog_path,
            left_routine_id,
            right_routine_id,
        } => {
            let resolved_catalog = catalog_path
                .as_deref()
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .unwrap_or(DEFAULT_CLONING_ROUTINE_CATALOG_PATH)
                .to_string();
            let catalog = load_cloning_routine_catalog(&resolved_catalog)?;
            let catalog_schema = catalog.schema.clone();

            let available_routine_ids = catalog
                .routines
                .iter()
                .map(|row| row.routine_id.clone())
                .collect::<Vec<_>>();
            let Some(left) = resolve_catalog_routine(&catalog, left_routine_id) else {
                return Err(format!(
                    "Routine '{}' was not found in catalog '{}'; available routine_id values: {}",
                    left_routine_id.trim(),
                    resolved_catalog,
                    available_routine_ids.join(", ")
                ));
            };
            let Some(right) = resolve_catalog_routine(&catalog, right_routine_id) else {
                return Err(format!(
                    "Routine '{}' was not found in catalog '{}'; available routine_id values: {}",
                    right_routine_id.trim(),
                    resolved_catalog,
                    available_routine_ids.join(", ")
                ));
            };
            let left = left.clone();
            let right = right.clone();

            let left_tag_set = left
                .vocabulary_tags
                .iter()
                .map(|value| value.to_ascii_lowercase())
                .collect::<BTreeSet<_>>();
            let right_tag_set = right
                .vocabulary_tags
                .iter()
                .map(|value| value.to_ascii_lowercase())
                .collect::<BTreeSet<_>>();
            let shared_tags = left
                .vocabulary_tags
                .iter()
                .filter(|value| right_tag_set.contains(&value.to_ascii_lowercase()))
                .cloned()
                .collect::<Vec<_>>();
            let left_only_tags = left
                .vocabulary_tags
                .iter()
                .filter(|value| !right_tag_set.contains(&value.to_ascii_lowercase()))
                .cloned()
                .collect::<Vec<_>>();
            let right_only_tags = right
                .vocabulary_tags
                .iter()
                .filter(|value| !left_tag_set.contains(&value.to_ascii_lowercase()))
                .cloned()
                .collect::<Vec<_>>();

            let left_axis_map = build_routine_axis_map(&left);
            let right_axis_map = build_routine_axis_map(&right);
            let mut keys = left_axis_map
                .keys()
                .chain(right_axis_map.keys())
                .cloned()
                .collect::<BTreeSet<_>>();
            if keys.is_empty() {
                keys.insert("core_family".to_string());
                keys.insert("status".to_string());
                keys.insert("template".to_string());
                keys.insert("required_inputs".to_string());
            }
            let mut axis_rows = vec![];
            for key in keys {
                let row = match key.as_str() {
                    "core_family" => json!({
                        "axis": "core_family",
                        "left": left.family.clone(),
                        "right": right.family.clone(),
                        "same": left.family.eq_ignore_ascii_case(right.family.as_str()),
                    }),
                    "status" => json!({
                        "axis": "status",
                        "left": left.status.clone(),
                        "right": right.status.clone(),
                        "same": left.status.eq_ignore_ascii_case(right.status.as_str()),
                    }),
                    "template" => json!({
                        "axis": "template",
                        "left": left.template_name.clone(),
                        "right": right.template_name.clone(),
                        "same": left.template_name.eq_ignore_ascii_case(right.template_name.as_str()),
                    }),
                    "required_inputs" => {
                        let left_inputs = build_default_routine_requirements(&left).join("; ");
                        let right_inputs = build_default_routine_requirements(&right).join("; ");
                        json!({
                            "axis": "required_inputs",
                            "left": left_inputs,
                            "right": right_inputs,
                            "same": left_inputs.eq_ignore_ascii_case(right_inputs.as_str()),
                        })
                    }
                    _ => {
                        let left_value = left_axis_map
                            .get(key.as_str())
                            .map(|(_, value)| value.clone())
                            .unwrap_or_else(|| "-".to_string());
                        let right_value = right_axis_map
                            .get(key.as_str())
                            .map(|(_, value)| value.clone())
                            .unwrap_or_else(|| "-".to_string());
                        let axis_label = left_axis_map
                            .get(key.as_str())
                            .map(|(axis, _)| axis.clone())
                            .or_else(|| {
                                right_axis_map
                                    .get(key.as_str())
                                    .map(|(axis, _)| axis.clone())
                            })
                            .unwrap_or_else(|| key.clone());
                        json!({
                            "axis": axis_label,
                            "left": left_value,
                            "right": right_value,
                            "same": left_value.eq_ignore_ascii_case(right_value.as_str()),
                        })
                    }
                };
                axis_rows.push(row);
            }

            let mut disambiguation_questions = left.disambiguation_questions.clone();
            for question in &right.disambiguation_questions {
                if !disambiguation_questions
                    .iter()
                    .any(|entry| entry.eq_ignore_ascii_case(question))
                {
                    disambiguation_questions.push(question.clone());
                }
            }

            let cross_referenced = left
                .confusing_alternatives
                .iter()
                .any(|entry| entry.eq_ignore_ascii_case(right.routine_id.as_str()))
                || right
                    .confusing_alternatives
                    .iter()
                    .any(|entry| entry.eq_ignore_ascii_case(left.routine_id.as_str()));
            let same_family = left.family.eq_ignore_ascii_case(right.family.as_str());
            let planning_enabled = engine.planning_meta_enabled();
            let planning_objective = engine.planning_objective();
            let planning_profile = engine.planning_effective_profile();
            let left_estimate = estimate_routine_planning(engine, &left);
            let right_estimate = estimate_routine_planning(engine, &right);
            axis_rows.push(json!({
                "axis": "estimated_time_hours",
                "left": format!("{:.2}", left_estimate.estimated_time_hours),
                "right": format!("{:.2}", right_estimate.estimated_time_hours),
                "same": (left_estimate.estimated_time_hours - right_estimate.estimated_time_hours).abs() < 1e-9,
            }));
            axis_rows.push(json!({
                "axis": "estimated_cost",
                "left": format!("{:.2}", left_estimate.estimated_cost),
                "right": format!("{:.2}", right_estimate.estimated_cost),
                "same": (left_estimate.estimated_cost - right_estimate.estimated_cost).abs() < 1e-9,
            }));
            axis_rows.push(json!({
                "axis": "local_fit_score",
                "left": format!("{:.3}", left_estimate.local_fit_score),
                "right": format!("{:.3}", right_estimate.local_fit_score),
                "same": (left_estimate.local_fit_score - right_estimate.local_fit_score).abs() < 1e-9,
            }));
            axis_rows.push(json!({
                "axis": "composite_meta_score",
                "left": format!("{:.3}", left_estimate.composite_meta_score),
                "right": format!("{:.3}", right_estimate.composite_meta_score),
                "same": (left_estimate.composite_meta_score - right_estimate.composite_meta_score).abs() < 1e-9,
            }));
            let preferred_routine_id =
                if left_estimate.composite_meta_score >= right_estimate.composite_meta_score {
                    left.routine_id.clone()
                } else {
                    right.routine_id.clone()
                };
            let left_payload = left.clone();
            let right_payload = right.clone();

            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": CLONING_ROUTINE_COMPARE_SCHEMA,
                    "catalog_path": resolved_catalog,
                    "catalog_schema": catalog_schema,
                    "left": left_payload,
                    "right": right_payload,
                    "comparison": {
                        "same_family": same_family,
                        "cross_referenced_as_alternatives": cross_referenced,
                        "shared_vocabulary_tags": shared_tags,
                        "left_only_tags": left_only_tags,
                        "right_only_tags": right_only_tags,
                        "difference_matrix": axis_rows,
                        "disambiguation_questions": disambiguation_questions,
                    },
                    "planning": {
                        "enabled": planning_enabled,
                        "profile_procurement_business_days_default": planning_profile.procurement_business_days_default,
                        "objective": planning_objective,
                        "left_estimate": left_estimate.clone(),
                        "right_estimate": right_estimate.clone(),
                        "preferred_routine_id": preferred_routine_id,
                        "estimate_schema": PLANNING_ESTIMATE_SCHEMA,
                    },
                }),
            }
        }
        ShellCommand::PlanningProfileShow { scope } => {
            let scope_profile = engine.planning_profile(*scope);
            let effective_profile = engine.planning_effective_profile();
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": "gentle.planning_profile_view.v1",
                    "profile_schema": PLANNING_PROFILE_SCHEMA,
                    "scope": scope.as_str(),
                    "profile": scope_profile,
                    "effective_profile": effective_profile,
                    "profile_merge_order": [
                        "global_profile",
                        "confirmed_agent_overlay",
                        "project_override"
                    ],
                }),
            }
        }
        ShellCommand::PlanningProfileSet {
            scope,
            payload_json,
        } => {
            let profile =
                parse_optional_json_payload::<PlanningProfile>(payload_json, "planning profile")?;
            engine
                .set_planning_profile(*scope, profile)
                .map_err(|e| e.to_string())?;
            let scope_profile = engine.planning_profile(*scope);
            let effective_profile = engine.planning_effective_profile();
            ShellRunResult {
                state_changed: true,
                output: json!({
                    "schema": "gentle.planning_profile_update.v1",
                    "profile_schema": PLANNING_PROFILE_SCHEMA,
                    "scope": scope.as_str(),
                    "profile": scope_profile,
                    "effective_profile": effective_profile,
                }),
            }
        }
        ShellCommand::PlanningObjectiveShow => {
            let objective = engine.planning_objective();
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": "gentle.planning_objective_view.v1",
                    "objective_schema": PLANNING_OBJECTIVE_SCHEMA,
                    "objective": objective,
                }),
            }
        }
        ShellCommand::PlanningObjectiveSet { payload_json } => {
            let objective = parse_optional_json_payload::<PlanningObjective>(
                payload_json,
                "planning objective",
            )?;
            engine
                .set_planning_objective(objective)
                .map_err(|e| e.to_string())?;
            let current = engine.planning_objective();
            ShellRunResult {
                state_changed: true,
                output: json!({
                    "schema": "gentle.planning_objective_update.v1",
                    "objective_schema": PLANNING_OBJECTIVE_SCHEMA,
                    "objective": current,
                }),
            }
        }
        ShellCommand::PlanningSuggestionsList { status } => {
            let suggestions = engine.list_planning_suggestions(*status);
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": "gentle.planning_suggestions_list.v1",
                    "suggestion_schema": PLANNING_SUGGESTION_SCHEMA,
                    "status_filter": status.map(|value| value.as_str()),
                    "suggestion_count": suggestions.len(),
                    "suggestions": suggestions,
                }),
            }
        }
        ShellCommand::PlanningSuggestionAccept { suggestion_id } => {
            let suggestion = engine
                .accept_planning_suggestion(suggestion_id)
                .map_err(|e| e.to_string())?;
            let sync_status = engine.planning_sync_status();
            ShellRunResult {
                state_changed: true,
                output: json!({
                    "schema": "gentle.planning_suggestion_resolution.v1",
                    "suggestion_schema": PLANNING_SUGGESTION_SCHEMA,
                    "sync_status_schema": PLANNING_SYNC_STATUS_SCHEMA,
                    "status": "accepted",
                    "suggestion": suggestion,
                    "effective_profile": engine.planning_effective_profile(),
                    "objective": engine.planning_objective(),
                    "sync_status": sync_status,
                }),
            }
        }
        ShellCommand::PlanningSuggestionReject {
            suggestion_id,
            reason,
        } => {
            let suggestion = engine
                .reject_planning_suggestion(suggestion_id, reason.as_deref())
                .map_err(|e| e.to_string())?;
            let sync_status = engine.planning_sync_status();
            ShellRunResult {
                state_changed: true,
                output: json!({
                    "schema": "gentle.planning_suggestion_resolution.v1",
                    "suggestion_schema": PLANNING_SUGGESTION_SCHEMA,
                    "sync_status_schema": PLANNING_SYNC_STATUS_SCHEMA,
                    "status": "rejected",
                    "suggestion": suggestion,
                    "effective_profile": engine.planning_effective_profile(),
                    "objective": engine.planning_objective(),
                    "sync_status": sync_status,
                }),
            }
        }
        ShellCommand::PlanningSyncStatus => {
            let status = engine.planning_sync_status();
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": "gentle.planning_sync_status_view.v1",
                    "sync_status_schema": PLANNING_SYNC_STATUS_SCHEMA,
                    "status": status,
                }),
            }
        }
        ShellCommand::PlanningSyncPull {
            payload_json,
            source,
            confidence,
            snapshot_id,
        } => {
            if let Some(value) = confidence {
                if !value.is_finite() || !(0.0..=1.0).contains(value) {
                    return Err(format!(
                        "Invalid planning sync pull confidence '{}'; expected 0.0..=1.0",
                        value
                    ));
                }
            }
            let payload = parse_required_json_payload::<PlanningSyncSuggestionPayload>(
                payload_json,
                "planning sync pull",
            )?;
            let source = source
                .as_deref()
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .unwrap_or("agent_pull");
            let suggestion = engine
                .propose_planning_suggestion(
                    "pull",
                    source,
                    *confidence,
                    snapshot_id.as_deref(),
                    payload.profile_patch,
                    payload.objective_patch,
                    payload.message.as_deref(),
                )
                .map_err(|e| e.to_string())?;
            let sync_status = engine.planning_sync_status();
            ShellRunResult {
                state_changed: true,
                output: json!({
                    "schema": "gentle.planning_sync_suggestion.v1",
                    "direction": "pull",
                    "suggestion_schema": PLANNING_SUGGESTION_SCHEMA,
                    "sync_status_schema": PLANNING_SYNC_STATUS_SCHEMA,
                    "suggestion": suggestion,
                    "sync_status": sync_status,
                }),
            }
        }
        ShellCommand::PlanningSyncPush {
            payload_json,
            source,
            confidence,
            snapshot_id,
        } => {
            if let Some(value) = confidence {
                if !value.is_finite() || !(0.0..=1.0).contains(value) {
                    return Err(format!(
                        "Invalid planning sync push confidence '{}'; expected 0.0..=1.0",
                        value
                    ));
                }
            }
            let payload = parse_required_json_payload::<PlanningSyncSuggestionPayload>(
                payload_json,
                "planning sync push",
            )?;
            let source = source
                .as_deref()
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .unwrap_or("agent_push");
            let suggestion = engine
                .propose_planning_suggestion(
                    "push",
                    source,
                    *confidence,
                    snapshot_id.as_deref(),
                    payload.profile_patch,
                    payload.objective_patch,
                    payload.message.as_deref(),
                )
                .map_err(|e| e.to_string())?;
            let sync_status = engine.planning_sync_status();
            ShellRunResult {
                state_changed: true,
                output: json!({
                    "schema": "gentle.planning_sync_suggestion.v1",
                    "direction": "push",
                    "suggestion_schema": PLANNING_SUGGESTION_SCHEMA,
                    "sync_status_schema": PLANNING_SYNC_STATUS_SCHEMA,
                    "suggestion": suggestion,
                    "sync_status": sync_status,
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
            max_hits_explicit,
            task,
            request_options_json,
            catalog_path,
            cache_dir,
        } => {
            let resolved_catalog = resolved_catalog_path(catalog_path, *helper_mode);
            let report = if *helper_mode {
                engine.blast_helper_genome_with_project_and_request_options(
                    genome_id,
                    query_sequence,
                    request_options_json.as_ref(),
                    task.as_deref(),
                    if *max_hits_explicit {
                        Some(*max_hits)
                    } else {
                        None
                    },
                    resolved_catalog,
                    cache_dir.as_deref(),
                )
            } else {
                engine.blast_reference_genome_with_project_and_request_options(
                    resolved_catalog,
                    genome_id,
                    query_sequence,
                    request_options_json.as_ref(),
                    task.as_deref(),
                    if *max_hits_explicit {
                        Some(*max_hits)
                    } else {
                        None
                    },
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
        ShellCommand::ReferenceBlastAsyncStart {
            helper_mode,
            genome_id,
            query_sequence,
            max_hits,
            max_hits_explicit,
            task,
            request_options_json,
            catalog_path,
            cache_dir,
        } => {
            let status = start_blast_async_job(
                engine,
                *helper_mode,
                genome_id,
                query_sequence,
                *max_hits,
                *max_hits_explicit,
                task.clone(),
                request_options_json.clone(),
                catalog_path.clone(),
                cache_dir.clone(),
            )?;
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": "gentle.blast_async_start.v1",
                    "job": status,
                }),
            }
        }
        ShellCommand::ReferenceBlastAsyncStatus {
            helper_mode: _,
            job_id,
            include_report,
        } => {
            let (status, report) = get_blast_async_job_snapshot(job_id, *include_report)?;
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": "gentle.blast_async_status.v1",
                    "job": status,
                    "report": report,
                }),
            }
        }
        ShellCommand::ReferenceBlastAsyncCancel {
            helper_mode: _,
            job_id,
        } => {
            let status = cancel_blast_async_job(job_id)?;
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": "gentle.blast_async_cancel.v1",
                    "job": status,
                }),
            }
        }
        ShellCommand::ReferenceBlastAsyncList { helper_mode } => {
            let jobs = collect_blast_async_job_snapshots(Some(*helper_mode))?;
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": "gentle.blast_async_list.v1",
                    "helper_mode": helper_mode,
                    "job_count": jobs.len(),
                    "jobs": jobs,
                }),
            }
        }
        ShellCommand::ReferenceBlastTrack {
            helper_mode,
            genome_id,
            query_sequence,
            target_seq_id,
            max_hits,
            max_hits_explicit,
            task,
            request_options_json,
            track_name,
            clear_existing,
            catalog_path,
            cache_dir,
        } => {
            let resolved_catalog = resolved_catalog_path(catalog_path, *helper_mode);
            let report = if *helper_mode {
                engine.blast_helper_genome_with_project_and_request_options(
                    genome_id,
                    query_sequence,
                    request_options_json.as_ref(),
                    task.as_deref(),
                    if *max_hits_explicit {
                        Some(*max_hits)
                    } else {
                        None
                    },
                    resolved_catalog,
                    cache_dir.as_deref(),
                )
            } else {
                engine.blast_reference_genome_with_project_and_request_options(
                    resolved_catalog,
                    genome_id,
                    query_sequence,
                    request_options_json.as_ref(),
                    task.as_deref(),
                    if *max_hits_explicit {
                        Some(*max_hits)
                    } else {
                        None
                    },
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
            let command_line = if report.command.is_empty() {
                report.blastn_executable.clone()
            } else {
                format!("{} {}", report.blastn_executable, report.command.join(" "))
            };
            let op_result = engine
                .apply(Operation::ImportBlastHitsTrack {
                    seq_id: target_seq_id.clone(),
                    hits: hit_inputs,
                    track_name: track_name.clone(),
                    clear_existing: Some(*clear_existing),
                    blast_provenance: Some(crate::engine::BlastInvocationProvenance {
                        genome_id: report.genome_id.clone(),
                        query_label: target_seq_id.clone(),
                        query_length: report.query_length,
                        max_hits: report.max_hits,
                        task: report.task.clone(),
                        blastn_executable: report.blastn_executable.clone(),
                        blast_db_prefix: report.blast_db_prefix.clone(),
                        command: report.command.clone(),
                        command_line,
                        catalog_path: resolved_catalog.map(str::to_string),
                        cache_dir: cache_dir.clone(),
                        options_override_json: report.options_override_json.clone(),
                        effective_options_json: report.effective_options_json.clone(),
                    }),
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
            annotation_scope,
            max_annotation_features,
            include_genomic_annotation,
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
                    annotation_scope: *annotation_scope,
                    max_annotation_features: *max_annotation_features,
                    include_genomic_annotation: *include_genomic_annotation,
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
            annotation_scope,
            max_annotation_features,
            include_genomic_annotation,
            catalog_path,
            cache_dir,
        } => {
            let op_result = engine
                .apply(Operation::ExtractGenomeGene {
                    genome_id: genome_id.clone(),
                    gene_query: gene_query.clone(),
                    occurrence: *occurrence,
                    output_id: output_id.clone(),
                    annotation_scope: *annotation_scope,
                    max_annotation_features: *max_annotation_features,
                    include_genomic_annotation: *include_genomic_annotation,
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
            prepared_genome_id,
        } => {
            let op_result = engine
                .apply(Operation::ExtendGenomeAnchor {
                    seq_id: seq_id.clone(),
                    side: *side,
                    length_bp: *length_bp,
                    output_id: output_id.clone(),
                    catalog_path: operation_catalog_path(catalog_path, *helper_mode),
                    cache_dir: cache_dir.clone(),
                    prepared_genome_id: prepared_genome_id.clone(),
                })
                .map_err(|e| e.to_string())?;
            let state_changed =
                !op_result.created_seq_ids.is_empty() || !op_result.changed_seq_ids.is_empty();
            ShellRunResult {
                state_changed,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::ReferenceVerifyAnchor {
            helper_mode,
            seq_id,
            catalog_path,
            cache_dir,
            prepared_genome_id,
        } => {
            let op_result = engine
                .apply(Operation::VerifyGenomeAnchor {
                    seq_id: seq_id.clone(),
                    catalog_path: operation_catalog_path(catalog_path, *helper_mode),
                    cache_dir: cache_dir.clone(),
                    prepared_genome_id: prepared_genome_id.clone(),
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
        } => {
            let start_operation_index = engine.operation_log().len();
            match run_workflow_macro(engine, script, *transactional, options) {
                Ok(mut run) => {
                    let macro_instance_id = record_macro_lineage_instance(
                        engine,
                        start_operation_index,
                        None,
                        None,
                        None,
                        vec![],
                        vec![],
                        MacroInstanceStatus::Ok,
                        None,
                    );
                    if let Some(map) = run.output.as_object_mut() {
                        map.insert("macro_instance_id".to_string(), json!(macro_instance_id));
                        map.insert("macro_recorded".to_string(), json!(true));
                    }
                    run.state_changed = true;
                    run
                }
                Err(err) => {
                    let status = macro_status_from_error(&err);
                    let macro_instance_id = record_macro_lineage_instance(
                        engine,
                        start_operation_index,
                        None,
                        None,
                        None,
                        vec![],
                        vec![],
                        status,
                        Some(err.clone()),
                    );
                    return Err(format!("{} [macro_instance_id={}]", err, macro_instance_id));
                }
            }
        }
        ShellCommand::MacrosInstanceList => {
            let instances = engine.lineage_macro_instances().to_vec();
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": "gentle.lineage_macro_instances.v1",
                    "count": instances.len(),
                    "instances": instances
                }),
            }
        }
        ShellCommand::MacrosInstanceShow { macro_instance_id } => {
            let instance = engine
                .lineage_macro_instances()
                .iter()
                .find(|instance| instance.macro_instance_id == *macro_instance_id)
                .cloned()
                .ok_or_else(|| format!("Macro instance '{}' was not found", macro_instance_id))?;
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": "gentle.lineage_macro_instance.v1",
                    "instance": instance
                }),
            }
        }
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
            input_ports,
            output_ports,
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
                    input_ports: input_ports.clone(),
                    output_ports: output_ports.clone(),
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
                    "input_port_count": input_ports.len(),
                    "output_port_count": output_ports.len(),
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
                    input_ports: template.input_ports.clone(),
                    output_ports: template.output_ports.clone(),
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
            validate_only,
        } => {
            let preflight = preflight_workflow_macro_template_run(engine, name, bindings)?;
            if *validate_only {
                ShellRunResult {
                    state_changed: false,
                    output: json!({
                        "template_name": name,
                        "bindings": bindings,
                        "validate_only": true,
                        "can_execute": preflight.can_execute(),
                        "preflight": preflight,
                    }),
                }
            } else {
                if !preflight.can_execute() {
                    let macro_instance_id = record_macro_lineage_instance(
                        engine,
                        engine.operation_log().len(),
                        Some(name.clone()),
                        preflight.routine_id.clone(),
                        preflight.routine_title.clone(),
                        macro_bindings_from_preflight(&preflight, RoutinePortDirection::Input),
                        macro_bindings_from_preflight(&preflight, RoutinePortDirection::Output),
                        MacroInstanceStatus::Failed,
                        Some(format!("preflight failed: {}", preflight.errors.join("; "))),
                    );
                    return Err(format!(
                        "Workflow macro template '{}' preflight failed: {} [macro_instance_id={}]",
                        name,
                        preflight.errors.join("; "),
                        macro_instance_id
                    ));
                }
                let script = engine
                    .render_workflow_macro_template_script(name, bindings)
                    .map_err(|e| e.to_string())?;
                let start_operation_index = engine.operation_log().len();
                match run_workflow_macro(engine, &script, *transactional, options) {
                    Ok(mut run) => {
                        let macro_instance_id = record_macro_lineage_instance(
                            engine,
                            start_operation_index,
                            Some(name.clone()),
                            preflight.routine_id.clone(),
                            preflight.routine_title.clone(),
                            macro_bindings_from_preflight(&preflight, RoutinePortDirection::Input),
                            macro_bindings_from_preflight(&preflight, RoutinePortDirection::Output),
                            MacroInstanceStatus::Ok,
                            None,
                        );
                        run.state_changed = true;
                        run.output = json!({
                            "template_name": name,
                            "bindings": bindings,
                            "expanded_script": script,
                            "macro_instance_id": macro_instance_id,
                            "preflight": preflight,
                            "run": run.output
                        });
                        run
                    }
                    Err(err) => {
                        let status = macro_status_from_error(&err);
                        let macro_instance_id = record_macro_lineage_instance(
                            engine,
                            start_operation_index,
                            Some(name.clone()),
                            preflight.routine_id.clone(),
                            preflight.routine_title.clone(),
                            macro_bindings_from_preflight(&preflight, RoutinePortDirection::Input),
                            macro_bindings_from_preflight(&preflight, RoutinePortDirection::Output),
                            status,
                            Some(err.clone()),
                        );
                        return Err(format!("{} [macro_instance_id={}]", err, macro_instance_id));
                    }
                }
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
        ShellCommand::PrimersSeedFromFeature { seq_id, feature_id } => {
            let dna = engine
                .state()
                .sequences
                .get(seq_id)
                .ok_or_else(|| format!("Sequence '{seq_id}' not found"))?;
            let (roi_start_0based, roi_end_0based_exclusive) =
                sequence_feature_roi_range_0based(dna, *feature_id)?;
            let primer_pairs = build_seeded_primer_pair_operation(
                seq_id,
                roi_start_0based,
                roi_end_0based_exclusive,
            );
            let qpcr =
                build_seeded_qpcr_operation(seq_id, roi_start_0based, roi_end_0based_exclusive);
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": "gentle.primer_seed_request.v1",
                    "template": seq_id,
                    "source": {
                        "kind": "feature",
                        "feature_id": feature_id,
                    },
                    "roi_start_0based": roi_start_0based,
                    "roi_end_0based_exclusive": roi_end_0based_exclusive,
                    "operations": {
                        "design_primer_pairs": primer_pairs,
                        "design_qpcr_assays": qpcr,
                    }
                }),
            }
        }
        ShellCommand::PrimersSeedFromSplicing { seq_id, feature_id } => {
            let expert = engine
                .inspect_feature_expert(
                    seq_id,
                    &FeatureExpertTarget::SplicingFeature {
                        feature_id: *feature_id,
                        scope: SplicingScopePreset::AllOverlappingBothStrands,
                    },
                )
                .map_err(|e| e.to_string())?;
            let splicing = match expert {
                FeatureExpertView::Splicing(view) => view,
                _ => {
                    return Err(format!(
                        "Feature n-{} on '{}' does not resolve to a splicing expert view",
                        feature_id, seq_id
                    ));
                }
            };
            if splicing.region_start_1based == 0
                || splicing.region_end_1based < splicing.region_start_1based
            {
                return Err(format!(
                    "Splicing region bounds are invalid for feature n-{} on '{}'",
                    feature_id, seq_id
                ));
            }
            let roi_start_0based = splicing.region_start_1based.saturating_sub(1);
            let roi_end_0based_exclusive = splicing.region_end_1based;
            let primer_pairs = build_seeded_primer_pair_operation(
                seq_id,
                roi_start_0based,
                roi_end_0based_exclusive,
            );
            let qpcr =
                build_seeded_qpcr_operation(seq_id, roi_start_0based, roi_end_0based_exclusive);
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": "gentle.primer_seed_request.v1",
                    "template": seq_id,
                    "source": {
                        "kind": "splicing",
                        "feature_id": feature_id,
                        "group_label": splicing.group_label,
                        "transcript_count": splicing.transcript_count,
                        "unique_exon_count": splicing.unique_exon_count,
                    },
                    "roi_start_0based": roi_start_0based,
                    "roi_end_0based_exclusive": roi_end_0based_exclusive,
                    "operations": {
                        "design_primer_pairs": primer_pairs,
                        "design_qpcr_assays": qpcr,
                    }
                }),
            }
        }
        ShellCommand::PrimersDesign {
            request_json,
            backend,
            primer3_executable,
        } => {
            let json_text = parse_json_payload(request_json)?;
            let op: Operation = serde_json::from_str(&json_text).map_err(|e| {
                format!(
                    "Invalid primers design request JSON: {} (expected Operation payload with DesignPrimerPairs)",
                    e
                )
            })?;
            let (template_id, requested_report_id) = match &op {
                Operation::DesignPrimerPairs {
                    template,
                    report_id,
                    ..
                } => (template.clone(), report_id.clone()),
                _ => {
                    return Err(
                        "primers design expects an Operation payload with DesignPrimerPairs"
                            .to_string(),
                    );
                }
            };
            let before = engine
                .state()
                .metadata
                .get(PRIMER_DESIGN_REPORTS_METADATA_KEY)
                .cloned();
            let previous_backend = engine.state().parameters.primer_design_backend;
            let previous_executable = engine.state().parameters.primer3_executable.clone();
            if let Some(override_backend) = backend {
                engine.state_mut().parameters.primer_design_backend = *override_backend;
            }
            if let Some(override_exec) = primer3_executable
                .as_deref()
                .map(str::trim)
                .filter(|v| !v.is_empty())
            {
                engine.state_mut().parameters.primer3_executable = override_exec.to_string();
            }
            let op_result = engine.apply(op).map_err(|e| e.to_string());
            engine.state_mut().parameters.primer_design_backend = previous_backend;
            engine.state_mut().parameters.primer3_executable = previous_executable;
            let op_result = op_result?;
            let after = engine
                .state()
                .metadata
                .get(PRIMER_DESIGN_REPORTS_METADATA_KEY)
                .cloned();
            let reports = engine.list_primer_design_reports();
            let selected_report = if let Some(report_id) = requested_report_id {
                engine.get_primer_design_report(&report_id).ok()
            } else {
                reports
                    .iter()
                    .filter(|summary| summary.template == template_id)
                    .max_by_key(|summary| summary.generated_at_unix_ms)
                    .and_then(|summary| engine.get_primer_design_report(&summary.report_id).ok())
            };
            let effective_backend = selected_report
                .as_ref()
                .map(|report| report.backend.clone());
            ShellRunResult {
                state_changed: before != after,
                output: json!({
                    "result": op_result,
                    "report": selected_report,
                    "effective_backend": effective_backend,
                }),
            }
        }
        ShellCommand::PrimersDesignQpcr {
            request_json,
            backend,
            primer3_executable,
        } => {
            let json_text = parse_json_payload(request_json)?;
            let op: Operation = serde_json::from_str(&json_text).map_err(|e| {
                format!(
                    "Invalid primers design-qpcr request JSON: {} (expected Operation payload with DesignQpcrAssays)",
                    e
                )
            })?;
            let (template_id, requested_report_id) = match &op {
                Operation::DesignQpcrAssays {
                    template,
                    report_id,
                    ..
                } => (template.clone(), report_id.clone()),
                _ => {
                    return Err(
                        "primers design-qpcr expects an Operation payload with DesignQpcrAssays"
                            .to_string(),
                    );
                }
            };
            let before = engine
                .state()
                .metadata
                .get(PRIMER_DESIGN_REPORTS_METADATA_KEY)
                .cloned();
            let previous_backend = engine.state().parameters.primer_design_backend;
            let previous_executable = engine.state().parameters.primer3_executable.clone();
            if let Some(override_backend) = backend {
                engine.state_mut().parameters.primer_design_backend = *override_backend;
            }
            if let Some(override_exec) = primer3_executable
                .as_deref()
                .map(str::trim)
                .filter(|v| !v.is_empty())
            {
                engine.state_mut().parameters.primer3_executable = override_exec.to_string();
            }
            let op_result = engine.apply(op).map_err(|e| e.to_string());
            engine.state_mut().parameters.primer_design_backend = previous_backend;
            engine.state_mut().parameters.primer3_executable = previous_executable;
            let op_result = op_result?;
            let after = engine
                .state()
                .metadata
                .get(PRIMER_DESIGN_REPORTS_METADATA_KEY)
                .cloned();
            let reports = engine.list_qpcr_design_reports();
            let selected_report = if let Some(report_id) = requested_report_id {
                engine.get_qpcr_design_report(&report_id).ok()
            } else {
                reports
                    .iter()
                    .filter(|summary| summary.template == template_id)
                    .max_by_key(|summary| summary.generated_at_unix_ms)
                    .and_then(|summary| engine.get_qpcr_design_report(&summary.report_id).ok())
            };
            let effective_backend = selected_report
                .as_ref()
                .map(|report| report.backend.clone());
            ShellRunResult {
                state_changed: before != after,
                output: json!({
                    "result": op_result,
                    "report": selected_report,
                    "effective_backend": effective_backend,
                }),
            }
        }
        ShellCommand::PrimersPreflight {
            backend,
            primer3_executable,
        } => {
            let report = engine.primer3_preflight_report(*backend, primer3_executable.as_deref());
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": "gentle.primer3_preflight.v1",
                    "preflight": report,
                }),
            }
        }
        ShellCommand::PrimersListReports => {
            let reports = engine.list_primer_design_reports();
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": "gentle.primer_design_report_list.v1",
                    "report_count": reports.len(),
                    "reports": reports,
                }),
            }
        }
        ShellCommand::PrimersShowReport { report_id } => {
            let report = engine
                .get_primer_design_report(report_id)
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "report": report,
                }),
            }
        }
        ShellCommand::PrimersExportReport { report_id, path } => {
            let report = engine
                .export_primer_design_report(report_id, path)
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": "gentle.primer_design_report_export.v1",
                    "report_id": report.report_id,
                    "path": path,
                    "pair_count": report.pair_count,
                }),
            }
        }
        ShellCommand::PrimersListQpcrReports => {
            let reports = engine.list_qpcr_design_reports();
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": "gentle.qpcr_design_report_list.v1",
                    "report_count": reports.len(),
                    "reports": reports,
                }),
            }
        }
        ShellCommand::PrimersShowQpcrReport { report_id } => {
            let report = engine
                .get_qpcr_design_report(report_id)
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "report": report,
                }),
            }
        }
        ShellCommand::PrimersExportQpcrReport { report_id, path } => {
            let report = engine
                .export_qpcr_design_report(report_id, path)
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": "gentle.qpcr_design_report_export.v1",
                    "report_id": report.report_id,
                    "path": path,
                    "assay_count": report.assay_count,
                }),
            }
        }
        ShellCommand::DotplotCompute {
            seq_id,
            span_start_0based,
            span_end_0based,
            mode,
            word_size,
            step_bp,
            max_mismatches,
            tile_bp,
            dotplot_id,
        } => {
            let before = engine
                .state()
                .metadata
                .get(DOTPLOT_ANALYSIS_METADATA_KEY)
                .cloned();
            let op_result = engine
                .apply(Operation::ComputeDotplot {
                    seq_id: seq_id.clone(),
                    span_start_0based: *span_start_0based,
                    span_end_0based: *span_end_0based,
                    mode: *mode,
                    word_size: *word_size,
                    step_bp: *step_bp,
                    max_mismatches: *max_mismatches,
                    tile_bp: *tile_bp,
                    store_as: dotplot_id.clone(),
                })
                .map_err(|e| e.to_string())?;
            let after = engine
                .state()
                .metadata
                .get(DOTPLOT_ANALYSIS_METADATA_KEY)
                .cloned();
            let selected = if let Some(id) = dotplot_id
                .as_deref()
                .map(str::trim)
                .filter(|value| !value.is_empty())
            {
                engine.get_dotplot_view(id).ok()
            } else {
                engine
                    .list_dotplot_views(Some(seq_id.as_str()))
                    .into_iter()
                    .max_by_key(|row| row.generated_at_unix_ms)
                    .and_then(|row| engine.get_dotplot_view(row.dotplot_id.as_str()).ok())
            };
            ShellRunResult {
                state_changed: before != after,
                output: json!({
                    "result": op_result,
                    "dotplot": selected,
                }),
            }
        }
        ShellCommand::DotplotList { seq_id } => {
            let rows = engine.list_dotplot_views(seq_id.as_deref());
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": "gentle.dotplot_view_list.v1",
                    "dotplot_count": rows.len(),
                    "dotplots": rows,
                }),
            }
        }
        ShellCommand::DotplotShow { dotplot_id } => {
            let view = engine
                .get_dotplot_view(dotplot_id)
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "dotplot": view,
                }),
            }
        }
        ShellCommand::FlexCompute {
            seq_id,
            span_start_0based,
            span_end_0based,
            model,
            bin_bp,
            smoothing_bp,
            track_id,
        } => {
            let before = engine
                .state()
                .metadata
                .get(DOTPLOT_ANALYSIS_METADATA_KEY)
                .cloned();
            let op_result = engine
                .apply(Operation::ComputeFlexibilityTrack {
                    seq_id: seq_id.clone(),
                    span_start_0based: *span_start_0based,
                    span_end_0based: *span_end_0based,
                    model: *model,
                    bin_bp: *bin_bp,
                    smoothing_bp: *smoothing_bp,
                    store_as: track_id.clone(),
                })
                .map_err(|e| e.to_string())?;
            let after = engine
                .state()
                .metadata
                .get(DOTPLOT_ANALYSIS_METADATA_KEY)
                .cloned();
            let selected = if let Some(id) = track_id
                .as_deref()
                .map(str::trim)
                .filter(|value| !value.is_empty())
            {
                engine.get_flexibility_track(id).ok()
            } else {
                engine
                    .list_flexibility_tracks(Some(seq_id.as_str()))
                    .into_iter()
                    .max_by_key(|row| row.generated_at_unix_ms)
                    .and_then(|row| engine.get_flexibility_track(row.track_id.as_str()).ok())
            };
            ShellRunResult {
                state_changed: before != after,
                output: json!({
                    "result": op_result,
                    "track": selected,
                }),
            }
        }
        ShellCommand::FlexList { seq_id } => {
            let rows = engine.list_flexibility_tracks(seq_id.as_deref());
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": "gentle.flexibility_track_list.v1",
                    "track_count": rows.len(),
                    "tracks": rows,
                }),
            }
        }
        ShellCommand::FlexShow { track_id } => {
            let track = engine
                .get_flexibility_track(track_id)
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "track": track,
                }),
            }
        }
        ShellCommand::RnaReadsInterpret {
            seq_id,
            seed_feature_id,
            input_path,
            profile,
            input_format,
            scope,
            origin_mode,
            target_gene_ids,
            roi_seed_capture_enabled,
            seed_filter,
            align_config,
            report_id,
            report_mode,
            checkpoint_path,
            checkpoint_every_reads,
            resume_from_checkpoint,
        } => {
            let op_result = engine
                .apply(Operation::InterpretRnaReads {
                    seq_id: seq_id.clone(),
                    seed_feature_id: *seed_feature_id,
                    profile: *profile,
                    input_path: input_path.clone(),
                    input_format: *input_format,
                    scope: *scope,
                    origin_mode: *origin_mode,
                    target_gene_ids: target_gene_ids.clone(),
                    roi_seed_capture_enabled: *roi_seed_capture_enabled,
                    seed_filter: seed_filter.clone(),
                    align_config: align_config.clone(),
                    report_id: report_id.clone(),
                    report_mode: *report_mode,
                    checkpoint_path: checkpoint_path.clone(),
                    checkpoint_every_reads: *checkpoint_every_reads,
                    resume_from_checkpoint: *resume_from_checkpoint,
                })
                .map_err(|e| e.to_string())?;
            let report = if let Some(id) = report_id
                .as_deref()
                .map(str::trim)
                .filter(|value| !value.is_empty())
            {
                engine.get_rna_read_report(id).ok()
            } else {
                engine
                    .list_rna_read_reports(Some(seq_id.as_str()))
                    .into_iter()
                    .max_by_key(|row| row.generated_at_unix_ms)
                    .and_then(|row| engine.get_rna_read_report(row.report_id.as_str()).ok())
            };
            ShellRunResult {
                state_changed: true,
                output: json!({
                    "result": op_result,
                    "report": report,
                }),
            }
        }
        ShellCommand::RnaReadsListReports { seq_id } => {
            let rows = engine.list_rna_read_reports(seq_id.as_deref());
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": "gentle.rna_read_report_list.v1",
                    "report_count": rows.len(),
                    "reports": rows,
                }),
            }
        }
        ShellCommand::RnaReadsShowReport { report_id } => {
            let report = engine
                .get_rna_read_report(report_id)
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "report": report,
                }),
            }
        }
        ShellCommand::RnaReadsExportReport { report_id, path } => {
            let report = engine
                .export_rna_read_report(report_id, path)
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": "gentle.rna_read_report_export.v1",
                    "report_id": report.report_id,
                    "path": path,
                    "read_count_total": report.read_count_total,
                    "read_count_seed_passed": report.read_count_seed_passed,
                    "read_count_aligned": report.read_count_aligned,
                }),
            }
        }
        ShellCommand::RnaReadsExportHitsFasta {
            report_id,
            path,
            selection,
        } => {
            let written = engine
                .export_rna_read_hits_fasta(report_id, path, *selection)
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": "gentle.rna_read_hits_fasta_export.v1",
                    "report_id": report_id,
                    "path": path,
                    "selection": selection.as_str(),
                    "written_records": written,
                }),
            }
        }
        ShellCommand::RnaReadsExportSampleSheet {
            path,
            seq_id,
            report_ids,
            append,
        } => {
            let export = engine
                .export_rna_read_sample_sheet(path, seq_id.as_deref(), report_ids, *append)
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": export.schema,
                    "path": export.path,
                    "report_count": export.report_count,
                    "append": export.appended,
                    "seq_id": seq_id,
                    "report_ids": report_ids,
                }),
            }
        }
        ShellCommand::RnaReadsExportExonPathsTsv {
            report_id,
            path,
            selection,
        } => {
            let export = engine
                .export_rna_read_exon_paths_tsv(report_id, path, *selection)
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": export.schema,
                    "report_id": export.report_id,
                    "path": export.path,
                    "selection": export.selection.as_str(),
                    "row_count": export.row_count,
                }),
            }
        }
        ShellCommand::RnaReadsExportExonAbundanceTsv {
            report_id,
            path,
            selection,
        } => {
            let export = engine
                .export_rna_read_exon_abundance_tsv(report_id, path, *selection)
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": export.schema,
                    "report_id": export.report_id,
                    "path": export.path,
                    "selection": export.selection.as_str(),
                    "selected_read_count": export.selected_read_count,
                    "exon_row_count": export.exon_row_count,
                    "transition_row_count": export.transition_row_count,
                }),
            }
        }
        ShellCommand::RnaReadsExportScoreDensitySvg {
            report_id,
            path,
            scale,
        } => {
            let export = engine
                .export_rna_read_score_density_svg(report_id, path, *scale)
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "schema": export.schema,
                    "report_id": export.report_id,
                    "path": export.path,
                    "scale": export.scale.as_str(),
                    "bin_count": export.bin_count,
                    "max_bin_count": export.max_bin_count,
                    "total_scored_reads": export.total_scored_reads,
                    "derived_from_report_hits_only": export.derived_from_report_hits_only,
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
            let workflow = parse_workflow_json_payload(&json_text)?;
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
    #[cfg(unix)]
    use std::os::unix::fs::PermissionsExt;
    use std::path::Path;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use tempfile::tempdir;

    static JASPAR_RELOAD_TEST_COUNTER: AtomicUsize = AtomicUsize::new(1);
    static BLAST_ASYNC_TEST_MUTEX: std::sync::Mutex<()> = std::sync::Mutex::new(());

    fn with_blast_async_test_overrides<R>(
        max_concurrent: usize,
        worker_delay_ms: u64,
        f: impl FnOnce() -> R,
    ) -> R {
        let _guard = BLAST_ASYNC_TEST_MUTEX
            .lock()
            .expect("blast async test mutex lock");
        let previous_max =
            BLAST_ASYNC_MAX_CONCURRENT_TEST_OVERRIDE.swap(max_concurrent, Ordering::Relaxed);
        let previous_delay =
            BLAST_ASYNC_WORKER_DELAY_MS_TEST_OVERRIDE.swap(worker_delay_ms, Ordering::Relaxed);
        {
            let mut jobs = BLAST_ASYNC_JOBS
                .lock()
                .expect("blast async registry lock before test");
            jobs.clear();
        }
        let result = f();
        {
            let mut jobs = BLAST_ASYNC_JOBS
                .lock()
                .expect("blast async registry lock after test");
            jobs.clear();
        }
        BLAST_ASYNC_WORKER_DELAY_MS_TEST_OVERRIDE.store(previous_delay, Ordering::Relaxed);
        BLAST_ASYNC_MAX_CONCURRENT_TEST_OVERRIDE.store(previous_max, Ordering::Relaxed);
        result
    }

    fn resource_fixture_path(name: &str) -> String {
        format!(
            "{}/test_files/fixtures/resources/{name}",
            env!("CARGO_MANIFEST_DIR")
        )
    }

    fn primer3_fixture_path(name: &str) -> String {
        format!(
            "{}/test_files/fixtures/primer3/{name}",
            env!("CARGO_MANIFEST_DIR")
        )
    }

    #[cfg(unix)]
    fn install_fake_primer3(path: &Path, fixture_path: &Path) -> String {
        let script_path = path.join("fake_primer3.sh");
        let script = format!(
            "#!/bin/sh\nif [ \"$1\" = \"--version\" ]; then\n  echo \"primer3_core synthetic-fixture 2.6.1\"\n  exit 0\nfi\ncat \"{}\"\n",
            fixture_path.display()
        );
        std::fs::write(&script_path, script).expect("write fake primer3");
        let mut perms = std::fs::metadata(&script_path)
            .expect("metadata fake primer3")
            .permissions();
        perms.set_mode(0o755);
        std::fs::set_permissions(&script_path, perms).expect("chmod fake primer3");
        script_path.display().to_string()
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
    fn parse_workflow_json_payload_accepts_raw_workflow() {
        let payload = r#"{ "run_id": "raw", "ops": [] }"#;
        let workflow = parse_workflow_json_payload(payload).expect("parse raw workflow");
        assert_eq!(workflow.run_id, "raw");
        assert!(workflow.ops.is_empty());
    }

    #[test]
    fn parse_workflow_json_payload_accepts_wrapped_example() {
        let payload = r#"{
            "schema": "gentle.workflow_example.v1",
            "id": "demo",
            "title": "Demo",
            "workflow": {
                "run_id": "wrapped",
                "ops": []
            }
        }"#;
        let workflow = parse_workflow_json_payload(payload).expect("parse wrapped workflow");
        assert_eq!(workflow.run_id, "wrapped");
        assert!(workflow.ops.is_empty());
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
    fn parse_export_run_bundle_with_run_id() {
        let cmd = parse_shell_line("export-run-bundle run_bundle.json --run-id demo_run")
            .expect("parse command");
        match cmd {
            ShellCommand::ExportRunBundle { output, run_id } => {
                assert_eq!(output, "run_bundle.json");
                assert_eq!(run_id.as_deref(), Some("demo_run"));
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
                max_hits_explicit,
                task,
                request_options_json,
                catalog_path,
                cache_dir,
            } => {
                assert!(!helper_mode);
                assert_eq!(genome_id, "ToyGenome");
                assert_eq!(query_sequence, "ACGTACGT");
                assert_eq!(max_hits, 12);
                assert!(max_hits_explicit);
                assert_eq!(task.as_deref(), Some("blastn"));
                assert!(request_options_json.is_none());
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
                max_hits_explicit,
                task,
                request_options_json,
                ..
            } => {
                assert!(helper_mode);
                assert_eq!(genome_id, "pUC19");
                assert_eq!(query_sequence, "ACGTAG");
                assert_eq!(max_hits, 25);
                assert!(!max_hits_explicit);
                assert!(task.is_none());
                assert!(request_options_json.is_none());
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_genomes_blast_start_with_options() {
        let cmd = parse_shell_line(
            "genomes blast-start ToyGenome ACGTACGT --max-hits 12 --task blastn --catalog c.json --cache-dir cache",
        )
        .expect("parse command");
        match cmd {
            ShellCommand::ReferenceBlastAsyncStart {
                helper_mode,
                genome_id,
                query_sequence,
                max_hits,
                max_hits_explicit,
                task,
                request_options_json,
                catalog_path,
                cache_dir,
            } => {
                assert!(!helper_mode);
                assert_eq!(genome_id, "ToyGenome");
                assert_eq!(query_sequence, "ACGTACGT");
                assert_eq!(max_hits, 12);
                assert!(max_hits_explicit);
                assert_eq!(task.as_deref(), Some("blastn"));
                assert!(request_options_json.is_none());
                assert_eq!(catalog_path.as_deref(), Some("c.json"));
                assert_eq!(cache_dir.as_deref(), Some("cache"));
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_genomes_blast_status_and_cancel() {
        let status = parse_shell_line("genomes blast-status blast-job-42 --with-report")
            .expect("parse status command");
        match status {
            ShellCommand::ReferenceBlastAsyncStatus {
                helper_mode,
                job_id,
                include_report,
            } => {
                assert!(!helper_mode);
                assert_eq!(job_id, "blast-job-42");
                assert!(include_report);
            }
            other => panic!("unexpected command: {other:?}"),
        }
        let cancel =
            parse_shell_line("helpers blast-cancel blast-job-42").expect("parse cancel command");
        match cancel {
            ShellCommand::ReferenceBlastAsyncCancel {
                helper_mode,
                job_id,
            } => {
                assert!(helper_mode);
                assert_eq!(job_id, "blast-job-42");
            }
            other => panic!("unexpected command: {other:?}"),
        }
        let listed = parse_shell_line("helpers blast-list").expect("parse list command");
        match listed {
            ShellCommand::ReferenceBlastAsyncList { helper_mode } => {
                assert!(helper_mode);
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_genomes_blast_with_options_json_override() {
        let cmd = parse_shell_line(
            "genomes blast ToyGenome ACGTACGT --options-json '{\"max_hits\":7,\"thresholds\":{\"min_identity_percent\":97.5}}'",
        )
        .expect("parse command");
        match cmd {
            ShellCommand::ReferenceBlast {
                helper_mode,
                genome_id,
                max_hits,
                max_hits_explicit,
                task,
                request_options_json,
                ..
            } => {
                assert!(!helper_mode);
                assert_eq!(genome_id, "ToyGenome");
                assert_eq!(max_hits, 25);
                assert!(!max_hits_explicit);
                assert!(task.is_none());
                let options = request_options_json.expect("options json");
                assert_eq!(options.get("max_hits").and_then(|v| v.as_u64()), Some(7));
                let min_ident = options
                    .get("thresholds")
                    .and_then(|v| v.get("min_identity_percent"))
                    .and_then(|v| v.as_f64());
                assert_eq!(min_ident, Some(97.5));
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_helpers_blast_with_options_file_override() {
        let td = tempdir().expect("tempdir");
        let options_path = td.path().join("blast_options.json");
        fs::write(
            &options_path,
            r#"{"max_hits":9,"thresholds":{"min_bit_score":55.0}}"#,
        )
        .expect("write options file");
        let cmd = parse_shell_line(&format!(
            "helpers blast pUC19 ACGTAG --options-file {}",
            options_path.to_string_lossy()
        ))
        .expect("parse command");
        match cmd {
            ShellCommand::ReferenceBlast {
                helper_mode,
                request_options_json,
                max_hits_explicit,
                ..
            } => {
                assert!(helper_mode);
                assert!(!max_hits_explicit);
                let options = request_options_json.expect("options json");
                assert_eq!(options.get("max_hits").and_then(|v| v.as_u64()), Some(9));
                let min_bit = options
                    .get("thresholds")
                    .and_then(|v| v.get("min_bit_score"))
                    .and_then(|v| v.as_f64());
                assert_eq!(min_bit, Some(55.0));
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
                max_hits_explicit,
                task,
                request_options_json,
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
                assert!(max_hits_explicit);
                assert_eq!(task.as_deref(), Some("blastn"));
                assert!(request_options_json.is_none());
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
                max_hits_explicit,
                task,
                request_options_json,
                track_name,
                clear_existing,
                ..
            } => {
                assert!(helper_mode);
                assert_eq!(genome_id, "pUC19");
                assert_eq!(query_sequence, "ACGTAG");
                assert_eq!(target_seq_id, "target_seq");
                assert_eq!(max_hits, 25);
                assert!(!max_hits_explicit);
                assert!(task.is_none());
                assert!(request_options_json.is_none());
                assert!(track_name.is_none());
                assert!(!clear_existing);
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_primers_design_with_backend_overrides() {
        let cmd = parse_shell_line(
            "primers design @request.json --backend primer3 --primer3-exec /opt/primer3/primer3_core",
        )
        .expect("parse command");
        match cmd {
            ShellCommand::PrimersDesign {
                request_json,
                backend,
                primer3_executable,
            } => {
                assert_eq!(request_json, "@request.json");
                assert_eq!(backend, Some(PrimerDesignBackend::Primer3));
                assert_eq!(
                    primer3_executable.as_deref(),
                    Some("/opt/primer3/primer3_core")
                );
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_primers_design_qpcr_with_backend_overrides() {
        let cmd = parse_shell_line(
            "primers design-qpcr @request.json --backend primer3 --primer3-exec /opt/primer3/primer3_core",
        )
        .expect("parse command");
        match cmd {
            ShellCommand::PrimersDesignQpcr {
                request_json,
                backend,
                primer3_executable,
            } => {
                assert_eq!(request_json, "@request.json");
                assert_eq!(backend, Some(PrimerDesignBackend::Primer3));
                assert_eq!(
                    primer3_executable.as_deref(),
                    Some("/opt/primer3/primer3_core")
                );
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_primers_preflight_with_backend_overrides() {
        let cmd = parse_shell_line(
            "primers preflight --backend primer3 --primer3-exec /opt/primer3/primer3_core",
        )
        .expect("parse command");
        match cmd {
            ShellCommand::PrimersPreflight {
                backend,
                primer3_executable,
            } => {
                assert_eq!(backend, Some(PrimerDesignBackend::Primer3));
                assert_eq!(
                    primer3_executable.as_deref(),
                    Some("/opt/primer3/primer3_core")
                );
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_primers_seed_from_feature_and_splicing() {
        let feature =
            parse_shell_line("primers seed-from-feature seq_a 7").expect("parse seed-from-feature");
        assert!(matches!(
            feature,
            ShellCommand::PrimersSeedFromFeature {
                seq_id,
                feature_id
            } if seq_id == "seq_a" && feature_id == 7
        ));

        let splicing = parse_shell_line("primers seed-from-splicing seq_a 11")
            .expect("parse seed-from-splicing");
        assert!(matches!(
            splicing,
            ShellCommand::PrimersSeedFromSplicing {
                seq_id,
                feature_id
            } if seq_id == "seq_a" && feature_id == 11
        ));
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
    fn parse_routines_explain_with_catalog() {
        let cmd = parse_shell_line(
            "routines explain golden_gate.type_iis_single_insert --catalog catalog.json",
        )
        .expect("parse routines explain");
        match cmd {
            ShellCommand::RoutinesExplain {
                catalog_path,
                routine_id,
            } => {
                assert_eq!(catalog_path.as_deref(), Some("catalog.json"));
                assert_eq!(routine_id, "golden_gate.type_iis_single_insert".to_string());
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_routines_compare_with_catalog() {
        let cmd = parse_shell_line(
            "routines compare golden_gate.type_iis_single_insert gibson.two_fragment_overlap_preview --catalog catalog.json",
        )
        .expect("parse routines compare");
        match cmd {
            ShellCommand::RoutinesCompare {
                catalog_path,
                left_routine_id,
                right_routine_id,
            } => {
                assert_eq!(catalog_path.as_deref(), Some("catalog.json"));
                assert_eq!(
                    left_routine_id,
                    "golden_gate.type_iis_single_insert".to_string()
                );
                assert_eq!(
                    right_routine_id,
                    "gibson.two_fragment_overlap_preview".to_string()
                );
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_planning_commands() {
        let profile = parse_shell_line("planning profile set @profile.json --scope global")
            .expect("parse planning profile set");
        assert!(matches!(
            profile,
            ShellCommand::PlanningProfileSet { scope, payload_json }
                if scope == PlanningProfileScope::Global && payload_json == "@profile.json"
        ));

        let objective =
            parse_shell_line("planning objective show").expect("parse planning objective show");
        assert!(matches!(objective, ShellCommand::PlanningObjectiveShow));

        let suggestions = parse_shell_line("planning suggestions list --status pending")
            .expect("parse planning suggestions list");
        assert!(matches!(
            suggestions,
            ShellCommand::PlanningSuggestionsList { status }
                if status == Some(PlanningSuggestionStatus::Pending)
        ));

        let sync = parse_shell_line(
            "planning sync pull @sync.json --source lab_manager --confidence 0.75 --snapshot-id snap_1",
        )
        .expect("parse planning sync pull");
        assert!(matches!(
            sync,
            ShellCommand::PlanningSyncPull {
                payload_json,
                source,
                confidence,
                snapshot_id
            } if payload_json == "@sync.json"
                && source.as_deref() == Some("lab_manager")
                && confidence == Some(0.75)
                && snapshot_id.as_deref() == Some("snap_1")
        ));
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
                validate_only,
            } => {
                assert_eq!(name, "clone");
                assert_eq!(bindings.get("seq_id"), Some(&"seqB".to_string()));
                assert_eq!(bindings.get("out_id"), Some(&"seqB_rev".to_string()));
                assert!(transactional);
                assert!(!validate_only);
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

        let instance_list =
            parse_shell_line("macros instance-list").expect("parse macros instance-list");
        match instance_list {
            ShellCommand::MacrosInstanceList => {}
            other => panic!("unexpected command: {other:?}"),
        }

        let instance_show =
            parse_shell_line("macros instance-show macro-1").expect("parse macros instance-show");
        match instance_show {
            ShellCommand::MacrosInstanceShow { macro_instance_id } => {
                assert_eq!(macro_instance_id, "macro-1");
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
    fn parse_macros_template_put_accepts_port_contracts() {
        let put = parse_shell_line(
            r#"macros template-put clone --script 'op {"Reverse":{"input":"${seq_id}","output_id":"${out_id}"}}' --param seq_id --param out_id=seqA_rev --input-port seq_id:sequence:one:required:Template --output-port out_id:sequence:one:optional:Derived"#,
        )
        .expect("parse macros template-put with ports");
        match put {
            ShellCommand::MacrosTemplateUpsert {
                input_ports,
                output_ports,
                ..
            } => {
                assert_eq!(input_ports.len(), 1);
                assert_eq!(input_ports[0].port_id, "seq_id");
                assert_eq!(input_ports[0].kind, "sequence");
                assert_eq!(input_ports[0].cardinality, "one");
                assert!(input_ports[0].required);
                assert_eq!(output_ports.len(), 1);
                assert_eq!(output_ports[0].port_id, "out_id");
                assert_eq!(output_ports[0].kind, "sequence");
                assert_eq!(output_ports[0].cardinality, "one");
                assert!(!output_ports[0].required);
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_macros_template_run_validate_only_flag() {
        let cmd = parse_shell_line("macros template-run clone --validate-only --bind seq_id=seqA")
            .expect("parse macros template-run validate-only");
        match cmd {
            ShellCommand::MacrosTemplateRun {
                name,
                bindings,
                transactional,
                validate_only,
            } => {
                assert_eq!(name, "clone");
                assert_eq!(bindings.get("seq_id").map(String::as_str), Some("seqA"));
                assert!(!transactional);
                assert!(validate_only);
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
    fn execute_planning_profile_merge_precedence_global_agent_project() {
        let mut engine = GentleEngine::default();

        execute_shell_command(
            &mut engine,
            &ShellCommand::PlanningProfileSet {
                scope: PlanningProfileScope::Global,
                payload_json: r#"{
                  "schema":"gentle.planning_profile.v1",
                  "capabilities":["gel_imager"],
                  "inventory":{"enzymes":{"available":false,"procurement_business_days":3}},
                  "machine_availability":{"thermocycler":{"available":true,"queue_business_days":0.5}}
                }"#
                .to_string(),
            },
        )
        .expect("set global profile");

        execute_shell_command(
            &mut engine,
            &ShellCommand::PlanningProfileSet {
                scope: PlanningProfileScope::ConfirmedAgentOverlay,
                payload_json: r#"{
                  "schema":"gentle.planning_profile.v1",
                  "capabilities":["thermocycler"],
                  "inventory":{"enzymes":{"available":true}}
                }"#
                .to_string(),
            },
        )
        .expect("set agent overlay");

        execute_shell_command(
            &mut engine,
            &ShellCommand::PlanningProfileSet {
                scope: PlanningProfileScope::ProjectOverride,
                payload_json: r#"{
                  "schema":"gentle.planning_profile.v1",
                  "inventory":{"enzymes":{"available":false,"procurement_business_days":2}}
                }"#
                .to_string(),
            },
        )
        .expect("set project override");

        let effective = engine.planning_effective_profile();
        let enzymes = effective
            .inventory
            .get("enzymes")
            .expect("enzymes inventory");
        assert!(
            !enzymes.available,
            "project override should win over agent/global"
        );
        assert_eq!(enzymes.procurement_business_days, Some(2.0));
        assert!(
            effective.capabilities.iter().any(|cap| cap == "gel_imager"),
            "global capability should remain in merged profile"
        );
        assert!(
            effective
                .capabilities
                .iter()
                .any(|cap| cap == "thermocycler"),
            "agent overlay capability should be merged"
        );
        let machine = effective
            .machine_availability
            .get("thermocycler")
            .expect("thermocycler machine entry");
        assert!(
            (machine.queue_business_days - 0.5).abs() < f64::EPSILON,
            "global machine queue assumption should survive merge"
        );
    }

    #[test]
    fn execute_planning_rejects_schema_mismatch() {
        let mut engine = GentleEngine::default();
        let err = execute_shell_command(
            &mut engine,
            &ShellCommand::PlanningProfileSet {
                scope: PlanningProfileScope::Global,
                payload_json: r#"{"schema":"gentle.planning_profile.v0"}"#.to_string(),
            },
        )
        .expect_err("profile schema mismatch must fail");
        assert!(err.contains("Unsupported planning profile schema"));

        let err = execute_shell_command(
            &mut engine,
            &ShellCommand::PlanningObjectiveSet {
                payload_json: r#"{"schema":"gentle.planning_objective.v0"}"#.to_string(),
            },
        )
        .expect_err("objective schema mismatch must fail");
        assert!(err.contains("Unsupported planning objective schema"));
    }

    #[test]
    fn execute_planning_sync_payload_validation_rejects_unknown_fields() {
        let mut engine = GentleEngine::default();
        let err = execute_shell_command(
            &mut engine,
            &ShellCommand::PlanningSyncPull {
                payload_json: r#"{"unknown_patch":{"x":1}}"#.to_string(),
                source: None,
                confidence: None,
                snapshot_id: None,
            },
        )
        .expect_err("unknown sync payload fields must fail");
        assert!(err.contains("Invalid planning sync pull JSON payload"));
    }

    #[test]
    fn execute_routines_list_without_planning_preserves_legacy_order() {
        let mut engine = GentleEngine::default();
        let tmp = tempdir().expect("tempdir");
        let catalog_path = tmp.path().join("routines_legacy_order.json");
        fs::write(
            &catalog_path,
            r#"{
  "schema": "gentle.cloning_routines.v1",
  "routines": [
    {
      "routine_id": "zeta.workflow",
      "title": "Zeta Routine",
      "family": "zeta",
      "status": "implemented",
      "vocabulary_tags": ["zeta"],
      "template_name": "zeta_template",
      "base_time_hours": 1.0,
      "required_material_classes": ["missing_item"],
      "input_ports": [{ "port_id": "seq_id", "kind": "sequence", "required": true, "cardinality": "one" }],
      "output_ports": []
    },
    {
      "routine_id": "alpha.workflow",
      "title": "Alpha Routine",
      "family": "alpha",
      "status": "implemented",
      "vocabulary_tags": ["alpha"],
      "template_name": "alpha_template",
      "base_time_hours": 10.0,
      "input_ports": [{ "port_id": "seq_id", "kind": "sequence", "required": true, "cardinality": "one" }],
      "output_ports": []
    }
  ]
}"#,
        )
        .expect("write legacy-order catalog");

        let out = execute_shell_command(
            &mut engine,
            &ShellCommand::RoutinesList {
                catalog_path: Some(catalog_path.to_string_lossy().to_string()),
                family: None,
                status: None,
                tag: None,
                query: None,
            },
        )
        .expect("list routines");
        assert_eq!(out.output["planning"]["enabled"].as_bool(), Some(false));
        let rows = out.output["routines"]
            .as_array()
            .cloned()
            .unwrap_or_default();
        assert_eq!(rows.len(), 2);
        assert_eq!(
            rows[0].get("routine_id").and_then(|v| v.as_str()),
            Some("alpha.workflow"),
            "legacy non-planning order should remain family/title based"
        );
    }

    #[test]
    fn execute_routines_list_applies_missing_material_penalty_default_business_days() {
        let mut engine = GentleEngine::default();
        let tmp = tempdir().expect("tempdir");
        let catalog_path = tmp.path().join("routines_planning_penalty.json");
        fs::write(
            &catalog_path,
            r#"{
  "schema": "gentle.cloning_routines.v1",
  "routines": [
    {
      "routine_id": "restriction.missing_mix",
      "title": "Restriction with Missing Mix",
      "family": "restriction",
      "status": "implemented",
      "vocabulary_tags": ["restriction"],
      "template_name": "restriction_missing_mix",
      "base_time_hours": 4.0,
      "base_cost": 8.0,
      "required_material_classes": ["custom_mix", "custom_mix"],
      "input_ports": [{ "port_id": "seq_id", "kind": "sequence", "required": true, "cardinality": "one" }],
      "output_ports": []
    }
  ]
}"#,
        )
        .expect("write planning catalog");

        execute_shell_command(
            &mut engine,
            &ShellCommand::PlanningProfileSet {
                scope: PlanningProfileScope::ProjectOverride,
                payload_json: r#"{
                  "schema":"gentle.planning_profile.v1",
                  "procurement_business_days_default":10,
                  "inventory":{"custom_mix":{"available":false}}
                }"#
                .to_string(),
            },
        )
        .expect("set planning profile");

        let out = execute_shell_command(
            &mut engine,
            &ShellCommand::RoutinesList {
                catalog_path: Some(catalog_path.to_string_lossy().to_string()),
                family: None,
                status: None,
                tag: None,
                query: None,
            },
        )
        .expect("list routines");
        let rows = out.output["routines"]
            .as_array()
            .cloned()
            .unwrap_or_default();
        assert_eq!(rows.len(), 1);
        let estimate = rows[0]["planning_estimate"].clone();
        let estimated_time_hours = estimate["estimated_time_hours"]
            .as_f64()
            .unwrap_or(f64::NAN);
        assert!(
            (estimated_time_hours - 340.0).abs() < 1e-9,
            "expected base 4h + 10 business days (weekend-aware elapsed-time) penalty"
        );
        let procurement_days = estimate["explanation"]["procurement_delay_business_days"]
            .as_f64()
            .unwrap_or(f64::NAN);
        assert!(
            (procurement_days - 10.0).abs() < 1e-9,
            "expected exactly one deduplicated missing material penalty"
        );
        let missing = estimate["explanation"]["missing_material_classes"]
            .as_array()
            .cloned()
            .unwrap_or_default();
        assert_eq!(missing.len(), 1);
        assert_eq!(missing[0].as_str(), Some("custom_mix"));
    }

    #[test]
    fn execute_planning_suggestion_lifecycle_pending_accept_reject() {
        let mut engine = GentleEngine::default();

        let pull = execute_shell_command(
            &mut engine,
            &ShellCommand::PlanningSyncPull {
                payload_json: r#"{
                  "profile_patch":{"schema":"gentle.planning_profile.v1","capabilities":["pcr_machine"]},
                  "message":"sync pull"
                }"#
                .to_string(),
                source: Some("lab_agent".to_string()),
                confidence: Some(0.8),
                snapshot_id: Some("snap_pull_1".to_string()),
            },
        )
        .expect("sync pull suggestion");
        let suggestion_id = pull.output["suggestion"]["suggestion_id"]
            .as_str()
            .expect("pull suggestion id")
            .to_string();

        let pending = execute_shell_command(
            &mut engine,
            &ShellCommand::PlanningSuggestionsList {
                status: Some(PlanningSuggestionStatus::Pending),
            },
        )
        .expect("list pending suggestions");
        assert_eq!(pending.output["suggestion_count"].as_u64(), Some(1));

        let accepted = execute_shell_command(
            &mut engine,
            &ShellCommand::PlanningSuggestionAccept {
                suggestion_id: suggestion_id.clone(),
            },
        )
        .expect("accept suggestion");
        assert_eq!(accepted.output["status"].as_str(), Some("accepted"));
        assert!(
            engine
                .planning_effective_profile()
                .capabilities
                .iter()
                .any(|cap| cap == "pcr_machine"),
            "accepted suggestion should update confirmed overlay profile"
        );

        let push = execute_shell_command(
            &mut engine,
            &ShellCommand::PlanningSyncPush {
                payload_json: r#"{
                  "objective_patch":{"schema":"gentle.planning_objective.v1","weight_time":2.5},
                  "message":"sync push"
                }"#
                .to_string(),
                source: Some("lab_agent".to_string()),
                confidence: Some(0.6),
                snapshot_id: Some("snap_push_1".to_string()),
            },
        )
        .expect("sync push suggestion");
        let reject_id = push.output["suggestion"]["suggestion_id"]
            .as_str()
            .expect("push suggestion id")
            .to_string();

        let rejected = execute_shell_command(
            &mut engine,
            &ShellCommand::PlanningSuggestionReject {
                suggestion_id: reject_id.clone(),
                reason: Some("manual_override".to_string()),
            },
        )
        .expect("reject suggestion");
        assert_eq!(rejected.output["status"].as_str(), Some("rejected"));
        assert_eq!(
            rejected.output["suggestion"]["rejection_reason"].as_str(),
            Some("manual_override")
        );

        let accepted_rows = execute_shell_command(
            &mut engine,
            &ShellCommand::PlanningSuggestionsList {
                status: Some(PlanningSuggestionStatus::Accepted),
            },
        )
        .expect("list accepted suggestions");
        assert_eq!(accepted_rows.output["suggestion_count"].as_u64(), Some(1));
        let rejected_rows = execute_shell_command(
            &mut engine,
            &ShellCommand::PlanningSuggestionsList {
                status: Some(PlanningSuggestionStatus::Rejected),
            },
        )
        .expect("list rejected suggestions");
        assert_eq!(rejected_rows.output["suggestion_count"].as_u64(), Some(1));
    }

    #[test]
    fn execute_routines_explain_with_explicit_alternative_metadata() {
        let mut engine = GentleEngine::default();
        let tmp = tempdir().expect("tempdir");
        let catalog_path = tmp.path().join("cloning_routines_explain.json");
        fs::write(
            &catalog_path,
            r#"{
  "schema": "gentle.cloning_routines.v1",
  "routines": [
    {
      "routine_id": "golden_gate.type_iis_single_insert",
      "title": "Golden Gate Type IIS Single Insert",
      "family": "golden_gate",
      "status": "implemented",
      "vocabulary_tags": ["golden_gate", "type_iis", "ligation"],
      "summary": "Type IIS one-pot assembly for one insert.",
      "purpose": "Assemble one insert into one destination vector with directional overhangs.",
      "mechanism": "Digest + ligate cycles with Type IIS enzymes and programmed overhang junctions.",
      "requires": ["Type IIS-compatible enzymes", "Defined non-conflicting overhang plan"],
      "contraindications": ["Ambiguous or duplicated internal overhangs"],
      "confusing_alternatives": ["gibson.two_fragment_overlap_preview"],
      "difference_matrix": [
        { "axis": "junction_constraint", "value": "explicit overhang grammar" },
        { "axis": "assembly_mode", "value": "restriction-ligation cycling" }
      ],
      "disambiguation_questions": ["Do you require Type IIS junction tokens?"],
      "failure_modes": ["duplicate junction overhang token"],
      "template_name": "golden_gate_single_insert",
      "input_ports": [
        { "port_id": "seq_id", "kind": "sequence", "required": true, "cardinality": "one" }
      ],
      "output_ports": [
        { "port_id": "output_id", "kind": "sequence", "required": false, "cardinality": "one" }
      ]
    },
    {
      "routine_id": "gibson.two_fragment_overlap_preview",
      "title": "Gibson Two-Fragment Overlap Preview",
      "family": "gibson",
      "status": "implemented",
      "vocabulary_tags": ["gibson", "overlap", "assembly"],
      "summary": "Overlap assembly preflight and preview path.",
      "template_name": "gibson_two_fragment_overlap_preview",
      "input_ports": [
        { "port_id": "left_seq_id", "kind": "sequence", "required": true, "cardinality": "one" },
        { "port_id": "right_seq_id", "kind": "sequence", "required": true, "cardinality": "one" }
      ],
      "output_ports": [
        { "port_id": "output_id", "kind": "sequence", "required": false, "cardinality": "one" }
      ]
    }
  ]
}"#,
        )
        .expect("write routines catalog");

        let run = execute_shell_command(
            &mut engine,
            &ShellCommand::RoutinesExplain {
                catalog_path: Some(catalog_path.to_string_lossy().to_string()),
                routine_id: "golden_gate.type_iis_single_insert".to_string(),
            },
        )
        .expect("routines explain");
        assert!(!run.state_changed);
        assert_eq!(
            run.output["schema"].as_str(),
            Some(CLONING_ROUTINE_EXPLAIN_SCHEMA)
        );
        assert_eq!(
            run.output["routine"]["routine_id"].as_str(),
            Some("golden_gate.type_iis_single_insert")
        );
        let alternatives = run
            .output
            .get("alternatives")
            .and_then(|value| value.as_array())
            .cloned()
            .unwrap_or_default();
        assert_eq!(alternatives.len(), 1);
        assert_eq!(
            alternatives[0]
                .get("routine_id")
                .and_then(|value| value.as_str()),
            Some("gibson.two_fragment_overlap_preview")
        );
        let explanation_requires = run
            .output
            .get("explanation")
            .and_then(|value| value.get("requires"))
            .and_then(|value| value.as_array())
            .cloned()
            .unwrap_or_default();
        assert_eq!(explanation_requires.len(), 2);
    }

    #[test]
    fn execute_routines_compare_returns_difference_matrix() {
        let mut engine = GentleEngine::default();
        let tmp = tempdir().expect("tempdir");
        let catalog_path = tmp.path().join("cloning_routines_compare.json");
        fs::write(
            &catalog_path,
            r#"{
  "schema": "gentle.cloning_routines.v1",
  "routines": [
    {
      "routine_id": "golden_gate.type_iis_single_insert",
      "title": "Golden Gate Type IIS Single Insert",
      "family": "golden_gate",
      "status": "implemented",
      "vocabulary_tags": ["golden_gate", "assembly"],
      "summary": "Type IIS one-pot assembly for one insert.",
      "difference_matrix": [
        { "axis": "assembly_mode", "value": "restriction-ligation cycling" },
        { "axis": "junction_constraint", "value": "explicit overhang grammar" }
      ],
      "confusing_alternatives": ["gibson.two_fragment_overlap_preview"],
      "template_name": "golden_gate_single_insert",
      "input_ports": [
        { "port_id": "seq_id", "kind": "sequence", "required": true, "cardinality": "one" }
      ],
      "output_ports": [
        { "port_id": "output_id", "kind": "sequence", "required": false, "cardinality": "one" }
      ]
    },
    {
      "routine_id": "gibson.two_fragment_overlap_preview",
      "title": "Gibson Two-Fragment Overlap Preview",
      "family": "gibson",
      "status": "implemented",
      "vocabulary_tags": ["gibson", "assembly"],
      "summary": "Overlap assembly preflight and preview path.",
      "difference_matrix": [
        { "axis": "assembly_mode", "value": "homology-overlap assembly" },
        { "axis": "junction_constraint", "value": "homology overlap sequence identity" }
      ],
      "confusing_alternatives": ["golden_gate.type_iis_single_insert"],
      "template_name": "gibson_two_fragment_overlap_preview",
      "input_ports": [
        { "port_id": "left_seq_id", "kind": "sequence", "required": true, "cardinality": "one" },
        { "port_id": "right_seq_id", "kind": "sequence", "required": true, "cardinality": "one" }
      ],
      "output_ports": [
        { "port_id": "output_id", "kind": "sequence", "required": false, "cardinality": "one" }
      ]
    }
  ]
}"#,
        )
        .expect("write routines catalog");

        let run = execute_shell_command(
            &mut engine,
            &ShellCommand::RoutinesCompare {
                catalog_path: Some(catalog_path.to_string_lossy().to_string()),
                left_routine_id: "golden_gate.type_iis_single_insert".to_string(),
                right_routine_id: "gibson.two_fragment_overlap_preview".to_string(),
            },
        )
        .expect("routines compare");
        assert!(!run.state_changed);
        assert_eq!(
            run.output["schema"].as_str(),
            Some(CLONING_ROUTINE_COMPARE_SCHEMA)
        );
        assert_eq!(
            run.output["comparison"]["cross_referenced_as_alternatives"].as_bool(),
            Some(true)
        );
        let matrix = run
            .output
            .get("comparison")
            .and_then(|value| value.get("difference_matrix"))
            .and_then(|value| value.as_array())
            .cloned()
            .unwrap_or_default();
        assert!(
            matrix.iter().any(
                |row| row.get("axis").and_then(|value| value.as_str()) == Some("assembly_mode")
            )
        );
    }

    #[test]
    fn load_cloning_routine_catalog_rejects_unknown_confusing_alternative() {
        let tmp = tempdir().expect("tempdir");
        let catalog_path = tmp.path().join("cloning_routines_bad_alt.json");
        fs::write(
            &catalog_path,
            r#"{
  "schema": "gentle.cloning_routines.v1",
  "routines": [
    {
      "routine_id": "golden_gate.type_iis_single_insert",
      "title": "Golden Gate Type IIS Single Insert",
      "family": "golden_gate",
      "status": "implemented",
      "vocabulary_tags": ["golden_gate"],
      "summary": "baseline",
      "confusing_alternatives": ["missing.routine.id"],
      "template_name": "golden_gate_single_insert",
      "input_ports": [
        { "port_id": "seq_id", "kind": "sequence", "required": true, "cardinality": "one" }
      ],
      "output_ports": [
        { "port_id": "output_id", "kind": "sequence", "required": false, "cardinality": "one" }
      ]
    }
  ]
}"#,
        )
        .expect("write routines catalog");

        let err = load_cloning_routine_catalog(catalog_path.to_string_lossy().as_ref())
            .expect_err("unknown alternatives must be rejected");
        assert!(err.contains("unknown confusing_alternative"));
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
                input_ports: vec![],
                output_ports: vec![],
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
                validate_only: false,
            },
        )
        .expect("run macros template");
        assert!(run.state_changed);
        assert_eq!(run.output["template_name"].as_str(), Some("clone"));
        assert!(engine.state().sequences.contains_key("seqA_rev"));
    }

    #[test]
    fn execute_macros_template_run_records_macro_instance() {
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
                input_ports: vec![WorkflowMacroTemplatePort {
                    port_id: "seq_id".to_string(),
                    kind: "sequence".to_string(),
                    required: true,
                    cardinality: "one".to_string(),
                    description: Some("Template input".to_string()),
                }],
                output_ports: vec![WorkflowMacroTemplatePort {
                    port_id: "out_id".to_string(),
                    kind: "sequence".to_string(),
                    required: false,
                    cardinality: "one".to_string(),
                    description: Some("Derived sequence id".to_string()),
                }],
                script: r#"op {"Reverse":{"input":"${seq_id}","output_id":"${out_id}"}}"#
                    .to_string(),
            },
        )
        .expect("upsert template");

        let run = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateRun {
                name: "clone".to_string(),
                bindings: HashMap::from([
                    ("seq_id".to_string(), "seqA".to_string()),
                    ("out_id".to_string(), "seqA_rev".to_string()),
                ]),
                transactional: false,
                validate_only: false,
            },
        )
        .expect("run template");
        assert!(run.state_changed);
        assert!(run.output["macro_instance_id"].as_str().is_some());
        assert_eq!(engine.state().lineage.macro_instances.len(), 1);
        let instance = &engine.state().lineage.macro_instances[0];
        assert_eq!(instance.template_name.as_deref(), Some("clone"));
        assert!(!instance.expanded_op_ids.is_empty());
        assert_eq!(instance.bound_inputs.len(), 1);
        assert_eq!(instance.bound_inputs[0].port_id, "seq_id");
        assert_eq!(instance.bound_outputs.len(), 1);
        assert_eq!(instance.bound_outputs[0].port_id, "out_id");
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
                validate_only: false,
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
                validate_only: false,
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
                validate_only: false,
            },
        )
        .expect("run imported branch template");
        assert!(run.state_changed);
        assert!(engine.state().sequences.contains_key("seqA_branch"));
        assert!(engine.state().sequences.contains_key("seqA_branch_rc"));
    }

    #[test]
    fn execute_macros_template_run_validate_only_reports_preflight_without_mutation() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "seqA".to_string(),
            DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
        );
        let mut engine = GentleEngine::from_state(state);

        let path = format!(
            "{}/assets/cloning_patterns_catalog/sequence/transform/branch_reverse_complement.json",
            env!("CARGO_MANIFEST_DIR")
        );
        execute_shell_command(&mut engine, &ShellCommand::MacrosTemplateImport { path })
            .expect("import branch template");
        assert!(!engine.state().sequences.contains_key("seqA_branch"));

        let run = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateRun {
                name: "branch_reverse_complement".to_string(),
                bindings: HashMap::new(),
                transactional: false,
                validate_only: true,
            },
        )
        .expect("validate-only run");
        assert!(!run.state_changed);
        assert_eq!(run.output["validate_only"].as_bool(), Some(true));
        assert_eq!(run.output["can_execute"].as_bool(), Some(false));
        assert!(
            run.output
                .get("preflight")
                .and_then(|v| v.get("errors"))
                .and_then(|v| v.as_array())
                .map(|rows| !rows.is_empty())
                .unwrap_or(false)
        );
        assert!(!engine.state().sequences.contains_key("seqA_branch"));
    }

    #[test]
    fn execute_macros_template_run_preflight_failure_records_macro_instance() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "seqA".to_string(),
            DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
        );
        let mut engine = GentleEngine::from_state(state);

        let path = format!(
            "{}/assets/cloning_patterns_catalog/sequence/transform/branch_reverse_complement.json",
            env!("CARGO_MANIFEST_DIR")
        );
        execute_shell_command(&mut engine, &ShellCommand::MacrosTemplateImport { path })
            .expect("import branch template");
        let err = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateRun {
                name: "branch_reverse_complement".to_string(),
                bindings: HashMap::new(),
                transactional: false,
                validate_only: false,
            },
        )
        .expect_err("preflight should fail");
        assert!(err.contains("macro_instance_id="));
        assert_eq!(engine.state().lineage.macro_instances.len(), 1);
        let instance = &engine.state().lineage.macro_instances[0];
        assert_eq!(instance.status, MacroInstanceStatus::Failed);
        assert_eq!(
            instance.template_name.as_deref(),
            Some("branch_reverse_complement")
        );
        assert!(instance.status_message.is_some());
    }

    #[test]
    fn execute_macros_template_validate_only_checks_anchor_semantics() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "seqA".to_string(),
            DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
        );
        let mut engine = GentleEngine::from_state(state);
        execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateUpsert {
                name: "anchor_semantics".to_string(),
                description: Some("anchor semantics test".to_string()),
                details_url: None,
                parameters: vec![
                    WorkflowMacroTemplateParam {
                        name: "seq_id".to_string(),
                        default_value: None,
                        required: true,
                    },
                    WorkflowMacroTemplateParam {
                        name: "anchor_a".to_string(),
                        default_value: None,
                        required: true,
                    },
                ],
                input_ports: vec![
                    WorkflowMacroTemplatePort {
                        port_id: "seq_id".to_string(),
                        kind: "sequence".to_string(),
                        required: true,
                        cardinality: "one".to_string(),
                        description: Some("sequence input".to_string()),
                    },
                    WorkflowMacroTemplatePort {
                        port_id: "anchor_a".to_string(),
                        kind: "sequence_anchor".to_string(),
                        required: true,
                        cardinality: "one".to_string(),
                        description: Some("anchor input".to_string()),
                    },
                ],
                output_ports: vec![],
                script: r#"op {"Reverse":{"input":"${seq_id}"}}"#.to_string(),
            },
        )
        .expect("upsert template");

        let run = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateRun {
                name: "anchor_semantics".to_string(),
                bindings: HashMap::from([
                    ("seq_id".to_string(), "seqA".to_string()),
                    ("anchor_a".to_string(), "999".to_string()),
                ]),
                transactional: false,
                validate_only: true,
            },
        )
        .expect("validate-only");
        assert!(!run.state_changed);
        assert_eq!(run.output["can_execute"].as_bool(), Some(false));
        let preflight_errors = run
            .output
            .get("preflight")
            .and_then(|preflight| preflight.get("errors"))
            .and_then(|rows| rows.as_array())
            .cloned()
            .unwrap_or_default()
            .into_iter()
            .filter_map(|row| row.as_str().map(|s| s.to_string()))
            .collect::<Vec<_>>();
        assert!(
            preflight_errors
                .iter()
                .any(|message| message.contains("anchor check failed")),
            "expected semantic anchor validation error, got: {:?}",
            preflight_errors
        );
    }

    #[test]
    fn execute_macros_template_validate_only_rejects_output_alias_collisions() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "seqA".to_string(),
            DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
        );
        let mut engine = GentleEngine::from_state(state);
        execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateUpsert {
                name: "alias_collision".to_string(),
                description: Some("output alias collision".to_string()),
                details_url: None,
                parameters: vec![
                    WorkflowMacroTemplateParam {
                        name: "seq_id".to_string(),
                        default_value: None,
                        required: true,
                    },
                    WorkflowMacroTemplateParam {
                        name: "out_a".to_string(),
                        default_value: None,
                        required: true,
                    },
                    WorkflowMacroTemplateParam {
                        name: "out_b".to_string(),
                        default_value: None,
                        required: true,
                    },
                ],
                input_ports: vec![WorkflowMacroTemplatePort {
                    port_id: "seq_id".to_string(),
                    kind: "sequence".to_string(),
                    required: true,
                    cardinality: "one".to_string(),
                    description: None,
                }],
                output_ports: vec![
                    WorkflowMacroTemplatePort {
                        port_id: "out_a".to_string(),
                        kind: "sequence".to_string(),
                        required: true,
                        cardinality: "one".to_string(),
                        description: None,
                    },
                    WorkflowMacroTemplatePort {
                        port_id: "out_b".to_string(),
                        kind: "sequence".to_string(),
                        required: true,
                        cardinality: "one".to_string(),
                        description: None,
                    },
                ],
                script: r#"op {"Reverse":{"input":"${seq_id}","output_id":"${out_a}"}}"#
                    .to_string(),
            },
        )
        .expect("upsert template");

        let run = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateRun {
                name: "alias_collision".to_string(),
                bindings: HashMap::from([
                    ("seq_id".to_string(), "seqA".to_string()),
                    ("out_a".to_string(), "same_out".to_string()),
                    ("out_b".to_string(), "same_out".to_string()),
                ]),
                transactional: false,
                validate_only: true,
            },
        )
        .expect("validate-only");
        assert_eq!(run.output["can_execute"].as_bool(), Some(false));
        let preflight_errors = run
            .output
            .get("preflight")
            .and_then(|preflight| preflight.get("errors"))
            .and_then(|rows| rows.as_array())
            .cloned()
            .unwrap_or_default()
            .into_iter()
            .filter_map(|row| row.as_str().map(|s| s.to_string()))
            .collect::<Vec<_>>();
        assert!(
            preflight_errors
                .iter()
                .any(|message| message.contains("Output alias conflict")),
            "expected output alias conflict error, got: {:?}",
            preflight_errors
        );
    }

    #[test]
    fn execute_macros_template_validate_only_applies_gibson_family_overlap_checks() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "left".to_string(),
            DNAsequence::from_sequence("ATGCCGTTACCGGTTAAACCCGGGTTT").expect("sequence"),
        );
        state.sequences.insert(
            "right_ok".to_string(),
            DNAsequence::from_sequence("AACCCGGGTTTGGGAAATTTCCCGGG").expect("sequence"),
        );
        state.sequences.insert(
            "right_bad".to_string(),
            DNAsequence::from_sequence("TTTAAACCCGGGGGGAAATTTCCCGGG").expect("sequence"),
        );
        let mut engine = GentleEngine::from_state(state);

        let path = format!(
            "{}/assets/cloning_patterns_catalog/gibson/overlap_assembly/gibson_two_fragment_overlap_preview.json",
            env!("CARGO_MANIFEST_DIR")
        );
        execute_shell_command(&mut engine, &ShellCommand::MacrosTemplateImport { path })
            .expect("import gibson template");

        let ok_run = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateRun {
                name: "gibson_two_fragment_overlap_preview".to_string(),
                bindings: HashMap::from([
                    ("left_seq_id".to_string(), "left".to_string()),
                    ("right_seq_id".to_string(), "right_ok".to_string()),
                    ("overlap_bp".to_string(), "11".to_string()),
                    ("assembly_prefix".to_string(), "gib_ok".to_string()),
                    ("output_id".to_string(), "gib_ok_forward".to_string()),
                ]),
                transactional: false,
                validate_only: true,
            },
        )
        .expect("gibson validate-only success");
        assert_eq!(ok_run.output["can_execute"].as_bool(), Some(true));
        let ok_errors = ok_run
            .output
            .get("preflight")
            .and_then(|preflight| preflight.get("errors"))
            .and_then(|rows| rows.as_array())
            .cloned()
            .unwrap_or_default();
        assert!(
            ok_errors.is_empty(),
            "did not expect Gibson preflight errors for matching overlap: {:?}",
            ok_errors
        );

        let bad_run = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateRun {
                name: "gibson_two_fragment_overlap_preview".to_string(),
                bindings: HashMap::from([
                    ("left_seq_id".to_string(), "left".to_string()),
                    ("right_seq_id".to_string(), "right_bad".to_string()),
                    ("overlap_bp".to_string(), "11".to_string()),
                ]),
                transactional: false,
                validate_only: true,
            },
        )
        .expect("gibson validate-only mismatch");
        assert_eq!(bad_run.output["can_execute"].as_bool(), Some(false));
        let bad_errors = bad_run
            .output
            .get("preflight")
            .and_then(|preflight| preflight.get("errors"))
            .and_then(|rows| rows.as_array())
            .cloned()
            .unwrap_or_default()
            .into_iter()
            .filter_map(|row| row.as_str().map(|s| s.to_string()))
            .collect::<Vec<_>>();
        assert!(
            bad_errors
                .iter()
                .any(|message| message.contains("Gibson overlap mismatch")),
            "expected Gibson overlap mismatch error, got: {:?}",
            bad_errors
        );
    }

    #[test]
    fn execute_macros_template_validate_only_applies_restriction_family_checks() {
        let enzymes = crate::enzymes::active_restriction_enzymes();
        assert!(
            enzymes.len() >= 2,
            "restriction enzyme catalog should provide at least two enzymes"
        );
        let enzyme_a = enzymes
            .iter()
            .find(|enzyme| !enzyme.sequence.is_empty())
            .expect("enzyme_a with sequence");
        let enzyme_b = enzymes
            .iter()
            .find(|enzyme| enzyme.name != enzyme_a.name && !enzyme.sequence.is_empty())
            .expect("enzyme_b with sequence");
        let sequence_text = format!("AAA{}TTT{}GGG", enzyme_a.sequence, enzyme_b.sequence);

        let mut state = ProjectState::default();
        state.sequences.insert(
            "vector".to_string(),
            DNAsequence::from_sequence(&sequence_text).expect("sequence"),
        );
        let mut engine = GentleEngine::from_state(state);

        let path = format!(
            "{}/assets/cloning_patterns_catalog/restriction/digest_ligation/digest_ligate_extract_sticky.json",
            env!("CARGO_MANIFEST_DIR")
        );
        execute_shell_command(&mut engine, &ShellCommand::MacrosTemplateImport { path })
            .expect("import restriction template");

        let ok_run = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateRun {
                name: "digest_ligate_extract_sticky".to_string(),
                bindings: HashMap::from([
                    ("seq_id".to_string(), "vector".to_string()),
                    ("enzyme_a".to_string(), enzyme_a.name.clone()),
                    ("enzyme_b".to_string(), enzyme_b.name.clone()),
                    ("left_fragment".to_string(), "1".to_string()),
                    ("right_fragment".to_string(), "2".to_string()),
                    ("extract_from".to_string(), "0".to_string()),
                    ("extract_to".to_string(), "20".to_string()),
                    ("digest_prefix".to_string(), "digest_ok".to_string()),
                    ("ligation_prefix".to_string(), "lig_ok".to_string()),
                    ("output_id".to_string(), "restriction_ok".to_string()),
                ]),
                transactional: false,
                validate_only: true,
            },
        )
        .expect("restriction validate-only success");
        assert_eq!(ok_run.output["can_execute"].as_bool(), Some(true));
        let ok_errors = ok_run
            .output
            .get("preflight")
            .and_then(|preflight| preflight.get("errors"))
            .and_then(|rows| rows.as_array())
            .cloned()
            .unwrap_or_default();
        assert!(
            ok_errors.is_empty(),
            "did not expect restriction preflight errors for valid inputs: {:?}",
            ok_errors
        );

        let duplicate_run = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateRun {
                name: "digest_ligate_extract_sticky".to_string(),
                bindings: HashMap::from([
                    ("seq_id".to_string(), "vector".to_string()),
                    ("enzyme_a".to_string(), enzyme_a.name.clone()),
                    ("enzyme_b".to_string(), enzyme_a.name.clone()),
                ]),
                transactional: false,
                validate_only: true,
            },
        )
        .expect("restriction validate-only duplicate enzyme");
        assert_eq!(duplicate_run.output["can_execute"].as_bool(), Some(false));
        let duplicate_errors = duplicate_run
            .output
            .get("preflight")
            .and_then(|preflight| preflight.get("errors"))
            .and_then(|rows| rows.as_array())
            .cloned()
            .unwrap_or_default()
            .into_iter()
            .filter_map(|row| row.as_str().map(|s| s.to_string()))
            .collect::<Vec<_>>();
        assert!(
            duplicate_errors
                .iter()
                .any(|message| message.contains("expects distinct enzyme inputs")),
            "expected duplicate-enzyme restriction error, got: {:?}",
            duplicate_errors
        );

        let unknown_run = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateRun {
                name: "digest_ligate_extract_sticky".to_string(),
                bindings: HashMap::from([
                    ("seq_id".to_string(), "vector".to_string()),
                    ("enzyme_a".to_string(), "__missing_enzyme__".to_string()),
                    ("enzyme_b".to_string(), enzyme_b.name.clone()),
                ]),
                transactional: false,
                validate_only: true,
            },
        )
        .expect("restriction validate-only unknown enzyme");
        assert_eq!(unknown_run.output["can_execute"].as_bool(), Some(false));
        let unknown_errors = unknown_run
            .output
            .get("preflight")
            .and_then(|preflight| preflight.get("errors"))
            .and_then(|rows| rows.as_array())
            .cloned()
            .unwrap_or_default()
            .into_iter()
            .filter_map(|row| row.as_str().map(|s| s.to_string()))
            .collect::<Vec<_>>();
        assert!(
            unknown_errors
                .iter()
                .any(|message| message.contains("unknown enzyme name")),
            "expected unknown-enzyme restriction error, got: {:?}",
            unknown_errors
        );
    }

    #[test]
    fn execute_macros_template_validate_only_applies_golden_gate_family_checks() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "gg_vector".to_string(),
            DNAsequence::from_sequence("AAAAGGTCTCTTTTGGTCTCAAAA").expect("sequence"),
        );
        state.sequences.insert(
            "gg_insert".to_string(),
            DNAsequence::from_sequence("CCCCGGTCTCGGGG").expect("sequence"),
        );
        state.sequences.insert(
            "plain_seq".to_string(),
            DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
        );
        let mut engine = GentleEngine::from_state(state);

        let base = env!("CARGO_MANIFEST_DIR");
        let single_path = format!(
            "{}/assets/cloning_patterns_catalog/golden_gate/type_iis/golden_gate_single_insert.json",
            base
        );
        let multi_path = format!(
            "{}/assets/cloning_patterns_catalog/golden_gate/type_iis/golden_gate_multi_insert.json",
            base
        );
        execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateImport { path: single_path },
        )
        .expect("import golden gate single");
        execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateImport { path: multi_path },
        )
        .expect("import golden gate multi");

        let ok_run = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateRun {
                name: "golden_gate_single_insert".to_string(),
                bindings: HashMap::from([
                    ("vector_seq_id".to_string(), "gg_vector".to_string()),
                    ("insert_seq_id".to_string(), "gg_insert".to_string()),
                    ("type_iis_enzyme".to_string(), "Eco31".to_string()),
                    ("junction_overhangs".to_string(), "AATG".to_string()),
                    ("vector_fragment".to_string(), "1".to_string()),
                    ("insert_fragment".to_string(), "1".to_string()),
                ]),
                transactional: false,
                validate_only: true,
            },
        )
        .expect("golden gate validate-only success");
        assert_eq!(ok_run.output["can_execute"].as_bool(), Some(true));

        let unknown_or_non_type_iis = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateRun {
                name: "golden_gate_single_insert".to_string(),
                bindings: HashMap::from([
                    ("vector_seq_id".to_string(), "gg_vector".to_string()),
                    ("insert_seq_id".to_string(), "gg_insert".to_string()),
                    ("type_iis_enzyme".to_string(), "EcoRI".to_string()),
                ]),
                transactional: false,
                validate_only: true,
            },
        )
        .expect("golden gate validate-only non-type-iis");
        assert_eq!(
            unknown_or_non_type_iis.output["can_execute"].as_bool(),
            Some(false)
        );
        let non_type_iis_errors = unknown_or_non_type_iis
            .output
            .get("preflight")
            .and_then(|preflight| preflight.get("errors"))
            .and_then(|rows| rows.as_array())
            .cloned()
            .unwrap_or_default()
            .into_iter()
            .filter_map(|row| row.as_str().map(|s| s.to_string()))
            .collect::<Vec<_>>();
        assert!(
            non_type_iis_errors
                .iter()
                .any(|message| message.contains("Type IIS-capable")),
            "expected Type IIS-capable error, got: {:?}",
            non_type_iis_errors
        );

        let missing_site = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateRun {
                name: "golden_gate_single_insert".to_string(),
                bindings: HashMap::from([
                    ("vector_seq_id".to_string(), "plain_seq".to_string()),
                    ("insert_seq_id".to_string(), "gg_insert".to_string()),
                    ("type_iis_enzyme".to_string(), "Eco31".to_string()),
                ]),
                transactional: false,
                validate_only: true,
            },
        )
        .expect("golden gate validate-only missing site");
        assert_eq!(missing_site.output["can_execute"].as_bool(), Some(false));
        let missing_site_errors = missing_site
            .output
            .get("preflight")
            .and_then(|preflight| preflight.get("errors"))
            .and_then(|rows| rows.as_array())
            .cloned()
            .unwrap_or_default()
            .into_iter()
            .filter_map(|row| row.as_str().map(|s| s.to_string()))
            .collect::<Vec<_>>();
        assert!(
            missing_site_errors
                .iter()
                .any(|message| message.contains("no recognition site")),
            "expected missing-site error, got: {:?}",
            missing_site_errors
        );

        let invalid_junction_or_fragment = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateRun {
                name: "golden_gate_multi_insert".to_string(),
                bindings: HashMap::from([
                    ("vector_seq_id".to_string(), "gg_vector".to_string()),
                    ("insert_a_seq_id".to_string(), "gg_insert".to_string()),
                    ("insert_b_seq_id".to_string(), "gg_insert".to_string()),
                    ("type_iis_enzyme".to_string(), "Eco31".to_string()),
                    ("junction_overhangs".to_string(), "AATG,,GCTT".to_string()),
                    ("vector_fragment".to_string(), "0".to_string()),
                ]),
                transactional: false,
                validate_only: true,
            },
        )
        .expect("golden gate validate-only invalid overhang/fragment");
        assert_eq!(
            invalid_junction_or_fragment.output["can_execute"].as_bool(),
            Some(false)
        );
        let invalid_errors = invalid_junction_or_fragment
            .output
            .get("preflight")
            .and_then(|preflight| preflight.get("errors"))
            .and_then(|rows| rows.as_array())
            .cloned()
            .unwrap_or_default()
            .into_iter()
            .filter_map(|row| row.as_str().map(|s| s.to_string()))
            .collect::<Vec<_>>();
        assert!(
            invalid_errors.iter().any(|message| {
                message.contains("contains an empty token")
                    || message.contains("expects vector_fragment >= 1")
            }),
            "expected invalid-junction/fragment errors, got: {:?}",
            invalid_errors
        );
    }

    #[test]
    fn execute_macros_template_validate_only_applies_gateway_topo_and_ta_gc_checks() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "gw_donor".to_string(),
            DNAsequence::from_sequence("ACGTACGTACGA").expect("sequence"),
        );
        state.sequences.insert(
            "gw_entry".to_string(),
            DNAsequence::from_sequence("TTTTACGTACGT").expect("sequence"),
        );
        state.sequences.insert(
            "topo_vector_t".to_string(),
            DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
        );
        state.sequences.insert(
            "topo_insert_a".to_string(),
            DNAsequence::from_sequence("ACGTACGTACGA").expect("sequence"),
        );
        state.sequences.insert(
            "topo_insert_bad".to_string(),
            DNAsequence::from_sequence("ACGTACGTACGG").expect("sequence"),
        );
        state.sequences.insert(
            "topo_insert_cacc".to_string(),
            DNAsequence::from_sequence("CACCACGTACGTACGA").expect("sequence"),
        );
        state.sequences.insert(
            "gc_vector_c".to_string(),
            DNAsequence::from_sequence("ACGTACGTACGC").expect("sequence"),
        );
        state.sequences.insert(
            "gc_insert_g".to_string(),
            DNAsequence::from_sequence("ACGTACGTACGG").expect("sequence"),
        );
        let mut engine = GentleEngine::from_state(state);

        let base = env!("CARGO_MANIFEST_DIR");
        for path in [
            format!(
                "{}/assets/cloning_patterns_catalog/gateway/bp_lr/gateway_bp_single_insert.json",
                base
            ),
            format!(
                "{}/assets/cloning_patterns_catalog/topo/entry/topo_ta_single_insert.json",
                base
            ),
            format!(
                "{}/assets/cloning_patterns_catalog/topo/entry/topo_directional_cacc_single_insert.json",
                base
            ),
            format!(
                "{}/assets/cloning_patterns_catalog/ta_gc/entry/ta_clone_single_insert.json",
                base
            ),
            format!(
                "{}/assets/cloning_patterns_catalog/ta_gc/entry/gc_clone_single_insert.json",
                base
            ),
        ] {
            execute_shell_command(&mut engine, &ShellCommand::MacrosTemplateImport { path })
                .expect("import family template");
        }

        let gateway_ok = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateRun {
                name: "gateway_bp_single_insert".to_string(),
                bindings: HashMap::from([
                    ("donor_seq_id".to_string(), "gw_donor".to_string()),
                    ("entry_vector_seq_id".to_string(), "gw_entry".to_string()),
                    ("gateway_phase".to_string(), "bp".to_string()),
                    ("att_tokens".to_string(), "attB,attP".to_string()),
                ]),
                transactional: false,
                validate_only: true,
            },
        )
        .expect("gateway bp validate-only ok");
        assert_eq!(gateway_ok.output["can_execute"].as_bool(), Some(true));

        let gateway_bad = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateRun {
                name: "gateway_bp_single_insert".to_string(),
                bindings: HashMap::from([
                    ("donor_seq_id".to_string(), "gw_donor".to_string()),
                    ("entry_vector_seq_id".to_string(), "gw_entry".to_string()),
                    ("gateway_phase".to_string(), "bp".to_string()),
                    ("att_tokens".to_string(), "attL,attR".to_string()),
                ]),
                transactional: false,
                validate_only: true,
            },
        )
        .expect("gateway bp validate-only bad");
        assert_eq!(gateway_bad.output["can_execute"].as_bool(), Some(false));

        let topo_ta_ok = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateRun {
                name: "topo_ta_single_insert".to_string(),
                bindings: HashMap::from([
                    ("vector_seq_id".to_string(), "topo_vector_t".to_string()),
                    ("insert_seq_id".to_string(), "topo_insert_a".to_string()),
                    ("topo_mode".to_string(), "ta".to_string()),
                ]),
                transactional: false,
                validate_only: true,
            },
        )
        .expect("topo ta validate-only ok");
        assert_eq!(topo_ta_ok.output["can_execute"].as_bool(), Some(true));

        let topo_dir_bad = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateRun {
                name: "topo_directional_cacc_single_insert".to_string(),
                bindings: HashMap::from([
                    ("vector_seq_id".to_string(), "topo_vector_t".to_string()),
                    ("insert_seq_id".to_string(), "topo_insert_bad".to_string()),
                    ("topo_mode".to_string(), "directional_cacc".to_string()),
                ]),
                transactional: false,
                validate_only: true,
            },
        )
        .expect("topo directional bad");
        assert_eq!(topo_dir_bad.output["can_execute"].as_bool(), Some(false));

        let ta_ok = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateRun {
                name: "ta_clone_single_insert".to_string(),
                bindings: HashMap::from([
                    ("vector_seq_id".to_string(), "topo_vector_t".to_string()),
                    ("insert_seq_id".to_string(), "topo_insert_a".to_string()),
                    ("tail_mode".to_string(), "ta".to_string()),
                ]),
                transactional: false,
                validate_only: true,
            },
        )
        .expect("ta clone validate-only ok");
        assert_eq!(ta_ok.output["can_execute"].as_bool(), Some(true));

        let gc_bad = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateRun {
                name: "gc_clone_single_insert".to_string(),
                bindings: HashMap::from([
                    ("vector_seq_id".to_string(), "topo_vector_t".to_string()),
                    ("insert_seq_id".to_string(), "topo_insert_a".to_string()),
                    ("tail_mode".to_string(), "gc".to_string()),
                ]),
                transactional: false,
                validate_only: true,
            },
        )
        .expect("gc clone bad");
        assert_eq!(gc_bad.output["can_execute"].as_bool(), Some(false));

        let gc_ok = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateRun {
                name: "gc_clone_single_insert".to_string(),
                bindings: HashMap::from([
                    ("vector_seq_id".to_string(), "gc_vector_c".to_string()),
                    ("insert_seq_id".to_string(), "gc_insert_g".to_string()),
                    ("tail_mode".to_string(), "gc".to_string()),
                ]),
                transactional: false,
                validate_only: true,
            },
        )
        .expect("gc clone ok");
        assert_eq!(gc_ok.output["can_execute"].as_bool(), Some(true));

        let topo_dir_ok = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateRun {
                name: "topo_directional_cacc_single_insert".to_string(),
                bindings: HashMap::from([
                    ("vector_seq_id".to_string(), "topo_vector_t".to_string()),
                    ("insert_seq_id".to_string(), "topo_insert_cacc".to_string()),
                    ("topo_mode".to_string(), "directional_cacc".to_string()),
                ]),
                transactional: false,
                validate_only: true,
            },
        )
        .expect("topo directional ok");
        assert_eq!(topo_dir_ok.output["can_execute"].as_bool(), Some(true));
    }

    #[test]
    fn execute_macros_template_validate_only_applies_infusion_and_nebuilder_checks() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "left".to_string(),
            DNAsequence::from_sequence("ACACACGGGGAAAA").expect("sequence"),
        );
        state.sequences.insert(
            "middle".to_string(),
            DNAsequence::from_sequence("GGGGAAAATTTTCCCC").expect("sequence"),
        );
        state.sequences.insert(
            "right".to_string(),
            DNAsequence::from_sequence("TTTTCCCCGAGAGA").expect("sequence"),
        );
        state.sequences.insert(
            "right_bad".to_string(),
            DNAsequence::from_sequence("CCCCAAAAGAGAGA").expect("sequence"),
        );
        let mut engine = GentleEngine::from_state(state);

        let base = env!("CARGO_MANIFEST_DIR");
        for path in [
            format!(
                "{}/assets/cloning_patterns_catalog/infusion/overlap/infusion_two_fragment_overlap.json",
                base
            ),
            format!(
                "{}/assets/cloning_patterns_catalog/nebuilder_hifi/overlap/nebuilder_multi_fragment_overlap.json",
                base
            ),
        ] {
            execute_shell_command(&mut engine, &ShellCommand::MacrosTemplateImport { path })
                .expect("import overlap template");
        }

        let infusion_ok = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateRun {
                name: "infusion_two_fragment_overlap".to_string(),
                bindings: HashMap::from([
                    ("left_seq_id".to_string(), "left".to_string()),
                    ("right_seq_id".to_string(), "middle".to_string()),
                    ("overlap_bp".to_string(), "8".to_string()),
                ]),
                transactional: false,
                validate_only: true,
            },
        )
        .expect("infusion validate-only ok");
        assert_eq!(infusion_ok.output["can_execute"].as_bool(), Some(true));

        let infusion_bad = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateRun {
                name: "infusion_two_fragment_overlap".to_string(),
                bindings: HashMap::from([
                    ("left_seq_id".to_string(), "left".to_string()),
                    ("right_seq_id".to_string(), "right_bad".to_string()),
                    ("overlap_bp".to_string(), "8".to_string()),
                ]),
                transactional: false,
                validate_only: true,
            },
        )
        .expect("infusion validate-only bad");
        assert_eq!(infusion_bad.output["can_execute"].as_bool(), Some(false));
        let infusion_bad_errors = infusion_bad
            .output
            .get("preflight")
            .and_then(|preflight| preflight.get("errors"))
            .and_then(|rows| rows.as_array())
            .cloned()
            .unwrap_or_default()
            .into_iter()
            .filter_map(|row| row.as_str().map(|s| s.to_string()))
            .collect::<Vec<_>>();
        assert!(
            infusion_bad_errors
                .iter()
                .any(|message| message.contains("In-Fusion overlap mismatch")),
            "expected In-Fusion overlap mismatch error, got: {:?}",
            infusion_bad_errors
        );

        let nebuilder_ok = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateRun {
                name: "nebuilder_multi_fragment_overlap".to_string(),
                bindings: HashMap::from([
                    ("left_seq_id".to_string(), "left".to_string()),
                    ("middle_seq_id".to_string(), "middle".to_string()),
                    ("right_seq_id".to_string(), "right".to_string()),
                    ("overlap_bp".to_string(), "8".to_string()),
                ]),
                transactional: false,
                validate_only: true,
            },
        )
        .expect("nebuilder validate-only ok");
        assert_eq!(nebuilder_ok.output["can_execute"].as_bool(), Some(true));

        let nebuilder_bad = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateRun {
                name: "nebuilder_multi_fragment_overlap".to_string(),
                bindings: HashMap::from([
                    ("left_seq_id".to_string(), "left".to_string()),
                    ("middle_seq_id".to_string(), "middle".to_string()),
                    ("right_seq_id".to_string(), "right".to_string()),
                    ("overlap_bp".to_string(), "40".to_string()),
                ]),
                transactional: false,
                validate_only: true,
            },
        )
        .expect("nebuilder validate-only bad");
        assert_eq!(nebuilder_bad.output["can_execute"].as_bool(), Some(false));
    }

    #[test]
    fn execute_macros_template_run_new_family_packs_transactionally() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "gg_vector".to_string(),
            DNAsequence::from_sequence("AAAAGGTCTCTTTTGGTCTCAAAA").expect("sequence"),
        );
        state.sequences.insert(
            "gg_insert".to_string(),
            DNAsequence::from_sequence("CCCCGGTCTCGGGG").expect("sequence"),
        );
        state.sequences.insert(
            "gw_donor".to_string(),
            DNAsequence::from_sequence("ACGTACGTACGA").expect("sequence"),
        );
        state.sequences.insert(
            "gw_entry".to_string(),
            DNAsequence::from_sequence("TTTTACGTACGT").expect("sequence"),
        );
        state.sequences.insert(
            "topo_vector_t".to_string(),
            DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
        );
        state.sequences.insert(
            "topo_insert_a".to_string(),
            DNAsequence::from_sequence("ACGTACGTACGA").expect("sequence"),
        );
        state.sequences.insert(
            "inf_left".to_string(),
            DNAsequence::from_sequence("ACACACGGGGAAAA").expect("sequence"),
        );
        state.sequences.insert(
            "inf_right".to_string(),
            DNAsequence::from_sequence("GGGGAAAATTTTCCCC").expect("sequence"),
        );
        state.sequences.insert(
            "neb_left".to_string(),
            DNAsequence::from_sequence("TTTTCCCCAAAAGGGG").expect("sequence"),
        );
        state.sequences.insert(
            "neb_right".to_string(),
            DNAsequence::from_sequence("AAAAGGGGCCCTTTAA").expect("sequence"),
        );
        let mut engine = GentleEngine::from_state(state);

        let base = env!("CARGO_MANIFEST_DIR");
        for path in [
            format!(
                "{}/assets/cloning_patterns_catalog/golden_gate/type_iis/golden_gate_single_insert.json",
                base
            ),
            format!(
                "{}/assets/cloning_patterns_catalog/gateway/bp_lr/gateway_bp_single_insert.json",
                base
            ),
            format!(
                "{}/assets/cloning_patterns_catalog/topo/entry/topo_ta_single_insert.json",
                base
            ),
            format!(
                "{}/assets/cloning_patterns_catalog/ta_gc/entry/ta_clone_single_insert.json",
                base
            ),
            format!(
                "{}/assets/cloning_patterns_catalog/infusion/overlap/infusion_two_fragment_overlap.json",
                base
            ),
            format!(
                "{}/assets/cloning_patterns_catalog/nebuilder_hifi/overlap/nebuilder_two_fragment_overlap.json",
                base
            ),
        ] {
            execute_shell_command(&mut engine, &ShellCommand::MacrosTemplateImport { path })
                .expect("import template");
        }

        // Golden Gate sticky-end compatibility depends on fragment pair selection.
        // Probe deterministic digest fragment index combinations and keep the first
        // transactional run that yields at least one ligation product.
        let mut gg_output_id: Option<String> = None;
        let mut gg_last_error: Option<String> = None;
        'gg_search: for vector_fragment in 1..=3 {
            for insert_fragment in 1..=2 {
                let output_id = format!("gg_tx_out_v{}_i{}", vector_fragment, insert_fragment);
                let ligation_prefix = format!("gg_run_v{}_i{}", vector_fragment, insert_fragment);
                let run = execute_shell_command(
                    &mut engine,
                    &ShellCommand::MacrosTemplateRun {
                        name: "golden_gate_single_insert".to_string(),
                        bindings: HashMap::from([
                            ("vector_seq_id".to_string(), "gg_vector".to_string()),
                            ("insert_seq_id".to_string(), "gg_insert".to_string()),
                            ("type_iis_enzyme".to_string(), "Eco31".to_string()),
                            ("vector_fragment".to_string(), vector_fragment.to_string()),
                            ("insert_fragment".to_string(), insert_fragment.to_string()),
                            ("ligation_prefix".to_string(), ligation_prefix),
                            ("output_id".to_string(), output_id.clone()),
                        ]),
                        transactional: true,
                        validate_only: false,
                    },
                );
                match run {
                    Ok(out) if out.state_changed => {
                        gg_output_id = Some(output_id);
                        break 'gg_search;
                    }
                    Ok(_) => {
                        gg_last_error = Some(format!(
                            "golden_gate_single_insert v{} i{} returned state_changed=false",
                            vector_fragment, insert_fragment
                        ));
                    }
                    Err(err) => {
                        gg_last_error = Some(format!(
                            "golden_gate_single_insert v{} i{} failed: {}",
                            vector_fragment, insert_fragment, err
                        ));
                    }
                }
            }
        }
        let gg_output_id = gg_output_id.unwrap_or_else(|| {
            panic!(
                "{}",
                gg_last_error.unwrap_or_else(|| {
                    "No Golden Gate fragment pair produced a ligation product".to_string()
                })
            )
        });

        let runs = vec![
            ShellCommand::MacrosTemplateRun {
                name: "gateway_bp_single_insert".to_string(),
                bindings: HashMap::from([
                    ("donor_seq_id".to_string(), "gw_donor".to_string()),
                    ("entry_vector_seq_id".to_string(), "gw_entry".to_string()),
                    ("gateway_phase".to_string(), "bp".to_string()),
                    ("att_tokens".to_string(), "attB,attP".to_string()),
                    ("assembly_prefix".to_string(), "gw_tx".to_string()),
                    ("output_id".to_string(), "gw_tx_out".to_string()),
                ]),
                transactional: true,
                validate_only: false,
            },
            ShellCommand::MacrosTemplateRun {
                name: "topo_ta_single_insert".to_string(),
                bindings: HashMap::from([
                    ("vector_seq_id".to_string(), "topo_vector_t".to_string()),
                    ("insert_seq_id".to_string(), "topo_insert_a".to_string()),
                    ("topo_mode".to_string(), "ta".to_string()),
                    ("assembly_prefix".to_string(), "topo_tx".to_string()),
                    ("output_id".to_string(), "topo_tx_out".to_string()),
                ]),
                transactional: true,
                validate_only: false,
            },
            ShellCommand::MacrosTemplateRun {
                name: "ta_clone_single_insert".to_string(),
                bindings: HashMap::from([
                    ("vector_seq_id".to_string(), "topo_vector_t".to_string()),
                    ("insert_seq_id".to_string(), "topo_insert_a".to_string()),
                    ("tail_mode".to_string(), "ta".to_string()),
                    ("assembly_prefix".to_string(), "ta_tx".to_string()),
                    ("output_id".to_string(), "ta_tx_out".to_string()),
                ]),
                transactional: true,
                validate_only: false,
            },
            ShellCommand::MacrosTemplateRun {
                name: "infusion_two_fragment_overlap".to_string(),
                bindings: HashMap::from([
                    ("left_seq_id".to_string(), "inf_left".to_string()),
                    ("right_seq_id".to_string(), "inf_right".to_string()),
                    ("overlap_bp".to_string(), "8".to_string()),
                    ("assembly_prefix".to_string(), "inf_tx".to_string()),
                    ("output_id".to_string(), "inf_tx_out".to_string()),
                ]),
                transactional: true,
                validate_only: false,
            },
            ShellCommand::MacrosTemplateRun {
                name: "nebuilder_two_fragment_overlap".to_string(),
                bindings: HashMap::from([
                    ("left_seq_id".to_string(), "neb_left".to_string()),
                    ("right_seq_id".to_string(), "neb_right".to_string()),
                    ("overlap_bp".to_string(), "8".to_string()),
                    ("assembly_prefix".to_string(), "neb_tx".to_string()),
                    ("output_id".to_string(), "neb_tx_out".to_string()),
                ]),
                transactional: true,
                validate_only: false,
            },
        ];
        for run_cmd in runs {
            let out = execute_shell_command(&mut engine, &run_cmd).expect("transactional run");
            assert!(out.state_changed);
        }

        assert!(engine.state().sequences.contains_key(&gg_output_id));
        for seq_id in [
            "gw_tx_out",
            "topo_tx_out",
            "ta_tx_out",
            "inf_tx_out",
            "neb_tx_out",
        ] {
            assert!(engine.state().sequences.contains_key(seq_id));
        }
    }

    #[test]
    fn execute_macros_template_run_gibson_preview_creates_deterministic_outputs() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "left".to_string(),
            DNAsequence::from_sequence("ATGCCGTTACCGGTTAAACCCGGGTTT").expect("sequence"),
        );
        state.sequences.insert(
            "right".to_string(),
            DNAsequence::from_sequence("AACCCGGGTTTGGGAAATTTCCCGGG").expect("sequence"),
        );
        let mut engine = GentleEngine::from_state(state);

        let path = format!(
            "{}/assets/cloning_patterns_catalog/gibson/overlap_assembly/gibson_two_fragment_overlap_preview.json",
            env!("CARGO_MANIFEST_DIR")
        );
        execute_shell_command(&mut engine, &ShellCommand::MacrosTemplateImport { path })
            .expect("import gibson template");

        let run = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateRun {
                name: "gibson_two_fragment_overlap_preview".to_string(),
                bindings: HashMap::from([
                    ("left_seq_id".to_string(), "left".to_string()),
                    ("right_seq_id".to_string(), "right".to_string()),
                    ("overlap_bp".to_string(), "11".to_string()),
                    ("assembly_prefix".to_string(), "gib_preview".to_string()),
                    ("output_id".to_string(), "gib_forward".to_string()),
                ]),
                transactional: true,
                validate_only: false,
            },
        )
        .expect("run gibson template");
        assert!(run.state_changed);
        assert!(engine.state().sequences.contains_key("gib_preview_1"));
        assert!(engine.state().sequences.contains_key("gib_preview_2"));
        assert!(engine.state().sequences.contains_key("gib_forward"));
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
    fn execute_primers_design_list_show_export() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "tpl".to_string(),
            DNAsequence::from_sequence(
                "GGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
            )
            .expect("sequence"),
        );
        let mut engine = GentleEngine::from_state(state);
        let tmp = tempdir().expect("tempdir");
        let export_path = tmp.path().join("primer_report.json");
        let request = serde_json::to_string(&Operation::DesignPrimerPairs {
            template: "tpl".to_string(),
            roi_start_0based: 30,
            roi_end_0based: 70,
            forward: crate::engine::PrimerDesignSideConstraint {
                min_length: 20,
                max_length: 20,
                location_0based: Some(5),
                start_0based: None,
                end_0based: None,
                min_tm_c: 40.0,
                max_tm_c: 90.0,
                min_gc_fraction: 0.0,
                max_gc_fraction: 1.0,
                max_anneal_hits: 10,
                ..Default::default()
            },
            reverse: crate::engine::PrimerDesignSideConstraint {
                min_length: 20,
                max_length: 20,
                location_0based: Some(60),
                start_0based: None,
                end_0based: None,
                min_tm_c: 40.0,
                max_tm_c: 90.0,
                min_gc_fraction: 0.0,
                max_gc_fraction: 1.0,
                max_anneal_hits: 10,
                ..Default::default()
            },
            min_amplicon_bp: 40,
            max_amplicon_bp: 130,
            pair_constraints: crate::engine::PrimerDesignPairConstraint::default(),
            max_tm_delta_c: Some(50.0),
            max_pairs: Some(10),
            report_id: Some("tp73_roi".to_string()),
        })
        .expect("serialize request");

        let design = execute_shell_command(
            &mut engine,
            &ShellCommand::PrimersDesign {
                request_json: request,
                backend: Some(PrimerDesignBackend::Internal),
                primer3_executable: None,
            },
        )
        .expect("primers design");
        assert!(design.state_changed);
        let report_id = design.output["report"]["report_id"]
            .as_str()
            .unwrap_or_default()
            .to_string();
        assert_eq!(report_id, "tp73_roi");

        let listed = execute_shell_command(&mut engine, &ShellCommand::PrimersListReports)
            .expect("primers list-reports");
        assert!(!listed.state_changed);
        assert_eq!(
            listed.output["schema"].as_str(),
            Some("gentle.primer_design_report_list.v1")
        );
        assert_eq!(listed.output["report_count"].as_u64(), Some(1));

        let shown = execute_shell_command(
            &mut engine,
            &ShellCommand::PrimersShowReport {
                report_id: report_id.clone(),
            },
        )
        .expect("primers show-report");
        assert!(!shown.state_changed);
        assert_eq!(
            shown.output["report"]["report_id"].as_str(),
            Some("tp73_roi")
        );

        let exported = execute_shell_command(
            &mut engine,
            &ShellCommand::PrimersExportReport {
                report_id,
                path: export_path.to_string_lossy().to_string(),
            },
        )
        .expect("primers export-report");
        assert!(!exported.state_changed);
        assert_eq!(
            exported.output["schema"].as_str(),
            Some("gentle.primer_design_report_export.v1")
        );
        let text = fs::read_to_string(&export_path).expect("read export");
        assert!(text.contains("gentle.primer_design_report.v1"));
    }

    #[cfg(unix)]
    #[test]
    fn execute_primers_design_internal_vs_primer3_fixture_normalization_parity() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "tpl".to_string(),
            DNAsequence::from_sequence(
                "ACGTTGCATGTCAGTACGATCGTACGTAGCTAGTCGATCGTACGATCGTAGCTAGCATCGATGCTAGCTAGTACGTAGCATCGATCGTAGCTAGCATGCTAGCTAGTCGATCGATCGTACGATCG",
            )
            .expect("sequence"),
        );
        let mut engine = GentleEngine::from_state(state);
        let tmp = tempdir().expect("tempdir");
        let fixture_path = Path::new(&primer3_fixture_path("pairs.location_5_60.kv")).to_path_buf();
        let fake_primer3 = install_fake_primer3(tmp.path(), &fixture_path);

        let forward = crate::engine::PrimerDesignSideConstraint {
            min_length: 20,
            max_length: 20,
            location_0based: Some(5),
            start_0based: None,
            end_0based: None,
            min_tm_c: 0.0,
            max_tm_c: 100.0,
            min_gc_fraction: 0.0,
            max_gc_fraction: 1.0,
            max_anneal_hits: 1000,
            ..Default::default()
        };
        let reverse = crate::engine::PrimerDesignSideConstraint {
            min_length: 20,
            max_length: 20,
            location_0based: Some(60),
            start_0based: None,
            end_0based: None,
            min_tm_c: 0.0,
            max_tm_c: 100.0,
            min_gc_fraction: 0.0,
            max_gc_fraction: 1.0,
            max_anneal_hits: 1000,
            ..Default::default()
        };

        let internal_request = serde_json::to_string(&Operation::DesignPrimerPairs {
            template: "tpl".to_string(),
            roi_start_0based: 40,
            roi_end_0based: 80,
            forward: forward.clone(),
            reverse: reverse.clone(),
            min_amplicon_bp: 40,
            max_amplicon_bp: 150,
            pair_constraints: crate::engine::PrimerDesignPairConstraint::default(),
            max_tm_delta_c: Some(100.0),
            max_pairs: Some(10),
            report_id: Some("internal_norm".to_string()),
        })
        .expect("serialize internal request");
        execute_shell_command(
            &mut engine,
            &ShellCommand::PrimersDesign {
                request_json: internal_request,
                backend: Some(PrimerDesignBackend::Internal),
                primer3_executable: None,
            },
        )
        .expect("primers design internal");
        let internal = engine
            .get_primer_design_report("internal_norm")
            .expect("internal report");

        let primer3_request = serde_json::to_string(&Operation::DesignPrimerPairs {
            template: "tpl".to_string(),
            roi_start_0based: 40,
            roi_end_0based: 80,
            forward,
            reverse,
            min_amplicon_bp: 40,
            max_amplicon_bp: 150,
            pair_constraints: crate::engine::PrimerDesignPairConstraint::default(),
            max_tm_delta_c: Some(100.0),
            max_pairs: Some(10),
            report_id: Some("primer3_norm".to_string()),
        })
        .expect("serialize primer3 request");
        execute_shell_command(
            &mut engine,
            &ShellCommand::PrimersDesign {
                request_json: primer3_request,
                backend: Some(PrimerDesignBackend::Primer3),
                primer3_executable: Some(fake_primer3.clone()),
            },
        )
        .expect("primers design primer3");
        let primer3 = engine
            .get_primer_design_report("primer3_norm")
            .expect("primer3 report");

        assert_eq!(internal.template, primer3.template);
        assert_eq!(internal.roi_start_0based, primer3.roi_start_0based);
        assert_eq!(internal.roi_end_0based, primer3.roi_end_0based);
        assert_eq!(internal.min_amplicon_bp, primer3.min_amplicon_bp);
        assert_eq!(internal.max_amplicon_bp, primer3.max_amplicon_bp);
        assert_eq!(internal.pair_count, 1);
        assert_eq!(primer3.pair_count, 1);
        let internal_pair = internal.pairs.first().expect("internal pair");
        let primer3_pair = primer3.pairs.first().expect("primer3 pair");
        assert_eq!(
            internal_pair.forward.start_0based,
            primer3_pair.forward.start_0based
        );
        assert_eq!(
            internal_pair.forward.end_0based_exclusive,
            primer3_pair.forward.end_0based_exclusive
        );
        assert_eq!(
            internal_pair.reverse.start_0based,
            primer3_pair.reverse.start_0based
        );
        assert_eq!(
            internal_pair.reverse.end_0based_exclusive,
            primer3_pair.reverse.end_0based_exclusive
        );
        assert_eq!(
            internal_pair.amplicon_length_bp,
            primer3_pair.amplicon_length_bp
        );
        assert_eq!(primer3.backend.requested, "primer3");
        assert_eq!(primer3.backend.used, "primer3");
        assert_eq!(
            primer3.backend.primer3_executable.as_deref(),
            Some(fake_primer3.as_str())
        );
        assert_eq!(
            primer3.backend.primer3_version.as_deref(),
            Some("primer3_core synthetic-fixture 2.6.1")
        );
    }

    #[cfg(unix)]
    #[test]
    fn execute_primers_preflight_reports_reachable_primer3() {
        let mut engine = GentleEngine::from_state(ProjectState::default());
        let tmp = tempdir().expect("tempdir");
        let fixture_path = Path::new(&primer3_fixture_path("pairs.location_5_60.kv")).to_path_buf();
        let fake_primer3 = install_fake_primer3(tmp.path(), &fixture_path);
        let out = execute_shell_command(
            &mut engine,
            &ShellCommand::PrimersPreflight {
                backend: Some(PrimerDesignBackend::Primer3),
                primer3_executable: Some(fake_primer3.clone()),
            },
        )
        .expect("primers preflight");
        assert!(!out.state_changed);
        assert_eq!(
            out.output["schema"].as_str(),
            Some("gentle.primer3_preflight.v1")
        );
        assert_eq!(out.output["preflight"]["reachable"].as_bool(), Some(true));
        assert_eq!(
            out.output["preflight"]["version_probe_ok"].as_bool(),
            Some(true)
        );
        assert_eq!(out.output["preflight"]["backend"].as_str(), Some("primer3"));
        assert_eq!(
            out.output["preflight"]["executable"].as_str(),
            Some(fake_primer3.as_str())
        );
    }

    #[test]
    fn execute_primers_design_qpcr_list_show_export() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "tpl".to_string(),
            DNAsequence::from_sequence(
                "GGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
            )
            .expect("sequence"),
        );
        let mut engine = GentleEngine::from_state(state);
        let tmp = tempdir().expect("tempdir");
        let export_path = tmp.path().join("qpcr_report.json");
        let request = serde_json::to_string(&Operation::DesignQpcrAssays {
            template: "tpl".to_string(),
            roi_start_0based: 30,
            roi_end_0based: 70,
            forward: crate::engine::PrimerDesignSideConstraint {
                min_length: 20,
                max_length: 20,
                location_0based: Some(5),
                start_0based: None,
                end_0based: None,
                min_tm_c: 40.0,
                max_tm_c: 90.0,
                min_gc_fraction: 0.0,
                max_gc_fraction: 1.0,
                max_anneal_hits: 100,
                ..Default::default()
            },
            reverse: crate::engine::PrimerDesignSideConstraint {
                min_length: 20,
                max_length: 20,
                location_0based: Some(60),
                start_0based: None,
                end_0based: None,
                min_tm_c: 40.0,
                max_tm_c: 90.0,
                min_gc_fraction: 0.0,
                max_gc_fraction: 1.0,
                max_anneal_hits: 100,
                ..Default::default()
            },
            probe: crate::engine::PrimerDesignSideConstraint {
                min_length: 20,
                max_length: 20,
                location_0based: Some(35),
                start_0based: None,
                end_0based: None,
                min_tm_c: 40.0,
                max_tm_c: 90.0,
                min_gc_fraction: 0.0,
                max_gc_fraction: 1.0,
                max_anneal_hits: 100,
                ..Default::default()
            },
            min_amplicon_bp: 40,
            max_amplicon_bp: 130,
            pair_constraints: crate::engine::PrimerDesignPairConstraint::default(),
            max_tm_delta_c: Some(50.0),
            max_probe_tm_delta_c: Some(50.0),
            max_assays: Some(10),
            report_id: Some("tp73_qpcr".to_string()),
        })
        .expect("serialize request");

        let design = execute_shell_command(
            &mut engine,
            &ShellCommand::PrimersDesignQpcr {
                request_json: request,
                backend: Some(PrimerDesignBackend::Internal),
                primer3_executable: None,
            },
        )
        .expect("primers design-qpcr");
        assert!(design.state_changed);
        let report_id = design.output["report"]["report_id"]
            .as_str()
            .unwrap_or_default()
            .to_string();
        assert_eq!(report_id, "tp73_qpcr");

        let listed = execute_shell_command(&mut engine, &ShellCommand::PrimersListQpcrReports)
            .expect("primers list-qpcr-reports");
        assert!(!listed.state_changed);
        assert_eq!(
            listed.output["schema"].as_str(),
            Some("gentle.qpcr_design_report_list.v1")
        );
        assert_eq!(listed.output["report_count"].as_u64(), Some(1));

        let shown = execute_shell_command(
            &mut engine,
            &ShellCommand::PrimersShowQpcrReport {
                report_id: report_id.clone(),
            },
        )
        .expect("primers show-qpcr-report");
        assert!(!shown.state_changed);
        assert_eq!(
            shown.output["report"]["report_id"].as_str(),
            Some("tp73_qpcr")
        );

        let exported = execute_shell_command(
            &mut engine,
            &ShellCommand::PrimersExportQpcrReport {
                report_id,
                path: export_path.to_string_lossy().to_string(),
            },
        )
        .expect("primers export-qpcr-report");
        assert!(!exported.state_changed);
        assert_eq!(
            exported.output["schema"].as_str(),
            Some("gentle.qpcr_design_report_export.v1")
        );
        let text = fs::read_to_string(&export_path).expect("read export");
        assert!(text.contains("gentle.qpcr_design_report.v1"));
    }

    #[test]
    fn execute_primers_seed_from_feature_and_splicing() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("seq_a".to_string(), tp53_isoform_test_sequence());
        let mut engine = GentleEngine::from_state(state);
        let feature_id = engine
            .state()
            .sequences
            .get("seq_a")
            .expect("sequence present")
            .features()
            .iter()
            .position(|feature| feature.kind.to_string().eq_ignore_ascii_case("mRNA"))
            .expect("mRNA feature id");

        let seeded_feature = execute_shell_command(
            &mut engine,
            &ShellCommand::PrimersSeedFromFeature {
                seq_id: "seq_a".to_string(),
                feature_id,
            },
        )
        .expect("seed from feature");
        assert!(!seeded_feature.state_changed);
        assert_eq!(
            seeded_feature.output["schema"].as_str(),
            Some("gentle.primer_seed_request.v1")
        );
        assert_eq!(
            seeded_feature.output["source"]["kind"].as_str(),
            Some("feature")
        );
        let feature_start = seeded_feature.output["roi_start_0based"]
            .as_u64()
            .expect("feature start roi") as usize;
        let feature_end = seeded_feature.output["roi_end_0based_exclusive"]
            .as_u64()
            .expect("feature end roi") as usize;
        assert!(feature_end > feature_start);
        assert_eq!(
            seeded_feature.output["operations"]["design_primer_pairs"]["DesignPrimerPairs"]
                ["template"]
                .as_str(),
            Some("seq_a")
        );
        assert_eq!(
            seeded_feature.output["operations"]["design_primer_pairs"]["DesignPrimerPairs"]
                ["roi_start_0based"]
                .as_u64(),
            Some(feature_start as u64)
        );
        assert_eq!(
            seeded_feature.output["operations"]["design_primer_pairs"]["DesignPrimerPairs"]
                ["roi_end_0based"]
                .as_u64(),
            Some(feature_end as u64)
        );

        let seeded_splicing = execute_shell_command(
            &mut engine,
            &ShellCommand::PrimersSeedFromSplicing {
                seq_id: "seq_a".to_string(),
                feature_id,
            },
        )
        .expect("seed from splicing");
        assert!(!seeded_splicing.state_changed);
        assert_eq!(
            seeded_splicing.output["schema"].as_str(),
            Some("gentle.primer_seed_request.v1")
        );
        assert_eq!(
            seeded_splicing.output["source"]["kind"].as_str(),
            Some("splicing")
        );
        let splicing_start = seeded_splicing.output["roi_start_0based"]
            .as_u64()
            .expect("splicing start roi") as usize;
        let splicing_end = seeded_splicing.output["roi_end_0based_exclusive"]
            .as_u64()
            .expect("splicing end roi") as usize;
        assert!(splicing_end > splicing_start);
        assert_eq!(splicing_start, feature_start);
        assert_eq!(splicing_end, feature_end);
        assert_eq!(
            seeded_splicing.output["operations"]["design_qpcr_assays"]["DesignQpcrAssays"]
                ["template"]
                .as_str(),
            Some("seq_a")
        );
    }

    #[test]
    fn execute_async_blast_start_and_status_reports_failure_for_missing_genome() {
        let _guard = BLAST_ASYNC_TEST_MUTEX.lock().expect("blast mutex");
        let mut engine = GentleEngine::new();
        let start = execute_shell_command(
            &mut engine,
            &ShellCommand::ReferenceBlastAsyncStart {
                helper_mode: true,
                genome_id: "missing_helper".to_string(),
                query_sequence: "ACGTACGT".to_string(),
                max_hits: 5,
                max_hits_explicit: true,
                task: Some("blastn-short".to_string()),
                request_options_json: None,
                catalog_path: None,
                cache_dir: None,
            },
        )
        .expect("start async blast");
        assert!(!start.state_changed);
        let job_id = start
            .output
            .get("job")
            .and_then(|job| job.get("job_id"))
            .and_then(|value| value.as_str())
            .unwrap_or_default()
            .to_string();
        assert!(!job_id.is_empty());

        let mut terminal_state = String::new();
        for _ in 0..40 {
            let status = execute_shell_command(
                &mut engine,
                &ShellCommand::ReferenceBlastAsyncStatus {
                    helper_mode: true,
                    job_id: job_id.clone(),
                    include_report: true,
                },
            )
            .expect("status");
            terminal_state = status
                .output
                .get("job")
                .and_then(|job| job.get("state"))
                .and_then(|value| value.as_str())
                .unwrap_or_default()
                .to_string();
            if matches!(
                terminal_state.as_str(),
                "completed" | "failed" | "cancelled"
            ) {
                break;
            }
            std::thread::sleep(std::time::Duration::from_millis(20));
        }
        assert!(
            matches!(terminal_state.as_str(), "failed" | "cancelled"),
            "unexpected async blast terminal state: {}",
            terminal_state
        );
    }

    #[test]
    fn execute_export_run_bundle_writes_schema_json() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "seqA".to_string(),
            DNAsequence::from_sequence("ACGTACGT").expect("sequence"),
        );
        let mut engine = GentleEngine::from_state(state);
        engine
            .apply(crate::engine::Operation::Reverse {
                input: "seqA".to_string(),
                output_id: Some("seqA_rev".to_string()),
            })
            .expect("reverse");

        let tmp = tempfile::NamedTempFile::new().expect("tmp file");
        let path = tmp.path().with_extension("shell.run_bundle.json");
        let path_text = path.display().to_string();
        let out = execute_shell_command(
            &mut engine,
            &ShellCommand::ExportRunBundle {
                output: path_text.clone(),
                run_id: Some("interactive".to_string()),
            },
        )
        .expect("export run bundle");
        assert!(!out.state_changed);
        let text = fs::read_to_string(path_text).expect("read bundle");
        let value: serde_json::Value = serde_json::from_str(&text).expect("parse bundle");
        assert_eq!(
            value.get("schema").and_then(|v| v.as_str()),
            Some("gentle.process_run_bundle.v1")
        );
        assert_eq!(
            value
                .get("run_id_filter")
                .and_then(|v| v.as_str())
                .unwrap_or_default(),
            "interactive"
        );
    }

    #[test]
    fn execute_async_blast_start_queues_when_capacity_is_reached() {
        with_blast_async_test_overrides(1, 200, || {
            let mut engine = GentleEngine::new();
            let start_one = execute_shell_command(
                &mut engine,
                &ShellCommand::ReferenceBlastAsyncStart {
                    helper_mode: true,
                    genome_id: "missing_helper".to_string(),
                    query_sequence: "ACGTACGT".to_string(),
                    max_hits: 5,
                    max_hits_explicit: true,
                    task: Some("blastn-short".to_string()),
                    request_options_json: None,
                    catalog_path: None,
                    cache_dir: None,
                },
            )
            .expect("start first async blast");
            let job_one = start_one.output["job"]["job_id"]
                .as_str()
                .unwrap_or_default()
                .to_string();
            assert!(!job_one.is_empty());
            assert_eq!(
                start_one.output["job"]["max_concurrent_jobs"].as_u64(),
                Some(1)
            );

            let start_two = execute_shell_command(
                &mut engine,
                &ShellCommand::ReferenceBlastAsyncStart {
                    helper_mode: true,
                    genome_id: "missing_helper".to_string(),
                    query_sequence: "TTTTAAAA".to_string(),
                    max_hits: 5,
                    max_hits_explicit: true,
                    task: Some("blastn-short".to_string()),
                    request_options_json: None,
                    catalog_path: None,
                    cache_dir: None,
                },
            )
            .expect("start second async blast");
            let job_two = start_two.output["job"]["job_id"]
                .as_str()
                .unwrap_or_default()
                .to_string();
            assert!(!job_two.is_empty());
            assert_eq!(start_two.output["job"]["state"].as_str(), Some("queued"));
            assert_eq!(start_two.output["job"]["queue_position"].as_u64(), Some(1));

            let _cancel = execute_shell_command(
                &mut engine,
                &ShellCommand::ReferenceBlastAsyncCancel {
                    helper_mode: true,
                    job_id: job_one,
                },
            )
            .expect("cancel first job");

            let mut observed_non_queued = false;
            for _ in 0..80 {
                let status = execute_shell_command(
                    &mut engine,
                    &ShellCommand::ReferenceBlastAsyncStatus {
                        helper_mode: true,
                        job_id: job_two.clone(),
                        include_report: false,
                    },
                )
                .expect("status second job");
                let state = status.output["job"]["state"]
                    .as_str()
                    .unwrap_or_default()
                    .to_string();
                if state != "queued" {
                    observed_non_queued = true;
                    break;
                }
                std::thread::sleep(std::time::Duration::from_millis(15));
            }
            assert!(
                observed_non_queued,
                "expected queued second job to be dispatched after first slot freed"
            );
        });
    }

    #[test]
    fn execute_async_blast_list_reports_scheduler_metadata() {
        with_blast_async_test_overrides(2, 0, || {
            let mut engine = GentleEngine::new();
            let start = execute_shell_command(
                &mut engine,
                &ShellCommand::ReferenceBlastAsyncStart {
                    helper_mode: true,
                    genome_id: "missing_helper".to_string(),
                    query_sequence: "ACGTACGT".to_string(),
                    max_hits: 5,
                    max_hits_explicit: true,
                    task: Some("blastn-short".to_string()),
                    request_options_json: None,
                    catalog_path: None,
                    cache_dir: None,
                },
            )
            .expect("start async blast");
            assert!(!start.state_changed);

            let listed = execute_shell_command(
                &mut engine,
                &ShellCommand::ReferenceBlastAsyncList { helper_mode: true },
            )
            .expect("list async blast jobs");
            assert!(!listed.state_changed);
            assert_eq!(listed.output["job_count"].as_u64(), Some(1));
            let jobs = listed.output["jobs"]
                .as_array()
                .expect("jobs should be array");
            assert_eq!(jobs.len(), 1);
            assert_eq!(jobs[0]["max_concurrent_jobs"].as_u64(), Some(2));
            assert!(jobs[0]["running_jobs"].as_u64().unwrap_or(0) <= 2);
            assert!(jobs[0]["queued_jobs"].as_u64().unwrap_or(0) <= 1);
        });
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
    fn execute_macros_run_records_inline_macro_instance() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "seqA".to_string(),
            DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
        );
        let mut engine = GentleEngine::from_state(state);
        let run = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosRun {
                script: r#"op {"Reverse":{"input":"seqA","output_id":"seqA_rev"}}"#.to_string(),
                transactional: false,
            },
        )
        .expect("run inline macro");
        assert!(run.state_changed);
        assert!(run.output["macro_instance_id"].as_str().is_some());
        assert_eq!(engine.state().lineage.macro_instances.len(), 1);
        let instance = &engine.state().lineage.macro_instances[0];
        assert!(instance.template_name.is_none());
        assert!(!instance.expanded_op_ids.is_empty());
    }

    #[test]
    fn execute_macros_run_failed_records_macro_instance() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "seqA".to_string(),
            DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
        );
        let mut engine = GentleEngine::from_state(state);
        let err = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosRun {
                script: r#"op {"Reverse":{"input":"missing","output_id":"bad"}}"#.to_string(),
                transactional: false,
            },
        )
        .expect_err("macro should fail");
        assert!(err.contains("macro_instance_id="));
        assert_eq!(engine.state().lineage.macro_instances.len(), 1);
        let instance = &engine.state().lineage.macro_instances[0];
        assert_eq!(instance.status, MacroInstanceStatus::Failed);
        assert!(instance.status_message.is_some());
        assert!(instance.expanded_op_ids.is_empty());
    }

    #[test]
    fn execute_macros_instance_list_and_show() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "seqA".to_string(),
            DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
        );
        let mut engine = GentleEngine::from_state(state);
        execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosRun {
                script: r#"op {"Reverse":{"input":"seqA","output_id":"seqA_rev"}}"#.to_string(),
                transactional: false,
            },
        )
        .expect("run macro");

        let listed = execute_shell_command(&mut engine, &ShellCommand::MacrosInstanceList)
            .expect("list macro instances");
        assert!(!listed.state_changed);
        assert_eq!(
            listed.output["schema"].as_str(),
            Some("gentle.lineage_macro_instances.v1")
        );
        let instance_id = listed
            .output
            .get("instances")
            .and_then(|rows| rows.as_array())
            .and_then(|rows| rows.first())
            .and_then(|row| row.get("macro_instance_id"))
            .and_then(|v| v.as_str())
            .unwrap_or_default()
            .to_string();
        assert!(!instance_id.is_empty());

        let shown = execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosInstanceShow {
                macro_instance_id: instance_id,
            },
        )
        .expect("show macro instance");
        assert!(!shown.state_changed);
        assert_eq!(
            shown.output["schema"].as_str(),
            Some("gentle.lineage_macro_instance.v1")
        );
        assert!(
            shown
                .output
                .get("instance")
                .and_then(|row| row.get("expanded_op_ids"))
                .and_then(|v| v.as_array())
                .map(|rows| !rows.is_empty())
                .unwrap_or(false)
        );
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
    fn execute_agent_suggestions_allows_blast_shell_route() {
        let mut engine = GentleEngine::from_state(ProjectState::default());
        let suggestions = vec![crate::agent_bridge::AgentSuggestedCommand {
            title: Some("BLAST".to_string()),
            rationale: Some("Check specificity".to_string()),
            command: "genomes blast __missing_genome__ ACGT".to_string(),
            execution: AgentExecutionIntent::Auto,
        }];
        let execute_indices = std::collections::BTreeSet::new();
        let (changed, reports) = execute_agent_suggested_commands(
            &mut engine,
            &suggestions,
            false,
            &execute_indices,
            true,
            &ShellExecutionOptions::default(),
        );
        assert!(!changed);
        assert_eq!(reports.len(), 1);
        let row = &reports[0];
        assert!(row.executed);
        assert!(!row.ok);
        let error = row
            .error
            .as_deref()
            .unwrap_or_default()
            .to_ascii_lowercase();
        assert!(
            !error.contains("agent-to-agent"),
            "BLAST route should not be blocked by recursion guardrail: {error}"
        );
        assert!(
            !error.contains("blocked in this context"),
            "BLAST route should be executable in suggested-command context: {error}"
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
    fn parse_genomes_extract_gene_with_scope_and_cap() {
        let cmd = parse_shell_line(
            "genomes extract-gene ToyGenome MYGENE --annotation-scope full --max-annotation-features 120 --output-id out",
        )
        .expect("parse genomes extract-gene with scope and cap");
        match cmd {
            ShellCommand::ReferenceExtractGene {
                helper_mode,
                genome_id,
                gene_query,
                output_id,
                annotation_scope,
                max_annotation_features,
                include_genomic_annotation,
                ..
            } => {
                assert!(!helper_mode);
                assert_eq!(genome_id, "ToyGenome".to_string());
                assert_eq!(gene_query, "MYGENE".to_string());
                assert_eq!(output_id, Some("out".to_string()));
                assert_eq!(annotation_scope, Some(GenomeAnnotationScope::Full));
                assert_eq!(max_annotation_features, Some(120));
                assert_eq!(include_genomic_annotation, None);
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_genomes_extract_gene_with_annotation_flag() {
        let cmd = parse_shell_line(
            "genomes extract-gene ToyGenome MYGENE --include-genomic-annotation --output-id out",
        )
        .expect("parse genomes extract-gene with annotation flag");
        match cmd {
            ShellCommand::ReferenceExtractGene {
                helper_mode,
                genome_id,
                gene_query,
                output_id,
                annotation_scope,
                include_genomic_annotation,
                ..
            } => {
                assert!(!helper_mode);
                assert_eq!(genome_id, "ToyGenome".to_string());
                assert_eq!(gene_query, "MYGENE".to_string());
                assert_eq!(output_id, Some("out".to_string()));
                assert_eq!(annotation_scope, Some(GenomeAnnotationScope::Core));
                assert_eq!(include_genomic_annotation, Some(true));
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_genomes_extract_region_with_annotation_flag() {
        let cmd = parse_shell_line(
            "genomes extract-region ToyGenome chr1 100 250 --include-genomic-annotation --output-id out",
        )
        .expect("parse genomes extract-region with annotation flag");
        match cmd {
            ShellCommand::ReferenceExtractRegion {
                helper_mode,
                genome_id,
                chromosome,
                start_1based,
                end_1based,
                output_id,
                annotation_scope,
                include_genomic_annotation,
                ..
            } => {
                assert!(!helper_mode);
                assert_eq!(genome_id, "ToyGenome".to_string());
                assert_eq!(chromosome, "chr1".to_string());
                assert_eq!(start_1based, 100);
                assert_eq!(end_1based, 250);
                assert_eq!(output_id, Some("out".to_string()));
                assert_eq!(annotation_scope, Some(GenomeAnnotationScope::Core));
                assert_eq!(include_genomic_annotation, Some(true));
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_genomes_extract_region_with_scope_and_cap() {
        let cmd = parse_shell_line(
            "genomes extract-region ToyGenome chr1 100 250 --annotation-scope full --max-annotation-features 1200 --output-id out",
        )
        .expect("parse genomes extract-region with scope and cap");
        match cmd {
            ShellCommand::ReferenceExtractRegion {
                helper_mode,
                genome_id,
                chromosome,
                start_1based,
                end_1based,
                output_id,
                annotation_scope,
                max_annotation_features,
                include_genomic_annotation,
                ..
            } => {
                assert!(!helper_mode);
                assert_eq!(genome_id, "ToyGenome".to_string());
                assert_eq!(chromosome, "chr1".to_string());
                assert_eq!(start_1based, 100);
                assert_eq!(end_1based, 250);
                assert_eq!(output_id, Some("out".to_string()));
                assert_eq!(annotation_scope, Some(GenomeAnnotationScope::Full));
                assert_eq!(max_annotation_features, Some(1200));
                assert_eq!(include_genomic_annotation, None);
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
                prepared_genome_id,
            } => {
                assert!(!helper_mode);
                assert_eq!(seq_id, "tp73".to_string());
                assert_eq!(side, GenomeAnchorSide::FivePrime);
                assert_eq!(length_bp, 150);
                assert_eq!(output_id, Some("tp73_ext".to_string()));
                assert_eq!(catalog_path, Some("c.json".to_string()));
                assert_eq!(cache_dir, Some("cache".to_string()));
                assert_eq!(prepared_genome_id, None);
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_genomes_extend_anchor_with_prepared_genome() {
        let cmd = parse_shell_line(
            "genomes extend-anchor tp73 3p 200 --prepared-genome Human_GRCh38 --catalog c.json",
        )
        .expect("parse genomes extend-anchor with prepared genome");
        match cmd {
            ShellCommand::ReferenceExtendAnchor {
                helper_mode,
                seq_id,
                side,
                length_bp,
                output_id,
                catalog_path,
                cache_dir,
                prepared_genome_id,
            } => {
                assert!(!helper_mode);
                assert_eq!(seq_id, "tp73".to_string());
                assert_eq!(side, GenomeAnchorSide::ThreePrime);
                assert_eq!(length_bp, 200);
                assert_eq!(output_id, None);
                assert_eq!(catalog_path, Some("c.json".to_string()));
                assert_eq!(cache_dir, None);
                assert_eq!(prepared_genome_id, Some("Human_GRCh38".to_string()));
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_genomes_verify_anchor_with_prepared_genome() {
        let cmd = parse_shell_line(
            "genomes verify-anchor tp73 --catalog c.json --cache-dir cache --prepared-genome Human_GRCh38",
        )
        .expect("parse genomes verify-anchor");
        match cmd {
            ShellCommand::ReferenceVerifyAnchor {
                helper_mode,
                seq_id,
                catalog_path,
                cache_dir,
                prepared_genome_id,
            } => {
                assert!(!helper_mode);
                assert_eq!(seq_id, "tp73".to_string());
                assert_eq!(catalog_path, Some("c.json".to_string()));
                assert_eq!(cache_dir, Some("cache".to_string()));
                assert_eq!(prepared_genome_id, Some("Human_GRCh38".to_string()));
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
                annotation_scope: None,
                max_annotation_features: None,
                include_genomic_annotation: None,
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
                prepared_genome_id: None,
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
    fn execute_genomes_extract_region_default_scope_core_with_telemetry() {
        let td = tempdir().expect("tempdir");
        let fasta = td.path().join("toy.fa");
        let gtf = td.path().join("toy.gtf");
        fs::write(&fasta, ">chr1\nACGTACGTACGTACGT\n").expect("write fasta");
        fs::write(
            &gtf,
            concat!(
                "chr1\tsrc\tgene\t1\t16\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n",
                "chr1\tsrc\ttranscript\t1\t16\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\"; transcript_id \"TX1\";\n",
                "chr1\tsrc\texon\t1\t8\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\"; transcript_id \"TX1\"; exon_number \"1\";\n",
            ),
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

        let command = parse_shell_line(&format!(
            "genomes extract-region ToyGenome chr1 1 16 --output-id shell_slice_default --catalog {}",
            catalog_path
        ))
        .expect("parse extract-region");
        let out = execute_shell_command(&mut engine, &command).expect("execute extract-region");
        assert!(out.state_changed);
        let telemetry = out
            .output
            .get("result")
            .and_then(|v| v.get("genome_annotation_projection"))
            .and_then(|v| v.as_object())
            .expect("annotation projection telemetry");
        assert_eq!(
            telemetry.get("requested_scope").and_then(|v| v.as_str()),
            Some("core")
        );
        assert_eq!(
            telemetry.get("effective_scope").and_then(|v| v.as_str()),
            Some("core")
        );
        assert!(
            telemetry
                .get("attached_feature_count")
                .and_then(|v| v.as_u64())
                .unwrap_or(0)
                >= 2
        );
        let seq = engine
            .state()
            .sequences
            .get("shell_slice_default")
            .expect("sequence created by extract-region");
        assert!(
            seq.features()
                .iter()
                .any(|f| f.kind.to_string().eq_ignore_ascii_case("gene"))
        );
    }

    #[test]
    fn execute_genomes_verify_anchor_updates_verification_status() {
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
                annotation_scope: None,
                max_annotation_features: None,
                include_genomic_annotation: None,
                catalog_path: Some(catalog_path.clone()),
                cache_dir: None,
            },
        )
        .expect("extract region");

        // Mutate sequence to force a mismatch and validate unverified status recording.
        engine.state_mut().sequences.insert(
            "slice".to_string(),
            DNAsequence::from_sequence("AAAAAAAA").expect("mutated sequence"),
        );

        let out = execute_shell_command(
            &mut engine,
            &ShellCommand::ReferenceVerifyAnchor {
                helper_mode: false,
                seq_id: "slice".to_string(),
                catalog_path: Some(catalog_path),
                cache_dir: None,
                prepared_genome_id: None,
            },
        )
        .expect("verify-anchor");

        assert!(out.state_changed);
        let changed = out
            .output
            .get("result")
            .and_then(|v| v.get("changed_seq_ids"))
            .and_then(|v| v.as_array())
            .expect("changed_seq_ids");
        assert!(
            changed
                .iter()
                .any(|v| v.as_str().map(|id| id == "slice").unwrap_or(false))
        );
        let anchor_status = engine
            .describe_sequence_genome_anchor("slice")
            .expect("anchor status");
        assert!(
            anchor_status.contains("unverified"),
            "expected unverified status, got: {}",
            anchor_status
        );
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
                    FeatureExpertTarget::SplicingFeature {
                        feature_id: 11,
                        scope: SplicingScopePreset::AllOverlappingBothStrands,
                    }
                );
            }
            other => panic!("unexpected command: {other:?}"),
        }

        let isoform = parse_shell_line("inspect-feature-expert s isoform tp53_isoforms_v1")
            .expect("parse isoform feature target");
        match isoform {
            ShellCommand::InspectFeatureExpert { seq_id, target } => {
                assert_eq!(seq_id, "s");
                assert_eq!(
                    target,
                    FeatureExpertTarget::IsoformArchitecture {
                        panel_id: "tp53_isoforms_v1".to_string()
                    }
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
    fn parse_panels_commands() {
        let import = parse_shell_line(
            "panels import-isoform seq_a assets/panels/tp53_isoforms_v1.json --panel-id tp53_isoforms_v1 --strict",
        )
        .expect("parse panels import-isoform");
        match import {
            ShellCommand::PanelsImportIsoform {
                seq_id,
                panel_path,
                panel_id,
                strict,
            } => {
                assert_eq!(seq_id, "seq_a");
                assert_eq!(panel_path, "assets/panels/tp53_isoforms_v1.json");
                assert_eq!(panel_id.as_deref(), Some("tp53_isoforms_v1"));
                assert!(strict);
            }
            other => panic!("unexpected command: {other:?}"),
        }

        let inspect = parse_shell_line("panels inspect-isoform seq_a tp53_isoforms_v1")
            .expect("parse panels inspect-isoform");
        assert!(matches!(inspect, ShellCommand::PanelsInspectIsoform { .. }));

        let render =
            parse_shell_line("panels render-isoform-svg seq_a tp53_isoforms_v1 tp53.panel.svg")
                .expect("parse panels render-isoform-svg");
        assert!(matches!(
            render,
            ShellCommand::PanelsRenderIsoformSvg { .. }
        ));

        let validate = parse_shell_line(
            "panels validate-isoform assets/panels/tp53_isoforms_v1.json --panel-id tp53_isoforms_v1",
        )
        .expect("parse panels validate-isoform");
        assert!(matches!(
            validate,
            ShellCommand::PanelsValidateIsoform { .. }
        ));
    }

    #[test]
    fn parse_uniprot_commands() {
        let genbank = parse_shell_line("genbank fetch NC_000001 --as-id tp73_fetch")
            .expect("parse genbank fetch");
        match genbank {
            ShellCommand::GenbankFetch { accession, as_id } => {
                assert_eq!(accession, "NC_000001");
                assert_eq!(as_id.as_deref(), Some("tp73_fetch"));
            }
            other => panic!("unexpected command: {other:?}"),
        }

        let fetch = parse_shell_line("uniprot fetch P04637 --entry-id TP53_UNIPROT")
            .expect("parse uniprot fetch");
        match fetch {
            ShellCommand::UniprotFetch { query, entry_id } => {
                assert_eq!(query, "P04637");
                assert_eq!(entry_id.as_deref(), Some("TP53_UNIPROT"));
            }
            other => panic!("unexpected command: {other:?}"),
        }

        let import = parse_shell_line("uniprot import-swissprot tp53.txt --entry-id TP53_FILE")
            .expect("parse uniprot import-swissprot");
        match import {
            ShellCommand::UniprotImportSwissProt { path, entry_id } => {
                assert_eq!(path, "tp53.txt");
                assert_eq!(entry_id.as_deref(), Some("TP53_FILE"));
            }
            other => panic!("unexpected command: {other:?}"),
        }

        let map = parse_shell_line(
            "uniprot map P04637 seq_a --projection-id tp53_map --transcript ENST00000269305.9",
        )
        .expect("parse uniprot map");
        match map {
            ShellCommand::UniprotMap {
                entry_id,
                seq_id,
                projection_id,
                transcript_id,
            } => {
                assert_eq!(entry_id, "P04637");
                assert_eq!(seq_id, "seq_a");
                assert_eq!(projection_id.as_deref(), Some("tp53_map"));
                assert_eq!(transcript_id.as_deref(), Some("ENST00000269305.9"));
            }
            other => panic!("unexpected command: {other:?}"),
        }

        let projections =
            parse_shell_line("uniprot projection-list --seq seq_a").expect("parse projection-list");
        assert!(matches!(
            projections,
            ShellCommand::UniprotProjectionList { .. }
        ));
    }

    fn tp53_isoform_test_sequence() -> DNAsequence {
        let mut dna = DNAsequence::from_sequence(&"ACGT".repeat(800)).expect("valid dna");
        dna.features_mut().push(Feature {
            kind: FeatureKind::from("mRNA"),
            location: Location::simple_range(99, 560),
            qualifiers: vec![
                ("gene".into(), Some("TP53".to_string())),
                (
                    "transcript_id".into(),
                    Some("ENST00000269305.9".to_string()),
                ),
                ("label".into(), Some("TP53-201".to_string())),
            ]
            .into_iter()
            .collect(),
        });
        dna
    }

    fn uniprot_projection_test_sequence() -> DNAsequence {
        let mut dna = DNAsequence::from_sequence(&"ACGT".repeat(300)).expect("valid dna");
        dna.features_mut().push(Feature {
            kind: FeatureKind::from("mRNA"),
            location: Location::simple_range(99, 360),
            qualifiers: vec![
                ("gene".into(), Some("TOY1".to_string())),
                ("transcript_id".into(), Some("TX1".to_string())),
                ("label".into(), Some("TX1".to_string())),
                (
                    "cds_ranges_1based".into(),
                    Some("100-180,300-360".to_string()),
                ),
            ]
            .into_iter()
            .collect(),
        });
        dna
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
    fn execute_panels_import_and_inspect_isoform_architecture() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("seq_a".to_string(), tp53_isoform_test_sequence());
        let mut engine = GentleEngine::from_state(state);

        let import = execute_shell_command(
            &mut engine,
            &ShellCommand::PanelsImportIsoform {
                seq_id: "seq_a".to_string(),
                panel_path: "assets/panels/tp53_isoforms_v1.json".to_string(),
                panel_id: Some("tp53_isoforms_v1".to_string()),
                strict: false,
            },
        )
        .expect("execute panels import");
        assert!(import.state_changed);

        let inspect = execute_shell_command(
            &mut engine,
            &ShellCommand::PanelsInspectIsoform {
                seq_id: "seq_a".to_string(),
                panel_id: "tp53_isoforms_v1".to_string(),
            },
        )
        .expect("execute panels inspect");
        assert!(!inspect.state_changed);
        assert_eq!(
            inspect.output["kind"].as_str(),
            Some("isoform_architecture")
        );
        assert_eq!(
            inspect.output["data"]["panel_id"].as_str(),
            Some("tp53_isoforms_v1")
        );
    }

    #[test]
    fn execute_uniprot_import_list_and_projection() {
        let td = tempdir().expect("tempdir");
        let swiss_path = td.path().join("toy_uniprot.txt");
        let swiss_text = r#"ID   TOY1_HUMAN              Reviewed;         30 AA.
AC   PTEST1;
DE   RecName: Full=Toy DNA-binding protein;
GN   Name=TOY1;
OS   Homo sapiens (Human).
DR   Ensembl; TX1; ENSPTOY1; ENSGTOY1.
FT   DOMAIN          2..8
FT                   /note="toy domain"
SQ   SEQUENCE   30 AA;  3333 MW;  0000000000000000 CRC64;
     MEEPQSDPSV EPPLSQETFSDLWKLLPEN
//
"#;
        fs::write(&swiss_path, swiss_text).expect("write swiss file");

        let mut state = ProjectState::default();
        state
            .sequences
            .insert("seq_u".to_string(), uniprot_projection_test_sequence());
        let mut engine = GentleEngine::from_state(state);

        let import = execute_shell_command(
            &mut engine,
            &ShellCommand::UniprotImportSwissProt {
                path: swiss_path.to_string_lossy().to_string(),
                entry_id: None,
            },
        )
        .expect("execute uniprot import");
        assert!(import.state_changed);

        let listed = execute_shell_command(&mut engine, &ShellCommand::UniprotList)
            .expect("execute uniprot list");
        assert!(!listed.state_changed);
        let rows = listed
            .output
            .as_array()
            .expect("uniprot list output must be array");
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0]["entry_id"].as_str(), Some("PTEST1"));

        let mapped = execute_shell_command(
            &mut engine,
            &ShellCommand::UniprotMap {
                entry_id: "PTEST1".to_string(),
                seq_id: "seq_u".to_string(),
                projection_id: None,
                transcript_id: None,
            },
        )
        .expect("execute uniprot map");
        assert!(mapped.state_changed);

        let projection_id = "PTEST1@seq_u".to_string();
        let projection = execute_shell_command(
            &mut engine,
            &ShellCommand::UniprotProjectionShow {
                projection_id: projection_id.clone(),
            },
        )
        .expect("execute uniprot projection-show");
        assert!(!projection.state_changed);
        assert_eq!(
            projection.output["projection_id"].as_str(),
            Some(projection_id.as_str())
        );
        assert_eq!(projection.output["entry_id"].as_str(), Some("PTEST1"));
        assert_eq!(projection.output["seq_id"].as_str(), Some("seq_u"));
        assert!(
            projection.output["transcript_projections"]
                .as_array()
                .map(|v| !v.is_empty())
                .unwrap_or(false)
        );
    }

    #[test]
    fn execute_panels_validate_isoform_returns_report() {
        let mut engine = GentleEngine::new();
        let out = execute_shell_command(
            &mut engine,
            &ShellCommand::PanelsValidateIsoform {
                panel_path: "assets/panels/tp53_isoforms_v1.json".to_string(),
                panel_id: Some("tp53_isoforms_v1".to_string()),
            },
        )
        .expect("validate panel");
        assert!(!out.state_changed);
        assert_eq!(
            out.output["schema"].as_str(),
            Some("gentle.isoform_panel_validation_report.v1")
        );
        assert_eq!(out.output["panel_id"].as_str(), Some("tp53_isoforms_v1"));
        assert!(matches!(
            out.output["status"].as_str(),
            Some("ok") | Some("warning")
        ));
        assert!(out.output["isoform_count"].as_u64().unwrap_or(0) >= 1);
    }

    #[test]
    fn execute_isoform_svg_routes_are_byte_identical() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("seq_a".to_string(), tp53_isoform_test_sequence());
        let mut engine = GentleEngine::from_state(state);
        engine
            .apply(Operation::ImportIsoformPanel {
                seq_id: "seq_a".to_string(),
                panel_path: "assets/panels/tp53_isoforms_v1.json".to_string(),
                panel_id: Some("tp53_isoforms_v1".to_string()),
                strict: false,
            })
            .expect("import panel");

        let tmp = tempdir().expect("temp dir");
        let op_svg = tmp.path().join("isoform.op.svg");
        let shell_svg = tmp.path().join("isoform.shell.svg");
        let expert_svg = tmp.path().join("isoform.expert.svg");
        let op_path = op_svg.display().to_string();
        let shell_path = shell_svg.display().to_string();
        let expert_path = expert_svg.display().to_string();

        engine
            .apply(Operation::RenderIsoformArchitectureSvg {
                seq_id: "seq_a".to_string(),
                panel_id: "tp53_isoforms_v1".to_string(),
                path: op_path.clone(),
            })
            .expect("render op route");

        execute_shell_command(
            &mut engine,
            &ShellCommand::PanelsRenderIsoformSvg {
                seq_id: "seq_a".to_string(),
                panel_id: "tp53_isoforms_v1".to_string(),
                output: shell_path.clone(),
            },
        )
        .expect("render shell route");

        execute_shell_command(
            &mut engine,
            &ShellCommand::RenderFeatureExpertSvg {
                seq_id: "seq_a".to_string(),
                target: FeatureExpertTarget::IsoformArchitecture {
                    panel_id: "tp53_isoforms_v1".to_string(),
                },
                output: expert_path.clone(),
            },
        )
        .expect("render expert route");

        let op_text = fs::read_to_string(op_path).expect("read op svg");
        let shell_text = fs::read_to_string(shell_path).expect("read shell svg");
        let expert_text = fs::read_to_string(expert_path).expect("read expert svg");
        assert_eq!(
            op_text, shell_text,
            "RenderIsoformArchitectureSvg and panels render-isoform-svg outputs must match"
        );
        assert_eq!(
            op_text, expert_text,
            "RenderIsoformArchitectureSvg and render-feature-expert-svg isoform outputs must match"
        );
    }

    #[test]
    fn execute_splicing_feature_expert_svg_shell_and_operation_routes_match() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("seq_a".to_string(), tp53_isoform_test_sequence());
        let mut engine = GentleEngine::from_state(state);
        let feature_id = engine
            .state()
            .sequences
            .get("seq_a")
            .expect("sequence present")
            .features()
            .iter()
            .position(|feature| feature.kind.to_string().eq_ignore_ascii_case("mRNA"))
            .expect("mRNA feature id");

        let tmp = tempdir().expect("temp dir");
        let op_svg = tmp.path().join("splicing.op.svg");
        let shell_svg = tmp.path().join("splicing.shell.svg");
        let op_path = op_svg.display().to_string();
        let shell_path = shell_svg.display().to_string();

        engine
            .apply(Operation::RenderFeatureExpertSvg {
                seq_id: "seq_a".to_string(),
                target: FeatureExpertTarget::SplicingFeature {
                    feature_id,
                    scope: SplicingScopePreset::AllOverlappingBothStrands,
                },
                path: op_path.clone(),
            })
            .expect("render splicing op route");

        execute_shell_command(
            &mut engine,
            &ShellCommand::RenderFeatureExpertSvg {
                seq_id: "seq_a".to_string(),
                target: FeatureExpertTarget::SplicingFeature {
                    feature_id,
                    scope: SplicingScopePreset::AllOverlappingBothStrands,
                },
                output: shell_path.clone(),
            },
        )
        .expect("render splicing shell route");

        let op_text = fs::read_to_string(op_path).expect("read op svg");
        let shell_text = fs::read_to_string(shell_path).expect("read shell svg");
        assert_eq!(
            op_text, shell_text,
            "RenderFeatureExpertSvg operation and shell routes must match for splicing"
        );
        assert!(op_text.contains("cell color intensity encodes exon support frequency"));
    }

    #[test]
    fn parse_dotplot_and_flex_commands() {
        let dotplot = parse_shell_line(
            "dotplot compute seq_a --start 10 --end 110 --mode self_reverse_complement --word-size 8 --step 4 --max-mismatches 1 --tile-bp 256 --id promoter_dp",
        )
        .expect("parse dotplot compute");
        match dotplot {
            ShellCommand::DotplotCompute {
                seq_id,
                span_start_0based,
                span_end_0based,
                mode,
                word_size,
                step_bp,
                max_mismatches,
                tile_bp,
                dotplot_id,
            } => {
                assert_eq!(seq_id, "seq_a");
                assert_eq!(span_start_0based, Some(10));
                assert_eq!(span_end_0based, Some(110));
                assert_eq!(mode, DotplotMode::SelfReverseComplement);
                assert_eq!(word_size, 8);
                assert_eq!(step_bp, 4);
                assert_eq!(max_mismatches, 1);
                assert_eq!(tile_bp, Some(256));
                assert_eq!(dotplot_id.as_deref(), Some("promoter_dp"));
            }
            other => panic!("expected DotplotCompute, got {other:?}"),
        }

        let flex = parse_shell_line(
            "flex compute seq_a --start 25 --end 325 --model at_skew --bin-bp 20 --smoothing-bp 60 --id promoter_flex",
        )
        .expect("parse flex compute");
        match flex {
            ShellCommand::FlexCompute {
                seq_id,
                span_start_0based,
                span_end_0based,
                model,
                bin_bp,
                smoothing_bp,
                track_id,
            } => {
                assert_eq!(seq_id, "seq_a");
                assert_eq!(span_start_0based, Some(25));
                assert_eq!(span_end_0based, Some(325));
                assert_eq!(model, FlexibilityModel::AtSkew);
                assert_eq!(bin_bp, 20);
                assert_eq!(smoothing_bp, Some(60));
                assert_eq!(track_id.as_deref(), Some("promoter_flex"));
            }
            other => panic!("expected FlexCompute, got {other:?}"),
        }
    }

    #[test]
    fn parse_rna_reads_commands() {
        let interpret = parse_shell_line(
            "rna-reads interpret seq_a 7 reads.fa --report-id tp73_reads --scope target_group_any_strand --kmer-len 9 --short-max-bp 420 --long-window-bp 140 --long-window-count 3 --min-seed-hit-fraction 0.30 --min-weighted-seed-hit-fraction 0.05 --min-unique-matched-kmers 12 --min-chain-consistency-fraction 0.55 --max-median-transcript-gap 4.0 --min-confirmed-transitions 1 --min-transition-support-fraction 0.20 --no-cdna-poly-t-flip --poly-t-prefix-min-bp 20 --align-band-bp 24 --align-min-identity 0.60 --max-secondary-mappings 2",
        )
        .expect("parse rna-reads interpret");
        match interpret {
            ShellCommand::RnaReadsInterpret {
                seq_id,
                seed_feature_id,
                input_path,
                scope,
                origin_mode,
                target_gene_ids,
                roi_seed_capture_enabled,
                seed_filter,
                align_config,
                report_id,
                report_mode,
                checkpoint_path,
                checkpoint_every_reads,
                resume_from_checkpoint,
                ..
            } => {
                assert_eq!(seq_id, "seq_a");
                assert_eq!(seed_feature_id, 7);
                assert_eq!(input_path, "reads.fa");
                assert_eq!(scope, SplicingScopePreset::TargetGroupAnyStrand);
                assert_eq!(origin_mode, RnaReadOriginMode::SingleGene);
                assert!(target_gene_ids.is_empty());
                assert!(!roi_seed_capture_enabled);
                assert_eq!(seed_filter.kmer_len, 9);
                assert_eq!(seed_filter.short_full_hash_max_bp, 420);
                assert_eq!(seed_filter.long_window_bp, 140);
                assert_eq!(seed_filter.long_window_count, 3);
                assert!((seed_filter.min_seed_hit_fraction - 0.30).abs() < f64::EPSILON);
                assert!((seed_filter.min_weighted_seed_hit_fraction - 0.05).abs() < f64::EPSILON);
                assert_eq!(seed_filter.min_unique_matched_kmers, 12);
                assert!((seed_filter.min_chain_consistency_fraction - 0.55).abs() < f64::EPSILON);
                assert!((seed_filter.max_median_transcript_gap - 4.0).abs() < f64::EPSILON);
                assert_eq!(seed_filter.min_confirmed_exon_transitions, 1);
                assert!((seed_filter.min_transition_support_fraction - 0.20).abs() < f64::EPSILON);
                assert!(!seed_filter.cdna_poly_t_flip_enabled);
                assert_eq!(seed_filter.poly_t_prefix_min_bp, 20);
                assert_eq!(align_config.band_width_bp, 24);
                assert!((align_config.min_identity_fraction - 0.60).abs() < f64::EPSILON);
                assert_eq!(align_config.max_secondary_mappings, 2);
                assert_eq!(report_id.as_deref(), Some("tp73_reads"));
                assert_eq!(report_mode, RnaReadReportMode::Full);
                assert!(checkpoint_path.is_none());
                assert_eq!(checkpoint_every_reads, 10_000);
                assert!(!resume_from_checkpoint);
            }
            other => panic!("expected RnaReadsInterpret, got {other:?}"),
        }

        let interpret_multi = parse_shell_line(
            "rna-reads interpret seq_a 7 reads.fa --origin-mode multi_gene_sparse --target-gene TP73 --target-gene TP53 --roi-seed-capture",
        )
        .expect("parse rna-reads interpret multi-gene scaffold");
        match interpret_multi {
            ShellCommand::RnaReadsInterpret {
                origin_mode,
                target_gene_ids,
                roi_seed_capture_enabled,
                ..
            } => {
                assert_eq!(origin_mode, RnaReadOriginMode::MultiGeneSparse);
                assert_eq!(
                    target_gene_ids,
                    vec!["TP73".to_string(), "TP53".to_string()]
                );
                assert!(roi_seed_capture_enabled);
            }
            other => panic!("expected RnaReadsInterpret, got {other:?}"),
        }

        let interpret_checkpoint = parse_shell_line(
            "rna-reads interpret seq_a 7 reads.fa --report-mode seed_passed_only --checkpoint-path /tmp/tp53.chk.json --checkpoint-every-reads 2500 --resume-from-checkpoint",
        )
        .expect("parse rna-reads interpret checkpoint options");
        match interpret_checkpoint {
            ShellCommand::RnaReadsInterpret {
                report_mode,
                checkpoint_path,
                checkpoint_every_reads,
                resume_from_checkpoint,
                ..
            } => {
                assert_eq!(report_mode, RnaReadReportMode::SeedPassedOnly);
                assert_eq!(
                    checkpoint_path.as_deref(),
                    Some("/tmp/tp53.chk.json")
                );
                assert_eq!(checkpoint_every_reads, 2500);
                assert!(resume_from_checkpoint);
            }
            other => panic!("expected RnaReadsInterpret, got {other:?}"),
        }

        let list =
            parse_shell_line("rna-reads list-reports seq_a").expect("parse rna-reads list-reports");
        assert!(matches!(
            list,
            ShellCommand::RnaReadsListReports { seq_id } if seq_id.as_deref() == Some("seq_a")
        ));

        let show =
            parse_shell_line("rna-reads show-report tp73_reads").expect("parse rna-reads show");
        assert!(matches!(
            show,
            ShellCommand::RnaReadsShowReport { report_id } if report_id == "tp73_reads"
        ));

        let export = parse_shell_line("rna-reads export-report tp73_reads out.json")
            .expect("parse rna-reads export-report");
        assert!(matches!(
            export,
            ShellCommand::RnaReadsExportReport { report_id, path } if report_id == "tp73_reads" && path == "out.json"
        ));

        let export_hits =
            parse_shell_line("rna-reads export-hits-fasta tp73_reads hits.fa --selection aligned")
                .expect("parse rna-reads export-hits-fasta");
        assert!(matches!(
            export_hits,
            ShellCommand::RnaReadsExportHitsFasta { report_id, path, selection } if report_id == "tp73_reads" && path == "hits.fa" && selection == RnaReadHitSelection::Aligned
        ));

        let export_sheet = parse_shell_line(
            "rna-reads export-sample-sheet samples.tsv --seq-id seq_a --report-id tp73_reads --append",
        )
        .expect("parse rna-reads export-sample-sheet");
        assert!(matches!(
            export_sheet,
            ShellCommand::RnaReadsExportSampleSheet { path, seq_id, report_ids, append }
                if path == "samples.tsv"
                    && seq_id.as_deref() == Some("seq_a")
                    && report_ids == vec!["tp73_reads".to_string()]
                    && append
        ));

        let export_paths = parse_shell_line(
            "rna-reads export-paths-tsv tp73_reads paths.tsv --selection seed_passed",
        )
        .expect("parse rna-reads export-paths-tsv");
        assert!(matches!(
            export_paths,
            ShellCommand::RnaReadsExportExonPathsTsv { report_id, path, selection }
                if report_id == "tp73_reads"
                    && path == "paths.tsv"
                    && selection == RnaReadHitSelection::SeedPassed
        ));

        let export_abundance = parse_shell_line(
            "rna-reads export-abundance-tsv tp73_reads abundance.tsv --selection aligned",
        )
        .expect("parse rna-reads export-abundance-tsv");
        assert!(matches!(
            export_abundance,
            ShellCommand::RnaReadsExportExonAbundanceTsv { report_id, path, selection }
                if report_id == "tp73_reads"
                    && path == "abundance.tsv"
                    && selection == RnaReadHitSelection::Aligned
        ));

        let export_score_density = parse_shell_line(
            "rna-reads export-score-density-svg tp73_reads score_density.svg --scale linear",
        )
        .expect("parse rna-reads export-score-density-svg");
        assert!(matches!(
            export_score_density,
            ShellCommand::RnaReadsExportScoreDensitySvg { report_id, path, scale }
                if report_id == "tp73_reads"
                    && path == "score_density.svg"
                    && scale == RnaReadScoreDensityScale::Linear
        ));
    }

    #[test]
    fn parse_rna_reads_interpret_defaults_to_engine_kmer_len() {
        let cmd = parse_shell_line(
            "rna-reads interpret seq_a 7 reads.fa --scope all_overlapping_both_strands",
        )
        .expect("parse rna-reads interpret defaults");
        match cmd {
            ShellCommand::RnaReadsInterpret {
                origin_mode,
                target_gene_ids,
                roi_seed_capture_enabled,
                seed_filter,
                report_mode,
                checkpoint_path,
                checkpoint_every_reads,
                resume_from_checkpoint,
                ..
            } => {
                assert_eq!(origin_mode, RnaReadOriginMode::SingleGene);
                assert!(target_gene_ids.is_empty());
                assert!(!roi_seed_capture_enabled);
                assert_eq!(
                    seed_filter.kmer_len,
                    RnaReadSeedFilterConfig::default().kmer_len
                );
                assert_eq!(seed_filter.kmer_len, 10);
                assert_eq!(report_mode, RnaReadReportMode::Full);
                assert!(checkpoint_path.is_none());
                assert_eq!(checkpoint_every_reads, 10_000);
                assert!(!resume_from_checkpoint);
            }
            other => panic!("expected RnaReadsInterpret, got {other:?}"),
        }
    }

    #[test]
    fn execute_dotplot_and_flex_commands_store_payloads() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "seq_a".to_string(),
            DNAsequence::from_sequence("ATGCATGCATGCATGCATGCATGC").expect("sequence"),
        );
        let mut engine = GentleEngine::from_state(state);

        let dotplot = execute_shell_command(
            &mut engine,
            &ShellCommand::DotplotCompute {
                seq_id: "seq_a".to_string(),
                span_start_0based: Some(0),
                span_end_0based: Some(24),
                mode: DotplotMode::SelfForward,
                word_size: 4,
                step_bp: 2,
                max_mismatches: 0,
                tile_bp: None,
                dotplot_id: Some("dp_1".to_string()),
            },
        )
        .expect("execute dotplot compute");
        assert!(dotplot.state_changed);
        assert_eq!(
            dotplot.output["dotplot"]["dotplot_id"].as_str(),
            Some("dp_1")
        );

        let flex = execute_shell_command(
            &mut engine,
            &ShellCommand::FlexCompute {
                seq_id: "seq_a".to_string(),
                span_start_0based: Some(0),
                span_end_0based: Some(24),
                model: FlexibilityModel::AtRichness,
                bin_bp: 6,
                smoothing_bp: Some(12),
                track_id: Some("flex_1".to_string()),
            },
        )
        .expect("execute flex compute");
        assert!(flex.state_changed);
        assert_eq!(flex.output["track"]["track_id"].as_str(), Some("flex_1"));

        let listed = execute_shell_command(
            &mut engine,
            &ShellCommand::DotplotList {
                seq_id: Some("seq_a".to_string()),
            },
        )
        .expect("list dotplots");
        assert_eq!(listed.output["dotplot_count"].as_u64(), Some(1));
        let flex_list = execute_shell_command(
            &mut engine,
            &ShellCommand::FlexList {
                seq_id: Some("seq_a".to_string()),
            },
        )
        .expect("list flex tracks");
        assert_eq!(flex_list.output["track_count"].as_u64(), Some(1));
    }

    #[test]
    fn execute_rna_reads_commands_store_and_export_reports() {
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("seq_a".to_string(), tp53_isoform_test_sequence());
        let mut engine = GentleEngine::from_state(state);
        let feature_id = engine
            .state()
            .sequences
            .get("seq_a")
            .expect("sequence present")
            .features()
            .iter()
            .position(|feature| feature.kind.to_string().eq_ignore_ascii_case("mRNA"))
            .expect("mRNA feature id");
        let fasta_dir = tempdir().expect("tempdir");
        let input_path = fasta_dir.path().join("reads.fa");
        fs::write(
            &input_path,
            ">read_1\nATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTTCCTAATGGGCCCGGATTCCTTTTCTCTGTGAACCTTCCCGATGATGATGGAGGTGGAATGGAGGAGCCGCAGTCA\n",
        )
        .expect("write input fasta");
        let report_id = "rna_reads_test".to_string();
        let run = execute_shell_command(
            &mut engine,
            &ShellCommand::RnaReadsInterpret {
                seq_id: "seq_a".to_string(),
                seed_feature_id: feature_id,
                input_path: input_path.display().to_string(),
                profile: RnaReadInterpretationProfile::NanoporeCdnaV1,
                input_format: RnaReadInputFormat::Fasta,
                scope: SplicingScopePreset::AllOverlappingBothStrands,
                origin_mode: RnaReadOriginMode::SingleGene,
                target_gene_ids: vec![],
                roi_seed_capture_enabled: false,
                seed_filter: RnaReadSeedFilterConfig::default(),
                align_config: RnaReadAlignConfig::default(),
                report_id: Some(report_id.clone()),
                report_mode: RnaReadReportMode::Full,
                checkpoint_path: None,
                checkpoint_every_reads: 10_000,
                resume_from_checkpoint: false,
            },
        )
        .expect("execute rna-reads interpret");
        assert!(run.state_changed);
        assert_eq!(
            run.output["report"]["report_id"].as_str(),
            Some(report_id.as_str())
        );
        assert_eq!(run.output["report"]["read_count_total"].as_u64(), Some(1));

        let listed = execute_shell_command(
            &mut engine,
            &ShellCommand::RnaReadsListReports {
                seq_id: Some("seq_a".to_string()),
            },
        )
        .expect("list rna-read reports");
        assert_eq!(listed.output["report_count"].as_u64(), Some(1));

        let shown = execute_shell_command(
            &mut engine,
            &ShellCommand::RnaReadsShowReport {
                report_id: report_id.clone(),
            },
        )
        .expect("show rna-read report");
        assert_eq!(
            shown.output["report"]["report_id"].as_str(),
            Some(report_id.as_str())
        );

        let exported_report = fasta_dir.path().join("report.json");
        let export_result = execute_shell_command(
            &mut engine,
            &ShellCommand::RnaReadsExportReport {
                report_id: report_id.clone(),
                path: exported_report.display().to_string(),
            },
        )
        .expect("export rna-read report");
        assert_eq!(
            export_result.output["report_id"].as_str(),
            Some(report_id.as_str())
        );
        assert!(exported_report.exists());

        let exported_hits = fasta_dir.path().join("hits.fa");
        let export_hits_result = execute_shell_command(
            &mut engine,
            &ShellCommand::RnaReadsExportHitsFasta {
                report_id,
                path: exported_hits.display().to_string(),
                selection: RnaReadHitSelection::All,
            },
        )
        .expect("export rna-read hits");
        assert_eq!(
            export_hits_result.output["written_records"].as_u64(),
            Some(1)
        );
        let fasta_text = fs::read_to_string(exported_hits).expect("read exported hits");
        assert!(fasta_text.contains(">read_1"));

        let exported_sheet = fasta_dir.path().join("samples.tsv");
        let export_sheet_result = execute_shell_command(
            &mut engine,
            &ShellCommand::RnaReadsExportSampleSheet {
                path: exported_sheet.display().to_string(),
                seq_id: Some("seq_a".to_string()),
                report_ids: vec![],
                append: false,
            },
        )
        .expect("export rna-read sample sheet");
        assert_eq!(export_sheet_result.output["report_count"].as_u64(), Some(1));
        let sheet_text = fs::read_to_string(exported_sheet).expect("read sample sheet");
        assert!(sheet_text.contains("sample_id"));
        assert!(sheet_text.contains("exon_support_frequencies_json"));

        let exported_paths = fasta_dir.path().join("paths.tsv");
        let export_paths_result = execute_shell_command(
            &mut engine,
            &ShellCommand::RnaReadsExportExonPathsTsv {
                report_id: "rna_reads_test".to_string(),
                path: exported_paths.display().to_string(),
                selection: RnaReadHitSelection::All,
            },
        )
        .expect("export rna-read exon paths");
        assert_eq!(export_paths_result.output["row_count"].as_u64(), Some(1));
        let paths_text = fs::read_to_string(exported_paths).expect("read path sheet");
        assert!(paths_text.contains("exon_path"));
        assert!(paths_text.contains("reverse_complement_applied"));

        let exported_abundance = fasta_dir.path().join("abundance.tsv");
        let export_abundance_result = execute_shell_command(
            &mut engine,
            &ShellCommand::RnaReadsExportExonAbundanceTsv {
                report_id: "rna_reads_test".to_string(),
                path: exported_abundance.display().to_string(),
                selection: RnaReadHitSelection::All,
            },
        )
        .expect("export rna-read abundance");
        assert_eq!(
            export_abundance_result.output["selected_read_count"].as_u64(),
            Some(1)
        );
        let abundance_text = fs::read_to_string(exported_abundance).expect("read abundance sheet");
        assert!(abundance_text.contains("row_kind"));

        let exported_density_svg = fasta_dir.path().join("score_density.svg");
        let export_density_result = execute_shell_command(
            &mut engine,
            &ShellCommand::RnaReadsExportScoreDensitySvg {
                report_id: "rna_reads_test".to_string(),
                path: exported_density_svg.display().to_string(),
                scale: RnaReadScoreDensityScale::Log,
            },
        )
        .expect("export rna-read score density svg");
        assert_eq!(export_density_result.output["scale"].as_str(), Some("log"));
        let density_text = fs::read_to_string(exported_density_svg).expect("read density svg");
        assert!(density_text.contains("<svg"));
        assert!(density_text.contains("seed-hit score density"));
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
