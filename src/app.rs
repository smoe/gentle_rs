use std::{
    collections::hash_map::DefaultHasher,
    collections::{BTreeMap, HashMap, HashSet},
    env, fs,
    hash::{Hash, Hasher},
    io::ErrorKind,
    panic::{AssertUnwindSafe, catch_unwind},
    path::{Path, PathBuf},
    process::Command,
    sync::{
        Arc, RwLock,
        atomic::{AtomicBool, AtomicU64, Ordering},
        mpsc,
    },
    time::{Duration, Instant, SystemTime, UNIX_EPOCH},
};

use crate::{
    TRANSLATIONS, about,
    agent_bridge::{
        AGENT_BASE_URL_ENV, AGENT_CONNECT_TIMEOUT_SECS_ENV, AGENT_MAX_RESPONSE_BYTES_ENV,
        AGENT_MAX_RETRIES_ENV, AGENT_MODEL_ENV, AGENT_READ_TIMEOUT_SECS_ENV,
        AGENT_TIMEOUT_SECS_ENV, AgentExecutionIntent, AgentInvocationOutcome, AgentResponse,
        AgentSystemSpec, AgentSystemTransport, DEFAULT_AGENT_SYSTEM_CATALOG_PATH,
        OPENAI_API_KEY_ENV, OPENAI_COMPAT_UNSPECIFIED_MODEL, agent_system_availability,
        discover_openai_models, invoke_agent_support_with_env_overrides, load_agent_system_catalog,
    },
    dna_sequence::{self, DNAsequence},
    engine::{
        BIGWIG_TO_BEDGRAPH_ENV_BIN, BlastHitFeatureInput, DEFAULT_BIGWIG_TO_BEDGRAPH_BIN,
        DisplaySettings, DisplayTarget, Engine, EngineError, ErrorCode, GenomeTrackImportProgress,
        GenomeTrackSource, GenomeTrackSubscription, GenomeTrackSyncReport, GentleEngine, OpResult,
        Operation, OperationProgress, ProjectState, SequenceGenomeAnchorSummary,
    },
    engine_shell::{
        ShellCommand, ShellExecutionOptions, UiIntentTarget, execute_shell_command_with_options,
        parse_shell_line,
    },
    enzymes,
    genomes::{
        BLASTN_ENV_BIN, DEFAULT_BLASTN_BIN, DEFAULT_GENOME_CACHE_DIR, DEFAULT_GENOME_CATALOG_PATH,
        DEFAULT_HELPER_GENOME_CATALOG_PATH, DEFAULT_MAKEBLASTDB_BIN, GenomeBlastReport,
        GenomeCatalog, GenomeChromosomeRecord, GenomeGeneRecord, GenomeSourcePlan,
        MAKEBLASTDB_ENV_BIN, PREPARE_GENOME_TIMEOUT_SECS_ENV, PrepareGenomeProgress,
        PreparedGenomeInspection,
    },
    icons::APP_ICON,
    resource_sync,
    shell_docs::{
        shell_help_markdown as render_shell_help_markdown,
        shell_help_markdown_from_glossary_json as render_shell_help_markdown_from_glossary_json,
    },
    tf_motifs, tool_overrides,
    window::Window,
    window_backdrop::{self, WindowBackdropKind, WindowBackdropSettings},
};
use anyhow::{Result, anyhow};
use eframe::egui::{self, Key, KeyboardShortcut, Modifiers, Pos2, Ui, Vec2, ViewportId, menu};
use egui_commonmark::{CommonMarkCache, CommonMarkViewer};
use pulldown_cmark::{Event, LinkType, Parser, Tag};
use regex::{Regex, RegexBuilder};
use serde::{Deserialize, Serialize};

const GUI_MANUAL_MD: &str = include_str!("../docs/gui.md");
const CLI_MANUAL_MD: &str = include_str!("../docs/cli.md");
const APP_CONFIGURATION_FILE_NAME: &str = ".gentle_gui_settings.json";
const MAX_RECENT_PROJECTS: usize = 12;
const LINEAGE_GRAPH_WORKSPACE_METADATA_KEY: &str = "gui.lineage_graph.workspace";
const LINEAGE_NODE_OFFSETS_METADATA_KEY: &str = "gui.lineage_graph.node_offsets";
const GUI_OPENAI_DEFAULT_BASE_URL: &str = "https://api.openai.com/v1";
const GUI_OPENAI_COMPAT_DEFAULT_BASE_URL: &str = "http://127.0.0.1:11434/v1";
static NATIVE_HELP_OPEN_REQUESTED: AtomicBool = AtomicBool::new(false);
static NATIVE_SETTINGS_OPEN_REQUESTED: AtomicBool = AtomicBool::new(false);
static NATIVE_WINDOWS_OPEN_REQUESTED: AtomicBool = AtomicBool::new(false);
const NATIVE_WINDOW_FOCUS_KEY_NONE: u64 = u64::MAX;
static NATIVE_WINDOWS_FOCUS_KEY_REQUESTED: AtomicU64 = AtomicU64::new(NATIVE_WINDOW_FOCUS_KEY_NONE);
static ACTIVE_VIEWPORT_KEY_REPORTED: AtomicU64 = AtomicU64::new(NATIVE_WINDOW_FOCUS_KEY_NONE);

fn viewport_native_menu_key(viewport_id: ViewportId) -> u64 {
    let mut hasher = DefaultHasher::new();
    viewport_id.hash(&mut hasher);
    hasher.finish() & 0x7fff_ffff_ffff_ffff
}

fn report_active_viewport_from_ui(viewport_id: ViewportId) {
    ACTIVE_VIEWPORT_KEY_REPORTED.store(viewport_native_menu_key(viewport_id), Ordering::SeqCst);
}

fn normalize_agent_model_name(raw: &str) -> Option<String> {
    let trimmed = raw.trim();
    if trimmed.is_empty() || trimmed.eq_ignore_ascii_case(OPENAI_COMPAT_UNSPECIFIED_MODEL) {
        None
    } else {
        Some(trimmed.to_string())
    }
}

const AGENT_PROMPT_TEMPLATE_DEFAULT_ID: &str = "structured";

fn agent_prompt_template_options() -> &'static [(&'static str, &'static str)] {
    &[
        ("structured", "Structured (recommended)"),
        ("candidate_anchors", "Candidate between anchors"),
        ("blast_specificity", "BLAST specificity check"),
        ("track_intersection", "Track import + prioritization"),
        ("macro_template", "Macro/template authoring"),
    ]
}

fn agent_prompt_template_label(id: &str) -> &'static str {
    agent_prompt_template_options()
        .iter()
        .find(|(value, _)| *value == id)
        .map(|(_, label)| *label)
        .unwrap_or("Structured (recommended)")
}

fn agent_prompt_template_text(id: &str) -> &'static str {
    match id {
        "candidate_anchors" => {
            r#"Objective:
Generate candidate windows between two local anchors and rank them.

Context:
Project sequence ID: <SEQ_ID>

Inputs:
- seq_id: <SEQ_ID>
- anchor A: <feature boundary or absolute position>
- anchor B: <feature boundary or absolute position>

Constraints:
- length: 20
- step: 1
- GC range: 40-80%
- additional constraints: <motifs/sites/strand>

Output wanted:
- exact `gentle_cli shell "candidates ..."` commands
- scoring + filter + top-k steps
- validation checklist

Execution policy:
ask-before-run"#
        }
        "blast_specificity" => {
            r#"Objective:
Run a specificity check for one sequence with BLAST.

Context:
Target catalog: genomes | helpers

Inputs:
- genome_id/helper_id: <ID>
- query_sequence: <ACGT...>

Constraints:
- max_hits: 20
- task: blastn-short

Output wanted:
- exact BLAST command
- concise interpretation checklist for top hits

Execution policy:
chat-only"#
        }
        "track_intersection" => {
            r#"Objective:
Import external track evidence and prioritize candidates near track features.

Context:
Anchored sequence ID: <SEQ_ID>

Inputs:
- seq_id: <SEQ_ID>
- track file path: <BED/BED.GZ/BIGWIG/VCF path>

Constraints:
- keep imported features in TRACK groups
- do not modify original sequence content

Output wanted:
- exact track-import commands
- candidate generation/scoring/filter commands near imported TRACK features
- validation checklist

Execution policy:
ask-before-run"#
        }
        "macro_template" => {
            r#"Objective:
Create or update a reusable candidate macro template and run it with bindings.

Context:
Template name: <NAME>

Inputs:
- template parameters: <param list>
- script intent: <generate/score/filter/top-k/...>

Constraints:
- transactional run enabled
- deterministic tie-break policy where applicable

Output wanted:
- `candidates template-put` and `candidates template-run` commands
- brief note on expected outputs and rollback behavior

Execution policy:
ask-before-run"#
        }
        _ => {
            r#"Objective:
<one clear goal>

Context:
<sequence/genome/helper IDs and short background>

Inputs:
- seq_id / genome_id / helper_id: ...
- anchors or coordinates: ...
- feature labels/kinds: ...

Constraints:
- length: ...
- GC range: ...
- motifs/sites to require or avoid: ...
- strand assumptions: ...

Output wanted:
- plan
- exact gentle_cli commands
- validation checklist

Execution policy:
chat-only | ask-before-run | allow-auto-exec"#
        }
    }
}

fn default_prepare_timeout_secs_string() -> String {
    env::var(PREPARE_GENOME_TIMEOUT_SECS_ENV)
        .ok()
        .map(|v| v.trim().to_string())
        .filter(|v| !v.is_empty())
        .unwrap_or_default()
}

fn default_agent_timeout_secs_string() -> String {
    env::var(AGENT_TIMEOUT_SECS_ENV)
        .ok()
        .map(|v| v.trim().to_string())
        .filter(|v| !v.is_empty())
        .unwrap_or_default()
}

fn default_agent_connect_timeout_secs_string() -> String {
    env::var(AGENT_CONNECT_TIMEOUT_SECS_ENV)
        .ok()
        .map(|v| v.trim().to_string())
        .filter(|v| !v.is_empty())
        .unwrap_or_default()
}

fn default_agent_read_timeout_secs_string() -> String {
    env::var(AGENT_READ_TIMEOUT_SECS_ENV)
        .ok()
        .map(|v| v.trim().to_string())
        .filter(|v| !v.is_empty())
        .unwrap_or_default()
}

fn default_agent_max_retries_string() -> String {
    env::var(AGENT_MAX_RETRIES_ENV)
        .ok()
        .map(|v| v.trim().to_string())
        .filter(|v| !v.is_empty())
        .unwrap_or_default()
}

fn default_agent_max_response_bytes_string() -> String {
    env::var(AGENT_MAX_RESPONSE_BYTES_ENV)
        .ok()
        .map(|v| v.trim().to_string())
        .filter(|v| !v.is_empty())
        .unwrap_or_default()
}

pub fn request_open_help_from_native_menu() {
    NATIVE_HELP_OPEN_REQUESTED.store(true, Ordering::SeqCst);
}

pub fn request_open_settings_from_native_menu() {
    NATIVE_SETTINGS_OPEN_REQUESTED.store(true, Ordering::SeqCst);
}

pub fn request_open_windows_from_native_menu() {
    NATIVE_WINDOWS_FOCUS_KEY_REQUESTED.store(NATIVE_WINDOW_FOCUS_KEY_NONE, Ordering::SeqCst);
    NATIVE_WINDOWS_OPEN_REQUESTED.store(true, Ordering::SeqCst);
}

pub fn request_focus_window_key_from_native_menu(key: u64) {
    NATIVE_WINDOWS_FOCUS_KEY_REQUESTED.store(key, Ordering::SeqCst);
    NATIVE_WINDOWS_OPEN_REQUESTED.store(true, Ordering::SeqCst);
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
struct PersistedConfiguration {
    rnapkin_executable: String,
    makeblastdb_executable: String,
    blastn_executable: String,
    bigwig_to_bedgraph_executable: String,
    graphics_defaults: DisplaySettings,
    window_backdrops: WindowBackdropSettings,
    recent_projects: Vec<String>,
}

impl Default for PersistedConfiguration {
    fn default() -> Self {
        Self {
            rnapkin_executable: String::new(),
            makeblastdb_executable: String::new(),
            blastn_executable: String::new(),
            bigwig_to_bedgraph_executable: String::new(),
            graphics_defaults: DisplaySettings::default(),
            window_backdrops: WindowBackdropSettings::default(),
            recent_projects: vec![],
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
enum GenomeTrackSourceSelection {
    #[default]
    Auto,
    Bed,
    BigWig,
    Vcf,
}

impl GenomeTrackSourceSelection {
    fn label(self) -> &'static str {
        match self {
            Self::Auto => "Auto (from extension)",
            Self::Bed => "BED",
            Self::BigWig => "BigWig",
            Self::Vcf => "VCF",
        }
    }

    fn resolve(self, path: &str) -> GenomeTrackSource {
        match self {
            Self::Auto => GenomeTrackSource::from_path(path),
            Self::Bed => GenomeTrackSource::Bed,
            Self::BigWig => GenomeTrackSource::BigWig,
            Self::Vcf => GenomeTrackSource::Vcf,
        }
    }
}

pub struct GENtleApp {
    engine: Arc<RwLock<GentleEngine>>,
    new_windows: Vec<Window>,
    windows: HashMap<ViewportId, Arc<RwLock<Window>>>,
    windows_to_close: Arc<RwLock<Vec<ViewportId>>>,
    pending_focus_viewports: Vec<ViewportId>,
    viewport_id_counter: usize,
    update_has_run_before: bool,
    show_about_dialog: bool,
    show_help_dialog: bool,
    help_doc: HelpDoc,
    help_markdown_cache: CommonMarkCache,
    help_gui_markdown: String,
    help_cli_markdown: String,
    help_shell_markdown: String,
    help_shell_interface: ShellHelpInterface,
    help_search_query: String,
    help_search_matches: Vec<HelpSearchMatch>,
    help_search_selected: usize,
    help_focus_search_box: bool,
    help_image_preview_path: Option<String>,
    help_image_preview_caption: String,
    show_configuration_dialog: bool,
    configuration_tab: ConfigurationTab,
    configuration_rnapkin_executable: String,
    configuration_makeblastdb_executable: String,
    configuration_blastn_executable: String,
    configuration_bigwig_to_bedgraph_executable: String,
    configuration_rnapkin_validation_ok: Option<bool>,
    configuration_rnapkin_validation_message: String,
    configuration_blast_validation_ok: Option<bool>,
    configuration_blast_validation_message: String,
    configuration_graphics: DisplaySettings,
    configuration_graphics_dirty: bool,
    window_backdrops: WindowBackdropSettings,
    configuration_window_backdrops: WindowBackdropSettings,
    configuration_window_backdrops_dirty: bool,
    configuration_status: String,
    current_project_path: Option<String>,
    recent_project_paths: Vec<String>,
    lineage_graph_view: bool,
    lineage_graph_zoom: f32,
    lineage_graph_area_height: f32,
    lineage_container_area_height: f32,
    lineage_graph_scroll_offset: Vec2,
    lineage_graph_pan_origin: Option<Vec2>,
    lineage_graph_compact_labels: bool,
    lineage_graph_selected_node_id: Option<String>,
    lineage_graph_node_offsets: HashMap<String, Vec2>,
    lineage_graph_drag_origin: Option<(String, Vec2)>,
    lineage_graph_offsets_synced_stamp: u64,
    lineage_cache_stamp: u64,
    lineage_cache_valid: bool,
    lineage_rows: Vec<LineageRow>,
    lineage_edges: Vec<(String, String, String)>,
    lineage_op_label_by_id: HashMap<String, String>,
    lineage_containers: Vec<ContainerRow>,
    lineage_arrangements: Vec<ArrangementRow>,
    clean_state_fingerprint: u64,
    dirty_cache_stamp: u64,
    dirty_cache_valid: bool,
    dirty_cache_value: bool,
    dirty_cache_last_deep_check: Instant,
    last_applied_window_title: String,
    last_native_window_entries: Vec<(u64, String)>,
    last_native_active_window_key: Option<u64>,
    native_window_key_to_viewport: HashMap<u64, ViewportId>,
    active_window_menu_key: Option<u64>,
    pending_project_action: Option<ProjectAction>,
    show_reference_genome_prepare_dialog: bool,
    show_reference_genome_retrieve_dialog: bool,
    show_reference_genome_blast_dialog: bool,
    show_reference_genome_inspector_dialog: bool,
    genome_catalog_path: String,
    genome_cache_dir: String,
    genome_id: String,
    genome_catalog_genomes: Vec<String>,
    genome_catalog_error: String,
    genome_prepare_task: Option<GenomePrepareTask>,
    genome_prepare_progress: Option<PrepareGenomeProgress>,
    genome_prepare_status: String,
    genome_prepare_timeout_secs: String,
    genome_gene_filter: String,
    genome_gene_filter_limit: usize,
    genome_gene_filter_page: usize,
    genome_biotype_filter: HashMap<String, bool>,
    genome_genes: Vec<GenomeGeneRecord>,
    genome_genes_loaded_key: Option<String>,
    genome_genes_error: String,
    genome_selected_gene: Option<usize>,
    genome_chromosome: String,
    genome_start_1based: String,
    genome_end_1based: String,
    genome_output_id: String,
    genome_retrieve_status: String,
    genome_blast_source_mode: GenomeBlastSourceMode,
    genome_blast_query_manual: String,
    genome_blast_query_seq_id: String,
    genome_blast_query_pool_id: String,
    genome_blast_max_hits: usize,
    genome_blast_task_name: String,
    genome_blast_task: Option<GenomeBlastTask>,
    genome_blast_progress_fraction: Option<f32>,
    genome_blast_progress_label: String,
    genome_blast_results: Vec<GenomeBlastQueryResult>,
    genome_blast_selected_result: usize,
    genome_blast_status: String,
    show_genome_bed_track_dialog: bool,
    show_agent_assistant_dialog: bool,
    genome_track_seq_id: String,
    genome_track_source_selection: GenomeTrackSourceSelection,
    genome_track_path: String,
    genome_track_name: String,
    genome_track_min_score: String,
    genome_track_max_score: String,
    genome_track_clear_existing: bool,
    genome_track_status: String,
    genome_track_import_task: Option<GenomeTrackImportTask>,
    genome_track_import_progress: Option<GenomeTrackImportProgress>,
    genome_bed_track_subscriptions: Vec<GenomeTrackSubscription>,
    genome_track_autosync_status: String,
    tracked_autosync_last_op_count: Option<usize>,
    genome_blast_import_track_name: String,
    genome_blast_import_clear_existing: bool,
    show_command_palette_dialog: bool,
    command_palette_query: String,
    command_palette_selected: usize,
    command_palette_focus_query: bool,
    show_jobs_panel: bool,
    show_history_panel: bool,
    hover_status_name: String,
    app_status: String,
    next_background_job_id: u64,
    job_event_log: Vec<BackgroundJobEvent>,
    genome_track_preflight_track_subscription: bool,
    agent_catalog_path: String,
    agent_catalog_loaded_path: String,
    agent_catalog_error: String,
    agent_systems: Vec<AgentSystemSpec>,
    agent_system_id: String,
    agent_openai_api_key: String,
    agent_base_url_override: String,
    agent_model_override: String,
    agent_timeout_secs: String,
    agent_connect_timeout_secs: String,
    agent_read_timeout_secs: String,
    agent_max_retries: String,
    agent_max_response_bytes: String,
    agent_discovered_models: Vec<String>,
    agent_discovered_model_pick: String,
    agent_model_discovery_status: String,
    agent_model_discovery_source_key: String,
    agent_model_discovery_task: Option<AgentModelDiscoveryTask>,
    agent_prompt_template_id: String,
    agent_prompt: String,
    agent_include_state_summary: bool,
    agent_allow_auto_exec: bool,
    agent_status: String,
    agent_task: Option<AgentAskTask>,
    agent_last_invocation: Option<AgentInvocationOutcome>,
    agent_execution_log: Vec<AgentCommandExecutionRecord>,
}

#[derive(Clone)]
enum ProjectAction {
    New,
    Open,
    OpenPath(String),
    Close,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum HelpDoc {
    Gui,
    Cli,
    Shell,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum ShellHelpInterface {
    All,
    GuiShell,
    CliShell,
    CliDirect,
    Js,
    Lua,
}

impl ShellHelpInterface {
    fn label(self) -> &'static str {
        match self {
            Self::All => "All",
            Self::GuiShell => "GUI shell",
            Self::CliShell => "CLI shell",
            Self::CliDirect => "CLI direct",
            Self::Js => "JS",
            Self::Lua => "Lua",
        }
    }

    fn glossary_filter(self) -> Option<&'static str> {
        match self {
            Self::All => None,
            Self::GuiShell => Some("gui-shell"),
            Self::CliShell => Some("cli-shell"),
            Self::CliDirect => Some("cli-direct"),
            Self::Js => Some("js"),
            Self::Lua => Some("lua"),
        }
    }
}

#[derive(Clone, Debug)]
struct HelpSearchMatch {
    line_number: usize,
    snippet: String,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum ConfigurationTab {
    ExternalApplications,
    Graphics,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum BackgroundJobKind {
    PrepareGenome,
    BlastGenome,
    TrackImport,
    AgentAssist,
}

impl BackgroundJobKind {
    fn label(self) -> &'static str {
        match self {
            Self::PrepareGenome => "PrepareGenome",
            Self::BlastGenome => "BlastGenome",
            Self::TrackImport => "TrackImport",
            Self::AgentAssist => "AgentAssist",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum BackgroundJobEventPhase {
    Started,
    CancelRequested,
    Completed,
    Failed,
    Retried,
    IgnoredStale,
}

impl BackgroundJobEventPhase {
    fn label(self) -> &'static str {
        match self {
            Self::Started => "started",
            Self::CancelRequested => "cancel-requested",
            Self::Completed => "completed",
            Self::Failed => "failed",
            Self::Retried => "retried",
            Self::IgnoredStale => "ignored-stale",
        }
    }
}

#[derive(Clone, Debug)]
struct BackgroundJobEvent {
    job_id: Option<u64>,
    kind: BackgroundJobKind,
    phase: BackgroundJobEventPhase,
    emitted_at_unix_ms: u128,
    summary: String,
}

impl BackgroundJobEvent {
    fn to_line(&self) -> String {
        let job_token = self
            .job_id
            .map(|id| format!("#{id}"))
            .unwrap_or_else(|| "#-".to_string());
        format!(
            "[{} {} {} @{}] {}",
            self.kind.label(),
            job_token,
            self.phase.label(),
            self.emitted_at_unix_ms,
            self.summary
        )
    }
}

struct GenomePrepareTask {
    job_id: u64,
    started: Instant,
    cancel_requested: Arc<AtomicBool>,
    timeout_seconds: Option<u64>,
    receiver: mpsc::Receiver<GenomePrepareTaskMessage>,
}

enum GenomePrepareTaskMessage {
    Progress {
        job_id: u64,
        progress: PrepareGenomeProgress,
    },
    Done {
        job_id: u64,
        result: Result<OpResult, EngineError>,
    },
}

struct GenomeTrackImportTask {
    job_id: u64,
    started: Instant,
    cancel_requested: Arc<AtomicBool>,
    receiver: mpsc::Receiver<GenomeTrackImportTaskMessage>,
}

enum GenomeTrackImportTaskMessage {
    Progress {
        job_id: u64,
        progress: GenomeTrackImportProgress,
    },
    Done {
        job_id: u64,
        result: Result<OpResult, EngineError>,
    },
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum GenomeBlastSourceMode {
    Manual,
    ProjectSequence,
    ProjectPool,
}

struct GenomeBlastTask {
    job_id: u64,
    started: Instant,
    receiver: mpsc::Receiver<GenomeBlastTaskMessage>,
}

enum GenomeBlastTaskMessage {
    Progress {
        job_id: u64,
        done_queries: usize,
        total_queries: usize,
        current_query_label: String,
    },
    Done {
        job_id: u64,
        result: Result<GenomeBlastBatchResult, String>,
    },
}

#[derive(Clone)]
struct GenomeBlastQueryResult {
    query_label: String,
    query_length: usize,
    report: GenomeBlastReport,
}

#[derive(Clone)]
struct GenomeBlastBatchResult {
    reports: Vec<GenomeBlastQueryResult>,
    failed_queries: Vec<String>,
}

struct AgentAskTask {
    job_id: u64,
    started: Instant,
    receiver: mpsc::Receiver<AgentAskTaskMessage>,
}

enum AgentAskTaskMessage {
    Done {
        job_id: u64,
        result: Result<AgentInvocationOutcome, String>,
    },
}

struct AgentModelDiscoveryTask {
    started: Instant,
    source_key: String,
    receiver: mpsc::Receiver<AgentModelDiscoveryTaskMessage>,
}

enum AgentModelDiscoveryTaskMessage {
    Done {
        source_key: String,
        result: Result<Vec<String>, String>,
    },
}

#[derive(Clone)]
struct AgentCommandExecutionRecord {
    index_1based: usize,
    command: String,
    trigger: String,
    ok: bool,
    state_changed: bool,
    summary: String,
    executed_at_unix_ms: u128,
}

#[derive(Clone)]
struct BlastPoolOption {
    container_id: String,
    label: String,
    members: Vec<String>,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum LineageNodeKind {
    Sequence,
    Arrangement,
}

#[derive(Clone)]
struct LineageRow {
    kind: LineageNodeKind,
    node_id: String,
    seq_id: String,
    display_name: String,
    origin: String,
    created_by_op: String,
    created_at: u128,
    parents: Vec<String>,
    length: usize,
    circular: bool,
    pool_size: usize,
    pool_members: Vec<String>,
    arrangement_id: Option<String>,
    arrangement_mode: Option<String>,
    lane_container_ids: Vec<String>,
    ladders: Vec<String>,
}

#[derive(Clone)]
struct ContainerRow {
    container_id: String,
    kind: String,
    member_count: usize,
    representative: String,
    members: Vec<String>,
}

#[derive(Clone)]
struct ArrangementRow {
    arrangement_id: String,
    mode: String,
    name: String,
    created_by_op: String,
    created_at: u128,
    lane_count: usize,
    lane_container_ids: Vec<String>,
    ladders: Vec<String>,
}

#[derive(Clone)]
struct OpenWindowEntry {
    native_menu_key: u64,
    viewport_id: ViewportId,
    title: String,
    detail: String,
}

#[derive(Clone, Default)]
struct HelpMarkdownImage {
    alt: String,
    title: String,
    path: String,
}

#[derive(Clone, Copy)]
enum CommandPaletteAction {
    NewProject,
    OpenProject,
    SaveProject,
    OpenSequence,
    OpenConfiguration,
    OpenPrepareGenome,
    OpenRetrieveGenome,
    OpenBlastGenome,
    OpenGenomeTracks,
    OpenAgentAssistant,
    OpenPreparedInspector,
    OpenGuiManual,
    OpenCliManual,
    OpenShellManual,
    ExportLineageSvg,
    ToggleJobsPanel,
    ToggleHistoryPanel,
    Undo,
    Redo,
    FocusViewport(ViewportId),
}

#[derive(Clone)]
struct CommandPaletteEntry {
    title: String,
    detail: String,
    keywords: String,
    action: CommandPaletteAction,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
struct PersistedLineageGraphWorkspace {
    zoom: f32,
    graph_area_height: f32,
    container_area_height: f32,
    scroll_offset: [f32; 2],
    compact_labels: bool,
    node_offsets: HashMap<String, [f32; 2]>,
}

impl Default for PersistedLineageGraphWorkspace {
    fn default() -> Self {
        Self {
            zoom: 1.0,
            graph_area_height: 420.0,
            container_area_height: 220.0,
            scroll_offset: [0.0, 0.0],
            compact_labels: true,
            node_offsets: HashMap::new(),
        }
    }
}

impl Default for GENtleApp {
    fn default() -> Self {
        Self {
            engine: Arc::new(RwLock::new(GentleEngine::new())),
            new_windows: vec![],
            windows: HashMap::new(),
            windows_to_close: Arc::new(RwLock::new(vec![])),
            pending_focus_viewports: vec![],
            viewport_id_counter: 0,
            update_has_run_before: false,
            show_about_dialog: false,
            show_help_dialog: false,
            help_doc: HelpDoc::Gui,
            help_markdown_cache: CommonMarkCache::default(),
            help_gui_markdown: GUI_MANUAL_MD.to_string(),
            help_cli_markdown: CLI_MANUAL_MD.to_string(),
            help_shell_markdown: Self::generate_shell_help_markdown(),
            help_shell_interface: ShellHelpInterface::GuiShell,
            help_search_query: String::new(),
            help_search_matches: vec![],
            help_search_selected: 0,
            help_focus_search_box: false,
            help_image_preview_path: None,
            help_image_preview_caption: String::new(),
            show_configuration_dialog: false,
            configuration_tab: ConfigurationTab::ExternalApplications,
            configuration_rnapkin_executable: env::var("GENTLE_RNAPKIN_BIN").unwrap_or_default(),
            configuration_makeblastdb_executable: env::var(MAKEBLASTDB_ENV_BIN).unwrap_or_default(),
            configuration_blastn_executable: env::var(BLASTN_ENV_BIN).unwrap_or_default(),
            configuration_bigwig_to_bedgraph_executable: env::var(BIGWIG_TO_BEDGRAPH_ENV_BIN)
                .unwrap_or_default(),
            configuration_rnapkin_validation_ok: None,
            configuration_rnapkin_validation_message: String::new(),
            configuration_blast_validation_ok: None,
            configuration_blast_validation_message: String::new(),
            configuration_graphics: DisplaySettings::default(),
            configuration_graphics_dirty: false,
            window_backdrops: WindowBackdropSettings::default(),
            configuration_window_backdrops: WindowBackdropSettings::default(),
            configuration_window_backdrops_dirty: false,
            configuration_status: String::new(),
            current_project_path: None,
            recent_project_paths: vec![],
            lineage_graph_view: false,
            lineage_graph_zoom: 1.0,
            lineage_graph_area_height: 420.0,
            lineage_container_area_height: 220.0,
            lineage_graph_scroll_offset: Vec2::ZERO,
            lineage_graph_pan_origin: None,
            lineage_graph_compact_labels: true,
            lineage_graph_selected_node_id: None,
            lineage_graph_node_offsets: HashMap::new(),
            lineage_graph_drag_origin: None,
            lineage_graph_offsets_synced_stamp: 0,
            lineage_cache_stamp: 0,
            lineage_cache_valid: false,
            lineage_rows: vec![],
            lineage_edges: vec![],
            lineage_op_label_by_id: HashMap::new(),
            lineage_containers: vec![],
            lineage_arrangements: vec![],
            clean_state_fingerprint: 0,
            dirty_cache_stamp: 0,
            dirty_cache_valid: false,
            dirty_cache_value: false,
            dirty_cache_last_deep_check: Instant::now(),
            last_applied_window_title: String::new(),
            last_native_window_entries: vec![],
            last_native_active_window_key: Some(viewport_native_menu_key(ViewportId::ROOT)),
            native_window_key_to_viewport: HashMap::new(),
            active_window_menu_key: Some(viewport_native_menu_key(ViewportId::ROOT)),
            pending_project_action: None,
            show_reference_genome_prepare_dialog: false,
            show_reference_genome_retrieve_dialog: false,
            show_reference_genome_blast_dialog: false,
            show_reference_genome_inspector_dialog: false,
            genome_catalog_path: DEFAULT_GENOME_CATALOG_PATH.to_string(),
            genome_cache_dir: DEFAULT_GENOME_CACHE_DIR.to_string(),
            genome_id: "Human GRCh38 Ensembl 113".to_string(),
            genome_catalog_genomes: vec![],
            genome_catalog_error: String::new(),
            genome_prepare_task: None,
            genome_prepare_progress: None,
            genome_prepare_status: String::new(),
            genome_prepare_timeout_secs: default_prepare_timeout_secs_string(),
            genome_gene_filter: String::new(),
            genome_gene_filter_limit: 5000,
            genome_gene_filter_page: 0,
            genome_biotype_filter: HashMap::new(),
            genome_genes: vec![],
            genome_genes_loaded_key: None,
            genome_genes_error: String::new(),
            genome_selected_gene: None,
            genome_chromosome: "1".to_string(),
            genome_start_1based: "1".to_string(),
            genome_end_1based: "1000".to_string(),
            genome_output_id: String::new(),
            genome_retrieve_status: String::new(),
            genome_blast_source_mode: GenomeBlastSourceMode::Manual,
            genome_blast_query_manual: String::new(),
            genome_blast_query_seq_id: String::new(),
            genome_blast_query_pool_id: String::new(),
            genome_blast_max_hits: 25,
            genome_blast_task_name: "blastn-short".to_string(),
            genome_blast_task: None,
            genome_blast_progress_fraction: None,
            genome_blast_progress_label: String::new(),
            genome_blast_results: vec![],
            genome_blast_selected_result: 0,
            genome_blast_status: String::new(),
            genome_blast_import_track_name: "blast_hits".to_string(),
            genome_blast_import_clear_existing: false,
            show_genome_bed_track_dialog: false,
            show_agent_assistant_dialog: false,
            genome_track_seq_id: String::new(),
            genome_track_source_selection: GenomeTrackSourceSelection::Auto,
            genome_track_path: String::new(),
            genome_track_name: String::new(),
            genome_track_min_score: String::new(),
            genome_track_max_score: String::new(),
            genome_track_clear_existing: false,
            genome_track_status: String::new(),
            genome_track_import_task: None,
            genome_track_import_progress: None,
            genome_bed_track_subscriptions: vec![],
            genome_track_autosync_status: String::new(),
            tracked_autosync_last_op_count: None,
            show_command_palette_dialog: false,
            command_palette_query: String::new(),
            command_palette_selected: 0,
            command_palette_focus_query: false,
            show_jobs_panel: false,
            show_history_panel: false,
            hover_status_name: String::new(),
            app_status: String::new(),
            next_background_job_id: 1,
            job_event_log: vec![],
            genome_track_preflight_track_subscription: true,
            agent_catalog_path: DEFAULT_AGENT_SYSTEM_CATALOG_PATH.to_string(),
            agent_catalog_loaded_path: String::new(),
            agent_catalog_error: String::new(),
            agent_systems: vec![],
            agent_system_id: "builtin_echo".to_string(),
            agent_openai_api_key: String::new(),
            agent_base_url_override: String::new(),
            agent_model_override: String::new(),
            agent_timeout_secs: default_agent_timeout_secs_string(),
            agent_connect_timeout_secs: default_agent_connect_timeout_secs_string(),
            agent_read_timeout_secs: default_agent_read_timeout_secs_string(),
            agent_max_retries: default_agent_max_retries_string(),
            agent_max_response_bytes: default_agent_max_response_bytes_string(),
            agent_discovered_models: vec![],
            agent_discovered_model_pick: String::new(),
            agent_model_discovery_status: String::new(),
            agent_model_discovery_source_key: String::new(),
            agent_model_discovery_task: None,
            agent_prompt_template_id: AGENT_PROMPT_TEMPLATE_DEFAULT_ID.to_string(),
            agent_prompt: String::new(),
            agent_include_state_summary: true,
            agent_allow_auto_exec: false,
            agent_status: String::new(),
            agent_task: None,
            agent_last_invocation: None,
            agent_execution_log: vec![],
        }
    }
}

impl GENtleApp {
    fn help_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("GENtle Help Viewport")
    }

    fn prepare_genome_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("GENtle Prepare Genome Viewport")
    }

    fn blast_genome_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("GENtle BLAST Genome Viewport")
    }

    fn configuration_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("GENtle Configuration Viewport")
    }

    fn bed_track_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("GENtle BED Tracks Viewport")
    }

    fn agent_assistant_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("GENtle Agent Assistant Viewport")
    }

    fn command_palette_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("GENtle Command Palette Viewport")
    }

    fn history_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("GENtle Operation History Viewport")
    }

    fn native_menu_key_for_viewport(viewport_id: ViewportId) -> u64 {
        viewport_native_menu_key(viewport_id)
    }

    fn set_active_window_viewport(&mut self, viewport_id: ViewportId) {
        self.active_window_menu_key = Some(Self::native_menu_key_for_viewport(viewport_id));
        report_active_viewport_from_ui(viewport_id);
    }

    fn note_viewport_focus_if_active(&mut self, ctx: &egui::Context, viewport_id: ViewportId) {
        if ctx.input(|i| i.viewport().focused.unwrap_or(false)) {
            self.set_active_window_viewport(viewport_id);
        }
    }

    fn configuration_store_path() -> PathBuf {
        let home = env::var_os("HOME")
            .map(PathBuf::from)
            .unwrap_or_else(|| PathBuf::from("."));
        home.join(APP_CONFIGURATION_FILE_NAME)
    }

    fn read_persisted_configuration_from_disk() -> Option<PersistedConfiguration> {
        let path = Self::configuration_store_path();
        let text = fs::read_to_string(&path).ok()?;
        serde_json::from_str::<PersistedConfiguration>(&text).ok()
    }

    fn write_persisted_configuration_to_disk(&self) -> std::result::Result<(), String> {
        let path = Self::configuration_store_path();
        if let Some(parent) = path.parent() {
            if !parent.as_os_str().is_empty() {
                fs::create_dir_all(parent).map_err(|e| {
                    format!(
                        "Could not create configuration directory '{}': {e}",
                        parent.display()
                    )
                })?;
            }
        }
        let payload = PersistedConfiguration {
            rnapkin_executable: self.configuration_rnapkin_executable.trim().to_string(),
            makeblastdb_executable: self.configuration_makeblastdb_executable.trim().to_string(),
            blastn_executable: self.configuration_blastn_executable.trim().to_string(),
            bigwig_to_bedgraph_executable: self
                .configuration_bigwig_to_bedgraph_executable
                .trim()
                .to_string(),
            graphics_defaults: self.configuration_graphics.clone(),
            window_backdrops: self.window_backdrops.clone(),
            recent_projects: self.recent_project_paths.clone(),
        };
        let json = serde_json::to_string_pretty(&payload)
            .map_err(|e| format!("Could not serialize GUI configuration: {e}"))?;
        fs::write(&path, json).map_err(|e| {
            format!(
                "Could not write GUI configuration '{}': {e}",
                path.display()
            )
        })
    }

    fn normalize_project_path(path: &str) -> String {
        let trimmed = path.trim();
        if trimmed.is_empty() {
            return String::new();
        }
        let input = PathBuf::from(trimmed);
        let absolute = if input.is_absolute() {
            input
        } else {
            env::current_dir()
                .unwrap_or_else(|_| PathBuf::from("."))
                .join(input)
        };
        absolute
            .canonicalize()
            .unwrap_or(absolute)
            .to_string_lossy()
            .to_string()
    }

    fn normalize_recent_project_paths(paths: Vec<String>) -> Vec<String> {
        let mut normalized = Vec::new();
        let mut seen = HashSet::new();
        for raw in paths {
            let normalized_path = Self::normalize_project_path(&raw);
            if normalized_path.is_empty() {
                continue;
            }
            if seen.insert(normalized_path.clone()) {
                normalized.push(normalized_path);
                if normalized.len() >= MAX_RECENT_PROJECTS {
                    break;
                }
            }
        }
        normalized
    }

    fn record_recent_project_path(&mut self, path: &str) {
        let normalized = Self::normalize_project_path(path);
        if normalized.is_empty() {
            return;
        }

        let old_paths = self.recent_project_paths.clone();
        self.recent_project_paths
            .retain(|existing| existing != &normalized);
        self.recent_project_paths.insert(0, normalized);
        if self.recent_project_paths.len() > MAX_RECENT_PROJECTS {
            self.recent_project_paths.truncate(MAX_RECENT_PROJECTS);
        }

        if self.recent_project_paths != old_paths {
            if let Err(err) = self.write_persisted_configuration_to_disk() {
                eprintln!("{err}");
            }
        }
    }

    fn remove_recent_project_path(&mut self, path: &str) {
        let normalized = Self::normalize_project_path(path);
        if normalized.is_empty() {
            return;
        }
        let len_before = self.recent_project_paths.len();
        self.recent_project_paths
            .retain(|existing| existing != &normalized);
        if self.recent_project_paths.len() != len_before {
            if let Err(err) = self.write_persisted_configuration_to_disk() {
                eprintln!("{err}");
            }
        }
    }

    fn clear_recent_project_paths(&mut self) {
        if self.recent_project_paths.is_empty() {
            return;
        }
        self.recent_project_paths.clear();
        if let Err(err) = self.write_persisted_configuration_to_disk() {
            eprintln!("{err}");
        }
    }

    fn recent_project_menu_label(path: &str) -> String {
        let parsed = Path::new(path);
        let name = parsed
            .file_name()
            .map(|s| s.to_string_lossy().to_string())
            .unwrap_or_else(|| path.to_string());
        let parent = parsed
            .parent()
            .map(|p| p.display().to_string())
            .unwrap_or_default();
        if parent.is_empty() {
            name
        } else {
            format!("{name} ({parent})")
        }
    }

    fn apply_graphics_settings_to_display(source: &DisplaySettings, target: &mut DisplaySettings) {
        target.show_sequence_panel = source.show_sequence_panel;
        target.show_map_panel = source.show_map_panel;
        target.show_features = source.show_features;
        target.show_cds_features = source.show_cds_features;
        target.show_gene_features = source.show_gene_features;
        target.show_mrna_features = source.show_mrna_features;
        target.show_tfbs = source.show_tfbs;
        target.regulatory_tracks_near_baseline = source.regulatory_tracks_near_baseline;
        target.regulatory_feature_max_view_span_bp = source.regulatory_feature_max_view_span_bp;
        target.tfbs_display_use_llr_bits = source.tfbs_display_use_llr_bits;
        target.tfbs_display_min_llr_bits = source.tfbs_display_min_llr_bits;
        target.tfbs_display_use_llr_quantile = source.tfbs_display_use_llr_quantile;
        target.tfbs_display_min_llr_quantile = source.tfbs_display_min_llr_quantile;
        target.tfbs_display_use_true_log_odds_bits = source.tfbs_display_use_true_log_odds_bits;
        target.tfbs_display_min_true_log_odds_bits = source.tfbs_display_min_true_log_odds_bits;
        target.tfbs_display_use_true_log_odds_quantile =
            source.tfbs_display_use_true_log_odds_quantile;
        target.tfbs_display_min_true_log_odds_quantile =
            source.tfbs_display_min_true_log_odds_quantile;
        target.vcf_display_show_snp = source.vcf_display_show_snp;
        target.vcf_display_show_ins = source.vcf_display_show_ins;
        target.vcf_display_show_del = source.vcf_display_show_del;
        target.vcf_display_show_sv = source.vcf_display_show_sv;
        target.vcf_display_show_other = source.vcf_display_show_other;
        target.vcf_display_pass_only = source.vcf_display_pass_only;
        target.vcf_display_use_min_qual = source.vcf_display_use_min_qual;
        target.vcf_display_min_qual = source.vcf_display_min_qual;
        target.vcf_display_use_max_qual = source.vcf_display_use_max_qual;
        target.vcf_display_max_qual = source.vcf_display_max_qual;
        target.vcf_display_required_info_keys = source.vcf_display_required_info_keys.clone();
        target.show_restriction_enzymes = source.show_restriction_enzymes;
        target.show_gc_contents = source.show_gc_contents;
        target.show_open_reading_frames = source.show_open_reading_frames;
        target.show_methylation_sites = source.show_methylation_sites;
        target.feature_details_font_size = source.feature_details_font_size;
        target.linear_view_start_bp = source.linear_view_start_bp;
        target.linear_view_span_bp = source.linear_view_span_bp;
    }

    fn apply_configuration_graphics_to_engine_state(&mut self) {
        let mut guard = self.engine.write().expect("Engine lock poisoned");
        let display = &mut guard.state_mut().display;
        Self::apply_graphics_settings_to_display(&self.configuration_graphics, display);
        self.configuration_graphics_dirty = false;
    }

    fn sync_runtime_tool_overrides_from_configuration(&self) {
        tool_overrides::set_tool_override(
            "GENTLE_RNAPKIN_BIN",
            &self.configuration_rnapkin_executable,
        );
        tool_overrides::set_tool_override(
            MAKEBLASTDB_ENV_BIN,
            &self.configuration_makeblastdb_executable,
        );
        tool_overrides::set_tool_override(BLASTN_ENV_BIN, &self.configuration_blastn_executable);
        tool_overrides::set_tool_override(
            BIGWIG_TO_BEDGRAPH_ENV_BIN,
            &self.configuration_bigwig_to_bedgraph_executable,
        );
    }

    fn load_persisted_configuration(&mut self, apply_graphics_to_current_project: bool) {
        let Some(saved) = Self::read_persisted_configuration_from_disk() else {
            return;
        };
        self.configuration_rnapkin_executable = saved.rnapkin_executable.trim().to_string();
        self.configuration_makeblastdb_executable = saved.makeblastdb_executable.trim().to_string();
        self.configuration_blastn_executable = saved.blastn_executable.trim().to_string();
        self.configuration_bigwig_to_bedgraph_executable =
            saved.bigwig_to_bedgraph_executable.trim().to_string();
        self.sync_runtime_tool_overrides_from_configuration();
        self.configuration_graphics = saved.graphics_defaults;
        let mut runtime_backdrops = saved.window_backdrops;
        runtime_backdrops.apply_runtime_defaults_if_legacy();
        self.window_backdrops = runtime_backdrops.clone();
        self.configuration_window_backdrops = runtime_backdrops;
        self.configuration_window_backdrops_dirty = false;
        window_backdrop::set_window_backdrop_settings(self.window_backdrops.clone());
        self.recent_project_paths = Self::normalize_recent_project_paths(saved.recent_projects);
        self.configuration_graphics_dirty = false;
        self.configuration_rnapkin_validation_ok = None;
        self.configuration_rnapkin_validation_message.clear();
        self.configuration_blast_validation_ok = None;
        self.configuration_blast_validation_message.clear();

        if apply_graphics_to_current_project {
            self.apply_configuration_graphics_to_engine_state();
        }
    }

    fn resolve_runtime_doc_path(path: &str) -> Option<PathBuf> {
        let direct = PathBuf::from(path);
        if direct.exists() {
            return Some(direct);
        }
        if let Ok(exe_path) = env::current_exe() {
            if let Some(exe_dir) = exe_path.parent() {
                let bundled = exe_dir.join("../Resources").join(path);
                if bundled.exists() {
                    return Some(bundled);
                }
            }
        }
        None
    }

    fn rewrite_markdown_relative_image_links(markdown: &str, base_dir: &Path) -> String {
        let mut replacements = Vec::new();
        for (event, range) in Parser::new(markdown).into_offset_iter() {
            let Event::Start(Tag::Image {
                link_type,
                dest_url,
                ..
            }) = event
            else {
                continue;
            };
            if link_type != LinkType::Inline {
                continue;
            }

            let Some(abs_path) = Self::resolve_relative_image_path(dest_url.as_ref(), base_dir)
            else {
                continue;
            };

            let span = &markdown[range.clone()];
            let Some(rewritten_span) =
                Self::rewrite_inline_image_destination(span, abs_path.as_path())
            else {
                continue;
            };
            replacements.push((range, rewritten_span));
        }

        if replacements.is_empty() {
            return markdown.to_string();
        }

        let mut rewritten = markdown.to_string();
        replacements.sort_by(|a, b| b.0.start.cmp(&a.0.start));
        for (range, replacement) in replacements {
            rewritten.replace_range(range, &replacement);
        }
        rewritten
    }

    fn rewrite_inline_image_destination(span: &str, absolute_dest: &Path) -> Option<String> {
        let (dest_start, dest_end) = Self::find_inline_image_destination(span)?;
        let mut rewritten = String::with_capacity(span.len() + 64);
        rewritten.push_str(&span[..dest_start]);
        rewritten.push_str(&absolute_dest.to_string_lossy());
        rewritten.push_str(&span[dest_end..]);
        Some(rewritten)
    }

    fn find_inline_image_destination(span: &str) -> Option<(usize, usize)> {
        let bytes = span.as_bytes();
        if bytes.len() < 4 || bytes[0] != b'!' || bytes[1] != b'[' {
            return None;
        }

        let mut i = 2usize;
        let mut bracket_depth = 1usize;
        while i < bytes.len() {
            match bytes[i] {
                b'\\' => {
                    i = (i + 2).min(bytes.len());
                }
                b'[' => {
                    bracket_depth += 1;
                    i += 1;
                }
                b']' => {
                    bracket_depth = bracket_depth.saturating_sub(1);
                    i += 1;
                    if bracket_depth == 0 {
                        break;
                    }
                }
                _ => i += 1,
            }
        }
        if bracket_depth != 0 {
            return None;
        }

        while i < bytes.len() && bytes[i].is_ascii_whitespace() {
            i += 1;
        }
        if i >= bytes.len() || bytes[i] != b'(' {
            return None;
        }
        i += 1;

        while i < bytes.len() && bytes[i].is_ascii_whitespace() {
            i += 1;
        }
        if i >= bytes.len() {
            return None;
        }

        if bytes[i] == b'<' {
            let dest_start = i + 1;
            i += 1;
            while i < bytes.len() {
                match bytes[i] {
                    b'\\' => {
                        i = (i + 2).min(bytes.len());
                    }
                    b'>' => return Some((dest_start, i)),
                    _ => i += 1,
                }
            }
            return None;
        }

        let dest_start = i;
        let mut paren_depth = 0usize;
        while i < bytes.len() {
            match bytes[i] {
                b'\\' => {
                    i = (i + 2).min(bytes.len());
                }
                b'(' => {
                    paren_depth += 1;
                    i += 1;
                }
                b')' => {
                    if paren_depth == 0 {
                        break;
                    }
                    paren_depth -= 1;
                    i += 1;
                }
                b if b.is_ascii_whitespace() && paren_depth == 0 => break,
                _ => i += 1,
            }
        }
        if dest_start == i {
            return None;
        }
        Some((dest_start, i))
    }

    fn resolve_relative_image_path(path: &str, base_dir: &Path) -> Option<PathBuf> {
        if path.is_empty()
            || path.starts_with('/')
            || path.starts_with('\\')
            || path.starts_with('#')
            || path.contains("://")
            || path.contains(':')
        {
            return None;
        }
        let joined = base_dir.join(path);
        Some(joined.canonicalize().unwrap_or(joined))
    }

    fn load_help_doc(path: &str, fallback: &'static str) -> String {
        if let Some(runtime_path) = Self::resolve_runtime_doc_path(path) {
            if let Ok(text) = fs::read_to_string(&runtime_path) {
                if let Some(base_dir) = runtime_path.parent() {
                    return Self::rewrite_markdown_relative_image_links(&text, base_dir);
                }
                return text;
            }
        }
        fallback.to_string()
    }

    fn generate_shell_help_markdown_for(interface: ShellHelpInterface) -> String {
        let interface_filter = interface.glossary_filter();
        if let Some(runtime_path) = Self::resolve_runtime_doc_path("docs/glossary.json") {
            if let Ok(raw) = fs::read_to_string(runtime_path) {
                if let Ok(markdown) =
                    render_shell_help_markdown_from_glossary_json(&raw, interface_filter)
                {
                    return markdown;
                }
            }
        }
        match render_shell_help_markdown(interface_filter) {
            Ok(markdown) => markdown,
            Err(err) => format!(
                "# GENtle Shell Command Reference\n\n\
Could not render command glossary from `docs/glossary.json`.\n\n\
Error: `{err}`"
            ),
        }
    }

    fn generate_shell_help_markdown() -> String {
        Self::generate_shell_help_markdown_for(ShellHelpInterface::GuiShell)
    }

    fn refresh_help_docs(&mut self) {
        self.help_gui_markdown = Self::load_help_doc("docs/gui.md", GUI_MANUAL_MD);
        self.help_cli_markdown = Self::load_help_doc("docs/cli.md", CLI_MANUAL_MD);
        self.help_shell_markdown =
            Self::generate_shell_help_markdown_for(self.help_shell_interface);
    }

    fn consume_native_help_request(&mut self) {
        if NATIVE_HELP_OPEN_REQUESTED.swap(false, Ordering::SeqCst) {
            self.open_help_doc(HelpDoc::Gui);
        }
    }

    fn consume_native_settings_request(&mut self) {
        if NATIVE_SETTINGS_OPEN_REQUESTED.swap(false, Ordering::SeqCst) {
            self.open_configuration_dialog();
        }
    }

    fn consume_native_windows_request(&mut self) {
        if NATIVE_WINDOWS_OPEN_REQUESTED.swap(false, Ordering::SeqCst) {
            let requested_key = NATIVE_WINDOWS_FOCUS_KEY_REQUESTED
                .swap(NATIVE_WINDOW_FOCUS_KEY_NONE, Ordering::SeqCst);
            if requested_key == NATIVE_WINDOW_FOCUS_KEY_NONE {
                self.queue_focus_viewport(ViewportId::ROOT);
                return;
            }
            let target_viewport = self
                .native_window_key_to_viewport
                .get(&requested_key)
                .copied()
                .unwrap_or(ViewportId::ROOT);
            self.queue_focus_viewport(target_viewport);
        }
    }

    fn consume_active_viewport_report(&mut self) {
        let reported_key =
            ACTIVE_VIEWPORT_KEY_REPORTED.swap(NATIVE_WINDOW_FOCUS_KEY_NONE, Ordering::SeqCst);
        if reported_key != NATIVE_WINDOW_FOCUS_KEY_NONE {
            self.active_window_menu_key = Some(reported_key);
        }
    }

    fn resolved_rnapkin_executable(&self) -> String {
        let configured = self.configuration_rnapkin_executable.trim();
        if configured.is_empty() {
            "rnapkin".to_string()
        } else {
            configured.to_string()
        }
    }

    fn resolved_makeblastdb_executable(&self) -> String {
        let configured = self.configuration_makeblastdb_executable.trim();
        if configured.is_empty() {
            DEFAULT_MAKEBLASTDB_BIN.to_string()
        } else {
            configured.to_string()
        }
    }

    fn resolved_blastn_executable(&self) -> String {
        let configured = self.configuration_blastn_executable.trim();
        if configured.is_empty() {
            DEFAULT_BLASTN_BIN.to_string()
        } else {
            configured.to_string()
        }
    }

    fn resolved_bigwig_to_bedgraph_executable(&self) -> String {
        let configured = self.configuration_bigwig_to_bedgraph_executable.trim();
        if configured.is_empty() {
            DEFAULT_BIGWIG_TO_BEDGRAPH_BIN.to_string()
        } else {
            configured.to_string()
        }
    }

    fn clear_rnapkin_validation(&mut self) {
        self.configuration_rnapkin_validation_ok = None;
        self.configuration_rnapkin_validation_message.clear();
    }

    fn clear_blast_validation(&mut self) {
        self.configuration_blast_validation_ok = None;
        self.configuration_blast_validation_message.clear();
    }

    fn first_non_empty_output_line(stdout: &[u8], stderr: &[u8]) -> String {
        String::from_utf8_lossy(stdout)
            .lines()
            .chain(String::from_utf8_lossy(stderr).lines())
            .find(|line| !line.trim().is_empty())
            .map(|s| s.trim().to_string())
            .unwrap_or_else(|| "No version output".to_string())
    }

    fn validate_rnapkin_executable(&mut self) {
        let executable = self.resolved_rnapkin_executable();
        match Command::new(&executable).arg("--version").output() {
            Ok(output) => {
                let first_non_empty =
                    Self::first_non_empty_output_line(&output.stdout, &output.stderr);
                if output.status.success() {
                    self.configuration_rnapkin_validation_ok = Some(true);
                    self.configuration_rnapkin_validation_message = format!(
                        "rnapkin validation succeeded using '{}': {}",
                        executable, first_non_empty
                    );
                } else {
                    self.configuration_rnapkin_validation_ok = Some(false);
                    self.configuration_rnapkin_validation_message = format!(
                        "rnapkin validation failed using '{}' (status {:?}): {}",
                        executable,
                        output.status.code(),
                        first_non_empty
                    );
                }
            }
            Err(e) if e.kind() == ErrorKind::NotFound => {
                self.configuration_rnapkin_validation_ok = Some(false);
                self.configuration_rnapkin_validation_message =
                    format!("rnapkin executable '{}' not found", executable);
            }
            Err(e) => {
                self.configuration_rnapkin_validation_ok = Some(false);
                self.configuration_rnapkin_validation_message =
                    format!("Could not execute '{}': {}", executable, e);
            }
        }
    }

    fn validate_blast_executables(&mut self) {
        let makeblastdb = self.resolved_makeblastdb_executable();
        let blastn = self.resolved_blastn_executable();

        let mut ok = true;
        let makeblastdb_message = match Command::new(&makeblastdb).arg("-version").output() {
            Ok(output) => {
                let first_non_empty =
                    Self::first_non_empty_output_line(&output.stdout, &output.stderr);
                if output.status.success() {
                    format!("makeblastdb OK ('{}'): {}", makeblastdb, first_non_empty)
                } else {
                    ok = false;
                    format!(
                        "makeblastdb failed ('{}', status {:?}): {}",
                        makeblastdb,
                        output.status.code(),
                        first_non_empty
                    )
                }
            }
            Err(e) if e.kind() == ErrorKind::NotFound => {
                ok = false;
                format!("makeblastdb executable '{}' not found", makeblastdb)
            }
            Err(e) => {
                ok = false;
                format!("Could not execute makeblastdb '{}': {}", makeblastdb, e)
            }
        };

        let blastn_message = match Command::new(&blastn).arg("-version").output() {
            Ok(output) => {
                let first_non_empty =
                    Self::first_non_empty_output_line(&output.stdout, &output.stderr);
                if output.status.success() {
                    format!("blastn OK ('{}'): {}", blastn, first_non_empty)
                } else {
                    ok = false;
                    format!(
                        "blastn failed ('{}', status {:?}): {}",
                        blastn,
                        output.status.code(),
                        first_non_empty
                    )
                }
            }
            Err(e) if e.kind() == ErrorKind::NotFound => {
                ok = false;
                format!("blastn executable '{}' not found", blastn)
            }
            Err(e) => {
                ok = false;
                format!("Could not execute blastn '{}': {}", blastn, e)
            }
        };

        self.configuration_blast_validation_ok = Some(ok);
        self.configuration_blast_validation_message =
            format!("{makeblastdb_message} | {blastn_message}");
    }

    fn sync_configuration_from_runtime(&mut self) {
        self.configuration_rnapkin_executable =
            tool_overrides::configured_or_env("GENTLE_RNAPKIN_BIN");
        self.configuration_makeblastdb_executable =
            tool_overrides::configured_or_env(MAKEBLASTDB_ENV_BIN);
        self.configuration_blastn_executable = tool_overrides::configured_or_env(BLASTN_ENV_BIN);
        self.configuration_bigwig_to_bedgraph_executable =
            tool_overrides::configured_or_env(BIGWIG_TO_BEDGRAPH_ENV_BIN);
        if let Ok(engine) = self.engine.read() {
            self.configuration_graphics = engine.state().display.clone();
            self.configuration_graphics_dirty = false;
        }
        self.configuration_window_backdrops = self.window_backdrops.clone();
        self.configuration_window_backdrops_dirty = false;
        self.clear_rnapkin_validation();
        self.clear_blast_validation();
    }

    fn open_configuration_dialog(&mut self) {
        self.sync_configuration_from_runtime();
        self.configuration_tab = ConfigurationTab::ExternalApplications;
        self.show_configuration_dialog = true;
        self.configuration_status.clear();
    }

    fn open_command_palette_dialog(&mut self) {
        if self.show_command_palette_dialog {
            self.queue_focus_viewport(Self::command_palette_viewport_id());
        }
        self.show_command_palette_dialog = true;
        self.command_palette_focus_query = true;
    }

    fn track_hover_status<S: Into<String>>(
        &mut self,
        response: egui::Response,
        stable_name: S,
    ) -> egui::Response {
        if response.hovered() {
            self.hover_status_name = stable_name.into();
        }
        response
    }

    fn alloc_background_job_id(&mut self) -> u64 {
        let next = self.next_background_job_id.max(1);
        self.next_background_job_id = next.saturating_add(1);
        next
    }

    fn now_unix_ms() -> u128 {
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map(|v| v.as_millis())
            .unwrap_or_default()
    }

    fn push_job_event<S: Into<String>>(
        &mut self,
        kind: BackgroundJobKind,
        phase: BackgroundJobEventPhase,
        job_id: Option<u64>,
        summary: S,
    ) {
        let summary = summary.into();
        let trimmed = summary.trim();
        if trimmed.is_empty() {
            return;
        }
        self.job_event_log.push(BackgroundJobEvent {
            job_id,
            kind,
            phase,
            emitted_at_unix_ms: Self::now_unix_ms(),
            summary: trimmed.to_string(),
        });
        if self.job_event_log.len() > 200 {
            let drain_len = self.job_event_log.len() - 200;
            self.job_event_log.drain(0..drain_len);
        }
    }

    fn request_prepare_task_cancel(&mut self, origin: &str) {
        let Some((job_id, already_requested)) = self.genome_prepare_task.as_ref().map(|task| {
            (
                task.job_id,
                task.cancel_requested.swap(true, Ordering::Relaxed),
            )
        }) else {
            self.genome_prepare_status = "No running genome preparation task to cancel".to_string();
            return;
        };

        if already_requested {
            self.genome_prepare_status =
                "Cancellation was already requested for genome preparation".to_string();
            return;
        }

        self.genome_prepare_status =
            "Cancellation requested for running genome preparation".to_string();
        self.push_job_event(
            BackgroundJobKind::PrepareGenome,
            BackgroundJobEventPhase::CancelRequested,
            Some(job_id),
            format!("Cancellation requested from {origin}"),
        );
    }

    fn request_track_import_task_cancel(&mut self, origin: &str) {
        let Some((job_id, already_requested)) =
            self.genome_track_import_task.as_ref().map(|task| {
                (
                    task.job_id,
                    task.cancel_requested.swap(true, Ordering::Relaxed),
                )
            })
        else {
            self.genome_track_status = "No running track-import job to cancel".to_string();
            return;
        };

        if already_requested {
            self.genome_track_status =
                "Cancellation was already requested for track import".to_string();
            return;
        }

        self.genome_track_status = "Cancellation requested for running track import".to_string();
        self.push_job_event(
            BackgroundJobKind::TrackImport,
            BackgroundJobEventPhase::CancelRequested,
            Some(job_id),
            format!("Cancellation requested from {origin}"),
        );
    }

    fn has_active_background_jobs(&self) -> bool {
        self.genome_prepare_task.is_some()
            || self.genome_blast_task.is_some()
            || self.genome_track_import_task.is_some()
            || self.agent_task.is_some()
    }

    fn refresh_sequence_windows_from_engine_state(&mut self) {
        for window in self.windows.values() {
            if let Ok(mut guard) = window.write() {
                guard.refresh_from_engine_settings();
            }
        }
    }

    fn handle_engine_state_after_history_transition(&mut self) {
        self.lineage_cache_valid = false;
        self.lineage_rows.clear();
        self.lineage_edges.clear();
        self.lineage_op_label_by_id.clear();
        self.lineage_containers.clear();
        self.lineage_arrangements.clear();
        self.load_bed_track_subscriptions_from_state();
        self.tracked_autosync_last_op_count = None;
        self.genome_track_autosync_status.clear();
        self.refresh_sequence_windows_from_engine_state();
    }

    fn undo_last_operation(&mut self) {
        if self.has_active_background_jobs() {
            self.app_status =
                "Undo is disabled while background jobs are active (prepare/blast/import)."
                    .to_string();
            return;
        }
        let outcome = self
            .engine
            .write()
            .map_err(|_| EngineError {
                code: ErrorCode::Internal,
                message: "Engine lock poisoned".to_string(),
            })
            .and_then(|mut engine| engine.undo_last_operation());
        match outcome {
            Ok(()) => {
                self.handle_engine_state_after_history_transition();
                self.app_status = "Undo applied".to_string();
            }
            Err(e) => {
                self.app_status = format!("Undo unavailable: {}", e.message);
            }
        }
    }

    fn redo_last_operation(&mut self) {
        if self.has_active_background_jobs() {
            self.app_status =
                "Redo is disabled while background jobs are active (prepare/blast/import)."
                    .to_string();
            return;
        }
        let outcome = self
            .engine
            .write()
            .map_err(|_| EngineError {
                code: ErrorCode::Internal,
                message: "Engine lock poisoned".to_string(),
            })
            .and_then(|mut engine| engine.redo_last_operation());
        match outcome {
            Ok(()) => {
                self.handle_engine_state_after_history_transition();
                self.app_status = "Redo applied".to_string();
            }
            Err(e) => {
                self.app_status = format!("Redo unavailable: {}", e.message);
            }
        }
    }

    fn collect_command_palette_entries(&self) -> Vec<CommandPaletteEntry> {
        let mut entries = vec![
            CommandPaletteEntry {
                title: "New Project".to_string(),
                detail: "Create a new empty project".to_string(),
                keywords: "file project new".to_string(),
                action: CommandPaletteAction::NewProject,
            },
            CommandPaletteEntry {
                title: "Open Project".to_string(),
                detail: "Open an existing project file".to_string(),
                keywords: "file project open".to_string(),
                action: CommandPaletteAction::OpenProject,
            },
            CommandPaletteEntry {
                title: "Save Project".to_string(),
                detail: "Save current project state".to_string(),
                keywords: "file project save".to_string(),
                action: CommandPaletteAction::SaveProject,
            },
            CommandPaletteEntry {
                title: "Open Sequence".to_string(),
                detail: "Import a sequence file into project".to_string(),
                keywords: "sequence import load".to_string(),
                action: CommandPaletteAction::OpenSequence,
            },
            CommandPaletteEntry {
                title: "Undo Last Operation".to_string(),
                detail: "Restore previous project state".to_string(),
                keywords: "undo history edit".to_string(),
                action: CommandPaletteAction::Undo,
            },
            CommandPaletteEntry {
                title: "Redo Last Operation".to_string(),
                detail: "Re-apply previously undone state".to_string(),
                keywords: "redo history edit".to_string(),
                action: CommandPaletteAction::Redo,
            },
            CommandPaletteEntry {
                title: "Configuration".to_string(),
                detail: "Open global app configuration".to_string(),
                keywords: "settings config preferences".to_string(),
                action: CommandPaletteAction::OpenConfiguration,
            },
            CommandPaletteEntry {
                title: "Prepare Reference Genome".to_string(),
                detail: "Download/index selected genome".to_string(),
                keywords: "genome prepare reference".to_string(),
                action: CommandPaletteAction::OpenPrepareGenome,
            },
            CommandPaletteEntry {
                title: "Retrieve Genomic Sequence".to_string(),
                detail: "Extract region/gene from prepared genome".to_string(),
                keywords: "genome retrieve extract anchor".to_string(),
                action: CommandPaletteAction::OpenRetrieveGenome,
            },
            CommandPaletteEntry {
                title: "BLAST Genome Sequence".to_string(),
                detail: "Run BLAST against prepared genome indices".to_string(),
                keywords: "genome blast".to_string(),
                action: CommandPaletteAction::OpenBlastGenome,
            },
            CommandPaletteEntry {
                title: "Import Genome Tracks".to_string(),
                detail: "Import BED/BigWig/VCF onto anchored sequences".to_string(),
                keywords: "genome tracks bed bigwig vcf".to_string(),
                action: CommandPaletteAction::OpenGenomeTracks,
            },
            CommandPaletteEntry {
                title: "Agent Assistant".to_string(),
                detail: "Ask configured agent systems and run suggested commands".to_string(),
                keywords: "agent assistant chat automation support".to_string(),
                action: CommandPaletteAction::OpenAgentAssistant,
            },
            CommandPaletteEntry {
                title: "Prepared References".to_string(),
                detail: "Inspect prepared genome installs".to_string(),
                keywords: "genome inspector prepared".to_string(),
                action: CommandPaletteAction::OpenPreparedInspector,
            },
            CommandPaletteEntry {
                title: "GUI Manual".to_string(),
                detail: "Open in-app GUI help".to_string(),
                keywords: "help gui manual docs".to_string(),
                action: CommandPaletteAction::OpenGuiManual,
            },
            CommandPaletteEntry {
                title: "CLI Manual".to_string(),
                detail: "Open in-app CLI help".to_string(),
                keywords: "help cli manual docs".to_string(),
                action: CommandPaletteAction::OpenCliManual,
            },
            CommandPaletteEntry {
                title: "Shell Commands".to_string(),
                detail: "Open in-app shell command reference".to_string(),
                keywords: "help shell commands glossary docs".to_string(),
                action: CommandPaletteAction::OpenShellManual,
            },
            CommandPaletteEntry {
                title: "Export DALG SVG".to_string(),
                detail: "Export lineage graph as SVG".to_string(),
                keywords: "export lineage svg".to_string(),
                action: CommandPaletteAction::ExportLineageSvg,
            },
            CommandPaletteEntry {
                title: "Toggle Background Jobs Panel".to_string(),
                detail: "Show or hide centralized jobs panel".to_string(),
                keywords: "jobs progress background".to_string(),
                action: CommandPaletteAction::ToggleJobsPanel,
            },
            CommandPaletteEntry {
                title: "Toggle History Panel".to_string(),
                detail: "Show or hide operation history".to_string(),
                keywords: "history undo redo".to_string(),
                action: CommandPaletteAction::ToggleHistoryPanel,
            },
        ];
        for entry in self.collect_open_window_entries() {
            entries.push(CommandPaletteEntry {
                title: format!("Focus: {}", entry.title),
                detail: entry.detail,
                keywords: "focus window".to_string(),
                action: CommandPaletteAction::FocusViewport(entry.viewport_id),
            });
        }
        entries
    }

    fn execute_command_palette_action(
        &mut self,
        ctx: &egui::Context,
        action: CommandPaletteAction,
    ) {
        match action {
            CommandPaletteAction::NewProject => self.request_project_action(ProjectAction::New),
            CommandPaletteAction::OpenProject => self.request_project_action(ProjectAction::Open),
            CommandPaletteAction::SaveProject => {
                let _ = self.save_current_project();
            }
            CommandPaletteAction::OpenSequence => self.prompt_open_sequence(),
            CommandPaletteAction::OpenConfiguration => self.open_configuration_dialog(),
            CommandPaletteAction::OpenPrepareGenome => self.open_reference_genome_prepare_dialog(),
            CommandPaletteAction::OpenRetrieveGenome => {
                self.open_reference_genome_retrieve_dialog()
            }
            CommandPaletteAction::OpenBlastGenome => self.open_reference_genome_blast_dialog(),
            CommandPaletteAction::OpenGenomeTracks => self.open_genome_bed_track_dialog(),
            CommandPaletteAction::OpenAgentAssistant => self.open_agent_assistant_dialog(),
            CommandPaletteAction::OpenPreparedInspector => {
                self.open_reference_genome_inspector_dialog()
            }
            CommandPaletteAction::OpenGuiManual => self.open_help_doc(HelpDoc::Gui),
            CommandPaletteAction::OpenCliManual => self.open_help_doc(HelpDoc::Cli),
            CommandPaletteAction::OpenShellManual => self.open_help_doc(HelpDoc::Shell),
            CommandPaletteAction::ExportLineageSvg => self.prompt_export_lineage_svg(),
            CommandPaletteAction::ToggleJobsPanel => {
                self.show_jobs_panel = !self.show_jobs_panel;
            }
            CommandPaletteAction::ToggleHistoryPanel => {
                self.show_history_panel = !self.show_history_panel;
            }
            CommandPaletteAction::Undo => self.undo_last_operation(),
            CommandPaletteAction::Redo => self.redo_last_operation(),
            CommandPaletteAction::FocusViewport(viewport_id) => {
                self.focus_window_viewport(ctx, viewport_id);
            }
        }
    }

    fn open_help_doc(&mut self, doc: HelpDoc) {
        self.refresh_help_docs();
        self.help_doc = doc;
        self.refresh_help_search_matches();
        self.help_focus_search_box = true;
        self.help_image_preview_path = None;
        self.help_image_preview_caption.clear();
        self.show_help_dialog = true;
    }

    fn active_help_title(&self) -> &'static str {
        match self.help_doc {
            HelpDoc::Gui => "GUI Manual",
            HelpDoc::Cli => "CLI Manual",
            HelpDoc::Shell => "Shell Commands",
        }
    }

    fn active_help_markdown(&self) -> &str {
        match self.help_doc {
            HelpDoc::Gui => &self.help_gui_markdown,
            HelpDoc::Cli => &self.help_cli_markdown,
            HelpDoc::Shell => &self.help_shell_markdown,
        }
    }

    fn collect_help_markdown_images(markdown: &str) -> Vec<HelpMarkdownImage> {
        let mut images = Vec::new();
        let mut current: Option<HelpMarkdownImage> = None;
        for event in Parser::new(markdown) {
            match event {
                Event::Start(Tag::Image {
                    link_type: LinkType::Inline,
                    dest_url,
                    title,
                    ..
                }) => {
                    current = Some(HelpMarkdownImage {
                        path: dest_url.to_string(),
                        title: title.to_string(),
                        ..HelpMarkdownImage::default()
                    });
                }
                Event::End(pulldown_cmark::TagEnd::Image) => {
                    if let Some(image) = current.take() {
                        if !image.path.trim().is_empty() {
                            images.push(image);
                        }
                    }
                }
                Event::Text(text) | Event::Code(text) => {
                    if let Some(image) = current.as_mut() {
                        if !image.alt.is_empty() {
                            image.alt.push(' ');
                        }
                        image.alt.push_str(text.as_ref());
                    }
                }
                Event::SoftBreak | Event::HardBreak => {
                    if let Some(image) = current.as_mut() {
                        if !image.alt.is_empty() {
                            image.alt.push(' ');
                        }
                    }
                }
                _ => {}
            }
        }
        images
    }

    fn help_image_caption(image: &HelpMarkdownImage) -> String {
        if !image.title.trim().is_empty() {
            image.title.trim().to_string()
        } else if !image.alt.trim().is_empty() {
            image.alt.trim().to_string()
        } else {
            Path::new(image.path.as_str())
                .file_name()
                .and_then(|name| name.to_str())
                .unwrap_or("Image")
                .to_string()
        }
    }

    fn refresh_help_search_matches(&mut self) {
        self.help_search_matches.clear();
        let needle = self.help_search_query.trim();
        if needle.is_empty() {
            self.help_search_selected = 0;
            return;
        }
        let needle_lower = needle.to_ascii_lowercase();
        let markdown = self.active_help_markdown().to_string();
        for (idx, line) in markdown.lines().enumerate() {
            if line.to_ascii_lowercase().contains(&needle_lower) {
                let snippet = line.trim();
                self.help_search_matches.push(HelpSearchMatch {
                    line_number: idx + 1,
                    snippet: if snippet.is_empty() {
                        "(empty line)".to_string()
                    } else {
                        snippet.chars().take(140).collect()
                    },
                });
            }
        }
        if self.help_search_matches.is_empty() {
            self.help_search_selected = 0;
        } else if self.help_search_selected >= self.help_search_matches.len() {
            self.help_search_selected = self.help_search_matches.len() - 1;
        }
    }

    fn select_next_help_match(&mut self) {
        if self.help_search_matches.is_empty() {
            return;
        }
        self.help_search_selected =
            (self.help_search_selected + 1) % self.help_search_matches.len();
    }

    fn select_prev_help_match(&mut self) {
        if self.help_search_matches.is_empty() {
            return;
        }
        if self.help_search_selected == 0 {
            self.help_search_selected = self.help_search_matches.len() - 1;
        } else {
            self.help_search_selected -= 1;
        }
    }

    pub fn new() -> Self {
        let mut app = Self::default();
        window_backdrop::set_window_backdrop_settings(app.window_backdrops.clone());
        app.refresh_help_docs();
        app.load_persisted_configuration(true);
        app.load_bed_track_subscriptions_from_state();
        app.mark_clean_snapshot();
        app
    }

    pub fn new_with_project(project_path: Option<&str>) -> Self {
        let mut ret = Self::new();
        if let Some(path) = project_path {
            let _ = ret.load_project_from_file(path);
        }
        ret
    }

    fn load_dna_from_genbank_file(filename: &str) -> Result<DNAsequence> {
        let dna = dna_sequence::DNAsequence::from_genbank_file(filename)?
            .pop()
            .ok_or_else(|| anyhow!("Could not read GenBank file {filename}"))?;
        Ok(dna)
    }

    fn load_dna_from_embl_file(filename: &str) -> Result<DNAsequence> {
        let dna = dna_sequence::DNAsequence::from_embl_file(filename)?
            .pop()
            .ok_or_else(|| anyhow!("Could not read EMBL file {filename}"))?;
        Ok(dna)
    }

    fn load_dna_from_fasta_file(filename: &str) -> Result<DNAsequence> {
        let dna = dna_sequence::DNAsequence::from_fasta_file(filename)?
            .pop()
            .ok_or_else(|| anyhow!("Could not read fasta file {filename}"))?;
        Ok(dna)
    }

    fn new_dna_window(&mut self, seq_id: String, dna: DNAsequence) {
        self.new_windows
            .push(Window::new_dna(dna, seq_id, self.engine.clone()));
    }

    fn find_open_sequence_viewport_id(&self, seq_id: &str) -> Option<ViewportId> {
        self.windows.iter().find_map(|(viewport_id, window)| {
            let Ok(guard) = window.read() else {
                return None;
            };
            if guard.sequence_id().as_deref() == Some(seq_id) {
                Some(*viewport_id)
            } else {
                None
            }
        })
    }

    fn find_pending_sequence_window_mut(&mut self, seq_id: &str) -> Option<&mut Window> {
        self.new_windows
            .iter_mut()
            .find(|window| window.sequence_id().as_deref() == Some(seq_id))
    }

    fn queue_focus_viewport(&mut self, viewport_id: ViewportId) {
        if !self.pending_focus_viewports.contains(&viewport_id) {
            self.pending_focus_viewports.push(viewport_id);
        }
    }

    fn open_sequence_window(&mut self, seq_id: &str) {
        if let Some(viewport_id) = self.find_open_sequence_viewport_id(seq_id) {
            self.queue_focus_viewport(viewport_id);
            return;
        }
        if self.find_pending_sequence_window_mut(seq_id).is_some() {
            return;
        }
        let dna = self
            .engine
            .read()
            .unwrap()
            .state()
            .sequences
            .get(seq_id)
            .cloned();
        if let Some(dna) = dna {
            self.new_dna_window(seq_id.to_string(), dna);
        }
    }

    fn open_pool_window(&mut self, representative_seq_id: &str, pool_seq_ids: Vec<String>) {
        if let Some(viewport_id) = self.find_open_sequence_viewport_id(representative_seq_id) {
            if let Some(window) = self.windows.get(&viewport_id) {
                if let Ok(mut window) = window.write() {
                    window.set_pool_context(pool_seq_ids);
                }
            }
            self.queue_focus_viewport(viewport_id);
            return;
        }
        if let Some(window) = self.find_pending_sequence_window_mut(representative_seq_id) {
            window.set_pool_context(pool_seq_ids);
            return;
        }
        let dna = self
            .engine
            .read()
            .unwrap()
            .state()
            .sequences
            .get(representative_seq_id)
            .cloned();
        if let Some(dna) = dna {
            let mut window =
                Window::new_dna(dna, representative_seq_id.to_string(), self.engine.clone());
            window.set_pool_context(pool_seq_ids);
            self.new_windows.push(window);
        }
    }

    pub fn load_from_file(path: &str) -> Result<DNAsequence> {
        let extension = Path::new(path)
            .extension()
            .and_then(|value| value.to_str())
            .map(|value| value.to_ascii_lowercase());
        match extension.as_deref() {
            Some("embl") | Some("emb") => return Self::load_dna_from_embl_file(path),
            Some("gb") | Some("gbk") | Some("genbank") => {
                return Self::load_dna_from_genbank_file(path);
            }
            Some("fa") | Some("fasta") | Some("fna") | Some("fas") => {
                return Self::load_dna_from_fasta_file(path);
            }
            _ => {}
        }

        for loader in [
            Self::load_dna_from_genbank_file,
            Self::load_dna_from_embl_file,
            Self::load_dna_from_fasta_file,
        ] {
            if let Ok(dna) = loader(path) {
                return Ok(dna);
            }
        }
        Err(anyhow!("Could not load file '{path}'"))
    }

    fn open_new_window_from_file(&mut self, path: &str) {
        let op = Operation::LoadFile {
            path: path.to_string(),
            as_id: None,
        };
        let load_result = {
            let mut engine = self.engine.write().unwrap();
            engine.apply(op)
        };
        if let Ok(result) = load_result {
            if let Some(seq_id) = result.created_seq_ids.first() {
                let seq_id = seq_id.to_string();
                let dna = self
                    .engine
                    .read()
                    .unwrap()
                    .state()
                    .sequences
                    .get(&seq_id)
                    .cloned();
                if let Some(dna) = dna {
                    self.new_dna_window(seq_id, dna);
                    return;
                }
            }
            // TODO warning: load op ran but no created sequence
        } else {
            // TODO error, could not load file through engine
        }
    }

    fn save_project_to_file(&self, path: &str) -> Result<()> {
        self.engine
            .read()
            .unwrap()
            .state()
            .save_to_path(path)
            .map_err(|e| anyhow!(e.to_string()))
    }

    fn current_state_fingerprint(&self) -> u64 {
        let state_json = {
            let engine = self.engine.read().unwrap();
            serde_json::to_vec(engine.state()).unwrap_or_default()
        };
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        state_json.hash(&mut hasher);
        hasher.finish()
    }

    fn current_state_change_stamp(&self) -> u64 {
        let engine = self.engine.read().unwrap();
        let state = engine.state();
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        engine.operation_log().len().hash(&mut hasher);
        state.sequences.len().hash(&mut hasher);
        state.lineage.nodes.len().hash(&mut hasher);
        state.lineage.edges.len().hash(&mut hasher);
        state.container_state.containers.len().hash(&mut hasher);
        state.container_state.arrangements.len().hash(&mut hasher);
        state.metadata.len().hash(&mut hasher);

        let mut metadata_keys: Vec<&String> = state.metadata.keys().collect();
        metadata_keys.sort_unstable();
        for key in metadata_keys {
            key.hash(&mut hasher);
        }

        let display = &state.display;
        display.show_sequence_panel.hash(&mut hasher);
        display.show_map_panel.hash(&mut hasher);
        display.show_cds_features.hash(&mut hasher);
        display.show_gene_features.hash(&mut hasher);
        display.show_mrna_features.hash(&mut hasher);
        display.show_tfbs.hash(&mut hasher);
        display.regulatory_tracks_near_baseline.hash(&mut hasher);
        display
            .regulatory_feature_max_view_span_bp
            .hash(&mut hasher);
        display.tfbs_display_use_llr_bits.hash(&mut hasher);
        display
            .tfbs_display_min_llr_bits
            .to_bits()
            .hash(&mut hasher);
        display.tfbs_display_use_llr_quantile.hash(&mut hasher);
        display
            .tfbs_display_min_llr_quantile
            .to_bits()
            .hash(&mut hasher);
        display
            .tfbs_display_use_true_log_odds_bits
            .hash(&mut hasher);
        display
            .tfbs_display_min_true_log_odds_bits
            .to_bits()
            .hash(&mut hasher);
        display
            .tfbs_display_use_true_log_odds_quantile
            .hash(&mut hasher);
        display
            .tfbs_display_min_true_log_odds_quantile
            .to_bits()
            .hash(&mut hasher);
        display.vcf_display_show_snp.hash(&mut hasher);
        display.vcf_display_show_ins.hash(&mut hasher);
        display.vcf_display_show_del.hash(&mut hasher);
        display.vcf_display_show_sv.hash(&mut hasher);
        display.vcf_display_show_other.hash(&mut hasher);
        display.vcf_display_pass_only.hash(&mut hasher);
        display.vcf_display_use_min_qual.hash(&mut hasher);
        display.vcf_display_min_qual.to_bits().hash(&mut hasher);
        display.vcf_display_use_max_qual.hash(&mut hasher);
        display.vcf_display_max_qual.to_bits().hash(&mut hasher);
        for key in &display.vcf_display_required_info_keys {
            key.hash(&mut hasher);
        }
        display.show_restriction_enzymes.hash(&mut hasher);
        display.show_gc_contents.hash(&mut hasher);
        display.show_open_reading_frames.hash(&mut hasher);
        display.show_methylation_sites.hash(&mut hasher);
        display
            .feature_details_font_size
            .to_bits()
            .hash(&mut hasher);
        display.linear_view_start_bp.hash(&mut hasher);
        display.linear_view_span_bp.hash(&mut hasher);

        hasher.finish()
    }

    fn current_lineage_change_stamp(&self) -> u64 {
        let engine = self.engine.read().unwrap();
        let state = engine.state();
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        engine.operation_log().len().hash(&mut hasher);
        state.sequences.len().hash(&mut hasher);
        state.lineage.nodes.len().hash(&mut hasher);
        state.lineage.edges.len().hash(&mut hasher);
        state.container_state.containers.len().hash(&mut hasher);
        state.container_state.arrangements.len().hash(&mut hasher);
        hasher.finish()
    }

    fn current_operation_count(&self) -> usize {
        self.engine.read().unwrap().operation_log().len()
    }

    fn project_has_user_content(&self) -> bool {
        let engine = self.engine.read().unwrap();
        let state = engine.state();
        !state.sequences.is_empty()
            || !state.lineage.nodes.is_empty()
            || !state.container_state.containers.is_empty()
            || !state.container_state.arrangements.is_empty()
    }

    fn mark_clean_snapshot(&mut self) {
        self.clean_state_fingerprint = self.current_state_fingerprint();
        self.dirty_cache_stamp = self.current_state_change_stamp();
        self.dirty_cache_valid = true;
        self.dirty_cache_value = false;
        self.dirty_cache_last_deep_check = Instant::now();
    }

    fn is_project_dirty(&mut self) -> bool {
        if !self.project_has_user_content() {
            return false;
        }
        const DIRTY_DEEP_CHECK_INTERVAL: Duration = Duration::from_secs(2);
        let now = Instant::now();
        let stamp = self.current_state_change_stamp();
        let stamp_changed = !self.dirty_cache_valid || stamp != self.dirty_cache_stamp;
        if stamp_changed
            || now.duration_since(self.dirty_cache_last_deep_check) >= DIRTY_DEEP_CHECK_INTERVAL
        {
            self.dirty_cache_value =
                self.current_state_fingerprint() != self.clean_state_fingerprint;
            self.dirty_cache_last_deep_check = now;
        }
        self.dirty_cache_stamp = stamp;
        self.dirty_cache_valid = true;
        self.dirty_cache_value
    }

    fn reset_to_empty_project(&mut self) {
        self.engine = Arc::new(RwLock::new(GentleEngine::new()));
        self.apply_configuration_graphics_to_engine_state();
        self.current_project_path = None;
        self.last_applied_window_title.clear();
        self.lineage_cache_valid = false;
        self.lineage_rows.clear();
        self.lineage_edges.clear();
        self.lineage_op_label_by_id.clear();
        self.lineage_containers.clear();
        self.lineage_arrangements.clear();
        self.lineage_graph_zoom = 1.0;
        self.lineage_graph_area_height = 420.0;
        self.lineage_container_area_height = 220.0;
        self.lineage_graph_scroll_offset = Vec2::ZERO;
        self.lineage_graph_pan_origin = None;
        self.lineage_graph_compact_labels = true;
        self.lineage_graph_selected_node_id = None;
        self.lineage_graph_node_offsets.clear();
        self.lineage_graph_drag_origin = None;
        self.lineage_graph_offsets_synced_stamp = 0;
        self.new_windows.clear();
        self.windows.clear();
        self.windows_to_close.write().unwrap().clear();
        self.pending_focus_viewports.clear();
        self.viewport_id_counter = 0;
        self.pending_project_action = None;
        self.show_reference_genome_blast_dialog = false;
        self.genome_blast_task = None;
        self.genome_blast_progress_fraction = None;
        self.genome_blast_progress_label.clear();
        self.genome_blast_results.clear();
        self.genome_blast_selected_result = 0;
        self.genome_blast_status.clear();
        self.show_genome_bed_track_dialog = false;
        self.show_agent_assistant_dialog = false;
        self.genome_track_status.clear();
        self.genome_track_import_task = None;
        self.genome_track_import_progress = None;
        self.genome_track_autosync_status.clear();
        self.tracked_autosync_last_op_count = None;
        self.agent_task = None;
        self.agent_model_discovery_task = None;
        self.agent_status.clear();
        self.agent_last_invocation = None;
        self.agent_execution_log.clear();
        self.agent_discovered_models.clear();
        self.agent_discovered_model_pick.clear();
        self.agent_model_discovery_status.clear();
        self.agent_model_discovery_source_key.clear();
        self.load_bed_track_subscriptions_from_state();
        self.load_lineage_graph_workspace_from_state();
        self.mark_clean_snapshot();
    }

    fn set_current_project_path_and_track_recent(&mut self, path: &str) {
        let normalized = Self::normalize_project_path(path);
        if normalized.is_empty() {
            self.current_project_path = None;
            return;
        }
        self.current_project_path = Some(normalized.clone());
        self.record_recent_project_path(&normalized);
    }

    fn open_project_path(&mut self, path: &str) {
        if let Err(err) = self.load_project_from_file(path) {
            eprintln!("Could not open project '{}': {err}", path);
            self.remove_recent_project_path(path);
        }
    }

    fn prompt_open_project(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("GENtle project", &["json"])
            .pick_file()
        {
            let path = path.display().to_string();
            self.open_project_path(&path);
        }
    }

    fn prompt_open_sequence(&mut self) {
        if let Some(path) = rfd::FileDialog::new().pick_file() {
            let path = path.display().to_string();
            self.open_new_window_from_file(&path);
        }
    }

    fn prompt_save_project(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .set_file_name("project.gentle.json")
            .add_filter("GENtle project", &["json"])
            .save_file()
        {
            let path = path.display().to_string();
            if self.save_project_to_file(&path).is_ok() {
                self.set_current_project_path_and_track_recent(&path);
                self.mark_clean_snapshot();
            }
        }
    }

    fn prompt_export_lineage_svg(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .set_file_name("lineage.svg")
            .add_filter("SVG", &["svg"])
            .save_file()
        {
            let path = path.display().to_string();
            let result = self
                .engine
                .write()
                .unwrap()
                .apply(Operation::RenderLineageSvg { path: path.clone() });
            if let Err(e) = result {
                eprintln!("Could not export lineage SVG to '{}': {}", path, e.message);
            }
        }
    }

    fn sanitize_file_stem(raw: &str, fallback: &str) -> String {
        let mut out: String = raw
            .trim()
            .chars()
            .map(|ch| {
                if ch.is_ascii_alphanumeric() || ch == '-' || ch == '_' {
                    ch
                } else {
                    '_'
                }
            })
            .collect();
        out = out.trim_matches('_').to_string();
        if out.is_empty() {
            fallback.to_string()
        } else {
            out
        }
    }

    fn prompt_export_serial_gel_svg(
        &mut self,
        default_stem: &str,
        container_ids: Option<Vec<String>>,
        arrangement_id: Option<String>,
        ladders: Option<Vec<String>>,
    ) {
        let stem = Self::sanitize_file_stem(default_stem, "serial_gel");
        let default_file_name = format!("{stem}.gel.svg");
        let path = rfd::FileDialog::new()
            .set_file_name(&default_file_name)
            .add_filter("SVG", &["svg"])
            .save_file();
        let Some(path) = path else {
            self.app_status = "Serial gel SVG export canceled".to_string();
            return;
        };
        let path_text = path.display().to_string();
        let result = self
            .engine
            .write()
            .unwrap()
            .apply(Operation::RenderPoolGelSvg {
                inputs: vec![],
                path: path_text.clone(),
                ladders,
                container_ids,
                arrangement_id,
            });
        match result {
            Ok(op_result) => {
                self.app_status = op_result
                    .messages
                    .first()
                    .cloned()
                    .unwrap_or_else(|| format!("Wrote serial gel SVG to '{path_text}'"));
            }
            Err(e) => {
                self.app_status = format!("Could not export serial gel SVG: {}", e.message);
            }
        }
    }

    fn open_reference_genome_prepare_dialog(&mut self) {
        if self.show_reference_genome_prepare_dialog {
            self.queue_focus_viewport(Self::prepare_genome_viewport_id());
        }
        self.show_reference_genome_prepare_dialog = true;
    }

    fn open_reference_genome_retrieve_dialog(&mut self) {
        self.show_reference_genome_retrieve_dialog = true;
    }

    fn open_reference_genome_blast_dialog(&mut self) {
        if self.show_reference_genome_blast_dialog {
            self.queue_focus_viewport(Self::blast_genome_viewport_id());
        }
        self.show_reference_genome_blast_dialog = true;
    }

    fn open_reference_genome_inspector_dialog(&mut self) {
        self.show_reference_genome_inspector_dialog = true;
    }

    fn open_genome_bed_track_dialog(&mut self) {
        self.load_bed_track_subscriptions_from_state();
        self.show_genome_bed_track_dialog = true;
    }

    fn open_agent_assistant_dialog(&mut self) {
        self.refresh_agent_system_catalog();
        self.show_agent_assistant_dialog = true;
    }

    fn track_name_default_from_path(path: &Path) -> String {
        path.file_name()
            .map(|value| value.to_string_lossy().to_string())
            .unwrap_or_else(|| path.display().to_string())
    }

    fn open_helper_genome_prepare_dialog(&mut self) {
        self.genome_catalog_path = DEFAULT_HELPER_GENOME_CATALOG_PATH.to_string();
        self.genome_cache_dir = "data/helper_genomes".to_string();
        self.invalidate_genome_genes();
        self.open_reference_genome_prepare_dialog();
    }

    fn open_helper_genome_retrieve_dialog(&mut self) {
        self.genome_catalog_path = DEFAULT_HELPER_GENOME_CATALOG_PATH.to_string();
        self.genome_cache_dir = "data/helper_genomes".to_string();
        self.invalidate_genome_genes();
        self.show_reference_genome_retrieve_dialog = true;
    }

    fn open_helper_genome_blast_dialog(&mut self) {
        self.genome_catalog_path = DEFAULT_HELPER_GENOME_CATALOG_PATH.to_string();
        self.genome_cache_dir = "data/helper_genomes".to_string();
        self.open_reference_genome_blast_dialog();
    }

    fn genome_catalog_path_opt(&self) -> Option<String> {
        let trimmed = self.genome_catalog_path.trim();
        if trimmed.is_empty() {
            None
        } else {
            Some(trimmed.to_string())
        }
    }

    fn genome_catalog_path_resolved(&self) -> String {
        self.genome_catalog_path_opt()
            .unwrap_or_else(|| DEFAULT_GENOME_CATALOG_PATH.to_string())
    }

    fn genome_cache_dir_opt(&self) -> Option<String> {
        let trimmed = self.genome_cache_dir.trim();
        if trimmed.is_empty() {
            None
        } else {
            Some(trimmed.to_string())
        }
    }

    fn anchored_sequence_anchor_summaries_for_tracks(&self) -> Vec<SequenceGenomeAnchorSummary> {
        self.engine
            .read()
            .unwrap()
            .list_sequence_genome_anchor_summaries()
    }

    fn describe_sequence_genome_anchor(&self, seq_id: &str) -> Option<String> {
        self.engine
            .read()
            .unwrap()
            .describe_sequence_genome_anchor(seq_id)
            .ok()
    }

    fn load_bed_track_subscriptions_from_state(&mut self) {
        let subscriptions = self
            .engine
            .read()
            .unwrap()
            .list_genome_track_subscriptions();
        self.genome_bed_track_subscriptions = subscriptions;
        self.tracked_autosync_last_op_count = None;
    }

    fn parse_optional_score_field(raw: &str, label: &str) -> Result<Option<f64>> {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            Ok(None)
        } else {
            trimmed
                .parse::<f64>()
                .map(Some)
                .map_err(|e| anyhow!("Invalid {label} value '{trimmed}': {e}"))
        }
    }

    fn anchors_share_mapping_group(
        left: &SequenceGenomeAnchorSummary,
        right: &SequenceGenomeAnchorSummary,
    ) -> bool {
        if left.genome_id != right.genome_id || left.chromosome != right.chromosome {
            return false;
        }
        match (left.strand, right.strand) {
            (Some(a), Some(b)) => a == b,
            _ => true,
        }
    }

    fn parse_bed_track_form(&self) -> Result<GenomeTrackSubscription> {
        let path = self.genome_track_path.trim().to_string();
        if path.is_empty() {
            return Err(anyhow!("Track path cannot be empty"));
        }
        let min_score =
            Self::parse_optional_score_field(&self.genome_track_min_score, "min score")?;
        let max_score =
            Self::parse_optional_score_field(&self.genome_track_max_score, "max score")?;
        if min_score
            .zip(max_score)
            .map(|(min, max)| min > max)
            .unwrap_or(false)
        {
            return Err(anyhow!("min score must be <= max score"));
        }
        let track_name = {
            let trimmed = self.genome_track_name.trim();
            if trimmed.is_empty() {
                None
            } else {
                Some(trimmed.to_string())
            }
        };
        Ok(GenomeTrackSubscription {
            source: self.genome_track_source_selection.resolve(&path),
            path,
            track_name,
            min_score,
            max_score,
            clear_existing: self.genome_track_clear_existing,
        })
    }

    fn genome_track_import_operation(
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

    fn start_genome_track_import_for_selected_sequence(
        &mut self,
        seq_id: String,
        subscription: GenomeTrackSubscription,
    ) {
        if self.genome_track_import_task.is_some() {
            self.genome_track_status = "Another genome track import is already running".to_string();
            return;
        }
        let op = Self::genome_track_import_operation(&seq_id, &subscription);
        let source_label = subscription.source.label().to_string();
        let path_label = subscription.path.clone();
        let job_id = self.alloc_background_job_id();
        let (tx, rx) = mpsc::channel::<GenomeTrackImportTaskMessage>();
        let cancel_requested = Arc::new(AtomicBool::new(false));
        self.genome_track_import_progress = None;
        self.genome_track_status = format!(
            "Importing {source_label} track in background: '{}' -> '{}'",
            path_label, seq_id
        );
        self.push_job_event(
            BackgroundJobKind::TrackImport,
            BackgroundJobEventPhase::Started,
            Some(job_id),
            format!(
                "Track import started: {} '{}' -> {}",
                source_label, path_label, seq_id
            ),
        );
        self.genome_track_import_task = Some(GenomeTrackImportTask {
            job_id,
            started: Instant::now(),
            cancel_requested: cancel_requested.clone(),
            receiver: rx,
        });

        let engine = self.engine.clone();
        std::thread::spawn(move || {
            let tx_progress = tx.clone();
            let cancel_flag = cancel_requested.clone();
            let outcome = {
                let mut guard = engine.write().expect("Engine lock poisoned");
                guard.apply_with_progress(op, move |progress| match progress {
                    OperationProgress::GenomeTrackImport(p) => {
                        let _ = tx_progress.send(GenomeTrackImportTaskMessage::Progress {
                            job_id,
                            progress: p,
                        });
                        !cancel_flag.load(Ordering::Relaxed)
                    }
                    _ => true,
                })
            };
            let _ = tx.send(GenomeTrackImportTaskMessage::Done {
                job_id,
                result: outcome,
            });
        });
    }

    fn validate_bigwig_converter_available(&self) -> Result<String> {
        let executable = self.resolved_bigwig_to_bedgraph_executable();
        match Command::new(&executable).arg("-version").output() {
            Ok(output) => {
                let detail = Self::first_non_empty_output_line(&output.stdout, &output.stderr);
                Ok(format!("{executable}: {detail}"))
            }
            Err(e) if e.kind() == ErrorKind::NotFound => Err(anyhow!(
                "bigWigToBedGraph executable '{}' not found",
                executable
            )),
            Err(e) => Err(anyhow!("Could not execute '{}': {}", executable, e)),
        }
    }

    fn ensure_bigwig_converter_ready(&mut self, subscription: &GenomeTrackSubscription) -> bool {
        if subscription.source != GenomeTrackSource::BigWig {
            return true;
        }
        match self.validate_bigwig_converter_available() {
            Ok(detail) => {
                self.genome_track_status = format!("BigWig converter check: {detail}");
                true
            }
            Err(e) => {
                self.genome_track_status = format!("BigWig import preflight failed: {e}");
                false
            }
        }
    }

    fn format_track_sync_status(prefix: &str, report: &GenomeTrackSyncReport) -> String {
        let mut text = format!(
            "{prefix}: subscriptions={}, targets={}, applied={}, failed={}, warnings={}",
            report.subscriptions_considered,
            report.target_sequences,
            report.applied_imports,
            report.failed_imports,
            report.warnings_count
        );
        if !report.errors.is_empty() {
            text.push_str(" | errors: ");
            text.push_str(
                &report
                    .errors
                    .iter()
                    .take(3)
                    .cloned()
                    .collect::<Vec<_>>()
                    .join(" | "),
            );
        }
        text
    }

    fn import_genome_bed_track_for_selected_sequence(&mut self) {
        if self.genome_track_import_task.is_some() {
            self.genome_track_status = "Another genome track import is already running".to_string();
            return;
        }
        let seq_id = self.genome_track_seq_id.trim().to_string();
        if seq_id.is_empty() {
            self.genome_track_status =
                "Select a genome-anchored sequence before importing a track".to_string();
            return;
        }
        let subscription = match self.parse_bed_track_form() {
            Ok(value) => value,
            Err(e) => {
                self.genome_track_status = format!("Import track failed: {e}");
                return;
            }
        };
        if !self.ensure_bigwig_converter_ready(&subscription) {
            return;
        }
        self.start_genome_track_import_for_selected_sequence(seq_id, subscription);
    }

    fn import_genome_bed_track_for_all_anchored_sequences(&mut self, track_subscription: bool) {
        if self.genome_track_import_task.is_some() {
            self.genome_track_status = "Another genome track import is already running".to_string();
            return;
        }
        let subscription = match self.parse_bed_track_form() {
            Ok(value) => value,
            Err(e) => {
                self.genome_track_status = format!("Import track failed: {e}");
                return;
            }
        };
        if !self.ensure_bigwig_converter_ready(&subscription) {
            return;
        }

        let sync_result = self
            .engine
            .write()
            .unwrap()
            .import_genome_track_to_all_anchored(subscription.clone(), track_subscription);
        match sync_result {
            Ok(report) => {
                let tracked_note = if track_subscription {
                    "tracked for auto-sync"
                } else {
                    "one-time import"
                };
                self.genome_track_status = format!(
                    "{} ({})",
                    Self::format_track_sync_status(
                        &format!("{} import (all anchored)", subscription.source.label()),
                        &report
                    ),
                    tracked_note
                );
                self.push_job_event(
                    BackgroundJobKind::TrackImport,
                    BackgroundJobEventPhase::Completed,
                    None,
                    format!(
                        "{} import (all anchored): applied={}, failed={}",
                        subscription.source.label(),
                        report.applied_imports,
                        report.failed_imports
                    ),
                );
                self.load_bed_track_subscriptions_from_state();
            }
            Err(e) => {
                self.genome_track_status = format!("Import track failed: {}", e.message);
                self.push_job_event(
                    BackgroundJobKind::TrackImport,
                    BackgroundJobEventPhase::Failed,
                    None,
                    format!(
                        "{} import (all anchored) failed: {}",
                        subscription.source.label(),
                        e.message
                    ),
                );
            }
        }
    }

    fn apply_tracked_bed_subscription_to_all_anchored(&mut self, index: usize) {
        if self.genome_track_import_task.is_some() {
            self.genome_track_status = "Another genome track import is already running".to_string();
            return;
        }
        let Some(subscription) = self.genome_bed_track_subscriptions.get(index).cloned() else {
            return;
        };
        if !self.ensure_bigwig_converter_ready(&subscription) {
            return;
        }

        let sync_result = self
            .engine
            .write()
            .unwrap()
            .apply_tracked_genome_track_subscription(index);
        match sync_result {
            Ok(report) => {
                self.genome_track_status = Self::format_track_sync_status(
                    &format!(
                        "Re-applied tracked {} '{}'",
                        subscription.source.label(),
                        subscription.path
                    ),
                    &report,
                );
                self.push_job_event(
                    BackgroundJobKind::TrackImport,
                    BackgroundJobEventPhase::Completed,
                    None,
                    format!(
                        "Re-applied tracked {} '{}': applied={}, failed={}",
                        subscription.source.label(),
                        subscription.path,
                        report.applied_imports,
                        report.failed_imports
                    ),
                );
            }
            Err(e) => {
                self.genome_track_status =
                    format!("Could not re-apply tracked subscription: {}", e.message);
                self.push_job_event(
                    BackgroundJobKind::TrackImport,
                    BackgroundJobEventPhase::Failed,
                    None,
                    format!("Tracked subscription apply failed: {}", e.message),
                );
            }
        }
    }

    fn sync_tracked_bed_tracks_for_new_anchors(&mut self) {
        if self.genome_track_import_task.is_some() {
            return;
        }
        if self.genome_bed_track_subscriptions.is_empty() {
            self.tracked_autosync_last_op_count = Some(self.current_operation_count());
            return;
        }
        let operation_count = self.current_operation_count();
        if self.tracked_autosync_last_op_count == Some(operation_count) {
            return;
        }
        self.tracked_autosync_last_op_count = Some(operation_count);
        let sync_result = {
            self.engine
                .write()
                .unwrap()
                .sync_tracked_genome_track_subscriptions(true)
        };
        match sync_result {
            Ok(report) => {
                if report.applied_imports > 0 || report.failed_imports > 0 {
                    self.genome_track_autosync_status =
                        Self::format_track_sync_status("Auto-sync", &report);
                    self.push_job_event(
                        BackgroundJobKind::TrackImport,
                        BackgroundJobEventPhase::Completed,
                        None,
                        format!(
                            "Auto-sync: applied={}, failed={}",
                            report.applied_imports, report.failed_imports
                        ),
                    );
                }
                self.tracked_autosync_last_op_count = Some(self.current_operation_count());
            }
            Err(e) => {
                self.genome_track_autosync_status =
                    format!("Auto-sync failed unexpectedly: {}", e.message);
                self.push_job_event(
                    BackgroundJobKind::TrackImport,
                    BackgroundJobEventPhase::Failed,
                    None,
                    format!("Auto-sync failed: {}", e.message),
                );
                self.tracked_autosync_last_op_count = Some(self.current_operation_count());
            }
        }
    }

    fn refresh_genome_catalog_list(&mut self) {
        let catalog_path = self.genome_catalog_path_resolved();
        let prev_genome = self.genome_id.clone();
        match GenomeCatalog::from_json_file(&catalog_path) {
            Ok(catalog) => {
                self.genome_catalog_genomes = catalog.list_genomes();
                self.genome_catalog_error.clear();
                if self.genome_catalog_genomes.is_empty() {
                    self.genome_id.clear();
                } else if !self.genome_catalog_genomes.contains(&self.genome_id) {
                    self.genome_id = self.genome_catalog_genomes[0].clone();
                }
            }
            Err(e) => {
                self.genome_catalog_genomes.clear();
                self.genome_catalog_error = e;
                self.genome_id.clear();
            }
        }
        if self.genome_id != prev_genome {
            self.invalidate_genome_genes();
        }
    }

    fn refresh_agent_system_catalog(&mut self) {
        let catalog_path = self.agent_catalog_path.trim().to_string();
        if !self.agent_systems.is_empty()
            && self.agent_catalog_loaded_path == catalog_path
            && self.agent_catalog_error.is_empty()
        {
            return;
        }
        self.agent_catalog_loaded_path = catalog_path.clone();
        match load_agent_system_catalog(Some(&catalog_path)) {
            Ok((_resolved, catalog)) => {
                self.agent_systems = catalog.systems;
                self.agent_catalog_error.clear();
                if self.agent_system_id.trim().is_empty()
                    || !self
                        .agent_systems
                        .iter()
                        .any(|system| system.id == self.agent_system_id)
                {
                    self.agent_system_id = self
                        .agent_systems
                        .first()
                        .map(|system| system.id.clone())
                        .unwrap_or_default();
                }
            }
            Err(err) => {
                self.agent_catalog_error = err;
                self.agent_systems.clear();
                self.agent_system_id.clear();
            }
        }
    }

    fn selected_agent_system(&self) -> Option<AgentSystemSpec> {
        self.agent_systems
            .iter()
            .find(|system| system.id == self.agent_system_id)
            .cloned()
    }

    fn selected_agent_system_with_session_overrides(
        &self,
        system: &AgentSystemSpec,
    ) -> AgentSystemSpec {
        let mut resolved = system.clone();
        let override_base_url = self.agent_base_url_override.trim();
        if !override_base_url.is_empty()
            && matches!(
                resolved.transport,
                AgentSystemTransport::NativeOpenai | AgentSystemTransport::NativeOpenaiCompat
            )
        {
            resolved.env.insert(
                AGENT_BASE_URL_ENV.to_string(),
                override_base_url.to_string(),
            );
        }
        let selected_discovered_model =
            normalize_agent_model_name(&self.agent_discovered_model_pick).filter(|picked| {
                self.agent_discovered_models
                    .iter()
                    .any(|item| item == picked)
            });
        let override_model = normalize_agent_model_name(self.agent_model_override.trim())
            .or(selected_discovered_model);
        if let Some(override_model) = override_model
            && matches!(
                resolved.transport,
                AgentSystemTransport::NativeOpenai | AgentSystemTransport::NativeOpenaiCompat
            )
        {
            resolved
                .env
                .insert(AGENT_MODEL_ENV.to_string(), override_model);
        }
        let timeout_override = self
            .agent_timeout_secs
            .trim()
            .parse::<u64>()
            .ok()
            .filter(|value| *value > 0);
        if let Some(timeout_override) = timeout_override
            && matches!(
                resolved.transport,
                AgentSystemTransport::ExternalJsonStdio
                    | AgentSystemTransport::NativeOpenai
                    | AgentSystemTransport::NativeOpenaiCompat
            )
        {
            resolved.env.insert(
                AGENT_TIMEOUT_SECS_ENV.to_string(),
                timeout_override.to_string(),
            );
        }
        let connect_timeout_override = self
            .agent_connect_timeout_secs
            .trim()
            .parse::<u64>()
            .ok()
            .filter(|value| *value > 0);
        if let Some(connect_timeout_override) = connect_timeout_override
            && matches!(
                resolved.transport,
                AgentSystemTransport::NativeOpenai | AgentSystemTransport::NativeOpenaiCompat
            )
        {
            resolved.env.insert(
                AGENT_CONNECT_TIMEOUT_SECS_ENV.to_string(),
                connect_timeout_override.to_string(),
            );
        }
        let read_timeout_override = self
            .agent_read_timeout_secs
            .trim()
            .parse::<u64>()
            .ok()
            .filter(|value| *value > 0);
        if let Some(read_timeout_override) = read_timeout_override
            && matches!(
                resolved.transport,
                AgentSystemTransport::ExternalJsonStdio
                    | AgentSystemTransport::NativeOpenai
                    | AgentSystemTransport::NativeOpenaiCompat
            )
        {
            resolved.env.insert(
                AGENT_READ_TIMEOUT_SECS_ENV.to_string(),
                read_timeout_override.to_string(),
            );
        }
        let max_retries_override = self.agent_max_retries.trim().parse::<usize>().ok();
        if let Some(max_retries_override) = max_retries_override
            && matches!(
                resolved.transport,
                AgentSystemTransport::ExternalJsonStdio
                    | AgentSystemTransport::NativeOpenai
                    | AgentSystemTransport::NativeOpenaiCompat
            )
        {
            resolved.env.insert(
                AGENT_MAX_RETRIES_ENV.to_string(),
                max_retries_override.to_string(),
            );
        }
        let max_response_bytes_override = self
            .agent_max_response_bytes
            .trim()
            .parse::<usize>()
            .ok()
            .filter(|value| *value > 0);
        if let Some(max_response_bytes_override) = max_response_bytes_override
            && matches!(
                resolved.transport,
                AgentSystemTransport::ExternalJsonStdio
                    | AgentSystemTransport::NativeOpenai
                    | AgentSystemTransport::NativeOpenaiCompat
            )
        {
            resolved.env.insert(
                AGENT_MAX_RESPONSE_BYTES_ENV.to_string(),
                max_response_bytes_override.to_string(),
            );
        }
        resolved
    }

    fn selected_agent_runtime_base_url(&self, system: &AgentSystemSpec) -> Option<String> {
        if !matches!(
            system.transport,
            AgentSystemTransport::NativeOpenai | AgentSystemTransport::NativeOpenaiCompat
        ) {
            return None;
        }
        let override_base_url = self.agent_base_url_override.trim();
        if !override_base_url.is_empty() {
            return Some(override_base_url.to_string());
        }
        if let Some(catalog_base_url) = system
            .base_url
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            return Some(catalog_base_url.to_string());
        }
        Some(match system.transport {
            AgentSystemTransport::NativeOpenai => GUI_OPENAI_DEFAULT_BASE_URL.to_string(),
            AgentSystemTransport::NativeOpenaiCompat => {
                GUI_OPENAI_COMPAT_DEFAULT_BASE_URL.to_string()
            }
            _ => return None,
        })
    }

    fn selected_agent_model_discovery_source_key(
        &self,
        system: &AgentSystemSpec,
    ) -> Option<String> {
        let base_url = self.selected_agent_runtime_base_url(system)?;
        let key_state = if self.agent_openai_api_key.trim().is_empty() {
            "nokey"
        } else {
            "key"
        };
        Some(format!(
            "{}|{}|{}|{}",
            system.id,
            system.transport.as_str(),
            base_url,
            key_state
        ))
    }

    fn parse_agent_timeout_seconds(&self) -> Result<Option<u64>, String> {
        let raw = self.agent_timeout_secs.trim();
        if raw.is_empty() {
            return Ok(None);
        }
        let parsed = raw
            .parse::<u64>()
            .map_err(|e| format!("Invalid timeout_sec '{}': {}", raw, e))?;
        if parsed == 0 {
            return Ok(None);
        }
        Ok(Some(parsed))
    }

    fn parse_agent_connect_timeout_seconds(&self) -> Result<Option<u64>, String> {
        let raw = self.agent_connect_timeout_secs.trim();
        if raw.is_empty() {
            return Ok(None);
        }
        let parsed = raw
            .parse::<u64>()
            .map_err(|e| format!("Invalid connect_timeout_sec '{}': {}", raw, e))?;
        if parsed == 0 {
            return Ok(None);
        }
        Ok(Some(parsed))
    }

    fn parse_agent_read_timeout_seconds(&self) -> Result<Option<u64>, String> {
        let raw = self.agent_read_timeout_secs.trim();
        if raw.is_empty() {
            return Ok(None);
        }
        let parsed = raw
            .parse::<u64>()
            .map_err(|e| format!("Invalid read_timeout_sec '{}': {}", raw, e))?;
        if parsed == 0 {
            return Ok(None);
        }
        Ok(Some(parsed))
    }

    fn parse_agent_max_retries(&self) -> Result<Option<usize>, String> {
        let raw = self.agent_max_retries.trim();
        if raw.is_empty() {
            return Ok(None);
        }
        let parsed = raw
            .parse::<usize>()
            .map_err(|e| format!("Invalid max_retries '{}': {}", raw, e))?;
        Ok(Some(parsed))
    }

    fn parse_agent_max_response_bytes(&self) -> Result<Option<usize>, String> {
        let raw = self.agent_max_response_bytes.trim();
        if raw.is_empty() {
            return Ok(None);
        }
        let parsed = raw
            .parse::<usize>()
            .map_err(|e| format!("Invalid max_response_bytes '{}': {}", raw, e))?;
        if parsed == 0 {
            return Ok(None);
        }
        Ok(Some(parsed))
    }

    fn start_agent_model_discovery_task(&mut self, system: &AgentSystemSpec, force: bool) {
        if !matches!(
            system.transport,
            AgentSystemTransport::NativeOpenai | AgentSystemTransport::NativeOpenaiCompat
        ) {
            return;
        }
        let Some(base_url) = self.selected_agent_runtime_base_url(system) else {
            return;
        };
        let Some(source_key) = self.selected_agent_model_discovery_source_key(system) else {
            return;
        };
        if !force {
            if let Some(task) = &self.agent_model_discovery_task {
                if task.source_key == source_key {
                    return;
                }
            }
            if self.agent_model_discovery_source_key == source_key
                && !self.agent_discovered_models.is_empty()
            {
                return;
            }
        }
        self.agent_model_discovery_source_key = source_key.clone();
        self.agent_model_discovery_status = format!("Discovering models at {base_url} ...");
        self.agent_model_discovery_task = None;
        let api_key = self.agent_openai_api_key.trim().to_string();
        let api_key = if api_key.is_empty() {
            None
        } else {
            Some(api_key)
        };
        let (tx, rx) = mpsc::channel::<AgentModelDiscoveryTaskMessage>();
        self.agent_model_discovery_task = Some(AgentModelDiscoveryTask {
            started: Instant::now(),
            source_key: source_key.clone(),
            receiver: rx,
        });
        std::thread::spawn(move || {
            let result = discover_openai_models(&base_url, api_key.as_deref());
            let _ = tx.send(AgentModelDiscoveryTaskMessage::Done { source_key, result });
        });
    }

    fn selected_agent_system_availability(
        &self,
        system: &AgentSystemSpec,
    ) -> (bool, Option<String>) {
        let resolved = self.selected_agent_system_with_session_overrides(system);
        let availability = agent_system_availability(&resolved);
        if availability.available {
            return (true, availability.reason);
        }
        if matches!(system.transport, AgentSystemTransport::NativeOpenai)
            && !self.agent_openai_api_key.trim().is_empty()
        {
            return (
                true,
                Some("using GUI-supplied OpenAI API key for this session".to_string()),
            );
        }
        (false, availability.reason)
    }

    fn start_agent_assistant_request(&mut self) {
        if self.agent_task.is_some() {
            self.agent_status = "Agent request is already running".to_string();
            return;
        }
        self.refresh_agent_system_catalog();
        if !self.agent_catalog_error.is_empty() {
            self.agent_status = format!("Agent catalog error: {}", self.agent_catalog_error);
            return;
        }
        let system_id = self.agent_system_id.trim().to_string();
        if system_id.is_empty() {
            self.agent_status = "Select an agent system first".to_string();
            return;
        }
        let Some(selected_system) = self.selected_agent_system() else {
            self.agent_status = "Selected agent system is not available in catalog".to_string();
            return;
        };
        let (available, reason) = self.selected_agent_system_availability(&selected_system);
        if !available {
            self.agent_status = format!(
                "Selected agent system is unavailable: {}",
                reason.unwrap_or_else(|| "unknown reason".to_string())
            );
            return;
        }
        let prompt = self.agent_prompt.trim().to_string();
        if prompt.is_empty() {
            self.agent_status = "Agent prompt cannot be empty".to_string();
            return;
        }
        let timeout_seconds = match self.parse_agent_timeout_seconds() {
            Ok(value) => value,
            Err(err) => {
                self.agent_status = err;
                return;
            }
        };
        let connect_timeout_seconds = match self.parse_agent_connect_timeout_seconds() {
            Ok(value) => value,
            Err(err) => {
                self.agent_status = err;
                return;
            }
        };
        let read_timeout_seconds = match self.parse_agent_read_timeout_seconds() {
            Ok(value) => value,
            Err(err) => {
                self.agent_status = err;
                return;
            }
        };
        let max_retries = match self.parse_agent_max_retries() {
            Ok(value) => value,
            Err(err) => {
                self.agent_status = err;
                return;
            }
        };
        let max_response_bytes = match self.parse_agent_max_response_bytes() {
            Ok(value) => value,
            Err(err) => {
                self.agent_status = err;
                return;
            }
        };
        let mut model_override = normalize_agent_model_name(self.agent_model_override.trim());
        let discovered_pick =
            normalize_agent_model_name(&self.agent_discovered_model_pick).filter(|picked| {
                self.agent_discovered_models
                    .iter()
                    .any(|item| item == picked)
            });
        if model_override.is_none() && discovered_pick.is_some() {
            model_override = discovered_pick.clone();
        }
        if matches!(
            selected_system.transport,
            AgentSystemTransport::NativeOpenaiCompat
        ) && model_override.is_none()
        {
            if discovered_pick.is_some() {
                model_override = discovered_pick;
            } else {
                let catalog_model = normalize_agent_model_name(
                    selected_system.model.as_deref().unwrap_or_default(),
                );
                if let Some(catalog_model) = catalog_model {
                    if !self.agent_discovered_models.is_empty()
                        && !self
                            .agent_discovered_models
                            .iter()
                            .any(|value| value == &catalog_model)
                    {
                        self.agent_status = format!(
                            "Catalog model '{catalog_model}' is not available on current endpoint. Select a discovered model or set Model override."
                        );
                        return;
                    }
                } else {
                    self.agent_status =
                        "Model is unspecified. Discover models and select one, or set Model override."
                            .to_string();
                    return;
                }
            }
        }

        let state_summary = if self.agent_include_state_summary {
            Some(self.engine.read().unwrap().summarize_state())
        } else {
            None
        };
        let catalog_path = self.agent_catalog_path.trim().to_string();
        let openai_api_key = self.agent_openai_api_key.trim().to_string();
        let base_url_override = self.agent_base_url_override.trim().to_string();
        let model_override = model_override.unwrap_or_default();
        let timeout_override = timeout_seconds
            .map(|value| value.to_string())
            .unwrap_or_default();
        let connect_timeout_override = connect_timeout_seconds
            .map(|value| value.to_string())
            .unwrap_or_default();
        let read_timeout_override = read_timeout_seconds
            .map(|value| value.to_string())
            .unwrap_or_default();
        let max_retries_override = max_retries
            .map(|value| value.to_string())
            .unwrap_or_default();
        let max_response_bytes_override = max_response_bytes
            .map(|value| value.to_string())
            .unwrap_or_default();
        let job_id = self.alloc_background_job_id();
        let (tx, rx) = mpsc::channel::<AgentAskTaskMessage>();
        self.agent_status = if let Some(timeout) = timeout_seconds {
            format!(
                "Asking agent '{}' in background (timeout={}s, retries={})",
                system_id,
                timeout,
                max_retries.unwrap_or(2)
            )
        } else {
            format!("Asking agent '{}' in background", system_id)
        };
        self.push_job_event(
            BackgroundJobKind::AgentAssist,
            BackgroundJobEventPhase::Started,
            Some(job_id),
            format!("Agent request started for system '{}'", system_id),
        );
        self.agent_task = Some(AgentAskTask {
            job_id,
            started: Instant::now(),
            receiver: rx,
        });
        std::thread::spawn(move || {
            let mut env_overrides = HashMap::new();
            if !openai_api_key.is_empty() {
                env_overrides.insert(OPENAI_API_KEY_ENV.to_string(), openai_api_key);
            }
            if !base_url_override.is_empty() {
                env_overrides.insert(AGENT_BASE_URL_ENV.to_string(), base_url_override);
            }
            if !model_override.is_empty() {
                env_overrides.insert(AGENT_MODEL_ENV.to_string(), model_override);
            }
            if !timeout_override.is_empty() {
                env_overrides.insert(AGENT_TIMEOUT_SECS_ENV.to_string(), timeout_override);
            }
            if !connect_timeout_override.is_empty() {
                env_overrides.insert(
                    AGENT_CONNECT_TIMEOUT_SECS_ENV.to_string(),
                    connect_timeout_override,
                );
            }
            if !read_timeout_override.is_empty() {
                env_overrides.insert(
                    AGENT_READ_TIMEOUT_SECS_ENV.to_string(),
                    read_timeout_override,
                );
            }
            if !max_retries_override.is_empty() {
                env_overrides.insert(AGENT_MAX_RETRIES_ENV.to_string(), max_retries_override);
            }
            if !max_response_bytes_override.is_empty() {
                env_overrides.insert(
                    AGENT_MAX_RESPONSE_BYTES_ENV.to_string(),
                    max_response_bytes_override,
                );
            }
            let result = invoke_agent_support_with_env_overrides(
                Some(catalog_path.as_str()),
                &system_id,
                &prompt,
                state_summary.as_ref(),
                if env_overrides.is_empty() {
                    None
                } else {
                    Some(&env_overrides)
                },
            );
            let _ = tx.send(AgentAskTaskMessage::Done { job_id, result });
        });
    }

    fn execute_agent_suggested_command(
        &mut self,
        index_1based: usize,
        command_text: &str,
        trigger: &str,
    ) {
        let trimmed = command_text.trim();
        if trimmed.is_empty() {
            self.agent_status = format!("Suggestion #{index_1based} is empty");
            return;
        }
        let command = match parse_shell_line(trimmed) {
            Ok(command) => command,
            Err(err) => {
                self.agent_status = format!("Suggestion #{index_1based} parse error: {err}");
                self.agent_execution_log.push(AgentCommandExecutionRecord {
                    index_1based,
                    command: trimmed.to_string(),
                    trigger: trigger.to_string(),
                    ok: false,
                    state_changed: false,
                    summary: format!("parse error: {err}"),
                    executed_at_unix_ms: Self::now_unix_ms(),
                });
                return;
            }
        };
        if matches!(command, ShellCommand::AgentsAsk { .. }) {
            self.agent_status = format!(
                "Suggestion #{index_1based} rejected: agent-to-agent 'agents ask' is blocked"
            );
            self.agent_execution_log.push(AgentCommandExecutionRecord {
                index_1based,
                command: trimmed.to_string(),
                trigger: trigger.to_string(),
                ok: false,
                state_changed: false,
                summary: "agent-to-agent agents ask blocked".to_string(),
                executed_at_unix_ms: Self::now_unix_ms(),
            });
            return;
        }
        if let Some(summary) = self.try_apply_shell_ui_intent(&command) {
            self.agent_status = format!("Suggestion #{index_1based}: {summary}");
            self.agent_execution_log.push(AgentCommandExecutionRecord {
                index_1based,
                command: trimmed.to_string(),
                trigger: trigger.to_string(),
                ok: true,
                state_changed: false,
                summary,
                executed_at_unix_ms: Self::now_unix_ms(),
            });
            if self.agent_execution_log.len() > 100 {
                let drain = self.agent_execution_log.len() - 100;
                self.agent_execution_log.drain(0..drain);
            }
            return;
        }
        let options = ShellExecutionOptions {
            allow_screenshots: false,
            allow_agent_commands: false,
        };
        let run = {
            let mut guard = self.engine.write().unwrap();
            execute_shell_command_with_options(&mut guard, &command, &options)
        };
        match run {
            Ok(run) => {
                if run.state_changed {
                    self.lineage_cache_valid = false;
                }
                let summary = if run.state_changed {
                    "executed (state changed)".to_string()
                } else {
                    "executed".to_string()
                };
                self.agent_status = format!("Suggestion #{index_1based}: {summary}");
                self.agent_execution_log.push(AgentCommandExecutionRecord {
                    index_1based,
                    command: trimmed.to_string(),
                    trigger: trigger.to_string(),
                    ok: true,
                    state_changed: run.state_changed,
                    summary,
                    executed_at_unix_ms: Self::now_unix_ms(),
                });
            }
            Err(err) => {
                self.agent_status = format!("Suggestion #{index_1based} failed: {err}");
                self.agent_execution_log.push(AgentCommandExecutionRecord {
                    index_1based,
                    command: trimmed.to_string(),
                    trigger: trigger.to_string(),
                    ok: false,
                    state_changed: false,
                    summary: err,
                    executed_at_unix_ms: Self::now_unix_ms(),
                });
            }
        }
        if self.agent_execution_log.len() > 100 {
            let drain = self.agent_execution_log.len() - 100;
            self.agent_execution_log.drain(0..drain);
        }
    }

    fn try_apply_shell_ui_intent(&mut self, command: &ShellCommand) -> Option<String> {
        let ShellCommand::UiIntent {
            action,
            target,
            genome_id,
            helper_mode,
            catalog_path,
            cache_dir,
            filter,
            species,
            latest,
        } = command
        else {
            return None;
        };
        let mut selected_genome_id = genome_id
            .as_deref()
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .map(str::to_string);
        if matches!(target, UiIntentTarget::PreparedReferences) {
            self.apply_prepared_reference_intent_scope(*helper_mode, catalog_path, cache_dir);
            if selected_genome_id.is_none() {
                match self.resolve_prepared_reference_intent_selection(
                    *helper_mode,
                    catalog_path.clone(),
                    cache_dir.clone(),
                    filter.clone(),
                    species.clone(),
                    *latest,
                ) {
                    Ok(Some(resolved)) => {
                        selected_genome_id = Some(resolved);
                    }
                    Ok(None) => {}
                    Err(err) => {
                        self.app_status = format!(
                            "Could not resolve prepared-reference selection for ui intent: {err}"
                        );
                    }
                }
            }
        }
        if let Some(genome_id) = selected_genome_id
            .as_deref()
            .map(str::trim)
            .filter(|v| !v.is_empty())
        {
            self.genome_id = genome_id.to_string();
            self.invalidate_genome_genes();
        }
        match target {
            UiIntentTarget::PreparedReferences => self.open_reference_genome_inspector_dialog(),
            UiIntentTarget::PrepareReferenceGenome => self.open_reference_genome_prepare_dialog(),
            UiIntentTarget::RetrieveGenomeSequence => self.open_reference_genome_retrieve_dialog(),
            UiIntentTarget::BlastGenomeSequence => self.open_reference_genome_blast_dialog(),
            UiIntentTarget::ImportGenomeTrack => self.open_genome_bed_track_dialog(),
            UiIntentTarget::AgentAssistant => self.open_agent_assistant_dialog(),
            UiIntentTarget::PrepareHelperGenome => self.open_helper_genome_prepare_dialog(),
            UiIntentTarget::RetrieveHelperSequence => self.open_helper_genome_retrieve_dialog(),
            UiIntentTarget::BlastHelperSequence => self.open_helper_genome_blast_dialog(),
        }
        let mut summary = format!("ui intent {} '{}'", action.as_str(), target.as_str());
        if let Some(genome_id) = selected_genome_id {
            summary.push_str(&format!(" (selected_genome_id={genome_id})"));
        }
        Some(summary)
    }

    fn apply_prepared_reference_intent_scope(
        &mut self,
        helper_mode: bool,
        catalog_path: &Option<String>,
        cache_dir: &Option<String>,
    ) {
        let normalized_catalog = catalog_path
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(str::to_string);
        let normalized_cache = cache_dir
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(str::to_string);
        let next_catalog = if helper_mode {
            normalized_catalog.unwrap_or_else(|| DEFAULT_HELPER_GENOME_CATALOG_PATH.to_string())
        } else {
            normalized_catalog.unwrap_or_else(|| self.genome_catalog_path.clone())
        };
        let next_cache = if helper_mode {
            normalized_cache.unwrap_or_else(|| "data/helper_genomes".to_string())
        } else {
            normalized_cache.unwrap_or_else(|| self.genome_cache_dir.clone())
        };
        let catalog_changed = self.genome_catalog_path != next_catalog;
        let cache_changed = self.genome_cache_dir != next_cache;
        if catalog_changed {
            self.genome_catalog_path = next_catalog;
        }
        if cache_changed {
            self.genome_cache_dir = next_cache;
        }
        if catalog_changed || cache_changed {
            self.invalidate_genome_genes();
        }
    }

    fn resolve_prepared_reference_intent_selection(
        &self,
        helper_mode: bool,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
        filter: Option<String>,
        species: Option<String>,
        latest: bool,
    ) -> Result<Option<String>, String> {
        let mut engine = self.engine.write().unwrap();
        let run = execute_shell_command_with_options(
            &mut engine,
            &ShellCommand::UiPreparedGenomes {
                helper_mode,
                catalog_path,
                cache_dir,
                filter,
                species,
                latest,
            },
            &ShellExecutionOptions::default(),
        )?;
        Ok(run
            .output
            .get("selected_genome_id")
            .and_then(|value| value.as_str())
            .map(str::to_string))
    }

    fn execute_agent_auto_suggestions(&mut self, response: &AgentResponse) {
        for (idx, suggestion) in response.suggested_commands.iter().enumerate() {
            if suggestion.execution == AgentExecutionIntent::Auto {
                self.execute_agent_suggested_command(idx + 1, &suggestion.command, "auto");
            }
        }
    }

    fn poll_agent_assistant_task(&mut self, ctx: &egui::Context) {
        if self.agent_task.is_none() {
            return;
        }
        ctx.request_repaint_after(Duration::from_millis(100));
        let mut done: Option<(u64, Result<AgentInvocationOutcome, String>)> = None;
        if let Some(task) = &self.agent_task {
            match task.receiver.try_recv() {
                Ok(AgentAskTaskMessage::Done { job_id, result }) => {
                    done = Some((job_id, result));
                }
                Err(mpsc::TryRecvError::Empty) => {}
                Err(mpsc::TryRecvError::Disconnected) => {
                    done = Some((task.job_id, Err("Agent worker disconnected".to_string())));
                }
            }
        }
        if let Some((job_id, outcome)) = done {
            let elapsed = self
                .agent_task
                .as_ref()
                .map(|task| task.started.elapsed().as_secs_f64())
                .unwrap_or(0.0);
            self.agent_task = None;
            match outcome {
                Ok(invocation) => {
                    let suggestion_count = invocation.response.suggested_commands.len();
                    self.agent_status = format!(
                        "Agent response received in {:.1}s (suggestions={})",
                        elapsed, suggestion_count
                    );
                    self.push_job_event(
                        BackgroundJobKind::AgentAssist,
                        BackgroundJobEventPhase::Completed,
                        Some(job_id),
                        format!(
                            "Agent '{}' completed in {:.1}s (suggestions={})",
                            invocation.system_id, elapsed, suggestion_count
                        ),
                    );
                    let response = invocation.response.clone();
                    self.agent_last_invocation = Some(invocation);
                    if self.agent_allow_auto_exec {
                        self.execute_agent_auto_suggestions(&response);
                    }
                }
                Err(err) => {
                    self.agent_status =
                        format!("Agent request failed after {:.1}s: {}", elapsed, err);
                    self.push_job_event(
                        BackgroundJobKind::AgentAssist,
                        BackgroundJobEventPhase::Failed,
                        Some(job_id),
                        format!("Agent request failed in {:.1}s: {}", elapsed, err),
                    );
                }
            }
        }
    }

    fn poll_agent_model_discovery_task(&mut self, ctx: &egui::Context) {
        if self.agent_model_discovery_task.is_none() {
            return;
        }
        ctx.request_repaint_after(Duration::from_millis(100));
        let mut done: Option<(String, Result<Vec<String>, String>)> = None;
        if let Some(task) = &self.agent_model_discovery_task {
            match task.receiver.try_recv() {
                Ok(AgentModelDiscoveryTaskMessage::Done { source_key, result }) => {
                    done = Some((source_key, result));
                }
                Err(mpsc::TryRecvError::Empty) => {}
                Err(mpsc::TryRecvError::Disconnected) => {
                    done = Some((
                        task.source_key.clone(),
                        Err("Model discovery worker disconnected".to_string()),
                    ));
                }
            }
        }
        if let Some((source_key, result)) = done {
            let elapsed = self
                .agent_model_discovery_task
                .as_ref()
                .map(|task| task.started.elapsed().as_secs_f64())
                .unwrap_or(0.0);
            self.agent_model_discovery_task = None;
            if source_key != self.agent_model_discovery_source_key {
                return;
            }
            match result {
                Ok(models) => {
                    self.agent_discovered_models = models;
                    if self.agent_discovered_models.is_empty() {
                        self.agent_model_discovery_status =
                            format!("Model discovery returned no models ({:.1}s)", elapsed);
                        self.agent_discovered_model_pick.clear();
                    } else {
                        if !self
                            .agent_discovered_models
                            .iter()
                            .any(|item| item == &self.agent_discovered_model_pick)
                        {
                            self.agent_discovered_model_pick.clear();
                        }
                        self.agent_model_discovery_status = format!(
                            "Discovered {} model(s) in {:.1}s",
                            self.agent_discovered_models.len(),
                            elapsed
                        );
                    }
                }
                Err(err) => {
                    self.agent_discovered_models.clear();
                    self.agent_discovered_model_pick.clear();
                    self.agent_model_discovery_status =
                        format!("Model discovery failed after {:.1}s: {}", elapsed, err);
                }
            }
        }
    }

    fn invalidate_genome_genes(&mut self) {
        self.genome_genes.clear();
        self.genome_genes_error.clear();
        self.genome_genes_loaded_key = None;
        self.genome_selected_gene = None;
        self.genome_gene_filter_page = 0;
        self.genome_biotype_filter.clear();
    }

    fn genome_genes_loaded_key(&self) -> String {
        format!(
            "{}|{}|{}",
            self.genome_catalog_path_resolved(),
            self.genome_cache_dir_opt().unwrap_or_default(),
            self.genome_id.trim()
        )
    }

    fn ensure_genome_genes_loaded(&mut self) {
        let genome_id = self.genome_id.trim().to_string();
        if genome_id.is_empty() {
            self.invalidate_genome_genes();
            return;
        }
        let key = self.genome_genes_loaded_key();
        if self
            .genome_genes_loaded_key
            .as_ref()
            .map(|v| v == &key)
            .unwrap_or(false)
        {
            return;
        }
        self.genome_genes.clear();
        self.genome_selected_gene = None;
        self.genome_genes_error.clear();
        let catalog_path = self.genome_catalog_path_resolved();
        let cache_dir = self.genome_cache_dir_opt();
        let loaded = GenomeCatalog::from_json_file(&catalog_path)
            .and_then(|catalog| catalog.list_gene_regions(&genome_id, cache_dir.as_deref()));
        match loaded {
            Ok(genes) => {
                self.genome_genes = genes;
                self.sync_genome_biotype_filter();
            }
            Err(e) => {
                self.genome_genes_error = e;
                self.genome_biotype_filter.clear();
            }
        }
        self.genome_genes_loaded_key = Some(key);
    }

    fn sync_genome_biotype_filter(&mut self) {
        let mut discovered: Vec<String> = self
            .genome_genes
            .iter()
            .filter_map(|gene| gene.biotype.as_ref())
            .map(|v| v.trim())
            .filter(|v| !v.is_empty())
            .map(str::to_string)
            .collect();
        discovered.sort();
        discovered.dedup();

        let mut next: HashMap<String, bool> = HashMap::with_capacity(discovered.len());
        for biotype in discovered {
            let enabled = self
                .genome_biotype_filter
                .get(&biotype)
                .copied()
                .unwrap_or(true);
            next.insert(biotype, enabled);
        }
        self.genome_biotype_filter = next;
    }

    fn unprepared_genomes_for_prepare_dialog(&self) -> Result<Vec<String>, String> {
        let catalog_path = self.genome_catalog_path_resolved();
        let catalog = GenomeCatalog::from_json_file(&catalog_path)?;
        let cache_dir = self.genome_cache_dir_opt();
        let mut names: Vec<String> = vec![];
        for name in catalog.list_genomes() {
            match catalog.is_prepared(&name, cache_dir.as_deref()) {
                Ok(false) => names.push(name),
                Ok(true) => {}
                Err(e) => {
                    return Err(format!(
                        "Could not check preparation status for '{name}': {e}"
                    ));
                }
            }
        }
        Ok(names)
    }

    fn prepared_genomes_for_retrieve_dialog(&self) -> Result<Vec<String>, String> {
        let catalog_path = self.genome_catalog_path_resolved();
        let catalog = GenomeCatalog::from_json_file(&catalog_path)?;
        let cache_dir = self.genome_cache_dir_opt();
        let mut names: Vec<String> = vec![];
        for name in catalog.list_genomes() {
            match catalog.is_prepared(&name, cache_dir.as_deref()) {
                Ok(true) => names.push(name),
                Ok(false) => {}
                Err(e) => {
                    return Err(format!(
                        "Could not check preparation status for '{name}': {e}"
                    ));
                }
            }
        }
        Ok(names)
    }

    fn project_sequence_ids_for_blast(&self) -> Vec<String> {
        let mut seq_ids: Vec<String> = self
            .engine
            .read()
            .unwrap()
            .state()
            .sequences
            .keys()
            .cloned()
            .collect();
        seq_ids.sort_unstable();
        seq_ids
    }

    fn blast_pool_options(&self) -> Vec<BlastPoolOption> {
        let engine = self.engine.read().unwrap();
        let state = engine.state();
        let mut options: Vec<BlastPoolOption> = state
            .container_state
            .containers
            .values()
            .filter_map(|container| {
                if container.members.len() <= 1 {
                    return None;
                }
                let mut members: Vec<String> = container
                    .members
                    .iter()
                    .filter(|seq_id| state.sequences.contains_key(*seq_id))
                    .cloned()
                    .collect();
                if members.len() <= 1 {
                    return None;
                }
                members.sort_unstable();
                members.dedup();
                let name = container
                    .name
                    .as_ref()
                    .map(|v| v.trim())
                    .filter(|v| !v.is_empty())
                    .unwrap_or("-");
                Some(BlastPoolOption {
                    container_id: container.container_id.clone(),
                    label: format!(
                        "{} (kind={:?}, name={}, n={})",
                        container.container_id,
                        container.kind,
                        name,
                        members.len()
                    ),
                    members,
                })
            })
            .collect();
        options.sort_by(|a, b| a.container_id.cmp(&b.container_id));
        options
    }

    fn selected_genome_source_plan(&self) -> Result<Option<GenomeSourcePlan>, String> {
        let genome_id = self.genome_id.trim();
        if genome_id.is_empty() {
            return Ok(None);
        }
        let catalog_path = self.genome_catalog_path_resolved();
        let catalog = GenomeCatalog::from_json_file(&catalog_path)?;
        let cache_dir = self.genome_cache_dir_opt();
        let plan = catalog.source_plan(genome_id, cache_dir.as_deref())?;
        Ok(Some(plan))
    }

    fn format_genome_source_plan_summary(plan: &GenomeSourcePlan) -> String {
        let mut parts = vec![format!(
            "sources: sequence={} | annotation={}",
            plan.sequence_source_type, plan.annotation_source_type
        )];
        if let Some(length_bp) = plan.nucleotide_length_bp {
            parts.push(format!("length={length_bp} bp"));
        }
        if let Some(mass_da) = plan.molecular_mass_da {
            let source = plan
                .molecular_mass_source
                .as_deref()
                .unwrap_or("unknown_source");
            parts.push(format!("mass={mass_da:.3e} Da ({source})"));
        }
        parts.join(" | ")
    }

    fn choose_genome_from_catalog(ui: &mut Ui, genome_id: &mut String, names: &[String]) -> bool {
        let mut changed = false;
        let selected_text = if names.iter().any(|name| name == genome_id) {
            genome_id.clone()
        } else if names.is_empty() {
            "(no genomes available)".to_string()
        } else {
            "(choose genome)".to_string()
        };
        egui::ComboBox::from_label("genome")
            .selected_text(selected_text)
            .show_ui(ui, |ui| {
                for name in names {
                    if ui.selectable_value(genome_id, name.clone(), name).changed() {
                        changed = true;
                    }
                }
            });
        changed
    }

    fn format_op_result_status(
        prefix: &str,
        created: &[String],
        warnings: &[String],
        messages: &[String],
    ) -> String {
        let created_text = if created.is_empty() {
            "-".to_string()
        } else {
            created.join(", ")
        };
        let warnings_text = if warnings.is_empty() {
            "-".to_string()
        } else {
            warnings.join(" | ")
        };
        let messages_text = if messages.is_empty() {
            "-".to_string()
        } else {
            messages.join(" | ")
        };
        format!(
            "{prefix}\ncreated: {created_text}\nwarnings: {warnings_text}\nmessages: {messages_text}"
        )
    }

    fn collect_prepared_genome_inspections(
        &self,
    ) -> Result<(Vec<PreparedGenomeInspection>, Vec<String>), String> {
        let catalog_path = self.genome_catalog_path_resolved();
        let catalog = GenomeCatalog::from_json_file(&catalog_path)?;
        let cache_dir = self.genome_cache_dir_opt();
        let mut inspections: Vec<PreparedGenomeInspection> = vec![];
        let mut errors: Vec<String> = vec![];
        for genome_id in catalog.list_genomes() {
            match catalog.inspect_prepared_genome(&genome_id, cache_dir.as_deref()) {
                Ok(Some(inspection)) => inspections.push(inspection),
                Ok(None) => {}
                Err(e) => errors.push(format!("{genome_id}: {e}")),
            }
        }
        inspections.sort_by(|a, b| a.genome_id.cmp(&b.genome_id));
        Ok((inspections, errors))
    }

    fn prepared_genome_chromosome_records(
        &self,
        genome_id: &str,
    ) -> Result<Vec<GenomeChromosomeRecord>, String> {
        let catalog_path = self.genome_catalog_path_resolved();
        let catalog = GenomeCatalog::from_json_file(&catalog_path)?;
        let cache_dir = self.genome_cache_dir_opt();
        catalog.list_chromosome_lengths(genome_id, cache_dir.as_deref())
    }

    fn render_chromosome_length_lines(ui: &mut Ui, chromosomes: &[GenomeChromosomeRecord]) {
        if chromosomes.is_empty() {
            ui.small("No chromosomes/contigs found in FASTA index.");
            return;
        }
        let longest_bp = chromosomes
            .iter()
            .map(|record| record.length_bp)
            .max()
            .unwrap_or(1)
            .max(1);
        let total_bp: u128 = chromosomes.iter().fold(0u128, |acc, record| {
            acc.saturating_add(record.length_bp as u128)
        });
        ui.small(format!(
            "contigs: {} | longest: {} ({} bp) | total span: {} bp",
            chromosomes.len(),
            chromosomes[0].chromosome,
            chromosomes[0].length_bp,
            total_bp
        ));
        egui::ScrollArea::vertical()
            .max_height(280.0)
            .show(ui, |ui| {
                for record in chromosomes {
                    ui.horizontal(|ui| {
                        ui.add_sized(
                            [220.0, 0.0],
                            egui::Label::new(
                                egui::RichText::new(record.chromosome.clone()).monospace(),
                            ),
                        );
                        ui.add_sized(
                            [120.0, 0.0],
                            egui::Label::new(format!("{} bp", record.length_bp)),
                        );
                        let width = ui.available_width().max(80.0);
                        let (rect, response) =
                            ui.allocate_exact_size(egui::vec2(width, 14.0), egui::Sense::hover());
                        let ratio = (record.length_bp as f32 / longest_bp as f32).clamp(0.0, 1.0);
                        let line_width = (rect.width() * ratio.max(0.005)).max(1.0);
                        let y = rect.center().y;
                        ui.painter().line_segment(
                            [
                                Pos2::new(rect.left(), y),
                                Pos2::new(rect.left() + line_width, y),
                            ],
                            egui::Stroke::new(2.0, egui::Color32::from_rgb(80, 130, 175)),
                        );
                        response.on_hover_text(format!(
                            "{} bp ({:.2}% of longest)",
                            record.length_bp,
                            ratio * 100.0
                        ));
                    });
                }
            });
    }

    fn format_bytes_compact(bytes: u64) -> String {
        const UNITS: [&str; 5] = ["B", "KB", "MB", "GB", "TB"];
        let mut value = bytes as f64;
        let mut unit = 0usize;
        while value >= 1024.0 && unit + 1 < UNITS.len() {
            value /= 1024.0;
            unit += 1;
        }
        if unit == 0 {
            format!("{bytes} {}", UNITS[unit])
        } else {
            format!("{value:.2} {}", UNITS[unit])
        }
    }

    fn format_short_sha1(value: &Option<String>) -> String {
        value
            .as_ref()
            .map(|v| v.trim())
            .filter(|v| !v.is_empty())
            .map(|v| {
                if v.len() > 12 {
                    v[..12].to_string()
                } else {
                    v.to_string()
                }
            })
            .unwrap_or_else(|| "-".to_string())
    }

    fn parse_prepare_timeout_seconds(&self) -> Result<Option<u64>, String> {
        let raw = self.genome_prepare_timeout_secs.trim();
        if raw.is_empty() {
            return Ok(None);
        }
        let parsed = raw
            .parse::<u64>()
            .map_err(|e| format!("Invalid timeout seconds '{}': {}", raw, e))?;
        if parsed == 0 {
            return Ok(None);
        }
        Ok(Some(parsed))
    }

    fn start_prepare_reference_genome(&mut self) {
        if self.genome_prepare_task.is_some() {
            self.genome_prepare_status = "Genome preparation is already running".to_string();
            return;
        }
        let genome_id = self.genome_id.trim().to_string();
        if genome_id.is_empty() {
            self.genome_prepare_status = "Select a genome first".to_string();
            return;
        }
        let timeout_seconds = match self.parse_prepare_timeout_seconds() {
            Ok(value) => value,
            Err(e) => {
                self.genome_prepare_status = e;
                return;
            }
        };
        let catalog_path = self.genome_catalog_path_opt();
        let cache_dir = self.genome_cache_dir_opt();
        let cancel_requested = Arc::new(AtomicBool::new(false));
        let job_id = self.alloc_background_job_id();
        let (tx, rx) = mpsc::channel::<GenomePrepareTaskMessage>();
        self.genome_prepare_progress = None;
        self.genome_prepare_status = if let Some(timeout) = timeout_seconds {
            format!(
                "Preparing genome '{genome_id}' in background (timeout: {} s). You can keep using the UI.",
                timeout
            )
        } else {
            format!("Preparing genome '{genome_id}' in background. You can keep using the UI.")
        };
        self.push_job_event(
            BackgroundJobKind::PrepareGenome,
            BackgroundJobEventPhase::Started,
            Some(job_id),
            format!("Prepare genome started: {}", genome_id),
        );
        self.genome_prepare_task = Some(GenomePrepareTask {
            job_id,
            started: Instant::now(),
            cancel_requested: cancel_requested.clone(),
            timeout_seconds,
            receiver: rx,
        });
        std::thread::spawn(move || {
            let tx_progress = tx.clone();
            let cancel_flag = cancel_requested.clone();
            let started = Instant::now();
            let mut last_phase = String::new();
            let mut last_percent_tenths: Option<i64> = None;
            let mut last_bytes_bucket: u64 = 0;
            let mut progress_forwarder = move |p: PrepareGenomeProgress| -> bool {
                if cancel_flag.load(Ordering::Relaxed) {
                    return false;
                }
                if let Some(limit) = timeout_seconds {
                    if started.elapsed() >= Duration::from_secs(limit) {
                        return false;
                    }
                }
                let phase_changed = p.phase != last_phase;
                let percent_tenths = p.percent.map(|v| (v * 10.0).floor() as i64);
                let percent_progressed = percent_tenths
                    .map(|value| last_percent_tenths.map(|prev| value > prev).unwrap_or(true))
                    .unwrap_or(false);
                let bytes_bucket = p.bytes_done / (8 * 1024 * 1024);
                let bytes_progressed = bytes_bucket > last_bytes_bucket;
                let done = p.phase == "ready"
                    || p.percent
                        .map(|v| (v - 100.0).abs() < f64::EPSILON)
                        .unwrap_or(false);
                if phase_changed || percent_progressed || bytes_progressed || done {
                    last_phase = p.phase.clone();
                    if let Some(v) = percent_tenths {
                        last_percent_tenths = Some(v);
                    }
                    last_bytes_bucket = bytes_bucket;
                    let _ = tx_progress.send(GenomePrepareTaskMessage::Progress {
                        job_id,
                        progress: p,
                    });
                }
                true
            };
            let outcome = GentleEngine::prepare_reference_genome_once(
                &genome_id,
                catalog_path.as_deref(),
                cache_dir.as_deref(),
                timeout_seconds,
                &mut progress_forwarder,
            )
            .map(|report| OpResult {
                op_id: "background-prepare-genome".to_string(),
                created_seq_ids: vec![],
                changed_seq_ids: vec![],
                warnings: report.warnings.clone(),
                messages: vec![GentleEngine::format_prepare_genome_message(
                    &genome_id,
                    cache_dir.as_deref(),
                    &report,
                )],
            });
            let _ = tx.send(GenomePrepareTaskMessage::Done {
                job_id,
                result: outcome,
            });
        });
    }

    fn poll_prepare_reference_genome_task(&mut self, ctx: &egui::Context) {
        if self.genome_prepare_task.is_none() {
            return;
        }
        ctx.request_repaint_after(Duration::from_millis(100));
        let mut done: Option<(u64, Result<OpResult, EngineError>)> = None;
        let mut stale_job_ids: Vec<u64> = vec![];
        if let Some(task) = &self.genome_prepare_task {
            let active_job_id = task.job_id;
            const MAX_MESSAGES_PER_TICK: usize = 128;
            for _ in 0..MAX_MESSAGES_PER_TICK {
                match task.receiver.try_recv() {
                    Ok(GenomePrepareTaskMessage::Progress { job_id, progress }) => {
                        if job_id != active_job_id {
                            if !stale_job_ids.contains(&job_id) {
                                stale_job_ids.push(job_id);
                            }
                            continue;
                        }
                        self.genome_prepare_progress = Some(progress.clone());
                        self.genome_prepare_status = format!(
                            "Preparing genome '{}': {} ({})",
                            progress.genome_id, progress.phase, progress.item
                        );
                    }
                    Ok(GenomePrepareTaskMessage::Done { job_id, result }) => {
                        if job_id != active_job_id {
                            if !stale_job_ids.contains(&job_id) {
                                stale_job_ids.push(job_id);
                            }
                            continue;
                        }
                        done = Some((job_id, result));
                        break;
                    }
                    Err(mpsc::TryRecvError::Empty) => break,
                    Err(mpsc::TryRecvError::Disconnected) => {
                        done = Some((
                            active_job_id,
                            Err(EngineError {
                                code: ErrorCode::Internal,
                                message: "Genome preparation worker disconnected".to_string(),
                            }),
                        ));
                        break;
                    }
                }
            }
        }
        for stale_job_id in stale_job_ids {
            self.push_job_event(
                BackgroundJobKind::PrepareGenome,
                BackgroundJobEventPhase::IgnoredStale,
                Some(stale_job_id),
                "Ignored stale prepare-genome worker message",
            );
        }
        if let Some((job_id, outcome)) = done {
            let elapsed = self
                .genome_prepare_task
                .as_ref()
                .map(|task| task.started.elapsed().as_secs_f64())
                .unwrap_or(0.0);
            self.genome_prepare_task = None;
            match outcome {
                Ok(result) => {
                    self.genome_prepare_status = format!(
                        "{}\nelapsed: {:.1}s",
                        Self::format_op_result_status(
                            "Prepare genome: ok",
                            &result.created_seq_ids,
                            &result.warnings,
                            &result.messages,
                        ),
                        elapsed
                    );
                    self.invalidate_genome_genes();
                    self.push_job_event(
                        BackgroundJobKind::PrepareGenome,
                        BackgroundJobEventPhase::Completed,
                        Some(job_id),
                        format!("Prepare genome completed in {:.1}s", elapsed),
                    );
                }
                Err(e) => {
                    let lower = e.message.to_ascii_lowercase();
                    if lower.contains("timed out") {
                        self.genome_prepare_status = format!(
                            "Prepare genome timed out after {:.1}s: {}",
                            elapsed, e.message
                        );
                    } else if lower.contains("cancelled") || lower.contains("canceled") {
                        self.genome_prepare_status = format!(
                            "Prepare genome cancelled after {:.1}s: {}",
                            elapsed, e.message
                        );
                    } else {
                        self.genome_prepare_status =
                            format!("Prepare genome failed after {:.1}s: {}", elapsed, e.message);
                    }
                    self.push_job_event(
                        BackgroundJobKind::PrepareGenome,
                        BackgroundJobEventPhase::Failed,
                        Some(job_id),
                        format!("Prepare genome ended in {:.1}s: {}", elapsed, e.message),
                    );
                }
            }
        }
    }

    fn poll_genome_track_import_task(&mut self, ctx: &egui::Context) {
        if self.genome_track_import_task.is_none() {
            return;
        }
        ctx.request_repaint_after(Duration::from_millis(100));
        let mut done: Option<(u64, Result<OpResult, EngineError>)> = None;
        let mut stale_job_ids: Vec<u64> = vec![];
        if let Some(task) = &self.genome_track_import_task {
            let active_job_id = task.job_id;
            const MAX_MESSAGES_PER_TICK: usize = 256;
            for _ in 0..MAX_MESSAGES_PER_TICK {
                match task.receiver.try_recv() {
                    Ok(GenomeTrackImportTaskMessage::Progress { job_id, progress }) => {
                        if job_id != active_job_id {
                            if !stale_job_ids.contains(&job_id) {
                                stale_job_ids.push(job_id);
                            }
                            continue;
                        }
                        self.genome_track_import_progress = Some(progress.clone());
                        let canceling = task.cancel_requested.load(Ordering::Relaxed);
                        self.genome_track_status = format!(
                            "Importing {} track for '{}': parsed={}, imported={}, skipped={}{}",
                            progress.source,
                            progress.seq_id,
                            progress.parsed_records,
                            progress.imported_features,
                            progress.skipped_records,
                            if canceling {
                                " (cancellation requested)"
                            } else {
                                ""
                            }
                        );
                    }
                    Ok(GenomeTrackImportTaskMessage::Done { job_id, result }) => {
                        if job_id != active_job_id {
                            if !stale_job_ids.contains(&job_id) {
                                stale_job_ids.push(job_id);
                            }
                            continue;
                        }
                        done = Some((job_id, result));
                        break;
                    }
                    Err(mpsc::TryRecvError::Empty) => break,
                    Err(mpsc::TryRecvError::Disconnected) => {
                        done = Some((
                            active_job_id,
                            Err(EngineError {
                                code: ErrorCode::Internal,
                                message: "Genome track import worker disconnected".to_string(),
                            }),
                        ));
                        break;
                    }
                }
            }
        }
        for stale_job_id in stale_job_ids {
            self.push_job_event(
                BackgroundJobKind::TrackImport,
                BackgroundJobEventPhase::IgnoredStale,
                Some(stale_job_id),
                "Ignored stale track-import worker message",
            );
        }

        if let Some((job_id, outcome)) = done {
            let elapsed = self
                .genome_track_import_task
                .as_ref()
                .map(|task| task.started.elapsed().as_secs_f64())
                .unwrap_or(0.0);
            let cancellation_requested = self
                .genome_track_import_task
                .as_ref()
                .map(|task| task.cancel_requested.load(Ordering::Relaxed))
                .unwrap_or(false);
            self.genome_track_import_task = None;
            self.genome_track_import_progress = None;
            match outcome {
                Ok(result) => {
                    let prefix = if cancellation_requested {
                        "Import track finished after cancellation request"
                    } else {
                        "Import track: ok"
                    };
                    self.genome_track_status = format!(
                        "{}\nelapsed: {:.1}s",
                        Self::format_op_result_status(
                            prefix,
                            &result.created_seq_ids,
                            &result.warnings,
                            &result.messages,
                        ),
                        elapsed
                    );
                    self.push_job_event(
                        BackgroundJobKind::TrackImport,
                        BackgroundJobEventPhase::Completed,
                        Some(job_id),
                        format!("{prefix} in {:.1}s", elapsed),
                    );
                }
                Err(e) => {
                    self.genome_track_status =
                        format!("Import track failed after {:.1}s: {}", elapsed, e.message);
                    self.push_job_event(
                        BackgroundJobKind::TrackImport,
                        BackgroundJobEventPhase::Failed,
                        Some(job_id),
                        format!("Track import failed in {:.1}s: {}", elapsed, e.message),
                    );
                }
            }
        }
    }

    fn start_reference_genome_blast(&mut self) {
        if self.genome_blast_task.is_some() {
            self.genome_blast_status = "A BLAST task is already running".to_string();
            return;
        }

        let genome_id = self.genome_id.trim().to_string();
        if genome_id.is_empty() {
            self.genome_blast_status = "Select a prepared genome first".to_string();
            return;
        }

        let blast_queries: Vec<(String, String)> = match self.genome_blast_source_mode {
            GenomeBlastSourceMode::Manual => {
                let query = self.genome_blast_query_manual.trim().to_string();
                if query.is_empty() {
                    self.genome_blast_status =
                        "Provide a query sequence for manual BLAST mode".to_string();
                    return;
                }
                vec![("manual_query".to_string(), query)]
            }
            GenomeBlastSourceMode::ProjectSequence => {
                let seq_id = self.genome_blast_query_seq_id.trim().to_string();
                if seq_id.is_empty() {
                    self.genome_blast_status = "Select a project sequence to blast".to_string();
                    return;
                }
                let sequence = {
                    let engine = self.engine.read().unwrap();
                    let Some(dna) = engine.state().sequences.get(&seq_id) else {
                        self.genome_blast_status =
                            format!("Selected sequence '{}' is no longer available", seq_id);
                        return;
                    };
                    dna.get_forward_string()
                };
                if sequence.is_empty() {
                    self.genome_blast_status = format!(
                        "Selected sequence '{}' has no bases and cannot be blasted",
                        seq_id
                    );
                    return;
                }
                vec![(seq_id, sequence)]
            }
            GenomeBlastSourceMode::ProjectPool => {
                let pool_id = self.genome_blast_query_pool_id.trim().to_string();
                if pool_id.is_empty() {
                    self.genome_blast_status =
                        "Select a project pool/container to blast".to_string();
                    return;
                }
                let mut queries: Vec<(String, String)> = vec![];
                {
                    let engine = self.engine.read().unwrap();
                    let state = engine.state();
                    let Some(container) = state.container_state.containers.get(&pool_id) else {
                        self.genome_blast_status =
                            format!("Selected pool/container '{}' does not exist", pool_id);
                        return;
                    };
                    for seq_id in &container.members {
                        let Some(dna) = state.sequences.get(seq_id) else {
                            continue;
                        };
                        let seq = dna.get_forward_string();
                        if seq.is_empty() {
                            continue;
                        }
                        queries.push((seq_id.clone(), seq));
                    }
                }
                if queries.is_empty() {
                    self.genome_blast_status = format!(
                        "Selected pool/container '{}' has no blastable sequence members",
                        pool_id
                    );
                    return;
                }
                queries
            }
        };

        let total_queries = blast_queries.len();
        let catalog_path = self.genome_catalog_path_opt();
        let cache_dir = self.genome_cache_dir_opt();
        let max_hits = self.genome_blast_max_hits.max(1);
        let blast_task_name = self.genome_blast_task_name.trim().to_string();
        let task_arg = if blast_task_name.is_empty() {
            None
        } else {
            Some(blast_task_name)
        };
        let job_id = self.alloc_background_job_id();
        let (tx, rx) = mpsc::channel::<GenomeBlastTaskMessage>();
        self.genome_blast_results.clear();
        self.genome_blast_selected_result = 0;
        self.genome_blast_progress_fraction = Some(0.0);
        self.genome_blast_progress_label = format!("0 / {total_queries}");
        self.genome_blast_status = format!(
            "Running BLAST for {} quer{} in background",
            total_queries,
            if total_queries == 1 { "y" } else { "ies" }
        );
        self.push_job_event(
            BackgroundJobKind::BlastGenome,
            BackgroundJobEventPhase::Started,
            Some(job_id),
            format!(
                "BLAST started: genome='{}', queries={}, max_hits={}",
                genome_id, total_queries, max_hits
            ),
        );
        self.genome_blast_task = Some(GenomeBlastTask {
            job_id,
            started: Instant::now(),
            receiver: rx,
        });
        std::thread::spawn(move || {
            let mut reports: Vec<GenomeBlastQueryResult> = vec![];
            let mut failed_queries: Vec<String> = vec![];

            for (idx, (label, query)) in blast_queries.into_iter().enumerate() {
                let _ = tx.send(GenomeBlastTaskMessage::Progress {
                    job_id,
                    done_queries: idx,
                    total_queries,
                    current_query_label: label.clone(),
                });

                let query_length = query.len();
                match GentleEngine::blast_reference_genome(
                    catalog_path.as_deref(),
                    &genome_id,
                    &query,
                    max_hits,
                    task_arg.as_deref(),
                    cache_dir.as_deref(),
                ) {
                    Ok(report) => {
                        reports.push(GenomeBlastQueryResult {
                            query_label: label.clone(),
                            query_length,
                            report,
                        });
                    }
                    Err(e) => {
                        failed_queries.push(format!("{label}: {}", e.message));
                    }
                }
                let _ = tx.send(GenomeBlastTaskMessage::Progress {
                    job_id,
                    done_queries: idx + 1,
                    total_queries,
                    current_query_label: label,
                });
            }

            let done_result = if reports.is_empty() && !failed_queries.is_empty() {
                Err(format!(
                    "BLAST failed for all queries: {}",
                    failed_queries.join(" | ")
                ))
            } else {
                Ok(GenomeBlastBatchResult {
                    reports,
                    failed_queries,
                })
            };
            let _ = tx.send(GenomeBlastTaskMessage::Done {
                job_id,
                result: done_result,
            });
        });
    }

    fn poll_reference_genome_blast_task(&mut self, ctx: &egui::Context) {
        if self.genome_blast_task.is_none() {
            return;
        }
        ctx.request_repaint_after(Duration::from_millis(100));
        let mut done: Option<(u64, Result<GenomeBlastBatchResult, String>)> = None;
        let mut stale_job_ids: Vec<u64> = vec![];
        if let Some(task) = &self.genome_blast_task {
            let active_job_id = task.job_id;
            const MAX_MESSAGES_PER_TICK: usize = 128;
            for _ in 0..MAX_MESSAGES_PER_TICK {
                match task.receiver.try_recv() {
                    Ok(GenomeBlastTaskMessage::Progress {
                        job_id,
                        done_queries,
                        total_queries,
                        current_query_label,
                    }) => {
                        if job_id != active_job_id {
                            if !stale_job_ids.contains(&job_id) {
                                stale_job_ids.push(job_id);
                            }
                            continue;
                        }
                        let fraction = if total_queries == 0 {
                            0.0
                        } else {
                            (done_queries as f32 / total_queries as f32).clamp(0.0, 1.0)
                        };
                        self.genome_blast_progress_fraction = Some(fraction);
                        self.genome_blast_progress_label =
                            format!("{done_queries} / {total_queries} ({current_query_label})");
                    }
                    Ok(GenomeBlastTaskMessage::Done { job_id, result }) => {
                        if job_id != active_job_id {
                            if !stale_job_ids.contains(&job_id) {
                                stale_job_ids.push(job_id);
                            }
                            continue;
                        }
                        done = Some((job_id, result));
                        break;
                    }
                    Err(mpsc::TryRecvError::Empty) => break,
                    Err(mpsc::TryRecvError::Disconnected) => {
                        done = Some((active_job_id, Err("BLAST worker disconnected".to_string())));
                        break;
                    }
                }
            }
        }
        for stale_job_id in stale_job_ids {
            self.push_job_event(
                BackgroundJobKind::BlastGenome,
                BackgroundJobEventPhase::IgnoredStale,
                Some(stale_job_id),
                "Ignored stale BLAST worker message",
            );
        }

        if let Some((job_id, outcome)) = done {
            let elapsed = self
                .genome_blast_task
                .as_ref()
                .map(|task| task.started.elapsed().as_secs_f64())
                .unwrap_or(0.0);
            self.genome_blast_task = None;
            self.genome_blast_progress_fraction = None;
            self.genome_blast_progress_label.clear();

            match outcome {
                Ok(batch) => {
                    self.genome_blast_results = batch.reports;
                    if self.genome_blast_results.is_empty() {
                        self.genome_blast_selected_result = 0;
                    } else if self.genome_blast_selected_result >= self.genome_blast_results.len() {
                        self.genome_blast_selected_result = 0;
                    }
                    let hit_total: usize = self
                        .genome_blast_results
                        .iter()
                        .map(|r| r.report.hit_count)
                        .sum();
                    if batch.failed_queries.is_empty() {
                        self.genome_blast_status = format!(
                            "BLAST finished in {:.1}s: {} query result(s), {} total hit(s)",
                            elapsed,
                            self.genome_blast_results.len(),
                            hit_total
                        );
                        self.push_job_event(
                            BackgroundJobKind::BlastGenome,
                            BackgroundJobEventPhase::Completed,
                            Some(job_id),
                            format!(
                                "BLAST completed in {:.1}s: {} result(s), {} hit(s)",
                                elapsed,
                                self.genome_blast_results.len(),
                                hit_total
                            ),
                        );
                    } else {
                        self.genome_blast_status = format!(
                            "BLAST finished in {:.1}s: {} query result(s), {} total hit(s), {} failed query(ies): {}",
                            elapsed,
                            self.genome_blast_results.len(),
                            hit_total,
                            batch.failed_queries.len(),
                            batch.failed_queries.join(" | ")
                        );
                        self.push_job_event(
                            BackgroundJobKind::BlastGenome,
                            BackgroundJobEventPhase::Failed,
                            Some(job_id),
                            format!(
                                "BLAST completed in {:.1}s with {} failed quer{}",
                                elapsed,
                                batch.failed_queries.len(),
                                if batch.failed_queries.len() == 1 {
                                    "y"
                                } else {
                                    "ies"
                                }
                            ),
                        );
                    }
                }
                Err(e) => {
                    self.genome_blast_results.clear();
                    self.genome_blast_selected_result = 0;
                    self.genome_blast_status = format!("BLAST failed after {:.1}s: {}", elapsed, e);
                    self.push_job_event(
                        BackgroundJobKind::BlastGenome,
                        BackgroundJobEventPhase::Failed,
                        Some(job_id),
                        format!("BLAST failed in {:.1}s: {}", elapsed, e),
                    );
                }
            }
        }
    }

    fn target_seq_id_for_blast_result(&self, result: &GenomeBlastQueryResult) -> Option<String> {
        let candidate = result.query_label.trim();
        if candidate.is_empty() {
            return None;
        }
        if self
            .engine
            .read()
            .unwrap()
            .state()
            .sequences
            .contains_key(candidate)
        {
            Some(candidate.to_string())
        } else {
            None
        }
    }

    fn import_selected_blast_hits_as_track(&mut self) {
        let Some(result) = self
            .genome_blast_results
            .get(self.genome_blast_selected_result)
            .cloned()
        else {
            self.genome_blast_status = "No BLAST result is selected".to_string();
            return;
        };
        let Some(seq_id) = self.target_seq_id_for_blast_result(&result) else {
            self.genome_blast_status = format!(
                "Cannot import BLAST hits for query '{}' because it is not a current project sequence ID",
                result.query_label
            );
            return;
        };
        if result.report.hits.is_empty() {
            self.genome_blast_status = format!(
                "No BLAST hits available for '{}' to import",
                result.query_label
            );
            return;
        }

        let hits = result
            .report
            .hits
            .iter()
            .map(|hit| BlastHitFeatureInput {
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
        let track_name = {
            let trimmed = self.genome_blast_import_track_name.trim();
            if trimmed.is_empty() {
                None
            } else {
                Some(trimmed.to_string())
            }
        };
        let op = Operation::ImportBlastHitsTrack {
            seq_id: seq_id.clone(),
            hits,
            track_name,
            clear_existing: Some(self.genome_blast_import_clear_existing),
        };
        match self.engine.write().unwrap().apply(op) {
            Ok(result) => {
                self.genome_blast_status = Self::format_op_result_status(
                    &format!("Import BLAST hits to '{}': ok", seq_id),
                    &result.created_seq_ids,
                    &result.warnings,
                    &result.messages,
                );
            }
            Err(e) => {
                self.genome_blast_status = format!("Import BLAST hits failed: {}", e.message);
            }
        }
    }

    fn select_gene_record(&mut self, gene_index: usize) {
        let Some(gene) = self.genome_genes.get(gene_index) else {
            return;
        };
        self.genome_selected_gene = Some(gene_index);
        self.genome_chromosome = gene.chromosome.clone();
        self.genome_start_1based = gene.start_1based.to_string();
        self.genome_end_1based = gene.end_1based.to_string();
        if self.genome_output_id.trim().is_empty() {
            let base = gene
                .gene_name
                .as_ref()
                .or(gene.gene_id.as_ref())
                .cloned()
                .unwrap_or_else(|| "gene_region".to_string());
            self.genome_output_id = base.replace(' ', "_");
        }
    }

    fn gene_record_matches_regex(gene: &GenomeGeneRecord, regex: &Regex) -> bool {
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

    fn gene_record_matches_filter(gene: &GenomeGeneRecord, regex: Option<&Regex>) -> bool {
        if let Some(re) = regex {
            return Self::gene_record_matches_regex(gene, re);
        }
        true
    }

    fn gene_record_matches_biotype_filter(&self, gene: &GenomeGeneRecord) -> bool {
        if self.genome_biotype_filter.is_empty() {
            return true;
        }
        let Some(raw_biotype) = gene.biotype.as_ref() else {
            return true;
        };
        let biotype = raw_biotype.trim();
        if biotype.is_empty() {
            return true;
        }
        self.genome_biotype_filter
            .get(biotype)
            .copied()
            .unwrap_or(true)
    }

    fn gene_record_label(gene: &GenomeGeneRecord) -> String {
        let mut name = gene
            .gene_name
            .clone()
            .or(gene.gene_id.clone())
            .unwrap_or_else(|| "unnamed_gene".to_string());
        if let (Some(gene_name), Some(gene_id)) = (&gene.gene_name, &gene.gene_id) {
            if gene_name != gene_id {
                name = format!("{gene_name} ({gene_id})");
            }
        }
        let strand = gene
            .strand
            .map(|s| s.to_string())
            .unwrap_or_else(|| ".".to_string());
        let biotype = gene
            .biotype
            .as_ref()
            .map(|v| v.trim())
            .filter(|v| !v.is_empty())
            .unwrap_or("-");
        format!(
            "{}: {}:{}-{} ({}, {})",
            name, gene.chromosome, gene.start_1based, gene.end_1based, strand, biotype
        )
    }

    fn selected_gene_query_and_occurrence(&self) -> Option<(String, usize)> {
        let selected_idx = self.genome_selected_gene?;
        let selected_gene = self.genome_genes.get(selected_idx)?;
        let query = selected_gene
            .gene_id
            .as_ref()
            .or(selected_gene.gene_name.as_ref())
            .cloned()?;
        let mut occurrence = 0usize;
        for (idx, gene) in self.genome_genes.iter().enumerate() {
            let matches_query = gene
                .gene_id
                .as_ref()
                .map(|id| id.eq_ignore_ascii_case(&query))
                .unwrap_or(false)
                || gene
                    .gene_name
                    .as_ref()
                    .map(|name| name.eq_ignore_ascii_case(&query))
                    .unwrap_or(false);
            if matches_query {
                occurrence += 1;
            }
            if idx == selected_idx {
                break;
            }
        }
        if occurrence == 0 {
            None
        } else {
            Some((query, occurrence))
        }
    }

    fn filtered_gene_candidate_indices(&self) -> (Vec<usize>, bool, Option<String>) {
        let limit = self.genome_gene_filter_limit.max(1);
        let mut indices: Vec<usize> = Vec::with_capacity(limit.min(2048));
        let mut overflow = false;
        let needle = self.genome_gene_filter.trim();
        let regex = if needle.is_empty() {
            None
        } else {
            match RegexBuilder::new(needle).case_insensitive(true).build() {
                Ok(re) => Some(re),
                Err(e) => {
                    return (vec![], false, Some(format!("Invalid regex: {e}")));
                }
            }
        };
        for (idx, gene) in self.genome_genes.iter().enumerate() {
            if Self::gene_record_matches_filter(gene, regex.as_ref())
                && self.gene_record_matches_biotype_filter(gene)
            {
                if indices.len() < limit {
                    indices.push(idx);
                } else {
                    overflow = true;
                    break;
                }
            }
        }
        (indices, overflow, None)
    }

    fn normalize_coordinate_field(value: &mut String) {
        value.retain(|ch| ch.is_ascii_digit());
        if value.len() > 10 {
            value.truncate(10);
        }
    }

    fn extract_reference_genome_gene(&mut self) {
        let genome_id = self.genome_id.trim().to_string();
        if genome_id.is_empty() {
            self.genome_retrieve_status = "Select a genome first".to_string();
            return;
        }
        let Some((gene_query, occurrence)) = self.selected_gene_query_and_occurrence() else {
            self.genome_retrieve_status =
                "Select a gene with a gene_id or gene_name first".to_string();
            return;
        };
        let output_id = if self.genome_output_id.trim().is_empty() {
            None
        } else {
            Some(self.genome_output_id.trim().to_string())
        };
        let op = Operation::ExtractGenomeGene {
            genome_id,
            gene_query,
            occurrence: Some(occurrence),
            output_id,
            catalog_path: self.genome_catalog_path_opt(),
            cache_dir: self.genome_cache_dir_opt(),
        };
        let result = { self.engine.write().unwrap().apply(op) };
        match result {
            Ok(r) => {
                for seq_id in &r.created_seq_ids {
                    self.open_sequence_window(seq_id);
                }
                self.genome_retrieve_status = Self::format_op_result_status(
                    "Extract gene: ok",
                    &r.created_seq_ids,
                    &r.warnings,
                    &r.messages,
                );
            }
            Err(e) => {
                self.genome_retrieve_status = format!("Extract gene failed: {}", e.message);
            }
        }
    }

    fn extract_reference_genome_region(&mut self) {
        let genome_id = self.genome_id.trim().to_string();
        if genome_id.is_empty() {
            self.genome_retrieve_status = "Select a genome first".to_string();
            return;
        }
        let chromosome = self.genome_chromosome.trim().to_string();
        if chromosome.is_empty() {
            self.genome_retrieve_status = "Chromosome cannot be empty".to_string();
            return;
        }
        let (start_1based, end_1based) = match (
            self.genome_start_1based.trim().parse::<usize>(),
            self.genome_end_1based.trim().parse::<usize>(),
        ) {
            (Ok(start), Ok(end)) if start > 0 && end >= start => (start, end),
            _ => {
                self.genome_retrieve_status =
                    "Invalid interval; require start >= 1 and end >= start".to_string();
                return;
            }
        };
        let output_id = if self.genome_output_id.trim().is_empty() {
            None
        } else {
            Some(self.genome_output_id.trim().to_string())
        };
        let catalog_path = self.genome_catalog_path_opt();
        let extract_op = Operation::ExtractGenomeRegion {
            genome_id: genome_id.clone(),
            chromosome,
            start_1based,
            end_1based,
            output_id,
            catalog_path,
            cache_dir: self.genome_cache_dir_opt(),
        };
        let result = { self.engine.write().unwrap().apply(extract_op) };
        match result {
            Ok(r) => {
                for seq_id in &r.created_seq_ids {
                    self.open_sequence_window(seq_id);
                }
                self.genome_retrieve_status = Self::format_op_result_status(
                    "Extract region: ok",
                    &r.created_seq_ids,
                    &r.warnings,
                    &r.messages,
                );
            }
            Err(e) => {
                self.genome_retrieve_status = format!("Extract region failed: {}", e.message);
            }
        }
    }

    fn render_reference_genome_prepare_contents(&mut self, ui: &mut Ui) {
        self.refresh_genome_catalog_list();
        self.render_specialist_window_nav(ui);
        ui.label("Download and index a reference genome once.");
        ui.horizontal(|ui| {
            ui.label("catalog");
            ui.text_edit_singleline(&mut self.genome_catalog_path);
            if ui
                .button("Browse...")
                .on_hover_text("Browse filesystem and fill this path")
                .clicked()
            {
                if let Some(path) = rfd::FileDialog::new()
                    .add_filter("JSON", &["json"])
                    .pick_file()
                {
                    self.genome_catalog_path = path.display().to_string();
                }
            }
        });
        ui.horizontal(|ui| {
            ui.label("cache_dir");
            ui.text_edit_singleline(&mut self.genome_cache_dir);
            if ui
                .button("Browse...")
                .on_hover_text("Browse filesystem and fill this path")
                .clicked()
            {
                if let Some(path) = rfd::FileDialog::new().pick_folder() {
                    self.genome_cache_dir = path.display().to_string();
                }
            }
        });
        ui.horizontal(|ui| {
            ui.label("timeout_sec");
            ui.add(
                egui::TextEdit::singleline(&mut self.genome_prepare_timeout_secs)
                    .desired_width(90.0),
            );
            ui.small("optional; empty or 0 means no timebox");
        });
        if !self.genome_catalog_error.is_empty() {
            ui.colored_label(
                egui::Color32::from_rgb(190, 70, 70),
                format!("Catalog error: {}", self.genome_catalog_error),
            );
        }
        let all_genomes = self.genome_catalog_genomes.clone();
        let preparable_genomes = match self.unprepared_genomes_for_prepare_dialog() {
            Ok(names) => names,
            Err(e) => {
                ui.colored_label(
                    egui::Color32::from_rgb(190, 70, 70),
                    format!("Prepared-state check error: {e}"),
                );
                vec![]
            }
        };
        let preparable_set: HashSet<String> = preparable_genomes.iter().cloned().collect();
        let mut selection_changed = false;
        let selected_text = if self.genome_id.trim().is_empty() {
            "Select genome".to_string()
        } else {
            self.genome_id.clone()
        };
        egui::ComboBox::from_label("genome")
            .selected_text(selected_text)
            .show_ui(ui, |ui| {
                for genome_name in &all_genomes {
                    let is_preparable = preparable_set.contains(genome_name);
                    let item_label = if is_preparable {
                        genome_name.clone()
                    } else {
                        format!("{genome_name} (already prepared)")
                    };
                    ui.add_enabled_ui(is_preparable, |ui| {
                        if ui
                            .selectable_label(self.genome_id == *genome_name, item_label)
                            .clicked()
                        {
                            self.genome_id = genome_name.clone();
                            selection_changed = true;
                            ui.close_menu();
                        }
                    });
                }
            });
        if selection_changed {
            self.invalidate_genome_genes();
        }
        match self.selected_genome_source_plan() {
            Ok(Some(plan)) => {
                ui.small(Self::format_genome_source_plan_summary(&plan));
            }
            Ok(None) => {}
            Err(e) => {
                ui.colored_label(
                    egui::Color32::from_rgb(190, 70, 70),
                    format!("Source-plan error: {e}"),
                );
            }
        }
        if preparable_genomes.is_empty() {
            ui.label("All genomes in this catalog are already prepared.");
        }
        let selected_preparable = preparable_genomes.iter().any(|n| n == &self.genome_id);
        if !self.genome_id.trim().is_empty() && !selected_preparable {
            ui.label("Selected genome is already prepared and cannot be selected here.");
        }
        let running = self.genome_prepare_task.is_some();
        if let Some(task) = &self.genome_prepare_task {
            ui.horizontal(|ui| {
                ui.add(egui::Spinner::new());
                let mut status = format!(
                    "Prepare task running ({:.1}s)",
                    task.started.elapsed().as_secs_f32()
                );
                if let Some(timeout) = task.timeout_seconds {
                    status.push_str(&format!(", timeout={}s", timeout));
                }
                if task.cancel_requested.load(Ordering::Relaxed) {
                    status.push_str(", cancellation requested");
                }
                ui.label(status);
            });
        }
        ui.horizontal(|ui| {
            if ui
                .add_enabled(
                    !running && selected_preparable,
                    egui::Button::new("Prepare Genome"),
                )
                .on_hover_text("Download and index the selected reference genome")
                .clicked()
            {
                self.start_prepare_reference_genome();
            }
            if ui
                .button("Close")
                .on_hover_text("Close this dialog")
                .clicked()
            {
                self.show_reference_genome_prepare_dialog = false;
            }
            if self.genome_prepare_task.is_some() {
                if ui
                    .button("Cancel Prepare")
                    .on_hover_text("Request cancellation of the running prepare task.")
                    .clicked()
                {
                    self.request_prepare_task_cancel("prepare dialog");
                }
            }
        });
        if let Some(progress) = &self.genome_prepare_progress {
            let fraction = progress
                .percent
                .map(|p| (p / 100.0) as f32)
                .or_else(|| {
                    progress.bytes_total.and_then(|total| {
                        if total == 0 {
                            None
                        } else {
                            Some((progress.bytes_done as f32 / total as f32).clamp(0.0, 1.0))
                        }
                    })
                })
                .unwrap_or(0.0);
            ui.add(
                egui::ProgressBar::new(fraction)
                    .show_percentage()
                    .text(format!("{}: {}", progress.phase, progress.item)),
            );
            let bytes_total = progress
                .bytes_total
                .map(|b| b.to_string())
                .unwrap_or_else(|| "?".to_string());
            ui.label(format!("bytes: {} / {}", progress.bytes_done, bytes_total));
        }
        if !self.genome_prepare_status.is_empty() {
            ui.separator();
            ui.monospace(&self.genome_prepare_status);
        }
    }

    fn render_reference_genome_prepare_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_reference_genome_prepare_dialog {
            return;
        }

        let builder = egui::ViewportBuilder::default()
            .with_title("Prepare Reference Genome")
            .with_inner_size([760.0, 560.0])
            .with_min_inner_size([520.0, 360.0]);
        ctx.show_viewport_immediate(Self::prepare_genome_viewport_id(), builder, |ctx, class| {
            self.note_viewport_focus_if_active(ctx, Self::prepare_genome_viewport_id());
            if class == egui::ViewportClass::Embedded {
                let mut open = self.show_reference_genome_prepare_dialog;
                egui::Window::new("Prepare Reference Genome")
                    .open(&mut open)
                    .collapsible(false)
                    .resizable(true)
                    .default_size(Vec2::new(760.0, 560.0))
                    .show(ctx, |ui| {
                        self.render_reference_genome_prepare_contents(ui);
                    });
                self.show_reference_genome_prepare_dialog = open;
                return;
            }

            egui::CentralPanel::default().show(ctx, |ui| {
                self.render_reference_genome_prepare_contents(ui);
            });

            if ctx.input(|i| i.viewport().close_requested()) {
                self.show_reference_genome_prepare_dialog = false;
            }
        });
    }

    fn render_reference_genome_retrieve_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_reference_genome_retrieve_dialog {
            return;
        }
        self.refresh_genome_catalog_list();
        let mut open = self.show_reference_genome_retrieve_dialog;
        egui::Window::new("Retrieve Genomic Sequence")
            .open(&mut open)
            .collapsible(false)
            .resizable(true)
            .default_size(Vec2::new(980.0, 760.0))
            .min_size(Vec2::new(860.0, 520.0))
            .show(ctx, |ui| {
                self.render_specialist_window_nav(ui);
                ui.label("Retrieve region sequences from a prepared genome.");
                ui.horizontal(|ui| {
                    ui.label("catalog");
                    ui.text_edit_singleline(&mut self.genome_catalog_path);
                    if ui
                        .button("Browse...")
                        .on_hover_text("Browse filesystem and fill this path")
                        .clicked()
                    {
                        if let Some(path) = rfd::FileDialog::new()
                            .add_filter("JSON", &["json"])
                            .pick_file()
                        {
                            self.genome_catalog_path = path.display().to_string();
                        }
                    }
                });
                ui.horizontal(|ui| {
                    ui.label("cache_dir");
                    ui.text_edit_singleline(&mut self.genome_cache_dir);
                    if ui
                        .button("Browse...")
                        .on_hover_text("Browse filesystem and fill this path")
                        .clicked()
                    {
                        if let Some(path) = rfd::FileDialog::new().pick_folder() {
                            self.genome_cache_dir = path.display().to_string();
                        }
                    }
                });
                if !self.genome_catalog_error.is_empty() {
                    ui.colored_label(
                        egui::Color32::from_rgb(190, 70, 70),
                        format!("Catalog error: {}", self.genome_catalog_error),
                    );
                }
                let prepared_genomes = match self.prepared_genomes_for_retrieve_dialog() {
                    Ok(names) => names,
                    Err(e) => {
                        ui.colored_label(
                            egui::Color32::from_rgb(190, 70, 70),
                            format!("Prepared-state check error: {e}"),
                        );
                        vec![]
                    }
                };
                if !prepared_genomes.is_empty() && !prepared_genomes.contains(&self.genome_id) {
                    self.genome_id = prepared_genomes[0].clone();
                    self.invalidate_genome_genes();
                }
                let selection_changed =
                    Self::choose_genome_from_catalog(ui, &mut self.genome_id, &prepared_genomes);
                if selection_changed {
                    self.invalidate_genome_genes();
                }
                match self.selected_genome_source_plan() {
                    Ok(Some(plan)) => {
                        ui.small(Self::format_genome_source_plan_summary(&plan));
                    }
                    Ok(None) => {}
                    Err(e) => {
                        ui.colored_label(
                            egui::Color32::from_rgb(190, 70, 70),
                            format!("Source-plan error: {e}"),
                        );
                    }
                }
                if prepared_genomes.is_empty() {
                    ui.label(
                        "No prepared genomes are available in this catalog/cache. Use 'Prepare Reference Genome...' first.",
                    );
                } else {
                    self.ensure_genome_genes_loaded();
                    if !self.genome_genes_error.is_empty() {
                        ui.colored_label(
                            egui::Color32::from_rgb(190, 70, 70),
                            format!("Gene index error: {}", self.genome_genes_error),
                        );
                    }
                    ui.separator();
                    ui.horizontal(|ui| {
                        ui.label("Gene filter");
                        let response = ui.text_edit_singleline(&mut self.genome_gene_filter);
                        if response.changed() {
                            self.genome_gene_filter_page = 0;
                        }
                        if ui
                            .button("Clear")
                            .on_hover_text("Clear the gene filter text")
                            .clicked()
                        {
                            self.genome_gene_filter.clear();
                            self.genome_gene_filter_page = 0;
                        }
                    });
                    ui.small("Supports case-insensitive regex (for example: ^TP53$).");
                    ui.horizontal(|ui| {
                        ui.label("Top matches");
                        if ui
                            .add(
                                egui::DragValue::new(&mut self.genome_gene_filter_limit)
                                    .speed(100.0)
                                    .range(100..=100_000),
                            )
                            .changed()
                        {
                            self.genome_gene_filter_page = 0;
                        }
                    });
                    if !self.genome_biotype_filter.is_empty() {
                        ui.separator();
                        ui.label("Biotype filter");
                        ui.horizontal(|ui| {
                            if ui
                                .button("All")
                                .on_hover_text("Enable all biotypes")
                                .clicked()
                            {
                                for enabled in self.genome_biotype_filter.values_mut() {
                                    *enabled = true;
                                }
                                self.genome_gene_filter_page = 0;
                            }
                            if ui
                                .button("None")
                                .on_hover_text("Disable all biotypes")
                                .clicked()
                            {
                                for enabled in self.genome_biotype_filter.values_mut() {
                                    *enabled = false;
                                }
                                self.genome_gene_filter_page = 0;
                            }
                        });
                        let mut biotypes: Vec<String> =
                            self.genome_biotype_filter.keys().cloned().collect();
                        biotypes.sort();
                        egui::ScrollArea::vertical()
                            .max_height(120.0)
                            .show(ui, |ui| {
                                for biotype in biotypes {
                                    if let Some(enabled) =
                                        self.genome_biotype_filter.get_mut(&biotype)
                                    {
                                        if ui.checkbox(enabled, &biotype).changed() {
                                            self.genome_gene_filter_page = 0;
                                        }
                                    }
                                }
                            });
                    }
                    let (filtered_indices, overflow, regex_error) =
                        self.filtered_gene_candidate_indices();
                    if let Some(error) = regex_error {
                        ui.colored_label(egui::Color32::from_rgb(190, 70, 70), error);
                    }
                    const GENE_PAGE_SIZE: usize = 100;
                    let page_count =
                        ((filtered_indices.len() + GENE_PAGE_SIZE - 1) / GENE_PAGE_SIZE).max(1);
                    if self.genome_gene_filter_page >= page_count {
                        self.genome_gene_filter_page = page_count.saturating_sub(1);
                    }
                    let page_start = self.genome_gene_filter_page * GENE_PAGE_SIZE;
                    let page_end = (page_start + GENE_PAGE_SIZE).min(filtered_indices.len());
                    let shown_note = if overflow {
                        format!(
                            " (limited to top {} matches)",
                            self.genome_gene_filter_limit
                        )
                    } else {
                        String::new()
                    };
                    ui.label(format!(
                        "Genes: {} candidate{} / {} total",
                        filtered_indices.len(),
                        shown_note,
                        self.genome_genes.len()
                    ));
                    ui.horizontal(|ui| {
                        if ui
                            .add_enabled(
                                self.genome_gene_filter_page > 0,
                                egui::Button::new("Prev"),
                            )
                            .on_hover_text("Show previous page of gene matches")
                            .clicked()
                        {
                            self.genome_gene_filter_page -= 1;
                        }
                        ui.label(format!(
                            "Page {}/{}",
                            self.genome_gene_filter_page + 1,
                            page_count
                        ));
                        if ui
                            .add_enabled(
                                self.genome_gene_filter_page + 1 < page_count,
                                egui::Button::new("Next"),
                            )
                            .on_hover_text("Show next page of gene matches")
                            .clicked()
                        {
                            self.genome_gene_filter_page += 1;
                        }
                    });
                    egui::ScrollArea::vertical()
                        .max_height(220.0)
                        .show(ui, |ui| {
                            for idx in filtered_indices
                                .iter()
                                .copied()
                                .skip(page_start)
                                .take(page_end.saturating_sub(page_start))
                            {
                                let label = Self::gene_record_label(&self.genome_genes[idx]);
                                let selected = self.genome_selected_gene == Some(idx);
                                if ui.selectable_label(selected, label).clicked() {
                                    self.select_gene_record(idx);
                                }
                            }
                        });
                    ui.separator();
                    ui.horizontal(|ui| {
                        ui.label("chr");
                        ui.text_edit_singleline(&mut self.genome_chromosome);
                        ui.label("start_1based");
                        let start_changed = ui
                            .add(
                                egui::TextEdit::singleline(&mut self.genome_start_1based)
                                    .desired_width(120.0)
                                    .char_limit(10),
                            )
                            .changed();
                        if start_changed {
                            Self::normalize_coordinate_field(&mut self.genome_start_1based);
                        }
                        ui.label("end_1based");
                        let end_changed = ui
                            .add(
                                egui::TextEdit::singleline(&mut self.genome_end_1based)
                                    .desired_width(120.0)
                                    .char_limit(10),
                            )
                            .changed();
                        if end_changed {
                            Self::normalize_coordinate_field(&mut self.genome_end_1based);
                        }
                    });
                    ui.horizontal(|ui| {
                        ui.label("output_id");
                        ui.text_edit_singleline(&mut self.genome_output_id);
                        if ui
                            .add_enabled(
                                self.genome_selected_gene.is_some(),
                                egui::Button::new("Extract Selected Gene"),
                            )
                            .on_hover_text(
                                "Extract the currently selected gene from the prepared reference",
                            )
                            .clicked()
                        {
                            self.extract_reference_genome_gene();
                        }
                        if ui
                            .button("Extract Region")
                            .on_hover_text("Extract the explicit chromosome start/end interval")
                            .clicked()
                        {
                            self.extract_reference_genome_region();
                        }
                    });
                }
                if !self.genome_retrieve_status.is_empty() {
                    ui.separator();
                    ui.monospace(&self.genome_retrieve_status);
                }
            });
        self.show_reference_genome_retrieve_dialog = open;
    }

    fn render_blast_query_result(ui: &mut Ui, result: &GenomeBlastQueryResult) {
        let report = &result.report;
        ui.label(format!(
            "Query: {} ({} bp) | Genome: {} | Hits: {}",
            result.query_label, result.query_length, report.genome_id, report.hit_count
        ));
        ui.label(format!(
            "Task: {} | max_hits: {}",
            report.task, report.max_hits
        ));
        ui.monospace(format!("blastn: {}", report.blastn_executable));
        ui.monospace(format!("db: {}", report.blast_db_prefix));
        if !report.command.is_empty() {
            ui.monospace(format!("command: {}", report.command.join(" ")));
        }
        if !report.warnings.is_empty() {
            ui.separator();
            ui.label("Warnings");
            for warning in &report.warnings {
                ui.monospace(format!("- {}", warning));
            }
        }
        if !report.stderr.trim().is_empty() {
            ui.separator();
            ui.label("BLAST stderr");
            ui.monospace(report.stderr.trim().to_string());
        }
        ui.separator();
        ui.label("Hits");
        let max_rows = 200usize;
        let shown_hits = report.hits.len().min(max_rows);
        if shown_hits < report.hits.len() {
            ui.small(format!(
                "Showing first {} of {} hits",
                shown_hits,
                report.hits.len()
            ));
        }
        egui::ScrollArea::both().max_height(280.0).show(ui, |ui| {
            egui::Grid::new(format!("blast_hits_{}", result.query_label))
                .striped(true)
                .show(ui, |ui| {
                    ui.strong("subject");
                    ui.strong("%identity");
                    ui.strong("aln_len");
                    ui.strong("q_range");
                    ui.strong("s_range");
                    ui.strong("evalue");
                    ui.strong("bit_score");
                    ui.strong("qcov%");
                    ui.end_row();
                    for hit in report.hits.iter().take(max_rows) {
                        ui.monospace(hit.subject_id.clone());
                        ui.label(format!("{:.2}", hit.identity_percent));
                        ui.label(hit.alignment_length.to_string());
                        ui.label(format!("{}-{}", hit.query_start, hit.query_end));
                        ui.label(format!("{}-{}", hit.subject_start, hit.subject_end));
                        ui.label(format!("{:.3e}", hit.evalue));
                        ui.label(format!("{:.2}", hit.bit_score));
                        ui.label(
                            hit.query_coverage_percent
                                .map(|v| format!("{v:.1}"))
                                .unwrap_or_else(|| "-".to_string()),
                        );
                        ui.end_row();
                    }
                });
        });
    }

    fn render_reference_genome_blast_contents(&mut self, ui: &mut Ui) {
        self.refresh_genome_catalog_list();
        self.render_specialist_window_nav(ui);
        let sequence_ids = self.project_sequence_ids_for_blast();
        if !sequence_ids.is_empty() && !sequence_ids.contains(&self.genome_blast_query_seq_id) {
            self.genome_blast_query_seq_id = sequence_ids[0].clone();
        }
        let pool_options = self.blast_pool_options();
        if !pool_options
            .iter()
            .any(|pool| pool.container_id == self.genome_blast_query_pool_id)
        {
            self.genome_blast_query_pool_id = pool_options
                .first()
                .map(|pool| pool.container_id.clone())
                .unwrap_or_default();
        }

        ui.label("Run BLAST against a prepared genome index.");
        ui.horizontal(|ui| {
            ui.label("catalog");
            ui.text_edit_singleline(&mut self.genome_catalog_path);
            if ui
                .button("Browse...")
                .on_hover_text("Browse filesystem and fill this path")
                .clicked()
            {
                if let Some(path) = rfd::FileDialog::new()
                    .add_filter("JSON", &["json"])
                    .pick_file()
                {
                    self.genome_catalog_path = path.display().to_string();
                }
            }
        });
        ui.horizontal(|ui| {
            ui.label("cache_dir");
            ui.text_edit_singleline(&mut self.genome_cache_dir);
            if ui
                .button("Browse...")
                .on_hover_text("Browse filesystem and fill this path")
                .clicked()
            {
                if let Some(path) = rfd::FileDialog::new().pick_folder() {
                    self.genome_cache_dir = path.display().to_string();
                }
            }
        });
        if !self.genome_catalog_error.is_empty() {
            ui.colored_label(
                egui::Color32::from_rgb(190, 70, 70),
                format!("Catalog error: {}", self.genome_catalog_error),
            );
        }

        let prepared_genomes = match self.prepared_genomes_for_retrieve_dialog() {
            Ok(names) => names,
            Err(e) => {
                ui.colored_label(
                    egui::Color32::from_rgb(190, 70, 70),
                    format!("Prepared-state check error: {e}"),
                );
                vec![]
            }
        };
        if !prepared_genomes.is_empty() && !prepared_genomes.contains(&self.genome_id) {
            self.genome_id = prepared_genomes[0].clone();
            self.invalidate_genome_genes();
        }
        let selection_changed =
            Self::choose_genome_from_catalog(ui, &mut self.genome_id, &prepared_genomes);
        if selection_changed {
            self.invalidate_genome_genes();
        }
        match self.selected_genome_source_plan() {
            Ok(Some(plan)) => {
                ui.small(Self::format_genome_source_plan_summary(&plan));
            }
            Ok(None) => {}
            Err(e) => {
                ui.colored_label(
                    egui::Color32::from_rgb(190, 70, 70),
                    format!("Source-plan error: {e}"),
                );
            }
        }
        if prepared_genomes.is_empty() {
            ui.label(
                "No prepared genomes are available in this catalog/cache. Use 'Prepare Reference Genome...' first.",
            );
        }

        ui.separator();
        ui.label("Query source");
        ui.horizontal(|ui| {
            ui.selectable_value(
                &mut self.genome_blast_source_mode,
                GenomeBlastSourceMode::Manual,
                "Manual sequence",
            );
            ui.selectable_value(
                &mut self.genome_blast_source_mode,
                GenomeBlastSourceMode::ProjectSequence,
                "Project sequence",
            );
            ui.selectable_value(
                &mut self.genome_blast_source_mode,
                GenomeBlastSourceMode::ProjectPool,
                "Project pool",
            );
        });
        match self.genome_blast_source_mode {
            GenomeBlastSourceMode::Manual => {
                ui.label("Query sequence (IUPAC DNA letters)");
                ui.add(
                    egui::TextEdit::multiline(&mut self.genome_blast_query_manual)
                        .desired_rows(4)
                        .hint_text("Paste query sequence here"),
                );
            }
            GenomeBlastSourceMode::ProjectSequence => {
                if sequence_ids.is_empty() {
                    ui.label("No sequences are loaded in this project.");
                } else {
                    egui::ComboBox::from_label("sequence")
                        .selected_text(if self.genome_blast_query_seq_id.trim().is_empty() {
                            "(choose sequence)"
                        } else {
                            self.genome_blast_query_seq_id.as_str()
                        })
                        .show_ui(ui, |ui| {
                            for seq_id in &sequence_ids {
                                ui.selectable_value(
                                    &mut self.genome_blast_query_seq_id,
                                    seq_id.clone(),
                                    seq_id,
                                );
                            }
                        });
                }
            }
            GenomeBlastSourceMode::ProjectPool => {
                if pool_options.is_empty() {
                    ui.label("No sequence pools/containers with >1 sequence are available.");
                } else {
                    egui::ComboBox::from_label("pool/container")
                        .selected_text(if self.genome_blast_query_pool_id.trim().is_empty() {
                            "(choose pool)"
                        } else {
                            self.genome_blast_query_pool_id.as_str()
                        })
                        .show_ui(ui, |ui| {
                            for pool in &pool_options {
                                ui.selectable_value(
                                    &mut self.genome_blast_query_pool_id,
                                    pool.container_id.clone(),
                                    pool.label.clone(),
                                );
                            }
                        });
                    if let Some(pool) = pool_options
                        .iter()
                        .find(|pool| pool.container_id == self.genome_blast_query_pool_id)
                    {
                        ui.small(format!(
                            "Selected pool will run {} query sequence(s).",
                            pool.members.len()
                        ));
                    }
                }
            }
        }

        ui.separator();
        ui.horizontal(|ui| {
            ui.label("max_hits");
            ui.add(
                egui::DragValue::new(&mut self.genome_blast_max_hits)
                    .range(1..=1000)
                    .speed(1.0),
            );
            ui.label("task");
            egui::ComboBox::from_id_salt("genome_blast_task_combo")
                .selected_text(self.genome_blast_task_name.clone())
                .show_ui(ui, |ui| {
                    ui.selectable_value(
                        &mut self.genome_blast_task_name,
                        "blastn-short".to_string(),
                        "blastn-short",
                    );
                    ui.selectable_value(
                        &mut self.genome_blast_task_name,
                        "blastn".to_string(),
                        "blastn",
                    );
                });
        });

        let running = self.genome_blast_task.is_some();
        ui.horizontal(|ui| {
            if ui
                .add_enabled(
                    !running && !prepared_genomes.is_empty(),
                    egui::Button::new("Run BLAST"),
                )
                .on_hover_text("Run BLAST query and collect hits for the selected prepared genome")
                .clicked()
            {
                self.start_reference_genome_blast();
            }
            if ui
                .button("Clear Results")
                .on_hover_text("Clear all BLAST query results from this dialog")
                .clicked()
            {
                self.genome_blast_results.clear();
                self.genome_blast_selected_result = 0;
            }
            if ui
                .button("Close")
                .on_hover_text("Close this dialog")
                .clicked()
            {
                self.show_reference_genome_blast_dialog = false;
            }
        });

        if let Some(fraction) = self.genome_blast_progress_fraction {
            ui.add(
                egui::ProgressBar::new(fraction.clamp(0.0, 1.0))
                    .show_percentage()
                    .text(self.genome_blast_progress_label.clone()),
            );
        }
        if !self.genome_blast_status.is_empty() {
            ui.separator();
            ui.monospace(self.genome_blast_status.clone());
        }
        if !self.genome_blast_results.is_empty() {
            ui.separator();
            if self.genome_blast_selected_result >= self.genome_blast_results.len() {
                self.genome_blast_selected_result = 0;
            }
            if self.genome_blast_results.len() > 1 {
                egui::ComboBox::from_label("Query result")
                    .selected_text(
                        self.genome_blast_results[self.genome_blast_selected_result]
                            .query_label
                            .clone(),
                    )
                    .show_ui(ui, |ui| {
                        for (idx, result) in self.genome_blast_results.iter().enumerate() {
                            let label = format!(
                                "{} (hits={})",
                                result.query_label, result.report.hit_count
                            );
                            ui.selectable_value(&mut self.genome_blast_selected_result, idx, label);
                        }
                    });
            }
            if let Some(result) = self
                .genome_blast_results
                .get(self.genome_blast_selected_result)
                .cloned()
            {
                let target_seq_id = self.target_seq_id_for_blast_result(&result);
                ui.separator();
                ui.label("Import BLAST hits as features");
                ui.horizontal(|ui| {
                    ui.label("track_name");
                    ui.text_edit_singleline(&mut self.genome_blast_import_track_name);
                    ui.checkbox(
                        &mut self.genome_blast_import_clear_existing,
                        "Clear existing BLAST-hit features first",
                    );
                });
                let can_import = target_seq_id.is_some() && !result.report.hits.is_empty();
                let button = ui
                    .add_enabled(can_import, egui::Button::new("Import Hits To Sequence"))
                    .on_hover_text(
                        "Import current BLAST hit intervals as features on the target sequence",
                    );
                if button.clicked() {
                    self.import_selected_blast_hits_as_track();
                }
                if let Some(seq_id) = target_seq_id {
                    ui.small(format!("Import target sequence: {seq_id}"));
                } else {
                    ui.small(
                        "Hit import requires query label to match an existing project sequence ID.",
                    );
                }
                Self::render_blast_query_result(ui, &result);
            }
        }
    }

    fn render_reference_genome_blast_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_reference_genome_blast_dialog {
            return;
        }

        let builder = egui::ViewportBuilder::default()
            .with_title("BLAST Genome Sequence")
            .with_inner_size([980.0, 700.0])
            .with_min_inner_size([640.0, 420.0]);
        ctx.show_viewport_immediate(Self::blast_genome_viewport_id(), builder, |ctx, class| {
            self.note_viewport_focus_if_active(ctx, Self::blast_genome_viewport_id());
            if class == egui::ViewportClass::Embedded {
                let mut open = self.show_reference_genome_blast_dialog;
                egui::Window::new("BLAST Genome Sequence")
                    .open(&mut open)
                    .collapsible(false)
                    .resizable(true)
                    .default_size(Vec2::new(980.0, 700.0))
                    .show(ctx, |ui| {
                        self.render_reference_genome_blast_contents(ui);
                    });
                self.show_reference_genome_blast_dialog = open;
                return;
            }

            egui::CentralPanel::default().show(ctx, |ui| {
                self.render_reference_genome_blast_contents(ui);
            });

            if ctx.input(|i| i.viewport().close_requested()) {
                self.show_reference_genome_blast_dialog = false;
            }
        });
    }

    fn render_reference_genome_inspector_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_reference_genome_inspector_dialog {
            return;
        }
        self.refresh_genome_catalog_list();
        let mut open = self.show_reference_genome_inspector_dialog;
        egui::Window::new("Prepared Genome References")
            .open(&mut open)
            .collapsible(false)
            .resizable(true)
            .show(ctx, |ui| {
                self.render_specialist_window_nav(ui);
                ui.label("Inspect prepared references and installation integrity metadata.");
                ui.horizontal(|ui| {
                    ui.label("catalog");
                    ui.text_edit_singleline(&mut self.genome_catalog_path);
                    if ui
                        .button("Browse...")
                        .on_hover_text("Browse filesystem and fill this path")
                        .clicked()
                    {
                        if let Some(path) = rfd::FileDialog::new()
                            .add_filter("JSON", &["json"])
                            .pick_file()
                        {
                            self.genome_catalog_path = path.display().to_string();
                        }
                    }
                });
                ui.horizontal(|ui| {
                    ui.label("cache_dir");
                    ui.text_edit_singleline(&mut self.genome_cache_dir);
                    if ui
                        .button("Browse...")
                        .on_hover_text("Browse filesystem and fill this path")
                        .clicked()
                    {
                        if let Some(path) = rfd::FileDialog::new().pick_folder() {
                            self.genome_cache_dir = path.display().to_string();
                        }
                    }
                });
                if !self.genome_catalog_error.is_empty() {
                    ui.colored_label(
                        egui::Color32::from_rgb(190, 70, 70),
                        format!("Catalog error: {}", self.genome_catalog_error),
                    );
                }
                match self.collect_prepared_genome_inspections() {
                    Ok((inspections, errors)) => {
                        let total_size: u64 = inspections.iter().map(|r| r.total_size_bytes).sum();
                        ui.label(format!(
                            "Prepared references: {} | total size: {}",
                            inspections.len(),
                            Self::format_bytes_compact(total_size)
                        ));
                        if inspections.is_empty() {
                            ui.label("No prepared references found for this catalog/cache.");
                        } else {
                            egui::ScrollArea::vertical()
                                .max_height(320.0)
                                .show(ui, |ui| {
                                    egui::Grid::new("prepared_genome_inspector_grid")
                                        .striped(true)
                                        .num_columns(8)
                                        .show(ui, |ui| {
                                            ui.strong("Genome");
                                            ui.strong("Size");
                                            ui.strong("Ready");
                                            ui.strong("Sources");
                                            ui.strong("SHA1 seq/ann");
                                            ui.strong("Installed");
                                            ui.strong("Path");
                                            ui.strong("");
                                            ui.end_row();
                                            for inspection in &inspections {
                                                ui.label(&inspection.genome_id);
                                                ui.label(Self::format_bytes_compact(
                                                    inspection.total_size_bytes,
                                                ));
                                                ui.label(format!(
                                                    "seq:{} ann:{} fai:{} gene:{}",
                                                    if inspection.sequence_present {
                                                        "y"
                                                    } else {
                                                        "n"
                                                    },
                                                    if inspection.annotation_present {
                                                        "y"
                                                    } else {
                                                        "n"
                                                    },
                                                    if inspection.fasta_index_ready {
                                                        "y"
                                                    } else {
                                                        "n"
                                                    },
                                                    if inspection.gene_index_ready {
                                                        "y"
                                                    } else {
                                                        "n"
                                                    }
                                                ));
                                                ui.label(format!(
                                                    "{}/{}",
                                                    inspection.sequence_source_type,
                                                    inspection.annotation_source_type
                                                ));
                                                ui.label(format!(
                                                    "{}/{}",
                                                    Self::format_short_sha1(
                                                        &inspection.sequence_sha1
                                                    ),
                                                    Self::format_short_sha1(
                                                        &inspection.annotation_sha1
                                                    )
                                                ));
                                                ui.label(
                                                    inspection.installed_at_unix_ms.to_string(),
                                                );
                                                let path_label = inspection.install_dir.clone();
                                                ui.label(
                                                    egui::RichText::new(path_label)
                                                        .monospace()
                                                        .small(),
                                                );
                                                if ui
                                                    .small_button("Retrieve")
                                                    .on_hover_text(
                                                        "Open retrieval dialog preselected for this prepared genome",
                                                    )
                                                    .clicked()
                                                {
                                                    self.genome_id = inspection.genome_id.clone();
                                                    self.invalidate_genome_genes();
                                                    self.open_reference_genome_retrieve_dialog();
                                                }
                                                ui.end_row();
                                            }
                                        });
                                });

                            ui.separator();
                            ui.label("Chromosome inspector");
                            let prepared_ids: Vec<String> =
                                inspections.iter().map(|i| i.genome_id.clone()).collect();
                            if !prepared_ids.is_empty()
                                && !prepared_ids.iter().any(|id| id == &self.genome_id)
                            {
                                self.genome_id = prepared_ids[0].clone();
                            }
                            ui.horizontal(|ui| {
                                ui.label("genome");
                                egui::ComboBox::from_id_salt("prepared_genome_chromosome_inspect")
                                    .selected_text(if self.genome_id.trim().is_empty() {
                                        "<select genome>"
                                    } else {
                                        self.genome_id.as_str()
                                    })
                                    .show_ui(ui, |ui| {
                                        for genome_id in &prepared_ids {
                                            ui.selectable_value(
                                                &mut self.genome_id,
                                                genome_id.clone(),
                                                genome_id,
                                            );
                                        }
                                    });
                                if ui
                                    .small_button("Use In Retrieve")
                                    .on_hover_text(
                                        "Open retrieval dialog preselected for this genome",
                                    )
                                    .clicked()
                                {
                                    self.invalidate_genome_genes();
                                    self.open_reference_genome_retrieve_dialog();
                                }
                            });
                            if !self.genome_id.trim().is_empty() {
                                match self.prepared_genome_chromosome_records(&self.genome_id) {
                                    Ok(chromosomes) => {
                                        Self::render_chromosome_length_lines(ui, &chromosomes);
                                    }
                                    Err(e) => {
                                        ui.colored_label(
                                            egui::Color32::from_rgb(190, 70, 70),
                                            format!("Chromosome inspector error: {e}"),
                                        );
                                    }
                                }
                            }
                        }
                        if !errors.is_empty() {
                            ui.separator();
                            ui.colored_label(
                                egui::Color32::from_rgb(190, 70, 70),
                                "Inspection errors:",
                            );
                            for err in errors {
                                ui.colored_label(egui::Color32::from_rgb(190, 70, 70), err);
                            }
                        }
                    }
                    Err(e) => {
                        ui.colored_label(
                            egui::Color32::from_rgb(190, 70, 70),
                            format!("Inspector error: {e}"),
                        );
                    }
                }
            });
        self.show_reference_genome_inspector_dialog = open;
    }

    fn render_genome_bed_track_contents(&mut self, ui: &mut Ui) {
        self.render_specialist_window_nav(ui);
        let anchor_summaries = self.anchored_sequence_anchor_summaries_for_tracks();
        let anchored_seq_ids = anchor_summaries
            .iter()
            .map(|summary| summary.seq_id.clone())
            .collect::<Vec<_>>();
        let import_running = self.genome_track_import_task.is_some();
        if !anchored_seq_ids.is_empty() && !anchored_seq_ids.contains(&self.genome_track_seq_id) {
            self.genome_track_seq_id = anchored_seq_ids[0].clone();
        }

        ui.label("Map genome signal tracks onto genome-anchored sequences.");
        ui.small("Supports .bed/.bed.gz and .bw/.bigWig input files.");
        ui.small(format!(
            "Anchored sequences available: {}",
            anchored_seq_ids.len()
        ));

        if anchored_seq_ids.is_empty() {
            ui.colored_label(
                egui::Color32::from_rgb(190, 70, 70),
                "No genome-anchored sequence is available. Extract a genome region or gene first.",
            );
        } else {
            ui.horizontal(|ui| {
                ui.label("sequence");
                egui::ComboBox::from_id_salt("genome_track_seq_id_combo")
                    .selected_text(if self.genome_track_seq_id.trim().is_empty() {
                        "<select sequence>"
                    } else {
                        self.genome_track_seq_id.as_str()
                    })
                    .show_ui(ui, |ui| {
                        for seq_id in &anchored_seq_ids {
                            ui.selectable_value(
                                &mut self.genome_track_seq_id,
                                seq_id.clone(),
                                seq_id,
                            );
                        }
                    });
                let selected_resp = self.track_hover_status(
                    ui.add_enabled(
                        !anchored_seq_ids.is_empty() && !import_running,
                        egui::Button::new("Import To Selected"),
                    )
                    .on_hover_text(
                        "Import this BED/BigWig/VCF signal file onto only the currently selected anchored sequence.",
                    ),
                    "Genome Tracks > Import Selected",
                );
                if selected_resp.clicked() {
                    self.import_genome_bed_track_for_selected_sequence();
                }
            });
            if let Some(anchor) = self.describe_sequence_genome_anchor(&self.genome_track_seq_id) {
                ui.small(format!("selected anchor: {anchor}"));
            }
            ui.small(
                "Tip: Keep this window open. Change 'sequence' and click 'Import To Selected' again to reuse the same track file for another anchored sequence.",
            );
        }

        let detected_source = GenomeTrackSource::from_path(&self.genome_track_path);
        ui.horizontal(|ui| {
            ui.label("source");
            egui::ComboBox::from_id_salt("genome_track_source_selection")
                .selected_text(self.genome_track_source_selection.label())
                .show_ui(ui, |ui| {
                    ui.selectable_value(
                        &mut self.genome_track_source_selection,
                        GenomeTrackSourceSelection::Auto,
                        GenomeTrackSourceSelection::Auto.label(),
                    );
                    ui.selectable_value(
                        &mut self.genome_track_source_selection,
                        GenomeTrackSourceSelection::Bed,
                        GenomeTrackSourceSelection::Bed.label(),
                    );
                    ui.selectable_value(
                        &mut self.genome_track_source_selection,
                        GenomeTrackSourceSelection::BigWig,
                        GenomeTrackSourceSelection::BigWig.label(),
                    );
                    ui.selectable_value(
                        &mut self.genome_track_source_selection,
                        GenomeTrackSourceSelection::Vcf,
                        GenomeTrackSourceSelection::Vcf.label(),
                    );
                });
            ui.small(format!("Detected extension: {}", detected_source.label()));
        });

        ui.horizontal(|ui| {
            ui.label("track_path");
            ui.text_edit_singleline(&mut self.genome_track_path);
            let browse_track_resp = self.track_hover_status(
                ui.button("Browse...")
                    .on_hover_text("Browse filesystem and fill this path"),
                "Genome Tracks > Browse Track Path",
            );
            if browse_track_resp.clicked() {
                if let Some(path) = rfd::FileDialog::new()
                    .add_filter("Signal tracks", &["bed", "gz", "bw", "bigwig", "vcf"])
                    .pick_file()
                {
                    self.genome_track_path = path.display().to_string();
                    if self.genome_track_name.trim().is_empty() {
                        self.genome_track_name = Self::track_name_default_from_path(&path);
                    }
                }
            }
        });
        let resolved_source = self
            .genome_track_source_selection
            .resolve(&self.genome_track_path);
        ui.small(format!(
            "Resolved source format for import: {}",
            resolved_source.label()
        ));
        ui.horizontal(|ui| {
            ui.label("track_name");
            ui.text_edit_singleline(&mut self.genome_track_name);
            ui.small("(optional)");
        });
        ui.horizontal(|ui| {
            ui.label("min_score");
            ui.add(
                egui::TextEdit::singleline(&mut self.genome_track_min_score).desired_width(90.0),
            );
            ui.label("max_score");
            ui.add(
                egui::TextEdit::singleline(&mut self.genome_track_max_score).desired_width(90.0),
            );
        });
        ui.checkbox(
            &mut self.genome_track_clear_existing,
            "Clear existing imported track features first",
        );
        if let Some(selected_anchor) = anchor_summaries
            .iter()
            .find(|summary| summary.seq_id == self.genome_track_seq_id)
        {
            let matching_anchor_count = anchor_summaries
                .iter()
                .filter(|candidate| Self::anchors_share_mapping_group(candidate, selected_anchor))
                .count();
            let projected_targets = anchored_seq_ids.len();
            let mapping_status = if matching_anchor_count == projected_targets {
                "all anchored sequences match selected genome/chromosome"
            } else {
                "mixed anchors detected; import still applies to all anchored sequences"
            };
            let track_path = self.genome_track_path.trim().to_string();
            let path_exists = !track_path.is_empty() && Path::new(&track_path).is_file();
            let projected_track_name = {
                let trimmed = self.genome_track_name.trim();
                if trimmed.is_empty() {
                    resolved_source.label().to_string()
                } else {
                    trimmed.to_string()
                }
            };

            ui.group(|ui| {
                ui.strong("Import Preflight");
                ui.small(format!(
                    "anchor: {}:{}:{}..{} (strand {})",
                    selected_anchor.genome_id,
                    selected_anchor.chromosome,
                    selected_anchor.start_1based,
                    selected_anchor.end_1based,
                    selected_anchor
                        .strand
                        .map(|s| s.to_string())
                        .unwrap_or_else(|| "?".to_string())
                ));
                ui.small(format!(
                    "matching status: {} / {} anchored sequences share this mapping group ({})",
                    matching_anchor_count, projected_targets, mapping_status
                ));
                ui.small(format!(
                    "projected tracks: {} target sequence(s) with track '{}'",
                    projected_targets, projected_track_name
                ));
                if path_exists {
                    ui.small(format!("track file: {}", track_path));
                } else {
                    ui.colored_label(
                        egui::Color32::from_rgb(190, 70, 70),
                        "track file: missing or unreadable path",
                    );
                }
                ui.horizontal(|ui| {
                    let toggle_resp = ui
                        .checkbox(
                            &mut self.genome_track_preflight_track_subscription,
                            "Track this file for auto-sync",
                        )
                        .on_hover_text(
                            "When enabled, this file is added to tracked subscriptions and auto-applied to newly anchored sequences.",
                        );
                    if toggle_resp.hovered() {
                        self.hover_status_name =
                            "Genome Tracks > Track Subscription Toggle".to_string();
                    }
                    let response = ui
                        .add_enabled(
                            !anchored_seq_ids.is_empty() && !import_running && path_exists,
                            egui::Button::new("Apply To All Anchored Now"),
                        )
                        .on_hover_text(
                            "One-click import to all anchored sequences using the current preflight settings.",
                        );
                    let response = self.track_hover_status(
                        response,
                        "Genome Tracks > Apply To All Anchored Now",
                    );
                    if response.clicked() {
                        self.import_genome_bed_track_for_all_anchored_sequences(
                            self.genome_track_preflight_track_subscription,
                        );
                    }
                });
            });
        }
        if self.genome_track_import_task.is_some() {
            let elapsed_secs = self
                .genome_track_import_task
                .as_ref()
                .map(|task| task.started.elapsed().as_secs_f32())
                .unwrap_or(0.0);
            let progress_snapshot = self.genome_track_import_progress.clone();
            ui.horizontal(|ui| {
                ui.add(egui::Spinner::new());
                if let Some(progress) = &progress_snapshot {
                    ui.label(format!(
                        "Running import: {} '{}' parsed={} imported={} skipped={} ({:.1}s)",
                        progress.source,
                        progress.path,
                        progress.parsed_records,
                        progress.imported_features,
                        progress.skipped_records,
                        elapsed_secs
                    ));
                } else {
                    ui.label(format!("Running import task... ({:.1}s)", elapsed_secs));
                }
                let cancel_import_resp = self.track_hover_status(
                    ui.button("Cancel Import").on_hover_text(
                        "Request cancellation. Imported features up to the cancellation point are kept.",
                    ),
                    "Genome Tracks > Cancel Import",
                );
                if cancel_import_resp.clicked() {
                    self.request_track_import_task_cancel("genome tracks dialog");
                }
            });
        }
        ui.horizontal(|ui| {
            let all_once_resp = self.track_hover_status(
                ui.add_enabled(
                    !anchored_seq_ids.is_empty() && !import_running,
                    egui::Button::new("Import To All Anchored (One-Time)"),
                )
                .on_hover_text(
                    "Import once onto every currently anchored sequence. No subscription is saved for future extracts.",
                ),
                "Genome Tracks > Import All Anchored One-Time",
            );
            if all_once_resp.clicked()
            {
                self.import_genome_bed_track_for_all_anchored_sequences(false);
            }
            let all_track_resp = self.track_hover_status(
                ui.add_enabled(
                    !anchored_seq_ids.is_empty() && !import_running,
                    egui::Button::new("Import To All Anchored + Track"),
                )
                .on_hover_text(
                    "Import onto all current anchored sequences and save this file as a tracked subscription for automatic future auto-sync.",
                ),
                "Genome Tracks > Import All Anchored And Track",
            );
            if all_track_resp.clicked()
            {
                self.import_genome_bed_track_for_all_anchored_sequences(true);
            }
        });

        ui.separator();
        ui.label("Tracked genome signal files for auto-sync to new anchored sequences");

        let mut apply_now_index: Option<usize> = None;
        let mut remove_index: Option<usize> = None;
        if self.genome_bed_track_subscriptions.is_empty() {
            ui.small("No tracked files yet.");
        } else {
            let subscription_rows = self.genome_bed_track_subscriptions.clone();
            egui::Grid::new("genome_bed_track_subscriptions_grid")
                .striped(true)
                .show(ui, |ui| {
                    ui.strong("Type");
                    ui.strong("Path");
                    ui.strong("Track");
                    ui.strong("Score Range");
                    ui.strong("Clear Existing");
                    ui.strong("Actions");
                    ui.end_row();
                    for (index, subscription) in subscription_rows.iter().enumerate() {
                        ui.label(subscription.source.label());
                        ui.monospace(subscription.path.as_str());
                        ui.label(
                            subscription
                                .track_name
                                .clone()
                                .unwrap_or_else(|| "-".to_string()),
                        );
                        ui.label(format!(
                            "{}..{}",
                            subscription
                                .min_score
                                .map(|v| format!("{v:.3}"))
                                .unwrap_or_else(|| "-".to_string()),
                            subscription
                                .max_score
                                .map(|v| format!("{v:.3}"))
                                .unwrap_or_else(|| "-".to_string())
                        ));
                        ui.label(if subscription.clear_existing {
                            "yes"
                        } else {
                            "no"
                        });
                        ui.horizontal(|ui| {
                            let apply_resp = self.track_hover_status(
                                ui.add_enabled(!import_running, egui::Button::new("Apply now"))
                                    .on_hover_text(
                                        "Re-apply this tracked file to all currently anchored sequences now.",
                                    ),
                                "Genome Tracks > Apply Tracked File",
                            );
                            if apply_resp.clicked()
                            {
                                apply_now_index = Some(index);
                            }
                            let remove_resp = self.track_hover_status(
                                ui.add_enabled(!import_running, egui::Button::new("Remove"))
                                    .on_hover_text(
                                        "Remove this tracked subscription (does not remove already imported features).",
                                    ),
                                "Genome Tracks > Remove Tracked File",
                            );
                            if remove_resp.clicked()
                            {
                                remove_index = Some(index);
                            }
                        });
                        ui.end_row();
                    }
                });
        }

        if let Some(index) = apply_now_index {
            self.apply_tracked_bed_subscription_to_all_anchored(index);
        }
        if let Some(index) = remove_index {
            let result = self
                .engine
                .write()
                .unwrap()
                .remove_genome_track_subscription(index);
            match result {
                Ok(removed) => {
                    self.load_bed_track_subscriptions_from_state();
                    self.genome_track_status = format!(
                        "Removed tracked {} '{}'",
                        removed.source.label(),
                        removed.path
                    );
                }
                Err(e) => {
                    self.genome_track_status =
                        format!("Could not remove tracked file: {}", e.message);
                }
            }
        }
        if !self.genome_bed_track_subscriptions.is_empty() {
            let clear_resp = self.track_hover_status(
                ui.add_enabled(!import_running, egui::Button::new("Clear Tracked Files"))
                    .on_hover_text(
                        "Remove all tracked subscriptions (does not remove already imported features).",
                    ),
                "Genome Tracks > Clear Tracked Files",
            );
            if clear_resp.clicked() {
                self.engine
                    .write()
                    .unwrap()
                    .clear_genome_track_subscriptions();
                self.load_bed_track_subscriptions_from_state();
                self.genome_track_status = "Cleared all tracked files".to_string();
            }
        }

        if !self.genome_track_autosync_status.is_empty() {
            ui.separator();
            ui.monospace(&self.genome_track_autosync_status);
        }
        if !self.genome_track_status.is_empty() {
            ui.separator();
            ui.monospace(&self.genome_track_status);
        }
    }

    fn render_genome_bed_track_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_genome_bed_track_dialog {
            return;
        }
        let builder = egui::ViewportBuilder::default()
            .with_title("Import Genome Tracks")
            .with_inner_size([980.0, 620.0])
            .with_min_inner_size([620.0, 320.0]);
        ctx.show_viewport_immediate(Self::bed_track_viewport_id(), builder, |ctx, class| {
            self.note_viewport_focus_if_active(ctx, Self::bed_track_viewport_id());
            if class == egui::ViewportClass::Embedded {
                let mut open = self.show_genome_bed_track_dialog;
                egui::Window::new("Import Genome Tracks")
                    .open(&mut open)
                    .collapsible(false)
                    .resizable(true)
                    .default_size(Vec2::new(980.0, 620.0))
                    .show(ctx, |ui| self.render_genome_bed_track_contents(ui));
                self.show_genome_bed_track_dialog = open;
                return;
            }

            egui::CentralPanel::default().show(ctx, |ui| {
                self.render_genome_bed_track_contents(ui);
            });

            if ctx.input(|i| i.viewport().close_requested()) {
                self.show_genome_bed_track_dialog = false;
            }
        });
    }

    fn render_agent_assistant_contents(&mut self, ui: &mut Ui) {
        self.refresh_agent_system_catalog();
        self.render_specialist_window_nav(ui);
        ui.label("Ask an agent system for project support, then execute suggested GENtle shell commands per reply.");
        ui.horizontal(|ui| {
            ui.label("catalog");
            ui.text_edit_singleline(&mut self.agent_catalog_path);
            if ui
                .button("Browse...")
                .on_hover_text("Browse filesystem and fill this path")
                .clicked()
            {
                if let Some(path) = rfd::FileDialog::new()
                    .add_filter("JSON", &["json"])
                    .pick_file()
                {
                    self.agent_catalog_path = path.display().to_string();
                    self.agent_catalog_loaded_path.clear();
                    self.refresh_agent_system_catalog();
                }
            }
        });
        if !self.agent_catalog_error.is_empty() {
            ui.colored_label(
                egui::Color32::from_rgb(190, 70, 70),
                format!("Catalog error: {}", self.agent_catalog_error),
            );
        }
        ui.horizontal(|ui| {
            ui.label("system");
            egui::ComboBox::from_id_salt("agent_system_combo")
                .selected_text(if self.agent_system_id.trim().is_empty() {
                    "(choose system)"
                } else {
                    self.agent_system_id.as_str()
                })
                .show_ui(ui, |ui| {
                    for system in &self.agent_systems {
                        let (available, reason) = self.selected_agent_system_availability(system);
                        let label = if available {
                            format!("{} ({})", system.label, system.id)
                        } else {
                            format!("{} ({}) [unavailable]", system.label, system.id)
                        };
                        let mut response = ui.add(
                            egui::Button::new(label).selected(self.agent_system_id == system.id),
                        );
                        if !available {
                            response = response.on_hover_text(
                                reason.unwrap_or_else(|| "agent system unavailable".to_string()),
                            );
                        }
                        if response.clicked() {
                            self.agent_system_id = system.id.clone();
                        }
                    }
                });
        });
        let mut selected_available = false;
        if let Some(system) = self.selected_agent_system() {
            let (available, reason) = self.selected_agent_system_availability(&system);
            selected_available = available;
            if let Some(description) = system.description.as_deref() {
                let trimmed = description.trim();
                if !trimmed.is_empty() {
                    ui.small(trimmed);
                }
            }
            ui.small(format!("transport: {}", system.transport.as_str()));
            if !system.command.is_empty() {
                ui.small(format!("command: {}", system.command.join(" ")));
            }
            if matches!(
                system.transport,
                AgentSystemTransport::NativeOpenai | AgentSystemTransport::NativeOpenaiCompat
            ) {
                if let Some(source_key) = self.selected_agent_model_discovery_source_key(&system) {
                    if self.agent_model_discovery_source_key != source_key {
                        self.agent_model_discovery_source_key = source_key;
                        self.agent_model_discovery_task = None;
                        self.agent_discovered_models.clear();
                        self.agent_discovered_model_pick.clear();
                        self.agent_model_discovery_status.clear();
                    }
                    if normalize_agent_model_name(self.agent_model_override.trim()).is_none()
                        && self.agent_model_discovery_task.is_none()
                        && self.agent_discovered_models.is_empty()
                    {
                        self.start_agent_model_discovery_task(&system, false);
                    }
                }
                let catalog_base_url = system
                    .base_url
                    .as_deref()
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .unwrap_or("(transport default)");
                if self.agent_base_url_override.trim().is_empty() {
                    ui.small(format!("base URL: {catalog_base_url}"));
                } else {
                    ui.small(format!(
                        "base URL: {} (session override)",
                        self.agent_base_url_override.trim()
                    ));
                }
                let catalog_model = system
                    .model
                    .as_deref()
                    .map(str::trim)
                    .and_then(normalize_agent_model_name)
                    .unwrap_or_else(|| OPENAI_COMPAT_UNSPECIFIED_MODEL.to_string());
                let model_override = normalize_agent_model_name(self.agent_model_override.trim());
                if let Some(model_override) = model_override {
                    ui.small(format!("model: {model_override} (session override)"));
                } else if matches!(system.transport, AgentSystemTransport::NativeOpenaiCompat)
                    && normalize_agent_model_name(&self.agent_discovered_model_pick).is_some()
                {
                    ui.small(format!(
                        "model: {} (selected discovered model)",
                        self.agent_discovered_model_pick.trim()
                    ));
                } else {
                    ui.small(format!("model: {catalog_model}"));
                }
            } else {
                self.agent_model_discovery_task = None;
                self.agent_discovered_models.clear();
                self.agent_discovered_model_pick.clear();
                self.agent_model_discovery_source_key.clear();
                self.agent_model_discovery_status.clear();
            }
            if !available {
                ui.colored_label(
                    egui::Color32::from_rgb(190, 70, 70),
                    format!(
                        "Unavailable: {}",
                        reason.unwrap_or_else(|| "unknown reason".to_string())
                    ),
                );
            }
        } else if self.agent_systems.is_empty() {
            ui.small("No systems loaded from this catalog.");
        }
        ui.horizontal(|ui| {
            ui.label("OpenAI API key");
            ui.add(
                egui::TextEdit::singleline(&mut self.agent_openai_api_key)
                    .password(true)
                    .hint_text("sk-..."),
            );
            if ui
                .button("Clear Key")
                .on_hover_text("Clear session-only API key override")
                .clicked()
            {
                self.agent_openai_api_key.clear();
            }
        });
        ui.horizontal(|ui| {
            ui.label("Base URL override");
            ui.add(
                egui::TextEdit::singleline(&mut self.agent_base_url_override)
                    .hint_text("http://localhost:11964/v1"),
            );
            if ui
                .button("Clear URL")
                .on_hover_text("Clear session-only base URL override")
                .clicked()
            {
                self.agent_base_url_override.clear();
            }
        });
        ui.horizontal(|ui| {
            ui.label("Model override");
            ui.add(
                egui::TextEdit::singleline(&mut self.agent_model_override).hint_text("unspecified"),
            );
            if ui
                .button("Clear Model")
                .on_hover_text("Clear session-only model override")
                .clicked()
            {
                self.agent_model_override.clear();
            }
        });
        ui.horizontal(|ui| {
            ui.label("timeout_sec");
            ui.add(
                egui::TextEdit::singleline(&mut self.agent_timeout_secs)
                    .desired_width(100.0)
                    .hint_text("default"),
            );
            if ui
                .button("Clear Timeout")
                .on_hover_text("Use default timeout")
                .clicked()
            {
                self.agent_timeout_secs.clear();
            }
        });
        ui.horizontal(|ui| {
            ui.label("connect_timeout_sec");
            ui.add(
                egui::TextEdit::singleline(&mut self.agent_connect_timeout_secs)
                    .desired_width(90.0)
                    .hint_text("default"),
            );
            ui.label("read_timeout_sec");
            ui.add(
                egui::TextEdit::singleline(&mut self.agent_read_timeout_secs)
                    .desired_width(90.0)
                    .hint_text("default"),
            );
            if ui
                .button("Clear HTTP Timeouts")
                .on_hover_text("Use default connect/read timeouts")
                .clicked()
            {
                self.agent_connect_timeout_secs.clear();
                self.agent_read_timeout_secs.clear();
            }
        });
        ui.horizontal(|ui| {
            ui.label("max_retries");
            ui.add(
                egui::TextEdit::singleline(&mut self.agent_max_retries)
                    .desired_width(90.0)
                    .hint_text("default"),
            );
            ui.label("max_response_bytes");
            ui.add(
                egui::TextEdit::singleline(&mut self.agent_max_response_bytes)
                    .desired_width(120.0)
                    .hint_text("default"),
            );
            if ui
                .button("Clear Limits")
                .on_hover_text("Use default retry/response-size limits")
                .clicked()
            {
                self.agent_max_retries.clear();
                self.agent_max_response_bytes.clear();
            }
        });
        if let Some(system) = self.selected_agent_system() {
            if matches!(
                system.transport,
                AgentSystemTransport::NativeOpenai | AgentSystemTransport::NativeOpenaiCompat
            ) {
                ui.horizontal(|ui| {
                    if ui
                        .button("Discover Models")
                        .on_hover_text("Query local/server model list from current base URL")
                        .clicked()
                    {
                        self.start_agent_model_discovery_task(&system, true);
                    }
                    if let Some(task) = &self.agent_model_discovery_task {
                        ui.add(egui::Spinner::new());
                        ui.small(format!(
                            "Discovering models ({:.1}s)",
                            task.started.elapsed().as_secs_f32()
                        ));
                    }
                });
                if !self.agent_discovered_models.is_empty() {
                    ui.horizontal(|ui| {
                        ui.label("Discovered model");
                        egui::ComboBox::from_id_salt("agent_discovered_model_combo")
                            .selected_text(if self.agent_discovered_model_pick.trim().is_empty() {
                                "(choose model)"
                            } else {
                                self.agent_discovered_model_pick.as_str()
                            })
                            .show_ui(ui, |ui| {
                                for model in &self.agent_discovered_models {
                                    ui.selectable_value(
                                        &mut self.agent_discovered_model_pick,
                                        model.clone(),
                                        model,
                                    );
                                }
                            });
                    });
                    ui.small(
                        "If Model override is unspecified, the selected discovered model is used.",
                    );
                }
                if !self.agent_model_discovery_status.trim().is_empty() {
                    ui.small(self.agent_model_discovery_status.clone());
                }
            }
        }
        ui.small(
            "Session only: if set, this key overrides OPENAI_API_KEY for agent requests started from this GUI window.",
        );
        ui.small(
            "Session only: Base URL override applies to native_openai/native_openai_compat. For local roots (e.g. http://localhost:11964), GENtle tries /chat/completions and /v1/chat/completions on that same base URL.",
        );
        ui.small(
            "Session only: Model override applies to native_openai/native_openai_compat and maps to GENTLE_AGENT_MODEL. Value 'unspecified' means no override.",
        );
        ui.small(
            "Session only: timeout_sec maps to GENTLE_AGENT_TIMEOUT_SECS and applies to agent requests (stdio and native transports).",
        );
        ui.small(
            "Session only: connect_timeout_sec/read_timeout_sec map to GENTLE_AGENT_CONNECT_TIMEOUT_SECS/GENTLE_AGENT_READ_TIMEOUT_SECS.",
        );
        ui.small(
            "Session only: max_retries/max_response_bytes map to GENTLE_AGENT_MAX_RETRIES/GENTLE_AGENT_MAX_RESPONSE_BYTES.",
        );
        ui.horizontal(|ui| {
            ui.checkbox(
                &mut self.agent_include_state_summary,
                "Include project state summary in request",
            );
            ui.checkbox(
                &mut self.agent_allow_auto_exec,
                "Auto-run suggestions marked as 'auto'",
            );
        });
        if !agent_prompt_template_options()
            .iter()
            .any(|(id, _)| *id == self.agent_prompt_template_id)
        {
            self.agent_prompt_template_id = AGENT_PROMPT_TEMPLATE_DEFAULT_ID.to_string();
        }
        ui.horizontal(|ui| {
            ui.label("Prompt template");
            egui::ComboBox::from_id_salt("agent_prompt_template_combo")
                .selected_text(agent_prompt_template_label(&self.agent_prompt_template_id))
                .show_ui(ui, |ui| {
                    for (id, label) in agent_prompt_template_options() {
                        ui.selectable_value(
                            &mut self.agent_prompt_template_id,
                            (*id).to_string(),
                            *label,
                        );
                    }
                });
            if ui
                .button("Insert")
                .on_hover_text("Replace current prompt with selected template")
                .clicked()
            {
                self.agent_prompt =
                    agent_prompt_template_text(&self.agent_prompt_template_id).to_string();
            }
            if ui
                .button("Append")
                .on_hover_text("Append selected template below current prompt")
                .clicked()
            {
                let template_text = agent_prompt_template_text(&self.agent_prompt_template_id);
                if self.agent_prompt.trim().is_empty() {
                    self.agent_prompt = template_text.to_string();
                } else {
                    if !self.agent_prompt.ends_with('\n') {
                        self.agent_prompt.push('\n');
                    }
                    self.agent_prompt.push('\n');
                    self.agent_prompt.push_str(template_text);
                }
            }
        });
        ui.label("Prompt");
        ui.add(
            egui::TextEdit::multiline(&mut self.agent_prompt)
                .desired_rows(6)
                .desired_width(f32::INFINITY),
        );
        let running = self.agent_task.is_some();
        ui.horizontal(|ui| {
            if ui
                .add_enabled(
                    !running && selected_available,
                    egui::Button::new("Ask Agent"),
                )
                .on_hover_text("Send prompt to selected agent system")
                .clicked()
            {
                self.start_agent_assistant_request();
            }
            if ui
                .button("Clear Response")
                .on_hover_text("Clear latest agent response and status")
                .clicked()
            {
                self.agent_last_invocation = None;
                self.agent_status.clear();
            }
            if ui
                .button("Clear Execution Log")
                .on_hover_text("Clear local execution history for agent suggestions")
                .clicked()
            {
                self.agent_execution_log.clear();
            }
            if ui
                .button("Close")
                .on_hover_text("Close this dialog")
                .clicked()
            {
                self.show_agent_assistant_dialog = false;
            }
        });
        if let Some(task) = &self.agent_task {
            ui.horizontal(|ui| {
                ui.add(egui::Spinner::new());
                ui.label(format!(
                    "Agent request running ({:.1}s)",
                    task.started.elapsed().as_secs_f32()
                ));
            });
        }
        if !self.agent_status.is_empty() {
            ui.separator();
            ui.monospace(self.agent_status.clone());
        }

        if let Some(invocation) = self.agent_last_invocation.clone() {
            ui.separator();
            ui.label(format!(
                "Latest response from {} ({})",
                invocation.system_label, invocation.system_id
            ));
            ui.small(format!(
                "elapsed={} ms | transport={} | exit_code={:?}",
                invocation.elapsed_ms, invocation.transport, invocation.exit_code
            ));
            ui.small(format!(
                "runtime: timeout={}s | connect={:?}s | read={:?}s | max_retries={} | max_response_bytes={}",
                invocation.runtime.timeout_secs,
                invocation.runtime.connect_timeout_secs,
                invocation.runtime.read_timeout_secs,
                invocation.runtime.max_retries,
                invocation.runtime.max_response_bytes
            ));
            if !invocation.runtime.endpoint_candidates.is_empty() {
                ui.small(format!(
                    "endpoint candidates: {}",
                    invocation.runtime.endpoint_candidates.join(" | ")
                ));
            }
            if !invocation.runtime.attempted_endpoints.is_empty() {
                ui.small(format!(
                    "attempted endpoints: {}",
                    invocation.runtime.attempted_endpoints.join(" | ")
                ));
            }
            if let Some(selected_endpoint) = invocation.runtime.selected_endpoint.as_deref() {
                ui.small(format!("selected endpoint: {}", selected_endpoint));
            }
            if !invocation.response.assistant_message.trim().is_empty() {
                ui.group(|ui| {
                    ui.strong("Agent message");
                    ui.label(invocation.response.assistant_message.trim());
                });
            }
            if !invocation.response.questions.is_empty() {
                ui.group(|ui| {
                    ui.strong("Agent questions");
                    for question in &invocation.response.questions {
                        ui.label(format!("- {}", question));
                    }
                });
            }
            if invocation.response.suggested_commands.is_empty() {
                ui.small("No executable suggestions in this reply.");
            } else {
                ui.separator();
                ui.strong("Suggested commands");
                let mut run_request: Option<(usize, String)> = None;
                egui::Grid::new("agent_suggested_commands_grid")
                    .striped(true)
                    .show(ui, |ui| {
                        ui.strong("#");
                        ui.strong("Intent");
                        ui.strong("Command");
                        ui.strong("Rationale");
                        ui.strong("Action");
                        ui.end_row();
                        for (idx, suggestion) in invocation
                            .response
                            .suggested_commands
                            .iter()
                            .enumerate()
                        {
                            let index_1based = idx + 1;
                            ui.label(index_1based.to_string());
                            ui.label(suggestion.execution.as_str());
                            ui.monospace(suggestion.command.as_str());
                            ui.label(suggestion.rationale.clone().unwrap_or_default());
                            let can_run = !suggestion.command.trim().is_empty()
                                && suggestion.execution != AgentExecutionIntent::Chat;
                            let run_resp = ui
                                .add_enabled(can_run, egui::Button::new("Run"))
                                .on_hover_text(
                                    "Execute this suggested command using GENtle shared shell parser/executor",
                                );
                            if run_resp.clicked() {
                                run_request = Some((index_1based, suggestion.command.clone()));
                            }
                            ui.end_row();
                        }
                    });
                if let Some((index_1based, command)) = run_request {
                    self.execute_agent_suggested_command(index_1based, &command, "manual");
                }
            }
            if !invocation.raw_stderr.trim().is_empty() {
                ui.separator();
                ui.strong("Agent stderr");
                let mut stderr = invocation.raw_stderr.clone();
                ui.add(
                    egui::TextEdit::multiline(&mut stderr)
                        .desired_rows(4)
                        .desired_width(f32::INFINITY),
                );
            }
        }

        if !self.agent_execution_log.is_empty() {
            ui.separator();
            ui.strong("Execution log");
            egui::ScrollArea::vertical()
                .max_height(180.0)
                .show(ui, |ui| {
                    for entry in self.agent_execution_log.iter().rev() {
                        ui.label(format!(
                            "#{} [{}] {} | {} | changed={} | t={}",
                            entry.index_1based,
                            entry.trigger,
                            if entry.ok { "ok" } else { "error" },
                            entry.command,
                            entry.state_changed,
                            entry.executed_at_unix_ms
                        ));
                        ui.small(entry.summary.clone());
                    }
                });
        }
    }

    fn render_agent_assistant_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_agent_assistant_dialog {
            return;
        }
        let builder = egui::ViewportBuilder::default()
            .with_title("Agent Assistant")
            .with_inner_size([980.0, 720.0])
            .with_min_inner_size([640.0, 420.0]);
        ctx.show_viewport_immediate(
            Self::agent_assistant_viewport_id(),
            builder,
            |ctx, class| {
                self.note_viewport_focus_if_active(ctx, Self::agent_assistant_viewport_id());
                if class == egui::ViewportClass::Embedded {
                    let mut open = self.show_agent_assistant_dialog;
                    egui::Window::new("Agent Assistant")
                        .open(&mut open)
                        .collapsible(false)
                        .resizable(true)
                        .default_size(Vec2::new(980.0, 720.0))
                        .show(ctx, |ui| self.render_agent_assistant_contents(ui));
                    self.show_agent_assistant_dialog = open;
                    return;
                }

                egui::CentralPanel::default().show(ctx, |ui| {
                    self.render_agent_assistant_contents(ui);
                });

                if ctx.input(|i| i.viewport().close_requested()) {
                    self.show_agent_assistant_dialog = false;
                }
            },
        );
    }

    fn refresh_project_restriction_enzymes(&mut self, resource_path: &str) -> Result<usize> {
        let enzymes = enzymes::load_restriction_enzymes_from_path(resource_path)?;
        let mut engine = self.engine.write().unwrap();
        for dna in engine.state_mut().sequences.values_mut() {
            enzymes.clone_into(dna.restriction_enzymes_mut());
            dna.update_computed_features();
        }
        Ok(enzymes.len())
    }

    fn prompt_import_rebase_resource(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("REBASE withrefm", &["txt", "withrefm", "f19", "f20", "f21"])
            .pick_file()
        {
            let input = path.display().to_string();
            match resource_sync::sync_rebase(
                &input,
                Some(resource_sync::DEFAULT_REBASE_RESOURCE_PATH),
                false,
            ) {
                Ok(report) => {
                    let seq_count = self.engine.read().unwrap().state().sequences.len();
                    let refresh_status = self.refresh_project_restriction_enzymes(&report.output);
                    match refresh_status {
                        Ok(loaded_count) => println!(
                            "Imported REBASE from '{}' ({} enzymes); active set now {} enzymes, refreshed {} sequence(s)",
                            report.source, report.item_count, loaded_count, seq_count
                        ),
                        Err(e) => eprintln!(
                            "Imported REBASE from '{}', but could not refresh loaded sequences: {}",
                            report.source, e
                        ),
                    }
                }
                Err(e) => eprintln!("Could not import REBASE from '{}': {}", input, e),
            }
        }
    }

    fn prompt_import_jaspar_resource(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("JASPAR motif matrix", &["txt", "jaspar", "pfm"])
            .pick_file()
        {
            let input = path.display().to_string();
            match resource_sync::sync_jaspar(
                &input,
                Some(resource_sync::DEFAULT_JASPAR_RESOURCE_PATH),
            ) {
                Ok(report) => {
                    tf_motifs::reload();
                    println!(
                        "Imported JASPAR from '{}' ({} motifs) to '{}'",
                        report.source, report.item_count, report.output
                    );
                }
                Err(e) => eprintln!("Could not import JASPAR from '{}': {}", input, e),
            }
        }
    }

    fn save_current_project(&mut self) -> bool {
        if let Some(path) = self.current_project_path.clone() {
            if self.save_project_to_file(&path).is_ok() {
                self.record_recent_project_path(&path);
                self.mark_clean_snapshot();
                return true;
            }
            return false;
        }

        if let Some(path) = rfd::FileDialog::new()
            .set_file_name("project.gentle.json")
            .add_filter("GENtle project", &["json"])
            .save_file()
        {
            let path = path.display().to_string();
            if self.save_project_to_file(&path).is_ok() {
                self.set_current_project_path_and_track_recent(&path);
                self.mark_clean_snapshot();
                return true;
            }
        }
        false
    }

    fn execute_project_action(&mut self, action: ProjectAction) {
        match action {
            ProjectAction::New => self.new_project(),
            ProjectAction::Open => self.prompt_open_project(),
            ProjectAction::OpenPath(path) => self.open_project_path(&path),
            ProjectAction::Close => self.close_project(),
        }
    }

    fn request_project_action(&mut self, action: ProjectAction) {
        if self.project_has_user_content() && self.is_project_dirty() {
            self.pending_project_action = Some(action);
        } else {
            self.execute_project_action(action);
        }
    }

    fn new_project(&mut self) {
        self.reset_to_empty_project();
    }

    fn close_project(&mut self) {
        self.reset_to_empty_project();
    }

    fn load_project_from_file(&mut self, path: &str) -> Result<()> {
        let state = ProjectState::load_from_path(path).map_err(|e| anyhow!(e.to_string()))?;

        self.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
        self.set_current_project_path_and_track_recent(path);
        self.lineage_cache_valid = false;
        self.lineage_rows.clear();
        self.lineage_edges.clear();
        self.lineage_op_label_by_id.clear();
        self.lineage_containers.clear();
        self.lineage_arrangements.clear();
        self.lineage_graph_zoom = 1.0;
        self.lineage_graph_area_height = 420.0;
        self.lineage_container_area_height = 220.0;
        self.lineage_graph_scroll_offset = Vec2::ZERO;
        self.lineage_graph_pan_origin = None;
        self.lineage_graph_compact_labels = true;
        self.lineage_graph_selected_node_id = None;
        self.lineage_graph_node_offsets.clear();
        self.lineage_graph_drag_origin = None;
        self.lineage_graph_offsets_synced_stamp = 0;
        self.tracked_autosync_last_op_count = None;
        self.new_windows.clear();
        self.windows.clear();
        self.windows_to_close.write().unwrap().clear();
        self.pending_focus_viewports.clear();
        self.agent_task = None;
        self.agent_model_discovery_task = None;
        self.agent_status.clear();
        self.agent_last_invocation = None;
        self.agent_execution_log.clear();
        self.agent_discovered_models.clear();
        self.agent_discovered_model_pick.clear();
        self.agent_model_discovery_status.clear();
        self.agent_model_discovery_source_key.clear();
        self.load_bed_track_subscriptions_from_state();
        self.load_lineage_graph_workspace_from_state();

        self.mark_clean_snapshot();
        Ok(())
    }

    fn current_project_name(&self) -> String {
        match &self.current_project_path {
            Some(path) => Path::new(path)
                .file_name()
                .map(|s| s.to_string_lossy().to_string())
                .unwrap_or_else(|| path.clone()),
            None => "Untitled Project".to_string(),
        }
    }

    fn load_lineage_graph_workspace_from_state(&mut self) {
        let (workspace_serialized, legacy_offsets_serialized) = {
            let engine = self.engine.read().unwrap();
            let metadata = &engine.state().metadata;
            (
                metadata.get(LINEAGE_GRAPH_WORKSPACE_METADATA_KEY).cloned(),
                metadata.get(LINEAGE_NODE_OFFSETS_METADATA_KEY).cloned(),
            )
        };

        self.lineage_graph_zoom = 1.0;
        self.lineage_graph_area_height = 420.0;
        self.lineage_container_area_height = 220.0;
        self.lineage_graph_scroll_offset = Vec2::ZERO;
        self.lineage_graph_pan_origin = None;
        self.lineage_graph_compact_labels = true;
        self.lineage_graph_selected_node_id = None;
        self.lineage_graph_node_offsets.clear();
        if let Some(serialized) = workspace_serialized {
            if let Ok(workspace) =
                serde_json::from_value::<PersistedLineageGraphWorkspace>(serialized)
            {
                if workspace.zoom.is_finite() {
                    self.lineage_graph_zoom = workspace.zoom.clamp(0.35, 4.0);
                }
                if workspace.graph_area_height.is_finite() {
                    self.lineage_graph_area_height =
                        workspace.graph_area_height.clamp(220.0, 2400.0);
                }
                if workspace.container_area_height.is_finite() {
                    self.lineage_container_area_height =
                        workspace.container_area_height.clamp(120.0, 1600.0);
                }
                if workspace.scroll_offset[0].is_finite() && workspace.scroll_offset[1].is_finite()
                {
                    self.lineage_graph_scroll_offset = Vec2::new(
                        workspace.scroll_offset[0].max(0.0),
                        workspace.scroll_offset[1].max(0.0),
                    );
                }
                self.lineage_graph_compact_labels = workspace.compact_labels;
                for (node_id, pair) in workspace.node_offsets {
                    if pair[0].is_finite() && pair[1].is_finite() {
                        self.lineage_graph_node_offsets
                            .insert(node_id, Vec2::new(pair[0], pair[1]));
                    }
                }
            }
        } else if let Some(serialized) = legacy_offsets_serialized {
            if let Ok(raw) = serde_json::from_value::<HashMap<String, [f32; 2]>>(serialized) {
                for (node_id, pair) in raw {
                    if pair[0].is_finite() && pair[1].is_finite() {
                        self.lineage_graph_node_offsets
                            .insert(node_id, Vec2::new(pair[0], pair[1]));
                    }
                }
            }
        }
        self.lineage_graph_drag_origin = None;
        self.lineage_graph_offsets_synced_stamp = 0;
    }

    fn persist_lineage_graph_workspace_to_state(&mut self) {
        let mut raw: HashMap<String, [f32; 2]> = HashMap::new();
        for (node_id, offset) in &self.lineage_graph_node_offsets {
            if offset.x.is_finite() && offset.y.is_finite() {
                raw.insert(node_id.clone(), [offset.x, offset.y]);
            }
        }

        let workspace = PersistedLineageGraphWorkspace {
            zoom: self.lineage_graph_zoom.clamp(0.35, 4.0),
            graph_area_height: self.lineage_graph_area_height.clamp(220.0, 2400.0),
            container_area_height: self.lineage_container_area_height.clamp(120.0, 1600.0),
            scroll_offset: [
                self.lineage_graph_scroll_offset.x.max(0.0),
                self.lineage_graph_scroll_offset.y.max(0.0),
            ],
            compact_labels: self.lineage_graph_compact_labels,
            node_offsets: raw.clone(),
        };

        let workspace_is_default = (workspace.zoom - 1.0).abs() <= 0.0001
            && (workspace.graph_area_height - 420.0).abs() <= 0.0001
            && (workspace.container_area_height - 220.0).abs() <= 0.0001
            && workspace.scroll_offset[0].abs() <= 0.0001
            && workspace.scroll_offset[1].abs() <= 0.0001
            && workspace.compact_labels
            && workspace.node_offsets.is_empty();

        let mut engine = self.engine.write().unwrap();
        let state = engine.state_mut();
        if workspace_is_default {
            state.metadata.remove(LINEAGE_GRAPH_WORKSPACE_METADATA_KEY);
        } else if let Ok(value) = serde_json::to_value(&workspace) {
            if state.metadata.get(LINEAGE_GRAPH_WORKSPACE_METADATA_KEY) != Some(&value) {
                state
                    .metadata
                    .insert(LINEAGE_GRAPH_WORKSPACE_METADATA_KEY.to_string(), value);
            }
        }

        if raw.is_empty() {
            state.metadata.remove(LINEAGE_NODE_OFFSETS_METADATA_KEY);
        } else if let Ok(value) = serde_json::to_value(raw) {
            if state.metadata.get(LINEAGE_NODE_OFFSETS_METADATA_KEY) != Some(&value) {
                state
                    .metadata
                    .insert(LINEAGE_NODE_OFFSETS_METADATA_KEY.to_string(), value);
            }
        }
    }

    fn collect_open_window_entries(&self) -> Vec<OpenWindowEntry> {
        let mut entries = vec![OpenWindowEntry {
            native_menu_key: Self::native_menu_key_for_viewport(ViewportId::ROOT),
            viewport_id: ViewportId::ROOT,
            title: format!("Main Window — {}", self.current_project_name()),
            detail: "Project workspace".to_string(),
        }];

        if self.show_configuration_dialog {
            entries.push(OpenWindowEntry {
                native_menu_key: Self::native_menu_key_for_viewport(
                    Self::configuration_viewport_id(),
                ),
                viewport_id: Self::configuration_viewport_id(),
                title: "Configuration".to_string(),
                detail: "External tools and graphics defaults".to_string(),
            });
        }
        if self.show_help_dialog {
            entries.push(OpenWindowEntry {
                native_menu_key: Self::native_menu_key_for_viewport(Self::help_viewport_id()),
                viewport_id: Self::help_viewport_id(),
                title: format!("Help — {}", self.active_help_title()),
                detail: "GUI/CLI manual".to_string(),
            });
        }
        if self.show_command_palette_dialog {
            entries.push(OpenWindowEntry {
                native_menu_key: Self::native_menu_key_for_viewport(
                    Self::command_palette_viewport_id(),
                ),
                viewport_id: Self::command_palette_viewport_id(),
                title: "Command Palette".to_string(),
                detail: "Action launcher".to_string(),
            });
        }
        if self.show_history_panel {
            entries.push(OpenWindowEntry {
                native_menu_key: Self::native_menu_key_for_viewport(Self::history_viewport_id()),
                viewport_id: Self::history_viewport_id(),
                title: "Operation History".to_string(),
                detail: "Undo/redo operation log".to_string(),
            });
        }
        if self.show_reference_genome_prepare_dialog {
            entries.push(OpenWindowEntry {
                native_menu_key: Self::native_menu_key_for_viewport(
                    Self::prepare_genome_viewport_id(),
                ),
                viewport_id: Self::prepare_genome_viewport_id(),
                title: "Prepare Reference Genome".to_string(),
                detail: "Reference/helper genome preparation".to_string(),
            });
        }
        if self.show_reference_genome_blast_dialog {
            entries.push(OpenWindowEntry {
                native_menu_key: Self::native_menu_key_for_viewport(
                    Self::blast_genome_viewport_id(),
                ),
                viewport_id: Self::blast_genome_viewport_id(),
                title: "BLAST Genome Sequence".to_string(),
                detail: "BLAST against prepared genomes/helpers".to_string(),
            });
        }
        if self.show_genome_bed_track_dialog {
            entries.push(OpenWindowEntry {
                native_menu_key: Self::native_menu_key_for_viewport(Self::bed_track_viewport_id()),
                viewport_id: Self::bed_track_viewport_id(),
                title: "Import Genome Tracks".to_string(),
                detail: "Import BED/BigWig/VCF track overlays".to_string(),
            });
        }
        if self.show_agent_assistant_dialog {
            entries.push(OpenWindowEntry {
                native_menu_key: Self::native_menu_key_for_viewport(
                    Self::agent_assistant_viewport_id(),
                ),
                viewport_id: Self::agent_assistant_viewport_id(),
                title: "Agent Assistant".to_string(),
                detail: "Agent chat and per-reply command execution".to_string(),
            });
        }

        let mut sequence_windows = self
            .windows
            .iter()
            .map(|(viewport_id, window)| {
                let title = window
                    .read()
                    .map(|w| w.name())
                    .unwrap_or_else(|_| "Sequence window".to_string());
                OpenWindowEntry {
                    native_menu_key: Self::native_menu_key_for_viewport(*viewport_id),
                    viewport_id: *viewport_id,
                    title,
                    detail: "Sequence map window".to_string(),
                }
            })
            .collect::<Vec<_>>();
        sequence_windows.sort_by(|left, right| left.title.cmp(&right.title));
        entries.extend(sequence_windows);
        entries
    }

    fn focus_window_viewport(&mut self, ctx: &egui::Context, viewport_id: ViewportId) {
        if viewport_id == Self::configuration_viewport_id() {
            self.show_configuration_dialog = true;
        } else if viewport_id == Self::help_viewport_id() {
            self.show_help_dialog = true;
        } else if viewport_id == Self::command_palette_viewport_id() {
            self.show_command_palette_dialog = true;
            self.command_palette_focus_query = true;
        } else if viewport_id == Self::history_viewport_id() {
            self.show_history_panel = true;
        } else if viewport_id == Self::prepare_genome_viewport_id() {
            self.show_reference_genome_prepare_dialog = true;
        } else if viewport_id == Self::blast_genome_viewport_id() {
            self.show_reference_genome_blast_dialog = true;
        } else if viewport_id == Self::bed_track_viewport_id() {
            self.show_genome_bed_track_dialog = true;
        } else if viewport_id == Self::agent_assistant_viewport_id() {
            self.show_agent_assistant_dialog = true;
        }

        ctx.send_viewport_cmd_to(viewport_id, egui::ViewportCommand::Visible(true));
        ctx.send_viewport_cmd_to(viewport_id, egui::ViewportCommand::Focus);
        self.set_active_window_viewport(viewport_id);
    }

    fn render_specialist_window_nav(&mut self, ui: &mut Ui) {
        ui.horizontal(|ui| {
            if ui
                .button("Main")
                .on_hover_text("Bring the main project window to front")
                .clicked()
            {
                self.queue_focus_viewport(ViewportId::ROOT);
            }
        });
        ui.separator();
    }

    pub fn render_menu_bar(&mut self, ui: &mut Ui) {
        let open_window_entries = self.collect_open_window_entries();
        self.native_window_key_to_viewport = open_window_entries
            .iter()
            .map(|entry| (entry.native_menu_key, entry.viewport_id))
            .collect();
        let native_window_entries = open_window_entries
            .iter()
            .map(|entry| (entry.native_menu_key, entry.title.clone()))
            .collect::<Vec<_>>();
        let active_window_key = self
            .active_window_menu_key
            .or_else(|| Some(Self::native_menu_key_for_viewport(ViewportId::ROOT)));
        if self.last_native_window_entries != native_window_entries
            || self.last_native_active_window_key != active_window_key
        {
            about::sync_native_open_windows_menu(&native_window_entries, active_window_key);
            self.last_native_window_entries = native_window_entries;
            self.last_native_active_window_key = active_window_key;
        }
        let (undo_count, redo_count) = {
            let engine = self.engine.read().expect("Engine lock poisoned");
            (engine.undo_available(), engine.redo_available())
        };
        let history_ops_enabled = !self.has_active_background_jobs();
        menu::bar(ui, |ui| {
            ui.menu_button(TRANSLATIONS.get("m_file"), |ui| {
                if ui
                    .button("New Project")
                    .on_hover_text("Create a new empty project state")
                    .clicked()
                {
                    self.request_project_action(ProjectAction::New);
                    ui.close_menu();
                }
                if ui
                    .button("Open Project...")
                    .on_hover_text("Open an existing .gentle.json project")
                    .clicked()
                {
                    self.request_project_action(ProjectAction::Open);
                    ui.close_menu();
                }
                ui.menu_button("Open Recent Project...", |ui| {
                    if self.recent_project_paths.is_empty() {
                        ui.add_enabled(false, egui::Button::new("No recent projects"));
                        return;
                    }

                    let recent_paths = self.recent_project_paths.clone();
                    let mut selected_recent_path: Option<String> = None;
                    for path in recent_paths {
                        let exists = Path::new(&path).is_file();
                        let mut label = Self::recent_project_menu_label(&path);
                        if !exists {
                            label.push_str(" [missing]");
                        }
                        let clicked = ui
                            .add_enabled(exists, egui::Button::new(label))
                            .on_hover_text(path.clone())
                            .clicked();
                        if clicked {
                            selected_recent_path = Some(path);
                        }
                    }

                    ui.separator();
                    if ui
                        .button("Clear Recent Projects")
                        .on_hover_text("Remove all items from Open Recent Project menu")
                        .clicked()
                    {
                        self.clear_recent_project_paths();
                        ui.close_menu();
                        return;
                    }

                    if let Some(path) = selected_recent_path {
                        self.request_project_action(ProjectAction::OpenPath(path));
                        ui.close_menu();
                    }
                });
                if ui
                    .button("Close Project")
                    .on_hover_text("Close the current project from the workspace")
                    .clicked()
                {
                    self.request_project_action(ProjectAction::Close);
                    ui.close_menu();
                }
                ui.separator();
                if ui
                    .button("Open Sequence...")
                    .on_hover_text("Import a sequence file into the current project")
                    .clicked()
                {
                    self.prompt_open_sequence();
                    ui.close_menu();
                }
                ui.separator();
                if ui
                    .button("Configuration...")
                    .on_hover_text("Open global app configuration and graphics defaults")
                    .clicked()
                {
                    self.open_configuration_dialog();
                    ui.close_menu();
                }
                ui.separator();
                if ui
                    .button("Prepare Reference Genome...")
                    .on_hover_text(
                        "Download/index a reference genome for local extraction and BLAST",
                    )
                    .clicked()
                {
                    self.open_reference_genome_prepare_dialog();
                    ui.close_menu();
                }
                if ui
                    .button("Prepared References...")
                    .on_hover_text("Inspect prepared reference/helper genome installations")
                    .clicked()
                {
                    self.open_reference_genome_inspector_dialog();
                    ui.close_menu();
                }
                if ui
                    .button("Retrieve Genomic Sequence...")
                    .on_hover_text(
                        "Extract anchored region/gene sequence from a prepared reference",
                    )
                    .clicked()
                {
                    self.open_reference_genome_retrieve_dialog();
                    ui.close_menu();
                }
                if ui
                    .button("BLAST Genome Sequence...")
                    .on_hover_text("Run BLAST against prepared reference genome indices")
                    .clicked()
                {
                    self.open_reference_genome_blast_dialog();
                    ui.close_menu();
                }
                if ui
                    .button("Import Genome Track...")
                    .on_hover_text("Import BED/BigWig/VCF tracks onto anchored sequences")
                    .clicked()
                {
                    self.open_genome_bed_track_dialog();
                    ui.close_menu();
                }
                if ui
                    .button("Agent Assistant...")
                    .on_hover_text(
                        "Ask configured agent systems and execute suggested shared-shell commands",
                    )
                    .clicked()
                {
                    self.open_agent_assistant_dialog();
                    ui.close_menu();
                }
                if ui
                    .button("Prepare Helper Genome...")
                    .on_hover_text("Prepare helper catalog genomes for extraction and BLAST")
                    .clicked()
                {
                    self.open_helper_genome_prepare_dialog();
                    ui.close_menu();
                }
                if ui
                    .button("Retrieve Helper Sequence...")
                    .on_hover_text("Extract sequence from prepared helper genomes")
                    .clicked()
                {
                    self.open_helper_genome_retrieve_dialog();
                    ui.close_menu();
                }
                if ui
                    .button("BLAST Helper Sequence...")
                    .on_hover_text("Run BLAST against prepared helper genome indices")
                    .clicked()
                {
                    self.open_helper_genome_blast_dialog();
                    ui.close_menu();
                }
                if ui
                    .button("Import REBASE Data...")
                    .on_hover_text("Import restriction-enzyme definitions from REBASE JSON/TXT")
                    .clicked()
                {
                    self.prompt_import_rebase_resource();
                    ui.close_menu();
                }
                if ui
                    .button("Import JASPAR Data...")
                    .on_hover_text("Import TF motif library from JASPAR JSON/TXT")
                    .clicked()
                {
                    self.prompt_import_jaspar_resource();
                    ui.close_menu();
                }
                if ui
                    .button("Save Project...")
                    .on_hover_text("Save current project state to disk")
                    .clicked()
                {
                    self.prompt_save_project();
                    ui.close_menu();
                }
                if ui
                    .button("Export DALG SVG...")
                    .on_hover_text("Export lineage graph as SVG")
                    .clicked()
                {
                    self.prompt_export_lineage_svg();
                    ui.close_menu();
                }
            });
            ui.menu_button("Edit", |ui| {
                let undo_resp = self.track_hover_status(
                    ui.add_enabled(
                        history_ops_enabled && undo_count > 0,
                        egui::Button::new("Undo"),
                    )
                    .on_hover_text("Undo the most recent operation-level state change"),
                    "Edit > Undo",
                );
                if undo_resp.clicked() {
                    self.undo_last_operation();
                    ui.close_menu();
                }
                let redo_resp = self.track_hover_status(
                    ui.add_enabled(
                        history_ops_enabled && redo_count > 0,
                        egui::Button::new("Redo"),
                    )
                    .on_hover_text("Redo the most recently undone operation"),
                    "Edit > Redo",
                );
                if redo_resp.clicked() {
                    self.redo_last_operation();
                    ui.close_menu();
                }
                if !history_ops_enabled {
                    ui.small("Undo/redo disabled while background jobs are running.");
                } else {
                    ui.small(format!("Undo {undo_count} | Redo {redo_count}"));
                }
                ui.separator();
                let palette_resp = self.track_hover_status(
                    ui.button("Command Palette...")
                        .on_hover_text("Open searchable command palette (Cmd/Ctrl+K)"),
                    "Edit > Command Palette",
                );
                if palette_resp.clicked() {
                    self.open_command_palette_dialog();
                    ui.close_menu();
                }
                let history_resp = self.track_hover_status(
                    ui.button("Operation History...")
                        .on_hover_text("Show operation history with undo/redo controls"),
                    "Edit > Operation History",
                );
                if history_resp.clicked() {
                    self.show_history_panel = true;
                    ui.close_menu();
                }
            });
            ui.menu_button("Settings", |ui| {
                if ui
                    .button("Configuration...")
                    .on_hover_text("Open global app configuration and graphics defaults")
                    .clicked()
                {
                    self.open_configuration_dialog();
                    ui.close_menu();
                }
            });
            ui.menu_button("Genome", |ui| {
                if ui
                    .button("Prepare Reference Genome...")
                    .on_hover_text(
                        "Download/index a reference genome for local extraction and BLAST",
                    )
                    .clicked()
                {
                    self.open_reference_genome_prepare_dialog();
                    ui.close_menu();
                }
                if ui
                    .button("Prepared References...")
                    .on_hover_text("Inspect prepared reference/helper genome installations")
                    .clicked()
                {
                    self.open_reference_genome_inspector_dialog();
                    ui.close_menu();
                }
                if ui
                    .button("Retrieve Genomic Sequence...")
                    .on_hover_text(
                        "Extract anchored region/gene sequence from a prepared reference",
                    )
                    .clicked()
                {
                    self.open_reference_genome_retrieve_dialog();
                    ui.close_menu();
                }
                if ui
                    .button("BLAST Genome Sequence...")
                    .on_hover_text("Run BLAST against prepared reference genome indices")
                    .clicked()
                {
                    self.open_reference_genome_blast_dialog();
                    ui.close_menu();
                }
                if ui
                    .button("Import Genome Track...")
                    .on_hover_text("Import BED/BigWig/VCF tracks onto anchored sequences")
                    .clicked()
                {
                    self.open_genome_bed_track_dialog();
                    ui.close_menu();
                }
                ui.separator();
                if ui
                    .button("Prepare Helper Genome...")
                    .on_hover_text("Prepare helper catalog genomes for extraction and BLAST")
                    .clicked()
                {
                    self.open_helper_genome_prepare_dialog();
                    ui.close_menu();
                }
                if ui
                    .button("Retrieve Helper Sequence...")
                    .on_hover_text("Extract sequence from prepared helper genomes")
                    .clicked()
                {
                    self.open_helper_genome_retrieve_dialog();
                    ui.close_menu();
                }
                if ui
                    .button("BLAST Helper Sequence...")
                    .on_hover_text("Run BLAST against prepared helper genome indices")
                    .clicked()
                {
                    self.open_helper_genome_blast_dialog();
                    ui.close_menu();
                }
            });
            ui.menu_button("Windows", |ui| {
                let jobs_panel_resp = self.track_hover_status(
                    ui.button(if self.show_jobs_panel {
                        "Hide Background Jobs"
                    } else {
                        "Show Background Jobs"
                    })
                    .on_hover_text("Toggle centralized panel for progress/cancel/retry summaries"),
                    "Window > Background Jobs Panel",
                );
                if jobs_panel_resp.clicked() {
                    self.show_jobs_panel = !self.show_jobs_panel;
                    ui.close_menu();
                }
                let history_panel_resp = self.track_hover_status(
                    ui.button(if self.show_history_panel {
                        "Hide Operation History"
                    } else {
                        "Show Operation History"
                    })
                    .on_hover_text("Toggle operation history panel with undo/redo"),
                    "Window > History Panel",
                );
                if history_panel_resp.clicked() {
                    self.show_history_panel = !self.show_history_panel;
                    ui.close_menu();
                }
                ui.separator();
                for entry in &open_window_entries {
                    let label = if entry.viewport_id == ViewportId::ROOT {
                        format!("Main: {}", entry.title)
                    } else {
                        entry.title.clone()
                    };
                    if ui
                        .button(label)
                        .on_hover_text(format!("Raise window: {}", entry.detail))
                        .clicked()
                    {
                        self.focus_window_viewport(ui.ctx(), entry.viewport_id);
                        ui.close_menu();
                    }
                }
            });
            ui.menu_button("Help", |ui| {
                if ui
                    .button("GUI Manual")
                    .on_hover_text("Open the GUI manual in the built-in help window")
                    .clicked()
                {
                    self.open_help_doc(HelpDoc::Gui);
                    ui.close_menu();
                }
                if ui
                    .button("CLI Manual")
                    .on_hover_text("Open the CLI manual in the built-in help window")
                    .clicked()
                {
                    self.open_help_doc(HelpDoc::Cli);
                    ui.close_menu();
                }
                if ui
                    .button("Shell Commands")
                    .on_hover_text("Open shell command reference generated from the glossary")
                    .clicked()
                {
                    self.open_help_doc(HelpDoc::Shell);
                    ui.close_menu();
                }
                ui.separator();
                if ui
                    .button("About GENtle")
                    .on_hover_text("Show version and build information")
                    .clicked()
                {
                    self.show_about_dialog = !about::show_native_about_panel();
                    ui.close_menu();
                }
            });
        });
    }

    fn show_window(&self, ctx: &egui::Context, id: ViewportId, window: Arc<RwLock<Window>>) {
        let windows_to_close = self.windows_to_close.clone();
        let window_number = self.get_window_number_from_id(id);
        let window_pos = Pos2 {
            x: window_number as f32 * 200.0,
            y: window_number as f32 * 200.0,
        };
        let window_title = window
            .read()
            .map(|w| w.name())
            .unwrap_or_else(|_| "GENtle".to_string());
        ctx.show_viewport_deferred(
            id,
            egui::ViewportBuilder::default()
                .with_title(window_title)
                // .with_maximized(true),
                .with_position(window_pos),
            move |ctx, class| {
                if class != egui::ViewportClass::Deferred {
                    eprintln!(
                        "W GENtleApp: unexpected viewport class, skipping deferred window update"
                    );
                    return;
                }
                if ctx.input(|i| i.viewport().focused.unwrap_or(false)) {
                    report_active_viewport_from_ui(id);
                }

                // Draw the window
                let update_result = catch_unwind(AssertUnwindSafe(|| {
                    if let Ok(mut w) = window.write() {
                        w.update(ctx);
                    } else {
                        eprintln!("W GENtleApp: window lock poisoned; skipping update");
                    }
                }));
                if update_result.is_err() {
                    eprintln!("E GENtleApp: recovered from panic while updating window");
                }

                // "Close window" action
                if ctx.input(|i| i.viewport().close_requested()) {
                    if let Ok(mut to_close) = windows_to_close.write() {
                        to_close.push(id);
                    } else {
                        eprintln!("W GENtleApp: close-queue lock poisoned");
                    }
                }
            },
        );
    }

    fn get_window_number_from_id(&self, id: ViewportId) -> usize {
        self.windows
            .keys()
            .enumerate()
            .find(|(_num, viewport_id)| **viewport_id == id)
            .map(|(num, _viewport_id)| num)
            .unwrap_or(0)
    }

    fn refresh_lineage_cache_if_needed(&mut self) {
        let stamp = self.current_lineage_change_stamp();
        if self.lineage_cache_valid && self.lineage_cache_stamp == stamp {
            return;
        }

        let (rows, lineage_edges, op_label_by_id, containers, arrangements) = {
            let engine = self.engine.read().unwrap();
            let state = engine.state();
            let mut op_created_count: HashMap<String, usize> = HashMap::new();
            let mut op_created_ids: HashMap<String, Vec<String>> = HashMap::new();
            let mut op_label_by_id: HashMap<String, String> = HashMap::new();
            for rec in engine.operation_log() {
                op_created_count.insert(rec.result.op_id.clone(), rec.result.created_seq_ids.len());
                op_created_ids.insert(rec.result.op_id.clone(), rec.result.created_seq_ids.clone());
                op_label_by_id.insert(rec.result.op_id.clone(), Self::summarize_operation(&rec.op));
            }

            let mut parents_by_node: HashMap<String, Vec<String>> = HashMap::new();
            for edge in &state.lineage.edges {
                let Some(parent_node) = state.lineage.nodes.get(&edge.from_node_id) else {
                    continue;
                };
                parents_by_node
                    .entry(edge.to_node_id.clone())
                    .or_default()
                    .push(parent_node.seq_id.clone());
            }

            let mut out: Vec<LineageRow> = state
                .lineage
                .nodes
                .values()
                .map(|node| {
                    let parents = parents_by_node
                        .get(&node.node_id)
                        .cloned()
                        .unwrap_or_default();
                    let (length, circular) = state
                        .sequences
                        .get(&node.seq_id)
                        .map(|dna| (dna.len(), dna.is_circular()))
                        .unwrap_or((0, false));
                    let display_name = state
                        .sequences
                        .get(&node.seq_id)
                        .and_then(|dna| dna.name().clone())
                        .unwrap_or_else(|| node.seq_id.clone());
                    let pool_size = node
                        .created_by_op
                        .as_ref()
                        .and_then(|op| op_created_count.get(op))
                        .cloned()
                        .unwrap_or(1);
                    let pool_members = node
                        .created_by_op
                        .as_ref()
                        .and_then(|op| op_created_ids.get(op))
                        .cloned()
                        .unwrap_or_else(|| vec![node.seq_id.clone()]);
                    LineageRow {
                        kind: LineageNodeKind::Sequence,
                        node_id: node.node_id.clone(),
                        seq_id: node.seq_id.clone(),
                        display_name,
                        origin: format!("{:?}", node.origin),
                        created_by_op: node
                            .created_by_op
                            .clone()
                            .unwrap_or_else(|| "-".to_string()),
                        created_at: node.created_at_unix_ms,
                        parents,
                        length,
                        circular,
                        pool_size,
                        pool_members,
                        arrangement_id: None,
                        arrangement_mode: None,
                        lane_container_ids: vec![],
                        ladders: vec![],
                    }
                })
                .collect();
            out.sort_by(|a, b| {
                a.created_at
                    .cmp(&b.created_at)
                    .then(a.node_id.cmp(&b.node_id))
            });
            let lineage_edges: Vec<(String, String, String)> = state
                .lineage
                .edges
                .iter()
                .filter_map(|e| {
                    let from = state.lineage.nodes.get(&e.from_node_id)?.node_id.clone();
                    let to = state.lineage.nodes.get(&e.to_node_id)?.node_id.clone();
                    Some((from, to, e.op_id.clone()))
                })
                .collect();
            let mut containers: Vec<ContainerRow> = state
                .container_state
                .containers
                .iter()
                .map(|(id, c)| ContainerRow {
                    container_id: id.clone(),
                    kind: format!("{:?}", c.kind),
                    member_count: c.members.len(),
                    representative: c.members.first().cloned().unwrap_or_default(),
                    members: c.members.clone(),
                })
                .collect();
            containers.sort_by(|a, b| a.container_id.cmp(&b.container_id));
            let mut arrangements: Vec<ArrangementRow> = state
                .container_state
                .arrangements
                .iter()
                .map(|(id, arrangement)| ArrangementRow {
                    arrangement_id: id.clone(),
                    mode: format!("{:?}", arrangement.mode),
                    name: arrangement.name.clone().unwrap_or_default(),
                    created_by_op: arrangement
                        .created_by_op
                        .clone()
                        .unwrap_or_else(|| "CreateArrangementSerial".to_string()),
                    created_at: arrangement.created_at_unix_ms,
                    lane_count: arrangement.lane_container_ids.len(),
                    lane_container_ids: arrangement.lane_container_ids.clone(),
                    ladders: arrangement.ladders.clone(),
                })
                .collect();
            arrangements.sort_by(|a, b| a.arrangement_id.cmp(&b.arrangement_id));
            (out, lineage_edges, op_label_by_id, containers, arrangements)
        };

        self.lineage_rows = rows;
        self.lineage_edges = lineage_edges;
        self.lineage_op_label_by_id = op_label_by_id;
        self.lineage_containers = containers;
        self.lineage_arrangements = arrangements;
        self.lineage_cache_stamp = stamp;
        self.lineage_cache_valid = true;
    }

    fn lineage_layout_positions(
        order_by_layer: &BTreeMap<usize, Vec<String>>,
    ) -> HashMap<String, usize> {
        let mut positions = HashMap::new();
        for nodes in order_by_layer.values() {
            for (index, node_id) in nodes.iter().enumerate() {
                positions.insert(node_id.clone(), index);
            }
        }
        positions
    }

    fn compute_lineage_dag_layout(
        rows: &[LineageRow],
        lineage_edges: &[(String, String, String)],
    ) -> (HashMap<String, (usize, usize)>, usize, usize) {
        if rows.is_empty() {
            return (HashMap::new(), 1, 1);
        }

        let mut row_index: HashMap<String, usize> = HashMap::new();
        for (index, row) in rows.iter().enumerate() {
            row_index.insert(row.node_id.clone(), index);
        }

        let mut parents_by_node: HashMap<String, Vec<String>> = HashMap::new();
        let mut children_by_node: HashMap<String, Vec<String>> = HashMap::new();
        let mut indegree_by_node: HashMap<String, usize> = HashMap::new();
        for row in rows {
            parents_by_node.insert(row.node_id.clone(), Vec::new());
            children_by_node.insert(row.node_id.clone(), Vec::new());
            indegree_by_node.insert(row.node_id.clone(), 0);
        }

        let mut seen_edges: HashSet<(String, String)> = HashSet::new();
        for (from_node, to_node, _op_id) in lineage_edges {
            if !row_index.contains_key(from_node) || !row_index.contains_key(to_node) {
                continue;
            }
            if !seen_edges.insert((from_node.clone(), to_node.clone())) {
                continue;
            }
            if let Some(children) = children_by_node.get_mut(from_node) {
                children.push(to_node.clone());
            }
            if let Some(parents) = parents_by_node.get_mut(to_node) {
                parents.push(from_node.clone());
            }
            if let Some(indegree) = indegree_by_node.get_mut(to_node) {
                *indegree = indegree.saturating_add(1);
            }
        }

        for parents in parents_by_node.values_mut() {
            parents.sort_by_key(|node_id| row_index.get(node_id).copied().unwrap_or(usize::MAX));
        }
        for children in children_by_node.values_mut() {
            children.sort_by_key(|node_id| row_index.get(node_id).copied().unwrap_or(usize::MAX));
        }

        let mut ready: Vec<String> = indegree_by_node
            .iter()
            .filter_map(|(node_id, indegree)| {
                if *indegree == 0 {
                    Some(node_id.clone())
                } else {
                    None
                }
            })
            .collect();
        ready.sort_by_key(|node_id| row_index.get(node_id).copied().unwrap_or(usize::MAX));

        let mut topo_order: Vec<String> = Vec::with_capacity(rows.len());
        let mut topo_seen: HashSet<String> = HashSet::with_capacity(rows.len());
        while !ready.is_empty() {
            let node_id = ready.remove(0);
            if !topo_seen.insert(node_id.clone()) {
                continue;
            }
            topo_order.push(node_id.clone());
            if let Some(children) = children_by_node.get(&node_id) {
                for child_id in children {
                    if let Some(indegree) = indegree_by_node.get_mut(child_id) {
                        if *indegree > 0 {
                            *indegree -= 1;
                            if *indegree == 0 {
                                ready.push(child_id.clone());
                            }
                        }
                    }
                }
            }
            ready.sort_by_key(|candidate| row_index.get(candidate).copied().unwrap_or(usize::MAX));
        }

        for row in rows {
            if topo_seen.insert(row.node_id.clone()) {
                topo_order.push(row.node_id.clone());
            }
        }

        let mut layer_by_node: HashMap<String, usize> = HashMap::new();
        for node_id in &topo_order {
            let layer = parents_by_node
                .get(node_id)
                .map(|parents| {
                    parents
                        .iter()
                        .filter_map(|parent_id| layer_by_node.get(parent_id).copied())
                        .max()
                        .map(|max_parent_layer| max_parent_layer + 1)
                        .unwrap_or(0)
                })
                .unwrap_or(0);
            layer_by_node.insert(node_id.clone(), layer);
        }

        let mut order_by_layer: BTreeMap<usize, Vec<String>> = BTreeMap::new();
        for node_id in &topo_order {
            let layer = layer_by_node.get(node_id).copied().unwrap_or(0);
            order_by_layer
                .entry(layer)
                .or_default()
                .push(node_id.clone());
        }
        for nodes in order_by_layer.values_mut() {
            nodes.sort_by_key(|node_id| row_index.get(node_id).copied().unwrap_or(usize::MAX));
        }

        let max_layer = order_by_layer.keys().copied().max().unwrap_or(0);
        let barycenter =
            |neighbors: &[String], positions: &HashMap<String, usize>| -> Option<f32> {
                let mut sum = 0.0f32;
                let mut count = 0usize;
                for node_id in neighbors {
                    if let Some(pos) = positions.get(node_id) {
                        sum += *pos as f32;
                        count += 1;
                    }
                }
                if count == 0 {
                    None
                } else {
                    Some(sum / count as f32)
                }
            };

        for _ in 0..6 {
            for layer in 1..=max_layer {
                let positions = Self::lineage_layout_positions(&order_by_layer);
                let Some(mut nodes) = order_by_layer.remove(&layer) else {
                    continue;
                };
                nodes.sort_by(|left, right| {
                    let left_score = parents_by_node
                        .get(left)
                        .and_then(|parents| barycenter(parents, &positions));
                    let right_score = parents_by_node
                        .get(right)
                        .and_then(|parents| barycenter(parents, &positions));
                    let fallback = row_index
                        .get(left)
                        .copied()
                        .unwrap_or(usize::MAX)
                        .cmp(&row_index.get(right).copied().unwrap_or(usize::MAX));
                    match (left_score, right_score) {
                        (Some(l), Some(r)) => l
                            .partial_cmp(&r)
                            .unwrap_or(std::cmp::Ordering::Equal)
                            .then(fallback),
                        (Some(_), None) => std::cmp::Ordering::Less,
                        (None, Some(_)) => std::cmp::Ordering::Greater,
                        (None, None) => fallback,
                    }
                });
                order_by_layer.insert(layer, nodes);
            }

            if max_layer == 0 {
                break;
            }

            for layer in (0..max_layer).rev() {
                let positions = Self::lineage_layout_positions(&order_by_layer);
                let Some(mut nodes) = order_by_layer.remove(&layer) else {
                    continue;
                };
                nodes.sort_by(|left, right| {
                    let left_score = children_by_node
                        .get(left)
                        .and_then(|children| barycenter(children, &positions));
                    let right_score = children_by_node
                        .get(right)
                        .and_then(|children| barycenter(children, &positions));
                    let fallback = row_index
                        .get(left)
                        .copied()
                        .unwrap_or(usize::MAX)
                        .cmp(&row_index.get(right).copied().unwrap_or(usize::MAX));
                    match (left_score, right_score) {
                        (Some(l), Some(r)) => l
                            .partial_cmp(&r)
                            .unwrap_or(std::cmp::Ordering::Equal)
                            .then(fallback),
                        (Some(_), None) => std::cmp::Ordering::Less,
                        (None, Some(_)) => std::cmp::Ordering::Greater,
                        (None, None) => fallback,
                    }
                });
                order_by_layer.insert(layer, nodes);
            }
        }

        let mut layout_by_node: HashMap<String, (usize, usize)> = HashMap::new();
        let mut max_rank_seen = 0usize;
        for (layer, nodes) in &order_by_layer {
            for node_id in nodes {
                let rank = row_index.get(node_id).copied().unwrap_or(0);
                max_rank_seen = max_rank_seen.max(rank);
                layout_by_node.insert(node_id.clone(), (*layer, rank));
            }
        }

        for row in rows {
            let fallback_rank = row_index.get(&row.node_id).copied().unwrap_or(0);
            max_rank_seen = max_rank_seen.max(fallback_rank);
            layout_by_node
                .entry(row.node_id.clone())
                .or_insert((0, fallback_rank));
        }

        let max_nodes_in_layer = max_rank_seen.saturating_add(1);
        (layout_by_node, max_layer + 1, max_nodes_in_layer)
    }

    fn compact_lineage_node_label(raw: &str, max_chars: usize) -> String {
        if max_chars == 0 {
            return String::new();
        }
        let char_count = raw.chars().count();
        if char_count <= max_chars {
            return raw.to_string();
        }
        let keep = max_chars.saturating_sub(3).max(1);
        let prefix: String = raw.chars().take(keep).collect();
        format!("{prefix}...")
    }

    fn compact_lineage_op_label(raw: &str) -> String {
        let mut head = raw.split(':').next().unwrap_or(raw).trim().to_string();
        if head.eq_ignore_ascii_case("Molecular weight filter (container)") {
            return "MW filter".to_string();
        }
        if head.eq_ignore_ascii_case("Merge containers by id") {
            return "Merge by id".to_string();
        }
        if head.len() > 20 {
            head.truncate(17);
            head.push_str("...");
        }
        head
    }

    fn render_main_lineage(&mut self, ui: &mut Ui) {
        self.refresh_lineage_cache_if_needed();

        ui.heading("Lineage Graph");
        ui.label(format!("Project: {}", self.current_project_name()));
        ui.label("Project-level sequence lineage (branch and merge aware)");
        ui.horizontal(|ui| {
            let prepare_resp = self.track_hover_status(
                ui.button("Prepare Reference Genome...")
                    .on_hover_text("Download and index reference genomes"),
                "Lineage > Prepare Reference Genome",
            );
            if prepare_resp.clicked() {
                self.open_reference_genome_prepare_dialog();
            }
            let retrieve_resp = self.track_hover_status(
                ui.button("Retrieve Genomic Sequence...")
                    .on_hover_text("Extract sequence regions from prepared genomes"),
                "Lineage > Retrieve Genomic Sequence",
            );
            if retrieve_resp.clicked() {
                self.open_reference_genome_retrieve_dialog();
            }
            let blast_resp = self.track_hover_status(
                ui.button("BLAST Genome Sequence...")
                    .on_hover_text("Run BLAST against prepared genome indices"),
                "Lineage > BLAST Genome Sequence",
            );
            if blast_resp.clicked() {
                self.open_reference_genome_blast_dialog();
            }
        });
        ui.horizontal(|ui| {
            let table = self.track_hover_status(
                ui.selectable_label(!self.lineage_graph_view, "Table")
                    .on_hover_text("Show lineage as table"),
                "Lineage > View Table",
            );
            if table.clicked() {
                self.lineage_graph_view = false;
            }
            let graph = self.track_hover_status(
                ui.selectable_label(self.lineage_graph_view, "Graph")
                    .on_hover_text("Show lineage as node graph"),
                "Lineage > View Graph",
            );
            if graph.clicked() {
                self.lineage_graph_view = true;
            }
        });
        ui.separator();

        let mut open_seq: Option<String> = None;
        let mut open_pool: Option<(String, Vec<String>)> = None;
        let mut open_lane_containers: Option<Vec<String>> = None;
        let mut export_container_gel: Option<(String, Vec<String>)> = None;
        let mut export_arrangement_gel: Option<(String, String)> = None;
        let mut select_candidate_from: Option<String> = None;
        let mut graph_zoom = self.lineage_graph_zoom.clamp(0.35, 4.0);
        let mut graph_area_height = self.lineage_graph_area_height.clamp(220.0, 2400.0);
        let mut container_area_height = self.lineage_container_area_height.clamp(120.0, 1600.0);
        let mut graph_scroll_offset = Vec2::new(
            if self.lineage_graph_scroll_offset.x.is_finite() {
                self.lineage_graph_scroll_offset.x.max(0.0)
            } else {
                0.0
            },
            if self.lineage_graph_scroll_offset.y.is_finite() {
                self.lineage_graph_scroll_offset.y.max(0.0)
            } else {
                0.0
            },
        );
        let mut graph_compact_labels = self.lineage_graph_compact_labels;
        let mut persist_workspace_after_frame = false;
        let mut graph_rows = self.lineage_rows.clone();
        let mut graph_edges = self.lineage_edges.clone();
        let mut graph_op_label_by_id = self.lineage_op_label_by_id.clone();
        let seq_node_by_seq_id: HashMap<String, String> = self
            .lineage_rows
            .iter()
            .map(|row| (row.seq_id.clone(), row.node_id.clone()))
            .collect();
        let container_members_by_id: HashMap<String, Vec<String>> = self
            .lineage_containers
            .iter()
            .map(|row| (row.container_id.clone(), row.members.clone()))
            .collect();
        for arrangement in &self.lineage_arrangements {
            let arrangement_node_id = format!("arr:{}", arrangement.arrangement_id);
            let mut source_node_ids: Vec<String> = vec![];
            let mut seen_sources: HashSet<String> = HashSet::new();
            for container_id in &arrangement.lane_container_ids {
                if let Some(members) = container_members_by_id.get(container_id) {
                    for seq_id in members {
                        if let Some(source_node_id) = seq_node_by_seq_id.get(seq_id).cloned() {
                            if seen_sources.insert(source_node_id.clone()) {
                                source_node_ids.push(source_node_id.clone());
                                graph_edges.push((
                                    source_node_id,
                                    arrangement_node_id.clone(),
                                    arrangement.created_by_op.clone(),
                                ));
                            }
                        }
                    }
                }
            }
            graph_op_label_by_id
                .entry(arrangement.created_by_op.clone())
                .or_insert_with(|| "Arrange serial lanes".to_string());
            graph_rows.push(LineageRow {
                kind: LineageNodeKind::Arrangement,
                node_id: arrangement_node_id,
                seq_id: arrangement.arrangement_id.clone(),
                display_name: if arrangement.name.trim().is_empty() {
                    arrangement.arrangement_id.clone()
                } else {
                    arrangement.name.clone()
                },
                origin: "Arrangement".to_string(),
                created_by_op: arrangement.created_by_op.clone(),
                created_at: arrangement.created_at,
                parents: source_node_ids,
                length: 0,
                circular: false,
                pool_size: 1,
                pool_members: vec![],
                arrangement_id: Some(arrangement.arrangement_id.clone()),
                arrangement_mode: Some(arrangement.mode.clone()),
                lane_container_ids: arrangement.lane_container_ids.clone(),
                ladders: arrangement.ladders.clone(),
            });
        }
        graph_rows.sort_by(|a, b| {
            a.created_at
                .cmp(&b.created_at)
                .then(a.node_id.cmp(&b.node_id))
        });
        if self.lineage_graph_view {
            if self.lineage_graph_offsets_synced_stamp != self.lineage_cache_stamp {
                let offset_count_before = self.lineage_graph_node_offsets.len();
                let active_node_ids = graph_rows
                    .iter()
                    .map(|row| row.node_id.clone())
                    .collect::<std::collections::HashSet<_>>();
                self.lineage_graph_node_offsets
                    .retain(|node_id, _| active_node_ids.contains(node_id));
                if self
                    .lineage_graph_drag_origin
                    .as_ref()
                    .is_some_and(|(node_id, _)| !active_node_ids.contains(node_id))
                {
                    self.lineage_graph_drag_origin = None;
                }
                if self
                    .lineage_graph_pan_origin
                    .as_ref()
                    .is_some_and(|_| active_node_ids.is_empty())
                {
                    self.lineage_graph_pan_origin = None;
                }
                if self
                    .lineage_graph_selected_node_id
                    .as_ref()
                    .is_some_and(|node_id| !active_node_ids.contains(node_id))
                {
                    self.lineage_graph_selected_node_id = None;
                }
                if self.lineage_graph_node_offsets.len() != offset_count_before {
                    persist_workspace_after_frame = true;
                }
                self.lineage_graph_offsets_synced_stamp = self.lineage_cache_stamp;
            }
            let mut request_fit_zoom = false;
            let mut request_fit_origin = false;
            ui.horizontal(|ui| {
                ui.label("Legend:");
                ui.colored_label(egui::Color32::from_rgb(90, 140, 210), "● single sequence");
                ui.colored_label(egui::Color32::from_rgb(180, 120, 70), "◆ pool");
                ui.colored_label(egui::Color32::from_rgb(108, 154, 122), "▭ arrangement");
                ui.separator();
                let zoom_out_resp = self.track_hover_status(
                    ui.button("−").on_hover_text("Zoom out"),
                    "Lineage Graph > Zoom Out",
                );
                if zoom_out_resp.clicked() {
                    graph_zoom = (graph_zoom / 1.15).clamp(0.35, 4.0);
                }
                let zoom_in_resp = self.track_hover_status(
                    ui.button("+").on_hover_text("Zoom in"),
                    "Lineage Graph > Zoom In",
                );
                if zoom_in_resp.clicked() {
                    graph_zoom = (graph_zoom * 1.15).clamp(0.35, 4.0);
                }
                let zoom_reset_resp = self.track_hover_status(
                    ui.button("Reset").on_hover_text("Reset zoom"),
                    "Lineage Graph > Zoom Reset",
                );
                if zoom_reset_resp.clicked() {
                    graph_zoom = 1.0;
                }
                let fit_resp = self.track_hover_status(
                    ui.button("Fit")
                        .on_hover_text("Fit full graph content into current graph area"),
                    "Lineage Graph > Fit",
                );
                if fit_resp.clicked()
                {
                    request_fit_zoom = true;
                    request_fit_origin = true;
                }
                let reset_layout_resp = self.track_hover_status(
                    ui.button("Reset Layout")
                        .on_hover_text("Reset manually moved node positions"),
                    "Lineage Graph > Reset Layout",
                );
                if reset_layout_resp.clicked()
                {
                    self.lineage_graph_node_offsets.clear();
                    self.lineage_graph_drag_origin = None;
                    self.lineage_graph_pan_origin = None;
                    persist_workspace_after_frame = true;
                }
                ui.separator();
                if ui
                    .checkbox(&mut graph_compact_labels, "Compact labels")
                    .on_hover_text(
                        "Reduce label density in crowded graphs by simplifying operation and node labels",
                    )
                    .changed()
                {
                    persist_workspace_after_frame = true;
                }
                ui.separator();
                ui.add(
                    egui::Slider::new(&mut graph_zoom, 0.35..=4.0)
                        .logarithmic(true)
                        .text("Zoom"),
                );
                ui.label(format!("{:.0}%", graph_zoom * 100.0));
                if ui
                    .add(
                        egui::DragValue::new(&mut graph_area_height)
                            .range(220.0..=2400.0)
                            .speed(2.0)
                            .prefix("Graph h "),
                    )
                    .on_hover_text("Persisted preferred height of the graph area")
                    .changed()
                {
                    persist_workspace_after_frame = true;
                }
            });
            ui.separator();
            let rows = &graph_rows;
            let lineage_edges = &graph_edges;
            let op_label_by_id = &graph_op_label_by_id;
            let graph_resize_max_height = ui.available_height().max(220.0);
            let graph_resize_width = ui.available_width().max(360.0);
            egui::Resize::default()
                .id_salt("lineage_graph_area_resize")
                .default_width(graph_resize_width)
                .min_width(graph_resize_width)
                .max_width(graph_resize_width)
                .default_height(graph_area_height.min(graph_resize_max_height))
                .min_height(220.0)
                .max_height(graph_resize_max_height)
                .resizable(egui::Vec2b::new(false, true))
                .show(ui, |ui| {
                    ui.set_min_width(graph_resize_width);
                    ui.set_max_width(graph_resize_width);
                    graph_area_height = ui.max_rect().height().clamp(220.0, 2400.0);
                    ui.small(
                        "Resize area by dragging lower edge. Hold Space + drag background to pan. Drag nodes to reposition. Cmd/Ctrl+scroll zooms.",
                    );
                    let graph_scroll_output = egui::ScrollArea::both()
                        .id_salt("lineage_graph_scroll")
                        .auto_shrink([false, false])
                        .scroll_offset(graph_scroll_offset)
                        .drag_to_scroll(false)
                        .max_height(ui.available_height())
                        .show(ui, |ui| {
                            let (layout_by_node, layer_count, max_nodes_in_layer) =
                                Self::compute_lineage_dag_layout(rows, lineage_edges);
                            let base_width = (layer_count.max(1) as f32) * 220.0 + 220.0;
                            let base_height = (max_nodes_in_layer.max(1) as f32) * 110.0 + 300.0;
                            if request_fit_zoom {
                                let available = ui.available_size();
                                let fit_x = ((available.x - 24.0).max(120.0) / base_width.max(1.0))
                                    .max(0.01);
                                let fit_y =
                                    ((available.y - 24.0).max(120.0) / base_height.max(1.0))
                                        .max(0.01);
                                graph_zoom = fit_x.min(fit_y).clamp(0.35, 4.0);
                                request_fit_zoom = false;
                            }
                            if request_fit_origin {
                                graph_scroll_offset = Vec2::ZERO;
                                self.lineage_graph_pan_origin = None;
                                request_fit_origin = false;
                            }
                            let width = (base_width * graph_zoom + 280.0 * graph_zoom)
                                .max(ui.available_width());
                            let height = base_height * graph_zoom;
                            let (resp, painter) = ui.allocate_painter(
                                Vec2::new(width, height),
                                egui::Sense::click_and_drag(),
                            );
                            if resp.hovered() {
                                let (zoom_modifier, scroll_y) = ui.input(|i| {
                                    (
                                        i.modifiers.ctrl || i.modifiers.command,
                                        i.smooth_scroll_delta.y,
                                    )
                                });
                                if zoom_modifier && scroll_y.abs() > f32::EPSILON {
                                    let factor = (1.0 + scroll_y * 0.0015).clamp(0.8, 1.25);
                                    graph_zoom = (graph_zoom * factor).clamp(0.35, 4.0);
                                }
                            }
                            let dense_graph = rows.len() >= 22 || lineage_edges.len() >= 30;
                            let simplify_labels = graph_compact_labels && dense_graph;
                            let edge_label_stride = if simplify_labels {
                                if lineage_edges.len() > 180 {
                                    4
                                } else if lineage_edges.len() > 120 {
                                    3
                                } else if lineage_edges.len() > 60 {
                                    2
                                } else {
                                    1
                                }
                            } else {
                                1
                            };
                            let rect = resp.rect;
                            let edge_stroke_width = (1.0 * graph_zoom).clamp(1.0, 2.5);
                            let node_id_font_size = (10.0 * graph_zoom).clamp(8.0, 15.0);
                            let op_font_size = if simplify_labels {
                                (9.0 * graph_zoom).clamp(8.0, 14.0)
                            } else {
                                (10.0 * graph_zoom).clamp(9.0, 16.0)
                            };
                            let name_font_size = if simplify_labels {
                                (11.0 * graph_zoom).clamp(9.0, 15.0)
                            } else {
                                (12.0 * graph_zoom).clamp(10.0, 18.0)
                            };
                            let details_font_size = (10.0 * graph_zoom).clamp(9.0, 15.0);
                            let node_radius = 16.0 * graph_zoom;
                            let mut pos_by_node: HashMap<String, Pos2> = HashMap::new();
                            for (fallback_rank, row) in rows.iter().enumerate() {
                                let (layer, rank) = layout_by_node
                                    .get(&row.node_id)
                                    .copied()
                                    .unwrap_or((0, fallback_rank));
                                let manual_offset = self
                                    .lineage_graph_node_offsets
                                    .get(&row.node_id)
                                    .copied()
                                    .unwrap_or(Vec2::ZERO);
                                let x = rect.left()
                                    + (120.0 + layer as f32 * 220.0) * graph_zoom
                                    + manual_offset.x;
                                let y = rect.top()
                                    + (120.0 + rank as f32 * 110.0) * graph_zoom
                                    + manual_offset.y;
                                pos_by_node.insert(row.node_id.clone(), Pos2::new(x, y));
                            }
                            let pointer = resp.interact_pointer_pos();
                            let is_node_hit = |row: &LineageRow, pointer: Pos2| -> bool {
                                let Some(pos) = pos_by_node.get(&row.node_id).copied() else {
                                    return false;
                                };
                                let glyph_hit = match row.kind {
                                    LineageNodeKind::Arrangement => egui::Rect::from_center_size(
                                        pos,
                                        Vec2::new(56.0 * graph_zoom, 28.0 * graph_zoom),
                                    )
                                    .contains(pointer),
                                    LineageNodeKind::Sequence if row.pool_size > 1 => {
                                        egui::Rect::from_center_size(
                                            pos,
                                            Vec2::new(34.0 * graph_zoom, 28.0 * graph_zoom),
                                        )
                                        .contains(pointer)
                                    }
                                    LineageNodeKind::Sequence => pointer.distance(pos) <= 19.0 * graph_zoom,
                                };
                                if glyph_hit {
                                    return true;
                                }
                                let display_name = if simplify_labels {
                                    Self::compact_lineage_node_label(&row.display_name, 26)
                                } else {
                                    row.display_name.clone()
                                };
                                let name_galley = painter.layout_no_wrap(
                                    display_name,
                                    egui::FontId::proportional(name_font_size),
                                    egui::Color32::BLACK,
                                );
                                let name_anchor = pos + Vec2::new(22.0 * graph_zoom, -4.0 * graph_zoom);
                                let name_rect = egui::Rect::from_min_size(
                                    Pos2::new(name_anchor.x, name_anchor.y - name_galley.size().y),
                                    name_galley.size(),
                                )
                                .expand(3.0 * graph_zoom);
                                if name_rect.contains(pointer) {
                                    return true;
                                }
                                if !simplify_labels {
                                    let detail_text = match row.kind {
                                        LineageNodeKind::Arrangement => {
                                            let ladders = if row.ladders.is_empty() {
                                                "auto".to_string()
                                            } else {
                                                row.ladders.join(" + ")
                                            };
                                            format!(
                                                "{} lane(s) | mode={} | ladders={}",
                                                row.lane_container_ids.len(),
                                                row.arrangement_mode
                                                    .as_deref()
                                                    .unwrap_or("Serial")
                                                    .to_lowercase(),
                                                ladders
                                            )
                                        }
                                        LineageNodeKind::Sequence if row.pool_size > 1 => {
                                            format!(
                                                "{} ({} bp) | pool={}",
                                                row.seq_id, row.length, row.pool_size
                                            )
                                        }
                                        LineageNodeKind::Sequence => {
                                            format!("{} ({} bp)", row.seq_id, row.length)
                                        }
                                    };
                                    let details_galley = painter.layout_no_wrap(
                                        detail_text,
                                        egui::FontId::proportional(details_font_size),
                                        egui::Color32::BLACK,
                                    );
                                    let details_anchor =
                                        pos + Vec2::new(22.0 * graph_zoom, 10.0 * graph_zoom);
                                    let details_rect =
                                        egui::Rect::from_min_size(details_anchor, details_galley.size())
                                            .expand(3.0 * graph_zoom);
                                    if details_rect.contains(pointer) {
                                        return true;
                                    }
                                }
                                false
                            };
                            let hovered_node_id = pointer.and_then(|pointer| {
                                rows.iter()
                                    .filter_map(|row| {
                                        if !is_node_hit(row, pointer) {
                                            return None;
                                        }
                                        pos_by_node
                                            .get(&row.node_id)
                                            .map(|pos| (row.node_id.clone(), pointer.distance_sq(*pos)))
                                    })
                                    .min_by(|a, b| {
                                        a.1.partial_cmp(&b.1)
                                            .unwrap_or(std::cmp::Ordering::Equal)
                                    })
                                    .map(|(node_id, _)| node_id)
                            });
                            let mut used_label_rects: Vec<egui::Rect> = Vec::new();
                            let mut op_label_galleys: HashMap<String, Arc<egui::Galley>> =
                                HashMap::new();
                            for (edge_idx, (from_node, to_node, op_id)) in
                                lineage_edges.iter().enumerate()
                            {
                                let Some(from) = pos_by_node.get(from_node) else {
                                    continue;
                                };
                                let Some(to) = pos_by_node.get(to_node) else {
                                    continue;
                                };
                                painter.line_segment(
                                    [*from, *to],
                                    egui::Stroke::new(edge_stroke_width, egui::Color32::GRAY),
                                );
                                if edge_label_stride > 1 && edge_idx % edge_label_stride != 0 {
                                    continue;
                                }
                                let mid = Pos2::new((from.x + to.x) * 0.5, (from.y + to.y) * 0.5);
                                let op_label = op_label_by_id
                                    .get(op_id)
                                    .cloned()
                                    .unwrap_or_else(|| op_id.clone());
                                let display = if simplify_labels {
                                    Self::compact_lineage_op_label(&op_label)
                                } else {
                                    op_label
                                };
                                let galley = op_label_galleys
                                    .entry(display.clone())
                                    .or_insert_with(|| {
                                        painter.layout_no_wrap(
                                            display.clone(),
                                            egui::FontId::proportional(op_font_size),
                                            egui::Color32::BLACK,
                                        )
                                    })
                                    .clone();
                                let edge = *to - *from;
                                let edge_len = edge.length();
                                if edge_len < 0.1 {
                                    continue;
                                }
                                let edge_dir = edge / edge_len;
                                let perp = Vec2::new(-edge_dir.y, edge_dir.x);
                                let mut placed = None;
                                for idx in 0..14 {
                                    let sign = if idx % 2 == 0 { 1.0 } else { -1.0 };
                                    let step = (idx / 2) as f32;
                                    let candidate_center =
                                        mid + perp * (12.0 + step * 12.0) * graph_zoom * sign;
                                    let candidate_rect = egui::Rect::from_center_size(
                                        candidate_center,
                                        Vec2::new(
                                            galley.size().x + 10.0 * graph_zoom,
                                            galley.size().y + 6.0 * graph_zoom,
                                        ),
                                    );
                                    if !used_label_rects
                                        .iter()
                                        .any(|r| r.intersects(candidate_rect))
                                    {
                                        placed = Some((candidate_center, candidate_rect));
                                        break;
                                    }
                                }
                                let (label_center, bg_rect) = placed.unwrap_or_else(|| {
                                    let fallback_center = mid + perp * 12.0 * graph_zoom;
                                    let rect = egui::Rect::from_center_size(
                                        fallback_center,
                                        Vec2::new(
                                            galley.size().x + 10.0 * graph_zoom,
                                            galley.size().y + 6.0 * graph_zoom,
                                        ),
                                    );
                                    (fallback_center, rect)
                                });
                                used_label_rects.push(bg_rect);
                                painter.rect_filled(
                                    bg_rect,
                                    3.0 * graph_zoom,
                                    egui::Color32::from_rgba_premultiplied(245, 245, 245, 235),
                                );
                                painter.text(
                                    label_center,
                                    egui::Align2::CENTER_CENTER,
                                    display,
                                    egui::FontId::proportional(op_font_size),
                                    egui::Color32::BLACK,
                                );
                            }
                            for row in rows {
                                let Some(pos) = pos_by_node.get(&row.node_id).cloned() else {
                                    continue;
                                };
                                let is_selected = self
                                    .lineage_graph_selected_node_id
                                    .as_ref()
                                    .is_some_and(|node_id| node_id == &row.node_id);
                                let is_hovered = hovered_node_id
                                    .as_ref()
                                    .is_some_and(|node_id| node_id == &row.node_id);
                                let highlight_stroke = if is_selected {
                                    egui::Stroke::new(
                                        (2.4 * graph_zoom).clamp(1.6, 4.0),
                                        egui::Color32::from_rgb(250, 220, 80),
                                    )
                                } else if is_hovered {
                                    egui::Stroke::new(
                                        (2.0 * graph_zoom).clamp(1.4, 3.2),
                                        egui::Color32::from_rgb(230, 230, 150),
                                    )
                                } else {
                                    egui::Stroke::new(
                                        edge_stroke_width,
                                        egui::Color32::from_rgb(235, 196, 150),
                                    )
                                };
                                match row.kind {
                                    LineageNodeKind::Arrangement => {
                                        let rect = egui::Rect::from_center_size(
                                            pos,
                                            Vec2::new(56.0 * graph_zoom, 28.0 * graph_zoom),
                                        );
                                        painter.rect_filled(
                                            rect,
                                            5.0 * graph_zoom,
                                            if is_selected {
                                                egui::Color32::from_rgb(98, 140, 112)
                                            } else {
                                                egui::Color32::from_rgb(108, 154, 122)
                                            },
                                        );
                                        painter.rect_stroke(
                                            rect,
                                            5.0 * graph_zoom,
                                            highlight_stroke,
                                        );
                                    }
                                    LineageNodeKind::Sequence if row.pool_size > 1 => {
                                        let points = vec![
                                            pos + Vec2::new(0.0, -16.0 * graph_zoom),
                                            pos + Vec2::new(16.0 * graph_zoom, 0.0),
                                            pos + Vec2::new(0.0, 16.0 * graph_zoom),
                                            pos + Vec2::new(-16.0 * graph_zoom, 0.0),
                                        ];
                                        painter.add(egui::Shape::convex_polygon(
                                            points,
                                            if is_selected {
                                                egui::Color32::from_rgb(205, 140, 80)
                                            } else {
                                                egui::Color32::from_rgb(180, 120, 70)
                                            },
                                            highlight_stroke,
                                        ));
                                    }
                                    LineageNodeKind::Sequence => {
                                        painter.circle_filled(
                                            pos,
                                            node_radius,
                                            if is_selected {
                                                egui::Color32::from_rgb(70, 125, 215)
                                            } else {
                                                egui::Color32::from_rgb(90, 140, 210)
                                            },
                                        );
                                        painter.circle_stroke(pos, node_radius, highlight_stroke);
                                    }
                                }
                                let node_id_label = match row.kind {
                                    LineageNodeKind::Arrangement => row
                                        .arrangement_id
                                        .as_ref()
                                        .map(|id| Self::compact_lineage_node_label(id, 12))
                                        .unwrap_or_else(|| {
                                            Self::compact_lineage_node_label(&row.node_id, 12)
                                        }),
                                    LineageNodeKind::Sequence => {
                                        Self::compact_lineage_node_label(&row.node_id, 10)
                                    }
                                };
                                painter.text(
                                    pos,
                                    egui::Align2::CENTER_CENTER,
                                    node_id_label,
                                    egui::FontId::monospace(node_id_font_size),
                                    egui::Color32::WHITE,
                                );
                                let display_name = if simplify_labels {
                                    Self::compact_lineage_node_label(&row.display_name, 26)
                                } else {
                                    row.display_name.clone()
                                };
                                painter.text(
                                    pos + Vec2::new(22.0 * graph_zoom, -4.0 * graph_zoom),
                                    egui::Align2::LEFT_BOTTOM,
                                    display_name,
                                    egui::FontId::proportional(name_font_size),
                                    egui::Color32::BLACK,
                                );
                                if !simplify_labels {
                                    let detail_text = match row.kind {
                                        LineageNodeKind::Arrangement => {
                                            let ladders = if row.ladders.is_empty() {
                                                "auto".to_string()
                                            } else {
                                                row.ladders.join(" + ")
                                            };
                                            format!(
                                                "{} lane(s) | mode={} | ladders={}",
                                                row.lane_container_ids.len(),
                                                row.arrangement_mode
                                                    .as_deref()
                                                    .unwrap_or("Serial")
                                                    .to_lowercase(),
                                                ladders
                                            )
                                        }
                                        LineageNodeKind::Sequence if row.pool_size > 1 => {
                                            format!(
                                                "{} ({} bp) | pool={}",
                                                row.seq_id, row.length, row.pool_size
                                            )
                                        }
                                        LineageNodeKind::Sequence => {
                                            format!("{} ({} bp)", row.seq_id, row.length)
                                        }
                                    };
                                    painter.text(
                                        pos + Vec2::new(22.0 * graph_zoom, 10.0 * graph_zoom),
                                        egui::Align2::LEFT_TOP,
                                        detail_text,
                                        egui::FontId::proportional(details_font_size),
                                        egui::Color32::BLACK,
                                    );
                                }
                            }

                            let hit_row = hovered_node_id
                                .as_ref()
                                .and_then(|node_id| rows.iter().find(|row| &row.node_id == node_id));
                            if let Some(row) = hit_row {
                                let mut hover_pool_range: Option<(usize, usize)> = None;
                                let mut hover_ladder_hint: Option<String> = None;
                                if row.kind == LineageNodeKind::Sequence && row.pool_size > 1 {
                                    let member_lengths: Vec<(String, usize)> = {
                                        let engine = self.engine.read().unwrap();
                                        row.pool_members
                                            .iter()
                                            .filter_map(|seq_id| {
                                                engine
                                                    .state()
                                                    .sequences
                                                    .get(seq_id)
                                                    .map(|dna| (seq_id.clone(), dna.len()))
                                            })
                                            .collect()
                                    };
                                    if !member_lengths.is_empty() {
                                        let min_bp = member_lengths
                                            .iter()
                                            .map(|(_, bp)| *bp)
                                            .min()
                                            .unwrap_or(0);
                                        let max_bp = member_lengths
                                            .iter()
                                            .map(|(_, bp)| *bp)
                                            .max()
                                            .unwrap_or(min_bp);
                                        hover_pool_range = Some((min_bp, max_bp));
                                        if let Ok(layout) =
                                            crate::pool_gel::build_pool_gel_layout(&member_lengths, &[])
                                        {
                                            if !layout.selected_ladders.is_empty() {
                                                hover_ladder_hint = Some(layout.selected_ladders.join(" + "));
                                            }
                                        }
                                    }
                                }
                                let tooltip_member_preview = if row.kind == LineageNodeKind::Sequence
                                    && row.pool_size > 1
                                {
                                    row.pool_members
                                        .iter()
                                        .take(6)
                                        .cloned()
                                        .collect::<Vec<_>>()
                                        .join(", ")
                                } else {
                                    String::new()
                                };
                                resp.clone().on_hover_ui_at_pointer(|ui| {
                                    ui.strong(&row.display_name);
                                    match row.kind {
                                        LineageNodeKind::Arrangement => {
                                            let ladders = if row.ladders.is_empty() {
                                                "auto".to_string()
                                            } else {
                                                row.ladders.join(" + ")
                                            };
                                            ui.monospace(format!(
                                                "{} | mode={} | lanes={} | ladders={}",
                                                row.arrangement_id
                                                    .as_deref()
                                                    .unwrap_or(&row.seq_id),
                                                row.arrangement_mode
                                                    .as_deref()
                                                    .unwrap_or("Serial")
                                                    .to_lowercase(),
                                                row.lane_container_ids.len(),
                                                ladders
                                            ));
                                        }
                                        LineageNodeKind::Sequence => {
                                            ui.monospace(format!(
                                                "{} | {} bp | {}",
                                                row.seq_id,
                                                row.length,
                                                if row.circular { "circular" } else { "linear" }
                                            ));
                                        }
                                    }
                                    ui.small(format!("node={} | origin={}", row.node_id, row.origin));
                                    ui.small(format!("parents={} | op={}", row.parents.len(), row.created_by_op));
                                    match row.kind {
                                        LineageNodeKind::Arrangement => {
                                            if !row.lane_container_ids.is_empty() {
                                                ui.separator();
                                                let preview = row
                                                    .lane_container_ids
                                                    .iter()
                                                    .take(6)
                                                    .cloned()
                                                    .collect::<Vec<_>>()
                                                    .join(", ");
                                                ui.small(format!("lane containers: {preview}"));
                                                if row.lane_container_ids.len() > 6 {
                                                    ui.small(format!(
                                                        "... and {} more",
                                                        row.lane_container_ids.len() - 6
                                                    ));
                                                }
                                            }
                                        }
                                        LineageNodeKind::Sequence if row.pool_size > 1 => {
                                            ui.separator();
                                            ui.small(format!("pool members={}", row.pool_size));
                                            if let Some((min_bp, max_bp)) = hover_pool_range {
                                                ui.small(format!(
                                                    "pool range={}..{} bp",
                                                    min_bp, max_bp
                                                ));
                                            }
                                            if let Some(ladders) = &hover_ladder_hint {
                                                ui.small(format!("suggested ladders={}", ladders));
                                            }
                                            if !tooltip_member_preview.is_empty() {
                                                ui.small(format!(
                                                    "members: {}",
                                                    tooltip_member_preview
                                                ));
                                            }
                                        }
                                        LineageNodeKind::Sequence => {}
                                    }
                                });
                            }
                            let space_pan_requested = ui.input(|i| i.key_down(Key::Space));
                            if resp.drag_started() {
                                if space_pan_requested && hit_row.is_none() {
                                    self.lineage_graph_pan_origin = Some(graph_scroll_offset);
                                    self.lineage_graph_drag_origin = None;
                                } else if let Some(row) = hit_row {
                                    let start_offset = self
                                        .lineage_graph_node_offsets
                                        .get(&row.node_id)
                                        .copied()
                                        .unwrap_or(Vec2::ZERO);
                                    self.lineage_graph_drag_origin =
                                        Some((row.node_id.clone(), start_offset));
                                    self.lineage_graph_pan_origin = None;
                                    self.lineage_graph_selected_node_id = Some(row.node_id.clone());
                                } else {
                                    self.lineage_graph_drag_origin = None;
                                    self.lineage_graph_pan_origin = None;
                                }
                            }
                            if let Some(start_scroll_offset) = self.lineage_graph_pan_origin {
                                if resp.dragged() {
                                    let next_offset = start_scroll_offset - resp.drag_delta();
                                    graph_scroll_offset = Vec2::new(
                                        next_offset.x.max(0.0),
                                        next_offset.y.max(0.0),
                                    );
                                }
                                if resp.drag_stopped() {
                                    self.lineage_graph_pan_origin = None;
                                    persist_workspace_after_frame = true;
                                }
                            }
                            if self.lineage_graph_pan_origin.is_none() {
                                if let Some((node_id, start_offset)) =
                                    self.lineage_graph_drag_origin.clone()
                                {
                                    if resp.dragged() {
                                        self.lineage_graph_node_offsets
                                            .insert(node_id.clone(), start_offset + resp.drag_delta());
                                    }
                                    if resp.drag_stopped() {
                                        self.lineage_graph_drag_origin = None;
                                        persist_workspace_after_frame = true;
                                    }
                                }
                            }
                            if self.lineage_graph_pan_origin.is_some() {
                                ui.output_mut(|o| o.cursor_icon = egui::CursorIcon::Grabbing);
                            } else if space_pan_requested && hit_row.is_none() {
                                ui.output_mut(|o| o.cursor_icon = egui::CursorIcon::Grab);
                            } else if self.lineage_graph_drag_origin.is_some() {
                                ui.output_mut(|o| o.cursor_icon = egui::CursorIcon::Grabbing);
                            }
                            if self.lineage_graph_pan_origin.is_none()
                                && self.lineage_graph_drag_origin.is_none()
                            {
                                if resp.clicked() {
                                    self.lineage_graph_selected_node_id =
                                        hit_row.map(|row| row.node_id.clone());
                                }
                                if let Some(row) = hit_row {
                                    if resp.double_clicked() {
                                        match row.kind {
                                            LineageNodeKind::Arrangement => {
                                                if !row.lane_container_ids.is_empty() {
                                                    open_lane_containers =
                                                        Some(row.lane_container_ids.clone());
                                                }
                                            }
                                            LineageNodeKind::Sequence if row.pool_size > 1 => {
                                                open_pool = Some((
                                                    row.seq_id.clone(),
                                                    row.pool_members.clone(),
                                                ));
                                            }
                                            LineageNodeKind::Sequence => {
                                                open_seq = Some(row.seq_id.clone());
                                            }
                                        }
                                    }
                                }
                            }
                            if self.lineage_graph_drag_origin.is_some()
                                && !ui.input(|i| i.pointer.primary_down())
                            {
                                self.lineage_graph_drag_origin = None;
                                persist_workspace_after_frame = true;
                            }
                            if self.lineage_graph_pan_origin.is_some()
                                && !ui.input(|i| i.pointer.primary_down())
                            {
                                self.lineage_graph_pan_origin = None;
                                persist_workspace_after_frame = true;
                            }
                        });
                    let max_scroll_x =
                        (graph_scroll_output.content_size.x - graph_scroll_output.inner_rect.width())
                            .max(0.0);
                    let max_scroll_y = (graph_scroll_output.content_size.y
                        - graph_scroll_output.inner_rect.height())
                    .max(0.0);
                    if self.lineage_graph_pan_origin.is_some() {
                        graph_scroll_offset.x = graph_scroll_offset.x.clamp(0.0, max_scroll_x);
                        graph_scroll_offset.y = graph_scroll_offset.y.clamp(0.0, max_scroll_y);
                    } else {
                        let mut measured_offset = graph_scroll_output.state.offset;
                        measured_offset.x = measured_offset.x.clamp(0.0, max_scroll_x);
                        measured_offset.y = measured_offset.y.clamp(0.0, max_scroll_y);
                        graph_scroll_offset = measured_offset;
                    }
                });
        } else {
            egui::ScrollArea::both()
                .id_salt("lineage_table_scroll")
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    egui::Grid::new("lineage_table_grid")
                        .striped(true)
                        .min_col_width(90.0)
                        .show(ui, |ui| {
                            ui.strong("Node");
                            ui.strong("Sequence");
                            ui.strong("Parents");
                            ui.strong("Origin");
                            ui.strong("Op");
                            ui.strong("Length");
                            ui.strong("Topology");
                            ui.strong("Action");
                            ui.end_row();
                            for row in &self.lineage_rows {
                                ui.monospace(&row.node_id);
                                if ui
                                    .button(&row.seq_id)
                                    .on_hover_text("Open this sequence in a dedicated window")
                                    .clicked()
                                {
                                    open_seq = Some(row.seq_id.clone());
                                }
                                ui.label(if row.parents.is_empty() {
                                    "-".to_string()
                                } else {
                                    row.parents.join(" + ")
                                });
                                ui.label(&row.origin);
                                ui.monospace(&row.created_by_op);
                                ui.monospace(format!("{} bp", row.length));
                                ui.label(if row.circular { "circular" } else { "linear" });
                                if ui
                                    .button("Select")
                                    .on_hover_text(
                                        "Run candidate selection operation using this sequence as input",
                                    )
                                    .clicked()
                                {
                                    select_candidate_from = Some(row.seq_id.clone());
                                }
                                ui.end_row();
                            }
                        });
                });
        }
        ui.separator();
        ui.heading("Containers");
        ui.label("Container-level view of candidate sequence sets");
        ui.horizontal(|ui| {
            if ui
                .add(
                    egui::DragValue::new(&mut container_area_height)
                        .range(120.0..=1600.0)
                        .speed(2.0)
                        .prefix("Container h "),
                )
                .on_hover_text("Persisted preferred height of the containers area")
                .changed()
            {
                persist_workspace_after_frame = true;
            }
        });
        let container_resize_max_height = ui.available_height().max(120.0);
        let container_resize_width = ui.available_width().max(360.0);
        egui::Resize::default()
            .id_salt("lineage_container_area_resize")
            .default_width(container_resize_width)
            .min_width(container_resize_width)
            .max_width(container_resize_width)
            .default_height(container_area_height.min(container_resize_max_height))
            .min_height(120.0)
            .max_height(container_resize_max_height)
            .resizable(egui::Vec2b::new(false, true))
            .show(ui, |ui| {
                ui.set_min_width(container_resize_width);
                ui.set_max_width(container_resize_width);
                container_area_height = ui.max_rect().height().clamp(120.0, 1600.0);
                egui::ScrollArea::both()
                    .id_salt("lineage_container_grid_scroll")
                    .auto_shrink([false, false])
                    .max_height(ui.available_height())
                    .show(ui, |ui| {
                        egui::Grid::new("container_grid")
                            .striped(true)
                            .show(ui, |ui| {
                                ui.strong("Container");
                                ui.strong("Kind");
                                ui.strong("Members");
                                ui.strong("Representative");
                                ui.strong("Action");
                                ui.end_row();
                                for c in &self.lineage_containers {
                                    ui.monospace(&c.container_id);
                                    ui.label(&c.kind);
                                    ui.monospace(format!("{}", c.member_count));
                                    ui.monospace(&c.representative);
                                    ui.horizontal(|ui| {
                                        if c.member_count > 1 {
                                            if ui
                                                .button("Open Pool")
                                                .on_hover_text("Open this container as a pool view")
                                                .clicked()
                                            {
                                                open_pool = Some((
                                                    c.representative.clone(),
                                                    c.members.clone(),
                                                ));
                                            }
                                        } else if !c.representative.is_empty() {
                                            if ui
                                                .button("Open Seq")
                                                .on_hover_text("Open this representative sequence")
                                                .clicked()
                                            {
                                                open_seq = Some(c.representative.clone());
                                            }
                                        } else {
                                            ui.label("-");
                                        }
                                        if c.member_count > 0
                                            && ui
                                                .button("Gel SVG")
                                                .on_hover_text(
                                                    "Export a serial gel lane for this container",
                                                )
                                                .clicked()
                                        {
                                            export_container_gel = Some((
                                                c.container_id.clone(),
                                                vec![c.container_id.clone()],
                                            ));
                                        }
                                    });
                                    ui.end_row();
                                }
                            });
                        ui.separator();
                        ui.heading("Arrangements");
                        ui.label("Serial lane setups from one or more containers");
                        if self.lineage_arrangements.is_empty() {
                            ui.label("No arrangements recorded");
                        } else {
                            egui::Grid::new("arrangement_grid")
                                .striped(true)
                                .show(ui, |ui| {
                                    ui.strong("Arrangement");
                                    ui.strong("Mode");
                                    ui.strong("Name");
                                    ui.strong("Lanes");
                                    ui.strong("Lane containers");
                                    ui.strong("Ladders");
                                    ui.strong("Action");
                                    ui.end_row();
                                    for arrangement in &self.lineage_arrangements {
                                        ui.monospace(&arrangement.arrangement_id);
                                        ui.label(&arrangement.mode);
                                        ui.label(if arrangement.name.trim().is_empty() {
                                            "-".to_string()
                                        } else {
                                            arrangement.name.clone()
                                        });
                                        ui.monospace(arrangement.lane_count.to_string());
                                        ui.label(if arrangement.lane_container_ids.is_empty() {
                                            "-".to_string()
                                        } else {
                                            arrangement.lane_container_ids.join(", ")
                                        });
                                        ui.label(if arrangement.ladders.is_empty() {
                                            "auto".to_string()
                                        } else {
                                            arrangement.ladders.join(", ")
                                        });
                                        ui.horizontal(|ui| {
                                            if ui
                                                .button("Export Gel")
                                                .on_hover_text(
                                                    "Export one serial gel using this arrangement",
                                                )
                                                .clicked()
                                            {
                                                let stem = if arrangement.name.trim().is_empty() {
                                                    arrangement.arrangement_id.clone()
                                                } else {
                                                    arrangement.name.clone()
                                                };
                                                export_arrangement_gel = Some((
                                                    stem,
                                                    arrangement.arrangement_id.clone(),
                                                ));
                                            }
                                            if !arrangement.lane_container_ids.is_empty()
                                                && ui
                                                    .button("Open Lanes")
                                                    .on_hover_text(
                                                        "Open windows for the lane containers in this arrangement",
                                                    )
                                                    .clicked()
                                            {
                                                open_lane_containers =
                                                    Some(arrangement.lane_container_ids.clone());
                                            }
                                        });
                                        ui.end_row();
                                    }
                                });
                        }
                    });
            });

        if let Some((stem, container_ids)) = export_container_gel.take() {
            self.prompt_export_serial_gel_svg(&stem, Some(container_ids), None, None);
        }
        if let Some((stem, arrangement_id)) = export_arrangement_gel.take() {
            self.prompt_export_serial_gel_svg(&stem, None, Some(arrangement_id), None);
        }
        if let Some(container_ids) = open_lane_containers.take() {
            for container_id in &container_ids {
                if let Some(container_row) = self
                    .lineage_containers
                    .iter()
                    .find(|row| row.container_id == *container_id)
                    .cloned()
                {
                    if container_row.member_count > 1 {
                        self.open_pool_window(
                            &container_row.representative,
                            container_row.members.clone(),
                        );
                    } else if !container_row.representative.is_empty() {
                        self.open_sequence_window(&container_row.representative);
                    }
                }
            }
        }

        if (self.lineage_graph_zoom - graph_zoom).abs() > 0.0001 {
            self.lineage_graph_zoom = graph_zoom;
            persist_workspace_after_frame = true;
        }
        if (self.lineage_graph_area_height - graph_area_height).abs() > 0.5 {
            self.lineage_graph_area_height = graph_area_height;
            persist_workspace_after_frame = true;
        }
        if (self.lineage_container_area_height - container_area_height).abs() > 0.5 {
            self.lineage_container_area_height = container_area_height;
            persist_workspace_after_frame = true;
        }
        if self.lineage_graph_compact_labels != graph_compact_labels {
            self.lineage_graph_compact_labels = graph_compact_labels;
            persist_workspace_after_frame = true;
        }
        if (self.lineage_graph_scroll_offset - graph_scroll_offset).length_sq() > 0.25 {
            self.lineage_graph_scroll_offset = graph_scroll_offset;
            persist_workspace_after_frame = true;
        }
        if persist_workspace_after_frame {
            self.persist_lineage_graph_workspace_to_state();
        }

        if let Some(input) = select_candidate_from {
            let criterion = format!("gui_lineage_select:{input}");
            let result = self
                .engine
                .write()
                .unwrap()
                .apply(Operation::SelectCandidate {
                    input: input.clone(),
                    criterion,
                    output_id: None,
                });
            if let Ok(op_result) = result {
                if let Some(seq_id) = op_result.created_seq_ids.first() {
                    open_seq = Some(seq_id.clone());
                }
            }
        }

        if let Some(seq_id) = open_seq {
            self.open_sequence_window(&seq_id);
        }
        if let Some((representative, pool_members)) = open_pool {
            self.open_pool_window(&representative, pool_members);
        }
    }

    fn apply_configuration_external_apps(&mut self) {
        self.configuration_rnapkin_executable =
            self.configuration_rnapkin_executable.trim().to_string();
        self.configuration_makeblastdb_executable =
            self.configuration_makeblastdb_executable.trim().to_string();
        self.configuration_blastn_executable =
            self.configuration_blastn_executable.trim().to_string();
        self.configuration_bigwig_to_bedgraph_executable = self
            .configuration_bigwig_to_bedgraph_executable
            .trim()
            .to_string();

        self.sync_runtime_tool_overrides_from_configuration();

        let rnapkin_status =
            tool_overrides::active_resolution_label("GENTLE_RNAPKIN_BIN", "rnapkin");
        let makeblastdb_status =
            tool_overrides::active_resolution_label(MAKEBLASTDB_ENV_BIN, DEFAULT_MAKEBLASTDB_BIN);
        let blastn_status =
            tool_overrides::active_resolution_label(BLASTN_ENV_BIN, DEFAULT_BLASTN_BIN);
        let bigwig_status = tool_overrides::active_resolution_label(
            BIGWIG_TO_BEDGRAPH_ENV_BIN,
            DEFAULT_BIGWIG_TO_BEDGRAPH_BIN,
        );
        self.configuration_status = format!(
            "External app settings applied (rnapkin: {}, makeblastdb: {}, blastn: {}, bigWigToBedGraph: {})",
            rnapkin_status, makeblastdb_status, blastn_status, bigwig_status
        );

        self.validate_rnapkin_executable();
        self.validate_blast_executables();
        match self.write_persisted_configuration_to_disk() {
            Ok(()) => {
                self.configuration_status
                    .push_str(" | persisted to app settings");
            }
            Err(e) => {
                self.configuration_status
                    .push_str(&format!(" | persistence failed: {e}"));
            }
        }
    }

    fn apply_configuration_graphics(&mut self) {
        self.apply_configuration_graphics_to_engine_state();
        self.configuration_status =
            "Graphics settings applied globally to the project display state".to_string();
        match self.write_persisted_configuration_to_disk() {
            Ok(()) => {
                self.configuration_status
                    .push_str(" | persisted to app settings");
            }
            Err(e) => {
                self.configuration_status
                    .push_str(&format!(" | persistence failed: {e}"));
            }
        }
    }

    fn reset_configuration_graphics_to_defaults(&mut self) {
        let defaults = DisplaySettings::default();
        self.configuration_graphics.show_sequence_panel = defaults.show_sequence_panel;
        self.configuration_graphics.show_map_panel = defaults.show_map_panel;
        self.configuration_graphics.show_features = defaults.show_features;
        self.configuration_graphics.show_cds_features = defaults.show_cds_features;
        self.configuration_graphics.show_gene_features = defaults.show_gene_features;
        self.configuration_graphics.show_mrna_features = defaults.show_mrna_features;
        self.configuration_graphics.show_tfbs = defaults.show_tfbs;
        self.configuration_graphics.regulatory_tracks_near_baseline =
            defaults.regulatory_tracks_near_baseline;
        self.configuration_graphics
            .regulatory_feature_max_view_span_bp = defaults.regulatory_feature_max_view_span_bp;
        self.configuration_graphics.tfbs_display_use_llr_bits = defaults.tfbs_display_use_llr_bits;
        self.configuration_graphics.tfbs_display_min_llr_bits = defaults.tfbs_display_min_llr_bits;
        self.configuration_graphics.tfbs_display_use_llr_quantile =
            defaults.tfbs_display_use_llr_quantile;
        self.configuration_graphics.tfbs_display_min_llr_quantile =
            defaults.tfbs_display_min_llr_quantile;
        self.configuration_graphics
            .tfbs_display_use_true_log_odds_bits = defaults.tfbs_display_use_true_log_odds_bits;
        self.configuration_graphics
            .tfbs_display_min_true_log_odds_bits = defaults.tfbs_display_min_true_log_odds_bits;
        self.configuration_graphics
            .tfbs_display_use_true_log_odds_quantile =
            defaults.tfbs_display_use_true_log_odds_quantile;
        self.configuration_graphics
            .tfbs_display_min_true_log_odds_quantile =
            defaults.tfbs_display_min_true_log_odds_quantile;
        self.configuration_graphics.vcf_display_show_snp = defaults.vcf_display_show_snp;
        self.configuration_graphics.vcf_display_show_ins = defaults.vcf_display_show_ins;
        self.configuration_graphics.vcf_display_show_del = defaults.vcf_display_show_del;
        self.configuration_graphics.vcf_display_show_sv = defaults.vcf_display_show_sv;
        self.configuration_graphics.vcf_display_show_other = defaults.vcf_display_show_other;
        self.configuration_graphics.vcf_display_pass_only = defaults.vcf_display_pass_only;
        self.configuration_graphics.vcf_display_use_min_qual = defaults.vcf_display_use_min_qual;
        self.configuration_graphics.vcf_display_min_qual = defaults.vcf_display_min_qual;
        self.configuration_graphics.vcf_display_use_max_qual = defaults.vcf_display_use_max_qual;
        self.configuration_graphics.vcf_display_max_qual = defaults.vcf_display_max_qual;
        self.configuration_graphics.vcf_display_required_info_keys =
            defaults.vcf_display_required_info_keys.clone();
        self.configuration_graphics.show_restriction_enzymes = defaults.show_restriction_enzymes;
        self.configuration_graphics.show_gc_contents = defaults.show_gc_contents;
        self.configuration_graphics.show_open_reading_frames = defaults.show_open_reading_frames;
        self.configuration_graphics.show_methylation_sites = defaults.show_methylation_sites;
        self.configuration_graphics.feature_details_font_size = defaults.feature_details_font_size;
        self.configuration_graphics.linear_view_start_bp = defaults.linear_view_start_bp;
        self.configuration_graphics.linear_view_span_bp = defaults.linear_view_span_bp;
        self.configuration_graphics_dirty = true;
    }

    fn apply_configuration_window_backdrops(&mut self) {
        self.window_backdrops = self.configuration_window_backdrops.clone();
        window_backdrop::set_window_backdrop_settings(self.window_backdrops.clone());
        self.configuration_window_backdrops_dirty = false;
        self.configuration_status = "Window styling settings applied".to_string();
        match self.write_persisted_configuration_to_disk() {
            Ok(()) => {
                self.configuration_status
                    .push_str(" | persisted to app settings");
            }
            Err(e) => {
                self.configuration_status
                    .push_str(&format!(" | persistence failed: {e}"));
            }
        }
    }

    fn reset_configuration_window_backdrops_to_defaults(&mut self) {
        self.configuration_window_backdrops = WindowBackdropSettings::default();
        self.configuration_window_backdrops_dirty = true;
    }

    fn render_window_backdrop_path_row(
        ui: &mut Ui,
        label: &str,
        value: &mut String,
        changed: &mut bool,
    ) {
        ui.vertical(|ui| {
            ui.horizontal(|ui| {
                ui.label(label);
                if ui
                    .add(
                        egui::TextEdit::singleline(value)
                            .desired_width(360.0)
                            .hint_text("/absolute/path/to/image.png"),
                    )
                    .on_hover_text(
                        "Optional absolute or working-directory-relative image path for this window type",
                    )
                    .changed()
                {
                    *changed = true;
                }
                if ui
                    .button("Browse...")
                    .on_hover_text("Pick an image file for this window backdrop")
                    .clicked()
                {
                    if let Some(path) = rfd::FileDialog::new()
                        .add_filter("Images", &["png", "jpg", "jpeg", "gif", "bmp", "webp"])
                        .pick_file()
                    {
                        *value = path.display().to_string();
                        *changed = true;
                    }
                }
                if ui
                    .button("Clear")
                    .on_hover_text("Clear custom image path for this window type")
                    .clicked()
                {
                    value.clear();
                    *changed = true;
                }
            });

            let trimmed = value.trim();
            if trimmed.is_empty() {
                ui.small("No image path configured (text watermark/tint fallback only)");
            } else {
                match window_backdrop::validate_window_backdrop_image_path(trimmed) {
                    Ok(resolved) => {
                        ui.colored_label(
                            egui::Color32::from_rgb(20, 140, 45),
                            format!("Resolved image: {resolved}"),
                        );
                    }
                    Err(message) => {
                        ui.colored_label(
                            egui::Color32::from_rgb(180, 50, 50),
                            format!("Path check failed: {message}"),
                        );
                    }
                }
            }
        });
    }

    fn refresh_open_sequence_windows(&mut self, ctx: &egui::Context) -> usize {
        let mut refreshed = 0usize;
        for window in self.windows.values() {
            if let Ok(mut guard) = window.write() {
                guard.refresh_from_engine_settings();
                refreshed += 1;
            }
        }
        ctx.request_repaint();
        refreshed
    }

    fn render_configuration_external_tab(&mut self, ui: &mut Ui) {
        ui.label("Configure external application integration used by shared engine features.");
        ui.separator();
        ui.label("rnapkin executable override");
        let rnapkin_edit_response = ui.add(
            egui::TextEdit::singleline(&mut self.configuration_rnapkin_executable)
                .hint_text("Leave empty to use PATH lookup for 'rnapkin'"),
        );
        if rnapkin_edit_response.changed() {
            self.clear_rnapkin_validation();
        }
        let active_rnapkin =
            tool_overrides::active_resolution_label("GENTLE_RNAPKIN_BIN", "rnapkin");
        ui.monospace(format!("Active resolution: {active_rnapkin}"));

        ui.separator();
        ui.label("makeblastdb executable override");
        let makeblastdb_edit_response = ui.add(
            egui::TextEdit::singleline(&mut self.configuration_makeblastdb_executable).hint_text(
                format!(
                    "Leave empty to use PATH lookup for '{}'",
                    DEFAULT_MAKEBLASTDB_BIN
                ),
            ),
        );
        if makeblastdb_edit_response.changed() {
            self.clear_blast_validation();
        }
        let active_makeblastdb =
            tool_overrides::active_resolution_label(MAKEBLASTDB_ENV_BIN, DEFAULT_MAKEBLASTDB_BIN);
        ui.monospace(format!("Active makeblastdb: {active_makeblastdb}"));

        ui.label("blastn executable override");
        let blastn_edit_response = ui.add(
            egui::TextEdit::singleline(&mut self.configuration_blastn_executable).hint_text(
                format!(
                    "Leave empty to use PATH lookup for '{}'",
                    DEFAULT_BLASTN_BIN
                ),
            ),
        );
        if blastn_edit_response.changed() {
            self.clear_blast_validation();
        }
        let active_blastn =
            tool_overrides::active_resolution_label(BLASTN_ENV_BIN, DEFAULT_BLASTN_BIN);
        ui.monospace(format!("Active blastn: {active_blastn}"));

        ui.separator();
        ui.label("bigWigToBedGraph executable override");
        ui.add(
            egui::TextEdit::singleline(&mut self.configuration_bigwig_to_bedgraph_executable)
                .hint_text(format!(
                    "Leave empty to use PATH lookup for '{}'",
                    DEFAULT_BIGWIG_TO_BEDGRAPH_BIN
                )),
        );
        let active_bigwig = tool_overrides::active_resolution_label(
            BIGWIG_TO_BEDGRAPH_ENV_BIN,
            DEFAULT_BIGWIG_TO_BEDGRAPH_BIN,
        );
        ui.monospace(format!("Active bigWigToBedGraph: {active_bigwig}"));

        ui.horizontal(|ui| {
            if ui
                .button("Use PATH")
                .on_hover_text("Clear rnapkin override and use PATH lookup")
                .clicked()
            {
                self.configuration_rnapkin_executable.clear();
                self.clear_rnapkin_validation();
            }
            if ui
                .button("Use PATH (BLAST)")
                .on_hover_text("Clear makeblastdb/blastn overrides and use PATH lookup")
                .clicked()
            {
                self.configuration_makeblastdb_executable.clear();
                self.configuration_blastn_executable.clear();
                self.clear_blast_validation();
            }
            if ui
                .button("Use PATH (BigWig)")
                .on_hover_text("Clear bigWigToBedGraph override and use PATH lookup")
                .clicked()
            {
                self.configuration_bigwig_to_bedgraph_executable.clear();
            }
            if ui
                .button("Validate rnapkin")
                .on_hover_text("Run rnapkin --version and capture validation status")
                .clicked()
            {
                self.validate_rnapkin_executable();
            }
            if ui
                .button("Validate BLAST tools")
                .on_hover_text("Run makeblastdb/blastn --version checks")
                .clicked()
            {
                self.validate_blast_executables();
            }
            if ui
                .button("Apply External Settings")
                .on_hover_text("Apply executable overrides to current runtime environment")
                .clicked()
            {
                self.apply_configuration_external_apps();
            }
        });
        if let Some(ok) = self.configuration_rnapkin_validation_ok {
            let color = if ok {
                egui::Color32::from_rgb(20, 140, 45)
            } else {
                egui::Color32::from_rgb(180, 50, 50)
            };
            ui.colored_label(color, self.configuration_rnapkin_validation_message.clone());
        } else if !self
            .configuration_rnapkin_validation_message
            .trim()
            .is_empty()
        {
            ui.monospace(self.configuration_rnapkin_validation_message.clone());
        }

        if let Some(ok) = self.configuration_blast_validation_ok {
            let color = if ok {
                egui::Color32::from_rgb(20, 140, 45)
            } else {
                egui::Color32::from_rgb(180, 50, 50)
            };
            ui.colored_label(color, self.configuration_blast_validation_message.clone());
        } else if !self
            .configuration_blast_validation_message
            .trim()
            .is_empty()
        {
            ui.monospace(self.configuration_blast_validation_message.clone());
        }
    }

    fn render_configuration_graphics_tab(&mut self, ui: &mut Ui) {
        ui.label("Configure project-level graphics visibility defaults.");
        ui.separator();
        let mut changed = false;
        let mut backdrop_changed = false;

        ui.heading("Panels");
        changed |= ui
            .checkbox(
                &mut self.configuration_graphics.show_sequence_panel,
                "Show sequence panel",
            )
            .changed();
        changed |= ui
            .checkbox(
                &mut self.configuration_graphics.show_map_panel,
                "Show map panel",
            )
            .changed();

        ui.separator();
        ui.heading("Feature Layers");
        changed |= ui
            .checkbox(
                &mut self.configuration_graphics.show_features,
                "Show feature overlays",
            )
            .changed();
        changed |= ui
            .checkbox(
                &mut self.configuration_graphics.show_cds_features,
                "Show CDS features",
            )
            .changed();
        changed |= ui
            .checkbox(
                &mut self.configuration_graphics.show_gene_features,
                "Show gene features",
            )
            .changed();
        changed |= ui
            .checkbox(
                &mut self.configuration_graphics.show_mrna_features,
                "Show mRNA features",
            )
            .changed();
        changed |= ui
            .checkbox(
                &mut self.configuration_graphics.show_tfbs,
                "Show TFBS features",
            )
            .changed();
        changed |= ui
            .checkbox(
                &mut self.configuration_graphics.regulatory_tracks_near_baseline,
                "Place regulatory features near DNA/GC strip",
            )
            .changed();
        ui.horizontal(|ui| {
            ui.label("Regulatory max view span");
            if ui
                .add(
                    egui::DragValue::new(
                        &mut self.configuration_graphics.regulatory_feature_max_view_span_bp,
                    )
                    .range(0..=5_000_000)
                    .speed(100.0)
                    .suffix(" bp"),
                )
                .on_hover_text(
                    "Regulatory features are hidden in linear view when current view span exceeds this threshold",
                )
                .changed()
            {
                changed = true;
            }
        });

        ui.separator();
        ui.heading("Overlays");
        changed |= ui
            .checkbox(
                &mut self.configuration_graphics.show_restriction_enzymes,
                "Show restriction enzymes",
            )
            .changed();
        changed |= ui
            .checkbox(
                &mut self.configuration_graphics.show_gc_contents,
                "Show GC contents",
            )
            .changed();
        changed |= ui
            .checkbox(
                &mut self.configuration_graphics.show_open_reading_frames,
                "Show ORFs",
            )
            .changed();
        changed |= ui
            .checkbox(
                &mut self.configuration_graphics.show_methylation_sites,
                "Show methylation sites",
            )
            .changed();

        ui.separator();
        ui.heading("Window Styling (experimental)");
        backdrop_changed |= ui
            .checkbox(
                &mut self.configuration_window_backdrops.enabled,
                "Enable themed window backdrops",
            )
            .on_hover_text("Apply subtle per-window-type color/image styling")
            .changed();
        backdrop_changed |= ui
            .checkbox(
                &mut self.configuration_window_backdrops.draw_images,
                "Use background images",
            )
            .on_hover_text("Render optional image watermark per window type")
            .changed();
        backdrop_changed |= ui
            .checkbox(
                &mut self.configuration_window_backdrops.show_text_watermark,
                "Show text watermark fallback",
            )
            .on_hover_text("Show low-contrast type text when no image is configured")
            .changed();
        ui.horizontal(|ui| {
            ui.label("Tint opacity");
            if ui
                .add(
                    egui::Slider::new(
                        &mut self.configuration_window_backdrops.tint_opacity,
                        0.0..=0.25,
                    )
                    .show_value(true),
                )
                .on_hover_text("Background color intensity behind UI content")
                .changed()
            {
                backdrop_changed = true;
            }
        });
        ui.horizontal(|ui| {
            ui.label("Image opacity");
            if ui
                .add(
                    egui::Slider::new(
                        &mut self.configuration_window_backdrops.image_opacity,
                        0.0..=0.25,
                    )
                    .show_value(true),
                )
                .on_hover_text("Image watermark intensity")
                .changed()
            {
                backdrop_changed = true;
            }
        });
        Self::render_window_backdrop_path_row(
            ui,
            "Main window image",
            &mut self.configuration_window_backdrops.main_image_path,
            &mut backdrop_changed,
        );
        Self::render_window_backdrop_path_row(
            ui,
            "Sequence window image",
            &mut self.configuration_window_backdrops.sequence_image_path,
            &mut backdrop_changed,
        );
        Self::render_window_backdrop_path_row(
            ui,
            "Pool window image",
            &mut self.configuration_window_backdrops.pool_image_path,
            &mut backdrop_changed,
        );
        Self::render_window_backdrop_path_row(
            ui,
            "Configuration window image",
            &mut self.configuration_window_backdrops.configuration_image_path,
            &mut backdrop_changed,
        );
        Self::render_window_backdrop_path_row(
            ui,
            "Help window image",
            &mut self.configuration_window_backdrops.help_image_path,
            &mut backdrop_changed,
        );

        if changed {
            self.configuration_graphics_dirty = true;
        }
        if backdrop_changed {
            self.configuration_window_backdrops_dirty = true;
        }

        ui.horizontal(|ui| {
            if ui
                .button("Reset Defaults")
                .on_hover_text("Reset graphics settings to built-in defaults")
                .clicked()
            {
                self.reset_configuration_graphics_to_defaults();
            }
            if ui
                .add_enabled(
                    self.configuration_graphics_dirty,
                    egui::Button::new("Apply Graphics Settings"),
                )
                .on_hover_text("Apply graphics settings to project state")
                .clicked()
            {
                self.apply_configuration_graphics();
            }
            if ui
                .button("Apply + Refresh Open Windows")
                .on_hover_text("Apply settings and refresh currently open sequence windows")
                .clicked()
            {
                self.apply_configuration_graphics();
                let refreshed = self.refresh_open_sequence_windows(ui.ctx());
                self.configuration_status.push_str(&format!(
                    " | refreshed {} open sequence window(s)",
                    refreshed
                ));
            }
            if ui
                .button("Reset Window Styling")
                .on_hover_text("Reset themed backdrop settings to defaults")
                .clicked()
            {
                self.reset_configuration_window_backdrops_to_defaults();
            }
            if ui
                .add_enabled(
                    self.configuration_window_backdrops_dirty,
                    egui::Button::new("Apply Window Styling"),
                )
                .on_hover_text("Apply themed backdrop settings to all windows")
                .clicked()
            {
                self.apply_configuration_window_backdrops();
            }
            if ui
                .button("Reload Current")
                .on_hover_text("Reload graphics settings from current project state")
                .clicked()
            {
                self.sync_configuration_from_runtime();
                self.configuration_status =
                    "Reloaded graphics settings from current project".to_string();
            }
        });
    }

    fn render_configuration_contents(&mut self, ui: &mut Ui) {
        window_backdrop::paint_window_backdrop(
            ui,
            WindowBackdropKind::Configuration,
            &self.window_backdrops,
        );
        self.render_specialist_window_nav(ui);
        ui.horizontal(|ui| {
            if ui
                .selectable_label(
                    self.configuration_tab == ConfigurationTab::ExternalApplications,
                    "External Applications",
                )
                .clicked()
            {
                self.configuration_tab = ConfigurationTab::ExternalApplications;
            }
            if ui
                .selectable_label(
                    self.configuration_tab == ConfigurationTab::Graphics,
                    "Graphics",
                )
                .clicked()
            {
                self.configuration_tab = ConfigurationTab::Graphics;
            }
            ui.separator();
            if ui
                .button("Close")
                .on_hover_text("Close configuration dialog")
                .clicked()
            {
                self.show_configuration_dialog = false;
            }
        });
        ui.separator();
        egui::ScrollArea::vertical()
            .auto_shrink([false, false])
            .show(ui, |ui| match self.configuration_tab {
                ConfigurationTab::ExternalApplications => {
                    self.render_configuration_external_tab(ui);
                }
                ConfigurationTab::Graphics => {
                    self.render_configuration_graphics_tab(ui);
                }
            });
        if !self.configuration_status.trim().is_empty() {
            ui.separator();
            ui.monospace(self.configuration_status.clone());
        }
    }

    fn render_configuration_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_configuration_dialog {
            return;
        }
        let builder = egui::ViewportBuilder::default()
            .with_title("Configuration")
            .with_inner_size([720.0, 540.0])
            .with_min_inner_size([460.0, 320.0]);
        ctx.show_viewport_immediate(Self::configuration_viewport_id(), builder, |ctx, class| {
            self.note_viewport_focus_if_active(ctx, Self::configuration_viewport_id());
            if class == egui::ViewportClass::Embedded {
                let mut open = self.show_configuration_dialog;
                egui::Window::new("Configuration")
                    .open(&mut open)
                    .resizable(true)
                    .default_size(Vec2::new(720.0, 540.0))
                    .show(ctx, |ui| self.render_configuration_contents(ui));
                self.show_configuration_dialog = open;
                return;
            }

            egui::CentralPanel::default().show(ctx, |ui| {
                self.render_configuration_contents(ui);
            });

            if ctx.input(|i| i.viewport().close_requested()) {
                self.show_configuration_dialog = false;
            }
        });
    }

    fn render_unsaved_changes_dialog(&mut self, ctx: &egui::Context) {
        if self.pending_project_action.is_none() {
            return;
        }
        egui::Window::new("Unsaved Changes")
            .collapsible(false)
            .resizable(false)
            .show(ctx, |ui| {
                ui.label("Save changes to the current project before continuing?");
                ui.horizontal(|ui| {
                    if ui
                        .button("Save")
                        .on_hover_text("Save current project, then continue")
                        .clicked()
                    {
                        if self.save_current_project() {
                            if let Some(action) = self.pending_project_action.take() {
                                self.execute_project_action(action);
                            }
                        }
                    }
                    if ui
                        .button("Don't Save")
                        .on_hover_text("Continue without saving current project changes")
                        .clicked()
                    {
                        if let Some(action) = self.pending_project_action.take() {
                            self.execute_project_action(action);
                        }
                    }
                    if ui
                        .button("Cancel")
                        .on_hover_text("Cancel and keep editing current project")
                        .clicked()
                    {
                        self.pending_project_action = None;
                    }
                });
            });
    }

    fn render_about_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_about_dialog {
            return;
        }
        egui::Window::new("About GENtle")
            .open(&mut self.show_about_dialog)
            .collapsible(false)
            .resizable(false)
            .show(ctx, |ui| {
                ui.vertical_centered(|ui| {
                    ui.add(APP_ICON.clone().fit_to_exact_size(Vec2::new(96.0, 96.0)));
                    for line in about::version_cli_text().lines() {
                        if line.starts_with("GENtle ") {
                            ui.heading(line);
                        } else {
                            ui.label(line);
                        }
                    }
                });
            });
    }

    fn render_command_palette_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_command_palette_dialog {
            return;
        }

        let mut entries = self.collect_command_palette_entries();
        let query = self.command_palette_query.trim().to_ascii_lowercase();
        if !query.is_empty() {
            entries.retain(|entry| {
                entry.title.to_ascii_lowercase().contains(&query)
                    || entry.detail.to_ascii_lowercase().contains(&query)
                    || entry.keywords.to_ascii_lowercase().contains(&query)
            });
        }
        if entries.is_empty() {
            self.command_palette_selected = 0;
        } else if self.command_palette_selected >= entries.len() {
            self.command_palette_selected = entries.len() - 1;
        }

        let mut open = self.show_command_palette_dialog;
        let mut execute_action: Option<CommandPaletteAction> = None;
        let builder = egui::ViewportBuilder::default()
            .with_title("Command Palette")
            .with_inner_size([760.0, 520.0])
            .with_min_inner_size([500.0, 320.0]);
        ctx.show_viewport_immediate(
            Self::command_palette_viewport_id(),
            builder,
            |ctx, class| {
                self.note_viewport_focus_if_active(ctx, Self::command_palette_viewport_id());
                let mut render_contents = |ui: &mut Ui| {
                    ui.label("Search actions, settings, and help topics");
                    let input_id = ui.make_persistent_id("gentle_command_palette_search");
                    let search_response = ui.add(
                        egui::TextEdit::singleline(&mut self.command_palette_query)
                            .id(input_id)
                            .desired_width(f32::INFINITY)
                            .hint_text("Type action name (Cmd/Ctrl+K)"),
                    );
                    if self.command_palette_focus_query {
                        search_response.request_focus();
                        self.command_palette_focus_query = false;
                    }
                    if ui.input(|i| i.key_pressed(Key::ArrowDown)) {
                        if !entries.is_empty() {
                            self.command_palette_selected =
                                (self.command_palette_selected + 1) % entries.len();
                        }
                    }
                    if ui.input(|i| i.key_pressed(Key::ArrowUp)) && !entries.is_empty() {
                        if self.command_palette_selected == 0 {
                            self.command_palette_selected = entries.len() - 1;
                        } else {
                            self.command_palette_selected -= 1;
                        }
                    }
                    if ui.input(|i| i.key_pressed(Key::Enter))
                        && !entries.is_empty()
                        && self.command_palette_selected < entries.len()
                    {
                        execute_action = Some(entries[self.command_palette_selected].action);
                    }

                    ui.separator();
                    if entries.is_empty() {
                        ui.small("No matching commands");
                    } else {
                        egui::ScrollArea::vertical()
                            .max_height(400.0)
                            .show(ui, |ui| {
                                for (idx, entry) in entries.iter().enumerate() {
                                    let selected = self.command_palette_selected == idx;
                                    let label = format!("{} — {}", entry.title, entry.detail);
                                    let response = ui.selectable_label(selected, label);
                                    if response.hovered() {
                                        self.command_palette_selected = idx;
                                        self.hover_status_name =
                                            format!("Command palette: {}", entry.title);
                                    }
                                    if response.clicked() {
                                        execute_action = Some(entry.action);
                                    }
                                }
                            });
                    }
                };

                if class == egui::ViewportClass::Embedded {
                    egui::Window::new("Command Palette")
                        .open(&mut open)
                        .collapsible(false)
                        .resizable(true)
                        .default_size(Vec2::new(760.0, 520.0))
                        .show(ctx, |ui| render_contents(ui));
                } else {
                    egui::CentralPanel::default().show(ctx, |ui| {
                        render_contents(ui);
                    });
                    if ctx.input(|i| i.viewport().close_requested()) {
                        open = false;
                    }
                }
            },
        );

        if let Some(action) = execute_action {
            self.execute_command_palette_action(ctx, action);
            open = false;
        }

        if ctx.input(|i| i.key_pressed(Key::Escape)) {
            open = false;
        }

        self.show_command_palette_dialog = open;
    }

    fn render_jobs_panel(&mut self, ctx: &egui::Context) {
        if !self.show_jobs_panel {
            return;
        }
        let mut open = self.show_jobs_panel;
        egui::Window::new("Background Jobs")
            .open(&mut open)
            .collapsible(false)
            .resizable(true)
            .default_size(Vec2::new(760.0, 480.0))
            .show(ctx, |ui| {
                ui.label("Centralized progress, cancellation, and completion summaries");
                ui.separator();

                ui.strong("Prepare Genome");
                let mut cancel_prepare_clicked = false;
                if self.genome_prepare_task.is_some() {
                    ui.horizontal(|ui| {
                        ui.add(egui::Spinner::new());
                        if let Some(progress) = &self.genome_prepare_progress {
                            ui.label(format!(
                                "{}: {} ({})",
                                progress.genome_id, progress.phase, progress.item
                            ));
                        } else {
                            ui.label("Running...");
                        }
                        if ui
                            .button("Cancel")
                            .on_hover_text("Request cancellation of genome prepare job")
                            .clicked()
                        {
                            cancel_prepare_clicked = true;
                        }
                    });
                } else {
                    ui.horizontal(|ui| {
                        ui.small("Idle");
                        if ui
                            .button("Retry")
                            .on_hover_text("Run prepare genome again using current dialog settings")
                            .clicked()
                        {
                            self.push_job_event(
                                BackgroundJobKind::PrepareGenome,
                                BackgroundJobEventPhase::Retried,
                                None,
                                "Retry requested from background jobs panel",
                            );
                            self.start_prepare_reference_genome();
                        }
                    });
                }
                if cancel_prepare_clicked {
                    self.request_prepare_task_cancel("background jobs panel");
                }
                if !self.genome_prepare_status.trim().is_empty() {
                    ui.small(self.genome_prepare_status.clone());
                }

                ui.separator();
                ui.strong("BLAST");
                if let Some(task) = &self.genome_blast_task {
                    ui.horizontal(|ui| {
                        ui.add(egui::Spinner::new());
                        ui.label(format!("Running ({:.1}s)", task.started.elapsed().as_secs_f32()));
                        ui.small("Cancellation is not yet available for BLAST jobs");
                    });
                    if let Some(fraction) = self.genome_blast_progress_fraction {
                        ui.add(
                            egui::ProgressBar::new(fraction.clamp(0.0, 1.0))
                                .show_percentage()
                                .text(self.genome_blast_progress_label.clone()),
                            );
                    }
                } else {
                    ui.horizontal(|ui| {
                        ui.small("Idle");
                        if ui
                            .button("Retry")
                            .on_hover_text("Run BLAST again using current BLAST dialog settings")
                            .clicked()
                        {
                            self.push_job_event(
                                BackgroundJobKind::BlastGenome,
                                BackgroundJobEventPhase::Retried,
                                None,
                                "Retry requested from background jobs panel",
                            );
                            self.start_reference_genome_blast();
                        }
                    });
                }
                if !self.genome_blast_status.trim().is_empty() {
                    ui.small(self.genome_blast_status.clone());
                }

                ui.separator();
                ui.strong("Genome Track Import");
                let mut cancel_track_import_clicked = false;
                if self.genome_track_import_task.is_some() {
                    ui.horizontal(|ui| {
                        ui.add(egui::Spinner::new());
                        if let Some(progress) = &self.genome_track_import_progress {
                            ui.label(format!(
                                "{} '{}' parsed={} imported={} skipped={}",
                                progress.source,
                                progress.seq_id,
                                progress.parsed_records,
                                progress.imported_features,
                                progress.skipped_records
                            ));
                        } else {
                            ui.label("Running...");
                        }
                        if ui
                            .button("Cancel")
                            .on_hover_text("Request cancellation of running track-import job")
                            .clicked()
                        {
                            cancel_track_import_clicked = true;
                        }
                    });
                } else {
                    ui.horizontal(|ui| {
                        ui.small("Idle");
                        if ui
                            .button("Retry")
                            .on_hover_text(
                                "Run track import again for the currently selected anchored sequence",
                            )
                            .clicked()
                        {
                            self.push_job_event(
                                BackgroundJobKind::TrackImport,
                                BackgroundJobEventPhase::Retried,
                                None,
                                "Retry requested from background jobs panel",
                            );
                            self.import_genome_bed_track_for_selected_sequence();
                        }
                    });
                }
                if cancel_track_import_clicked {
                    self.request_track_import_task_cancel("background jobs panel");
                }
                if !self.genome_track_status.trim().is_empty() {
                    ui.small(self.genome_track_status.clone());
                }

                ui.separator();
                ui.strong("Agent Assistant");
                if let Some(task) = &self.agent_task {
                    ui.horizontal(|ui| {
                        ui.add(egui::Spinner::new());
                        ui.label(format!("Running ({:.1}s)", task.started.elapsed().as_secs_f32()));
                        ui.small("Cancellation is not yet available for agent requests");
                    });
                } else {
                    ui.horizontal(|ui| {
                        ui.small("Idle");
                        if ui
                            .button("Retry")
                            .on_hover_text("Run the agent assistant request again with current prompt/settings")
                            .clicked()
                        {
                            self.push_job_event(
                                BackgroundJobKind::AgentAssist,
                                BackgroundJobEventPhase::Retried,
                                None,
                                "Retry requested from background jobs panel",
                            );
                            self.start_agent_assistant_request();
                        }
                    });
                }
                if !self.agent_status.trim().is_empty() {
                    ui.small(self.agent_status.clone());
                }

                ui.separator();
                ui.strong("Recent job events");
                egui::ScrollArea::vertical()
                    .max_height(180.0)
                    .show(ui, |ui| {
                        for event in self.job_event_log.iter().rev().take(40) {
                            ui.small(event.to_line());
                        }
                    });
            });
        self.show_jobs_panel = open;
    }

    fn render_history_panel(&mut self, ctx: &egui::Context) {
        if !self.show_history_panel {
            return;
        }
        let mut open = self.show_history_panel;
        let (undo_count, redo_count, history_rows) = {
            let engine = self.engine.read().unwrap();
            let rows = engine
                .operation_log()
                .iter()
                .rev()
                .take(120)
                .map(|record| {
                    (
                        record.result.op_id.clone(),
                        record.run_id.clone(),
                        Self::summarize_operation(&record.op),
                    )
                })
                .collect::<Vec<_>>();
            (engine.undo_available(), engine.redo_available(), rows)
        };

        let builder = egui::ViewportBuilder::default()
            .with_title("Operation History")
            .with_inner_size([820.0, 520.0])
            .with_min_inner_size([560.0, 320.0]);
        ctx.show_viewport_immediate(Self::history_viewport_id(), builder, |ctx, class| {
            self.note_viewport_focus_if_active(ctx, Self::history_viewport_id());
            let mut render_contents = |ui: &mut Ui| {
                ui.label("Operation-level history with undo/redo");
                ui.horizontal(|ui| {
                    if ui
                        .add_enabled(undo_count > 0, egui::Button::new("Undo"))
                        .on_hover_text("Undo the most recent operation-level state transition")
                        .clicked()
                    {
                        self.undo_last_operation();
                    }
                    if ui
                        .add_enabled(redo_count > 0, egui::Button::new("Redo"))
                        .on_hover_text("Redo the most recently undone operation-level transition")
                        .clicked()
                    {
                        self.redo_last_operation();
                    }
                    ui.small(format!(
                        "undo available: {} | redo available: {}",
                        undo_count, redo_count
                    ));
                });
                ui.separator();
                if history_rows.is_empty() {
                    ui.small("No operations recorded yet.");
                } else {
                    egui::ScrollArea::vertical().show(ui, |ui| {
                        for (op_id, run_id, summary) in &history_rows {
                            ui.monospace(format!("[{op_id}] run={run_id}"));
                            ui.small(summary);
                            ui.separator();
                        }
                    });
                }
            };

            if class == egui::ViewportClass::Embedded {
                egui::Window::new("Operation History")
                    .open(&mut open)
                    .collapsible(false)
                    .resizable(true)
                    .default_size(Vec2::new(820.0, 520.0))
                    .show(ctx, |ui| render_contents(ui));
            } else {
                egui::CentralPanel::default().show(ctx, |ui| {
                    render_contents(ui);
                });
                if ctx.input(|i| i.viewport().close_requested()) {
                    open = false;
                }
            }
        });
        self.show_history_panel = open;
    }

    fn render_status_bar(&mut self, ctx: &egui::Context) {
        let (undo_count, redo_count) = {
            let engine = self.engine.read().unwrap();
            (engine.undo_available(), engine.redo_available())
        };
        egui::TopBottomPanel::bottom("gentle_status_bar").show(ctx, |ui| {
            ui.horizontal_wrapped(|ui| {
                let hover_text = if self.hover_status_name.trim().is_empty() {
                    "Hover: -".to_string()
                } else {
                    format!("Hover: {}", self.hover_status_name)
                };
                ui.monospace(hover_text);
                ui.separator();
                ui.small(format!("Undo {undo_count} / Redo {redo_count}"));
                if !self.app_status.trim().is_empty() {
                    ui.separator();
                    ui.small(self.app_status.clone());
                }
            });
        });
    }

    fn render_help_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_help_dialog {
            return;
        }
        let title = format!("Help - {}", self.active_help_title());
        let builder = egui::ViewportBuilder::default()
            .with_title(title.clone())
            .with_inner_size([860.0, 680.0])
            .with_min_inner_size([420.0, 320.0]);
        ctx.show_viewport_immediate(Self::help_viewport_id(), builder, |ctx, class| {
            self.note_viewport_focus_if_active(ctx, Self::help_viewport_id());
            if class == egui::ViewportClass::Embedded {
                let mut open = self.show_help_dialog;
                egui::Window::new(title.clone())
                    .open(&mut open)
                    .collapsible(false)
                    .resizable(true)
                    .default_size(Vec2::new(860.0, 680.0))
                    .show(ctx, |ui| {
                        self.render_help_contents(ui);
                    });
                self.show_help_dialog = open;
                return;
            }

            egui::CentralPanel::default().show(ctx, |ui| {
                self.render_help_contents(ui);
            });

            if ctx.input(|i| i.viewport().close_requested()) {
                self.show_help_dialog = false;
            }
        });
    }

    fn render_help_contents(&mut self, ui: &mut Ui) {
        window_backdrop::paint_window_backdrop(
            ui,
            WindowBackdropKind::Help,
            &self.window_backdrops,
        );
        let find_shortcut = KeyboardShortcut::new(Modifiers::COMMAND, Key::F);
        if ui.ctx().input_mut(|i| i.consume_shortcut(&find_shortcut)) {
            self.help_focus_search_box = true;
        }

        let mut active_doc_changed = false;
        ui.horizontal(|ui| {
            if ui
                .selectable_label(self.help_doc == HelpDoc::Gui, "GUI Manual")
                .clicked()
            {
                self.help_doc = HelpDoc::Gui;
                active_doc_changed = true;
            }
            if ui
                .selectable_label(self.help_doc == HelpDoc::Cli, "CLI Manual")
                .clicked()
            {
                self.help_doc = HelpDoc::Cli;
                active_doc_changed = true;
            }
            if ui
                .selectable_label(self.help_doc == HelpDoc::Shell, "Shell Commands")
                .on_hover_text("Generated from docs/glossary.json")
                .clicked()
            {
                self.help_doc = HelpDoc::Shell;
                active_doc_changed = true;
            }
            if self.help_doc == HelpDoc::Shell {
                ui.separator();
                ui.label("Interface:");
                let options = [
                    ShellHelpInterface::All,
                    ShellHelpInterface::GuiShell,
                    ShellHelpInterface::CliShell,
                    ShellHelpInterface::CliDirect,
                    ShellHelpInterface::Js,
                    ShellHelpInterface::Lua,
                ];
                let mut shell_filter_changed = false;
                for option in options {
                    if ui
                        .selectable_label(self.help_shell_interface == option, option.label())
                        .clicked()
                    {
                        self.help_shell_interface = option;
                        shell_filter_changed = true;
                    }
                }
                if shell_filter_changed {
                    self.help_shell_markdown =
                        Self::generate_shell_help_markdown_for(self.help_shell_interface);
                    self.help_markdown_cache = CommonMarkCache::default();
                    active_doc_changed = true;
                }
            }
            ui.separator();
            if ui
                .button("Reload")
                .on_hover_text("Reload help markdown files from disk")
                .clicked()
            {
                self.refresh_help_docs();
                self.help_markdown_cache = CommonMarkCache::default();
                active_doc_changed = true;
            }
            if ui
                .button("Close")
                .on_hover_text("Close help window")
                .clicked()
            {
                self.show_help_dialog = false;
            }
        });

        ui.horizontal(|ui| {
            ui.label("Find:");
            let search_id = ui.make_persistent_id("gentle_help_search_input");
            let search_response = ui.add(
                egui::TextEdit::singleline(&mut self.help_search_query)
                    .id(search_id)
                    .hint_text("Search help text (Cmd/Ctrl+F)")
                    .desired_width(260.0),
            );
            if self.help_focus_search_box {
                search_response.request_focus();
                self.help_focus_search_box = false;
            }
            if search_response.changed() {
                self.help_search_selected = 0;
                self.refresh_help_search_matches();
            }
            if search_response.has_focus() && ui.input(|i| i.key_pressed(Key::Enter)) {
                self.select_next_help_match();
            }

            let has_matches = !self.help_search_matches.is_empty();
            if ui
                .add_enabled(has_matches, egui::Button::new("Prev"))
                .on_hover_text("Jump to previous help search match")
                .clicked()
            {
                self.select_prev_help_match();
            }
            if ui
                .add_enabled(has_matches, egui::Button::new("Next"))
                .on_hover_text("Jump to next help search match")
                .clicked()
            {
                self.select_next_help_match();
            }
            if ui
                .button("Clear")
                .on_hover_text("Clear search query and match list")
                .clicked()
            {
                self.help_search_query.clear();
                self.help_search_selected = 0;
                self.help_search_matches.clear();
            }

            if self.help_search_query.trim().is_empty() {
                ui.small("No active search");
            } else if self.help_search_matches.is_empty() {
                ui.small("No matches");
            } else {
                ui.small(format!(
                    "Match {}/{}",
                    self.help_search_selected + 1,
                    self.help_search_matches.len()
                ));
            }
        });

        if active_doc_changed {
            self.refresh_help_search_matches();
        }
        if let Some(current) = self.help_search_matches.get(self.help_search_selected) {
            ui.small(format!("Line {}: {}", current.line_number, current.snippet));
        }

        ui.separator();
        let markdown = self.active_help_markdown().to_string();
        let help_images = Self::collect_help_markdown_images(markdown.as_str());
        egui::ScrollArea::vertical()
            .auto_shrink([false, false])
            .show(ui, |ui| {
                let max_image_width = (ui.available_width() * 0.75).round().clamp(220.0, 1600.0);
                CommonMarkViewer::new()
                    .max_image_width(Some(max_image_width as usize))
                    .show(ui, &mut self.help_markdown_cache, markdown.as_str());
                if !help_images.is_empty() {
                    ui.separator();
                    ui.heading("Image Captions");
                    ui.small("Click any preview to open an enlarged view.");
                    for image in &help_images {
                        let caption = Self::help_image_caption(image);
                        let response = ui.add(
                            egui::Image::new(image.path.clone())
                                .max_width(max_image_width)
                                .maintain_aspect_ratio(true)
                                .sense(egui::Sense::click()),
                        );
                        if response.clicked() {
                            self.help_image_preview_path = Some(image.path.clone());
                            self.help_image_preview_caption = caption.clone();
                        }
                        ui.small(caption);
                        ui.separator();
                    }
                }
            });

        if let Some(path) = self.help_image_preview_path.clone() {
            let title = if self.help_image_preview_caption.trim().is_empty() {
                "Help image".to_string()
            } else {
                self.help_image_preview_caption.clone()
            };
            let mut open = true;
            egui::Window::new(title)
                .open(&mut open)
                .resizable(true)
                .vscroll(true)
                .default_size(Vec2::new(960.0, 720.0))
                .show(ui.ctx(), |ui| {
                    ui.add(egui::Image::new(path.clone()).shrink_to_fit());
                    ui.horizontal(|ui| {
                        if ui.button("Close").clicked() {
                            self.help_image_preview_path = None;
                            self.help_image_preview_caption.clear();
                        }
                    });
                });
            if !open {
                self.help_image_preview_path = None;
                self.help_image_preview_caption.clear();
            }
        }
    }

    fn summarize_operation(op: &Operation) -> String {
        match op {
            Operation::LoadFile { path, as_id } => match as_id {
                Some(id) => format!("Load file: path={path}, as_id={id}"),
                None => format!("Load file: path={path}"),
            },
            Operation::DigestContainer {
                container_id,
                enzymes,
                output_prefix,
            } => format!(
                "Digest container: container_id={container_id}, enzymes=[{}], output_prefix={}",
                enzymes.join(", "),
                output_prefix.clone().unwrap_or_else(|| "-".to_string())
            ),
            Operation::MergeContainersById {
                container_ids,
                output_prefix,
            } => format!(
                "Merge containers by id: container_ids={}, output_prefix={}",
                container_ids.join(", "),
                output_prefix.clone().unwrap_or_else(|| "-".to_string())
            ),
            Operation::LigationContainer {
                container_id,
                circularize_if_possible,
                protocol,
                output_prefix,
                unique,
                ..
            } => format!(
                "Ligation container: container_id={container_id}, protocol={:?}, circularize_if_possible={}, output_prefix={}, unique={}",
                protocol,
                circularize_if_possible,
                output_prefix.clone().unwrap_or_else(|| "-".to_string()),
                unique.unwrap_or(false)
            ),
            Operation::FilterContainerByMolecularWeight {
                container_id,
                min_bp,
                max_bp,
                error,
                unique,
                ..
            } => format!(
                "Molecular weight filter (container): container_id={}, min_bp={}, max_bp={}, error={:.2}, unique={}",
                container_id, min_bp, max_bp, error, unique
            ),
            Operation::Digest {
                input,
                enzymes,
                output_prefix,
            } => {
                format!(
                    "Digest: input={input}, enzymes=[{}], output_prefix={}",
                    enzymes.join(", "),
                    output_prefix.clone().unwrap_or_else(|| "-".to_string())
                )
            }
            Operation::MergeContainers {
                inputs,
                output_prefix,
            } => format!(
                "Merge containers: inputs={}, output_prefix={}",
                inputs.join(", "),
                output_prefix.clone().unwrap_or_else(|| "-".to_string())
            ),
            Operation::Ligation {
                inputs,
                circularize_if_possible,
                protocol,
                output_prefix,
                unique,
                ..
            } => format!(
                "Ligation: inputs={}, protocol={:?}, circularize_if_possible={}, output_prefix={}, unique={}",
                inputs.join(", "),
                protocol,
                circularize_if_possible,
                output_prefix.clone().unwrap_or_else(|| "-".to_string()),
                unique.unwrap_or(false)
            ),
            Operation::Pcr {
                template,
                forward_primer,
                reverse_primer,
                unique,
                ..
            } => format!(
                "PCR: template={template}, forward_primer={}, reverse_primer={}, unique={}",
                forward_primer,
                reverse_primer,
                unique.unwrap_or(false)
            ),
            Operation::PcrAdvanced {
                template,
                forward_primer,
                reverse_primer,
                unique,
                ..
            } => format!(
                "PCR advanced: template={template}, fwd_len={}, fwd_anneal={:?}, fwd_mm={:?}, rev_len={}, rev_anneal={:?}, rev_mm={:?}, unique={}",
                forward_primer.sequence.len(),
                forward_primer.anneal_len,
                forward_primer.max_mismatches,
                reverse_primer.sequence.len(),
                reverse_primer.anneal_len,
                reverse_primer.max_mismatches,
                unique.unwrap_or(false)
            ),
            Operation::PcrMutagenesis {
                template,
                forward_primer,
                reverse_primer,
                mutations,
                require_all_mutations,
                unique,
                ..
            } => format!(
                "PCR mutagenesis: template={template}, fwd_len={}, rev_len={}, mutations={}, require_all_mutations={}, unique={}",
                forward_primer.sequence.len(),
                reverse_primer.sequence.len(),
                mutations.len(),
                require_all_mutations.unwrap_or(false),
                unique.unwrap_or(false)
            ),
            Operation::FilterByMolecularWeight {
                inputs,
                min_bp,
                max_bp,
                error,
                unique,
                ..
            } => format!(
                "Molecular weight filter: inputs={}, min_bp={}, max_bp={}, error={:.2}, unique={}",
                inputs.join(", "),
                min_bp,
                max_bp,
                error,
                unique
            ),
            Operation::FilterByDesignConstraints {
                inputs,
                gc_min,
                gc_max,
                max_homopolymer_run,
                reject_ambiguous_bases,
                avoid_u6_terminator_tttt,
                forbidden_motifs,
                unique,
                ..
            } => format!(
                "Design-constraint filter: inputs={}, gc_min={:?}, gc_max={:?}, max_homopolymer_run={:?}, reject_ambiguous_bases={}, avoid_u6_terminator_tttt={}, forbidden_motifs={}, unique={}",
                inputs.join(", "),
                gc_min,
                gc_max,
                max_homopolymer_run,
                reject_ambiguous_bases.unwrap_or(true),
                avoid_u6_terminator_tttt.unwrap_or(true),
                forbidden_motifs.len(),
                unique
            ),
            Operation::SelectCandidate {
                input,
                criterion,
                output_id,
            } => format!(
                "Select candidate: input={input}, criterion={criterion}, output_id={}",
                output_id.clone().unwrap_or_else(|| "-".to_string())
            ),
            Operation::Reverse { input, output_id } => format!(
                "Reverse: input={input}, output_id={}",
                output_id.clone().unwrap_or_else(|| "-".to_string())
            ),
            Operation::Complement { input, output_id } => format!(
                "Complement: input={input}, output_id={}",
                output_id.clone().unwrap_or_else(|| "-".to_string())
            ),
            Operation::ReverseComplement { input, output_id } => format!(
                "Reverse complement: input={input}, output_id={}",
                output_id.clone().unwrap_or_else(|| "-".to_string())
            ),
            Operation::Branch { input, output_id } => format!(
                "Branch: input={input}, output_id={}",
                output_id.clone().unwrap_or_else(|| "-".to_string())
            ),
            Operation::SetDisplayVisibility { target, visible } => {
                let target_name = match target {
                    DisplayTarget::SequencePanel => "Sequence panel",
                    DisplayTarget::MapPanel => "Map panel",
                    DisplayTarget::Features => "Features",
                    DisplayTarget::CdsFeatures => "CDS features",
                    DisplayTarget::GeneFeatures => "Gene features",
                    DisplayTarget::MrnaFeatures => "mRNA features",
                    DisplayTarget::Tfbs => "TFBS",
                    DisplayTarget::RestrictionEnzymes => "Restriction enzymes",
                    DisplayTarget::GcContents => "GC contents",
                    DisplayTarget::OpenReadingFrames => "Open reading frames",
                    DisplayTarget::MethylationSites => "Methylation sites",
                };
                format!("Set display visibility: target={target_name}, visible={visible}")
            }
            Operation::SetLinearViewport { start_bp, span_bp } => {
                format!("Set linear viewport: start_bp={start_bp}, span_bp={span_bp}")
            }
            Operation::SetTopology { seq_id, circular } => {
                if *circular {
                    format!("Set topology: seq_id={seq_id}, topology=circular")
                } else {
                    format!("Set topology: seq_id={seq_id}, topology=linear")
                }
            }
            Operation::RecomputeFeatures { seq_id } => {
                format!("Recompute features: seq_id={seq_id}")
            }
            Operation::SetParameter { name, value } => {
                format!("Set parameter: name={name}, value={value}")
            }
            Operation::AnnotateTfbs {
                seq_id,
                motifs,
                min_llr_bits,
                min_llr_quantile,
                per_tf_thresholds,
                clear_existing,
                max_hits,
            } => format!(
                "Annotate TFBS: seq_id={seq_id}, motifs=[{}], min_llr_bits={:?}, min_llr_quantile={:?}, per_tf_overrides={}, clear_existing={}, max_hits={}",
                motifs.join(", "),
                min_llr_bits,
                min_llr_quantile,
                per_tf_thresholds.len(),
                clear_existing.unwrap_or(true),
                max_hits
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "default(500)".to_string())
            ),
            Operation::ExtractRegion {
                input,
                from,
                to,
                output_id,
            } => format!(
                "Extract region: input={input}, from={from}, to={to}, output_id={}",
                output_id.clone().unwrap_or_else(|| "-".to_string())
            ),
            Operation::ExtractAnchoredRegion {
                input,
                anchor,
                direction,
                target_length_bp,
                length_tolerance_bp,
                required_re_sites,
                required_tf_motifs,
                output_prefix,
                unique,
                max_candidates,
                ..
            } => format!(
                "Extract anchored region: input={input}, anchor={anchor:?}, direction={direction:?}, target_length_bp={}, tolerance_bp={}, re_sites=[{}], tf_motifs=[{}], output_prefix={}, unique={}, max_candidates={}",
                target_length_bp,
                length_tolerance_bp,
                required_re_sites.join(", "),
                required_tf_motifs.join(", "),
                output_prefix.clone().unwrap_or_else(|| "-".to_string()),
                unique.unwrap_or(false),
                max_candidates
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string())
            ),
            Operation::SaveFile {
                seq_id,
                path,
                format,
            } => format!("Save file: seq_id={seq_id}, path={path}, format={format:?}"),
            Operation::RenderSequenceSvg { seq_id, mode, path } => {
                format!("Render sequence SVG: seq_id={seq_id}, mode={mode:?}, path={path}")
            }
            Operation::RenderRnaStructureSvg { seq_id, path } => {
                format!("Render RNA structure SVG: seq_id={seq_id}, path={path}")
            }
            Operation::RenderLineageSvg { path } => {
                format!("Render lineage SVG: path={path}")
            }
            Operation::RenderPoolGelSvg {
                inputs,
                path,
                ladders,
                container_ids,
                arrangement_id,
            } => format!(
                "Render serial gel SVG: inputs={}, container_ids={}, arrangement_id={}, path={}, ladders={}",
                inputs.join(", "),
                container_ids
                    .as_ref()
                    .map(|v| v.join(", "))
                    .filter(|s| !s.is_empty())
                    .unwrap_or_else(|| "-".to_string()),
                arrangement_id
                    .as_deref()
                    .map(str::trim)
                    .filter(|v| !v.is_empty())
                    .unwrap_or("-"),
                path,
                ladders
                    .as_ref()
                    .map(|v| v.join(", "))
                    .filter(|s| !s.is_empty())
                    .unwrap_or_else(|| "auto".to_string())
            ),
            Operation::CreateArrangementSerial {
                container_ids,
                arrangement_id,
                name,
                ladders,
            } => format!(
                "Create serial arrangement: containers={}, arrangement_id={}, name={}, ladders={}",
                container_ids.join(", "),
                arrangement_id
                    .as_deref()
                    .map(str::trim)
                    .filter(|v| !v.is_empty())
                    .unwrap_or("auto"),
                name.as_deref()
                    .map(str::trim)
                    .filter(|v| !v.is_empty())
                    .unwrap_or("-"),
                ladders
                    .as_ref()
                    .map(|v| v.join(", "))
                    .filter(|s| !s.is_empty())
                    .unwrap_or_else(|| "auto".to_string())
            ),
            Operation::ExportDnaLadders { path, name_filter } => format!(
                "Export DNA ladders: path={}, filter={}",
                path,
                name_filter
                    .as_deref()
                    .map(str::trim)
                    .filter(|v| !v.is_empty())
                    .unwrap_or("-")
            ),
            Operation::ExportRnaLadders { path, name_filter } => format!(
                "Export RNA ladders: path={}, filter={}",
                path,
                name_filter
                    .as_deref()
                    .map(str::trim)
                    .filter(|v| !v.is_empty())
                    .unwrap_or("-")
            ),
            Operation::ExportPool {
                inputs,
                path,
                pool_id,
                human_id,
            } => format!(
                "Export pool: inputs={}, path={}, pool_id={}, human_id={}",
                inputs.join(", "),
                path,
                pool_id.clone().unwrap_or_else(|| "-".to_string()),
                human_id.clone().unwrap_or_else(|| "-".to_string())
            ),
            Operation::PrepareGenome {
                genome_id,
                catalog_path,
                cache_dir,
                timeout_seconds,
            } => format!(
                "Prepare genome: genome_id={}, catalog_path={}, cache_dir={}, timeout_seconds={}",
                genome_id,
                catalog_path.clone().unwrap_or_else(|| "-".to_string()),
                cache_dir.clone().unwrap_or_else(|| "-".to_string()),
                timeout_seconds
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string())
            ),
            Operation::ExtractGenomeRegion {
                genome_id,
                chromosome,
                start_1based,
                end_1based,
                output_id,
                catalog_path,
                cache_dir,
            } => format!(
                "Extract genome region: genome_id={}, chromosome={}, start={}, end={}, output_id={}, catalog_path={}, cache_dir={}",
                genome_id,
                chromosome,
                start_1based,
                end_1based,
                output_id.clone().unwrap_or_else(|| "-".to_string()),
                catalog_path.clone().unwrap_or_else(|| "-".to_string()),
                cache_dir.clone().unwrap_or_else(|| "-".to_string())
            ),
            Operation::ExtractGenomeGene {
                genome_id,
                gene_query,
                occurrence,
                output_id,
                catalog_path,
                cache_dir,
            } => format!(
                "Extract genome gene: genome_id={}, gene_query={}, occurrence={}, output_id={}, catalog_path={}, cache_dir={}",
                genome_id,
                gene_query,
                occurrence
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                output_id.clone().unwrap_or_else(|| "-".to_string()),
                catalog_path.clone().unwrap_or_else(|| "-".to_string()),
                cache_dir.clone().unwrap_or_else(|| "-".to_string())
            ),
            Operation::ExtendGenomeAnchor {
                seq_id,
                side,
                length_bp,
                output_id,
                catalog_path,
                cache_dir,
            } => format!(
                "Extend genome anchor: seq_id={}, side={:?}, length_bp={}, output_id={}, catalog_path={}, cache_dir={}",
                seq_id,
                side,
                length_bp,
                output_id.clone().unwrap_or_else(|| "-".to_string()),
                catalog_path.clone().unwrap_or_else(|| "-".to_string()),
                cache_dir.clone().unwrap_or_else(|| "-".to_string())
            ),
            Operation::ImportGenomeBedTrack {
                seq_id,
                path,
                track_name,
                min_score,
                max_score,
                clear_existing,
            } => format!(
                "Import genome BED track: seq_id={}, path={}, track_name={}, min_score={:?}, max_score={:?}, clear_existing={}",
                seq_id,
                path,
                track_name.clone().unwrap_or_else(|| "-".to_string()),
                min_score,
                max_score,
                clear_existing.unwrap_or(false)
            ),
            Operation::ImportGenomeBigWigTrack {
                seq_id,
                path,
                track_name,
                min_score,
                max_score,
                clear_existing,
            } => format!(
                "Import genome BigWig track: seq_id={}, path={}, track_name={}, min_score={:?}, max_score={:?}, clear_existing={}",
                seq_id,
                path,
                track_name.clone().unwrap_or_else(|| "-".to_string()),
                min_score,
                max_score,
                clear_existing.unwrap_or(false)
            ),
            Operation::ImportGenomeVcfTrack {
                seq_id,
                path,
                track_name,
                min_score,
                max_score,
                clear_existing,
            } => format!(
                "Import genome VCF track: seq_id={}, path={}, track_name={}, min_score={:?}, max_score={:?}, clear_existing={}",
                seq_id,
                path,
                track_name.clone().unwrap_or_else(|| "-".to_string()),
                min_score,
                max_score,
                clear_existing.unwrap_or(false)
            ),
            Operation::ImportBlastHitsTrack {
                seq_id,
                hits,
                track_name,
                clear_existing,
            } => format!(
                "Import BLAST hits track: seq_id={}, track_name={}, hits={}, clear_existing={}",
                seq_id,
                track_name.clone().unwrap_or_else(|| "-".to_string()),
                hits.len(),
                clear_existing.unwrap_or(false)
            ),
            other => format!("{other:?}"),
        }
    }
}

impl eframe::App for GENtleApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        let update_result = catch_unwind(AssertUnwindSafe(|| {
            if !self.update_has_run_before {
                egui_extras::install_image_loaders(ctx);
                self.update_has_run_before = true;
            }
            about::install_native_help_menu_bridge();
            about::install_native_settings_menu_bridge();
            about::install_native_windows_menu_bridge();
            about::install_native_app_windows_menu_bridge();
            self.consume_native_help_request();
            self.consume_native_settings_request();
            self.consume_native_windows_request();
            self.consume_active_viewport_report();
            if ctx.input(|i| i.viewport().focused.unwrap_or(false)) {
                self.set_active_window_viewport(ViewportId::ROOT);
            }
            self.hover_status_name.clear();
            let project_dirty = self.is_project_dirty();
            let dirty_marker = if project_dirty { " *" } else { "" };
            let window_title = format!(
                "GENtle - {}{} (v{})",
                self.current_project_name(),
                dirty_marker,
                about::GENTLE_DISPLAY_VERSION
            );
            if self.last_applied_window_title != window_title {
                ctx.send_viewport_cmd(egui::ViewportCommand::Title(window_title.clone()));
                self.last_applied_window_title = window_title;
            }

            let open_project = KeyboardShortcut::new(Modifiers::COMMAND, Key::O);
            let new_project = KeyboardShortcut::new(Modifiers::COMMAND, Key::N);
            let open_sequence =
                KeyboardShortcut::new(Modifiers::COMMAND | Modifiers::SHIFT, Key::O);
            let open_retrieve_genome =
                KeyboardShortcut::new(Modifiers::COMMAND | Modifiers::SHIFT, Key::G);
            let open_prepare_genome =
                KeyboardShortcut::new(Modifiers::COMMAND | Modifiers::SHIFT, Key::P);
            let open_blast_genome =
                KeyboardShortcut::new(Modifiers::COMMAND | Modifiers::SHIFT, Key::L);
            let open_import_bed_track =
                KeyboardShortcut::new(Modifiers::COMMAND | Modifiers::SHIFT, Key::B);
            let open_agent_assistant =
                KeyboardShortcut::new(Modifiers::COMMAND | Modifiers::SHIFT, Key::A);
            let save_project = KeyboardShortcut::new(Modifiers::COMMAND, Key::S);
            let close_project =
                KeyboardShortcut::new(Modifiers::COMMAND | Modifiers::SHIFT, Key::W);
            let open_configuration = KeyboardShortcut::new(Modifiers::COMMAND, Key::Comma);
            let focus_main_window = KeyboardShortcut::new(Modifiers::COMMAND, Key::Backtick);
            let open_command_palette = KeyboardShortcut::new(Modifiers::COMMAND, Key::K);
            let undo_shortcut = KeyboardShortcut::new(Modifiers::COMMAND, Key::Z);
            let redo_shortcut_shift =
                KeyboardShortcut::new(Modifiers::COMMAND | Modifiers::SHIFT, Key::Z);
            let redo_shortcut_y = KeyboardShortcut::new(Modifiers::COMMAND, Key::Y);
            if ctx.input_mut(|i| i.consume_shortcut(&new_project)) {
                self.request_project_action(ProjectAction::New);
            }
            if ctx.input_mut(|i| i.consume_shortcut(&open_project)) {
                self.request_project_action(ProjectAction::Open);
            }
            if ctx.input_mut(|i| i.consume_shortcut(&open_sequence)) {
                self.prompt_open_sequence();
            }
            if ctx.input_mut(|i| i.consume_shortcut(&open_retrieve_genome)) {
                self.open_reference_genome_retrieve_dialog();
            }
            if ctx.input_mut(|i| i.consume_shortcut(&open_prepare_genome)) {
                self.open_reference_genome_prepare_dialog();
            }
            if ctx.input_mut(|i| i.consume_shortcut(&open_blast_genome)) {
                self.open_reference_genome_blast_dialog();
            }
            if ctx.input_mut(|i| i.consume_shortcut(&open_import_bed_track)) {
                self.open_genome_bed_track_dialog();
            }
            if ctx.input_mut(|i| i.consume_shortcut(&open_agent_assistant)) {
                self.open_agent_assistant_dialog();
            }
            if ctx.input_mut(|i| i.consume_shortcut(&save_project)) {
                let _ = self.save_current_project();
            }
            if ctx.input_mut(|i| i.consume_shortcut(&close_project)) {
                self.request_project_action(ProjectAction::Close);
            }
            if ctx.input_mut(|i| i.consume_shortcut(&open_configuration)) {
                self.open_configuration_dialog();
            }
            if ctx.input_mut(|i| i.consume_shortcut(&focus_main_window)) {
                self.queue_focus_viewport(ViewportId::ROOT);
            }
            if ctx.input_mut(|i| i.consume_shortcut(&open_command_palette)) {
                self.open_command_palette_dialog();
            }
            if ctx.input_mut(|i| i.consume_shortcut(&undo_shortcut)) {
                self.undo_last_operation();
            }
            if ctx.input_mut(|i| i.consume_shortcut(&redo_shortcut_shift))
                || ctx.input_mut(|i| i.consume_shortcut(&redo_shortcut_y))
            {
                self.redo_last_operation();
            }

            self.poll_prepare_reference_genome_task(ctx);
            self.poll_reference_genome_blast_task(ctx);
            self.poll_genome_track_import_task(ctx);
            self.poll_agent_assistant_task(ctx);
            self.poll_agent_model_discovery_task(ctx);
            self.sync_tracked_bed_tracks_for_new_anchors();

            // Show menu bar
            egui::TopBottomPanel::top("top").show(ctx, |ui| {
                self.render_menu_bar(ui);
            });

            // Show main window
            egui::CentralPanel::default().show(ctx, |ui| {
                window_backdrop::paint_window_backdrop(
                    ui,
                    WindowBackdropKind::Main,
                    &self.window_backdrops,
                );
                self.render_main_lineage(ui);
                if project_dirty {
                    ui.label("Status: unsaved changes");
                } else {
                    ui.label("Status: saved");
                }
            });
            self.render_reference_genome_prepare_dialog(ctx);
            self.render_reference_genome_retrieve_dialog(ctx);
            self.render_reference_genome_blast_dialog(ctx);
            self.render_reference_genome_inspector_dialog(ctx);
            self.render_genome_bed_track_dialog(ctx);
            self.render_agent_assistant_dialog(ctx);
            self.render_configuration_dialog(ctx);
            self.render_help_dialog(ctx);
            self.render_about_dialog(ctx);
            self.render_command_palette_dialog(ctx);
            self.render_jobs_panel(ctx);
            self.render_history_panel(ctx);
            self.render_unsaved_changes_dialog(ctx);
            self.render_status_bar(ctx);

            // Open new windows
            let mut new_windows: Vec<Window> = self.new_windows.drain(..).collect();
            for window in new_windows.drain(..) {
                let id = format!("Viewport {}", self.viewport_id_counter);
                let id = ViewportId::from_hash_of(id);
                self.windows.insert(id, Arc::new(RwLock::new(window)));
                self.viewport_id_counter += 1;
            }

            // Close windows
            if let Ok(mut to_close) = self.windows_to_close.write() {
                for id in to_close.drain(..) {
                    self.windows.remove(&id);
                }
            } else {
                eprintln!("W GENtleApp: close-queue lock poisoned");
            }

            // Show windows
            for (id, window) in self.windows.iter() {
                let id = id.to_owned();
                self.show_window(ctx, id, window.clone());
            }

            if !self.pending_focus_viewports.is_empty() {
                let to_focus: Vec<ViewportId> = self.pending_focus_viewports.drain(..).collect();
                for viewport_id in to_focus {
                    self.focus_window_viewport(ctx, viewport_id);
                }
            }
        }));
        if update_result.is_err() {
            eprintln!("E GENtleApp: recovered from panic in app update");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{
        BackgroundJobEventPhase, BackgroundJobKind, EngineError, ErrorCode, GENtleApp,
        GenomePrepareTask, GenomePrepareTaskMessage, MAX_RECENT_PROJECTS,
    };
    use eframe::egui;
    use std::{
        fs,
        sync::{
            Arc,
            atomic::{AtomicBool, Ordering},
            mpsc,
        },
        time::Instant,
    };
    use tempfile::tempdir;

    #[test]
    fn load_help_doc_returns_fallback_when_missing() {
        let loaded = GENtleApp::load_help_doc(
            "/definitely/missing/gentle/help-doc.md",
            "fallback help markdown",
        );
        assert_eq!(loaded, "fallback help markdown");
    }

    #[test]
    fn load_help_doc_uses_runtime_file_and_rewrites_relative_images() {
        let temp = tempdir().unwrap();
        let docs_dir = temp.path().join("docs");
        let images_dir = docs_dir.join("images");
        fs::create_dir_all(&images_dir).unwrap();

        let image_path = images_dir.join("gui.png");
        fs::write(&image_path, b"fake image").unwrap();

        let markdown_path = docs_dir.join("gui.md");
        fs::write(
            &markdown_path,
            "# Help\n\n![GUI](<images/gui.png> \"GUI Screenshot\")\n",
        )
        .unwrap();

        let loaded = GENtleApp::load_help_doc(markdown_path.to_str().unwrap(), "fallback");
        let abs_image = image_path.canonicalize().unwrap();
        let expected = format!(
            "![GUI](<{}> \"GUI Screenshot\")",
            abs_image.to_string_lossy()
        );
        assert!(loaded.contains(&expected), "{loaded}");
    }

    #[test]
    fn rewrite_markdown_relative_image_links_handles_paths_with_parentheses() {
        let temp = tempdir().unwrap();
        let docs_dir = temp.path().join("docs");
        let images_dir = docs_dir.join("images");
        fs::create_dir_all(&images_dir).unwrap();

        let image_path = images_dir.join("gui(1).png");
        fs::write(&image_path, b"fake image").unwrap();

        let markdown = "![Shot](images/gui(1).png)\n";
        let rewritten = GENtleApp::rewrite_markdown_relative_image_links(markdown, &docs_dir);
        let abs_image = image_path.canonicalize().unwrap();
        let expected = format!("![Shot]({})", abs_image.to_string_lossy());
        assert!(rewritten.contains(&expected), "{rewritten}");
    }

    #[test]
    fn rewrite_markdown_relative_image_links_keeps_absolute_and_reference_images() {
        let markdown = "![Web](https://example.com/gui.png)\n![Root](/tmp/gui.png)\n![ByRef][img]\n\n[img]: images/gui.png\n";
        let temp = tempdir().unwrap();
        let rewritten = GENtleApp::rewrite_markdown_relative_image_links(markdown, temp.path());
        assert_eq!(rewritten, markdown);
    }

    #[test]
    fn collect_help_markdown_images_extracts_inline_targets_and_labels() {
        let markdown = "# Help\n\n![Main](<docs/screenshots/main.png> \"Main window\")\n\n![Alt only](images/alt.png)\n";
        let images = GENtleApp::collect_help_markdown_images(markdown);
        assert_eq!(images.len(), 2);
        assert_eq!(images[0].path, "docs/screenshots/main.png");
        assert_eq!(images[0].title, "Main window");
        assert_eq!(images[0].alt, "Main");
        assert_eq!(images[1].path, "images/alt.png");
        assert_eq!(images[1].title, "");
        assert_eq!(images[1].alt, "Alt only");
    }

    #[test]
    fn help_image_caption_prefers_title_then_alt_then_filename() {
        let title_first = super::HelpMarkdownImage {
            alt: "Alt".to_string(),
            title: "Caption".to_string(),
            path: "/tmp/ignored.png".to_string(),
        };
        assert_eq!(GENtleApp::help_image_caption(&title_first), "Caption");

        let alt_second = super::HelpMarkdownImage {
            alt: "Fallback alt".to_string(),
            title: String::new(),
            path: "/tmp/ignored.png".to_string(),
        };
        assert_eq!(GENtleApp::help_image_caption(&alt_second), "Fallback alt");

        let filename_last = super::HelpMarkdownImage {
            alt: String::new(),
            title: String::new(),
            path: "/tmp/final-name.png".to_string(),
        };
        assert_eq!(
            GENtleApp::help_image_caption(&filename_last),
            "final-name.png"
        );
    }

    #[test]
    fn normalize_recent_project_paths_deduplicates_and_limits() {
        let mut paths = Vec::new();
        for idx in 0..(MAX_RECENT_PROJECTS + 4) {
            paths.push(format!("/tmp/gentle_project_{idx}.gentle.json"));
        }
        paths.insert(2, "/tmp/gentle_project_1.gentle.json".to_string());
        paths.insert(4, "   ".to_string());

        let normalized = GENtleApp::normalize_recent_project_paths(paths);
        assert_eq!(normalized.len(), MAX_RECENT_PROJECTS);
        assert_eq!(
            normalized[0],
            GENtleApp::normalize_project_path("/tmp/gentle_project_0.gentle.json")
        );
        assert_eq!(
            normalized[1],
            GENtleApp::normalize_project_path("/tmp/gentle_project_1.gentle.json")
        );
        let needle = GENtleApp::normalize_project_path("/tmp/gentle_project_1.gentle.json");
        assert_eq!(normalized.iter().filter(|p| *p == &needle).count(), 1);
    }

    #[test]
    fn recent_project_menu_label_shows_name_and_parent() {
        let temp = tempdir().unwrap();
        let project_path = temp.path().join("my_project.gentle.json");
        let label = GENtleApp::recent_project_menu_label(project_path.to_string_lossy().as_ref());
        assert!(label.starts_with("my_project.gentle.json ("));
    }

    #[test]
    fn request_prepare_cancel_is_idempotent() {
        let mut app = GENtleApp::default();
        let (_tx, rx) = mpsc::channel::<GenomePrepareTaskMessage>();
        app.genome_prepare_task = Some(GenomePrepareTask {
            job_id: 42,
            started: Instant::now(),
            cancel_requested: Arc::new(AtomicBool::new(false)),
            timeout_seconds: None,
            receiver: rx,
        });

        app.request_prepare_task_cancel("test");
        app.request_prepare_task_cancel("test");

        let cancel_events = app
            .job_event_log
            .iter()
            .filter(|event| {
                event.kind == BackgroundJobKind::PrepareGenome
                    && event.phase == BackgroundJobEventPhase::CancelRequested
                    && event.job_id == Some(42)
            })
            .count();
        assert_eq!(cancel_events, 1);
        assert!(
            app.genome_prepare_task
                .as_ref()
                .unwrap()
                .cancel_requested
                .load(Ordering::Relaxed)
        );
    }

    #[test]
    fn poll_prepare_ignores_stale_job_messages() {
        let mut app = GENtleApp::default();
        let (tx, rx) = mpsc::channel::<GenomePrepareTaskMessage>();
        app.genome_prepare_task = Some(GenomePrepareTask {
            job_id: 7,
            started: Instant::now(),
            cancel_requested: Arc::new(AtomicBool::new(false)),
            timeout_seconds: None,
            receiver: rx,
        });

        tx.send(GenomePrepareTaskMessage::Done {
            job_id: 6,
            result: Err(EngineError {
                code: ErrorCode::Internal,
                message: "stale".to_string(),
            }),
        })
        .unwrap();
        tx.send(GenomePrepareTaskMessage::Done {
            job_id: 7,
            result: Err(EngineError {
                code: ErrorCode::Internal,
                message: "actual".to_string(),
            }),
        })
        .unwrap();

        app.poll_prepare_reference_genome_task(&egui::Context::default());

        assert!(app.genome_prepare_task.is_none());
        assert!(app.job_event_log.iter().any(|event| {
            event.kind == BackgroundJobKind::PrepareGenome
                && event.phase == BackgroundJobEventPhase::IgnoredStale
                && event.job_id == Some(6)
        }));
        assert!(app.job_event_log.iter().any(|event| {
            event.kind == BackgroundJobKind::PrepareGenome
                && event.phase == BackgroundJobEventPhase::Failed
                && event.job_id == Some(7)
        }));
    }
}
