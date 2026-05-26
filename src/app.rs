//! Top-level GUI application wiring and event loop state.
//!
//! `GENtleApp` owns application-wide GUI state that spans projects, windows,
//! background jobs, help/settings surfaces, lineage views, and specialist
//! dialogs. It is the GUI-side coordinator above individual `MainAreaDna`
//! sequence windows.
//!
//! File map (source order):
//! - extracted helper modules first:
//!   - `app/agent_assistant_config.rs` for Agent Assistant templates/defaults
//!   - `app/routine_and_agent_assistant_ui.rs` for Routine Assistant and
//!     Agent Assistant dialog/rendering helpers
//!   - `app/sequence_ingress_dialogs_ui.rs` for GenBank, dbSNP, UniProt,
//!     Ensembl protein, reverse-translation, and protease dialogs
//!   - `app/gibson_ui.rs` for the Gibson specialist dialog and preview/apply
//!     helpers
//!   - `app/genome_catalog_ui.rs` for reference-genome catalog dialogs
//!     and cache/track rendering helpers
//!   - `app/main_lineage_ui.rs` for project overview and main lineage
//!     workspace rendering helpers
//!   - `app/help_docs.rs` for markdown/help loading and rewrite helpers
//!   - `app/rack_workspace_ui.rs` for rack, arrangement-gel, and rack-label
//!     workspace helpers
//!   - `app/window_registry.rs` for open-window listing/focus helpers and
//!     deferred child-window placement
//! - top-level imports/constants plus native-menu open/focus request helpers
//! - persisted configuration and background-job state records
//! - genome-preparation, BLAST, lineage, Gibson, planning, and command-palette
//!   UI structs/enums
//! - `GENtleApp` is the main application state object
//! - the main `impl GENtleApp` blocks then handle startup/persistence, menu and
//!   dialog orchestration, background-task polling, lineage/planning/routine
//!   windows, and cross-window coordination
//! - the `eframe::App` impl near the end is the event-loop entry point
//! - tests and large embedded markdown fixtures live at the tail of the file
//!
//! Start here when...
//! - changing app-wide menus, dialogs, or background jobs: search for
//!   `pub struct GENtleApp`
//! - changing help/tutorial markdown loading: look in `app/help_docs.rs` first
//! - debugging per-sequence editing/rendering: jump to `src/main_area_dna.rs`
//!   instead of adding more logic here

#[path = "app/help_docs.rs"]
mod help_docs;

#[path = "app/window_registry.rs"]
mod window_registry;

pub use window_registry::{GuiProminentGlossaryEntry, gui_prominent_glossary_entries};

#[path = "app/jaspar_expert.rs"]
mod jaspar_expert;

#[path = "app/history_ui.rs"]
mod history_ui;

#[path = "app/root_ui.rs"]
mod root_ui;

#[path = "app/splash_ui.rs"]
mod splash_ui;

#[path = "app/configuration_ui.rs"]
mod configuration_ui;

#[path = "app/external_services_ui.rs"]
mod external_services_ui;

#[path = "app/clawbio_bridge.rs"]
mod clawbio_bridge;

#[path = "app/clawbio_ui.rs"]
mod clawbio_ui;

#[path = "app/agent_assistant_config.rs"]
mod agent_assistant_config;

#[path = "app/routine_and_agent_assistant_ui.rs"]
mod routine_and_agent_assistant_ui;

#[path = "app/genome_catalog_ui.rs"]
mod genome_catalog_ui;

#[path = "app/sequence_ingress_dialogs_ui.rs"]
mod sequence_ingress_dialogs_ui;

#[path = "app/main_lineage_ui.rs"]
mod main_lineage_ui;

#[path = "app/rack_workspace_ui.rs"]
mod rack_workspace_ui;

#[path = "app/gibson_ui.rs"]
mod gibson_ui;

use std::{
    collections::hash_map::DefaultHasher,
    collections::{BTreeMap, BTreeSet, HashMap, HashSet},
    env, fs,
    hash::{Hash, Hasher},
    io::ErrorKind,
    panic::{AssertUnwindSafe, catch_unwind},
    path::{Path, PathBuf},
    process::Command,
    sync::{
        Arc, Mutex, RwLock,
        atomic::{AtomicBool, AtomicU64, Ordering},
        mpsc,
    },
    time::{Duration, Instant, SystemTime, UNIX_EPOCH},
};

use crate::{
    about,
    agent_bridge::{
        AGENT_BASE_URL_ENV, AGENT_CONNECT_TIMEOUT_SECS_ENV, AGENT_MAX_RESPONSE_BYTES_ENV,
        AGENT_MAX_RETRIES_ENV, AGENT_MODEL_ENV, AGENT_READ_TIMEOUT_SECS_ENV,
        AGENT_TIMEOUT_SECS_ENV, ANTHROPIC_API_KEY_AUTH_HINT, ANTHROPIC_API_KEY_ENV,
        AgentExecutionIntent, AgentInvocationOutcome, AgentResponse, AgentSystemSpec,
        AgentSystemTransport, DEFAULT_AGENT_SYSTEM_CATALOG_PATH, MISTRAL_API_KEY_AUTH_HINT,
        MISTRAL_API_KEY_ENV, OPENAI_API_KEY_ENV, OPENAI_BILLING_URL,
        OPENAI_COMPAT_UNSPECIFIED_MODEL, OPENAI_USAGE_URL, agent_system_availability,
        invoke_agent_support_with_env_overrides, load_agent_system_catalog,
    },
    agent_transport::{
        AgentLiveProbeStatusClass, AgentSystemPreflight, build_agent_system_preflight_with_live,
        discover_models_for_agent_system,
    },
    dna_sequence::{self, DNAsequence},
    engine::{
        BIGWIG_TO_BEDGRAPH_ENV_BIN, BlastHitFeatureInput, BlastInvocationProvenance,
        ConstructReasoningGraph, DEFAULT_BIGWIG_TO_BEDGRAPH_BIN, DEFAULT_HOST_PROFILE_CATALOG_PATH,
        DEFAULT_JASPAR_PRESENTATION_RANDOM_SEED,
        DEFAULT_JASPAR_PRESENTATION_RANDOM_SEQUENCE_LENGTH_BP, DbSnpFetchProgress, DbSnpFetchStage,
        DisplaySettings, DisplayTarget, Engine, EngineError, ErrorCode, FeatureExpertTarget,
        GenomeAnnotationScope, GenomeGeneExtractMode, GenomeTrackImportProgress, GenomeTrackSource,
        GenomeTrackSubscription, GenomeTrackSyncReport, GentleEngine, HostProfileRecord,
        JasparCatalogReport, JasparEntryExpertView, LabAssistantInstructionsFormat,
        LineageMacroPortBinding, LinearSequenceLetterLayoutMode, MacroTemplateSuggestion, OpResult,
        Operation, OperationProgress, PlanningEstimate, PlanningObjective, PlanningProfile,
        PlanningProfileScope, PlanningSuggestionStatus, ProjectState, ProteaseDigestReport,
        ProteinToDnaHandoffRankingGoal, ROUTINE_DECISION_TRACE_SCHEMA,
        ROUTINE_DECISION_TRACE_STORE_SCHEMA, ROUTINE_DECISION_TRACES_METADATA_KEY, Rack,
        RackAuthoringTemplate, RackCarrierLabelPreset, RackFillDirection, RackLabelSheetPreset,
        RackOccupant, RackPhysicalTemplateKind, RackProfileKind, RenderSvgMode,
        RestrictionEnzymeDisplayMode, ReverseTranslationReport, RoutineDecisionTrace,
        RoutineDecisionTraceCandidateScore, RoutineDecisionTraceComparison,
        RoutineDecisionTraceDisambiguationAnswer, RoutineDecisionTraceDisambiguationQuestion,
        RoutineDecisionTraceExportEvent, RoutineDecisionTracePreflightSnapshot,
        RoutineDecisionTraceStore, RoutinePreferenceContextRecord, SequenceGenomeAnchorSummary,
        SequenceScanTarget, TranslationSpeedMark, TranslationSpeedProfile,
        UniprotFeatureCodingDnaQueryMode, UniprotFeatureCodingDnaQueryReport,
        UniprotProjectionAuditParityReport, UniprotProjectionAuditReport,
    },
    engine_shell::{
        ShellCommand, ShellExecutionOptions, UiIntentAction, UiIntentTarget,
        execute_shell_command_with_options, normalize_pasted_iupac_sequence, parse_shell_line,
    },
    ensembl_protein::{EnsemblProteinEntry, EnsemblProteinEntrySummary},
    enzymes,
    genomes::{
        BLASTN_ENV_BIN, DEFAULT_BLASTN_BIN, DEFAULT_GENOME_CACHE_DIR, DEFAULT_GENOME_CATALOG_PATH,
        DEFAULT_HELPER_GENOME_CATALOG_PATH, DEFAULT_MAKEBLASTDB_BIN, EnsemblCatalogUpdatePreview,
        EnsemblInstallableGenomeCatalog, EnsemblQuickInstallPreview, GenomeBlastReport,
        GenomeCatalog, GenomeCatalogListEntry, GenomeChromosomeRecord, GenomeGeneRecord,
        GenomeSourcePlan, HelperConstructInterpretation, HelperVectorCard,
        HelperVectorCatalogDoctorIssue, MAKEBLASTDB_ENV_BIN, PREPARE_GENOME_TIMEOUT_SECS_ENV,
        PrepareGenomePlan, PrepareGenomePlanStep, PrepareGenomeProgress, PrepareGenomeStepId,
        PreparedCacheArtifactGroup, PreparedCacheCleanupMode, PreparedCacheCleanupRequest,
        PreparedCacheInspectionEntry, PreparedCacheInspectionReport, PreparedGenomeInspection,
        configured_helper_genome_cache_dir, configured_reference_genome_cache_dir,
    },
    gibson_planning::{
        GibsonAssemblyPlan, GibsonAssemblyPreview, GibsonDesignAdjustmentTarget,
        GibsonDestinationOpeningSuggestion, GibsonPlanAssemblyMember, GibsonPlanDesignTargets,
        GibsonPlanDestination, GibsonPlanEndStrategy, GibsonPlanFragment, GibsonPlanJunction,
        GibsonPlanOpening, GibsonPlanOverlapPartition, GibsonPlanProduct,
        GibsonPlanUniquenessChecks, GibsonPlanValidationPolicy, GibsonSuggestedDesignAdjustment,
        suggest_gibson_destination_openings,
    },
    i18n::{I18n, UiLanguage},
    icons::APP_ICON,
    lineage_export::{
        LineageSvgEdge, LineageSvgNode, LineageSvgNodeKind, export_projected_lineage_svg,
    },
    mcp_server::DEFAULT_MCP_STATE_PATH,
    protocol_cartoon::{
        ProtocolCartoonKind, protocol_cartoon_template_for_kind, render_protocol_cartoon_spec_svg,
        resolve_protocol_cartoon_template_with_bindings,
    },
    resource_sync,
    scroll_input_policy::{self, WheelIntent, ZoomDirection},
    shell_docs::{
        shell_help_markdown as render_shell_help_markdown,
        shell_help_markdown_from_glossary_json as render_shell_help_markdown_from_glossary_json,
    },
    tf_motifs, tool_overrides,
    uniprot::{UniprotEntrySummary, UniprotGenomeProjectionSummary},
    window::Window,
    window_backdrop::{
        self, WindowBackdropKind, WindowBackdropSettings, with_window_content_inset,
    },
    workflow_examples::{
        DEFAULT_TUTORIAL_CATALOG_PATH, DEFAULT_TUTORIAL_MANIFEST_PATH, DEFAULT_TUTORIAL_OUTPUT_DIR,
        DEFAULT_TUTORIAL_REVIEW_MANIFEST_PATH, DEFAULT_TUTORIAL_SOURCE_DIR,
        DEFAULT_WORKFLOW_EXAMPLE_DIR, ExampleTestMode, TutorialTier, WorkflowExample,
        build_tutorial_feedback_context_text, load_tutorial_catalog, load_tutorial_manifest,
        load_tutorial_review_manifest, load_workflow_examples,
        run_example_workflow_for_project_state_with_progress, validate_example_required_files,
        validate_tutorial_manifest_against_examples,
    },
};
use agent_assistant_config::{
    AGENT_PROMPT_TEMPLATE_DEFAULT_ID, agent_prompt_template_includes_state_summary_by_default,
    agent_prompt_template_label, agent_prompt_template_options, agent_prompt_template_text,
    default_agent_connect_timeout_secs_string, default_agent_max_response_bytes_string,
    default_agent_max_retries_string, default_agent_read_timeout_secs_string,
    default_agent_timeout_secs_string, normalize_agent_model_name,
    preferred_anthropic_agent_system_id, preferred_local_agent_system_id,
    preferred_mistral_agent_system_id, preferred_openai_agent_system_id,
};
use anyhow::{Result, anyhow};
use eframe::egui::{self, Key, KeyboardShortcut, Modifiers, Pos2, Ui, Vec2, ViewportId};
use egui_commonmark::{CommonMarkCache, CommonMarkViewer};
use gentle_gui::theme::{self, SciencePalette};
use gibson_ui::*;
use pulldown_cmark::{Event, LinkType, Parser, Tag};
use rack_workspace_ui::*;
use regex::{Regex, RegexBuilder};
use resvg::{self, tiny_skia, usvg};
use serde::{Deserialize, Serialize};
use sha1::{Digest, Sha1};

const GUI_MANUAL_MD: &str = include_str!("../docs/gui.md");
const CLI_MANUAL_MD: &str = include_str!("../docs/cli.md");
const AGENT_INTERFACE_MD: &str = include_str!("../docs/agent_interface.md");
const REVIEWER_PREVIEW_MD: &str = include_str!("../docs/reviewer_preview.md");
const APP_CONFIGURATION_FILE_NAME: &str = ".gentle_gui_settings.json";
const APP_CONFIGURATION_SCHEMA_VERSION: u32 = 3;
const MAX_RECENT_PROJECTS: usize = 12;
const MAX_BACKGROUND_JOB_EVENTS: usize = 200;
const MAX_BACKGROUND_JOB_RETRY_SNAPSHOTS: usize = 120;
const MAX_RETRY_SNAPSHOT_CLEANUP_AUDIT_ENTRIES: usize = 240;
const BACKGROUND_JOBS_RECENT_JOB_EVENTS_SCROLL_ID: &str =
    "background_jobs_recent_job_events_scroll";
const BACKGROUND_JOBS_RETRY_SNAPSHOTS_REMOVED_PREVIEW_SCROLL_ID: &str =
    "background_jobs_retry_snapshots_removed_preview_scroll";
const BACKGROUND_JOBS_RETRY_SNAPSHOTS_RETAINED_PREVIEW_SCROLL_ID: &str =
    "background_jobs_retry_snapshots_retained_preview_scroll";
const BACKGROUND_JOBS_RETRY_SNAPSHOTS_SCROLL_ID: &str = "background_jobs_retry_snapshots_scroll";
const BACKGROUND_JOBS_RETRY_CLEANUP_AUDIT_SCROLL_ID: &str =
    "background_jobs_retry_cleanup_audit_scroll";
const OPERATION_HISTORY_SCROLL_ID: &str = "operation_history_scroll";
const LINEAGE_GRAPH_WORKSPACE_METADATA_KEY: &str = "gui.lineage_graph.workspace";
const RACK_WORKSPACE_METADATA_KEY: &str = "gui.rack.workspace";
const LINEAGE_NODE_OFFSETS_METADATA_KEY: &str = "gui.lineage_graph.node_offsets";
const LINEAGE_NODE_GROUPS_METADATA_KEY: &str = "gui.lineage_graph.node_groups";
const LINEAGE_MAIN_TOP_PANEL_MIN_HEIGHT: f32 = 180.0;
const BACKGROUND_JOB_HISTORY_METADATA_KEY: &str = "gui.background_job_history";
const BACKGROUND_JOB_HISTORY_SCHEMA: &str = "gentle.gui_background_job_history.v1";
const DEFAULT_LINEAGE_MAIN_SPLIT_FRACTION: f32 = 420.0 / (420.0 + 220.0);
const DEFAULT_LINEAGE_CONTAINER_ARRANGEMENT_SPLIT_FRACTION: f32 = 0.58;
const LAB_ASSISTANT_REPORT_EXTENSIONS: &[&str] = &["odt", "docx", "md"];
const DEFAULT_CLONING_PATTERN_CATALOG_DIR: &str = "assets/cloning_patterns_catalog";
const DEFAULT_CLONING_PATTERN_PACK_PATH: &str = "assets/cloning_patterns.json";
const DEFAULT_CLONING_ROUTINE_CATALOG_PATH: &str = "assets/cloning_routines.json";
const DEFAULT_HELPER_GENOME_CACHE_DIR: &str = "data/helper_genomes";
const DEFAULT_DBSNP_TUTORIAL_RS_ID: &str = "rs9923231";
const RACK_HELP_AUTO_MINIMIZE_MOVE_THRESHOLD: u32 = 3;
const GUI_OPENAI_DEFAULT_BASE_URL: &str = "https://api.openai.com/v1";
const GUI_ANTHROPIC_DEFAULT_MODEL: &str = "claude-sonnet-4-6";
const GUI_ANTHROPIC_DEFAULT_BASE_URL: &str = "https://api.anthropic.com/v1";
const GUI_MISTRAL_DEFAULT_MODEL: &str = "mistral-large-latest";
const GUI_MISTRAL_DEFAULT_BASE_URL: &str = "https://api.mistral.ai/v1";
const GUI_OPENAI_COMPAT_DEFAULT_BASE_URL: &str = "http://127.0.0.1:11434/v1";
const WINDOW_OPEN_SLOW_THRESHOLD_MS: u128 = 400;
// This margin is subtracted on both sides from the hosted window's inner
// max-size, so keep it small; large values cap DNA viewers far inside the
// root window, while zero lets inner content fight the outer chrome.
const EMBEDDED_SEQUENCE_WINDOW_DRAG_MARGIN_X_PX: f32 = 32.0;
const EMBEDDED_SEQUENCE_WINDOW_DRAG_MARGIN_Y_PX: f32 = 32.0;
const HELP_MARKDOWN_REFLOW_DELTA_PX: f32 = 8.0;
const MACOS_NATIVE_CHILD_VIEWPORTS_ENV: &str = "GENTLE_MACOS_NATIVE_CHILD_VIEWPORTS";
static NATIVE_HELP_OPEN_REQUESTED: AtomicBool = AtomicBool::new(false);
static NATIVE_HELP_OPEN_REQUESTED_AT_MS: AtomicU64 = AtomicU64::new(0);
static NATIVE_SETTINGS_OPEN_REQUESTED: AtomicBool = AtomicBool::new(false);
static NATIVE_SETTINGS_OPEN_GRAPHICS_TAB_REQUESTED: AtomicBool = AtomicBool::new(false);
static NATIVE_PCR_DESIGN_OPEN_REQUESTED: AtomicBool = AtomicBool::new(false);
static NATIVE_PCR_DESIGN_SEQ_ID_REQUESTED: Mutex<Option<String>> = Mutex::new(None);
static NATIVE_JASPAR_EXPERT_OPEN_REQUESTED: AtomicBool = AtomicBool::new(false);
static NATIVE_JASPAR_EXPERT_MOTIF_ID_REQUESTED: Mutex<Option<String>> = Mutex::new(None);
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

fn now_unix_ms_u64() -> u64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|v| v.as_millis() as u64)
        .unwrap_or(0)
}

fn default_prepare_timeout_secs_string() -> String {
    env::var(PREPARE_GENOME_TIMEOUT_SECS_ENV)
        .ok()
        .map(|v| v.trim().to_string())
        .filter(|v| !v.is_empty())
        .unwrap_or_default()
}

pub fn request_open_help_from_native_menu() {
    NATIVE_HELP_OPEN_REQUESTED_AT_MS.store(now_unix_ms_u64(), Ordering::SeqCst);
    NATIVE_HELP_OPEN_REQUESTED.store(true, Ordering::SeqCst);
}

pub fn request_open_settings_from_native_menu() {
    NATIVE_SETTINGS_OPEN_GRAPHICS_TAB_REQUESTED.store(false, Ordering::SeqCst);
    NATIVE_SETTINGS_OPEN_REQUESTED.store(true, Ordering::SeqCst);
}

pub fn request_open_graphics_settings_from_native_menu() {
    NATIVE_SETTINGS_OPEN_GRAPHICS_TAB_REQUESTED.store(true, Ordering::SeqCst);
    NATIVE_SETTINGS_OPEN_REQUESTED.store(true, Ordering::SeqCst);
}

pub fn request_open_pcr_design_from_native_menu() {
    if let Ok(mut requested_seq_id) = NATIVE_PCR_DESIGN_SEQ_ID_REQUESTED.lock() {
        requested_seq_id.take();
    }
    NATIVE_PCR_DESIGN_OPEN_REQUESTED.store(true, Ordering::SeqCst);
}

pub fn request_open_pcr_design_for_sequence_from_native_menu(seq_id: &str) {
    if let Ok(mut requested_seq_id) = NATIVE_PCR_DESIGN_SEQ_ID_REQUESTED.lock() {
        let trimmed = seq_id.trim();
        *requested_seq_id = if trimmed.is_empty() {
            None
        } else {
            Some(trimmed.to_string())
        };
    }
    NATIVE_PCR_DESIGN_OPEN_REQUESTED.store(true, Ordering::SeqCst);
}

pub fn request_open_jaspar_expert_for_motif_from_native_menu(motif_id: &str) {
    if let Ok(mut requested_motif_id) = NATIVE_JASPAR_EXPERT_MOTIF_ID_REQUESTED.lock() {
        let trimmed = motif_id.trim();
        *requested_motif_id = if trimmed.is_empty() {
            None
        } else {
            Some(trimmed.to_string())
        };
    }
    NATIVE_JASPAR_EXPERT_OPEN_REQUESTED.store(true, Ordering::SeqCst);
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
    #[serde(default)]
    schema_version: u32,
    rnapkin_executable: String,
    makeblastdb_executable: String,
    blastn_executable: String,
    bigwig_to_bedgraph_executable: String,
    graphics_defaults: DisplaySettings,
    window_backdrops: WindowBackdropSettings,
    ui_language: UiLanguage,
    recent_projects: Vec<String>,
}

impl Default for PersistedConfiguration {
    fn default() -> Self {
        Self {
            schema_version: APP_CONFIGURATION_SCHEMA_VERSION,
            rnapkin_executable: String::new(),
            makeblastdb_executable: String::new(),
            blastn_executable: String::new(),
            bigwig_to_bedgraph_executable: String::new(),
            graphics_defaults: DisplaySettings::default(),
            window_backdrops: WindowBackdropSettings::default(),
            ui_language: UiLanguage::default(),
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

    fn as_str(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Bed => "bed",
            Self::BigWig => "bigwig",
            Self::Vcf => "vcf",
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
enum GenomeDialogScope {
    #[default]
    Reference,
    Helper,
}

impl GenomeDialogScope {
    fn prepare_title(self) -> &'static str {
        match self {
            Self::Reference => "Prepare Reference Genome",
            Self::Helper => "Prepare Helper Genome",
        }
    }

    fn retrieve_title(self) -> &'static str {
        match self {
            Self::Reference => "Retrieve Genomic Sequence",
            Self::Helper => "Retrieve Helper Sequence",
        }
    }

    fn blast_title(self) -> &'static str {
        match self {
            Self::Reference => "BLAST Genome Sequence",
            Self::Helper => "BLAST Helper Sequence",
        }
    }

    fn description(self) -> &'static str {
        match self {
            Self::Reference => "reference genome",
            Self::Helper => "helper genome",
        }
    }

    fn default_catalog_path(self) -> &'static str {
        match self {
            Self::Reference => DEFAULT_GENOME_CATALOG_PATH,
            Self::Helper => DEFAULT_HELPER_GENOME_CATALOG_PATH,
        }
    }

    fn default_cache_dir(self) -> &'static str {
        match self {
            Self::Reference => DEFAULT_GENOME_CACHE_DIR,
            Self::Helper => DEFAULT_HELPER_GENOME_CACHE_DIR,
        }
    }

    fn empty_prepared_hint(self) -> String {
        format!(
            "No prepared genomes are available in this catalog/cache. Use '{}...' first.",
            self.prepare_title()
        )
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PrepareGenomeDialogPrimaryAction {
    None,
    Prepare,
    Reindex,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum GenomePrepareLaunchMode {
    Prepare,
    ReindexCachedFiles,
    RefreshFromSources,
}

impl GenomePrepareLaunchMode {
    fn progress_label(self) -> &'static str {
        match self {
            Self::Prepare => "Preparing",
            Self::ReindexCachedFiles => "Reindexing",
            Self::RefreshFromSources => "Refreshing",
        }
    }

    fn job_summary(self) -> &'static str {
        match self {
            Self::Prepare => "Prepare genome started",
            Self::ReindexCachedFiles => "Reindex genome started",
            Self::RefreshFromSources => "Refresh genome started",
        }
    }

    fn result_prefix(self, cancellation_requested: bool) -> &'static str {
        match (self, cancellation_requested) {
            (Self::Prepare, false) => "Prepare genome: ok",
            (Self::Prepare, true) => "Prepare genome finished after cancellation request",
            (Self::ReindexCachedFiles, false) => "Reindex genome: ok",
            (Self::ReindexCachedFiles, true) => {
                "Reindex genome finished after cancellation request"
            }
            (Self::RefreshFromSources, false) => "Refresh genome: ok",
            (Self::RefreshFromSources, true) => {
                "Refresh genome finished after cancellation request"
            }
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PrepareGenomeUiStepStatus {
    Pending,
    Running,
    Completed,
    Failed,
    Cancelled,
}

#[derive(Debug, Clone)]
struct PrepareGenomeUiStepState {
    step_id: PrepareGenomeStepId,
    label: String,
    operation_summary: String,
    determinate_hint: bool,
    status: PrepareGenomeUiStepStatus,
    progress_fraction: Option<f32>,
    detail: String,
    raw_phase: Option<String>,
    bytes_done: Option<u64>,
    bytes_total: Option<u64>,
    eta_remaining: Option<Duration>,
}

impl PrepareGenomeUiStepState {
    fn from_plan_step(step: PrepareGenomePlanStep) -> Self {
        Self {
            step_id: step.step_id,
            label: step.label,
            operation_summary: step.operation_summary,
            determinate_hint: step.determinate_hint,
            status: PrepareGenomeUiStepStatus::Pending,
            progress_fraction: None,
            detail: String::new(),
            raw_phase: None,
            bytes_done: None,
            bytes_total: None,
            eta_remaining: None,
        }
    }
}

#[derive(Debug, Clone)]
struct PrepareGenomeEtaBaseline {
    step_id: PrepareGenomeStepId,
    observed_at: Instant,
    bytes_done: u64,
    bytes_total: u64,
}

#[derive(Debug, Clone)]
struct PrepareGenomeFailureRecovery {
    genome_id: String,
    scope: GenomeDialogScope,
    catalog_path: String,
    cache_dir: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PreparedGenomeReinstallDialogHost {
    Root,
    PrepareDialog,
}

pub struct GENtleApp {
    engine: Arc<RwLock<GentleEngine>>,
    new_windows: Vec<Window>,
    windows: HashMap<ViewportId, Arc<RwLock<Window>>>,
    detached_auxiliary_window_hosts: HashMap<ViewportId, Arc<RwLock<Window>>>,
    windows_to_close: Arc<RwLock<Vec<ViewportId>>>,
    pending_focus_viewports: Vec<ViewportId>,
    pending_window_initial_positions: HashMap<ViewportId, Pos2>,
    viewport_id_counter: usize,
    update_has_run_before: bool,
    splash_started_at: Instant,
    splash_dismissed: bool,
    show_about_dialog: bool,
    show_help_dialog: bool,
    help_doc: HelpDoc,
    help_markdown_cache: CommonMarkCache,
    help_rendered_markdown_cache: HashMap<String, RenderedHelpMarkdownEntry>,
    help_gui_markdown: String,
    help_cli_markdown: String,
    help_agent_interface_markdown: String,
    help_reviewer_preview_markdown: String,
    help_shell_markdown: String,
    help_tutorial_markdown: String,
    help_tutorial_title: String,
    help_tutorial_entries: Vec<HelpTutorialDocEntry>,
    help_tutorial_selected: usize,
    help_shell_interface: ShellHelpInterface,
    help_search_query: String,
    help_search_matches: Vec<HelpSearchMatch>,
    help_search_selected: usize,
    help_focus_search_box: bool,
    help_last_content_width: f32,
    help_selectable_text_mode: bool,
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
    ui_language: UiLanguage,
    configuration_ui_language: UiLanguage,
    configuration_language_dirty: bool,
    i18n: I18n,
    configuration_status: String,
    window_backdrop_path_status_cache: HashMap<String, (egui::Color32, String)>,
    current_project_path: Option<String>,
    recent_project_paths: Vec<String>,
    lineage_graph_view: bool,
    lineage_main_split_fraction: f32,
    lineage_graph_zoom: f32,
    lineage_graph_area_height: f32,
    lineage_container_area_height: f32,
    lineage_container_arrangement_split_fraction: f32,
    lineage_graph_scroll_offset: Vec2,
    lineage_graph_pan_origin: Option<Vec2>,
    lineage_graph_compact_labels: bool,
    lineage_main_split_drag_origin: Option<f32>,
    lineage_main_split_drag_start_y: Option<f32>,
    lineage_container_arrangement_split_drag_origin: Option<f32>,
    lineage_container_arrangement_split_drag_start_y: Option<f32>,
    lineage_graph_selected_node_id: Option<String>,
    lineage_graph_node_offsets: HashMap<String, Vec2>,
    lineage_graph_drag_origin: Option<(String, Vec2)>,
    lineage_graph_offsets_synced_stamp: u64,
    lineage_cache_stamp: u64,
    lineage_cache_valid: bool,
    lineage_rows: Vec<LineageRow>,
    lineage_edges: Vec<(String, String, String)>,
    lineage_op_label_by_id: HashMap<String, String>,
    lineage_reopenable_gibson_op_ids: HashSet<String>,
    lineage_reopenable_pcr_op_seq_ids: HashMap<String, String>,
    lineage_node_groups: Vec<PersistedLineageNodeGroup>,
    lineage_group_form_label: String,
    lineage_group_form_members: String,
    lineage_group_form_representative: String,
    lineage_group_form_editing_id: Option<String>,
    lineage_group_status: String,
    lineage_group_marked_nodes: HashSet<String>,
    lineage_node_rename_target_id: Option<String>,
    lineage_node_rename_text: String,
    lineage_node_remove_target_id: Option<String>,
    lineage_containers: Vec<ContainerRow>,
    lineage_arrangements: Vec<ArrangementRow>,
    lineage_racks: Vec<RackRow>,
    clean_state_fingerprint: u64,
    dirty_cache_stamp: u64,
    last_display_sync_stamp: u64,
    dirty_cache_valid: bool,
    dirty_cache_value: bool,
    dirty_cache_last_deep_check: Instant,
    last_applied_window_title: String,
    last_native_window_entries: Vec<(u64, String)>,
    last_native_active_window_key: Option<u64>,
    native_window_key_to_viewport: HashMap<u64, ViewportId>,
    active_window_menu_key: Option<u64>,
    pending_window_open_timestamps: HashMap<ViewportId, Instant>,
    pending_viewport_focus_timestamps: HashMap<ViewportId, Instant>,
    pending_project_action: Option<ProjectAction>,
    pending_app_quit: bool,
    pending_prepared_genome_reinstall: Option<PreparedGenomeReinstallRequest>,
    pending_ensembl_catalog_update: Option<PendingEnsemblCatalogUpdateDialog>,
    pending_ensembl_installable_genomes: Option<PendingEnsemblInstallableGenomeDialog>,
    pending_ensembl_quick_install: Option<PendingEnsemblQuickInstallDialog>,
    pending_prepared_genome_removal: Option<PendingPreparedGenomeRemovalRequest>,
    pending_genome_catalog_entry_removal: Option<PendingGenomeCatalogEntryRemovalRequest>,
    show_reference_genome_prepare_dialog: bool,
    show_reference_genome_retrieve_dialog: bool,
    show_reference_genome_blast_dialog: bool,
    show_reference_genome_inspector_dialog: bool,
    show_cache_cleanup_dialog: bool,
    show_new_sequence_dialog: bool,
    show_uniprot_dialog: bool,
    show_genbank_dialog: bool,
    genome_dialog_scope: GenomeDialogScope,
    reference_genome_catalog_path: String,
    reference_genome_cache_dir: String,
    helper_genome_catalog_path: String,
    helper_genome_cache_dir: String,
    genome_catalog_path: String,
    genome_cache_dir: String,
    genome_id: String,
    genome_catalog_genomes: Vec<String>,
    genome_catalog_error: String,
    genome_prepare_task: Option<GenomePrepareTask>,
    genome_prepare_progress: Option<PrepareGenomeProgress>,
    genome_prepare_steps: Vec<PrepareGenomeUiStepState>,
    genome_prepare_eta_baseline: Option<PrepareGenomeEtaBaseline>,
    genome_prepare_failure_recovery: Option<PrepareGenomeFailureRecovery>,
    genome_prepare_status: String,
    tutorial_project_task: Option<TutorialProjectTask>,
    tutorial_project_progress_fraction: Option<f32>,
    tutorial_project_progress_label: String,
    tutorial_project_status: String,
    cache_cleanup_scope: CacheCleanupScope,
    cache_cleanup_reference_cache_dir: String,
    cache_cleanup_helper_cache_dir: String,
    cache_cleanup_mode: PreparedCacheCleanupMode,
    cache_cleanup_include_orphans: bool,
    cache_cleanup_inspection: Option<PreparedCacheInspectionReport>,
    cache_cleanup_selected_paths: HashSet<String>,
    cache_cleanup_rebuild_candidates: Vec<PreparedGenomeReinstallRequest>,
    cache_cleanup_status: String,
    cache_cleanup_confirm_pending: bool,
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
    genome_output_id_autofilled: bool,
    genome_gene_extract_mode: GenomeGeneExtractMode,
    genome_gene_promoter_upstream_bp: String,
    genome_annotation_scope: GenomeAnnotationScope,
    genome_max_annotation_features: String,
    genome_include_genomic_annotation: bool,
    genome_retrieve_status: String,
    genome_retrieve_contig_suggestions: Vec<String>,
    new_sequence_text: String,
    new_sequence_output_id: String,
    new_sequence_name: String,
    new_sequence_circular: bool,
    new_sequence_status: String,
    uniprot_query: String,
    uniprot_entry_id: String,
    uniprot_swiss_path: String,
    ensembl_protein_query: String,
    ensembl_protein_entry_id: String,
    ensembl_protein_output_id: String,
    protein_feature_key_include: String,
    protein_feature_key_exclude: String,
    uniprot_linked_accession: String,
    uniprot_linked_as_id: String,
    uniprot_map_seq_id: String,
    uniprot_map_projection_id: String,
    uniprot_map_transcript_id: String,
    uniprot_feature_query: String,
    uniprot_feature_transcript_id: String,
    uniprot_feature_query_mode: UniprotFeatureCodingDnaQueryMode,
    uniprot_feature_speed_profile: Option<TranslationSpeedProfile>,
    uniprot_feature_report: Option<UniprotFeatureCodingDnaQueryReport>,
    uniprot_audit_transcript_id: String,
    uniprot_audit_report_id: String,
    uniprot_audit_parity_report_id: String,
    uniprot_audit_report: Option<UniprotProjectionAuditReport>,
    uniprot_audit_parity_report: Option<UniprotProjectionAuditParityReport>,
    reverse_translate_protein_seq_id: String,
    reverse_translate_output_id: String,
    reverse_translate_speed_profile: Option<TranslationSpeedProfile>,
    reverse_translate_speed_mark: Option<TranslationSpeedMark>,
    reverse_translate_translation_table: String,
    reverse_translate_target_anneal_tm_c: String,
    reverse_translate_anneal_window_bp: String,
    reverse_translation_report: Option<ReverseTranslationReport>,
    protease_digest_names: String,
    protease_digest_output_prefix: String,
    protease_digest_min_length_aa: String,
    protease_digest_materialize: bool,
    protease_digest_report: Option<ProteaseDigestReport>,
    protein_handoff_ranking_goal: ProteinToDnaHandoffRankingGoal,
    protein_handoff_graph: Option<ConstructReasoningGraph>,
    protein_handoff_selected_candidate_id: String,
    protein_handoff_status: String,
    uniprot_status: String,
    genbank_accession: String,
    genbank_as_id: String,
    genbank_status: String,
    dbsnp_rs_id: String,
    dbsnp_genome_id: String,
    dbsnp_flank_bp: String,
    dbsnp_output_id: String,
    dbsnp_status: String,
    dbsnp_fetch_task: Option<DbSnpFetchTask>,
    genome_blast_source_mode: GenomeBlastSourceMode,
    genome_blast_query_manual: String,
    genome_blast_query_seq_id: String,
    genome_blast_query_pool_id: String,
    genome_blast_max_hits: usize,
    genome_blast_task_name: String,
    genome_blast_options_preset: GenomeBlastOptionsPreset,
    genome_blast_threshold_use_max_evalue: bool,
    genome_blast_threshold_max_evalue: String,
    genome_blast_threshold_use_min_identity_percent: bool,
    genome_blast_threshold_min_identity_percent: String,
    genome_blast_threshold_use_min_query_coverage_percent: bool,
    genome_blast_threshold_min_query_coverage_percent: String,
    genome_blast_threshold_use_min_alignment_length_bp: bool,
    genome_blast_threshold_min_alignment_length_bp: String,
    genome_blast_threshold_use_min_bit_score: bool,
    genome_blast_threshold_min_bit_score: String,
    genome_blast_threshold_unique_best_hit: bool,
    genome_blast_options_json: String,
    genome_blast_task: Option<GenomeBlastTask>,
    genome_blast_progress_fraction: Option<f32>,
    genome_blast_progress_label: String,
    genome_blast_results: Vec<GenomeBlastQueryResult>,
    genome_blast_selected_result: usize,
    genome_blast_status: String,
    show_genome_bed_track_dialog: bool,
    show_gibson_dialog: bool,
    show_arrangement_gel_preview_dialog: bool,
    show_rack_labels_preview_dialog: bool,
    show_rack_dialog: bool,
    show_place_arrangement_rack_dialog: bool,
    show_pcr_design_dialog: bool,
    show_sequencing_confirmation_dialog: bool,
    show_planning_dialog: bool,
    show_routine_assistant_dialog: bool,
    show_agent_assistant_dialog: bool,
    show_jaspar_expert_dialog: bool,
    pcr_design_seq_id: String,
    sequencing_confirmation_seq_id: String,
    jaspar_expert_filter: String,
    jaspar_expert_selected_motif_id: String,
    jaspar_expert_random_length_bp: String,
    jaspar_expert_random_seed: String,
    jaspar_expert_fetch_remote_metadata: bool,
    jaspar_expert_status: String,
    jaspar_catalog_report: Option<JasparCatalogReport>,
    jaspar_expert_view: Option<JasparEntryExpertView>,
    gibson_destination_seq_id: String,
    gibson_opening_mode: GibsonUiOpeningMode,
    gibson_opening_start_0based: String,
    gibson_opening_end_0based_exclusive: String,
    gibson_insert_seq_id: String,
    gibson_insert_orientation: GibsonUiInsertOrientation,
    gibson_extra_inserts: Vec<GibsonUiInsertRow>,
    gibson_overlap_bp_min: String,
    gibson_overlap_bp_max: String,
    gibson_minimum_overlap_tm_celsius: String,
    gibson_priming_tm_min_celsius: String,
    gibson_priming_tm_max_celsius: String,
    gibson_priming_length_min_bp: String,
    gibson_priming_length_max_bp: String,
    gibson_unique_restriction_site_enzyme_name: String,
    gibson_output_id_hint: String,
    gibson_show_all_unique_cutters: bool,
    gibson_status: String,
    gibson_preview_output: Option<GibsonAssemblyPreview>,
    gibson_preview_svg_uri: String,
    arrangement_gel_preview: ArrangementGelPreviewState,
    rack_labels_preview: RackLabelsPreviewState,
    rack_view_rack_id: String,
    rack_view_status: String,
    rack_profile_editor_kind: RackProfileKind,
    rack_authoring_template_editor: RackAuthoringTemplate,
    rack_fill_direction_editor: RackFillDirection,
    rack_custom_profile_rows: String,
    rack_custom_profile_columns: String,
    rack_blocked_coordinates_text: String,
    rack_label_sheet_preset: RackLabelSheetPreset,
    rack_carrier_label_preset: RackCarrierLabelPreset,
    rack_physical_template_kind: RackPhysicalTemplateKind,
    rack_view_scroll_offset: Vec2,
    rack_view_selected_coordinates: BTreeSet<String>,
    rack_view_selected_arrangement_ids: BTreeSet<String>,
    rack_view_drag_state: Option<RackDragState>,
    rack_view_hover_target_coordinate: Option<String>,
    rack_view_recent_drop_ghost: Option<RackDropGhostState>,
    rack_help_strip_collapsed: bool,
    rack_help_strip_pinned_open: bool,
    rack_help_strip_successful_move_count: u32,
    rack_help_strip_auto_minimized: bool,
    place_arrangement_source_id: String,
    place_arrangement_target_rack_id: String,
    place_arrangement_status: String,
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
    genome_track_subscription_filter: String,
    genome_track_autosync_status: String,
    tracked_autosync_last_op_count: Option<usize>,
    genome_blast_import_track_name: String,
    genome_blast_import_clear_existing: bool,
    show_command_palette_dialog: bool,
    command_palette_query: String,
    command_palette_selected: usize,
    command_palette_focus_query: bool,
    external_services_ui: ExternalServicesUiState,
    clawbio_panel: clawbio_ui::ClawBioPanelState,
    show_jobs_panel: bool,
    history_ui: HistoryUiState,
    hover_status_name: String,
    app_status: String,
    routine_assistant_stage: RoutineAssistantStage,
    routine_assistant_goal: String,
    routine_assistant_query: String,
    routine_assistant_candidates: Vec<CloningRoutineCatalogRow>,
    routine_assistant_preference_context: Option<RoutinePreferenceContextRecord>,
    routine_assistant_macro_suggestions: Vec<MacroTemplateSuggestion>,
    routine_assistant_selected_routine_id: String,
    routine_assistant_compare_routine_id: String,
    routine_assistant_bindings: BTreeMap<String, String>,
    routine_assistant_disambiguation_answers: BTreeMap<String, String>,
    routine_assistant_explain_output: Option<serde_json::Value>,
    routine_assistant_compare_output: Option<serde_json::Value>,
    routine_assistant_preflight_output: Option<serde_json::Value>,
    routine_assistant_execute_output: Option<serde_json::Value>,
    routine_assistant_status: String,
    routine_assistant_decision_trace: Option<RoutineDecisionTrace>,
    routine_assistant_trace_counter: u64,
    planning_profile_global_json: String,
    planning_profile_overlay_json: String,
    planning_profile_project_json: String,
    planning_objective_json: String,
    planning_sync_payload_json: String,
    planning_sync_source: String,
    planning_sync_confidence: String,
    planning_sync_snapshot_id: String,
    planning_sync_message: String,
    planning_rejection_reason: String,
    planning_suggestions_filter: Option<PlanningSuggestionStatus>,
    planning_host_filter: String,
    planning_host_selected_id: String,
    planning_helper_filter: String,
    planning_helper_selected_id: String,
    planning_status: String,
    next_background_job_id: u64,
    job_event_log: Vec<BackgroundJobEvent>,
    next_retry_snapshot_id: u64,
    retry_argument_snapshots: Vec<BackgroundJobRetrySnapshot>,
    next_retry_cleanup_audit_id: u64,
    retry_snapshot_cleanup_audit: Vec<RetrySnapshotCleanupAuditEntry>,
    retry_snapshot_retain_count: usize,
    retry_snapshot_kind_filter: RetrySnapshotKindFilter,
    retry_snapshot_text_filter: String,
    retry_cleanup_audit_retain_count: usize,
    retry_cleanup_audit_action_filter: RetryCleanupAuditActionFilter,
    retry_cleanup_audit_text_filter: String,
    retry_cleanup_audit_pending_clear_all: bool,
    retry_cleanup_audit_clear_confirm_text: String,
    retry_snapshot_pending_cleanup_action: Option<RetrySnapshotPendingCleanupAction>,
    retry_snapshot_cleanup_confirm_text: String,
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
    agent_model_discovery_failed_source_key: String,
    agent_model_discovery_task: Option<AgentModelDiscoveryTask>,
    agent_prompt_template_id: String,
    agent_prompt: String,
    agent_include_state_summary: bool,
    agent_allow_auto_exec: bool,
    agent_status: String,
    agent_preflight_output: Option<AgentSystemPreflight>,
    agent_task: Option<AgentAskTask>,
    agent_last_invocation: Option<AgentInvocationOutcome>,
    agent_execution_log: Vec<AgentCommandExecutionRecord>,
}

#[derive(Clone)]
enum ProjectAction {
    New,
    Open,
    OpenPath(String),
    OpenTutorialChapter(String),
    Close,
    Quit,
}

#[derive(Clone)]
struct TutorialProjectEntry {
    chapter_id: String,
    chapter_order: usize,
    chapter_title: String,
    chapter_summary: String,
    chapter_guide_path: Option<String>,
    group_label: Option<String>,
    group_order: Option<usize>,
    group_position: Option<usize>,
    decimal_id: Option<String>,
    review_status: Option<String>,
    codex_reviewed_at: Option<String>,
    human_reviewed_at: Option<String>,
    human_reviewer: Option<String>,
    review_stale: bool,
    tier: TutorialTier,
    example: WorkflowExample,
    repo_root: PathBuf,
}

#[derive(Clone)]
struct TutorialCatalogRuntimePaths {
    repo_root: PathBuf,
    manifest_path: PathBuf,
    examples_dir: PathBuf,
}

#[derive(Clone, Debug)]
struct HelpTutorialDocEntry {
    title: String,
    path: String,
    summary: String,
    audiences: Vec<String>,
    group_label: Option<String>,
    group_order: Option<usize>,
    group_position: Option<usize>,
    decimal_id: Option<String>,
    review_status: Option<String>,
    codex_reviewed_at: Option<String>,
    human_reviewed_at: Option<String>,
    human_reviewer: Option<String>,
    review_stale: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum HelpDoc {
    Gui,
    Cli,
    AgentInterface,
    ReviewerPreview,
    Shell,
    Tutorial,
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

#[derive(Clone, Debug)]
struct RenderedHelpMarkdownEntry {
    source_hash: u64,
    rendered: Arc<str>,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum ConfigurationTab {
    ExternalApplications,
    Graphics,
    Language,
}

#[derive(Clone, Debug, Default)]
struct HistoryUiState {
    show_panel: bool,
}

#[derive(Clone, Debug)]
struct ExternalServicesUiState {
    show_panel: bool,
    selected_provider: String,
    selected_service_kind: String,
    request_json: String,
    status: String,
    provider_catalog_output: Option<serde_json::Value>,
    provider_doctor_output: Option<serde_json::Value>,
    preflight_output: Option<serde_json::Value>,
    quote_output: Option<serde_json::Value>,
    selected_quote_payload: usize,
    quote_output_dir: String,
}

impl Default for ExternalServicesUiState {
    fn default() -> Self {
        Self {
            show_panel: false,
            selected_provider: "metabion".to_string(),
            selected_service_kind: "dna_oligo_single_tube".to_string(),
            request_json: String::new(),
            status: String::new(),
            provider_catalog_output: None,
            provider_doctor_output: None,
            preflight_output: None,
            quote_output: None,
            selected_quote_payload: 0,
            quote_output_dir: "artifacts/external_services/metabion_dna_oligo_single_tube_handoff"
                .to_string(),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
enum BackgroundJobKind {
    PrepareGenome,
    BlastGenome,
    TrackImport,
    OpenTutorialProject,
    AgentAssist,
}

impl BackgroundJobKind {
    fn label(self) -> &'static str {
        match self {
            Self::PrepareGenome => "PrepareGenome",
            Self::BlastGenome => "BlastGenome",
            Self::TrackImport => "TrackImport",
            Self::OpenTutorialProject => "OpenTutorialProject",
            Self::AgentAssist => "AgentAssist",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum RetrySnapshotKindFilter {
    All,
    PrepareGenome,
    BlastGenome,
    TrackImport,
    OpenTutorialProject,
    AgentAssist,
}

impl RetrySnapshotKindFilter {
    fn label(self) -> &'static str {
        match self {
            Self::All => "all",
            Self::PrepareGenome => "PrepareGenome",
            Self::BlastGenome => "BlastGenome",
            Self::TrackImport => "TrackImport",
            Self::OpenTutorialProject => "OpenTutorialProject",
            Self::AgentAssist => "AgentAssist",
        }
    }

    fn export_token(self) -> &'static str {
        match self {
            Self::All => "all",
            Self::PrepareGenome => "prepare_genome",
            Self::BlastGenome => "blast_genome",
            Self::TrackImport => "track_import",
            Self::OpenTutorialProject => "open_tutorial_project",
            Self::AgentAssist => "agent_assist",
        }
    }

    fn matches_kind(self, kind: BackgroundJobKind) -> bool {
        match self {
            Self::All => true,
            Self::PrepareGenome => kind == BackgroundJobKind::PrepareGenome,
            Self::BlastGenome => kind == BackgroundJobKind::BlastGenome,
            Self::TrackImport => kind == BackgroundJobKind::TrackImport,
            Self::OpenTutorialProject => kind == BackgroundJobKind::OpenTutorialProject,
            Self::AgentAssist => kind == BackgroundJobKind::AgentAssist,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum RetryCleanupAuditActionFilter {
    All,
    DeleteFiltered,
    ArchiveDeleteFiltered,
    PruneOldest,
    ClearAll,
}

impl RetryCleanupAuditActionFilter {
    fn label(self) -> &'static str {
        match self {
            Self::All => "all",
            Self::DeleteFiltered => "delete_filtered",
            Self::ArchiveDeleteFiltered => "archive_delete_filtered",
            Self::PruneOldest => "prune_oldest",
            Self::ClearAll => "clear_all",
        }
    }

    fn matches_action(self, action: &str) -> bool {
        match self {
            Self::All => true,
            Self::DeleteFiltered => action.eq_ignore_ascii_case("delete_filtered"),
            Self::ArchiveDeleteFiltered => action.eq_ignore_ascii_case("archive_delete_filtered"),
            Self::PruneOldest => action.eq_ignore_ascii_case("prune_oldest"),
            Self::ClearAll => action.eq_ignore_ascii_case("clear_all"),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum RetrySnapshotPendingCleanupAction {
    DeleteFiltered,
    ArchiveAndDeleteFiltered,
}

impl RetrySnapshotPendingCleanupAction {
    fn label(self) -> &'static str {
        match self {
            Self::DeleteFiltered => "Delete filtered",
            Self::ArchiveAndDeleteFiltered => "Archive & delete filtered",
        }
    }

    fn confirm_phrase(self, target_count: usize) -> String {
        match self {
            Self::DeleteFiltered => format!("delete {target_count}"),
            Self::ArchiveAndDeleteFiltered => format!("archive-delete {target_count}"),
        }
    }
}

#[derive(Clone, Debug, Default)]
struct RetrySnapshotDryRunDiff {
    removed: Vec<BackgroundJobRetrySnapshot>,
    retained: Vec<BackgroundJobRetrySnapshot>,
}

impl RetrySnapshotDryRunDiff {
    fn summary_line(&self) -> String {
        format!(
            "Dry-run diff: remove {} snapshot(s), retain {}",
            self.removed.len(),
            self.retained.len()
        )
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
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

#[derive(Clone, Debug, Serialize, Deserialize)]
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

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default)]
struct BackgroundJobRetrySnapshot {
    snapshot_id: u64,
    kind: BackgroundJobKind,
    origin: String,
    captured_at_unix_ms: u128,
    arguments: serde_json::Value,
}

impl Default for BackgroundJobRetrySnapshot {
    fn default() -> Self {
        Self {
            snapshot_id: 0,
            kind: BackgroundJobKind::PrepareGenome,
            origin: String::new(),
            captured_at_unix_ms: 0,
            arguments: serde_json::Value::Null,
        }
    }
}

impl BackgroundJobRetrySnapshot {
    fn to_line(&self) -> String {
        let args = serde_json::to_string(&self.arguments)
            .unwrap_or_else(|_| "{\"error\":\"could not format retry args\"}".to_string());
        let mut compact = args.replace('\n', " ");
        if compact.len() > 220 {
            compact.truncate(220);
            compact.push_str("...");
        }
        format!(
            "[{} retry#{} @{} from {}] {}",
            self.kind.label(),
            self.snapshot_id,
            self.captured_at_unix_ms,
            if self.origin.trim().is_empty() {
                "-"
            } else {
                self.origin.trim()
            },
            compact
        )
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default)]
#[derive(Default)]
struct RetrySnapshotCleanupAuditEntry {
    audit_id: u64,
    action: String,
    performed_at_unix_ms: u128,
    target_count: usize,
    removed_count: usize,
    retained_before: usize,
    retained_after: usize,
    filter_kind: String,
    filter_text: String,
    archive_path: Option<String>,
    summary: String,
}

impl RetrySnapshotCleanupAuditEntry {
    fn to_line(&self) -> String {
        let mut summary = self.summary.replace('\n', " ");
        if summary.len() > 180 {
            summary.truncate(180);
            summary.push_str("...");
        }
        let filter_text = if self.filter_text.trim().is_empty() {
            "-".to_string()
        } else {
            self.filter_text.trim().to_string()
        };
        let archive_token = self
            .archive_path
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .unwrap_or("-");
        format!(
            "[audit#{} {} @{} removed={} retained {}->{} target={} filter(kind={}, text={}) archive={}] {}",
            self.audit_id,
            self.action,
            self.performed_at_unix_ms,
            self.removed_count,
            self.retained_before,
            self.retained_after,
            self.target_count,
            self.filter_kind,
            filter_text,
            archive_token,
            summary
        )
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
struct PersistedBackgroundJobHistory {
    schema: String,
    next_background_job_id: u64,
    events: Vec<BackgroundJobEvent>,
    next_retry_snapshot_id: u64,
    retry_snapshots: Vec<BackgroundJobRetrySnapshot>,
    next_retry_cleanup_audit_id: u64,
    retry_cleanup_audit: Vec<RetrySnapshotCleanupAuditEntry>,
}

impl Default for PersistedBackgroundJobHistory {
    fn default() -> Self {
        Self {
            schema: BACKGROUND_JOB_HISTORY_SCHEMA.to_string(),
            next_background_job_id: 1,
            events: vec![],
            next_retry_snapshot_id: 1,
            retry_snapshots: vec![],
            next_retry_cleanup_audit_id: 1,
            retry_cleanup_audit: vec![],
        }
    }
}

struct GenomePrepareTask {
    job_id: u64,
    started: Instant,
    cancel_requested: Arc<AtomicBool>,
    timeout_seconds: Option<u64>,
    mode: GenomePrepareLaunchMode,
    genome_id: String,
    scope: GenomeDialogScope,
    catalog_path: String,
    cache_dir: String,
    receiver: mpsc::Receiver<GenomePrepareTaskMessage>,
}

struct TutorialProjectTask {
    job_id: u64,
    started: Instant,
    cancel_requested: Arc<AtomicBool>,
    chapter_id: String,
    chapter_title: String,
    receiver: mpsc::Receiver<TutorialProjectTaskMessage>,
}

#[derive(Debug, Clone)]
struct TutorialProjectTaskProgress {
    chapter_id: String,
    chapter_title: String,
    phase: String,
    item: String,
    percent: Option<f32>,
}

#[derive(Debug, Clone)]
struct TutorialProjectOpenOutcome {
    chapter_id: String,
    chapter_title: String,
    project_path: String,
    guide_path: String,
    guide_summary: String,
}

enum TutorialProjectTaskMessage {
    Progress {
        job_id: u64,
        progress: TutorialProjectTaskProgress,
    },
    Done {
        job_id: u64,
        result: Result<TutorialProjectOpenOutcome, String>,
    },
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

#[derive(Debug, Clone)]
struct PreparedGenomeReinstallRequest {
    genome_id: String,
    scope: GenomeDialogScope,
    catalog_path: String,
    cache_dir: String,
    dialog_host: PreparedGenomeReinstallDialogHost,
}

#[derive(Debug, Clone)]
struct PendingEnsemblCatalogUpdateDialog {
    scope: GenomeDialogScope,
    catalog_path: String,
    preview: EnsemblCatalogUpdatePreview,
    output_catalog_path: String,
}

#[derive(Debug, Clone)]
struct PendingEnsemblInstallableGenomeDialog {
    scope: GenomeDialogScope,
    collection_filter: String,
    filter: String,
    report: EnsemblInstallableGenomeCatalog,
}

#[derive(Debug, Clone)]
struct PendingEnsemblQuickInstallDialog {
    scope: GenomeDialogScope,
    collection: String,
    species_dir: String,
    preview: EnsemblQuickInstallPreview,
    genome_id: String,
    output_catalog_path: String,
}

#[derive(Debug, Clone)]
struct PendingPreparedGenomeRemovalRequest {
    genome_id: String,
    scope: GenomeDialogScope,
    catalog_path: String,
    cache_dir: Option<String>,
    install_dir: String,
}

#[derive(Debug, Clone)]
struct PendingGenomeCatalogEntryRemovalRequest {
    genome_id: String,
    scope: GenomeDialogScope,
    catalog_path: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CacheCleanupScope {
    References,
    Helpers,
    Both,
}

impl CacheCleanupScope {
    fn label(self) -> &'static str {
        match self {
            Self::References => "References",
            Self::Helpers => "Helpers",
            Self::Both => "Both",
        }
    }

    fn roots(self, reference_cache_dir: &str, helper_cache_dir: &str) -> Vec<String> {
        let reference = if reference_cache_dir.trim().is_empty() {
            configured_reference_genome_cache_dir()
        } else {
            reference_cache_dir.trim().to_string()
        };
        let helper = if helper_cache_dir.trim().is_empty() {
            configured_helper_genome_cache_dir()
        } else {
            helper_cache_dir.trim().to_string()
        };
        match self {
            Self::References => vec![reference],
            Self::Helpers => vec![helper],
            Self::Both => vec![reference, helper],
        }
    }
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

struct DbSnpFetchTask {
    started: Instant,
    receiver: mpsc::Receiver<DbSnpFetchTaskMessage>,
}

enum DbSnpFetchTaskMessage {
    Progress(DbSnpFetchProgress),
    Done(Result<OpResult, EngineError>),
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum GenomeBlastSourceMode {
    Manual,
    ProjectSequence,
    ProjectPool,
}

impl GenomeBlastSourceMode {
    fn as_str(self) -> &'static str {
        match self {
            Self::Manual => "manual",
            Self::ProjectSequence => "project_sequence",
            Self::ProjectPool => "project_pool",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
enum GenomeBlastOptionsPreset {
    #[default]
    None,
    StrictIdentityCoverage,
    UniqueBestHit,
    HighStringency,
}

impl GenomeBlastOptionsPreset {
    fn label(self) -> &'static str {
        match self {
            Self::None => "None",
            Self::StrictIdentityCoverage => "Strict identity+coverage",
            Self::UniqueBestHit => "Unique best hit",
            Self::HighStringency => "High stringency",
        }
    }

    fn as_request_override_json(self) -> Option<serde_json::Value> {
        match self {
            Self::None => None,
            Self::StrictIdentityCoverage => Some(serde_json::json!({
                "thresholds": {
                    "min_identity_percent": 97.0,
                    "min_query_coverage_percent": 80.0
                }
            })),
            Self::UniqueBestHit => Some(serde_json::json!({
                "thresholds": {
                    "unique_best_hit": true
                }
            })),
            Self::HighStringency => Some(serde_json::json!({
                "task": "blastn",
                "max_hits": 50,
                "thresholds": {
                    "max_evalue": 1e-6,
                    "min_identity_percent": 99.0,
                    "min_query_coverage_percent": 90.0,
                    "min_bit_score": 50.0
                }
            })),
        }
    }
}

struct GenomeBlastTask {
    job_id: u64,
    started: Instant,
    cancel_requested: Arc<AtomicBool>,
    receiver: mpsc::Receiver<GenomeBlastTaskMessage>,
}

enum GenomeBlastTaskMessage {
    Progress {
        job_id: u64,
        done_queries: usize,
        total_queries: usize,
        current_query_label: String,
    },
    Status {
        job_id: u64,
        status: String,
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

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum LineageNodeKind {
    Sequence,
    Arrangement,
    Macro,
    Analysis,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum ProjectOverviewTarget {
    Lineage,
    Containers,
    Arrangements,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
struct ProjectOverviewMetric {
    label: &'static str,
    count: usize,
    target: ProjectOverviewTarget,
    hover: &'static str,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum LineageCopyPayloadKind {
    NodeId,
    PrimaryId,
    DisplayLabel,
    RowSummary,
}

impl LineageCopyPayloadKind {
    fn menu_label(self) -> &'static str {
        match self {
            Self::NodeId => "Copy node ID",
            Self::PrimaryId => "Copy sequence/report ID",
            Self::DisplayLabel => "Copy display label",
            Self::RowSummary => "Copy row summary",
        }
    }

    fn hover_text(self) -> &'static str {
        match self {
            Self::NodeId => "Copy the stable lineage node identifier",
            Self::PrimaryId => "Copy the most relevant sequence, report, macro, or arrangement ID",
            Self::DisplayLabel => "Copy the label shown for this lineage row",
            Self::RowSummary => "Copy a compact multi-line summary of this lineage row",
        }
    }

    fn status_label(self) -> &'static str {
        match self {
            Self::NodeId => "node ID",
            Self::PrimaryId => "sequence/report ID",
            Self::DisplayLabel => "display label",
            Self::RowSummary => "row summary",
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum LineageAnalysisKind {
    Dotplot,
    FlexibilityTrack,
    RnaReadInterpretation,
    PrimerDesign,
    QpcrDesign,
    RestrictionCloningPcrHandoff,
    ProteinDerivation,
    ReverseTranslation,
    ConstructReasoning,
    UniprotProjection,
    SequencingConfirmation,
}

impl LineageAnalysisKind {
    fn as_str(self) -> &'static str {
        match self {
            Self::Dotplot => "dotplot",
            Self::FlexibilityTrack => "flexibility_track",
            Self::RnaReadInterpretation => "rna_read_interpretation",
            Self::PrimerDesign => "primer_design",
            Self::QpcrDesign => "qpcr_design",
            Self::RestrictionCloningPcrHandoff => "restriction_cloning_pcr_handoff",
            Self::ProteinDerivation => "protein_derivation",
            Self::ReverseTranslation => "reverse_translation",
            Self::ConstructReasoning => "construct_reasoning",
            Self::UniprotProjection => "uniprot_projection",
            Self::SequencingConfirmation => "sequencing_confirmation",
        }
    }
}

#[derive(Clone, Debug)]
enum LineageRetrievalDescriptor {
    GenomeGene {
        scope: GenomeDialogScope,
        genome_id: String,
        gene_query: String,
        occurrence: Option<usize>,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
    },
    GenomeRegion {
        scope: GenomeDialogScope,
        genome_id: String,
        chromosome: String,
        start_1based: usize,
        end_1based: usize,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
    },
    GenBankAccession {
        accession: String,
        as_id: Option<String>,
    },
    DbSnpRegion {
        genome_id: String,
        rs_id: String,
        flank_bp: usize,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
    },
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
    genome_anchor_summary: Option<SequenceGenomeAnchorSummary>,
    genome_anchor_display: Option<String>,
    is_full_genome_sequence: bool,
    retrieval_descriptor: Option<LineageRetrievalDescriptor>,
    analysis_kind: Option<LineageAnalysisKind>,
    analysis_artifact_id: Option<String>,
    analysis_reference_seq_id: Option<String>,
    analysis_mode: Option<String>,
    analysis_status: Option<String>,
    analysis_point_count: Option<usize>,
    analysis_bin_count: Option<usize>,
    analysis_read_count: Option<usize>,
    analysis_trace_count: Option<usize>,
    analysis_target_count: Option<usize>,
    analysis_variant_count: Option<usize>,
    macro_instance_id: Option<String>,
    macro_routine_id: Option<String>,
    macro_template_name: Option<String>,
    macro_status: Option<String>,
    macro_status_message: Option<String>,
    macro_op_ids: Vec<String>,
    macro_inputs: Vec<LineageMacroPortBinding>,
    macro_outputs: Vec<LineageMacroPortBinding>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(default)]
#[derive(Default)]
struct PersistedLineageNodeGroup {
    group_id: String,
    label: String,
    representative_node_id: String,
    member_node_ids: Vec<String>,
    collapsed: bool,
}

#[derive(Clone)]
struct LineageTableEntry {
    row: LineageRow,
    indent_level: usize,
    group_id: Option<String>,
    group_label: Option<String>,
    is_group_representative: bool,
    hidden_group_member_count: usize,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct LineageGroupHiddenOpBadge {
    total_ops: usize,
    families: Vec<(String, usize)>,
}

#[derive(Clone, Debug)]
struct AnchorProvenanceSnapshot {
    catalog_path: Option<String>,
    cache_dir: Option<String>,
    recorded_at_unix_ms: u128,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct GenomeLengthCacheKey {
    catalog_path: String,
    genome_id: String,
    cache_dir: Option<String>,
}

#[derive(Clone)]
struct ContainerRow {
    container_id: String,
    kind: String,
    declared_contents_exclusive: bool,
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
    default_rack_id: Option<String>,
}

#[derive(Clone)]
struct RackRow {
    rack_id: String,
    name: String,
    profile: String,
    occupied_positions: usize,
    arrangement_ids: Vec<String>,
}

#[derive(Clone)]
struct CloningPatternCatalogEntry {
    label: String,
    path: String,
    is_file: bool,
    children: Vec<CloningPatternCatalogEntry>,
}

#[derive(Clone, Deserialize, Default)]
#[serde(default)]
struct CloningRoutineDifferenceAxisRow {
    axis: String,
    value: String,
}

#[derive(Clone, Deserialize, Default)]
#[serde(default)]
struct CloningRoutineCatalogRow {
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
    difference_matrix: Vec<CloningRoutineDifferenceAxisRow>,
    disambiguation_questions: Vec<String>,
    failure_modes: Vec<String>,
    details_url: Option<String>,
    template_name: String,
    template_path: Option<String>,
    input_ports: Vec<serde_json::Value>,
    output_ports: Vec<serde_json::Value>,
    estimated_time_hours: Option<f64>,
    estimated_cost: Option<f64>,
    local_fit_score: Option<f64>,
    composite_meta_score: Option<f64>,
    planning_estimate: Option<PlanningEstimate>,
}

#[derive(Clone, Debug)]
struct RoutineAssistantBoundSequenceTopology {
    port_id: String,
    seq_id: String,
    circular: bool,
    length_bp: usize,
}

#[derive(Clone, Copy, PartialEq, Eq, Default)]
enum RoutineAssistantStage {
    #[default]
    GoalAndCandidates,
    Compare,
    Parameters,
    Preflight,
    ExecuteAndExport,
}

impl RoutineAssistantStage {
    fn label(self) -> &'static str {
        match self {
            Self::GoalAndCandidates => "1. Goal",
            Self::Compare => "2. Compare",
            Self::Parameters => "3. Parameters",
            Self::Preflight => "4. Preflight",
            Self::ExecuteAndExport => "5. Run/Export",
        }
    }
}

#[derive(Debug, Clone, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
struct PlanningSyncSuggestionPayloadUi {
    profile_patch: Option<PlanningProfile>,
    objective_patch: Option<PlanningObjective>,
}

#[derive(Clone, Copy)]
enum CommandPaletteAction {
    NewProject,
    OpenProject,
    SaveProject,
    NewSequence,
    NewSequenceFromClipboard,
    OpenSequence,
    OpenUniprot,
    OpenGenbank,
    OpenConfiguration,
    OpenExternalServices,
    OpenGibson,
    OpenPlanning,
    OpenRoutineAssistant,
    UiIntent(UiIntentTarget),
    OpenGuiManual,
    OpenCliManual,
    OpenAgentInterfaceManual,
    OpenReviewerPreviewManual,
    OpenShellManual,
    ExportLineageSvg,
    ExportLabAssistantReport,
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

impl CommandPaletteEntry {
    fn ui_intent(target: UiIntentTarget) -> Self {
        Self {
            title: target.discoverability_title().to_string(),
            detail: target.discoverability_detail().to_string(),
            keywords: format!(
                "{} {} {}",
                target.discoverability_keywords(),
                target.menu_path().to_ascii_lowercase(),
                target.as_str()
            ),
            action: CommandPaletteAction::UiIntent(target),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
struct PersistedLineageGraphWorkspace {
    graph_view: Option<bool>,
    main_split_fraction: Option<f32>,
    container_arrangement_split_fraction: Option<f32>,
    zoom: f32,
    graph_area_height: f32,
    container_area_height: f32,
    scroll_offset: [f32; 2],
    compact_labels: bool,
    node_offsets: HashMap<String, [f32; 2]>,
    node_groups: Vec<PersistedLineageNodeGroup>,
}

impl Default for PersistedLineageGraphWorkspace {
    fn default() -> Self {
        Self {
            graph_view: Some(true),
            main_split_fraction: Some(DEFAULT_LINEAGE_MAIN_SPLIT_FRACTION),
            container_arrangement_split_fraction: Some(
                DEFAULT_LINEAGE_CONTAINER_ARRANGEMENT_SPLIT_FRACTION,
            ),
            zoom: 1.0,
            graph_area_height: 420.0,
            container_area_height: 220.0,
            scroll_offset: [0.0, 0.0],
            compact_labels: true,
            node_offsets: HashMap::new(),
            node_groups: vec![],
        }
    }
}

impl Default for GENtleApp {
    fn default() -> Self {
        Self {
            engine: Arc::new(RwLock::new(GentleEngine::new())),
            new_windows: vec![],
            windows: HashMap::new(),
            detached_auxiliary_window_hosts: HashMap::new(),
            windows_to_close: Arc::new(RwLock::new(vec![])),
            pending_focus_viewports: vec![],
            pending_window_initial_positions: HashMap::new(),
            viewport_id_counter: 0,
            update_has_run_before: false,
            splash_started_at: Instant::now(),
            splash_dismissed: false,
            show_about_dialog: false,
            show_help_dialog: false,
            help_doc: HelpDoc::Gui,
            help_markdown_cache: CommonMarkCache::default(),
            help_rendered_markdown_cache: HashMap::new(),
            help_gui_markdown: GUI_MANUAL_MD.to_string(),
            help_cli_markdown: CLI_MANUAL_MD.to_string(),
            help_agent_interface_markdown: AGENT_INTERFACE_MD.to_string(),
            help_reviewer_preview_markdown: REVIEWER_PREVIEW_MD.to_string(),
            help_shell_markdown: Self::generate_shell_help_markdown(),
            help_tutorial_markdown: String::new(),
            help_tutorial_title: "Tutorial".to_string(),
            help_tutorial_entries: vec![],
            help_tutorial_selected: 0,
            help_shell_interface: ShellHelpInterface::GuiShell,
            help_search_query: String::new(),
            help_search_matches: vec![],
            help_search_selected: 0,
            help_focus_search_box: false,
            help_last_content_width: 0.0,
            help_selectable_text_mode: false,
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
            ui_language: UiLanguage::default(),
            configuration_ui_language: UiLanguage::default(),
            configuration_language_dirty: false,
            i18n: I18n::default(),
            configuration_status: String::new(),
            window_backdrop_path_status_cache: HashMap::new(),
            current_project_path: None,
            recent_project_paths: vec![],
            lineage_graph_view: true,
            lineage_main_split_fraction: DEFAULT_LINEAGE_MAIN_SPLIT_FRACTION,
            lineage_graph_zoom: 1.0,
            lineage_graph_area_height: 420.0,
            lineage_container_area_height: 220.0,
            lineage_container_arrangement_split_fraction:
                DEFAULT_LINEAGE_CONTAINER_ARRANGEMENT_SPLIT_FRACTION,
            lineage_graph_scroll_offset: Vec2::ZERO,
            lineage_graph_pan_origin: None,
            lineage_graph_compact_labels: true,
            lineage_main_split_drag_origin: None,
            lineage_main_split_drag_start_y: None,
            lineage_container_arrangement_split_drag_origin: None,
            lineage_container_arrangement_split_drag_start_y: None,
            lineage_graph_selected_node_id: None,
            lineage_graph_node_offsets: HashMap::new(),
            lineage_graph_drag_origin: None,
            lineage_graph_offsets_synced_stamp: 0,
            lineage_cache_stamp: 0,
            lineage_cache_valid: false,
            lineage_rows: vec![],
            lineage_edges: vec![],
            lineage_op_label_by_id: HashMap::new(),
            lineage_reopenable_gibson_op_ids: HashSet::new(),
            lineage_reopenable_pcr_op_seq_ids: HashMap::new(),
            lineage_node_groups: vec![],
            lineage_group_form_label: String::new(),
            lineage_group_form_members: String::new(),
            lineage_group_form_representative: String::new(),
            lineage_group_form_editing_id: None,
            lineage_group_status: String::new(),
            lineage_group_marked_nodes: HashSet::new(),
            lineage_node_rename_target_id: None,
            lineage_node_rename_text: String::new(),
            lineage_node_remove_target_id: None,
            lineage_containers: vec![],
            lineage_arrangements: vec![],
            lineage_racks: vec![],
            clean_state_fingerprint: 0,
            dirty_cache_stamp: 0,
            last_display_sync_stamp: 0,
            dirty_cache_valid: false,
            dirty_cache_value: false,
            dirty_cache_last_deep_check: Instant::now(),
            last_applied_window_title: String::new(),
            last_native_window_entries: vec![],
            last_native_active_window_key: Some(viewport_native_menu_key(ViewportId::ROOT)),
            native_window_key_to_viewport: HashMap::new(),
            active_window_menu_key: Some(viewport_native_menu_key(ViewportId::ROOT)),
            pending_window_open_timestamps: HashMap::new(),
            pending_viewport_focus_timestamps: HashMap::new(),
            pending_project_action: None,
            pending_app_quit: false,
            pending_prepared_genome_reinstall: None,
            pending_ensembl_catalog_update: None,
            pending_ensembl_installable_genomes: None,
            pending_ensembl_quick_install: None,
            pending_prepared_genome_removal: None,
            pending_genome_catalog_entry_removal: None,
            show_reference_genome_prepare_dialog: false,
            show_reference_genome_retrieve_dialog: false,
            show_reference_genome_blast_dialog: false,
            show_reference_genome_inspector_dialog: false,
            show_cache_cleanup_dialog: false,
            show_new_sequence_dialog: false,
            show_uniprot_dialog: false,
            show_genbank_dialog: false,
            genome_dialog_scope: GenomeDialogScope::Reference,
            reference_genome_catalog_path: DEFAULT_GENOME_CATALOG_PATH.to_string(),
            reference_genome_cache_dir: configured_reference_genome_cache_dir(),
            helper_genome_catalog_path: DEFAULT_HELPER_GENOME_CATALOG_PATH.to_string(),
            helper_genome_cache_dir: configured_helper_genome_cache_dir(),
            genome_catalog_path: DEFAULT_GENOME_CATALOG_PATH.to_string(),
            genome_cache_dir: configured_reference_genome_cache_dir(),
            genome_id: "Human GRCh38 Ensembl 113".to_string(),
            genome_catalog_genomes: vec![],
            genome_catalog_error: String::new(),
            genome_prepare_task: None,
            genome_prepare_progress: None,
            genome_prepare_steps: vec![],
            genome_prepare_eta_baseline: None,
            genome_prepare_failure_recovery: None,
            genome_prepare_status: String::new(),
            tutorial_project_task: None,
            tutorial_project_progress_fraction: None,
            tutorial_project_progress_label: String::new(),
            tutorial_project_status: String::new(),
            cache_cleanup_scope: CacheCleanupScope::References,
            cache_cleanup_reference_cache_dir: configured_reference_genome_cache_dir(),
            cache_cleanup_helper_cache_dir: configured_helper_genome_cache_dir(),
            cache_cleanup_mode: PreparedCacheCleanupMode::DerivedIndexesOnly,
            cache_cleanup_include_orphans: false,
            cache_cleanup_inspection: None,
            cache_cleanup_selected_paths: HashSet::new(),
            cache_cleanup_rebuild_candidates: vec![],
            cache_cleanup_status: String::new(),
            cache_cleanup_confirm_pending: false,
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
            genome_output_id_autofilled: false,
            genome_gene_extract_mode: GenomeGeneExtractMode::Gene,
            genome_gene_promoter_upstream_bp: "1000".to_string(),
            genome_annotation_scope: GenomeAnnotationScope::Core,
            genome_max_annotation_features: String::new(),
            genome_include_genomic_annotation: true,
            genome_retrieve_status: String::new(),
            genome_retrieve_contig_suggestions: vec![],
            new_sequence_text: String::new(),
            new_sequence_output_id: String::new(),
            new_sequence_name: String::new(),
            new_sequence_circular: false,
            new_sequence_status: String::new(),
            uniprot_query: String::new(),
            uniprot_entry_id: String::new(),
            uniprot_swiss_path: String::new(),
            ensembl_protein_query: String::new(),
            ensembl_protein_entry_id: String::new(),
            ensembl_protein_output_id: String::new(),
            protein_feature_key_include: String::new(),
            protein_feature_key_exclude: String::new(),
            uniprot_linked_accession: String::new(),
            uniprot_linked_as_id: String::new(),
            uniprot_map_seq_id: String::new(),
            uniprot_map_projection_id: String::new(),
            uniprot_map_transcript_id: String::new(),
            uniprot_feature_query: String::new(),
            uniprot_feature_transcript_id: String::new(),
            uniprot_feature_query_mode: UniprotFeatureCodingDnaQueryMode::Both,
            uniprot_feature_speed_profile: None,
            uniprot_feature_report: None,
            uniprot_audit_transcript_id: String::new(),
            uniprot_audit_report_id: String::new(),
            uniprot_audit_parity_report_id: String::new(),
            uniprot_audit_report: None,
            uniprot_audit_parity_report: None,
            reverse_translate_protein_seq_id: String::new(),
            reverse_translate_output_id: String::new(),
            reverse_translate_speed_profile: None,
            reverse_translate_speed_mark: None,
            reverse_translate_translation_table: String::new(),
            reverse_translate_target_anneal_tm_c: String::new(),
            reverse_translate_anneal_window_bp: "20".to_string(),
            reverse_translation_report: None,
            protease_digest_names: "Trypsin".to_string(),
            protease_digest_output_prefix: String::new(),
            protease_digest_min_length_aa: "1".to_string(),
            protease_digest_materialize: true,
            protease_digest_report: None,
            protein_handoff_ranking_goal: ProteinToDnaHandoffRankingGoal::BalancedProvenance,
            protein_handoff_graph: None,
            protein_handoff_selected_candidate_id: String::new(),
            protein_handoff_status: String::new(),
            uniprot_status: String::new(),
            genbank_accession: String::new(),
            genbank_as_id: String::new(),
            genbank_status: String::new(),
            dbsnp_rs_id: DEFAULT_DBSNP_TUTORIAL_RS_ID.to_string(),
            dbsnp_genome_id: "Human GRCh38 Ensembl 113".to_string(),
            dbsnp_flank_bp: "3000".to_string(),
            dbsnp_output_id: String::new(),
            dbsnp_status: String::new(),
            dbsnp_fetch_task: None,
            genome_blast_source_mode: GenomeBlastSourceMode::Manual,
            genome_blast_query_manual: String::new(),
            genome_blast_query_seq_id: String::new(),
            genome_blast_query_pool_id: String::new(),
            genome_blast_max_hits: 25,
            genome_blast_task_name: "blastn-short".to_string(),
            genome_blast_options_preset: GenomeBlastOptionsPreset::None,
            genome_blast_threshold_use_max_evalue: false,
            genome_blast_threshold_max_evalue: String::new(),
            genome_blast_threshold_use_min_identity_percent: false,
            genome_blast_threshold_min_identity_percent: String::new(),
            genome_blast_threshold_use_min_query_coverage_percent: false,
            genome_blast_threshold_min_query_coverage_percent: String::new(),
            genome_blast_threshold_use_min_alignment_length_bp: false,
            genome_blast_threshold_min_alignment_length_bp: String::new(),
            genome_blast_threshold_use_min_bit_score: false,
            genome_blast_threshold_min_bit_score: String::new(),
            genome_blast_threshold_unique_best_hit: false,
            genome_blast_options_json: String::new(),
            genome_blast_task: None,
            genome_blast_progress_fraction: None,
            genome_blast_progress_label: String::new(),
            genome_blast_results: vec![],
            genome_blast_selected_result: 0,
            genome_blast_status: String::new(),
            genome_blast_import_track_name: "blast_hits".to_string(),
            genome_blast_import_clear_existing: false,
            show_genome_bed_track_dialog: false,
            show_gibson_dialog: false,
            show_arrangement_gel_preview_dialog: false,
            show_rack_labels_preview_dialog: false,
            show_rack_dialog: false,
            show_place_arrangement_rack_dialog: false,
            show_pcr_design_dialog: false,
            show_sequencing_confirmation_dialog: false,
            show_planning_dialog: false,
            show_routine_assistant_dialog: false,
            show_agent_assistant_dialog: false,
            show_jaspar_expert_dialog: false,
            pcr_design_seq_id: String::new(),
            sequencing_confirmation_seq_id: String::new(),
            jaspar_expert_filter: String::new(),
            jaspar_expert_selected_motif_id: tf_motifs::list_motif_summaries()
                .first()
                .map(|row| row.id.clone())
                .unwrap_or_default(),
            jaspar_expert_random_length_bp: DEFAULT_JASPAR_PRESENTATION_RANDOM_SEQUENCE_LENGTH_BP
                .to_string(),
            jaspar_expert_random_seed: DEFAULT_JASPAR_PRESENTATION_RANDOM_SEED.to_string(),
            jaspar_expert_fetch_remote_metadata: false,
            jaspar_expert_status: String::new(),
            jaspar_catalog_report: None,
            jaspar_expert_view: None,
            gibson_destination_seq_id: String::new(),
            gibson_opening_mode: GibsonUiOpeningMode::DefinedSite,
            gibson_opening_start_0based: String::new(),
            gibson_opening_end_0based_exclusive: String::new(),
            gibson_insert_seq_id: String::new(),
            gibson_insert_orientation: GibsonUiInsertOrientation::Forward,
            gibson_extra_inserts: vec![],
            gibson_overlap_bp_min: "20".to_string(),
            gibson_overlap_bp_max: "40".to_string(),
            gibson_minimum_overlap_tm_celsius: "60.0".to_string(),
            gibson_priming_tm_min_celsius: "58.0".to_string(),
            gibson_priming_tm_max_celsius: "68.0".to_string(),
            gibson_priming_length_min_bp: "18".to_string(),
            gibson_priming_length_max_bp: "35".to_string(),
            gibson_unique_restriction_site_enzyme_name: String::new(),
            gibson_output_id_hint: String::new(),
            gibson_show_all_unique_cutters: false,
            gibson_status: String::new(),
            gibson_preview_output: None,
            gibson_preview_svg_uri: String::new(),
            arrangement_gel_preview: ArrangementGelPreviewState::default(),
            rack_labels_preview: RackLabelsPreviewState::default(),
            rack_view_rack_id: String::new(),
            rack_view_status: String::new(),
            rack_profile_editor_kind: RackProfileKind::SmallTube4x6,
            rack_authoring_template_editor: RackAuthoringTemplate::BenchRows,
            rack_fill_direction_editor: RackFillDirection::RowMajor,
            rack_custom_profile_rows: String::new(),
            rack_custom_profile_columns: String::new(),
            rack_blocked_coordinates_text: String::new(),
            rack_label_sheet_preset: RackLabelSheetPreset::default(),
            rack_carrier_label_preset: RackCarrierLabelPreset::default(),
            rack_physical_template_kind: RackPhysicalTemplateKind::default(),
            rack_view_scroll_offset: Vec2::ZERO,
            rack_view_selected_coordinates: BTreeSet::new(),
            rack_view_selected_arrangement_ids: BTreeSet::new(),
            rack_view_drag_state: None,
            rack_view_hover_target_coordinate: None,
            rack_view_recent_drop_ghost: None,
            rack_help_strip_collapsed: false,
            rack_help_strip_pinned_open: false,
            rack_help_strip_successful_move_count: 0,
            rack_help_strip_auto_minimized: false,
            place_arrangement_source_id: String::new(),
            place_arrangement_target_rack_id: String::new(),
            place_arrangement_status: String::new(),
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
            genome_track_subscription_filter: String::new(),
            genome_track_autosync_status: String::new(),
            tracked_autosync_last_op_count: None,
            show_command_palette_dialog: false,
            command_palette_query: String::new(),
            command_palette_selected: 0,
            command_palette_focus_query: false,
            external_services_ui: ExternalServicesUiState::default(),
            clawbio_panel: clawbio_ui::ClawBioPanelState::default(),
            show_jobs_panel: false,
            history_ui: HistoryUiState::default(),
            hover_status_name: String::new(),
            app_status: String::new(),
            routine_assistant_stage: RoutineAssistantStage::GoalAndCandidates,
            routine_assistant_goal: String::new(),
            routine_assistant_query: String::new(),
            routine_assistant_candidates: vec![],
            routine_assistant_preference_context: None,
            routine_assistant_macro_suggestions: vec![],
            routine_assistant_selected_routine_id: String::new(),
            routine_assistant_compare_routine_id: String::new(),
            routine_assistant_bindings: BTreeMap::new(),
            routine_assistant_disambiguation_answers: BTreeMap::new(),
            routine_assistant_explain_output: None,
            routine_assistant_compare_output: None,
            routine_assistant_preflight_output: None,
            routine_assistant_execute_output: None,
            routine_assistant_status: String::new(),
            routine_assistant_decision_trace: None,
            routine_assistant_trace_counter: 1,
            planning_profile_global_json: String::new(),
            planning_profile_overlay_json: String::new(),
            planning_profile_project_json: String::new(),
            planning_objective_json: String::new(),
            planning_sync_payload_json: String::new(),
            planning_sync_source: "lab_manager".to_string(),
            planning_sync_confidence: String::new(),
            planning_sync_snapshot_id: String::new(),
            planning_sync_message: String::new(),
            planning_rejection_reason: String::new(),
            planning_suggestions_filter: Some(PlanningSuggestionStatus::Pending),
            planning_host_filter: String::new(),
            planning_host_selected_id: String::new(),
            planning_helper_filter: String::new(),
            planning_helper_selected_id: String::new(),
            planning_status: String::new(),
            next_background_job_id: 1,
            job_event_log: vec![],
            next_retry_snapshot_id: 1,
            retry_argument_snapshots: vec![],
            next_retry_cleanup_audit_id: 1,
            retry_snapshot_cleanup_audit: vec![],
            retry_snapshot_retain_count: MAX_BACKGROUND_JOB_RETRY_SNAPSHOTS,
            retry_snapshot_kind_filter: RetrySnapshotKindFilter::All,
            retry_snapshot_text_filter: String::new(),
            retry_cleanup_audit_retain_count: MAX_RETRY_SNAPSHOT_CLEANUP_AUDIT_ENTRIES,
            retry_cleanup_audit_action_filter: RetryCleanupAuditActionFilter::All,
            retry_cleanup_audit_text_filter: String::new(),
            retry_cleanup_audit_pending_clear_all: false,
            retry_cleanup_audit_clear_confirm_text: String::new(),
            retry_snapshot_pending_cleanup_action: None,
            retry_snapshot_cleanup_confirm_text: String::new(),
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
            agent_model_discovery_failed_source_key: String::new(),
            agent_model_discovery_task: None,
            agent_prompt_template_id: AGENT_PROMPT_TEMPLATE_DEFAULT_ID.to_string(),
            agent_prompt: String::new(),
            agent_include_state_summary: true,
            agent_allow_auto_exec: false,
            agent_status: String::new(),
            agent_preflight_output: None,
            agent_task: None,
            agent_last_invocation: None,
            agent_execution_log: vec![],
        }
    }
}

impl GENtleApp {
    fn main_workspace_hosted_window_id() -> egui::Id {
        egui::Id::new("main_workspace_hosted_window")
    }

    fn help_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("GENtle Help Viewport")
    }

    fn hosted_help_window_id() -> egui::Id {
        egui::Id::new(("hosted_help_window", Self::help_viewport_id()))
    }

    fn hosted_routine_assistant_window_id() -> egui::Id {
        egui::Id::new((
            "hosted_routine_assistant_window",
            Self::routine_assistant_viewport_id(),
        ))
    }

    fn hosted_agent_assistant_window_id() -> egui::Id {
        egui::Id::new((
            "hosted_agent_assistant_window",
            Self::agent_assistant_viewport_id(),
        ))
    }

    fn unsaved_changes_dialog_id() -> egui::Id {
        egui::Id::new("unsaved_changes_dialog")
    }

    fn legacy_root_help_layer_ids(title: &str) -> [egui::LayerId; 2] {
        [
            egui::LayerId::new(egui::Order::Middle, egui::Id::new(title.to_string())),
            egui::LayerId::new(egui::Order::Middle, egui::Id::new(Self::help_viewport_id())),
        ]
    }

    fn stale_hosted_window_title_layer_id(title: &str) -> egui::LayerId {
        crate::egui_compat::hosted_window_title_layer_id(title)
    }

    #[cfg(test)]
    fn stale_help_title_layer_id(title: &str) -> egui::LayerId {
        Self::stale_hosted_window_title_layer_id(title)
    }

    fn reset_root_help_areas_if_legacy_layers_visible(ctx: &egui::Context, title: &str) -> bool {
        let stale_hosted_help_layers = Self::legacy_root_help_layer_ids(title);
        if stale_hosted_help_layers
            .iter()
            .any(|layer_id| ctx.memory(|mem| mem.areas().is_visible(layer_id)))
        {
            ctx.memory_mut(|mem| mem.reset_areas());
            true
        } else {
            false
        }
    }

    fn prepare_genome_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("GENtle Prepare Genome Viewport")
    }

    fn retrieve_genome_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("GENtle Retrieve Genome Viewport")
    }

    fn blast_genome_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("GENtle BLAST Genome Viewport")
    }

    fn configuration_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("GENtle Configuration Viewport")
    }

    fn hosted_configuration_window_id() -> egui::Id {
        egui::Id::new((
            "hosted_configuration_window",
            Self::configuration_viewport_id(),
        ))
    }

    fn bed_track_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("GENtle BED Tracks Viewport")
    }

    fn gibson_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("GENtle Gibson Viewport")
    }

    fn arrangement_gel_preview_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("GENtle Arrangement Gel Preview Viewport")
    }

    fn rack_labels_preview_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("GENtle Rack Labels Preview Viewport")
    }

    fn pcr_design_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("GENtle PCR Designer Viewport")
    }

    fn sequencing_confirmation_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("GENtle Sequencing Confirmation Viewport")
    }

    fn planning_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("GENtle Planning Viewport")
    }

    fn routine_assistant_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("GENtle Routine Assistant Viewport")
    }

    fn agent_assistant_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("GENtle Agent Assistant Viewport")
    }

    fn jaspar_expert_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("GENtle JASPAR Expert Viewport")
    }

    fn command_palette_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("GENtle Command Palette Viewport")
    }

    fn background_jobs_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("GENtle Background Jobs Viewport")
    }

    fn external_services_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("GENtle External Services Viewport")
    }

    fn clawbio_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("GENtle ClawBio Viewport")
    }

    fn hosted_clawbio_window_id() -> egui::Id {
        egui::Id::new(("hosted_clawbio_window", Self::clawbio_viewport_id()))
    }

    fn history_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("GENtle Operation History Viewport")
    }

    fn uniprot_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("GENtle UniProt Viewport")
    }

    fn new_sequence_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("GENtle New Sequence Viewport")
    }

    fn genbank_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("GENtle GenBank Viewport")
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
            self.finalize_viewport_focus_probe(viewport_id);
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

    fn tr(&self, key: &str) -> String {
        self.i18n.t(key)
    }

    fn write_persisted_configuration_to_disk(&self) -> std::result::Result<(), String> {
        let path = Self::configuration_store_path();
        if let Some(parent) = path.parent()
            && !parent.as_os_str().is_empty()
        {
            fs::create_dir_all(parent).map_err(|e| {
                format!(
                    "Could not create configuration directory '{}': {e}",
                    parent.display()
                )
            })?;
        }
        let payload = PersistedConfiguration {
            schema_version: APP_CONFIGURATION_SCHEMA_VERSION,
            rnapkin_executable: self.configuration_rnapkin_executable.trim().to_string(),
            makeblastdb_executable: self.configuration_makeblastdb_executable.trim().to_string(),
            blastn_executable: self.configuration_blastn_executable.trim().to_string(),
            bigwig_to_bedgraph_executable: self
                .configuration_bigwig_to_bedgraph_executable
                .trim()
                .to_string(),
            graphics_defaults: self.configuration_graphics.clone(),
            window_backdrops: self.window_backdrops.clone(),
            ui_language: self.ui_language,
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

        if self.recent_project_paths != old_paths
            && let Err(err) = self.write_persisted_configuration_to_disk()
        {
            eprintln!("{err}");
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
        if self.recent_project_paths.len() != len_before
            && let Err(err) = self.write_persisted_configuration_to_disk()
        {
            eprintln!("{err}");
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

    fn clamp_feature_details_font_size(value: f32) -> f32 {
        if value.is_finite() {
            value.clamp(8.0, 24.0)
        } else {
            9.0
        }
    }

    fn clamp_external_feature_label_font_size(value: f32) -> f32 {
        if value.is_finite() {
            value.clamp(8.0, 24.0)
        } else {
            11.0
        }
    }

    fn clamp_external_feature_label_background_opacity(value: f32) -> f32 {
        if value.is_finite() {
            value.clamp(0.0, 1.0)
        } else {
            0.9
        }
    }

    fn clamp_reverse_strand_visual_opacity(value: f32) -> f32 {
        if value.is_finite() {
            value.clamp(0.2, 1.0)
        } else {
            DisplaySettings::default_reverse_strand_visual_opacity()
        }
    }

    fn clamp_sequence_panel_max_text_length_bp(value: usize) -> usize {
        if value == 0 { 0 } else { value.min(5_000_000) }
    }

    fn clamp_linear_helical_phase_offset_bp(value: usize) -> usize {
        value % 10
    }

    fn apply_graphics_settings_to_display(source: &DisplaySettings, target: &mut DisplaySettings) {
        target.show_sequence_panel = source.show_sequence_panel;
        target.show_linear_sequence_panel = source.show_linear_sequence_panel;
        target.sequence_panel_max_text_length_bp =
            Self::clamp_sequence_panel_max_text_length_bp(source.sequence_panel_max_text_length_bp);
        target.auto_hide_sequence_panel_when_linear_bases_visible =
            source.auto_hide_sequence_panel_when_linear_bases_visible;
        target.show_map_panel = source.show_map_panel;
        target.show_features = source.show_features;
        target.show_cds_features = source.show_cds_features;
        target.show_gene_features = source.show_gene_features;
        target.show_mrna_features = source.show_mrna_features;
        target.show_repeat_features = source.show_repeat_features;
        target.show_array_features = source.show_array_features;
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
        target.restriction_enzyme_display_mode = source.restriction_enzyme_display_mode;
        target.preferred_restriction_enzymes =
            crate::dna_display::DnaDisplay::normalize_preferred_restriction_enzymes(
                &source.preferred_restriction_enzymes,
            );
        target.show_gc_contents = source.show_gc_contents;
        target.gc_content_bin_size_bp = source.gc_content_bin_size_bp;
        target.show_open_reading_frames = source.show_open_reading_frames;
        target.show_methylation_sites = source.show_methylation_sites;
        target.feature_details_font_size =
            Self::clamp_feature_details_font_size(source.feature_details_font_size);
        target.linear_external_feature_label_font_size =
            Self::clamp_external_feature_label_font_size(
                source.linear_external_feature_label_font_size,
            );
        target.linear_external_feature_label_background_opacity =
            Self::clamp_external_feature_label_background_opacity(
                source.linear_external_feature_label_background_opacity,
            );
        target.linear_view_start_bp = source.linear_view_start_bp;
        target.linear_view_span_bp = source.linear_view_span_bp;
        target.linear_view_vertical_offset_px = source.linear_view_vertical_offset_px;
        target.linear_sequence_base_text_max_view_span_bp =
            source.linear_sequence_base_text_max_view_span_bp;
        target.linear_sequence_helical_letters_enabled =
            source.linear_sequence_helical_letters_enabled;
        target.linear_sequence_helical_max_view_span_bp =
            source.linear_sequence_helical_max_view_span_bp;
        target.linear_sequence_condensed_max_view_span_bp =
            source.linear_sequence_condensed_max_view_span_bp;
        target.linear_sequence_letter_layout_mode = source.linear_sequence_letter_layout_mode;
        target.linear_sequence_helical_phase_offset_bp = Self::clamp_linear_helical_phase_offset_bp(
            source.linear_sequence_helical_phase_offset_bp,
        );
        target.linear_show_double_strand_bases = source.linear_show_double_strand_bases;
        target.linear_helical_parallel_strands = source.linear_helical_parallel_strands;
        target.linear_hide_backbone_when_sequence_bases_visible =
            source.linear_hide_backbone_when_sequence_bases_visible;
        target.linear_reverse_strand_use_upside_down_letters =
            source.linear_reverse_strand_use_upside_down_letters;
        target.reverse_strand_visual_opacity =
            Self::clamp_reverse_strand_visual_opacity(source.reverse_strand_visual_opacity);
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
        let Some(mut saved) = Self::read_persisted_configuration_from_disk() else {
            return;
        };
        let upgraded = Self::upgrade_persisted_configuration(&mut saved);
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
        self.window_backdrop_path_status_cache.clear();
        self.ui_language = saved.ui_language;
        self.configuration_ui_language = saved.ui_language;
        self.configuration_language_dirty = false;
        self.i18n.set_language(saved.ui_language);
        self.recent_project_paths = Self::normalize_recent_project_paths(saved.recent_projects);
        self.configuration_graphics_dirty = false;
        self.configuration_rnapkin_validation_ok = None;
        self.configuration_rnapkin_validation_message.clear();
        self.configuration_blast_validation_ok = None;
        self.configuration_blast_validation_message.clear();

        if apply_graphics_to_current_project {
            self.apply_configuration_graphics_to_engine_state();
        }
        if upgraded {
            let _ = self.write_persisted_configuration_to_disk();
        }
    }

    fn upgrade_persisted_configuration(saved: &mut PersistedConfiguration) -> bool {
        if saved.schema_version >= APP_CONFIGURATION_SCHEMA_VERSION {
            return false;
        }
        if !saved
            .graphics_defaults
            .linear_sequence_helical_letters_enabled
        {
            saved
                .graphics_defaults
                .linear_sequence_helical_letters_enabled = true;
        }
        if saved.graphics_defaults.linear_sequence_letter_layout_mode
            != LinearSequenceLetterLayoutMode::AutoAdaptive
        {
            saved.graphics_defaults.linear_sequence_letter_layout_mode =
                LinearSequenceLetterLayoutMode::AutoAdaptive;
        }
        saved.schema_version = APP_CONFIGURATION_SCHEMA_VERSION;
        true
    }

    fn resolve_runtime_doc_path(path: &str) -> Option<PathBuf> {
        let direct = PathBuf::from(path);
        if direct.exists() {
            return Some(direct);
        }
        if let Ok(exe_path) = env::current_exe()
            && let Some(exe_dir) = exe_path.parent()
        {
            let bundled = exe_dir.join("../Resources").join(path);
            if bundled.exists() {
                return Some(bundled);
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
            let image_path = Self::help_image_render_path(&abs_path);
            let Some(rewritten_span) =
                Self::rewrite_inline_image_destination(span, image_path.as_path())
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

    fn help_display_markdown(markdown: &str) -> String {
        let markdown = Self::summarize_markdown_front_matter_for_help(markdown);
        Self::rewrite_markdown_inline_code_soft_breaks(&markdown)
    }

    fn help_markdown_hash(markdown: &str) -> u64 {
        let mut hasher = DefaultHasher::new();
        markdown.hash(&mut hasher);
        hasher.finish()
    }

    fn help_render_cache_key(&self, doc: HelpDoc) -> String {
        match doc {
            HelpDoc::Gui => "gui".to_string(),
            HelpDoc::Cli => "cli".to_string(),
            HelpDoc::AgentInterface => "agent_interface".to_string(),
            HelpDoc::ReviewerPreview => "reviewer_preview".to_string(),
            HelpDoc::Shell => format!("shell:{}", self.help_shell_interface.label()),
            HelpDoc::Tutorial => format!(
                "tutorial:{}:{}",
                self.help_tutorial_selected, self.help_tutorial_title
            ),
        }
    }

    fn help_source_markdown(&self, doc: HelpDoc) -> &str {
        match doc {
            HelpDoc::Gui => self.help_gui_markdown.as_str(),
            HelpDoc::Cli => self.help_cli_markdown.as_str(),
            HelpDoc::AgentInterface => self.help_agent_interface_markdown.as_str(),
            HelpDoc::ReviewerPreview => self.help_reviewer_preview_markdown.as_str(),
            HelpDoc::Shell => self.help_shell_markdown.as_str(),
            HelpDoc::Tutorial => self.help_tutorial_markdown.as_str(),
        }
    }

    fn rendered_help_markdown_for(&mut self, doc: HelpDoc) -> Arc<str> {
        let cache_key = self.help_render_cache_key(doc);
        let source_hash = Self::help_markdown_hash(self.help_source_markdown(doc));
        let needs_refresh = self
            .help_rendered_markdown_cache
            .get(&cache_key)
            .map(|entry| entry.source_hash != source_hash)
            .unwrap_or(true);
        if needs_refresh {
            let rendered =
                Arc::<str>::from(Self::help_display_markdown(self.help_source_markdown(doc)));
            self.help_rendered_markdown_cache.insert(
                cache_key.clone(),
                RenderedHelpMarkdownEntry {
                    source_hash,
                    rendered,
                },
            );
        }
        self.help_rendered_markdown_cache
            .get(&cache_key)
            .map(|entry| entry.rendered.clone())
            .unwrap_or_else(|| Arc::<str>::from(self.help_source_markdown(doc)))
    }

    fn generate_shell_help_markdown_for(interface: ShellHelpInterface) -> String {
        let interface_filter = interface.glossary_filter();
        if let Some(runtime_path) = Self::resolve_runtime_doc_path("docs/glossary.json")
            && let Ok(raw) = fs::read_to_string(runtime_path)
            && let Ok(markdown) =
                render_shell_help_markdown_from_glossary_json(&raw, interface_filter)
        {
            return markdown;
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
        self.help_agent_interface_markdown =
            Self::load_help_doc("docs/agent_interface.md", AGENT_INTERFACE_MD);
        self.help_reviewer_preview_markdown =
            Self::load_help_doc("docs/reviewer_preview.md", REVIEWER_PREVIEW_MD);
        self.help_shell_markdown =
            Self::generate_shell_help_markdown_for(self.help_shell_interface);
        self.help_rendered_markdown_cache.clear();
        self.help_tutorial_entries = Self::discover_help_tutorial_entries();
        self.set_help_tutorial_selected(self.help_tutorial_selected);
    }

    fn ensure_help_docs_loaded(&mut self) {
        if self.help_gui_markdown.is_empty()
            || self.help_cli_markdown.is_empty()
            || self.help_agent_interface_markdown.is_empty()
            || self.help_reviewer_preview_markdown.is_empty()
            || self.help_shell_markdown.is_empty()
        {
            self.refresh_help_docs();
        }
    }

    fn consume_native_help_request(&mut self) {
        if NATIVE_HELP_OPEN_REQUESTED.swap(false, Ordering::SeqCst) {
            let requested_at_ms = NATIVE_HELP_OPEN_REQUESTED_AT_MS.swap(0, Ordering::SeqCst);
            if requested_at_ms > 0 {
                let now_ms = now_unix_ms_u64();
                let dispatch_elapsed_ms = now_ms.saturating_sub(requested_at_ms) as u128;
                self.note_slow_phase("Native Help menu dispatch", dispatch_elapsed_ms);
            }
            self.open_help_doc(HelpDoc::Gui);
        }
    }

    fn consume_native_settings_request(&mut self) {
        if NATIVE_SETTINGS_OPEN_REQUESTED.swap(false, Ordering::SeqCst) {
            let graphics_tab_requested =
                NATIVE_SETTINGS_OPEN_GRAPHICS_TAB_REQUESTED.swap(false, Ordering::SeqCst);
            if graphics_tab_requested {
                self.open_configuration_graphics_dialog();
            } else {
                self.open_configuration_dialog();
            }
        }
    }

    fn consume_native_pcr_design_request(&mut self) {
        if NATIVE_PCR_DESIGN_OPEN_REQUESTED.swap(false, Ordering::SeqCst) {
            let requested_seq_id = NATIVE_PCR_DESIGN_SEQ_ID_REQUESTED
                .lock()
                .ok()
                .and_then(|mut guard| guard.take());
            if let Some(seq_id) = requested_seq_id.filter(|seq_id| !seq_id.trim().is_empty()) {
                if let Err(err) = self.open_pcr_design_dialog_for_seq_id(&seq_id) {
                    self.app_status = format!("Cannot open PCR Designer for '{seq_id}': {err}");
                }
            } else {
                self.open_pcr_design_dialog();
            }
        }
    }

    fn consume_native_jaspar_expert_request(&mut self) {
        if NATIVE_JASPAR_EXPERT_OPEN_REQUESTED.swap(false, Ordering::SeqCst) {
            let requested_motif_id = NATIVE_JASPAR_EXPERT_MOTIF_ID_REQUESTED
                .lock()
                .ok()
                .and_then(|mut guard| guard.take());
            if let Some(motif_id) =
                requested_motif_id.filter(|motif_id| !motif_id.trim().is_empty())
            {
                self.open_jaspar_expert_dialog_for_motif(&motif_id);
            } else {
                self.open_jaspar_expert_dialog();
            }
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
        self.configuration_ui_language = self.ui_language;
        self.configuration_language_dirty = false;
        self.window_backdrop_path_status_cache.clear();
        self.clear_rnapkin_validation();
        self.clear_blast_validation();
    }

    fn external_apps_configuration_dirty(&self) -> bool {
        let configured_rnapkin = tool_overrides::configured_or_env("GENTLE_RNAPKIN_BIN");
        let configured_makeblastdb = tool_overrides::configured_or_env(MAKEBLASTDB_ENV_BIN);
        let configured_blastn = tool_overrides::configured_or_env(BLASTN_ENV_BIN);
        let configured_bigwig = tool_overrides::configured_or_env(BIGWIG_TO_BEDGRAPH_ENV_BIN);
        self.configuration_rnapkin_executable.trim() != configured_rnapkin.trim()
            || self.configuration_makeblastdb_executable.trim() != configured_makeblastdb.trim()
            || self.configuration_blastn_executable.trim() != configured_blastn.trim()
            || self.configuration_bigwig_to_bedgraph_executable.trim() != configured_bigwig.trim()
    }

    fn configuration_has_unapplied_changes(&self) -> bool {
        self.external_apps_configuration_dirty()
            || self.configuration_graphics_dirty
            || self.configuration_window_backdrops_dirty
            || self.configuration_language_dirty
    }

    fn apply_pending_configuration_changes(&mut self) {
        let external_dirty = self.external_apps_configuration_dirty();
        let graphics_dirty = self.configuration_graphics_dirty;
        let backdrop_dirty = self.configuration_window_backdrops_dirty;
        let language_dirty = self.configuration_language_dirty;
        if !external_dirty && !graphics_dirty && !backdrop_dirty && !language_dirty {
            self.configuration_status = "No unapplied configuration changes".to_string();
            return;
        }
        if external_dirty {
            self.apply_configuration_external_apps();
        }
        if graphics_dirty {
            self.apply_configuration_graphics();
        }
        if backdrop_dirty {
            self.apply_configuration_window_backdrops();
        }
        if language_dirty {
            self.apply_configuration_language();
        }
    }

    fn open_configuration_dialog_for_tab(&mut self, tab: ConfigurationTab) {
        if self.show_configuration_dialog {
            self.configuration_tab = tab;
            self.mark_window_open_or_focus(Self::configuration_viewport_id(), true);
            return;
        }
        let sync_started = Instant::now();
        self.sync_configuration_from_runtime();
        self.note_slow_phase(
            "Configuration runtime sync",
            sync_started.elapsed().as_millis(),
        );
        self.configuration_tab = tab;
        self.show_configuration_dialog = true;
        self.configuration_status.clear();
        self.mark_window_open_or_focus(Self::configuration_viewport_id(), false);
    }

    fn open_configuration_dialog(&mut self) {
        self.open_configuration_dialog_for_tab(ConfigurationTab::ExternalApplications);
    }

    fn open_configuration_graphics_dialog(&mut self) {
        self.open_configuration_dialog_for_tab(ConfigurationTab::Graphics);
    }

    pub(crate) fn consume_command_or_ctrl_shortcut(ctx: &egui::Context, key: Key) -> bool {
        let command_shortcut = KeyboardShortcut::new(Modifiers::COMMAND, key);
        let ctrl_shortcut = KeyboardShortcut::new(Modifiers::CTRL, key);
        ctx.input_mut(|i| i.consume_shortcut(&command_shortcut))
            || ctx.input_mut(|i| i.consume_shortcut(&ctrl_shortcut))
    }

    pub(crate) fn viewport_close_requested_or_shortcut(ctx: &egui::Context) -> bool {
        let shortcut_triggered = Self::consume_command_or_ctrl_shortcut(ctx, Key::W);
        shortcut_triggered || ctx.input(|i| i.viewport().close_requested())
    }

    fn quit_application(&mut self) {
        self.pending_app_quit = true;
    }

    fn open_command_palette_dialog(&mut self) {
        let was_open = self.show_command_palette_dialog;
        self.show_command_palette_dialog = true;
        self.command_palette_focus_query = true;
        self.mark_window_open_or_focus(Self::command_palette_viewport_id(), was_open);
    }

    fn open_background_jobs_panel(&mut self) {
        let was_open = self.show_jobs_panel;
        self.show_jobs_panel = true;
        self.mark_window_open_or_focus(Self::background_jobs_viewport_id(), was_open);
    }

    fn toggle_background_jobs_panel(&mut self) {
        if self.show_jobs_panel {
            self.show_jobs_panel = false;
        } else {
            self.open_background_jobs_panel();
        }
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
        self.persist_background_job_history_to_state();
        next
    }

    fn now_unix_ms() -> u128 {
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map(|v| v.as_millis())
            .unwrap_or_default()
    }

    fn sha1_hex(value: &str) -> String {
        let mut hasher = Sha1::new();
        hasher.update(value.as_bytes());
        format!("{:x}", hasher.finalize())
    }

    fn current_prepare_retry_arguments(&self) -> serde_json::Value {
        let timeout_parse = self.parse_prepare_timeout_seconds();
        let (timeout_seconds, timeout_parse_error) = match timeout_parse {
            Ok(value) => (value, Option::<String>::None),
            Err(err) => (None, Some(err)),
        };
        serde_json::json!({
            "schema": "gentle.gui_retry_snapshot.prepare.v1",
            "genome_id": self.genome_id.trim(),
            "catalog_path": self.genome_catalog_path_opt(),
            "cache_dir": self.genome_cache_dir_opt(),
            "timeout_seconds_raw": self.genome_prepare_timeout_secs.trim(),
            "timeout_seconds": timeout_seconds,
            "timeout_parse_error": timeout_parse_error
        })
    }

    fn current_blast_retry_arguments(&self) -> serde_json::Value {
        let query_context = match self.genome_blast_source_mode {
            GenomeBlastSourceMode::Manual => {
                let query = self.genome_blast_query_manual.trim();
                serde_json::json!({
                    "manual_query_length_bp": query.len(),
                    "manual_query_sha1": if query.is_empty() {
                        None
                    } else {
                        Some(Self::sha1_hex(query))
                    }
                })
            }
            GenomeBlastSourceMode::ProjectSequence => serde_json::json!({
                "project_sequence_id": self.genome_blast_query_seq_id.trim()
            }),
            GenomeBlastSourceMode::ProjectPool => serde_json::json!({
                "project_pool_id": self.genome_blast_query_pool_id.trim()
            }),
        };
        let request_override = self.build_genome_blast_request_override_json();
        let (request_override_json, request_override_error) = match request_override {
            Ok(value) => (value, Option::<String>::None),
            Err(err) => (None, Some(err)),
        };
        serde_json::json!({
            "schema": "gentle.gui_retry_snapshot.blast.v1",
            "genome_id": self.genome_id.trim(),
            "source_mode": self.genome_blast_source_mode.as_str(),
            "query_context": query_context,
            "catalog_path": self.genome_catalog_path_opt(),
            "cache_dir": self.genome_cache_dir_opt(),
            "max_hits_legacy": self.genome_blast_max_hits.max(1),
            "task_name_raw": self.genome_blast_task_name.trim(),
            "options_preset": self.genome_blast_options_preset.label(),
            "thresholds": {
                "use_max_evalue": self.genome_blast_threshold_use_max_evalue,
                "max_evalue": self.genome_blast_threshold_max_evalue.trim(),
                "use_min_identity_percent": self.genome_blast_threshold_use_min_identity_percent,
                "min_identity_percent": self.genome_blast_threshold_min_identity_percent.trim(),
                "use_min_query_coverage_percent": self.genome_blast_threshold_use_min_query_coverage_percent,
                "min_query_coverage_percent": self.genome_blast_threshold_min_query_coverage_percent.trim(),
                "use_min_alignment_length_bp": self.genome_blast_threshold_use_min_alignment_length_bp,
                "min_alignment_length_bp": self.genome_blast_threshold_min_alignment_length_bp.trim(),
                "use_min_bit_score": self.genome_blast_threshold_use_min_bit_score,
                "min_bit_score": self.genome_blast_threshold_min_bit_score.trim(),
                "unique_best_hit": self.genome_blast_threshold_unique_best_hit
            },
            "options_json_raw": self.genome_blast_options_json.trim(),
            "request_override_json": request_override_json,
            "request_override_error": request_override_error
        })
    }

    fn current_track_import_retry_arguments(&self) -> serde_json::Value {
        let parsed = self.parse_bed_track_form();
        let (subscription, parse_error) = match parsed {
            Ok(value) => (Some(value), Option::<String>::None),
            Err(err) => (None, Some(err.to_string())),
        };
        let subscription_value = subscription.map(|item| {
            serde_json::json!({
                "source": item.source.label(),
                "path": item.path,
                "track_name": item.track_name,
                "min_score": item.min_score,
                "max_score": item.max_score,
                "clear_existing": item.clear_existing
            })
        });
        serde_json::json!({
            "schema": "gentle.gui_retry_snapshot.track_import.v1",
            "selected_sequence_id": self.genome_track_seq_id.trim(),
            "source_selection": self.genome_track_source_selection.as_str(),
            "path_raw": self.genome_track_path.trim(),
            "track_name_raw": self.genome_track_name.trim(),
            "min_score_raw": self.genome_track_min_score.trim(),
            "max_score_raw": self.genome_track_max_score.trim(),
            "clear_existing": self.genome_track_clear_existing,
            "parsed_subscription": subscription_value,
            "parse_error": parse_error
        })
    }

    fn current_agent_retry_arguments(&self) -> serde_json::Value {
        let prompt = self.agent_prompt.trim();
        let prompt_preview = if prompt.chars().count() > 200 {
            format!("{}...", prompt.chars().take(200).collect::<String>())
        } else {
            prompt.to_string()
        };
        let timeout_parse = self.parse_agent_timeout_seconds();
        let (timeout_secs, timeout_error) = match timeout_parse {
            Ok(value) => (value, Option::<String>::None),
            Err(err) => (None, Some(err)),
        };
        let connect_parse = self.parse_agent_connect_timeout_seconds();
        let (connect_timeout_secs, connect_timeout_error) = match connect_parse {
            Ok(value) => (value, Option::<String>::None),
            Err(err) => (None, Some(err)),
        };
        let read_parse = self.parse_agent_read_timeout_seconds();
        let (read_timeout_secs, read_timeout_error) = match read_parse {
            Ok(value) => (value, Option::<String>::None),
            Err(err) => (None, Some(err)),
        };
        let retries_parse = self.parse_agent_max_retries();
        let (max_retries, max_retries_error) = match retries_parse {
            Ok(value) => (value, Option::<String>::None),
            Err(err) => (None, Some(err)),
        };
        let max_response_parse = self.parse_agent_max_response_bytes();
        let (max_response_bytes, max_response_bytes_error) = match max_response_parse {
            Ok(value) => (value, Option::<String>::None),
            Err(err) => (None, Some(err)),
        };
        serde_json::json!({
            "schema": "gentle.gui_retry_snapshot.agent_assist.v1",
            "catalog_path": self.agent_catalog_path.trim(),
            "system_id": self.agent_system_id.trim(),
            "include_state_summary": self.agent_include_state_summary,
            "allow_auto_exec": self.agent_allow_auto_exec,
            "prompt_length_chars": prompt.len(),
            "prompt_sha1": if prompt.is_empty() {
                None
            } else {
                Some(Self::sha1_hex(prompt))
            },
            "prompt_preview": prompt_preview,
            "base_url_override": self.agent_base_url_override.trim(),
            "model_override": self.agent_model_override.trim(),
            "discovered_model_pick": self.agent_discovered_model_pick.trim(),
            "openai_api_key_override_set": !self.agent_openai_api_key.trim().is_empty(),
            "timeout_seconds_raw": self.agent_timeout_secs.trim(),
            "timeout_seconds": timeout_secs,
            "timeout_error": timeout_error,
            "connect_timeout_seconds_raw": self.agent_connect_timeout_secs.trim(),
            "connect_timeout_seconds": connect_timeout_secs,
            "connect_timeout_error": connect_timeout_error,
            "read_timeout_seconds_raw": self.agent_read_timeout_secs.trim(),
            "read_timeout_seconds": read_timeout_secs,
            "read_timeout_error": read_timeout_error,
            "max_retries_raw": self.agent_max_retries.trim(),
            "max_retries": max_retries,
            "max_retries_error": max_retries_error,
            "max_response_bytes_raw": self.agent_max_response_bytes.trim(),
            "max_response_bytes": max_response_bytes,
            "max_response_bytes_error": max_response_bytes_error
        })
    }

    fn normalize_background_job_retry_snapshots(
        snapshots: Vec<BackgroundJobRetrySnapshot>,
    ) -> Vec<BackgroundJobRetrySnapshot> {
        let mut normalized = snapshots
            .into_iter()
            .filter_map(|snapshot| {
                if snapshot.snapshot_id == 0 {
                    return None;
                }
                let origin = snapshot.origin.trim().to_string();
                Some(BackgroundJobRetrySnapshot { origin, ..snapshot })
            })
            .collect::<Vec<_>>();
        if normalized.len() > MAX_BACKGROUND_JOB_RETRY_SNAPSHOTS {
            let drain_len = normalized.len() - MAX_BACKGROUND_JOB_RETRY_SNAPSHOTS;
            normalized.drain(0..drain_len);
        }
        normalized
    }

    fn normalize_retry_snapshot_cleanup_audit_entries(
        entries: Vec<RetrySnapshotCleanupAuditEntry>,
    ) -> Vec<RetrySnapshotCleanupAuditEntry> {
        let mut normalized = entries
            .into_iter()
            .filter_map(|entry| {
                if entry.audit_id == 0 {
                    return None;
                }
                let action = entry.action.trim().to_string();
                if action.is_empty() {
                    return None;
                }
                let filter_kind = entry.filter_kind.trim().to_string();
                if filter_kind.is_empty() {
                    return None;
                }
                let summary = entry.summary.trim().to_string();
                let filter_text = entry.filter_text.trim().to_string();
                let archive_path = entry
                    .archive_path
                    .as_deref()
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .map(ToString::to_string);
                Some(RetrySnapshotCleanupAuditEntry {
                    action,
                    filter_kind,
                    summary,
                    filter_text,
                    archive_path,
                    ..entry
                })
            })
            .collect::<Vec<_>>();
        if normalized.len() > MAX_RETRY_SNAPSHOT_CLEANUP_AUDIT_ENTRIES {
            let drain_len = normalized.len() - MAX_RETRY_SNAPSHOT_CLEANUP_AUDIT_ENTRIES;
            normalized.drain(0..drain_len);
        }
        normalized
    }

    fn capture_retry_snapshot(
        &mut self,
        kind: BackgroundJobKind,
        origin: &str,
        arguments: serde_json::Value,
    ) -> u64 {
        let snapshot_id = self.next_retry_snapshot_id.max(1);
        self.next_retry_snapshot_id = snapshot_id.saturating_add(1);
        self.retry_argument_snapshots
            .push(BackgroundJobRetrySnapshot {
                snapshot_id,
                kind,
                origin: origin.trim().to_string(),
                captured_at_unix_ms: Self::now_unix_ms(),
                arguments,
            });
        let retention_limit = self.retry_snapshot_retention_limit();
        if self.retry_argument_snapshots.len() > retention_limit {
            let drain_len = self.retry_argument_snapshots.len() - retention_limit;
            self.retry_argument_snapshots.drain(0..drain_len);
        }
        self.persist_background_job_history_to_state();
        snapshot_id
    }

    fn record_retry_snapshot_cleanup_audit(
        &mut self,
        action: &str,
        summary: &str,
        target_count: usize,
        removed_count: usize,
        retained_before: usize,
        retained_after: usize,
        archive_path: Option<&str>,
    ) -> u64 {
        let action = action.trim();
        let action = if action.is_empty() {
            "unknown_cleanup_action"
        } else {
            action
        };
        let audit_id = self.next_retry_cleanup_audit_id.max(1);
        self.next_retry_cleanup_audit_id = audit_id.saturating_add(1);
        self.retry_snapshot_cleanup_audit
            .push(RetrySnapshotCleanupAuditEntry {
                audit_id,
                action: action.to_string(),
                performed_at_unix_ms: Self::now_unix_ms(),
                target_count,
                removed_count,
                retained_before,
                retained_after,
                filter_kind: self.retry_snapshot_kind_filter.export_token().to_string(),
                filter_text: self.retry_snapshot_text_filter.trim().to_string(),
                archive_path: archive_path
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .map(ToString::to_string),
                summary: summary.trim().to_string(),
            });
        let retention_limit = self.retry_cleanup_audit_retention_limit();
        if self.retry_snapshot_cleanup_audit.len() > retention_limit {
            let drain_len = self.retry_snapshot_cleanup_audit.len() - retention_limit;
            self.retry_snapshot_cleanup_audit.drain(0..drain_len);
        }
        self.persist_background_job_history_to_state();
        audit_id
    }

    fn retry_snapshot_retention_limit(&self) -> usize {
        self.retry_snapshot_retain_count
            .clamp(1, MAX_BACKGROUND_JOB_RETRY_SNAPSHOTS)
    }

    fn retry_cleanup_audit_retention_limit(&self) -> usize {
        self.retry_cleanup_audit_retain_count
            .clamp(1, MAX_RETRY_SNAPSHOT_CLEANUP_AUDIT_ENTRIES)
    }

    fn prune_retry_snapshots_to_limit(&mut self, retain_count: usize) -> usize {
        let retain_count = retain_count.clamp(1, MAX_BACKGROUND_JOB_RETRY_SNAPSHOTS);
        self.retry_snapshot_retain_count = retain_count;
        let before = self.retry_argument_snapshots.len();
        if before > retain_count {
            let drain_len = before - retain_count;
            self.retry_argument_snapshots.drain(0..drain_len);
            self.persist_background_job_history_to_state();
        }
        before.saturating_sub(self.retry_argument_snapshots.len())
    }

    fn prune_retry_cleanup_audit_to_limit(&mut self, retain_count: usize) -> usize {
        let retain_count = retain_count.clamp(1, MAX_RETRY_SNAPSHOT_CLEANUP_AUDIT_ENTRIES);
        self.retry_cleanup_audit_retain_count = retain_count;
        let before = self.retry_snapshot_cleanup_audit.len();
        if before > retain_count {
            let drain_len = before - retain_count;
            self.retry_snapshot_cleanup_audit.drain(0..drain_len);
            self.persist_background_job_history_to_state();
        }
        before.saturating_sub(self.retry_snapshot_cleanup_audit.len())
    }

    fn clear_retry_snapshots(&mut self) -> usize {
        let removed = self.retry_argument_snapshots.len();
        if removed > 0 {
            self.retry_argument_snapshots.clear();
            self.persist_background_job_history_to_state();
        }
        removed
    }

    fn clear_retry_cleanup_audit(&mut self) -> usize {
        let removed = self.retry_snapshot_cleanup_audit.len();
        if removed > 0 {
            self.retry_snapshot_cleanup_audit.clear();
            self.persist_background_job_history_to_state();
        }
        removed
    }

    fn retry_snapshot_filter_help_text() -> &'static str {
        "Free text matches kind/origin/arguments. Scoped terms: kind:blast origin:panel args:hg38"
    }

    fn retry_cleanup_audit_filter_help_text() -> &'static str {
        "Free text matches action/filter/summary/archive. Scoped terms: action:delete kind:blast summary:pruned archive:json"
    }

    fn retry_snapshot_matches_filter(
        snapshot: &BackgroundJobRetrySnapshot,
        kind_filter: RetrySnapshotKindFilter,
        filter_text: &str,
    ) -> bool {
        if !kind_filter.matches_kind(snapshot.kind) {
            return false;
        }
        let terms = Self::parse_track_filter_terms(filter_text);
        if terms.is_empty() {
            return true;
        }
        let kind_terms = vec![
            snapshot.kind.label().to_ascii_lowercase(),
            serde_json::to_string(&snapshot.kind)
                .unwrap_or_else(|_| String::new())
                .to_ascii_lowercase(),
        ];
        let origin_terms = vec![snapshot.origin.to_ascii_lowercase()];
        let args_terms = vec![snapshot.arguments.to_string().to_ascii_lowercase()];
        let all_terms = kind_terms
            .iter()
            .chain(origin_terms.iter())
            .chain(args_terms.iter())
            .cloned()
            .collect::<Vec<_>>();
        terms.iter().all(|(scope, needle)| {
            let search_space: &Vec<String> = match scope.as_deref() {
                Some("kind") | Some("job") | Some("type") => &kind_terms,
                Some("origin") | Some("from") => &origin_terms,
                Some("args") | Some("json") => &args_terms,
                _ => &all_terms,
            };
            search_space.iter().any(|value| value.contains(needle))
        })
    }

    fn filtered_retry_snapshots(&self) -> Vec<BackgroundJobRetrySnapshot> {
        self.retry_argument_snapshots
            .iter()
            .filter(|snapshot| {
                Self::retry_snapshot_matches_filter(
                    snapshot,
                    self.retry_snapshot_kind_filter,
                    &self.retry_snapshot_text_filter,
                )
            })
            .cloned()
            .collect()
    }

    fn retry_cleanup_audit_entry_matches_filter(
        entry: &RetrySnapshotCleanupAuditEntry,
        action_filter: RetryCleanupAuditActionFilter,
        filter_text: &str,
    ) -> bool {
        if !action_filter.matches_action(entry.action.as_str()) {
            return false;
        }
        let terms = Self::parse_track_filter_terms(filter_text);
        if terms.is_empty() {
            return true;
        }
        let action_terms = vec![entry.action.to_ascii_lowercase()];
        let filter_kind_terms = vec![entry.filter_kind.to_ascii_lowercase()];
        let filter_text_terms = if entry.filter_text.trim().is_empty() {
            Vec::<String>::new()
        } else {
            vec![entry.filter_text.to_ascii_lowercase()]
        };
        let summary_terms = if entry.summary.trim().is_empty() {
            Vec::<String>::new()
        } else {
            vec![entry.summary.to_ascii_lowercase()]
        };
        let archive_terms = {
            let mut values = Vec::<String>::new();
            if let Some(path) = entry.archive_path.as_deref() {
                values.push(path.to_ascii_lowercase());
                if let Some(file_name) =
                    Path::new(path).file_name().and_then(|value| value.to_str())
                {
                    values.push(file_name.to_ascii_lowercase());
                }
            }
            values
        };
        let target_terms = vec![entry.target_count.to_string()];
        let removed_terms = vec![entry.removed_count.to_string()];
        let before_terms = vec![entry.retained_before.to_string()];
        let after_terms = vec![entry.retained_after.to_string()];
        let count_terms = vec![
            entry.target_count.to_string(),
            entry.removed_count.to_string(),
            entry.retained_before.to_string(),
            entry.retained_after.to_string(),
        ];
        let all_terms = action_terms
            .iter()
            .chain(filter_kind_terms.iter())
            .chain(filter_text_terms.iter())
            .chain(summary_terms.iter())
            .chain(archive_terms.iter())
            .chain(count_terms.iter())
            .cloned()
            .collect::<Vec<_>>();
        terms.iter().all(|(scope, needle)| {
            let search_space: &Vec<String> = match scope.as_deref() {
                Some("action") | Some("audit") => &action_terms,
                Some("kind") | Some("filter_kind") => &filter_kind_terms,
                Some("filter") | Some("text") => &filter_text_terms,
                Some("summary") => &summary_terms,
                Some("archive") | Some("path") | Some("file") => &archive_terms,
                Some("target") => &target_terms,
                Some("removed") => &removed_terms,
                Some("before") | Some("retained_before") => &before_terms,
                Some("after") | Some("retained_after") => &after_terms,
                Some("count") | Some("counts") => &count_terms,
                _ => &all_terms,
            };
            search_space.iter().any(|value| value.contains(needle))
        })
    }

    fn filtered_retry_cleanup_audit_entries(&self) -> Vec<RetrySnapshotCleanupAuditEntry> {
        self.retry_snapshot_cleanup_audit
            .iter()
            .filter(|entry| {
                Self::retry_cleanup_audit_entry_matches_filter(
                    entry,
                    self.retry_cleanup_audit_action_filter,
                    &self.retry_cleanup_audit_text_filter,
                )
            })
            .cloned()
            .collect()
    }

    fn retry_snapshot_dry_run_diff(&self) -> RetrySnapshotDryRunDiff {
        let mut diff = RetrySnapshotDryRunDiff::default();
        let kind_filter = self.retry_snapshot_kind_filter;
        let filter_text = self.retry_snapshot_text_filter.clone();
        for snapshot in self.retry_argument_snapshots.iter().cloned() {
            if Self::retry_snapshot_matches_filter(&snapshot, kind_filter, &filter_text) {
                diff.removed.push(snapshot);
            } else {
                diff.retained.push(snapshot);
            }
        }
        diff
    }

    fn build_retry_snapshot_export_payload(
        &self,
        snapshots: Vec<BackgroundJobRetrySnapshot>,
    ) -> serde_json::Value {
        let snapshot_count = snapshots.len();
        serde_json::json!({
            "schema": "gentle.gui_retry_snapshot_export.v1",
            "exported_at_unix_ms": Self::now_unix_ms(),
            "filters": {
                "kind": self.retry_snapshot_kind_filter.export_token(),
                "text": self.retry_snapshot_text_filter.trim()
            },
            "snapshot_count": snapshot_count,
            "snapshots": snapshots
        })
    }

    fn export_filtered_retry_snapshots_to_path(
        &self,
        path: &Path,
    ) -> std::result::Result<usize, String> {
        let snapshots = self.filtered_retry_snapshots();
        if snapshots.is_empty() {
            return Err("No retry snapshots match current filters".to_string());
        }
        let snapshot_count = snapshots.len();
        let payload = self.build_retry_snapshot_export_payload(snapshots);
        let json = serde_json::to_string_pretty(&payload)
            .map_err(|e| format!("Could not serialize retry snapshot export JSON: {e}"))?;
        fs::write(path, json).map_err(|e| {
            format!(
                "Could not write retry snapshot export '{}': {e}",
                path.display()
            )
        })?;
        Ok(snapshot_count)
    }

    fn remove_filtered_retry_snapshots(&mut self) -> usize {
        let before = self.retry_argument_snapshots.len();
        if before == 0 {
            return 0;
        }
        let kind_filter = self.retry_snapshot_kind_filter;
        let filter_text = self.retry_snapshot_text_filter.clone();
        self.retry_argument_snapshots.retain(|snapshot| {
            !Self::retry_snapshot_matches_filter(snapshot, kind_filter, &filter_text)
        });
        let removed = before.saturating_sub(self.retry_argument_snapshots.len());
        if removed > 0 {
            self.persist_background_job_history_to_state();
        }
        removed
    }

    fn summarize_retry_snapshot_cleanup_targets(
        snapshots: &[BackgroundJobRetrySnapshot],
    ) -> String {
        if snapshots.is_empty() {
            return "No matching retry snapshots".to_string();
        }
        let mut prepare = 0usize;
        let mut blast = 0usize;
        let mut track_import = 0usize;
        let mut tutorial = 0usize;
        let mut agent = 0usize;
        let mut origin_counts: BTreeMap<String, usize> = BTreeMap::new();
        let mut min_snapshot_id = u64::MAX;
        let mut max_snapshot_id = 0u64;
        let mut min_ts = u128::MAX;
        let mut max_ts = 0u128;
        for snapshot in snapshots {
            match snapshot.kind {
                BackgroundJobKind::PrepareGenome => prepare += 1,
                BackgroundJobKind::BlastGenome => blast += 1,
                BackgroundJobKind::TrackImport => track_import += 1,
                BackgroundJobKind::OpenTutorialProject => tutorial += 1,
                BackgroundJobKind::AgentAssist => agent += 1,
            }
            let origin = if snapshot.origin.trim().is_empty() {
                "-".to_string()
            } else {
                snapshot.origin.trim().to_string()
            };
            *origin_counts.entry(origin).or_insert(0usize) += 1;
            min_snapshot_id = min_snapshot_id.min(snapshot.snapshot_id);
            max_snapshot_id = max_snapshot_id.max(snapshot.snapshot_id);
            min_ts = min_ts.min(snapshot.captured_at_unix_ms);
            max_ts = max_ts.max(snapshot.captured_at_unix_ms);
        }
        let mut origin_entries = origin_counts.into_iter().collect::<Vec<_>>();
        origin_entries.sort_by(|(left_origin, left_count), (right_origin, right_count)| {
            right_count
                .cmp(left_count)
                .then_with(|| left_origin.cmp(right_origin))
        });
        let mut origin_preview = origin_entries
            .iter()
            .take(3)
            .map(|(origin, count)| format!("{origin}={count}"))
            .collect::<Vec<_>>()
            .join(", ");
        if origin_entries.len() > 3 {
            if !origin_preview.is_empty() {
                origin_preview.push_str(", ...");
            } else {
                origin_preview = "...".to_string();
            }
        }
        format!(
            "{} match(es) · kinds: PrepareGenome={} BlastGenome={} TrackImport={} OpenTutorialProject={} AgentAssist={} · origins: {} · ids #{}..#{} · captured_at {}..{}",
            snapshots.len(),
            prepare,
            blast,
            track_import,
            tutorial,
            agent,
            if origin_preview.is_empty() {
                "-".to_string()
            } else {
                origin_preview
            },
            min_snapshot_id,
            max_snapshot_id,
            min_ts,
            max_ts
        )
    }

    fn format_retry_snapshot_cleanup_status(
        action: RetrySnapshotPendingCleanupAction,
        removed: usize,
        before: usize,
        after: usize,
    ) -> String {
        format!(
            "{} confirmed: removed {} snapshot(s), retained {} (was {})",
            action.label(),
            removed,
            after,
            before
        )
    }

    fn retry_snapshot_cleanup_confirm_input_matches(
        action: RetrySnapshotPendingCleanupAction,
        target_count: usize,
        input: &str,
    ) -> bool {
        let expected = action.confirm_phrase(target_count);
        input.trim().eq_ignore_ascii_case(expected.as_str())
    }

    fn retry_cleanup_audit_clear_confirm_phrase(target_count: usize) -> String {
        format!("clear-audit {target_count}")
    }

    fn retry_cleanup_audit_clear_confirm_input_matches(target_count: usize, input: &str) -> bool {
        let expected = Self::retry_cleanup_audit_clear_confirm_phrase(target_count);
        input.trim().eq_ignore_ascii_case(expected.as_str())
    }

    fn build_retry_cleanup_audit_report_payload(
        &self,
        entries: Vec<RetrySnapshotCleanupAuditEntry>,
    ) -> serde_json::Value {
        let entry_count = entries.len();
        let mut action_counts: BTreeMap<String, usize> = BTreeMap::new();
        for entry in entries.iter() {
            let key = entry.action.trim().to_string();
            *action_counts.entry(key).or_insert(0usize) += 1;
        }
        serde_json::json!({
            "schema": "gentle.gui_retry_cleanup_audit_report.v1",
            "exported_at_unix_ms": Self::now_unix_ms(),
            "filters": {
                "action": self.retry_cleanup_audit_action_filter.label(),
                "text": self.retry_cleanup_audit_text_filter.trim()
            },
            "entry_count": entry_count,
            "retained_entry_count": self.retry_snapshot_cleanup_audit.len(),
            "action_counts": action_counts,
            "entries": entries
        })
    }

    fn export_retry_cleanup_audit_report_to_path(
        &self,
        path: &Path,
    ) -> std::result::Result<usize, String> {
        let entries = self.filtered_retry_cleanup_audit_entries();
        if entries.is_empty() {
            return Err("No retry cleanup audit entries match current filters".to_string());
        }
        let entry_count = entries.len();
        let payload = self.build_retry_cleanup_audit_report_payload(entries);
        let json = serde_json::to_string_pretty(&payload)
            .map_err(|e| format!("Could not serialize cleanup audit report JSON: {e}"))?;
        fs::write(path, json).map_err(|e| {
            format!(
                "Could not write cleanup audit report '{}': {e}",
                path.display()
            )
        })?;
        Ok(entry_count)
    }

    fn build_retry_snapshot_archive_payload(
        &self,
        snapshots: Vec<BackgroundJobRetrySnapshot>,
    ) -> serde_json::Value {
        let archived_count = snapshots.len();
        serde_json::json!({
            "schema": "gentle.gui_retry_snapshot_archive.v1",
            "archived_at_unix_ms": Self::now_unix_ms(),
            "archive_mode": "archive_and_delete_filtered",
            "filters": {
                "kind": self.retry_snapshot_kind_filter.export_token(),
                "text": self.retry_snapshot_text_filter.trim()
            },
            "archived_count": archived_count,
            "snapshots": snapshots
        })
    }

    fn archive_and_delete_filtered_retry_snapshots_to_path(
        &mut self,
        path: &Path,
    ) -> std::result::Result<usize, String> {
        let snapshots = self.filtered_retry_snapshots();
        if snapshots.is_empty() {
            return Err("No retry snapshots match current filters".to_string());
        }
        let payload = self.build_retry_snapshot_archive_payload(snapshots.clone());
        let json = serde_json::to_string_pretty(&payload)
            .map_err(|e| format!("Could not serialize retry snapshot archive JSON: {e}"))?;
        fs::write(path, json).map_err(|e| {
            format!(
                "Could not write retry snapshot archive '{}': {e}",
                path.display()
            )
        })?;
        let removed = self.remove_filtered_retry_snapshots();
        Ok(removed)
    }

    fn normalize_background_job_events(events: Vec<BackgroundJobEvent>) -> Vec<BackgroundJobEvent> {
        let mut normalized = events
            .into_iter()
            .filter_map(|event| {
                let summary = event.summary.trim();
                if summary.is_empty() {
                    None
                } else {
                    Some(BackgroundJobEvent {
                        summary: summary.to_string(),
                        ..event
                    })
                }
            })
            .collect::<Vec<_>>();
        if normalized.len() > MAX_BACKGROUND_JOB_EVENTS {
            let drain_len = normalized.len() - MAX_BACKGROUND_JOB_EVENTS;
            normalized.drain(0..drain_len);
        }
        normalized
    }

    fn load_background_job_history_from_state(&mut self) {
        let raw = self
            .engine
            .read()
            .unwrap()
            .state()
            .metadata
            .get(BACKGROUND_JOB_HISTORY_METADATA_KEY)
            .cloned();
        let Some(raw) = raw else {
            self.next_background_job_id = 1;
            self.job_event_log.clear();
            self.next_retry_snapshot_id = 1;
            self.retry_argument_snapshots.clear();
            self.next_retry_cleanup_audit_id = 1;
            self.retry_snapshot_cleanup_audit.clear();
            return;
        };
        let Ok(mut parsed) = serde_json::from_value::<PersistedBackgroundJobHistory>(raw) else {
            self.next_background_job_id = 1;
            self.job_event_log.clear();
            self.next_retry_snapshot_id = 1;
            self.retry_argument_snapshots.clear();
            self.next_retry_cleanup_audit_id = 1;
            self.retry_snapshot_cleanup_audit.clear();
            return;
        };
        if parsed.schema.trim().is_empty() {
            parsed.schema = BACKGROUND_JOB_HISTORY_SCHEMA.to_string();
        }
        if parsed.schema != BACKGROUND_JOB_HISTORY_SCHEMA {
            self.next_background_job_id = 1;
            self.job_event_log.clear();
            self.next_retry_snapshot_id = 1;
            self.retry_argument_snapshots.clear();
            self.next_retry_cleanup_audit_id = 1;
            self.retry_snapshot_cleanup_audit.clear();
            return;
        }
        self.next_background_job_id = parsed.next_background_job_id.max(1);
        self.job_event_log = Self::normalize_background_job_events(parsed.events);
        self.next_retry_snapshot_id = parsed.next_retry_snapshot_id.max(1);
        self.retry_argument_snapshots =
            Self::normalize_background_job_retry_snapshots(parsed.retry_snapshots);
        self.next_retry_cleanup_audit_id = parsed.next_retry_cleanup_audit_id.max(1);
        self.retry_snapshot_cleanup_audit =
            Self::normalize_retry_snapshot_cleanup_audit_entries(parsed.retry_cleanup_audit);
    }

    fn persist_background_job_history_to_state(&mut self) {
        let payload = PersistedBackgroundJobHistory {
            schema: BACKGROUND_JOB_HISTORY_SCHEMA.to_string(),
            next_background_job_id: self.next_background_job_id.max(1),
            events: Self::normalize_background_job_events(self.job_event_log.clone()),
            next_retry_snapshot_id: self.next_retry_snapshot_id.max(1),
            retry_snapshots: Self::normalize_background_job_retry_snapshots(
                self.retry_argument_snapshots.clone(),
            ),
            next_retry_cleanup_audit_id: self.next_retry_cleanup_audit_id.max(1),
            retry_cleanup_audit: Self::normalize_retry_snapshot_cleanup_audit_entries(
                self.retry_snapshot_cleanup_audit.clone(),
            ),
        };
        let Ok(value) = serde_json::to_value(payload) else {
            return;
        };
        let mut engine = self.engine.write().unwrap();
        let state = engine.state_mut();
        if state.metadata.get(BACKGROUND_JOB_HISTORY_METADATA_KEY) == Some(&value) {
            return;
        }
        state
            .metadata
            .insert(BACKGROUND_JOB_HISTORY_METADATA_KEY.to_string(), value);
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
        if self.job_event_log.len() > MAX_BACKGROUND_JOB_EVENTS {
            let drain_len = self.job_event_log.len() - MAX_BACKGROUND_JOB_EVENTS;
            self.job_event_log.drain(0..drain_len);
        }
        self.persist_background_job_history_to_state();
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

    fn request_tutorial_project_task_cancel(&mut self, origin: &str) {
        let Some((job_id, chapter_title, already_requested)) =
            self.tutorial_project_task.as_ref().map(|task| {
                (
                    task.job_id,
                    task.chapter_title.clone(),
                    task.cancel_requested.swap(true, Ordering::Relaxed),
                )
            })
        else {
            self.tutorial_project_status =
                "No running tutorial project build to cancel".to_string();
            self.app_status = self.tutorial_project_status.clone();
            return;
        };

        if already_requested {
            self.tutorial_project_status =
                format!("Cancellation was already requested for tutorial '{chapter_title}'");
            self.app_status = self.tutorial_project_status.clone();
            return;
        }

        self.tutorial_project_status =
            format!("Cancellation requested for tutorial '{chapter_title}'");
        self.app_status = self.tutorial_project_status.clone();
        self.push_job_event(
            BackgroundJobKind::OpenTutorialProject,
            BackgroundJobEventPhase::CancelRequested,
            Some(job_id),
            format!("Cancellation requested from {origin}"),
        );
    }

    fn request_blast_task_cancel(&mut self, origin: &str) {
        let Some((job_id, already_requested)) = self.genome_blast_task.as_ref().map(|task| {
            (
                task.job_id,
                task.cancel_requested.swap(true, Ordering::Relaxed),
            )
        }) else {
            self.genome_blast_status = "No running BLAST job to cancel".to_string();
            return;
        };

        if already_requested {
            self.genome_blast_status = "Cancellation was already requested for BLAST".to_string();
            return;
        }

        self.genome_blast_status = "Cancellation requested for running BLAST job".to_string();
        self.push_job_event(
            BackgroundJobKind::BlastGenome,
            BackgroundJobEventPhase::CancelRequested,
            Some(job_id),
            format!("Cancellation requested from {origin}"),
        );
    }

    fn has_active_background_jobs(&self) -> bool {
        self.genome_prepare_task.is_some()
            || self.genome_blast_task.is_some()
            || self.genome_track_import_task.is_some()
            || self.tutorial_project_task.is_some()
            || self.agent_task.is_some()
    }

    fn refresh_sequence_windows_for_seq_ids(&mut self, seq_ids: &[String]) -> usize {
        if seq_ids.is_empty() {
            return 0;
        }
        let target_ids: HashSet<&str> = seq_ids.iter().map(|id| id.as_str()).collect();
        let mut refreshed = 0usize;
        for window in self.windows.values() {
            if let Ok(mut guard) = window.write() {
                let should_refresh = guard
                    .sequence_id()
                    .as_deref()
                    .map(|seq_id| target_ids.contains(seq_id))
                    .unwrap_or(false);
                if should_refresh {
                    guard.refresh_from_engine_state();
                    refreshed += 1;
                }
            }
        }
        refreshed
    }

    fn refresh_sequence_windows_from_engine_state(&mut self) -> usize {
        let mut refreshed = 0usize;
        for window in self.windows.values() {
            if let Ok(mut guard) = window.write() {
                guard.refresh_from_engine_state();
                refreshed += 1;
            }
        }
        refreshed
    }

    fn handle_engine_state_after_history_transition(&mut self) {
        self.lineage_cache_valid = false;
        self.lineage_rows.clear();
        self.lineage_edges.clear();
        self.lineage_op_label_by_id.clear();
        self.lineage_containers.clear();
        self.lineage_arrangements.clear();
        self.lineage_racks.clear();
        self.load_bed_track_subscriptions_from_state();
        self.tracked_autosync_last_op_count = None;
        self.genome_track_autosync_status.clear();
        self.refresh_sequence_windows_from_engine_state();
    }

    fn undo_last_operation(&mut self) {
        if self.has_active_background_jobs() {
            self.app_status = "Undo is disabled while background jobs are active.".to_string();
            return;
        }
        let outcome = self
            .engine
            .write()
            .map_err(|_| EngineError {
                code: ErrorCode::Internal,
                message: "Engine lock poisoned".to_string(),

                cause_chain: vec![],
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
            self.app_status = "Redo is disabled while background jobs are active.".to_string();
            return;
        }
        let outcome = self
            .engine
            .write()
            .map_err(|_| EngineError {
                code: ErrorCode::Internal,
                message: "Engine lock poisoned".to_string(),

                cause_chain: vec![],
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
                title: "New Sequence".to_string(),
                detail: "Type or paste IUPAC DNA into a new project sequence".to_string(),
                keywords: "sequence new paste clipboard type create".to_string(),
                action: CommandPaletteAction::NewSequence,
            },
            CommandPaletteEntry {
                title: "New Sequence from Clipboard".to_string(),
                detail: "Read clipboard text into the new-sequence dialog".to_string(),
                keywords: "sequence new paste clipboard create".to_string(),
                action: CommandPaletteAction::NewSequenceFromClipboard,
            },
            CommandPaletteEntry {
                title: "Open Sequence".to_string(),
                detail: "Import one or more sequence files into project".to_string(),
                keywords: "sequence import load".to_string(),
                action: CommandPaletteAction::OpenSequence,
            },
            CommandPaletteEntry {
                title: "Protein Evidence".to_string(),
                detail: "Fetch/import UniProt or Ensembl protein evidence and compare it against transcript-native products".to_string(),
                keywords: "uniprot swiss-prot ensembl protein map compare evidence".to_string(),
                action: CommandPaletteAction::OpenUniprot,
            },
            CommandPaletteEntry {
                title: "GenBank / dbSNP Fetch".to_string(),
                detail: "Fetch GenBank accessions or rsID-anchored genome regions from NCBI"
                    .to_string(),
                keywords: "genbank dbsnp rsid snp ncbi accession fetch import genome".to_string(),
                action: CommandPaletteAction::OpenGenbank,
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
                title: "External Services".to_string(),
                detail:
                    "Inspect providers, preflight service requests, and prepare quote handoff bundles"
                        .to_string(),
                keywords: "services provider metabion geneart quote handoff wop email excel oligo m-block external vendor cro"
                    .to_string(),
                action: CommandPaletteAction::OpenExternalServices,
            },
            CommandPaletteEntry {
                title: "Gibson".to_string(),
                detail: "Plan one destination-first Gibson assembly with cartoon preview."
                    .to_string(),
                keywords: "gibson cloning overlaps primers cartoon".to_string(),
                action: CommandPaletteAction::OpenGibson,
            },
            CommandPaletteEntry {
                title: "Planning".to_string(),
                detail: "Edit planning profiles/objectives and resolve suggestions".to_string(),
                keywords: "planning profile objective suggestions sync meta-layer".to_string(),
                action: CommandPaletteAction::OpenPlanning,
            },
            CommandPaletteEntry {
                title: "Routine Assistant".to_string(),
                detail: "Compare routine families and run semantic preflight".to_string(),
                keywords: "routine assistant cloning preflight".to_string(),
                action: CommandPaletteAction::OpenRoutineAssistant,
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
                title: "Agent Interface".to_string(),
                detail: "Open in-app agent integration guide".to_string(),
                keywords: "help agent interface mcp cli assistant prompt docs".to_string(),
                action: CommandPaletteAction::OpenAgentInterfaceManual,
            },
            CommandPaletteEntry {
                title: "Reviewer Quickstart".to_string(),
                detail: "Open internal preview quickstart and known limitations".to_string(),
                keywords: "help reviewer quickstart preview internal limitations".to_string(),
                action: CommandPaletteAction::OpenReviewerPreviewManual,
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
                title: "Export Lab Assistant Report".to_string(),
                detail: "Export bench-facing cloning handoff as ODT, DOCX, or Markdown"
                    .to_string(),
                keywords: "export lab assistant report handoff odt docx markdown".to_string(),
                action: CommandPaletteAction::ExportLabAssistantReport,
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
        entries.extend(
            UiIntentTarget::all()
                .iter()
                .copied()
                .map(CommandPaletteEntry::ui_intent),
        );
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
            CommandPaletteAction::NewSequence => self.open_new_sequence_dialog(),
            CommandPaletteAction::NewSequenceFromClipboard => {
                self.open_new_sequence_from_clipboard()
            }
            CommandPaletteAction::OpenSequence => self.prompt_open_sequence(),
            CommandPaletteAction::OpenUniprot => self.open_uniprot_dialog(),
            CommandPaletteAction::OpenGenbank => self.open_genbank_dialog(),
            CommandPaletteAction::OpenConfiguration => self.open_configuration_dialog(),
            CommandPaletteAction::OpenExternalServices => self.open_external_services_dialog(),
            CommandPaletteAction::OpenGibson => self.open_gibson_dialog(),
            CommandPaletteAction::OpenPlanning => self.open_planning_dialog(),
            CommandPaletteAction::OpenRoutineAssistant => self.open_routine_assistant_dialog(),
            CommandPaletteAction::UiIntent(target) => {
                let command = ShellCommand::UiIntent {
                    action: UiIntentAction::Open,
                    target,
                    genome_id: None,
                    helper_mode: false,
                    catalog_path: None,
                    cache_dir: None,
                    filter: None,
                    species: None,
                    latest: false,
                };
                if let Some(summary) = self.try_apply_shell_ui_intent(&command) {
                    self.app_status = summary;
                }
            }
            CommandPaletteAction::OpenGuiManual => self.open_help_doc(HelpDoc::Gui),
            CommandPaletteAction::OpenCliManual => self.open_help_doc(HelpDoc::Cli),
            CommandPaletteAction::OpenAgentInterfaceManual => {
                self.open_help_doc(HelpDoc::AgentInterface)
            }
            CommandPaletteAction::OpenReviewerPreviewManual => {
                self.open_help_doc(HelpDoc::ReviewerPreview)
            }
            CommandPaletteAction::OpenShellManual => self.open_help_doc(HelpDoc::Shell),
            CommandPaletteAction::ExportLineageSvg => self.prompt_export_lineage_svg(),
            CommandPaletteAction::ExportLabAssistantReport => {
                self.prompt_export_lab_assistant_report()
            }
            CommandPaletteAction::ToggleJobsPanel => {
                self.toggle_background_jobs_panel();
            }
            CommandPaletteAction::ToggleHistoryPanel => {
                self.history_ui.show_panel = !self.history_ui.show_panel;
            }
            CommandPaletteAction::Undo => self.undo_last_operation(),
            CommandPaletteAction::Redo => self.redo_last_operation(),
            CommandPaletteAction::FocusViewport(viewport_id) => {
                self.focus_window_viewport(ctx, viewport_id);
            }
        }
    }

    fn open_help_doc(&mut self, doc: HelpDoc) {
        if doc == HelpDoc::Tutorial {
            self.open_help_tutorial_doc(self.help_tutorial_selected);
            return;
        }
        if self.show_help_dialog && self.help_doc == doc {
            self.queue_focus_viewport(Self::help_viewport_id());
            return;
        }
        let opening_from_closed_state = !self.show_help_dialog;
        if !self.show_help_dialog {
            self.mark_viewport_open_requested(Self::help_viewport_id());
        }
        let help_load_started = Instant::now();
        self.ensure_help_docs_loaded();
        self.note_slow_phase("Help payload load", help_load_started.elapsed().as_millis());
        let help_doc_changed = self.help_doc != doc;
        self.help_doc = doc;
        if opening_from_closed_state || help_doc_changed {
            self.refresh_help_search_matches();
        }
        self.help_focus_search_box = true;
        self.show_help_dialog = true;
        self.queue_focus_viewport(Self::help_viewport_id());
    }

    fn active_help_title(&self) -> &str {
        match self.help_doc {
            HelpDoc::Gui => "GUI Manual",
            HelpDoc::Cli => "CLI Manual",
            HelpDoc::AgentInterface => "Agent Interface",
            HelpDoc::ReviewerPreview => "Reviewer Quickstart",
            HelpDoc::Shell => "Shell Commands",
            HelpDoc::Tutorial => self.help_tutorial_title.as_str(),
        }
    }

    fn active_help_markdown(&self) -> &str {
        match self.help_doc {
            HelpDoc::Gui => &self.help_gui_markdown,
            HelpDoc::Cli => &self.help_cli_markdown,
            HelpDoc::AgentInterface => &self.help_agent_interface_markdown,
            HelpDoc::ReviewerPreview => &self.help_reviewer_preview_markdown,
            HelpDoc::Shell => &self.help_shell_markdown,
            HelpDoc::Tutorial => &self.help_tutorial_markdown,
        }
    }

    fn active_help_copyable_text(&self) -> String {
        self.active_help_markdown().to_string()
    }

    fn tutorial_audience_group_label(entry: &HelpTutorialDocEntry) -> String {
        if let Some(label) = entry
            .group_label
            .as_deref()
            .map(str::trim)
            .filter(|label| !label.is_empty())
        {
            return label.to_string();
        }
        let has = |needle: &str| entry.audiences.iter().any(|value| value == needle);
        if has("orientation_interfaces")
            || has("agent_users")
            || has("mcp_users")
            || has("new_users")
        {
            return "Orientation And Interfaces".to_string();
        }
        if has("pcr_qpcr_sequence_inspection") || has("primer_design") {
            return "PCR, qPCR, And Direct Sequence Inspection".to_string();
        }
        if has("protein_rna_projection") || has("protein_workflows") {
            return "Protein, RNA, And Projection Audits".to_string();
        }
        if has("sequencing_confirmation") || entry.path.contains("sequencing_confirmation") {
            return "Sequencing Confirmation".to_string();
        }
        if has("gibson_assembly_physical_layout") || entry.path.contains("gibson") {
            return "Gibson Assembly And Physical Layout".to_string();
        }
        if has("cloning_reporter_external_handoff")
            || has("cloning_planning")
            || has("synthetic_biology")
        {
            return "Cloning, Reporter Design, And External Handoffs".to_string();
        }
        "Executable Reference Chapters".to_string()
    }

    fn tutorial_audience_group_rank(entry: &HelpTutorialDocEntry, label: &str) -> usize {
        if let Some(order) = entry.group_order {
            return order.saturating_sub(1);
        }
        match label {
            "Orientation And Interfaces" => 0,
            "PCR, qPCR, And Direct Sequence Inspection" => 1,
            "Cloning, Reporter Design, And External Handoffs" => 2,
            "Gibson Assembly And Physical Layout" => 3,
            "Protein, RNA, And Projection Audits" => 4,
            "Sequencing Confirmation" => 5,
            _ => 6,
        }
    }

    fn sort_help_tutorial_entries_by_audience_group(entries: &mut [HelpTutorialDocEntry]) {
        entries.sort_by(|left, right| {
            let left_label = Self::tutorial_audience_group_label(left);
            let right_label = Self::tutorial_audience_group_label(right);
            Self::tutorial_audience_group_rank(left, &left_label)
                .cmp(&Self::tutorial_audience_group_rank(right, &right_label))
                .then_with(|| {
                    left.group_position
                        .unwrap_or(usize::MAX)
                        .cmp(&right.group_position.unwrap_or(usize::MAX))
                })
                .then_with(|| left.title.cmp(&right.title))
        });
    }

    fn tutorial_display_label(
        decimal_id: Option<&str>,
        order: Option<usize>,
        title: &str,
    ) -> String {
        if let Some(decimal_id) = decimal_id.map(str::trim).filter(|value| !value.is_empty()) {
            return format!("{decimal_id} {title}");
        }
        if let Some(order) = order {
            return format!("{order:02}. {title}");
        }
        title.to_string()
    }

    fn help_tutorial_review_label(entry: &HelpTutorialDocEntry) -> String {
        crate::workflow_examples::tutorial_review_badge_label(
            entry.review_status.as_deref(),
            entry.review_stale,
            entry.codex_reviewed_at.as_deref(),
            entry.human_reviewed_at.as_deref(),
            entry.human_reviewer.as_deref(),
        )
    }

    fn active_help_search_context(&self) -> Option<String> {
        if self.help_search_query.trim().is_empty() {
            return None;
        }
        self.help_search_matches
            .get(self.help_search_selected)
            .map(|hit| format!("line {}: {}", hit.line_number, hit.snippet))
    }

    fn selected_tutorial_feedback_context_text(&self) -> String {
        let Some(entry) = self.help_tutorial_entries.get(self.help_tutorial_selected) else {
            return "Tutorial feedback context\n-------------------------\nNo tutorial is selected.\n"
                .to_string();
        };
        let search_context = self.active_help_search_context();
        Self::tutorial_feedback_context_text_for_entry(
            entry,
            search_context.as_deref(),
            env!("CARGO_PKG_VERSION"),
            std::env::consts::OS,
        )
        .unwrap_or_else(|error| {
            format!(
                "Tutorial feedback context\n-------------------------\nTutorial title: {}\nTutorial path: {}\nCould not resolve catalog metadata: {}\n",
                entry.title, entry.path, error
            )
        })
    }

    fn tutorial_feedback_context_text_for_entry(
        entry: &HelpTutorialDocEntry,
        current_search_section: Option<&str>,
        gentle_version: &str,
        platform: &str,
    ) -> std::result::Result<String, String> {
        let repo_root = Self::detect_tutorial_repo_root().ok_or_else(|| {
            format!(
                "Could not locate '{}' from current runtime paths",
                DEFAULT_TUTORIAL_CATALOG_PATH
            )
        })?;
        let catalog_path = repo_root.join(DEFAULT_TUTORIAL_CATALOG_PATH);
        let catalog = load_tutorial_catalog(&catalog_path)?;
        let catalog_entry = catalog
            .entries
            .iter()
            .find(|candidate| {
                let resolved = repo_root
                    .join(&candidate.path)
                    .to_string_lossy()
                    .to_string();
                candidate.path == entry.path || resolved == entry.path
            })
            .ok_or_else(|| {
                format!(
                    "Selected tutorial path '{}' is not listed in '{}'",
                    entry.path,
                    catalog_path.display()
                )
            })?;
        let manifest_path = repo_root.join(&catalog.generated_runtime.manifest_path);
        let manifest = load_tutorial_manifest(&manifest_path).ok();
        let review_path = repo_root.join(DEFAULT_TUTORIAL_REVIEW_MANIFEST_PATH);
        let review_manifest = load_tutorial_review_manifest(&review_path).ok();
        let tutorial_source_dir = repo_root.join(DEFAULT_TUTORIAL_SOURCE_DIR);
        let generated_output_dir = repo_root.join(&catalog.generated_runtime.generated_readme);
        let generated_output_dir = generated_output_dir
            .parent()
            .map(Path::to_path_buf)
            .unwrap_or_else(|| repo_root.join(DEFAULT_TUTORIAL_OUTPUT_DIR));
        let workflow_example_dir = repo_root.join(DEFAULT_WORKFLOW_EXAMPLE_DIR);
        Ok(build_tutorial_feedback_context_text(
            catalog_entry,
            manifest.as_ref(),
            review_manifest.as_ref(),
            &tutorial_source_dir,
            &generated_output_dir,
            &workflow_example_dir,
            current_search_section,
            gentle_version,
            platform,
        ))
    }

    fn refresh_help_search_matches(&mut self) {
        let needle = self.help_search_query.trim();
        if needle.is_empty() {
            self.help_search_matches.clear();
            self.help_search_selected = 0;
            return;
        }
        let needle_lower = needle.to_ascii_lowercase();
        let matches = self
            .active_help_markdown()
            .lines()
            .enumerate()
            .filter_map(|(idx, line)| {
                if !line.to_ascii_lowercase().contains(&needle_lower) {
                    return None;
                }
                let snippet = line.trim();
                Some(HelpSearchMatch {
                    line_number: idx + 1,
                    snippet: if snippet.is_empty() {
                        "(empty line)".to_string()
                    } else {
                        snippet.chars().take(140).collect()
                    },
                })
            })
            .collect::<Vec<_>>();
        self.help_search_matches = matches;
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

    fn format_open_sequence_import_status(
        imported_seq_ids: &[String],
        failures: &[String],
    ) -> String {
        match (imported_seq_ids.len(), failures.len()) {
            (0, 0) => "Open sequence canceled".to_string(),
            (1, 0) => format!("Open sequence: loaded '{}'", imported_seq_ids[0]),
            (loaded, 0) => format!(
                "Open sequence: loaded {loaded} sequences ({})",
                imported_seq_ids.join(", ")
            ),
            (0, 1) => failures[0].clone(),
            (0, failed) => {
                let first = failures.first().cloned().unwrap_or_default();
                format!("Open sequence failed for {failed} files; first error: {first}")
            }
            (loaded, 1) => format!(
                "Open sequence: loaded {loaded} sequence{} ({}); 1 import failed: {}",
                if loaded == 1 { "" } else { "s" },
                imported_seq_ids.join(", "),
                failures[0]
            ),
            (loaded, failed) => {
                let first = failures.first().cloned().unwrap_or_default();
                format!(
                    "Open sequence: loaded {loaded} sequence{} ({}); {failed} imports failed (first error: {first})",
                    if loaded == 1 { "" } else { "s" },
                    imported_seq_ids.join(", ")
                )
            }
        }
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
            self.pending_viewport_focus_timestamps
                .entry(viewport_id)
                .or_insert_with(Instant::now);
        }
    }

    fn mark_viewport_open_requested(&mut self, viewport_id: ViewportId) {
        self.pending_window_open_timestamps
            .insert(viewport_id, Instant::now());
    }

    fn mark_window_open_or_focus(&mut self, viewport_id: ViewportId, was_open: bool) {
        if !was_open {
            self.mark_viewport_open_requested(viewport_id);
        }
        self.queue_focus_viewport(viewport_id);
    }

    fn viewport_foreground_requested(&self, viewport_id: ViewportId) -> bool {
        self.pending_focus_viewports.contains(&viewport_id)
            || self
                .pending_viewport_focus_timestamps
                .contains_key(&viewport_id)
    }

    fn hosted_window_spec_for_viewport(
        &self,
        title: impl Into<String>,
        stable_id: egui::Id,
        viewport_id: ViewportId,
        default_size: Vec2,
        min_size: Vec2,
    ) -> crate::egui_compat::HostedWindowSpec {
        crate::egui_compat::HostedWindowSpec::new(title, stable_id, default_size, min_size)
            .foreground(self.viewport_foreground_requested(viewport_id))
    }

    /// Clear the one-frame foreground request after an embedded hosted window
    /// has rendered. This mirrors `note_viewport_focus_if_active` for native
    /// viewports while letting the hosted window spec observe the request first.
    fn clear_viewport_foreground_request_after_render(&mut self, viewport_id: ViewportId) {
        if self.viewport_foreground_requested(viewport_id) {
            self.set_active_window_viewport(viewport_id);
            self.pending_focus_viewports.retain(|id| *id != viewport_id);
            self.finalize_viewport_focus_probe(viewport_id);
        }
    }

    fn note_slow_phase(&mut self, label: &str, elapsed_ms: u128) {
        if elapsed_ms >= WINDOW_OPEN_SLOW_THRESHOLD_MS {
            self.app_status = format!("{label} took {elapsed_ms} ms");
        }
    }

    fn note_slow_open_phase(&mut self, viewport_id: ViewportId, label: &str, elapsed_ms: u128) {
        if self
            .pending_window_open_timestamps
            .contains_key(&viewport_id)
        {
            self.note_slow_phase(label, elapsed_ms);
        }
    }

    fn viewport_focus_label(viewport_id: ViewportId) -> &'static str {
        if viewport_id == ViewportId::ROOT {
            "Main window focus acquisition"
        } else if viewport_id == Self::help_viewport_id() {
            "Help focus acquisition"
        } else if viewport_id == Self::configuration_viewport_id() {
            "Configuration focus acquisition"
        } else if viewport_id == Self::command_palette_viewport_id() {
            "Command Palette focus acquisition"
        } else if viewport_id == Self::history_viewport_id() {
            "History focus acquisition"
        } else if viewport_id == Self::prepare_genome_viewport_id() {
            "Prepare Genome focus acquisition"
        } else if viewport_id == Self::retrieve_genome_viewport_id() {
            "Retrieve Genome focus acquisition"
        } else if viewport_id == Self::blast_genome_viewport_id() {
            "BLAST Genome focus acquisition"
        } else if viewport_id == Self::bed_track_viewport_id() {
            "Track Import focus acquisition"
        } else if viewport_id == Self::gibson_viewport_id() {
            "Gibson focus acquisition"
        } else if viewport_id == Self::arrangement_gel_preview_viewport_id() {
            "Arrangement Gel focus acquisition"
        } else if viewport_id == Self::rack_labels_preview_viewport_id() {
            "Rack Labels focus acquisition"
        } else if viewport_id == Self::pcr_design_viewport_id() {
            "PCR Designer focus acquisition"
        } else if viewport_id == Self::sequencing_confirmation_viewport_id() {
            "Sequencing Confirmation focus acquisition"
        } else if viewport_id == Self::planning_viewport_id() {
            "Planning focus acquisition"
        } else if viewport_id == Self::routine_assistant_viewport_id() {
            "Routine Assistant focus acquisition"
        } else if viewport_id == Self::agent_assistant_viewport_id() {
            "Agent Assistant focus acquisition"
        } else if viewport_id == Self::jaspar_expert_viewport_id() {
            "JASPAR Expert focus acquisition"
        } else if viewport_id == Self::new_sequence_viewport_id() {
            "New Sequence focus acquisition"
        } else if viewport_id == Self::uniprot_viewport_id() {
            "UniProt focus acquisition"
        } else if viewport_id == Self::genbank_viewport_id() {
            "GenBank / dbSNP focus acquisition"
        } else {
            "Window focus acquisition"
        }
    }

    fn finalize_viewport_focus_probe(&mut self, viewport_id: ViewportId) {
        let Some(started) = self.pending_viewport_focus_timestamps.remove(&viewport_id) else {
            return;
        };
        let elapsed_ms = started.elapsed().as_millis();
        self.note_slow_phase(Self::viewport_focus_label(viewport_id), elapsed_ms);
    }

    fn native_window_menu_sync_blocked_by_open_probe(&self) -> bool {
        !self.pending_window_open_timestamps.is_empty()
    }

    fn finalize_viewport_open_probe(&mut self, viewport_id: ViewportId, label: &str) {
        let Some(started) = self.pending_window_open_timestamps.remove(&viewport_id) else {
            return;
        };
        let elapsed_ms = started.elapsed().as_millis();
        self.note_slow_phase(&format!("{label} window open"), elapsed_ms);
    }

    fn open_sequence_window_with_compact_lane_layout(
        &mut self,
        seq_id: &str,
        compact_lane_layout: bool,
    ) {
        if let Some(viewport_id) = self.find_open_sequence_viewport_id(seq_id) {
            if compact_lane_layout
                && let Some(window) = self.windows.get(&viewport_id)
                && let Ok(mut window) = window.write()
            {
                window.enable_compact_lane_layout();
            }
            self.queue_focus_viewport(viewport_id);
            return;
        }
        if let Some(window) = self.find_pending_sequence_window_mut(seq_id) {
            if compact_lane_layout {
                window.enable_compact_lane_layout();
            }
            return;
        }
        let exists = self
            .engine
            .read()
            .unwrap()
            .state()
            .sequences
            .contains_key(seq_id);
        if exists {
            let mut window = Window::new_dna_lazy(seq_id.to_string(), self.engine.clone());
            if compact_lane_layout {
                window.enable_compact_lane_layout();
            }
            self.new_windows.push(window);
        }
    }

    fn open_sequence_window(&mut self, seq_id: &str) {
        self.open_sequence_window_with_compact_lane_layout(seq_id, false);
    }

    fn open_sequence_window_compact(&mut self, seq_id: &str) {
        self.open_sequence_window_with_compact_lane_layout(seq_id, true);
    }

    fn open_sequence_window_for_dotplot_analysis(&mut self, seq_id: &str, dotplot_id: &str) {
        if dotplot_id.trim().is_empty() {
            self.open_sequence_window(seq_id);
            return;
        }
        if let Some(viewport_id) = self.find_open_sequence_viewport_id(seq_id) {
            if let Some(window) = self.windows.get(&viewport_id)
                && let Ok(mut window) = window.write()
            {
                window.focus_dotplot_analysis(dotplot_id);
            }
            self.queue_focus_viewport(viewport_id);
            return;
        }
        if let Some(window) = self.find_pending_sequence_window_mut(seq_id) {
            window.focus_dotplot_analysis(dotplot_id);
            return;
        }
        let exists = self
            .engine
            .read()
            .unwrap()
            .state()
            .sequences
            .contains_key(seq_id);
        if exists {
            let mut window = Window::new_dna_lazy(seq_id.to_string(), self.engine.clone());
            window.focus_dotplot_analysis(dotplot_id);
            self.new_windows.push(window);
        }
    }

    fn open_sequence_window_for_flexibility_track_analysis(
        &mut self,
        seq_id: &str,
        track_id: &str,
    ) {
        if track_id.trim().is_empty() {
            self.open_sequence_window(seq_id);
            return;
        }
        if let Some(viewport_id) = self.find_open_sequence_viewport_id(seq_id) {
            if let Some(window) = self.windows.get(&viewport_id)
                && let Ok(mut window) = window.write()
            {
                window.focus_flexibility_track_analysis(track_id);
            }
            self.queue_focus_viewport(viewport_id);
            return;
        }
        if let Some(window) = self.find_pending_sequence_window_mut(seq_id) {
            window.focus_flexibility_track_analysis(track_id);
            return;
        }
        let exists = self
            .engine
            .read()
            .unwrap()
            .state()
            .sequences
            .contains_key(seq_id);
        if exists {
            let mut window = Window::new_dna_lazy(seq_id.to_string(), self.engine.clone());
            window.focus_flexibility_track_analysis(track_id);
            self.new_windows.push(window);
        }
    }

    fn open_sequence_window_for_sequencing_confirmation_report(
        &mut self,
        seq_id: &str,
        report_id: &str,
    ) {
        if !report_id.trim().is_empty() {
            if let Some(viewport_id) = self.find_open_sequence_viewport_id(seq_id) {
                if let Some(window) = self.windows.get(&viewport_id)
                    && let Ok(mut window) = window.write()
                {
                    window.focus_sequencing_confirmation_report(report_id);
                }
            } else if let Some(window) = self.find_pending_sequence_window_mut(seq_id) {
                window.focus_sequencing_confirmation_report(report_id);
            } else {
                let exists = self
                    .engine
                    .read()
                    .unwrap()
                    .state()
                    .sequences
                    .contains_key(seq_id);
                if exists {
                    let mut window = Window::new_dna_lazy(seq_id.to_string(), self.engine.clone());
                    window.focus_sequencing_confirmation_report(report_id);
                    self.new_windows.push(window);
                }
            }
        } else {
            self.open_sequence_window(seq_id);
        }
        self.sequencing_confirmation_seq_id = seq_id.to_string();
        self.show_sequencing_confirmation_dialog = true;
        self.queue_focus_viewport(Self::sequencing_confirmation_viewport_id());
    }

    fn open_sequence_window_for_protein_derivation_report(
        &mut self,
        seq_id: &str,
        report_id: &str,
    ) {
        let transcript_filter = if report_id.trim().is_empty() {
            None
        } else {
            {
                self.engine
                    .read()
                    .unwrap()
                    .get_protein_derivation_report(report_id)
                    .ok()
                    .and_then(|report| {
                        if report.seq_id.eq_ignore_ascii_case(seq_id) && report.rows.len() == 1 {
                            report
                                .rows
                                .first()
                                .map(|row| row.derivation.transcript_id.clone())
                        } else {
                            None
                        }
                    })
            }
        };
        self.open_sequence_window_for_transcript_protein_expert(
            seq_id,
            transcript_filter.as_deref(),
        );
    }

    fn open_sequence_window_for_reverse_translation_report(
        &mut self,
        seq_id: &str,
        report_id: &str,
    ) {
        let coding_seq_id = if report_id.trim().is_empty() {
            None
        } else {
            {
                self.engine
                    .read()
                    .unwrap()
                    .get_reverse_translation_report(report_id)
                    .ok()
                    .and_then(|report| {
                        if report.protein_seq_id.eq_ignore_ascii_case(seq_id) {
                            Some(report.coding_seq_id)
                        } else {
                            None
                        }
                    })
            }
        };
        if let Some(coding_seq_id) = coding_seq_id.as_deref() {
            self.open_sequence_window(coding_seq_id);
        } else {
            self.open_sequence_window(seq_id);
        }
    }

    fn open_sequence_window_for_construct_reasoning_graph(&mut self, seq_id: &str, graph_id: &str) {
        if !graph_id.trim().is_empty() {
            if let Some(viewport_id) = self.find_open_sequence_viewport_id(seq_id) {
                if let Some(window) = self.windows.get(&viewport_id)
                    && let Ok(mut window) = window.write()
                {
                    window.focus_construct_reasoning_graph(graph_id);
                }
            } else if let Some(window) = self.find_pending_sequence_window_mut(seq_id) {
                window.focus_construct_reasoning_graph(graph_id);
            } else {
                let exists = self
                    .engine
                    .read()
                    .unwrap()
                    .state()
                    .sequences
                    .contains_key(seq_id);
                if exists {
                    let mut window = Window::new_dna_lazy(seq_id.to_string(), self.engine.clone());
                    window.focus_construct_reasoning_graph(graph_id);
                    self.new_windows.push(window);
                }
            }
        } else {
            self.open_sequence_window(seq_id);
        }
    }

    fn open_protein_to_dna_handoff_dialog_for_graph(&mut self, seq_id: &str, graph_id: &str) {
        let graph = match self
            .engine
            .read()
            .unwrap()
            .construct_reasoning_graph(graph_id)
        {
            Ok(graph) => graph,
            Err(err) => {
                self.protein_handoff_status = format!(
                    "Could not reopen protein-to-DNA handoff graph '{}': {}",
                    graph_id, err.message
                );
                self.show_uniprot_dialog = true;
                return;
            }
        };
        if !GentleEngine::construct_reasoning_graph_has_protein_to_dna_handoff(&graph) {
            self.open_sequence_window_for_construct_reasoning_graph(seq_id, graph_id);
            return;
        }
        self.show_uniprot_dialog = true;
        self.uniprot_map_seq_id = seq_id.to_string();
        self.protein_handoff_graph = Some(graph.clone());
        self.sync_protein_handoff_candidate_selection();
        if let Some(candidate) = graph
            .candidates
            .iter()
            .find(|candidate| candidate.candidate_id == self.protein_handoff_selected_candidate_id)
            .or_else(|| {
                graph
                    .candidates
                    .iter()
                    .find(|candidate| candidate.protein_to_dna_handoff.is_some())
            })
            && let Some(detail) = candidate.protein_to_dna_handoff.as_ref()
        {
            self.reverse_translate_protein_seq_id = detail.source_protein_seq_id.clone();
            self.protein_handoff_ranking_goal = detail.ranking_goal;
            if let Some(transcript_id) = detail.transcript_id.as_ref() {
                self.uniprot_map_transcript_id = transcript_id.clone();
                self.uniprot_feature_transcript_id = transcript_id.clone();
            }
            if let Some(projection_id) = detail.projection_id.as_ref() {
                self.uniprot_map_projection_id = projection_id.clone();
            }
            if let Some(feature_query) = detail.feature_query.as_ref() {
                self.uniprot_feature_query = feature_query.clone();
            }
            if let Some(entry_id) = detail.ensembl_entry_id.as_ref() {
                self.ensembl_protein_entry_id = entry_id.clone();
            }
        }
        self.protein_handoff_status = format!(
            "Opened protein-to-DNA handoff reasoning graph '{}' in Protein Evidence",
            graph_id
        );
    }

    fn open_sequence_window_for_primer_design_report(&mut self, seq_id: &str, report_id: &str) {
        if !report_id.trim().is_empty() {
            if let Some(viewport_id) = self.find_open_sequence_viewport_id(seq_id) {
                if let Some(window) = self.windows.get(&viewport_id)
                    && let Ok(mut window) = window.write()
                {
                    window.focus_primer_design_report(report_id);
                }
            } else if let Some(window) = self.find_pending_sequence_window_mut(seq_id) {
                window.focus_primer_design_report(report_id);
            } else {
                let exists = self
                    .engine
                    .read()
                    .unwrap()
                    .state()
                    .sequences
                    .contains_key(seq_id);
                if exists {
                    let mut window = Window::new_dna_lazy(seq_id.to_string(), self.engine.clone());
                    window.focus_primer_design_report(report_id);
                    self.new_windows.push(window);
                }
            }
        } else {
            self.open_sequence_window(seq_id);
        }
        if let Err(err) = self.open_pcr_design_dialog_for_seq_id(seq_id) {
            self.app_status = format!(
                "Cannot open PCR Designer for primer report '{}' on '{}': {}",
                report_id, seq_id, err
            );
        }
    }

    fn open_sequence_window_for_qpcr_design_report(&mut self, seq_id: &str, report_id: &str) {
        if !report_id.trim().is_empty() {
            if let Some(viewport_id) = self.find_open_sequence_viewport_id(seq_id) {
                if let Some(window) = self.windows.get(&viewport_id)
                    && let Ok(mut window) = window.write()
                {
                    window.focus_qpcr_design_report(report_id);
                }
            } else if let Some(window) = self.find_pending_sequence_window_mut(seq_id) {
                window.focus_qpcr_design_report(report_id);
            } else {
                let exists = self
                    .engine
                    .read()
                    .unwrap()
                    .state()
                    .sequences
                    .contains_key(seq_id);
                if exists {
                    let mut window = Window::new_dna_lazy(seq_id.to_string(), self.engine.clone());
                    window.focus_qpcr_design_report(report_id);
                    self.new_windows.push(window);
                }
            }
        } else {
            self.open_sequence_window(seq_id);
        }
        if let Err(err) = self.open_pcr_design_dialog_for_seq_id(seq_id) {
            self.app_status = format!(
                "Cannot open PCR Designer for qPCR report '{}' on '{}': {}",
                report_id, seq_id, err
            );
        }
    }

    fn open_sequence_window_for_restriction_cloning_pcr_handoff(
        &mut self,
        seq_id: &str,
        report_id: &str,
    ) {
        if !report_id.trim().is_empty() {
            if let Some(viewport_id) = self.find_open_sequence_viewport_id(seq_id) {
                if let Some(window) = self.windows.get(&viewport_id)
                    && let Ok(mut window) = window.write()
                {
                    window.focus_restriction_cloning_handoff_report(report_id);
                }
            } else if let Some(window) = self.find_pending_sequence_window_mut(seq_id) {
                window.focus_restriction_cloning_handoff_report(report_id);
            } else {
                let exists = self
                    .engine
                    .read()
                    .unwrap()
                    .state()
                    .sequences
                    .contains_key(seq_id);
                if exists {
                    let mut window = Window::new_dna_lazy(seq_id.to_string(), self.engine.clone());
                    window.focus_restriction_cloning_handoff_report(report_id);
                    self.new_windows.push(window);
                }
            }
        } else {
            self.open_sequence_window(seq_id);
        }
        if let Err(err) = self.open_pcr_design_dialog_for_seq_id(seq_id) {
            self.app_status = format!(
                "Cannot open PCR Designer for restriction-cloning handoff '{}' on '{}': {}",
                report_id, seq_id, err
            );
        }
    }

    fn open_sequence_window_for_rna_read_report(&mut self, seq_id: &str, report_id: &str) {
        if !report_id.trim().is_empty() {
            if let Some(viewport_id) = self.find_open_sequence_viewport_id(seq_id) {
                if let Some(window) = self.windows.get(&viewport_id)
                    && let Ok(mut window) = window.write()
                {
                    window.focus_rna_read_report(report_id);
                }
            } else if let Some(window) = self.find_pending_sequence_window_mut(seq_id) {
                window.focus_rna_read_report(report_id);
            } else {
                let exists = self
                    .engine
                    .read()
                    .unwrap()
                    .state()
                    .sequences
                    .contains_key(seq_id);
                if exists {
                    let mut window = Window::new_dna_lazy(seq_id.to_string(), self.engine.clone());
                    window.focus_rna_read_report(report_id);
                    self.new_windows.push(window);
                }
            }
        } else {
            self.open_sequence_window(seq_id);
        }
    }

    fn open_lineage_analysis_artifact(
        &mut self,
        kind: LineageAnalysisKind,
        seq_id: &str,
        artifact_id: &str,
    ) {
        match kind {
            LineageAnalysisKind::Dotplot => {
                self.open_sequence_window_for_dotplot_analysis(seq_id, artifact_id);
            }
            LineageAnalysisKind::FlexibilityTrack => {
                self.open_sequence_window_for_flexibility_track_analysis(seq_id, artifact_id);
            }
            LineageAnalysisKind::RnaReadInterpretation => {
                self.open_sequence_window_for_rna_read_report(seq_id, artifact_id);
            }
            LineageAnalysisKind::PrimerDesign => {
                self.open_sequence_window_for_primer_design_report(seq_id, artifact_id);
            }
            LineageAnalysisKind::QpcrDesign => {
                self.open_sequence_window_for_qpcr_design_report(seq_id, artifact_id);
            }
            LineageAnalysisKind::RestrictionCloningPcrHandoff => {
                self.open_sequence_window_for_restriction_cloning_pcr_handoff(seq_id, artifact_id);
            }
            LineageAnalysisKind::ProteinDerivation => {
                self.open_sequence_window_for_protein_derivation_report(seq_id, artifact_id);
            }
            LineageAnalysisKind::ReverseTranslation => {
                self.open_sequence_window_for_reverse_translation_report(seq_id, artifact_id);
            }
            LineageAnalysisKind::ConstructReasoning => {
                self.open_protein_to_dna_handoff_dialog_for_graph(seq_id, artifact_id);
            }
            LineageAnalysisKind::UniprotProjection => {
                self.open_sequence_window_for_uniprot_projection_expert(
                    seq_id,
                    artifact_id,
                    Default::default(),
                );
            }
            LineageAnalysisKind::SequencingConfirmation => {
                self.open_sequence_window_for_sequencing_confirmation_report(seq_id, artifact_id);
            }
        }
    }

    fn open_pool_window_with_compact_lane_layout(
        &mut self,
        representative_seq_id: &str,
        pool_seq_ids: Vec<String>,
        compact_lane_layout: bool,
    ) {
        if let Some(viewport_id) = self.find_open_sequence_viewport_id(representative_seq_id) {
            if let Some(window) = self.windows.get(&viewport_id)
                && let Ok(mut window) = window.write()
            {
                window.set_pool_context(pool_seq_ids);
                if compact_lane_layout {
                    window.enable_compact_lane_layout();
                }
            }
            self.queue_focus_viewport(viewport_id);
            return;
        }
        if let Some(window) = self.find_pending_sequence_window_mut(representative_seq_id) {
            window.set_pool_context(pool_seq_ids);
            if compact_lane_layout {
                window.enable_compact_lane_layout();
            }
            return;
        }
        let exists = self
            .engine
            .read()
            .unwrap()
            .state()
            .sequences
            .contains_key(representative_seq_id);
        if exists {
            let mut window =
                Window::new_dna_lazy(representative_seq_id.to_string(), self.engine.clone());
            window.set_pool_context(pool_seq_ids);
            if compact_lane_layout {
                window.enable_compact_lane_layout();
            }
            self.new_windows.push(window);
        }
    }

    fn open_pool_window(&mut self, representative_seq_id: &str, pool_seq_ids: Vec<String>) {
        self.open_pool_window_with_compact_lane_layout(representative_seq_id, pool_seq_ids, false);
    }

    fn open_pool_window_compact(&mut self, representative_seq_id: &str, pool_seq_ids: Vec<String>) {
        self.open_pool_window_with_compact_lane_layout(representative_seq_id, pool_seq_ids, true);
    }

    fn arrangement_lane_container_rows<'a>(
        &'a self,
        lane_container_ids: &[String],
    ) -> Vec<&'a ContainerRow> {
        lane_container_ids
            .iter()
            .filter_map(|container_id| {
                self.lineage_containers
                    .iter()
                    .find(|row| row.container_id == *container_id)
            })
            .collect()
    }

    fn arrangement_lane_menu_label(row: &ContainerRow) -> String {
        if row.member_count > 1 {
            format!(
                "{} ({}, {} seqs)",
                row.container_id, row.kind, row.member_count
            )
        } else if row.representative.trim().is_empty() {
            row.container_id.clone()
        } else {
            format!("{} ({})", row.container_id, row.representative)
        }
    }

    fn container_contents_mode_label(exclusive: bool) -> &'static str {
        if exclusive {
            "Declared only"
        } else {
            "Known subset"
        }
    }

    fn container_contents_mode_hover(exclusive: bool) -> &'static str {
        if exclusive {
            "This container is declared exhaustive: the listed members are intended to be the whole known contents."
        } else {
            "This container is non-exclusive: the listed members are known/measured constituents, but other compounds may also be present."
        }
    }

    fn arrangement_lane_menu_hover(row: &ContainerRow) -> String {
        if row.member_count > 1 {
            format!(
                "Open this pooled lane as a compact DNA window with the sequence text panel hidden by default.\nRepresentative: {}\nMembers: {}\nContents: {}",
                row.representative,
                row.members.join(", "),
                Self::container_contents_mode_label(row.declared_contents_exclusive)
            )
        } else {
            format!(
                "Open this lane as a compact DNA window with the sequence text panel hidden by default.\nSequence: {}\nContents: {}",
                row.representative,
                Self::container_contents_mode_label(row.declared_contents_exclusive)
            )
        }
    }

    fn render_open_lanes_menu(
        &self,
        ui: &mut egui::Ui,
        lane_container_ids: &[String],
        open_lane_containers: &mut Option<Vec<String>>,
    ) {
        let lane_rows = self.arrangement_lane_container_rows(lane_container_ids);
        if lane_rows.is_empty() {
            ui.label("-");
            return;
        }
        ui.menu_button("Open Lanes", |ui| {
            if ui
                .button(format!("All ({})", lane_rows.len()))
                .on_hover_text(
                    "Open all arrangement lanes as compact DNA windows with the sequence text panel hidden by default.",
                )
                .clicked()
            {
                *open_lane_containers =
                    Some(lane_rows.iter().map(|row| row.container_id.clone()).collect());
                ui.close();
            }
            ui.separator();
            for row in &lane_rows {
                if ui
                    .button(Self::arrangement_lane_menu_label(row))
                    .on_hover_text(Self::arrangement_lane_menu_hover(row))
                    .clicked()
                {
                    *open_lane_containers = Some(vec![row.container_id.clone()]);
                    ui.close();
                }
            }
        })
        .response
        .on_hover_text(
            "Choose one arrangement lane to open, or open all lanes at once as compact DNA windows.",
        );
    }

    fn sanitize_fasta_header(raw: &str) -> String {
        let mut header: String = raw
            .chars()
            .map(|ch| if ch.is_whitespace() { '_' } else { ch })
            .collect();
        if header.trim().is_empty() {
            header = "sequence".to_string();
        }
        header
    }

    fn build_fasta_payload(header: &str, sequence: &str) -> String {
        let mut output = format!(">{header}\n");
        for chunk in sequence.as_bytes().chunks(80) {
            output.push_str(std::str::from_utf8(chunk).unwrap_or_default());
            output.push('\n');
        }
        output
    }

    fn sequence_copy_payload(
        &self,
        seq_id: &str,
        fasta: bool,
    ) -> std::result::Result<(String, usize), String> {
        let guard = self
            .engine
            .read()
            .map_err(|_| "Engine lock poisoned".to_string())?;
        let dna = guard
            .state()
            .sequences
            .get(seq_id)
            .ok_or_else(|| format!("Sequence '{seq_id}' is not available in project state"))?;
        let plain = dna.get_forward_string();
        let bp = plain.len();
        if fasta {
            let header = Self::sanitize_fasta_header(seq_id);
            Ok((Self::build_fasta_payload(&header, &plain), bp))
        } else {
            Ok((plain, bp))
        }
    }

    fn copy_project_sequence_to_clipboard(
        &mut self,
        ctx: &egui::Context,
        seq_id: &str,
        fasta: bool,
    ) {
        match self.sequence_copy_payload(seq_id, fasta) {
            Ok((payload, bp)) => {
                ctx.copy_text(payload);
                if fasta {
                    self.app_status = format!("Copied FASTA for '{seq_id}' ({bp} bp)");
                } else {
                    self.app_status = format!("Copied sequence for '{seq_id}' ({bp} bp)");
                }
            }
            Err(err) => {
                self.app_status = format!("Could not copy '{seq_id}': {err}");
            }
        }
    }

    pub fn load_from_file(path: &str) -> Result<DNAsequence> {
        dna_sequence::load_from_file(path)
    }

    fn import_sequence_file_into_project(
        &mut self,
        path: &str,
    ) -> std::result::Result<Vec<String>, String> {
        let op = Operation::LoadFile {
            path: path.to_string(),
            as_id: None,
        };
        let load_result = {
            let mut engine = self.engine.write().unwrap();
            engine.apply(op)
        };
        match load_result {
            Ok(result) => {
                let loaded_sequences = {
                    let engine = self.engine.read().unwrap();
                    result
                        .created_seq_ids
                        .iter()
                        .filter_map(|seq_id| {
                            engine
                                .state()
                                .sequences
                                .get(seq_id)
                                .cloned()
                                .map(|dna| (seq_id.clone(), dna))
                        })
                        .collect::<Vec<_>>()
                };
                if loaded_sequences.is_empty() {
                    return Err(Self::format_op_result_status(
                        "Open sequence: operation finished but no sequence window was created",
                        &result.created_seq_ids,
                        &result.warnings,
                        &result.messages,
                    ));
                }
                let imported_seq_ids = loaded_sequences
                    .iter()
                    .map(|(seq_id, _)| seq_id.clone())
                    .collect::<Vec<_>>();
                for (seq_id, dna) in loaded_sequences {
                    self.new_dna_window(seq_id, dna);
                }
                Ok(imported_seq_ids)
            }
            Err(err) => Err(format!("Open sequence failed: {}", err.message)),
        }
    }

    fn open_new_window_from_file(&mut self, path: &str) {
        self.app_status = match self.import_sequence_file_into_project(path) {
            Ok(imported_seq_ids) => {
                Self::format_open_sequence_import_status(&imported_seq_ids, &[])
            }
            Err(err) => err,
        };
    }

    fn open_new_windows_from_files(&mut self, paths: &[PathBuf]) {
        if paths.is_empty() {
            self.app_status = "Open sequence canceled".to_string();
            return;
        }
        if paths.len() == 1 {
            let path_text = paths[0].display().to_string();
            self.open_new_window_from_file(&path_text);
            return;
        }
        let mut imported_seq_ids: Vec<String> = vec![];
        let mut failures: Vec<String> = vec![];
        for path in paths {
            let path_text = path.display().to_string();
            match self.import_sequence_file_into_project(&path_text) {
                Ok(mut seq_ids) => imported_seq_ids.append(&mut seq_ids),
                Err(err) => failures.push(err),
            }
        }
        self.app_status = Self::format_open_sequence_import_status(&imported_seq_ids, &failures);
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
        state.container_state.racks.len().hash(&mut hasher);
        state.metadata.len().hash(&mut hasher);

        let mut metadata_keys: Vec<&String> = state.metadata.keys().collect();
        metadata_keys.sort_unstable();
        for key in metadata_keys {
            key.hash(&mut hasher);
        }

        let display = &state.display;
        display.show_sequence_panel.hash(&mut hasher);
        display.show_linear_sequence_panel.hash(&mut hasher);
        display.sequence_panel_max_text_length_bp.hash(&mut hasher);
        display.show_map_panel.hash(&mut hasher);
        display.show_cds_features.hash(&mut hasher);
        display.show_gene_features.hash(&mut hasher);
        display.show_mrna_features.hash(&mut hasher);
        display.show_repeat_features.hash(&mut hasher);
        display.show_array_features.hash(&mut hasher);
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
        display.restriction_enzyme_display_mode.hash(&mut hasher);
        for name in &display.preferred_restriction_enzymes {
            name.hash(&mut hasher);
        }
        display.show_gc_contents.hash(&mut hasher);
        display.gc_content_bin_size_bp.hash(&mut hasher);
        display.show_open_reading_frames.hash(&mut hasher);
        display.show_methylation_sites.hash(&mut hasher);
        display
            .feature_details_font_size
            .to_bits()
            .hash(&mut hasher);
        display
            .linear_external_feature_label_font_size
            .to_bits()
            .hash(&mut hasher);
        display
            .linear_external_feature_label_background_opacity
            .to_bits()
            .hash(&mut hasher);
        display
            .auto_hide_sequence_panel_when_linear_bases_visible
            .hash(&mut hasher);
        display.linear_view_start_bp.hash(&mut hasher);
        display.linear_view_span_bp.hash(&mut hasher);
        display
            .linear_view_vertical_offset_px
            .to_bits()
            .hash(&mut hasher);
        display
            .linear_sequence_base_text_max_view_span_bp
            .hash(&mut hasher);
        display
            .linear_sequence_helical_letters_enabled
            .hash(&mut hasher);
        display
            .linear_sequence_helical_max_view_span_bp
            .hash(&mut hasher);
        display
            .linear_sequence_condensed_max_view_span_bp
            .hash(&mut hasher);
        display.linear_sequence_letter_layout_mode.hash(&mut hasher);
        display
            .linear_sequence_helical_phase_offset_bp
            .hash(&mut hasher);
        display.linear_show_double_strand_bases.hash(&mut hasher);
        display.linear_helical_parallel_strands.hash(&mut hasher);
        display
            .linear_hide_backbone_when_sequence_bases_visible
            .hash(&mut hasher);
        display
            .linear_reverse_strand_use_upside_down_letters
            .hash(&mut hasher);
        display
            .reverse_strand_visual_opacity
            .to_bits()
            .hash(&mut hasher);

        hasher.finish()
    }

    fn current_display_change_stamp(&self) -> u64 {
        let display_bytes = {
            let engine = self.engine.read().unwrap();
            serde_json::to_vec(&engine.state().display).unwrap_or_default()
        };
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        display_bytes.hash(&mut hasher);
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
        state.lineage.macro_instances.len().hash(&mut hasher);
        state.container_state.containers.len().hash(&mut hasher);
        state.container_state.arrangements.len().hash(&mut hasher);
        state.container_state.racks.len().hash(&mut hasher);
        for report in GentleEngine::sequencing_confirmation_reports_from_state(state) {
            report.report_id.hash(&mut hasher);
            report.expected_seq_id.hash(&mut hasher);
            report.baseline_seq_id.hash(&mut hasher);
            report.overall_status.as_str().hash(&mut hasher);
            report.generated_at_unix_ms.hash(&mut hasher);
        }
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
            || !state.container_state.racks.is_empty()
    }

    fn can_close_project(&self) -> bool {
        self.current_project_path.is_some() || self.project_has_user_content()
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
        self.lineage_racks.clear();
        self.lineage_graph_view = true;
        self.lineage_main_split_fraction = DEFAULT_LINEAGE_MAIN_SPLIT_FRACTION;
        self.lineage_graph_zoom = 1.0;
        self.lineage_graph_area_height = 420.0;
        self.lineage_container_area_height = 220.0;
        self.lineage_container_arrangement_split_fraction =
            DEFAULT_LINEAGE_CONTAINER_ARRANGEMENT_SPLIT_FRACTION;
        self.lineage_graph_scroll_offset = Vec2::ZERO;
        self.lineage_graph_pan_origin = None;
        self.lineage_graph_compact_labels = true;
        self.lineage_main_split_drag_origin = None;
        self.lineage_main_split_drag_start_y = None;
        self.lineage_container_arrangement_split_drag_origin = None;
        self.lineage_container_arrangement_split_drag_start_y = None;
        self.lineage_graph_selected_node_id = None;
        self.lineage_graph_node_offsets.clear();
        self.lineage_graph_drag_origin = None;
        self.lineage_graph_offsets_synced_stamp = 0;
        self.lineage_node_rename_target_id = None;
        self.lineage_node_rename_text.clear();
        self.lineage_node_remove_target_id = None;
        self.new_windows.clear();
        self.windows.clear();
        self.detached_auxiliary_window_hosts.clear();
        self.windows_to_close.write().unwrap().clear();
        self.pending_focus_viewports.clear();
        self.viewport_id_counter = 0;
        self.pending_project_action = None;
        self.pending_app_quit = false;
        self.show_reference_genome_blast_dialog = false;
        self.genome_blast_task = None;
        self.genome_blast_progress_fraction = None;
        self.genome_blast_progress_label.clear();
        self.genome_blast_results.clear();
        self.genome_blast_selected_result = 0;
        self.genome_blast_status.clear();
        self.show_genome_bed_track_dialog = false;
        self.show_gibson_dialog = false;
        self.show_arrangement_gel_preview_dialog = false;
        self.show_rack_labels_preview_dialog = false;
        self.show_rack_dialog = false;
        self.show_place_arrangement_rack_dialog = false;
        self.show_pcr_design_dialog = false;
        self.show_sequencing_confirmation_dialog = false;
        self.show_routine_assistant_dialog = false;
        self.show_agent_assistant_dialog = false;
        self.show_jaspar_expert_dialog = false;
        self.show_uniprot_dialog = false;
        self.show_genbank_dialog = false;
        self.pcr_design_seq_id.clear();
        self.sequencing_confirmation_seq_id.clear();
        self.gibson_destination_seq_id.clear();
        self.gibson_opening_start_0based.clear();
        self.gibson_opening_end_0based_exclusive.clear();
        self.gibson_insert_seq_id.clear();
        self.gibson_extra_inserts.clear();
        self.gibson_output_id_hint.clear();
        self.gibson_show_all_unique_cutters = false;
        self.rack_labels_preview = RackLabelsPreviewState::default();
        self.gibson_status.clear();
        self.gibson_preview_output = None;
        self.gibson_preview_svg_uri.clear();
        self.arrangement_gel_preview = ArrangementGelPreviewState::default();
        self.rack_view_rack_id.clear();
        self.rack_view_status.clear();
        self.rack_view_scroll_offset = Vec2::ZERO;
        self.rack_authoring_template_editor = RackAuthoringTemplate::BenchRows;
        self.rack_fill_direction_editor = RackFillDirection::RowMajor;
        self.rack_blocked_coordinates_text.clear();
        self.rack_view_selected_coordinates.clear();
        self.rack_view_selected_arrangement_ids.clear();
        self.rack_view_drag_state = None;
        self.rack_view_hover_target_coordinate = None;
        self.rack_view_recent_drop_ghost = None;
        self.place_arrangement_source_id.clear();
        self.place_arrangement_target_rack_id.clear();
        self.place_arrangement_status.clear();
        self.routine_assistant_stage = RoutineAssistantStage::GoalAndCandidates;
        self.routine_assistant_goal.clear();
        self.routine_assistant_query.clear();
        self.routine_assistant_candidates.clear();
        self.routine_assistant_preference_context = None;
        self.routine_assistant_macro_suggestions.clear();
        self.routine_assistant_selected_routine_id.clear();
        self.routine_assistant_compare_routine_id.clear();
        self.routine_assistant_bindings.clear();
        self.routine_assistant_disambiguation_answers.clear();
        self.routine_assistant_explain_output = None;
        self.routine_assistant_compare_output = None;
        self.routine_assistant_preflight_output = None;
        self.routine_assistant_execute_output = None;
        self.routine_assistant_status.clear();
        self.routine_assistant_decision_trace = None;
        self.routine_assistant_trace_counter = 1;
        self.genome_track_status.clear();
        self.genome_track_import_task = None;
        self.genome_track_import_progress = None;
        self.genome_track_autosync_status.clear();
        self.tracked_autosync_last_op_count = None;
        self.tutorial_project_task = None;
        self.tutorial_project_progress_fraction = None;
        self.tutorial_project_progress_label.clear();
        self.tutorial_project_status.clear();
        self.agent_task = None;
        self.agent_model_discovery_task = None;
        self.agent_status.clear();
        self.agent_preflight_output = None;
        self.agent_last_invocation = None;
        self.agent_execution_log.clear();
        self.agent_discovered_models.clear();
        self.agent_discovered_model_pick.clear();
        self.agent_model_discovery_status.clear();
        self.agent_model_discovery_source_key.clear();
        self.agent_model_discovery_failed_source_key.clear();
        self.uniprot_status.clear();
        self.genbank_status.clear();
        self.load_bed_track_subscriptions_from_state();
        self.load_lineage_graph_workspace_from_state();
        self.load_rack_workspace_from_state();
        self.load_background_job_history_from_state();
        self.mark_clean_snapshot();
    }

    fn set_current_project_path(&mut self, path: &str, track_recent: bool) {
        let normalized = Self::normalize_project_path(path);
        if normalized.is_empty() {
            self.current_project_path = None;
            return;
        }
        self.current_project_path = Some(normalized.clone());
        if track_recent {
            self.record_recent_project_path(&normalized);
        }
    }

    fn set_current_project_path_and_track_recent(&mut self, path: &str) {
        self.set_current_project_path(path, true);
    }

    fn open_project_path(&mut self, path: &str) {
        if let Err(err) = self.load_project_from_file(path) {
            eprintln!("Could not open project '{}': {err}", path);
            self.remove_recent_project_path(path);
        }
    }

    fn tutorial_cache_root_dir() -> PathBuf {
        env::temp_dir().join("gentle_tutorial_projects")
    }

    fn detect_tutorial_repo_root() -> Option<PathBuf> {
        let mut candidates: Vec<PathBuf> = vec![];
        if let Ok(cwd) = env::current_dir() {
            candidates.push(cwd);
        }
        if let Ok(exe) = env::current_exe() {
            let mut cursor = exe.parent().map(|p| p.to_path_buf());
            while let Some(path) = cursor {
                candidates.push(path.clone());
                cursor = path.parent().map(|p| p.to_path_buf());
            }
        }
        for candidate in candidates {
            let catalog = candidate.join(DEFAULT_TUTORIAL_CATALOG_PATH);
            let manifest = candidate.join(DEFAULT_TUTORIAL_MANIFEST_PATH);
            let examples = candidate.join(DEFAULT_WORKFLOW_EXAMPLE_DIR);
            if catalog.is_file() || (manifest.is_file() && examples.is_dir()) {
                return Some(candidate);
            }
        }
        None
    }

    fn tutorial_catalog_runtime_paths() -> std::result::Result<TutorialCatalogRuntimePaths, String>
    {
        let repo_root = Self::detect_tutorial_repo_root().ok_or_else(|| {
            format!(
                "Could not locate '{}' or '{}' from current runtime paths",
                DEFAULT_TUTORIAL_CATALOG_PATH, DEFAULT_TUTORIAL_MANIFEST_PATH
            )
        })?;
        let catalog_path = repo_root.join(DEFAULT_TUTORIAL_CATALOG_PATH);
        let examples_dir = repo_root.join(DEFAULT_WORKFLOW_EXAMPLE_DIR);
        if !examples_dir.is_dir() {
            return Err(format!(
                "Could not locate workflow examples directory '{}'",
                examples_dir.display()
            ));
        }
        if catalog_path.is_file() {
            let catalog = load_tutorial_catalog(&catalog_path)?;
            let manifest_path = repo_root.join(&catalog.generated_runtime.manifest_path);
            if !manifest_path.is_file() {
                return Err(format!(
                    "Tutorial catalog references missing manifest '{}'",
                    manifest_path.display()
                ));
            }
            return Ok(TutorialCatalogRuntimePaths {
                repo_root,
                manifest_path,
                examples_dir,
            });
        }
        let manifest_path = repo_root.join(DEFAULT_TUTORIAL_MANIFEST_PATH);
        if !manifest_path.is_file() {
            return Err(format!(
                "Could not locate tutorial manifest '{}'",
                manifest_path.display()
            ));
        }
        Ok(TutorialCatalogRuntimePaths {
            repo_root,
            manifest_path,
            examples_dir,
        })
    }

    fn load_tutorial_project_entries() -> std::result::Result<Vec<TutorialProjectEntry>, String> {
        let runtime_paths = Self::tutorial_catalog_runtime_paths()?;
        let repo_root = runtime_paths.repo_root;
        let manifest_path = runtime_paths.manifest_path;
        let examples_dir = runtime_paths.examples_dir;
        let manifest = load_tutorial_manifest(&manifest_path)?;
        let loaded_examples = load_workflow_examples(&examples_dir)?;
        validate_tutorial_manifest_against_examples(&manifest, &loaded_examples)?;
        let review_by_id = Self::resolve_runtime_doc_path(DEFAULT_TUTORIAL_CATALOG_PATH)
            .and_then(|path| load_tutorial_catalog(&path).ok())
            .map(|catalog| {
                catalog
                    .entries
                    .into_iter()
                    .map(|entry| (entry.id.clone(), entry))
                    .collect::<HashMap<_, _>>()
            })
            .unwrap_or_default();
        let mut by_id: HashMap<String, WorkflowExample> = HashMap::new();
        for loaded in loaded_examples {
            by_id.insert(loaded.example.id.clone(), loaded.example);
        }
        let mut entries = vec![];
        for chapter in manifest.chapters {
            let example = by_id.get(&chapter.example_id).cloned().ok_or_else(|| {
                format!(
                    "Tutorial chapter '{}' references missing example id '{}'",
                    chapter.id, chapter.example_id
                )
            })?;
            let review_entry = review_by_id.get(&chapter.id);
            entries.push(TutorialProjectEntry {
                chapter_id: chapter.id,
                chapter_order: chapter.order,
                chapter_title: chapter.title,
                chapter_summary: chapter.summary,
                chapter_guide_path: chapter.guide_path,
                group_label: chapter.group_label,
                group_order: chapter.group_order,
                group_position: chapter.group_position,
                decimal_id: chapter.decimal_id,
                review_status: review_entry.and_then(|entry| entry.review_status.clone()),
                codex_reviewed_at: review_entry.and_then(|entry| entry.codex_reviewed_at.clone()),
                human_reviewed_at: review_entry.and_then(|entry| entry.human_reviewed_at.clone()),
                human_reviewer: review_entry.and_then(|entry| entry.human_reviewer.clone()),
                review_stale: review_entry
                    .map(|entry| entry.review_stale)
                    .unwrap_or(false),
                tier: chapter.tier,
                example,
                repo_root: repo_root.clone(),
            });
        }
        entries.sort_by(|left, right| {
            left.group_order
                .unwrap_or(usize::MAX)
                .cmp(&right.group_order.unwrap_or(usize::MAX))
                .then_with(|| {
                    left.group_position
                        .unwrap_or(usize::MAX)
                        .cmp(&right.group_position.unwrap_or(usize::MAX))
                })
                .then_with(|| left.chapter_order.cmp(&right.chapter_order))
        });
        Ok(entries)
    }

    fn tutorial_project_guided_walkthrough_entries() -> Vec<HelpTutorialDocEntry> {
        Self::discover_guided_walkthrough_entries()
    }

    fn tutorial_project_generated_chapter_path(entry: &TutorialProjectEntry) -> String {
        if let Some(decimal_id) = entry
            .decimal_id
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            return format!(
                "docs/tutorial/generated/chapters/{}_{}.md",
                decimal_id.replace('.', "-"),
                entry.chapter_id
            );
        }
        format!(
            "docs/tutorial/generated/chapters/{:02}_{}.md",
            entry.chapter_order, entry.chapter_id
        )
    }

    fn tutorial_project_progress_message(
        chapter_id: &str,
        chapter_title: &str,
        phase: &str,
        item: &str,
        percent: Option<f32>,
    ) -> TutorialProjectTaskProgress {
        TutorialProjectTaskProgress {
            chapter_id: chapter_id.to_string(),
            chapter_title: chapter_title.to_string(),
            phase: phase.to_string(),
            item: item.to_string(),
            percent,
        }
    }

    fn tutorial_project_scale_workflow_percent(percent: Option<f32>) -> Option<f32> {
        percent.map(|value| 0.15 + value.clamp(0.0, 1.0) * 0.75)
    }

    fn tutorial_project_progress_from_operation(
        chapter_id: &str,
        chapter_title: &str,
        progress: &OperationProgress,
    ) -> TutorialProjectTaskProgress {
        let mut message = match progress {
            OperationProgress::PrimerDesign(p) => {
                let percent = match (p.pair_evaluated, p.pair_evaluation_limit) {
                    (Some(done), Some(total)) if total > 0 => Some(done as f32 / total as f32),
                    _ if p.done => Some(1.0),
                    _ => None,
                };
                Self::tutorial_project_progress_message(
                    chapter_id,
                    chapter_title,
                    "execute_workflow",
                    &format!(
                        "Primer design ({}) via {}: {} ({})",
                        p.design_kind, p.backend_used, p.stage, p.detail
                    ),
                    percent,
                )
            }
            OperationProgress::Tfbs(p) => Self::tutorial_project_progress_message(
                chapter_id,
                chapter_title,
                "execute_workflow",
                &format!(
                    "TFBS motif {} ({}/{})",
                    p.motif_id, p.motif_index, p.motif_count
                ),
                Some((p.total_percent as f32 / 100.0).clamp(0.0, 1.0)),
            ),
            OperationProgress::GenomePrepare(p) => Self::tutorial_project_progress_message(
                chapter_id,
                chapter_title,
                "execute_workflow",
                &format!("Prepare genome: {} ({})", p.phase, p.item),
                p.percent
                    .map(|value| (value as f32 / 100.0).clamp(0.0, 1.0)),
            ),
            OperationProgress::GenomeTrackImport(p) => Self::tutorial_project_progress_message(
                chapter_id,
                chapter_title,
                "execute_workflow",
                &format!(
                    "Track import '{}' parsed={} imported={} skipped={}",
                    p.source, p.parsed_records, p.imported_features, p.skipped_records
                ),
                if p.done { Some(1.0) } else { None },
            ),
            OperationProgress::DbSnpFetch(p) => Self::tutorial_project_progress_message(
                chapter_id,
                chapter_title,
                "execute_workflow",
                &format!("dbSNP {}: {}", p.stage.as_str(), p.detail),
                None,
            ),
            OperationProgress::ReadAcquisition(p) => {
                let phase = p.phase.as_deref().unwrap_or(&p.lifecycle_status);
                let item = p.item.as_deref().unwrap_or(&p.display_name);
                Self::tutorial_project_progress_message(
                    chapter_id,
                    chapter_title,
                    "execute_workflow",
                    &format!("Read acquisition: {phase} ({item})"),
                    p.percent
                        .map(|value| (value as f32 / 100.0).clamp(0.0, 1.0)),
                )
            }
            OperationProgress::RnaReadInterpret(p) => {
                let percent = if p.reads_total > 0 {
                    Some((p.reads_processed as f32 / p.reads_total as f32).clamp(0.0, 1.0))
                } else {
                    None
                };
                Self::tutorial_project_progress_message(
                    chapter_id,
                    chapter_title,
                    "execute_workflow",
                    &format!(
                        "RNA-read mapping: reads {} / {}, aligned {}",
                        p.reads_processed, p.reads_total, p.aligned
                    ),
                    percent,
                )
            }
        };
        if message.phase == "execute_workflow" {
            message.percent = Self::tutorial_project_scale_workflow_percent(message.percent);
        }
        message
    }

    fn tutorial_project_cancelled_message(context: &str) -> String {
        format!("Tutorial project opening cancelled by caller ({context})")
    }

    fn tutorial_project_status_from_progress(progress: &TutorialProjectTaskProgress) -> String {
        format!(
            "Opening tutorial project '{}' ({}): {} ({})",
            progress.chapter_title, progress.chapter_id, progress.phase, progress.item
        )
    }

    fn open_tutorial_project_chapter(&mut self, chapter_id: &str) {
        if self.tutorial_project_task.is_some() {
            self.tutorial_project_status =
                "A tutorial project is already opening in the background".to_string();
            self.app_status = self.tutorial_project_status.clone();
            self.open_background_jobs_panel();
            return;
        }

        let job_id = self.alloc_background_job_id();
        let cancel_requested = Arc::new(AtomicBool::new(false));
        let (tx, rx) = mpsc::channel::<TutorialProjectTaskMessage>();
        let chapter_id_owned = chapter_id.trim().to_string();
        self.tutorial_project_progress_fraction = Some(0.0);
        self.tutorial_project_progress_label = "Queued tutorial build".to_string();
        self.tutorial_project_status = format!(
            "Opening tutorial project '{}' in background. You can keep using the UI.",
            chapter_id_owned
        );
        self.app_status = self.tutorial_project_status.clone();
        self.open_background_jobs_panel();
        self.push_job_event(
            BackgroundJobKind::OpenTutorialProject,
            BackgroundJobEventPhase::Started,
            Some(job_id),
            format!(
                "Opening tutorial chapter '{}' in background",
                chapter_id_owned
            ),
        );
        self.tutorial_project_task = Some(TutorialProjectTask {
            job_id,
            started: Instant::now(),
            cancel_requested: cancel_requested.clone(),
            chapter_id: chapter_id_owned.clone(),
            chapter_title: chapter_id_owned.clone(),
            receiver: rx,
        });

        std::thread::spawn(move || {
            let send_progress = |progress: TutorialProjectTaskProgress| {
                let _ = tx.send(TutorialProjectTaskMessage::Progress { job_id, progress });
            };
            send_progress(Self::tutorial_project_progress_message(
                &chapter_id_owned,
                &chapter_id_owned,
                "load_catalog",
                "Loading tutorial catalog and chapter metadata",
                Some(0.0),
            ));
            let entries = match Self::load_tutorial_project_entries() {
                Ok(entries) => entries,
                Err(err) => {
                    let _ = tx.send(TutorialProjectTaskMessage::Done {
                        job_id,
                        result: Err(format!("Tutorial catalog unavailable: {err}")),
                    });
                    return;
                }
            };
            let entry = match entries
                .into_iter()
                .find(|entry| entry.chapter_id == chapter_id_owned)
            {
                Some(entry) => entry,
                None => {
                    let _ = tx.send(TutorialProjectTaskMessage::Done {
                        job_id,
                        result: Err(format!(
                            "Unknown tutorial chapter id '{}'",
                            chapter_id_owned
                        )),
                    });
                    return;
                }
            };
            if cancel_requested.load(Ordering::Relaxed) {
                let _ = tx.send(TutorialProjectTaskMessage::Done {
                    job_id,
                    result: Err(Self::tutorial_project_cancelled_message(
                        "before validation",
                    )),
                });
                return;
            }
            send_progress(Self::tutorial_project_progress_message(
                &entry.chapter_id,
                &entry.chapter_title,
                "validate_inputs",
                "Checking required local files for the tutorial workflow",
                Some(0.05),
            ));
            if let Err(err) = validate_example_required_files(&entry.example, &entry.repo_root) {
                let _ = tx.send(TutorialProjectTaskMessage::Done {
                    job_id,
                    result: Err(format!(
                        "Tutorial '{}' missing required files: {err}",
                        entry.chapter_id
                    )),
                });
                return;
            }
            if cancel_requested.load(Ordering::Relaxed) {
                let _ = tx.send(TutorialProjectTaskMessage::Done {
                    job_id,
                    result: Err(Self::tutorial_project_cancelled_message("after validation")),
                });
                return;
            }
            let slug = format!(
                "{:02}_{}",
                entry.chapter_order,
                Self::sanitize_file_stem(&entry.chapter_id, "tutorial")
            );
            let cache_root = Self::tutorial_cache_root_dir();
            let run_dir = cache_root.join(&slug).join("run");
            send_progress(Self::tutorial_project_progress_message(
                &entry.chapter_id,
                &entry.chapter_title,
                "prepare_workspace",
                &format!("Preparing temporary run directory '{}'", run_dir.display()),
                Some(0.1),
            ));
            if run_dir.exists() {
                let _ = fs::remove_dir_all(&run_dir);
            }
            if let Err(err) = fs::create_dir_all(&run_dir) {
                let _ = tx.send(TutorialProjectTaskMessage::Done {
                    job_id,
                    result: Err(format!(
                        "Could not create tutorial run directory '{}': {err}",
                        run_dir.display()
                    )),
                });
                return;
            }
            if cancel_requested.load(Ordering::Relaxed) {
                let _ = tx.send(TutorialProjectTaskMessage::Done {
                    job_id,
                    result: Err(Self::tutorial_project_cancelled_message(
                        "before workflow execution",
                    )),
                });
                return;
            }
            send_progress(Self::tutorial_project_progress_message(
                &entry.chapter_id,
                &entry.chapter_title,
                "execute_workflow",
                "Running canonical tutorial workflow example",
                Some(0.15),
            ));
            let mut on_progress = |progress: OperationProgress| -> bool {
                if cancel_requested.load(Ordering::Relaxed) {
                    return false;
                }
                let progress = Self::tutorial_project_progress_from_operation(
                    &entry.chapter_id,
                    &entry.chapter_title,
                    &progress,
                );
                let _ = tx.send(TutorialProjectTaskMessage::Progress { job_id, progress });
                true
            };
            let state = match run_example_workflow_for_project_state_with_progress(
                &entry.example,
                &entry.repo_root,
                &run_dir,
                &mut on_progress,
            ) {
                Ok(state) => state,
                Err(err) => {
                    let _ = tx.send(TutorialProjectTaskMessage::Done {
                        job_id,
                        result: Err(format!(
                            "Could not build tutorial project '{}': {err}",
                            entry.chapter_id
                        )),
                    });
                    return;
                }
            };
            let project_path = cache_root.join(format!("{slug}.project.gentle.json"));
            send_progress(Self::tutorial_project_progress_message(
                &entry.chapter_id,
                &entry.chapter_title,
                "persist_project",
                &format!(
                    "Writing temporary tutorial project '{}'",
                    project_path.display()
                ),
                Some(0.95),
            ));
            if let Err(err) = state.save_to_path(project_path.to_string_lossy().as_ref()) {
                let _ = tx.send(TutorialProjectTaskMessage::Done {
                    job_id,
                    result: Err(format!(
                        "Could not persist tutorial project '{}': {err}",
                        project_path.display()
                    )),
                });
                return;
            }
            let guide_path = entry
                .chapter_guide_path
                .clone()
                .unwrap_or_else(|| Self::tutorial_project_generated_chapter_path(&entry));
            let guide_summary = if entry.chapter_guide_path.is_some() {
                format!("tutorial project guide: {}", entry.chapter_title)
            } else {
                format!("generated tutorial chapter: {}", entry.chapter_title)
            };
            send_progress(Self::tutorial_project_progress_message(
                &entry.chapter_id,
                &entry.chapter_title,
                "open_project",
                "Tutorial project built; opening the project and guide",
                Some(0.98),
            ));
            let _ = tx.send(TutorialProjectTaskMessage::Done {
                job_id,
                result: Ok(TutorialProjectOpenOutcome {
                    chapter_id: entry.chapter_id,
                    chapter_title: entry.chapter_title,
                    project_path: project_path.to_string_lossy().to_string(),
                    guide_path,
                    guide_summary,
                }),
            });
        });
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
        if let Some(paths) = rfd::FileDialog::new()
            .add_filter(
                "Supported sequence files",
                &[
                    "dna", "gb", "gbk", "gbff", "genbank", "embl", "emb", "fa", "fasta", "fna",
                    "fas", "xml",
                ],
            )
            .add_filter("SnapGene DNA", &["dna"])
            .add_filter("GenBank", &["gb", "gbk", "gbff", "genbank"])
            .add_filter("EMBL", &["embl", "emb"])
            .add_filter("FASTA", &["fa", "fasta", "fna", "fas"])
            .add_filter("NCBI Sequence XML (GBSet/GBSeq, INSDSet/INSDSeq)", &["xml"])
            .pick_files()
        {
            self.open_new_windows_from_files(&paths);
        }
    }

    fn prompt_save_project(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .set_file_name(self.default_project_save_file_name())
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
        let Some(path) = rfd::FileDialog::new()
            .set_file_name(self.default_lineage_svg_file_name())
            .add_filter("SVG", &["svg"])
            .save_file()
        else {
            self.app_status = "Lineage SVG export canceled".to_string();
            return;
        };
        let path = path.display().to_string();
        match self.export_visible_lineage_svg_to_path(&path) {
            Ok(status) => self.app_status = status,
            Err(err) => self.app_status = err,
        }
    }

    fn prompt_export_lab_assistant_report(&mut self) {
        let Some(mut path) = rfd::FileDialog::new()
            .set_file_name(self.default_lab_assistant_report_file_name())
            .add_filter("Lab assistant report", LAB_ASSISTANT_REPORT_EXTENSIONS)
            .add_filter("OpenDocument Text", &["odt"])
            .add_filter("Word document", &["docx"])
            .add_filter("Markdown", &["md"])
            .save_file()
        else {
            self.app_status = "Lab assistant report export canceled".to_string();
            return;
        };
        if path.extension().and_then(|value| value.to_str()).is_none() {
            path.set_extension(LabAssistantInstructionsFormat::Odt.extension());
        }
        let path_text = path.display().to_string();
        let result =
            self.engine
                .write()
                .unwrap()
                .apply(Operation::ExportLabAssistantInstructions {
                    path: path_text.clone(),
                    run_id: None,
                    title: None,
                    audience: Some("Lab assistant".to_string()),
                    format: None,
                });
        match result {
            Ok(result) => {
                self.app_status = Self::format_op_result_status(
                    &format!("Export lab assistant report to '{}': ok", path_text),
                    &result.created_seq_ids,
                    &result.warnings,
                    &result.messages,
                );
            }
            Err(err) => {
                self.app_status = format!(
                    "Export lab assistant report to '{}' failed: {}",
                    path_text, err
                );
            }
        }
    }

    fn prompt_export_dotplot_svg_from_lineage(
        &mut self,
        seq_id: &str,
        dotplot_id: &str,
        reference_seq_id: Option<&str>,
    ) {
        let default_stem = if let Some(reference_seq_id) = reference_seq_id
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            format!(
                "{}_vs_{}_{}",
                seq_id.trim(),
                reference_seq_id,
                dotplot_id.trim()
            )
        } else {
            format!("{}_{}", seq_id.trim(), dotplot_id.trim())
        };
        let stem = Self::sanitize_file_stem(&default_stem, "dotplot");
        let default_file_name = format!("{stem}.dotplot.svg");
        let path = rfd::FileDialog::new()
            .set_file_name(&default_file_name)
            .add_filter("SVG", &["svg"])
            .save_file();
        let Some(path) = path else {
            self.app_status = "Dotplot SVG export canceled".to_string();
            return;
        };
        let path_text = path.display().to_string();
        let result = self
            .engine
            .write()
            .unwrap()
            .apply(Operation::RenderDotplotSvg {
                seq_id: seq_id.to_string(),
                dotplot_id: dotplot_id.to_string(),
                path: path_text.clone(),
                flex_track_id: None,
                display_density_threshold: None,
                display_intensity_gain: None,
                overlay_x_axis_mode: Default::default(),
                overlay_anchor_exon: None,
            });
        match result {
            Ok(op_result) => {
                self.app_status = op_result
                    .messages
                    .first()
                    .cloned()
                    .unwrap_or_else(|| format!("Wrote dotplot SVG to '{path_text}'"));
            }
            Err(e) => {
                self.app_status = format!("Could not export dotplot SVG: {}", e.message);
            }
        }
    }

    fn prompt_export_filtered_retry_snapshots(&mut self) {
        if self.filtered_retry_snapshots().is_empty() {
            self.app_status = "No retry snapshots match current filters".to_string();
            return;
        }
        let default_file_name = format!(
            "retry_snapshots_{}.json",
            self.retry_snapshot_kind_filter.export_token()
        );
        let path = rfd::FileDialog::new()
            .set_file_name(&default_file_name)
            .add_filter("JSON", &["json"])
            .save_file();
        let Some(path) = path else {
            self.app_status = "Retry snapshot export canceled".to_string();
            return;
        };
        let path_text = path.display().to_string();
        match self.export_filtered_retry_snapshots_to_path(path.as_path()) {
            Ok(count) => {
                self.app_status = format!("Exported {count} retry snapshot(s) to '{path_text}'");
            }
            Err(err) => {
                self.app_status =
                    format!("Could not export retry snapshots '{}': {err}", path_text);
            }
        }
    }

    fn prompt_export_retry_cleanup_audit_report(&mut self) {
        if self.filtered_retry_cleanup_audit_entries().is_empty() {
            self.app_status = "No retry cleanup audit entries match current filters".to_string();
            return;
        }
        let path = rfd::FileDialog::new()
            .set_file_name("retry_cleanup_audit_report_filtered.json")
            .add_filter("JSON", &["json"])
            .save_file();
        let Some(path) = path else {
            self.app_status = "Retry cleanup audit export canceled".to_string();
            return;
        };
        let path_text = path.display().to_string();
        match self.export_retry_cleanup_audit_report_to_path(path.as_path()) {
            Ok(count) => {
                // Exporting the audit report is read-only and intentionally not self-audited.
                self.app_status = format!(
                    "Exported {count} filtered retry cleanup audit entries to '{path_text}'"
                );
            }
            Err(err) => {
                self.app_status = format!(
                    "Could not export retry cleanup audit report '{}': {err}",
                    path_text
                );
            }
        }
    }

    fn prompt_archive_and_delete_filtered_retry_snapshots(&mut self) -> Option<(usize, String)> {
        if self.filtered_retry_snapshots().is_empty() {
            self.app_status = "No retry snapshots match current filters".to_string();
            return None;
        }
        let default_file_name = format!(
            "retry_snapshots_archive_{}.json",
            self.retry_snapshot_kind_filter.export_token()
        );
        let path = rfd::FileDialog::new()
            .set_file_name(&default_file_name)
            .add_filter("JSON", &["json"])
            .save_file();
        let Some(path) = path else {
            self.app_status = "Retry snapshot archive canceled".to_string();
            return None;
        };
        let path_text = path.display().to_string();
        match self.archive_and_delete_filtered_retry_snapshots_to_path(path.as_path()) {
            Ok(count) => {
                self.app_status =
                    format!("Archived and removed {count} retry snapshot(s) to '{path_text}'");
                Some((count, path_text))
            }
            Err(err) => {
                self.app_status = format!(
                    "Could not archive/delete retry snapshots '{}': {err}",
                    path_text
                );
                None
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

    fn set_scope_genome_paths(&mut self, scope: GenomeDialogScope, catalog: String, cache: String) {
        match scope {
            GenomeDialogScope::Reference => {
                self.reference_genome_catalog_path = catalog;
                self.reference_genome_cache_dir = cache;
            }
            GenomeDialogScope::Helper => {
                self.helper_genome_catalog_path = catalog;
                self.helper_genome_cache_dir = cache;
            }
        }
    }

    fn scope_genome_paths_resolved(&self, scope: GenomeDialogScope) -> (String, String) {
        let (catalog, cache) = match scope {
            GenomeDialogScope::Reference => (
                self.reference_genome_catalog_path.trim(),
                self.reference_genome_cache_dir.trim(),
            ),
            GenomeDialogScope::Helper => (
                self.helper_genome_catalog_path.trim(),
                self.helper_genome_cache_dir.trim(),
            ),
        };
        let resolved_catalog = if catalog.is_empty() {
            scope.default_catalog_path().to_string()
        } else {
            catalog.to_string()
        };
        let resolved_cache = if cache.is_empty() {
            scope.default_cache_dir().to_string()
        } else {
            cache.to_string()
        };
        (resolved_catalog, resolved_cache)
    }

    fn sync_active_genome_scope_paths_from_fields(&mut self) {
        self.set_scope_genome_paths(
            self.genome_dialog_scope,
            self.genome_catalog_path.clone(),
            self.genome_cache_dir.clone(),
        );
    }

    fn set_genome_dialog_scope(&mut self, scope: GenomeDialogScope) {
        self.sync_active_genome_scope_paths_from_fields();
        self.genome_dialog_scope = scope;
        let (next_catalog, next_cache) = self.scope_genome_paths_resolved(scope);
        let catalog_changed = self.genome_catalog_path != next_catalog;
        let cache_changed = self.genome_cache_dir != next_cache;
        self.genome_catalog_path = next_catalog;
        self.genome_cache_dir = next_cache;
        if catalog_changed || cache_changed {
            self.invalidate_genome_genes();
        }
    }

    fn queue_prepared_genome_reinstall(
        &mut self,
        genome_id: impl Into<String>,
        scope: GenomeDialogScope,
        dialog_host: PreparedGenomeReinstallDialogHost,
    ) {
        self.pending_prepared_genome_reinstall = Some(PreparedGenomeReinstallRequest {
            genome_id: genome_id.into(),
            scope,
            catalog_path: self.genome_catalog_path.clone(),
            cache_dir: self.genome_cache_dir.clone(),
            dialog_host,
        });
    }

    fn pending_prepared_genome_reinstall_targets_host(
        &self,
        dialog_host: PreparedGenomeReinstallDialogHost,
    ) -> bool {
        self.pending_prepared_genome_reinstall
            .as_ref()
            .is_some_and(|request| request.dialog_host == dialog_host)
    }

    fn dismiss_pending_prepared_genome_reinstall_for_host(
        &mut self,
        dialog_host: PreparedGenomeReinstallDialogHost,
    ) {
        if self.pending_prepared_genome_reinstall_targets_host(dialog_host) {
            self.pending_prepared_genome_reinstall = None;
        }
    }

    fn apply_prepared_genome_reinstall_request(&mut self, request: PreparedGenomeReinstallRequest) {
        self.sync_active_genome_scope_paths_from_fields();
        let scope_changed = self.genome_dialog_scope != request.scope;
        let catalog_changed = self.genome_catalog_path != request.catalog_path;
        let cache_changed = self.genome_cache_dir != request.cache_dir;
        let genome_changed = self.genome_id != request.genome_id;
        self.genome_dialog_scope = request.scope;
        self.genome_catalog_path = request.catalog_path.clone();
        self.genome_cache_dir = request.cache_dir.clone();
        self.set_scope_genome_paths(request.scope, request.catalog_path, request.cache_dir);
        self.genome_id = request.genome_id;
        if scope_changed || catalog_changed || cache_changed || genome_changed {
            self.invalidate_genome_genes();
        }
    }

    fn confirm_prepared_genome_reinstall(&mut self, mode: GenomePrepareLaunchMode) {
        if let Some(request) = self.pending_prepared_genome_reinstall.take() {
            self.apply_prepared_genome_reinstall_request(request);
            self.start_prepare_reference_genome_with_mode(mode);
        }
    }

    fn open_prepare_genome_dialog_for_scope(&mut self, scope: GenomeDialogScope) {
        self.set_genome_dialog_scope(scope);
        let was_open = self.show_reference_genome_prepare_dialog;
        self.show_reference_genome_prepare_dialog = true;
        self.mark_window_open_or_focus(Self::prepare_genome_viewport_id(), was_open);
    }

    fn open_retrieve_genome_dialog_for_scope(&mut self, scope: GenomeDialogScope) {
        self.set_genome_dialog_scope(scope);
        let was_open = self.show_reference_genome_retrieve_dialog;
        self.show_reference_genome_retrieve_dialog = true;
        self.mark_window_open_or_focus(Self::retrieve_genome_viewport_id(), was_open);
    }

    fn open_blast_genome_dialog_for_scope(&mut self, scope: GenomeDialogScope) {
        self.set_genome_dialog_scope(scope);
        let was_open = self.show_reference_genome_blast_dialog;
        self.show_reference_genome_blast_dialog = true;
        self.mark_window_open_or_focus(Self::blast_genome_viewport_id(), was_open);
    }

    fn open_reference_genome_prepare_dialog(&mut self) {
        self.open_prepare_genome_dialog_for_scope(GenomeDialogScope::Reference);
    }

    fn open_reference_genome_retrieve_dialog(&mut self) {
        self.open_retrieve_genome_dialog_for_scope(GenomeDialogScope::Reference);
    }

    fn open_reference_genome_blast_dialog(&mut self) {
        self.open_blast_genome_dialog_for_scope(GenomeDialogScope::Reference);
    }

    fn open_reference_genome_inspector_dialog(&mut self) {
        self.show_reference_genome_inspector_dialog = true;
    }

    fn open_cache_cleanup_dialog(&mut self) {
        self.sync_active_genome_scope_paths_from_fields();
        self.cache_cleanup_scope = match self.genome_dialog_scope {
            GenomeDialogScope::Reference => CacheCleanupScope::References,
            GenomeDialogScope::Helper => CacheCleanupScope::Helpers,
        };
        self.cache_cleanup_reference_cache_dir =
            if self.reference_genome_cache_dir.trim().is_empty() {
                configured_reference_genome_cache_dir()
            } else {
                self.reference_genome_cache_dir.trim().to_string()
            };
        self.cache_cleanup_helper_cache_dir = if self.helper_genome_cache_dir.trim().is_empty() {
            configured_helper_genome_cache_dir()
        } else {
            self.helper_genome_cache_dir.trim().to_string()
        };
        self.cache_cleanup_mode = PreparedCacheCleanupMode::DerivedIndexesOnly;
        self.cache_cleanup_include_orphans = false;
        self.cache_cleanup_confirm_pending = false;
        self.cache_cleanup_selected_paths.clear();
        self.cache_cleanup_rebuild_candidates.clear();
        self.show_cache_cleanup_dialog = true;
        self.refresh_cache_cleanup_inspection();
    }

    fn cache_cleanup_selected_roots(&self) -> Vec<String> {
        self.cache_cleanup_scope.roots(
            &self.cache_cleanup_reference_cache_dir,
            &self.cache_cleanup_helper_cache_dir,
        )
    }

    fn refresh_cache_cleanup_inspection(&mut self) {
        self.cache_cleanup_confirm_pending = false;
        self.cache_cleanup_rebuild_candidates.clear();
        let roots = self.cache_cleanup_selected_roots();
        match GentleEngine::inspect_prepared_cache_roots(&roots) {
            Ok(report) => {
                let valid_paths = report
                    .entries
                    .iter()
                    .map(|entry| entry.path.clone())
                    .collect::<HashSet<_>>();
                self.cache_cleanup_selected_paths
                    .retain(|path| valid_paths.contains(path));
                self.cache_cleanup_inspection = Some(report);
                self.cache_cleanup_status = if self
                    .cache_cleanup_inspection
                    .as_ref()
                    .is_some_and(|report| report.entry_count == 0)
                {
                    "Nothing to clean in selected cache roots".to_string()
                } else {
                    format!("Inspected {} selected cache root(s).", roots.len())
                };
            }
            Err(e) => {
                self.cache_cleanup_inspection = None;
                self.cache_cleanup_status = format!("Could not inspect prepared caches: {e}");
            }
        }
    }

    fn cache_cleanup_group_label(group: PreparedCacheArtifactGroup) -> &'static str {
        match group {
            PreparedCacheArtifactGroup::CachedSources => "sources",
            PreparedCacheArtifactGroup::DerivedIndexes => "indexes",
            PreparedCacheArtifactGroup::BlastDb => "blast",
        }
    }

    fn format_cache_cleanup_groups(entry: &PreparedCacheInspectionEntry) -> String {
        let mut labels = entry
            .artifact_stats
            .iter()
            .map(|stat| Self::cache_cleanup_group_label(stat.group).to_string())
            .collect::<Vec<_>>();
        labels.sort();
        labels.join(", ")
    }

    fn cache_cleanup_group_bytes(
        entry: &PreparedCacheInspectionEntry,
        group: PreparedCacheArtifactGroup,
    ) -> u64 {
        entry
            .artifact_stats
            .iter()
            .find(|stat| stat.group == group)
            .map(|stat| stat.total_size_bytes)
            .unwrap_or(0)
    }

    fn cache_cleanup_preview_targets(
        &self,
    ) -> (Vec<&PreparedCacheInspectionEntry>, u64, usize, String) {
        let Some(report) = &self.cache_cleanup_inspection else {
            return (
                vec![],
                0,
                0,
                "Inspect selected cache roots first.".to_string(),
            );
        };
        let entries = report
            .entries
            .iter()
            .filter(|entry| match self.cache_cleanup_mode {
                PreparedCacheCleanupMode::AllPreparedInCache => {
                    entry.classification == crate::genomes::PreparedCacheEntryKind::PreparedInstall
                        || (self.cache_cleanup_include_orphans
                            && entry.classification
                                == crate::genomes::PreparedCacheEntryKind::OrphanedRemnant)
                }
                PreparedCacheCleanupMode::SelectedPreparedInstalls => {
                    self.cache_cleanup_selected_paths.contains(&entry.path)
                }
                PreparedCacheCleanupMode::BlastDbOnly
                | PreparedCacheCleanupMode::DerivedIndexesOnly => {
                    entry.classification == crate::genomes::PreparedCacheEntryKind::PreparedInstall
                        && self.cache_cleanup_selected_paths.contains(&entry.path)
                }
            })
            .collect::<Vec<_>>();
        let bytes = entries
            .iter()
            .map(|entry| match self.cache_cleanup_mode {
                PreparedCacheCleanupMode::BlastDbOnly => {
                    Self::cache_cleanup_group_bytes(entry, PreparedCacheArtifactGroup::BlastDb)
                }
                PreparedCacheCleanupMode::DerivedIndexesOnly => {
                    Self::cache_cleanup_group_bytes(
                        entry,
                        PreparedCacheArtifactGroup::DerivedIndexes,
                    ) + Self::cache_cleanup_group_bytes(entry, PreparedCacheArtifactGroup::BlastDb)
                }
                PreparedCacheCleanupMode::SelectedPreparedInstalls
                | PreparedCacheCleanupMode::AllPreparedInCache => entry.total_size_bytes,
            })
            .sum::<u64>();
        let explanation = match self.cache_cleanup_mode {
            PreparedCacheCleanupMode::BlastDbOnly => {
                "This removes BLAST databases only; cached FASTA/annotation sources and rebuildable indexes remain."
            }
            PreparedCacheCleanupMode::DerivedIndexesOnly => {
                "This removes rebuildable indexes (FASTA index, gene index, and BLAST DB) while keeping cached FASTA/annotation sources and manifests."
            }
            PreparedCacheCleanupMode::SelectedPreparedInstalls => {
                "This removes the selected prepared install directories entirely. Cached sources inside those installs will be deleted too."
            }
            PreparedCacheCleanupMode::AllPreparedInCache => {
                "This removes all prepared installs under the selected cache roots. Orphaned remnants are removed only when that checkbox is enabled."
            }
        }
        .to_string();
        (entries, bytes, report.cache_roots.len(), explanation)
    }

    fn cache_cleanup_root_matches(selected_root: &str, entry_root: &str) -> bool {
        crate::genomes::paths_refer_to_same_location(
            Path::new(selected_root.trim()),
            Path::new(entry_root.trim()),
        )
    }

    fn cache_cleanup_scope_for_root(&self, cache_root: &str) -> Option<GenomeDialogScope> {
        let reference_root = if self.cache_cleanup_reference_cache_dir.trim().is_empty() {
            configured_reference_genome_cache_dir()
        } else {
            self.cache_cleanup_reference_cache_dir.trim().to_string()
        };
        let helper_root = if self.cache_cleanup_helper_cache_dir.trim().is_empty() {
            configured_helper_genome_cache_dir()
        } else {
            self.cache_cleanup_helper_cache_dir.trim().to_string()
        };
        let reference_matches = !matches!(self.cache_cleanup_scope, CacheCleanupScope::Helpers)
            && Self::cache_cleanup_root_matches(&reference_root, cache_root);
        let helper_matches = !matches!(self.cache_cleanup_scope, CacheCleanupScope::References)
            && Self::cache_cleanup_root_matches(&helper_root, cache_root);
        match (reference_matches, helper_matches) {
            (true, false) => Some(GenomeDialogScope::Reference),
            (false, true) => Some(GenomeDialogScope::Helper),
            _ => None,
        }
    }

    fn cache_cleanup_rebuild_candidates_from_report(
        &self,
        report: &crate::genomes::PreparedCacheCleanupReport,
    ) -> Vec<PreparedGenomeReinstallRequest> {
        if !matches!(
            report.mode,
            PreparedCacheCleanupMode::BlastDbOnly | PreparedCacheCleanupMode::DerivedIndexesOnly
        ) {
            return vec![];
        }
        let mut requests = report
            .results
            .iter()
            .filter(|result| {
                result.removed
                    && matches!(
                        result.classification,
                        crate::genomes::PreparedCacheEntryKind::PreparedInstall
                    )
            })
            .filter_map(|result| {
                let scope = self.cache_cleanup_scope_for_root(&result.cache_root)?;
                let (catalog_path, _) = self.scope_genome_paths_resolved(scope);
                Some(PreparedGenomeReinstallRequest {
                    genome_id: result.entry_id.clone(),
                    scope,
                    catalog_path,
                    cache_dir: result.cache_root.clone(),
                    dialog_host: PreparedGenomeReinstallDialogHost::Root,
                })
            })
            .collect::<Vec<_>>();
        Self::sort_cache_cleanup_rebuild_candidates(&mut requests);
        requests
    }

    fn sort_cache_cleanup_rebuild_candidates(requests: &mut Vec<PreparedGenomeReinstallRequest>) {
        requests.sort_by(|left, right| {
            (
                left.scope.description(),
                left.genome_id.as_str(),
                left.cache_dir.as_str(),
            )
                .cmp(&(
                    right.scope.description(),
                    right.genome_id.as_str(),
                    right.cache_dir.as_str(),
                ))
        });
        requests.dedup_by(|left, right| {
            left.scope == right.scope
                && left.genome_id == right.genome_id
                && left.cache_dir == right.cache_dir
        });
    }

    fn queue_cache_cleanup_rebuild_candidate(&mut self, index: usize) {
        if let Some(request) = self.cache_cleanup_rebuild_candidates.get(index).cloned() {
            self.pending_prepared_genome_reinstall = Some(request);
        }
    }

    fn apply_cache_cleanup(&mut self) {
        let roots = self.cache_cleanup_selected_roots();
        let request = PreparedCacheCleanupRequest {
            mode: self.cache_cleanup_mode,
            cache_roots: roots.clone(),
            prepared_ids: vec![],
            prepared_paths: self.cache_cleanup_selected_paths.iter().cloned().collect(),
            include_orphaned_remnants: self.cache_cleanup_include_orphans,
        };
        match GentleEngine::clear_prepared_cache_roots(&request) {
            Ok(report) => {
                let rebuild_candidates = self.cache_cleanup_rebuild_candidates_from_report(&report);
                let status = format!(
                    "Cache cleanup removed {} item(s), reclaiming {}.",
                    report.removed_item_count,
                    Self::format_bytes_compact(report.removed_bytes)
                );
                self.cache_cleanup_confirm_pending = false;
                self.cache_cleanup_selected_paths.clear();
                self.refresh_cache_cleanup_inspection();
                self.cache_cleanup_rebuild_candidates = rebuild_candidates;
                self.cache_cleanup_status = status;
            }
            Err(e) => {
                self.cache_cleanup_status = format!("Could not clear caches: {e}");
            }
        }
        self.cache_cleanup_confirm_pending = false;
    }

    fn open_genome_bed_track_dialog(&mut self) {
        self.load_bed_track_subscriptions_from_state();
        let was_open = self.show_genome_bed_track_dialog;
        self.show_genome_bed_track_dialog = true;
        self.mark_window_open_or_focus(Self::bed_track_viewport_id(), was_open);
    }

    fn active_dna_window_context(&self) -> Option<(String, Option<(usize, usize)>)> {
        let active_key = self.active_window_menu_key?;
        let viewport_id = *self.native_window_key_to_viewport.get(&active_key)?;
        let window = self.windows.get(&viewport_id)?;
        let guard = window.read().ok()?;
        Some((guard.sequence_id()?, guard.selection_range_0based()))
    }

    fn gibson_ui_insert_rows(&self) -> Vec<GibsonUiInsertRow> {
        let mut rows = vec![];
        if !self.gibson_insert_seq_id.trim().is_empty() {
            rows.push(GibsonUiInsertRow {
                seq_id: self.gibson_insert_seq_id.trim().to_string(),
                orientation: self.gibson_insert_orientation,
            });
        }
        rows.extend(
            self.gibson_extra_inserts
                .iter()
                .filter(|row| !row.seq_id.trim().is_empty())
                .cloned(),
        );
        rows
    }

    fn gibson_is_multi_insert_draft(&self) -> bool {
        self.gibson_ui_insert_rows().len() > 1
    }

    fn gibson_multi_insert_defined_opening_note() -> &'static str {
        "Multi-insert Gibson currently requires a defined destination opening; 'existing_termini' remains the single-fragment handoff path."
    }

    fn gibson_apply_blocked_status(preview: &GibsonAssemblyPreview) -> String {
        if preview
            .errors
            .iter()
            .any(|err| err == Self::gibson_multi_insert_defined_opening_note())
        {
            format!(
                "Gibson apply blocked: {}",
                Self::gibson_multi_insert_defined_opening_note()
            )
        } else if let Some(first_error) = preview.errors.first() {
            format!(
                "Gibson apply blocked: {} ({} blocking error(s) total). Resolve them first and rerun preview.",
                first_error,
                preview.errors.len()
            )
        } else {
            format!(
                "Gibson apply blocked: preview still has {} blocking error(s). Resolve them first and rerun preview.",
                preview.errors.len()
            )
        }
    }

    fn refresh_gibson_output_id_hint_default(&mut self) {
        if !self.gibson_output_id_hint.trim().is_empty() {
            return;
        }
        let insert_ids = self
            .gibson_ui_insert_rows()
            .into_iter()
            .map(|row| row.seq_id)
            .collect::<Vec<_>>();
        if self.gibson_destination_seq_id.trim().is_empty() || insert_ids.is_empty() {
            return;
        }
        self.gibson_output_id_hint = format!(
            "{}_with_{}",
            self.gibson_destination_seq_id,
            insert_ids.join("_")
        );
    }

    fn prefill_gibson_from_active_context(&mut self) {
        let active_context = self.active_dna_window_context();
        if self.gibson_destination_seq_id.trim().is_empty() {
            if let Some((seq_id, _)) = active_context.as_ref() {
                self.gibson_destination_seq_id = seq_id.clone();
            } else if let Some(seq_id) = self.project_sequence_ids_for_blast().first() {
                self.gibson_destination_seq_id = seq_id.clone();
            }
        }
        if self.gibson_insert_seq_id.trim().is_empty() {
            self.gibson_insert_seq_id = self
                .project_sequence_ids_for_blast()
                .into_iter()
                .find(|seq_id| seq_id != &self.gibson_destination_seq_id)
                .unwrap_or_default();
        }
        if let Some((active_seq_id, Some((start, end)))) = active_context
            && active_seq_id == self.gibson_destination_seq_id
            && self.gibson_opening_start_0based.trim().is_empty()
            && self.gibson_opening_end_0based_exclusive.trim().is_empty()
        {
            self.gibson_opening_mode = GibsonUiOpeningMode::DefinedSite;
            self.gibson_opening_start_0based = start.to_string();
            self.gibson_opening_end_0based_exclusive = end.to_string();
        }
        if self.gibson_output_id_hint.trim().is_empty() {
            self.refresh_gibson_output_id_hint_default();
        }
    }

    fn open_pcr_design_dialog_for_seq_id(
        &mut self,
        seq_id: &str,
    ) -> std::result::Result<(), String> {
        let seq_id = seq_id.trim();
        if seq_id.is_empty() {
            return Err("sequence id is empty".to_string());
        }
        let exists = self
            .engine
            .read()
            .map_err(|_| "could not read engine state".to_string())?
            .state()
            .sequences
            .contains_key(seq_id);
        if !exists {
            return Err(format!(
                "sequence '{seq_id}' was not found in the current project"
            ));
        }
        if self.find_open_sequence_viewport_id(seq_id).is_none()
            && self.find_pending_sequence_window_mut(seq_id).is_none()
        {
            self.open_sequence_window(seq_id);
        }
        self.pcr_design_seq_id = seq_id.to_string();
        let was_open = self.show_pcr_design_dialog;
        self.show_pcr_design_dialog = true;
        self.mark_window_open_or_focus(Self::pcr_design_viewport_id(), was_open);
        Ok(())
    }

    fn open_pcr_design_dialog(&mut self) {
        if self.show_pcr_design_dialog && !self.pcr_design_seq_id.trim().is_empty() {
            self.mark_window_open_or_focus(Self::pcr_design_viewport_id(), true);
            return;
        }
        let target_seq_id = self
            .active_dna_window_context()
            .map(|(seq_id, _)| seq_id)
            .or_else(|| self.project_sequence_ids_for_blast().first().cloned());
        let Some(seq_id) = target_seq_id else {
            self.app_status =
                "Cannot open PCR Designer: no active sequence window or project sequence"
                    .to_string();
            return;
        };
        if let Err(err) = self.open_pcr_design_dialog_for_seq_id(&seq_id) {
            self.app_status = format!("Cannot open PCR Designer for '{seq_id}': {err}");
        }
    }

    fn open_sequencing_confirmation_dialog(&mut self) {
        if self.show_sequencing_confirmation_dialog {
            self.mark_window_open_or_focus(Self::sequencing_confirmation_viewport_id(), true);
            return;
        }
        let target_seq_id = self
            .active_dna_window_context()
            .map(|(seq_id, _)| seq_id)
            .or_else(|| self.project_sequence_ids_for_blast().first().cloned());
        let Some(seq_id) = target_seq_id else {
            self.app_status = "Cannot open Sequencing Confirmation: no active sequence window or project sequence".to_string();
            return;
        };
        if self.find_open_sequence_viewport_id(&seq_id).is_none() {
            self.open_sequence_window(&seq_id);
        }
        self.sequencing_confirmation_seq_id = seq_id;
        let was_open = self.show_sequencing_confirmation_dialog;
        self.show_sequencing_confirmation_dialog = true;
        self.mark_window_open_or_focus(Self::sequencing_confirmation_viewport_id(), was_open);
    }

    fn open_gibson_dialog(&mut self) {
        if self.show_gibson_dialog {
            self.mark_window_open_or_focus(Self::gibson_viewport_id(), true);
            return;
        }
        self.prefill_gibson_from_active_context();
        self.show_gibson_dialog = true;
        self.mark_window_open_or_focus(Self::gibson_viewport_id(), false);
    }

    fn open_planning_dialog(&mut self) {
        if self.show_planning_dialog {
            self.mark_window_open_or_focus(Self::planning_viewport_id(), true);
            return;
        }
        self.refresh_planning_editor_buffers_from_engine();
        self.show_planning_dialog = true;
        self.mark_window_open_or_focus(Self::planning_viewport_id(), false);
    }

    fn open_jaspar_expert_dialog(&mut self) {
        let was_open = self.show_jaspar_expert_dialog;
        if self.show_jaspar_expert_dialog {
            self.mark_window_open_or_focus(Self::jaspar_expert_viewport_id(), true);
            return;
        }
        if self.jaspar_expert_selected_motif_id.trim().is_empty() {
            self.jaspar_expert_selected_motif_id = tf_motifs::list_motif_summaries()
                .first()
                .map(|row| row.id.clone())
                .unwrap_or_default();
        }
        self.show_jaspar_expert_dialog = true;
        self.mark_window_open_or_focus(Self::jaspar_expert_viewport_id(), was_open);
    }

    fn open_jaspar_expert_dialog_for_motif(&mut self, motif_id: &str) {
        let trimmed = motif_id.trim();
        if !trimmed.is_empty() {
            self.jaspar_expert_selected_motif_id = trimmed.to_string();
            self.jaspar_expert_view = None;
        }
        self.open_jaspar_expert_dialog();
        if !trimmed.is_empty() {
            self.refresh_jaspar_expert_view();
        }
    }

    fn open_uniprot_dialog(&mut self) {
        let was_open = self.show_uniprot_dialog;
        if self.uniprot_map_seq_id.trim().is_empty()
            && let Some(seq_id) = self.project_sequence_ids_for_blast().first()
        {
            self.uniprot_map_seq_id = seq_id.clone();
        }
        self.show_uniprot_dialog = true;
        self.mark_window_open_or_focus(Self::uniprot_viewport_id(), was_open);
    }

    fn open_new_sequence_dialog(&mut self) {
        let was_open = self.show_new_sequence_dialog;
        self.show_new_sequence_dialog = true;
        self.mark_window_open_or_focus(Self::new_sequence_viewport_id(), was_open);
    }

    fn open_new_sequence_from_clipboard(&mut self) {
        match Self::read_text_clipboard() {
            Ok(text) => {
                self.new_sequence_text = text;
                let pasted_chars = self.new_sequence_text.chars().count();
                self.new_sequence_status =
                    format!("Loaded {pasted_chars} clipboard character(s); review before create.");
            }
            Err(err) => {
                self.new_sequence_status = err;
            }
        }
        self.open_new_sequence_dialog();
    }

    fn open_genbank_dialog(&mut self) {
        if self.show_genbank_dialog {
            self.mark_window_open_or_focus(Self::genbank_viewport_id(), true);
            return;
        }
        self.seed_dbsnp_defaults_if_needed();
        self.show_genbank_dialog = true;
        self.mark_window_open_or_focus(Self::genbank_viewport_id(), false);
    }

    fn infer_genome_dialog_scope_from_paths(
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
    ) -> GenomeDialogScope {
        let catalog = catalog_path
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .unwrap_or_default();
        let cache = cache_dir
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .unwrap_or_default();
        if catalog.eq_ignore_ascii_case(DEFAULT_HELPER_GENOME_CATALOG_PATH)
            || cache.eq_ignore_ascii_case(DEFAULT_HELPER_GENOME_CACHE_DIR)
        {
            GenomeDialogScope::Helper
        } else {
            GenomeDialogScope::Reference
        }
    }

    fn lineage_retrieval_descriptor_from_operation(
        op: &Operation,
    ) -> Option<LineageRetrievalDescriptor> {
        match op {
            Operation::ExtractGenomeGene {
                genome_id,
                gene_query,
                occurrence,
                catalog_path,
                cache_dir,
                ..
            } => Some(LineageRetrievalDescriptor::GenomeGene {
                scope: Self::infer_genome_dialog_scope_from_paths(
                    catalog_path.as_deref(),
                    cache_dir.as_deref(),
                ),
                genome_id: genome_id.clone(),
                gene_query: gene_query.clone(),
                occurrence: *occurrence,
                catalog_path: catalog_path.clone(),
                cache_dir: cache_dir.clone(),
            }),
            Operation::ExtractGenomeRegion {
                genome_id,
                chromosome,
                start_1based,
                end_1based,
                catalog_path,
                cache_dir,
                ..
            } => Some(LineageRetrievalDescriptor::GenomeRegion {
                scope: Self::infer_genome_dialog_scope_from_paths(
                    catalog_path.as_deref(),
                    cache_dir.as_deref(),
                ),
                genome_id: genome_id.clone(),
                chromosome: chromosome.clone(),
                start_1based: *start_1based,
                end_1based: *end_1based,
                catalog_path: catalog_path.clone(),
                cache_dir: cache_dir.clone(),
            }),
            Operation::FetchGenBankAccession { accession, as_id } => {
                Some(LineageRetrievalDescriptor::GenBankAccession {
                    accession: accession.clone(),
                    as_id: as_id.clone(),
                })
            }
            Operation::FetchDbSnpRegion {
                rs_id,
                genome_id,
                flank_bp,
                catalog_path,
                cache_dir,
                ..
            } => Some(LineageRetrievalDescriptor::DbSnpRegion {
                genome_id: genome_id.clone(),
                rs_id: rs_id.clone(),
                flank_bp: flank_bp.unwrap_or(3000),
                catalog_path: catalog_path.clone(),
                cache_dir: cache_dir.clone(),
            }),
            _ => None,
        }
    }

    fn lineage_retrieval_pattern_label(descriptor: &LineageRetrievalDescriptor) -> &'static str {
        match descriptor {
            LineageRetrievalDescriptor::GenomeGene { .. } => "GENE",
            LineageRetrievalDescriptor::GenomeRegion { .. } => "REGION",
            LineageRetrievalDescriptor::GenBankAccession { .. } => "GB",
            LineageRetrievalDescriptor::DbSnpRegion { .. } => "RS",
        }
    }

    fn lineage_retrieval_pattern_tooltip(
        descriptor: &LineageRetrievalDescriptor,
        anchor: Option<&SequenceGenomeAnchorSummary>,
    ) -> String {
        let mut lines: Vec<String> = Vec::with_capacity(4);
        match descriptor {
            LineageRetrievalDescriptor::GenomeGene {
                scope,
                genome_id,
                gene_query,
                occurrence,
                ..
            } => {
                let route = match scope {
                    GenomeDialogScope::Reference => "genome catalog",
                    GenomeDialogScope::Helper => "helper catalog",
                };
                lines.push(format!(
                    "Retrieval pattern: gene query '{}' on '{}' ({route})",
                    gene_query, genome_id
                ));
                if let Some(idx) = occurrence {
                    lines.push(format!("Occurrence: {}", idx));
                }
            }
            LineageRetrievalDescriptor::GenomeRegion {
                scope,
                genome_id,
                chromosome,
                start_1based,
                end_1based,
                ..
            } => {
                let route = match scope {
                    GenomeDialogScope::Reference => "genome catalog",
                    GenomeDialogScope::Helper => "helper catalog",
                };
                lines.push(format!(
                    "Retrieval pattern: region {}:{}-{} on '{}' ({route})",
                    chromosome, start_1based, end_1based, genome_id
                ));
            }
            LineageRetrievalDescriptor::GenBankAccession { accession, .. } => {
                lines.push(format!(
                    "Retrieval pattern: GenBank accession '{}'",
                    accession
                ));
                lines.push(
                    "GenBank retrieval and genome anchoring are complementary (no source conflict)"
                        .to_string(),
                );
            }
            LineageRetrievalDescriptor::DbSnpRegion {
                genome_id,
                rs_id,
                flank_bp,
                ..
            } => {
                lines.push(format!(
                    "Retrieval pattern: dbSNP '{}' +/-{} bp on '{}'",
                    rs_id, flank_bp, genome_id
                ));
                lines.push(
                    "dbSNP resolution uses NCBI Variation, then extracts the annotated locus from the selected prepared genome."
                        .to_string(),
                );
            }
        }
        if let Some(anchor) = anchor {
            let strand = anchor.strand.unwrap_or('+');
            lines.push(format!(
                "Genome anchor: {}:{}-{} (strand {})",
                anchor.chromosome, anchor.start_1based, anchor.end_1based, strand
            ));
        }
        lines.push("Click to reopen the matching retrieval dialog".to_string());
        lines.join("\n")
    }

    fn open_lineage_retrieval_descriptor(&mut self, descriptor: &LineageRetrievalDescriptor) {
        match descriptor {
            LineageRetrievalDescriptor::GenomeGene {
                scope,
                genome_id,
                gene_query,
                occurrence,
                catalog_path,
                cache_dir,
            } => {
                self.set_genome_dialog_scope(*scope);
                if let Some(path) = catalog_path
                    .as_deref()
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                {
                    self.genome_catalog_path = path.to_string();
                }
                if let Some(path) = cache_dir
                    .as_deref()
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                {
                    self.genome_cache_dir = path.to_string();
                }
                self.sync_active_genome_scope_paths_from_fields();
                if self.genome_id != *genome_id {
                    self.genome_id = genome_id.clone();
                    self.invalidate_genome_genes();
                }
                self.genome_gene_filter = gene_query.clone();
                self.genome_gene_filter_page = 0;
                self.open_retrieve_genome_dialog_for_scope(*scope);
                self.genome_retrieve_status = match occurrence {
                    Some(value) => format!(
                        "Reopened retrieval pattern: gene '{}' (occurrence {}) on '{}'",
                        gene_query, value, genome_id
                    ),
                    None => format!(
                        "Reopened retrieval pattern: gene '{}' on '{}'",
                        gene_query, genome_id
                    ),
                };
                self.genome_retrieve_contig_suggestions.clear();
            }
            LineageRetrievalDescriptor::GenomeRegion {
                scope,
                genome_id,
                chromosome,
                start_1based,
                end_1based,
                catalog_path,
                cache_dir,
            } => {
                self.set_genome_dialog_scope(*scope);
                if let Some(path) = catalog_path
                    .as_deref()
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                {
                    self.genome_catalog_path = path.to_string();
                }
                if let Some(path) = cache_dir
                    .as_deref()
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                {
                    self.genome_cache_dir = path.to_string();
                }
                self.sync_active_genome_scope_paths_from_fields();
                if self.genome_id != *genome_id {
                    self.genome_id = genome_id.clone();
                    self.invalidate_genome_genes();
                }
                self.genome_chromosome = chromosome.clone();
                self.genome_start_1based = start_1based.to_string();
                self.genome_end_1based = end_1based.to_string();
                self.open_retrieve_genome_dialog_for_scope(*scope);
                self.genome_retrieve_status = format!(
                    "Reopened retrieval pattern: {}:{}-{} on '{}'",
                    chromosome, start_1based, end_1based, genome_id
                );
                self.genome_retrieve_contig_suggestions.clear();
            }
            LineageRetrievalDescriptor::GenBankAccession { accession, as_id } => {
                self.genbank_accession = accession.clone();
                self.genbank_as_id = as_id.clone().unwrap_or_default();
                self.open_genbank_dialog();
                self.genbank_status = format!(
                    "Reopened retrieval pattern: GenBank accession '{}'",
                    accession
                );
            }
            LineageRetrievalDescriptor::DbSnpRegion {
                genome_id,
                rs_id,
                flank_bp,
                catalog_path,
                cache_dir,
            } => {
                if let Some(path) = catalog_path
                    .as_deref()
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                {
                    self.reference_genome_catalog_path = path.to_string();
                }
                if let Some(path) = cache_dir
                    .as_deref()
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                {
                    self.reference_genome_cache_dir = path.to_string();
                }
                self.dbsnp_genome_id = genome_id.clone();
                self.dbsnp_rs_id = rs_id.clone();
                self.dbsnp_flank_bp = flank_bp.to_string();
                self.dbsnp_output_id.clear();
                self.open_genbank_dialog();
                self.dbsnp_status = format!(
                    "Reopened retrieval pattern: dbSNP '{}' +/-{} bp on '{}'",
                    rs_id, flank_bp, genome_id
                );
            }
        }
    }

    fn pretty_json_or_fallback<T: Serialize>(value: &T) -> String {
        serde_json::to_string_pretty(value).unwrap_or_else(|_| "{}".to_string())
    }

    fn refresh_planning_editor_buffers_from_engine(&mut self) {
        let (global, overlay, project_override, objective) = {
            let engine = self.engine.read().unwrap();
            (
                engine.planning_profile(PlanningProfileScope::Global),
                engine.planning_profile(PlanningProfileScope::ConfirmedAgentOverlay),
                engine.planning_profile(PlanningProfileScope::ProjectOverride),
                engine.planning_objective(),
            )
        };
        self.planning_profile_global_json = global
            .as_ref()
            .map(Self::pretty_json_or_fallback)
            .unwrap_or_default();
        self.planning_profile_overlay_json = overlay
            .as_ref()
            .map(Self::pretty_json_or_fallback)
            .unwrap_or_default();
        self.planning_profile_project_json = project_override
            .as_ref()
            .map(Self::pretty_json_or_fallback)
            .unwrap_or_default();
        self.planning_objective_json = Self::pretty_json_or_fallback(&objective);
    }

    fn parse_optional_planning_profile_json(
        raw_json: &str,
    ) -> Result<Option<PlanningProfile>, String> {
        let trimmed = raw_json.trim();
        if trimmed.is_empty() {
            return Ok(None);
        }
        serde_json::from_str::<PlanningProfile>(trimmed)
            .map(Some)
            .map_err(|e| format!("Could not parse planning profile JSON: {e}"))
    }

    fn parse_optional_planning_objective_json(
        raw_json: &str,
    ) -> Result<Option<PlanningObjective>, String> {
        let trimmed = raw_json.trim();
        if trimmed.is_empty() {
            return Ok(None);
        }
        serde_json::from_str::<PlanningObjective>(trimmed)
            .map(Some)
            .map_err(|e| format!("Could not parse planning objective JSON: {e}"))
    }

    fn apply_planning_profile_json(&mut self, scope: PlanningProfileScope, raw_json: &str) {
        match Self::parse_optional_planning_profile_json(raw_json) {
            Ok(profile) => {
                let outcome = self
                    .engine
                    .write()
                    .unwrap()
                    .set_planning_profile(scope, profile);
                match outcome {
                    Ok(()) => {
                        self.refresh_planning_editor_buffers_from_engine();
                        self.planning_status =
                            format!("Planning profile '{}' updated", scope.as_str());
                    }
                    Err(e) => {
                        self.planning_status = format!(
                            "Could not update planning profile '{}': {}",
                            scope.as_str(),
                            e.message
                        );
                    }
                }
            }
            Err(err) => {
                self.planning_status = err;
            }
        }
    }

    fn clear_planning_profile_scope(&mut self, scope: PlanningProfileScope) {
        let outcome = self
            .engine
            .write()
            .unwrap()
            .set_planning_profile(scope, None);
        match outcome {
            Ok(()) => {
                self.refresh_planning_editor_buffers_from_engine();
                self.planning_status = format!("Planning profile '{}' cleared", scope.as_str());
            }
            Err(e) => {
                self.planning_status = format!(
                    "Could not clear planning profile '{}': {}",
                    scope.as_str(),
                    e.message
                );
            }
        }
    }

    fn apply_planning_objective_json(&mut self, raw_json: &str) {
        match Self::parse_optional_planning_objective_json(raw_json) {
            Ok(objective) => {
                let outcome = self
                    .engine
                    .write()
                    .unwrap()
                    .set_planning_objective(objective);
                match outcome {
                    Ok(()) => {
                        self.refresh_planning_editor_buffers_from_engine();
                        self.planning_status = "Planning objective updated".to_string();
                    }
                    Err(e) => {
                        self.planning_status =
                            format!("Could not update planning objective: {}", e.message);
                    }
                }
            }
            Err(err) => {
                self.planning_status = err;
            }
        }
    }

    fn clear_planning_objective(&mut self) {
        let outcome = self.engine.write().unwrap().set_planning_objective(None);
        match outcome {
            Ok(()) => {
                self.refresh_planning_editor_buffers_from_engine();
                self.planning_status =
                    "Planning objective cleared (engine defaults active)".to_string();
            }
            Err(e) => {
                self.planning_status = format!("Could not clear planning objective: {}", e.message);
            }
        }
    }

    fn register_planning_sync_suggestion(&mut self, direction: &str) {
        let payload_text = self.planning_sync_payload_json.trim().to_string();
        if payload_text.is_empty() {
            self.planning_status =
                "Planning sync payload is empty; provide JSON with profile_patch/objective_patch"
                    .to_string();
            return;
        }
        let payload = match serde_json::from_str::<PlanningSyncSuggestionPayloadUi>(&payload_text) {
            Ok(value) => value,
            Err(e) => {
                self.planning_status = format!("Could not parse planning sync payload JSON: {e}");
                return;
            }
        };
        let source = {
            let trimmed = self.planning_sync_source.trim();
            if trimmed.is_empty() {
                "lab_manager".to_string()
            } else {
                trimmed.to_string()
            }
        };
        let confidence = if self.planning_sync_confidence.trim().is_empty() {
            None
        } else {
            match self.planning_sync_confidence.trim().parse::<f64>() {
                Ok(value) if value.is_finite() && (0.0..=1.0).contains(&value) => Some(value),
                _ => {
                    self.planning_status =
                        "Planning sync confidence must be a number in range 0.0..=1.0".to_string();
                    return;
                }
            }
        };
        let snapshot_id = self.planning_sync_snapshot_id.trim();
        let snapshot_id = if snapshot_id.is_empty() {
            None
        } else {
            Some(snapshot_id.to_string())
        };
        let message = self.planning_sync_message.trim().to_string();
        let message = if message.is_empty() {
            None
        } else {
            Some(message)
        };
        let outcome = self.engine.write().unwrap().propose_planning_suggestion(
            direction,
            &source,
            confidence,
            snapshot_id.as_deref(),
            payload.profile_patch,
            payload.objective_patch,
            message.as_deref(),
        );
        match outcome {
            Ok(suggestion) => {
                self.refresh_planning_editor_buffers_from_engine();
                self.planning_status = format!(
                    "Registered planning {} suggestion '{}'",
                    direction, suggestion.suggestion_id
                );
            }
            Err(e) => {
                self.planning_status = format!(
                    "Could not register planning {} suggestion: {}",
                    direction, e.message
                );
            }
        }
    }

    fn accept_planning_suggestion(&mut self, suggestion_id: &str) {
        let outcome = self
            .engine
            .write()
            .unwrap()
            .accept_planning_suggestion(suggestion_id);
        match outcome {
            Ok(suggestion) => {
                self.refresh_planning_editor_buffers_from_engine();
                self.planning_status = format!(
                    "Accepted planning suggestion '{}'",
                    suggestion.suggestion_id
                );
            }
            Err(e) => {
                self.planning_status = format!(
                    "Could not accept planning suggestion '{}': {}",
                    suggestion_id, e.message
                );
            }
        }
    }

    fn reject_planning_suggestion(&mut self, suggestion_id: &str) {
        let rejection_reason = self.planning_rejection_reason.trim().to_string();
        let rejection_reason = if rejection_reason.is_empty() {
            None
        } else {
            Some(rejection_reason)
        };
        let outcome = self
            .engine
            .write()
            .unwrap()
            .reject_planning_suggestion(suggestion_id, rejection_reason.as_deref());
        match outcome {
            Ok(suggestion) => {
                self.refresh_planning_editor_buffers_from_engine();
                self.planning_status = format!(
                    "Rejected planning suggestion '{}'",
                    suggestion.suggestion_id
                );
            }
            Err(e) => {
                self.planning_status = format!(
                    "Could not reject planning suggestion '{}': {}",
                    suggestion_id, e.message
                );
            }
        }
    }

    fn planning_filter_label(filter: Option<PlanningSuggestionStatus>) -> &'static str {
        match filter {
            None => "all",
            Some(PlanningSuggestionStatus::Pending) => "pending",
            Some(PlanningSuggestionStatus::Accepted) => "accepted",
            Some(PlanningSuggestionStatus::Rejected) => "rejected",
        }
    }

    fn track_name_default_from_path(path: &Path) -> String {
        path.file_name()
            .map(|value| value.to_string_lossy().to_string())
            .unwrap_or_else(|| path.display().to_string())
    }

    fn subscription_matches_filter(
        subscription: &GenomeTrackSubscription,
        filter_text: &str,
    ) -> bool {
        let terms = Self::parse_track_filter_terms(filter_text);
        if terms.is_empty() {
            return true;
        }
        let source_terms = vec![subscription.source.label().to_ascii_lowercase()];
        let path_terms = {
            let mut values = vec![subscription.path.to_ascii_lowercase()];
            if let Some(file_name) = Path::new(subscription.path.as_str())
                .file_name()
                .and_then(|value| value.to_str())
            {
                values.push(file_name.to_ascii_lowercase());
            }
            values
        };
        let track_terms = subscription
            .track_name
            .as_ref()
            .map(|value| vec![value.to_ascii_lowercase()])
            .unwrap_or_default();
        let all_terms = source_terms
            .iter()
            .chain(path_terms.iter())
            .chain(track_terms.iter())
            .cloned()
            .collect::<Vec<_>>();

        terms.iter().all(|(scope, needle)| {
            let search_space: &Vec<String> = match scope.as_deref() {
                Some("source") | Some("type") => &source_terms,
                Some("path") | Some("file") => &path_terms,
                Some("track") | Some("name") => &track_terms,
                _ => &all_terms,
            };
            search_space.iter().any(|value| value.contains(needle))
        })
    }

    fn parse_track_filter_terms(filter_text: &str) -> Vec<(Option<String>, String)> {
        filter_text
            .split_whitespace()
            .filter_map(|raw_term| {
                let term = raw_term.trim();
                if term.is_empty() {
                    return None;
                }
                if let Some((raw_scope, raw_needle)) = term.split_once(':') {
                    let scope = raw_scope.trim().to_ascii_lowercase();
                    let needle = raw_needle.trim().to_ascii_lowercase();
                    if !scope.is_empty() && !needle.is_empty() {
                        return Some((Some(scope), needle));
                    }
                }
                Some((None, term.to_ascii_lowercase()))
            })
            .collect()
    }

    fn tracked_file_filter_help_text() -> &'static str {
        "Free text matches source/path/track. Scoped terms: source:BED source:VCF path:peaks.bed track:chip"
    }

    fn append_filter_term(filter_text: &mut String, term: &str) {
        let normalized = term.trim();
        if normalized.is_empty() {
            return;
        }
        if filter_text
            .split_whitespace()
            .any(|existing| existing.eq_ignore_ascii_case(normalized))
        {
            return;
        }
        if !filter_text.trim().is_empty() {
            filter_text.push(' ');
        }
        filter_text.push_str(normalized);
    }

    fn open_helper_genome_prepare_dialog(&mut self) {
        self.open_prepare_genome_dialog_for_scope(GenomeDialogScope::Helper);
    }

    fn open_helper_genome_retrieve_dialog(&mut self) {
        self.open_retrieve_genome_dialog_for_scope(GenomeDialogScope::Helper);
    }

    fn open_helper_genome_blast_dialog(&mut self) {
        self.open_blast_genome_dialog_for_scope(GenomeDialogScope::Helper);
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
            .unwrap_or_else(|| self.genome_dialog_scope.default_catalog_path().to_string())
    }

    fn genome_cache_dir_opt(&self) -> Option<String> {
        let trimmed = self.genome_cache_dir.trim();
        if trimmed.is_empty() {
            None
        } else {
            Some(trimmed.to_string())
        }
    }

    fn dbsnp_catalog_path_opt(&self) -> Option<String> {
        let trimmed = self.reference_genome_catalog_path.trim();
        if trimmed.is_empty() {
            None
        } else {
            Some(trimmed.to_string())
        }
    }

    fn dbsnp_cache_dir_opt(&self) -> Option<String> {
        let trimmed = self.reference_genome_cache_dir.trim();
        if trimmed.is_empty() {
            None
        } else {
            Some(trimmed.to_string())
        }
    }

    fn seed_dbsnp_defaults_if_needed(&mut self) {
        if self.dbsnp_rs_id.trim().is_empty() {
            self.dbsnp_rs_id = DEFAULT_DBSNP_TUTORIAL_RS_ID.to_string();
        }
        if self.dbsnp_genome_id.trim().is_empty() {
            if self.genome_id.trim().is_empty() {
                self.dbsnp_genome_id = "Human GRCh38 Ensembl 113".to_string();
            } else {
                self.dbsnp_genome_id = self.genome_id.trim().to_string();
            }
        }
        if self.dbsnp_flank_bp.trim().is_empty() {
            self.dbsnp_flank_bp = "3000".to_string();
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

    fn latest_genome_anchor_provenance_by_seq(
        state: &ProjectState,
    ) -> HashMap<String, AnchorProvenanceSnapshot> {
        let mut out: HashMap<String, AnchorProvenanceSnapshot> = HashMap::new();
        let Some(entries) = state
            .metadata
            .get("provenance")
            .and_then(|v| v.get("genome_extractions"))
            .and_then(|v| v.as_array())
        else {
            return out;
        };
        for entry in entries {
            let Some(seq_id) = entry
                .get("seq_id")
                .and_then(|v| v.as_str())
                .map(str::trim)
                .filter(|v| !v.is_empty())
            else {
                continue;
            };
            let recorded_at_unix_ms = entry
                .get("recorded_at_unix_ms")
                .and_then(|v| v.as_u64())
                .map(|v| v as u128)
                .unwrap_or(0);
            let snapshot = AnchorProvenanceSnapshot {
                catalog_path: entry
                    .get("catalog_path")
                    .and_then(|v| v.as_str())
                    .map(str::trim)
                    .filter(|v| !v.is_empty())
                    .map(|v| v.to_string()),
                cache_dir: entry
                    .get("cache_dir")
                    .and_then(|v| v.as_str())
                    .map(str::trim)
                    .filter(|v| !v.is_empty())
                    .map(|v| v.to_string()),
                recorded_at_unix_ms,
            };
            let replace = out
                .get(seq_id)
                .map(|existing| recorded_at_unix_ms >= existing.recorded_at_unix_ms)
                .unwrap_or(true);
            if replace {
                out.insert(seq_id.to_string(), snapshot);
            }
        }
        out
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

    fn suggest_chromosome_names(query: &str, available: &[String], limit: usize) -> Vec<String> {
        if limit == 0 {
            return vec![];
        }
        let query_trimmed = query.trim();
        if query_trimmed.is_empty() {
            return vec![];
        }
        let query_lower = query_trimmed.to_ascii_lowercase();
        let query_normalized = Self::normalize_chromosome_alias(query_trimmed);
        let mut scored: Vec<(u8, String)> = available
            .iter()
            .filter_map(|candidate| {
                let candidate_trimmed = candidate.trim();
                if candidate_trimmed.is_empty() {
                    return None;
                }
                let candidate_lower = candidate_trimmed.to_ascii_lowercase();
                let candidate_normalized = Self::normalize_chromosome_alias(candidate_trimmed);
                let query_numeric = !query_normalized.is_empty()
                    && query_normalized.chars().all(|ch| ch.is_ascii_digit());
                let candidate_numeric = !candidate_normalized.is_empty()
                    && candidate_normalized.chars().all(|ch| ch.is_ascii_digit());
                let allow_shorter_partial = !(query_numeric && candidate_numeric);
                let score =
                    if !query_normalized.is_empty() && candidate_normalized == query_normalized {
                        0u8
                    } else if !query_normalized.is_empty()
                        && (candidate_normalized.starts_with(&query_normalized)
                            || (allow_shorter_partial
                                && query_normalized.starts_with(&candidate_normalized)))
                    {
                        1u8
                    } else if !query_lower.is_empty()
                        && (candidate_lower.contains(&query_lower)
                            || (allow_shorter_partial && query_lower.contains(&candidate_lower)))
                    {
                        2u8
                    } else {
                        return None;
                    };
                Some((score, candidate_trimmed.to_string()))
            })
            .collect();
        scored.sort_by(|left, right| {
            left.0
                .cmp(&right.0)
                .then(left.1.len().cmp(&right.1.len()))
                .then(left.1.cmp(&right.1))
        });
        let mut out: Vec<String> = vec![];
        for (_, candidate) in scored {
            if out.iter().any(|existing| existing == &candidate) {
                continue;
            }
            out.push(candidate);
            if out.len() >= limit {
                break;
            }
        }
        out
    }

    fn chromosome_aliases_match(left: &str, right: &str) -> bool {
        Self::normalize_chromosome_alias(left) == Self::normalize_chromosome_alias(right)
    }

    fn refresh_retrieve_contig_suggestions(&mut self, query_chromosome: &str) {
        self.genome_retrieve_contig_suggestions.clear();
        let genome_id = self.genome_id.trim();
        if genome_id.is_empty() {
            return;
        }
        let query = query_chromosome.trim();
        if query.is_empty() {
            return;
        }
        let available = match self.prepared_genome_chromosome_records(genome_id) {
            Ok(records) => records
                .into_iter()
                .map(|record| record.chromosome)
                .collect::<Vec<_>>(),
            Err(_) => return,
        };
        self.genome_retrieve_contig_suggestions =
            Self::suggest_chromosome_names(query, &available, 8);
    }

    fn update_retrieve_contig_suggestions_from_error(
        &mut self,
        query_chromosome: &str,
        error_message: &str,
    ) {
        self.genome_retrieve_contig_suggestions.clear();
        if error_message.contains("Chromosome/contig")
            && error_message.contains("not found in genome")
        {
            self.refresh_retrieve_contig_suggestions(query_chromosome);
        }
    }

    fn format_genome_anchor_summary(anchor: &SequenceGenomeAnchorSummary) -> String {
        let strand = anchor.strand.unwrap_or('+');
        format!(
            "{} | {}:{}-{} (strand {})",
            anchor.genome_id, anchor.chromosome, anchor.start_1based, anchor.end_1based, strand
        )
    }

    fn chromosome_length_for_anchor(
        anchor: &SequenceGenomeAnchorSummary,
        provenance: Option<&AnchorProvenanceSnapshot>,
        chromosome_length_cache: &mut HashMap<
            GenomeLengthCacheKey,
            Option<Vec<GenomeChromosomeRecord>>,
        >,
    ) -> Option<usize> {
        let provenance = provenance?;
        let catalog_path = provenance.catalog_path.as_deref()?.trim();
        if catalog_path.is_empty() {
            return None;
        }
        let cache_dir = provenance
            .cache_dir
            .as_ref()
            .map(|v| v.trim().to_string())
            .filter(|v| !v.is_empty());
        let key = GenomeLengthCacheKey {
            catalog_path: catalog_path.to_string(),
            genome_id: anchor.genome_id.clone(),
            cache_dir: cache_dir.clone(),
        };
        let lengths = chromosome_length_cache
            .entry(key.clone())
            .or_insert_with(|| {
                let catalog = GenomeCatalog::from_json_file(&key.catalog_path).ok()?;
                catalog
                    .list_chromosome_lengths(&key.genome_id, key.cache_dir.as_deref())
                    .ok()
            })
            .as_ref()?;
        lengths
            .iter()
            .find(|record| Self::chromosome_aliases_match(&record.chromosome, &anchor.chromosome))
            .map(|record| record.length_bp)
    }

    fn sequence_represents_full_genome(
        anchor: &SequenceGenomeAnchorSummary,
        seq_len: usize,
        provenance: Option<&AnchorProvenanceSnapshot>,
        chromosome_length_cache: &mut HashMap<
            GenomeLengthCacheKey,
            Option<Vec<GenomeChromosomeRecord>>,
        >,
    ) -> bool {
        if anchor.start_1based != 1 {
            return false;
        }
        let Some(span_len) = anchor
            .end_1based
            .checked_sub(anchor.start_1based)
            .and_then(|delta| delta.checked_add(1))
        else {
            return false;
        };
        if span_len != seq_len {
            return false;
        }
        let Some(chrom_length) =
            Self::chromosome_length_for_anchor(anchor, provenance, chromosome_length_cache)
        else {
            return false;
        };
        chrom_length == seq_len && anchor.end_1based == chrom_length
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
                let refreshed = self.refresh_sequence_windows_from_engine_state();
                if refreshed > 0 {
                    self.genome_track_status.push_str(&format!(
                        " | refreshed {} open sequence window(s)",
                        refreshed
                    ));
                }
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
                let refreshed = self.refresh_sequence_windows_from_engine_state();
                if refreshed > 0 {
                    self.genome_track_status.push_str(&format!(
                        " | refreshed {} open sequence window(s)",
                        refreshed
                    ));
                }
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

    fn invalidate_genome_genes(&mut self) {
        self.genome_genes.clear();
        self.genome_genes_error.clear();
        self.genome_genes_loaded_key = None;
        self.genome_selected_gene = None;
        self.genome_gene_filter_page = 0;
        self.genome_biotype_filter.clear();
        self.genome_retrieve_contig_suggestions.clear();
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

    fn prepare_dialog_primary_action(
        genome_id: &str,
        all_genomes: &[String],
        preparable_set: &HashSet<String>,
    ) -> PrepareGenomeDialogPrimaryAction {
        let trimmed = genome_id.trim();
        if trimmed.is_empty() {
            PrepareGenomeDialogPrimaryAction::None
        } else if preparable_set.contains(trimmed) {
            PrepareGenomeDialogPrimaryAction::Prepare
        } else if all_genomes.iter().any(|name| name == trimmed) {
            PrepareGenomeDialogPrimaryAction::Reindex
        } else {
            PrepareGenomeDialogPrimaryAction::None
        }
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

    fn selected_genome_catalog_is_writable(&self) -> bool {
        let catalog_path = self.genome_catalog_path_resolved();
        GenomeCatalog::from_json_file(&catalog_path)
            .map(|catalog| catalog.catalog_file_is_writable())
            .unwrap_or(false)
    }

    fn selected_genome_catalog_has_ensembl_templates(&self) -> bool {
        let catalog_path = self.genome_catalog_path_resolved();
        GenomeCatalog::from_json_file(&catalog_path)
            .map(|catalog| catalog.has_ensembl_updatable_entries())
            .unwrap_or(false)
    }

    fn queue_ensembl_catalog_update_preview(&mut self) {
        let catalog_path = self.genome_catalog_path_resolved();
        match GentleEngine::preview_reference_genome_ensembl_catalog_updates(Some(&catalog_path)) {
            Ok(preview) => {
                self.pending_ensembl_catalog_update = Some(PendingEnsemblCatalogUpdateDialog {
                    scope: self.genome_dialog_scope,
                    catalog_path: catalog_path.clone(),
                    preview,
                    output_catalog_path: String::new(),
                });
                self.genome_prepare_status = format!(
                    "Loaded Ensembl catalog update preview for '{}'.",
                    catalog_path
                );
            }
            Err(e) => {
                self.genome_prepare_status = format!(
                    "Could not preview Ensembl catalog updates for '{}': {}",
                    catalog_path, e
                );
                self.pending_ensembl_catalog_update = None;
            }
        }
    }

    fn queue_ensembl_installable_genome_discovery(
        &mut self,
        collection_filter: Option<&str>,
        filter: Option<&str>,
    ) {
        let collection_filter = collection_filter
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .unwrap_or("all");
        let filter = filter
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .unwrap_or("");
        match GentleEngine::discover_ensembl_installable_genomes(
            Some(collection_filter),
            if filter.is_empty() {
                None
            } else {
                Some(filter)
            },
        ) {
            Ok(report) => {
                self.pending_ensembl_installable_genomes =
                    Some(PendingEnsemblInstallableGenomeDialog {
                        scope: self.genome_dialog_scope,
                        collection_filter: collection_filter.to_string(),
                        filter: filter.to_string(),
                        report,
                    });
                self.genome_prepare_status = format!(
                    "Loaded Ensembl installable-genome candidates for collection '{}'.",
                    collection_filter
                );
            }
            Err(e) => {
                self.genome_prepare_status = format!(
                    "Could not inspect current Ensembl installable-genome candidates: {}",
                    e
                );
                self.pending_ensembl_installable_genomes = None;
            }
        }
    }

    fn queue_ensembl_quick_install_preview(
        &mut self,
        scope: GenomeDialogScope,
        collection: &str,
        species_dir: &str,
    ) {
        let catalog_path = self.genome_catalog_path_opt();
        let preview = match scope {
            GenomeDialogScope::Reference => {
                GentleEngine::preview_reference_genome_ensembl_quick_install(
                    catalog_path.as_deref(),
                    collection,
                    species_dir,
                    None,
                    None,
                )
            }
            GenomeDialogScope::Helper => GentleEngine::preview_helper_genome_ensembl_quick_install(
                catalog_path.as_deref(),
                collection,
                species_dir,
                None,
                None,
            ),
        };
        match preview {
            Ok(preview) => {
                self.pending_ensembl_quick_install = Some(PendingEnsemblQuickInstallDialog {
                    scope,
                    collection: collection.trim().to_string(),
                    species_dir: species_dir.trim().to_string(),
                    genome_id: preview.genome_id.clone(),
                    output_catalog_path: preview.output_catalog_path.clone(),
                    preview,
                });
                self.genome_prepare_status = format!(
                    "Prepared Ensembl quick-install preview for '{}' ({})",
                    species_dir, collection
                );
            }
            Err(e) => {
                self.pending_ensembl_quick_install = None;
                self.genome_prepare_status = format!(
                    "Could not preview Ensembl quick install for '{}' ({}): {}",
                    species_dir, collection, e
                );
            }
        }
    }

    fn apply_pending_ensembl_quick_install(&mut self) {
        let Some(dialog) = self.pending_ensembl_quick_install.clone() else {
            return;
        };
        let genome_id = dialog.genome_id.trim().to_string();
        if genome_id.is_empty() {
            self.genome_prepare_status =
                "Choose a non-empty genome id for the Ensembl quick install.".to_string();
            return;
        }
        let output_catalog_path = dialog.output_catalog_path.trim().to_string();
        let catalog_path = self.genome_catalog_path_opt();
        let write_report = match dialog.scope {
            GenomeDialogScope::Reference => {
                GentleEngine::apply_reference_genome_ensembl_quick_install(
                    catalog_path.as_deref(),
                    &dialog.collection,
                    &dialog.species_dir,
                    if output_catalog_path.is_empty() {
                        None
                    } else {
                        Some(output_catalog_path.as_str())
                    },
                    Some(genome_id.as_str()),
                )
            }
            GenomeDialogScope::Helper => GentleEngine::apply_helper_genome_ensembl_quick_install(
                catalog_path.as_deref(),
                &dialog.collection,
                &dialog.species_dir,
                if output_catalog_path.is_empty() {
                    None
                } else {
                    Some(output_catalog_path.as_str())
                },
                Some(genome_id.as_str()),
            ),
        };
        match write_report {
            Ok(report) => {
                self.pending_ensembl_quick_install = None;
                self.genome_dialog_scope = dialog.scope;
                self.genome_catalog_path = report.preview.output_catalog_path.clone();
                self.genome_id = report.preview.genome_id.clone();
                self.genome_prepare_status = format!(
                    "Wrote Ensembl quick-install catalog entry '{}' to '{}'. Starting prepare.",
                    report.preview.genome_id, report.preview.output_catalog_path
                );
                self.start_prepare_reference_genome_with_mode(GenomePrepareLaunchMode::Prepare);
            }
            Err(e) => {
                self.genome_prepare_status = format!(
                    "Could not apply Ensembl quick install for '{}' ({}): {}",
                    dialog.species_dir, dialog.collection, e
                );
            }
        }
    }

    fn clear_prepare_dialog_ephemeral_state(&mut self) {
        self.pending_ensembl_catalog_update = None;
        self.pending_ensembl_installable_genomes = None;
        self.pending_ensembl_quick_install = None;
        self.clear_prepare_step_state();
    }

    fn apply_pending_ensembl_catalog_update(&mut self) {
        let Some(dialog) = self.pending_ensembl_catalog_update.clone() else {
            return;
        };
        let output_catalog_path = dialog.output_catalog_path.trim().to_string();
        let output_catalog_path = if output_catalog_path.is_empty() {
            None
        } else {
            Some(output_catalog_path.as_str())
        };
        match GentleEngine::apply_reference_genome_ensembl_catalog_updates(
            Some(&dialog.catalog_path),
            output_catalog_path,
        ) {
            Ok(report) => {
                self.pending_ensembl_catalog_update = None;
                self.genome_prepare_status = format!(
                    "Ensembl catalog specs updated: {} update(s) written to '{}'.",
                    report.updates.len(),
                    report.output_catalog_path
                );
                if report.output_catalog_path == dialog.catalog_path {
                    self.refresh_genome_catalog_list();
                }
            }
            Err(e) => {
                self.genome_prepare_status = format!(
                    "Could not apply Ensembl catalog updates for '{}': {}",
                    dialog.catalog_path, e
                );
            }
        }
    }

    fn queue_prepared_genome_removal(
        &mut self,
        genome_id: String,
        scope: GenomeDialogScope,
        install_dir: String,
    ) {
        self.pending_prepared_genome_removal = Some(PendingPreparedGenomeRemovalRequest {
            genome_id,
            scope,
            catalog_path: self.genome_catalog_path_resolved(),
            cache_dir: self.genome_cache_dir_opt(),
            install_dir,
        });
    }

    fn confirm_prepared_genome_removal(&mut self) {
        let Some(request) = self.pending_prepared_genome_removal.take() else {
            return;
        };
        match GentleEngine::remove_prepared_reference_genome(
            Some(&request.catalog_path),
            &request.genome_id,
            request.cache_dir.as_deref(),
        ) {
            Ok(report) => {
                self.genome_prepare_status = if report.removed {
                    format!(
                        "Removed prepared {} '{}' from '{}'.",
                        request.scope.description(),
                        request.genome_id,
                        report.install_dir
                    )
                } else {
                    format!(
                        "No prepared {} files were present for '{}'.",
                        request.scope.description(),
                        request.genome_id
                    )
                };
            }
            Err(e) => {
                self.genome_prepare_status = format!(
                    "Could not remove prepared {} '{}': {}",
                    request.scope.description(),
                    request.genome_id,
                    e
                );
            }
        }
    }

    fn queue_catalog_entry_removal(&mut self, genome_id: String, scope: GenomeDialogScope) {
        self.pending_genome_catalog_entry_removal = Some(PendingGenomeCatalogEntryRemovalRequest {
            genome_id,
            scope,
            catalog_path: self.genome_catalog_path_resolved(),
        });
    }

    fn confirm_catalog_entry_removal(&mut self) {
        let Some(request) = self.pending_genome_catalog_entry_removal.take() else {
            return;
        };
        match GentleEngine::remove_reference_genome_catalog_entry(
            Some(&request.catalog_path),
            &request.genome_id,
            None,
        ) {
            Ok(report) => {
                self.genome_prepare_status = format!(
                    "Removed {} catalog entry '{}' from '{}'. Prepared cache files, if any, were left unchanged.",
                    request.scope.description(),
                    request.genome_id,
                    report.output_catalog_path
                );
                self.refresh_genome_catalog_list();
            }
            Err(e) => {
                self.genome_prepare_status = format!(
                    "Could not remove {} catalog entry '{}': {}",
                    request.scope.description(),
                    request.genome_id,
                    e
                );
            }
        }
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

    fn host_profile_catalog_entries(
        &self,
        catalog_path: Option<&str>,
        filter: Option<&str>,
    ) -> Result<Vec<HostProfileRecord>, String> {
        let resolved_catalog_path = catalog_path
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .unwrap_or(DEFAULT_HOST_PROFILE_CATALOG_PATH);
        GentleEngine::list_host_profile_catalog_entries(catalog_path, filter).map_err(|e| {
            format!(
                "Could not load host profile catalog entries from '{}': {}",
                resolved_catalog_path, e.message
            )
        })
    }

    fn host_profile_matches_selection(profile: &HostProfileRecord, selected: &str) -> bool {
        let needle = selected.trim();
        if needle.is_empty() {
            return false;
        }
        profile.profile_id.eq_ignore_ascii_case(needle)
            || profile
                .aliases
                .iter()
                .any(|alias| alias.eq_ignore_ascii_case(needle))
            || profile.strain.eq_ignore_ascii_case(needle)
    }

    fn helper_catalog_entries_for_catalog_path(
        &self,
        catalog_path: &str,
        filter: Option<&str>,
    ) -> Result<Vec<GenomeCatalogListEntry>, String> {
        GentleEngine::list_helper_catalog_entries(Some(catalog_path), filter).map_err(|e| {
            format!(
                "Could not load helper catalog entries from '{}': {}",
                catalog_path, e.message
            )
        })
    }

    fn helper_vector_cards_for_catalog_path(
        &self,
        catalog_path: &str,
        filter: Option<&str>,
    ) -> Result<Vec<HelperVectorCard>, String> {
        GentleEngine::list_helper_vector_cards(Some(catalog_path), filter)
            .map(|report| report.cards)
            .map_err(|e| {
                format!(
                    "Could not load helper vector cards from '{}': {}",
                    catalog_path, e.message
                )
            })
    }

    fn helper_vector_doctor_issues_for_catalog_path(
        &self,
        catalog_path: &str,
    ) -> Result<Vec<HelperVectorCatalogDoctorIssue>, String> {
        GentleEngine::doctor_helper_vector_catalog(Some(catalog_path))
            .map(|report| report.issues)
            .map_err(|e| {
                format!(
                    "Could not doctor helper catalog '{}': {}",
                    catalog_path, e.message
                )
            })
    }

    fn helper_catalog_entries_for_scope(
        &self,
        scope: GenomeDialogScope,
        filter: Option<&str>,
    ) -> Result<Vec<GenomeCatalogListEntry>, String> {
        if !matches!(scope, GenomeDialogScope::Helper) {
            return Ok(vec![]);
        }
        let catalog_path = if scope == self.genome_dialog_scope {
            self.genome_catalog_path_resolved()
        } else {
            self.scope_genome_paths_resolved(scope).0
        };
        self.helper_catalog_entries_for_catalog_path(&catalog_path, filter)
    }

    fn helper_catalog_entry_matches_selection(
        entry: &GenomeCatalogListEntry,
        selected: &str,
    ) -> bool {
        let needle = selected.trim();
        if needle.is_empty() {
            return false;
        }
        entry.genome_id.eq_ignore_ascii_case(needle)
            || entry
                .aliases
                .iter()
                .any(|alias| alias.eq_ignore_ascii_case(needle))
    }

    fn selected_helper_catalog_entry_for_scope(
        &self,
        scope: GenomeDialogScope,
    ) -> Result<Option<GenomeCatalogListEntry>, String> {
        if !matches!(scope, GenomeDialogScope::Helper) {
            return Ok(None);
        }
        let selected = self.genome_id.trim();
        if selected.is_empty() {
            return Ok(None);
        }
        let entries = self.helper_catalog_entries_for_scope(scope, None)?;
        Ok(entries
            .into_iter()
            .find(|entry| Self::helper_catalog_entry_matches_selection(entry, selected)))
    }

    fn format_helper_construct_interpretation_detail_lines(
        interpretation: &HelperConstructInterpretation,
    ) -> Vec<String> {
        let mut lines = vec![format!("helper id: {}", interpretation.helper_id)];
        if !interpretation.aliases.is_empty() {
            lines.push(format!("aliases: {}", interpretation.aliases.join(", ")));
        }
        if !interpretation.helper_kinds.is_empty() {
            lines.push(format!(
                "helper kind: {}",
                interpretation.helper_kinds.join(", ")
            ));
        }
        if !interpretation.host_systems.is_empty() {
            lines.push(format!(
                "host system: {}",
                interpretation.host_systems.join(", ")
            ));
        }
        if !interpretation.offered_functions.is_empty() {
            lines.push(format!(
                "offered functions: {}",
                interpretation.offered_functions.join(", ")
            ));
        }
        if !interpretation.constraints.is_empty() {
            lines.push(format!(
                "constraints: {}",
                interpretation.constraints.join(", ")
            ));
        }
        if !interpretation.procurement_channels.is_empty() {
            lines.push(format!(
                "procurement channels: {}",
                interpretation.procurement_channels.join(", ")
            ));
        }
        if !interpretation.routine_hints.is_empty() {
            let families = interpretation
                .routine_hints
                .iter()
                .map(|hint| hint.family.as_str())
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect::<Vec<_>>();
            lines.push(format!("routine hints: {}", families.join(", ")));
        }
        if !interpretation.normalized_terms.is_empty() {
            let preview = interpretation
                .normalized_terms
                .iter()
                .take(6)
                .map(|term| {
                    let suffix = term
                        .vocabulary_label
                        .as_deref()
                        .filter(|label| *label != term.value)
                        .map(|label| format!(" ({label})"))
                        .unwrap_or_default();
                    format!("{}={}{}", term.axis, term.value, suffix)
                })
                .collect::<Vec<_>>();
            let suffix = if interpretation.normalized_terms.len() > preview.len() {
                format!(
                    " (+{} more)",
                    interpretation.normalized_terms.len() - preview.len()
                )
            } else {
                String::new()
            };
            lines.push(format!(
                "normalized terms: {}{}",
                preview.join(", "),
                suffix
            ));
        }
        if interpretation.local_variant_unpublished {
            lines.push("local variant: unpublished".to_string());
        }
        lines
    }

    fn format_helper_vector_doctor_issue(issue: &HelperVectorCatalogDoctorIssue) -> String {
        let field = issue
            .field
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|field| format!(" field={field}"))
            .unwrap_or_default();
        format!(
            "{} [{}] {}{}: {}",
            issue.helper_id, issue.severity, issue.code, field, issue.message
        )
    }

    fn format_helper_component_summary(
        component: &crate::genomes::HelperConstructInterpretationComponent,
    ) -> String {
        let title = component
            .label
            .as_deref()
            .filter(|value| !value.trim().is_empty())
            .unwrap_or(component.id.as_str());
        let mut parts = vec![format!("{title} [{}]", component.kind)];
        if title != component.id {
            parts.push(format!("id={}", component.id));
        }
        if !component.tags.is_empty() {
            parts.push(format!("tags={}", component.tags.join(",")));
        }
        if !component.attributes.is_empty() {
            let attrs = component
                .attributes
                .iter()
                .map(|(key, value)| format!("{key}={value}"))
                .collect::<Vec<_>>()
                .join(", ");
            parts.push(format!("attrs={attrs}"));
        }
        if let Some(description) = component
            .description
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            parts.push(description.to_string());
        }
        parts.join(" | ")
    }

    fn format_helper_relationship_summary(
        relationship: &crate::genomes::HelperConstructRelationship,
    ) -> String {
        let mut text = format!(
            "{} --{}--> {}",
            relationship.subject, relationship.predicate, relationship.object
        );
        if let Some(note) = relationship
            .note
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            text.push_str(&format!(" ({note})"));
        }
        text
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

    fn format_extract_region_like_status(prefix: &str, result: &OpResult) -> String {
        let mut text = Self::format_op_result_status(
            prefix,
            &result.created_seq_ids,
            &result.warnings,
            &result.messages,
        );
        if let Some(telemetry) = result.genome_annotation_projection.as_ref() {
            let cap = telemetry
                .max_features_cap
                .map(|value| value.to_string())
                .unwrap_or_else(|| "-".to_string());
            text.push_str(&format!(
                "\nannotation: requested={} effective={} attached={} dropped={} cap={}",
                telemetry.requested_scope,
                telemetry.effective_scope,
                telemetry.attached_feature_count,
                telemetry.dropped_feature_count,
                cap
            ));
            text.push_str(&format!(
                "\nannotation kinds: genes={} transcripts={} exons={} cds={}",
                telemetry.genes_attached,
                telemetry.transcripts_attached,
                telemetry.exons_attached,
                telemetry.cds_attached
            ));
            if let Some(reason) = telemetry
                .fallback_reason
                .as_deref()
                .map(str::trim)
                .filter(|value| !value.is_empty())
            {
                text.push_str(&format!("\nannotation fallback reason: {reason}"));
            }
        }
        text
    }

    fn format_extract_region_status(result: &OpResult) -> String {
        Self::format_extract_region_like_status("Extract region: ok", result)
    }

    fn format_dbsnp_fetch_status(result: &OpResult) -> String {
        Self::format_extract_region_like_status("dbSNP fetch: ok", result)
    }

    fn format_dbsnp_progress_status(progress: &DbSnpFetchProgress) -> String {
        match progress.stage {
            DbSnpFetchStage::ValidateInput => progress.detail.clone(),
            DbSnpFetchStage::InspectPreparedGenome => progress.detail.clone(),
            DbSnpFetchStage::ContactServer => progress.detail.clone(),
            DbSnpFetchStage::WaitResponse => progress.detail.clone(),
            DbSnpFetchStage::ParseResponse => progress.detail.clone(),
            DbSnpFetchStage::ResolvePlacement => progress.detail.clone(),
            DbSnpFetchStage::ExtractRegion => progress.detail.clone(),
            DbSnpFetchStage::AttachVariantMarker => progress.detail.clone(),
        }
    }

    fn poll_dbsnp_fetch_task(&mut self, ctx: &egui::Context) {
        if self.dbsnp_fetch_task.is_none() {
            return;
        }
        ctx.request_repaint_after(Duration::from_millis(100));
        let mut done: Option<Result<OpResult, EngineError>> = None;
        if let Some(task) = &self.dbsnp_fetch_task {
            const MAX_MESSAGES_PER_TICK: usize = 64;
            for _ in 0..MAX_MESSAGES_PER_TICK {
                match task.receiver.try_recv() {
                    Ok(DbSnpFetchTaskMessage::Progress(progress)) => {
                        self.dbsnp_status = Self::format_dbsnp_progress_status(&progress);
                    }
                    Ok(DbSnpFetchTaskMessage::Done(result)) => {
                        done = Some(result);
                        break;
                    }
                    Err(mpsc::TryRecvError::Empty) => break,
                    Err(mpsc::TryRecvError::Disconnected) => {
                        done = Some(Err(EngineError {
                            code: ErrorCode::Io,
                            message: "dbSNP fetch worker disconnected".to_string(),

                            cause_chain: vec![],
                        }));
                        break;
                    }
                }
            }
        }

        if let Some(result) = done {
            let elapsed = self
                .dbsnp_fetch_task
                .as_ref()
                .map(|task| task.started.elapsed().as_secs_f64())
                .unwrap_or(0.0);
            self.dbsnp_fetch_task = None;
            match result {
                Ok(result) => {
                    for seq_id in &result.created_seq_ids {
                        self.open_sequence_window(seq_id);
                    }
                    self.dbsnp_status = format!(
                        "{}\ncompleted in {:.1}s",
                        Self::format_dbsnp_fetch_status(&result),
                        elapsed
                    );
                }
                Err(err) => {
                    self.dbsnp_status =
                        format!("dbSNP fetch failed after {:.1}s: {}", elapsed, err.message);
                }
            }
        }
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

    fn summarize_binary_probe(
        probe: &crate::genomes::ExternalBinaryPreflightProbe,
        include_version: bool,
    ) -> String {
        let version = if include_version {
            probe
                .version
                .as_ref()
                .map(|v| format!(" version={}", v))
                .unwrap_or_default()
        } else {
            String::new()
        };
        if probe.found {
            let path = probe
                .resolved_path
                .as_deref()
                .unwrap_or(probe.executable.as_str());
            let probe_note = if probe.version_probe_ok {
                ""
            } else {
                " (version probe failed)"
            };
            format!("found path={}{}{}", path, version, probe_note)
        } else {
            let error = probe.error.as_deref().unwrap_or("not found");
            format!("missing executable={} ({})", probe.executable, error)
        }
    }

    fn format_prepare_failure_status(
        message: &str,
        elapsed_secs: f64,
        cancellation_requested: bool,
        timeout_seconds: Option<u64>,
        mode: GenomePrepareLaunchMode,
    ) -> String {
        let action = match mode {
            GenomePrepareLaunchMode::Prepare => "Prepare genome",
            GenomePrepareLaunchMode::ReindexCachedFiles => "Reindex genome",
            GenomePrepareLaunchMode::RefreshFromSources => "Refresh genome",
        };
        let lower = message.to_ascii_lowercase();
        let timed_out_hint = lower.contains("timed out") || lower.contains("timeout");
        let cancelled_hint = lower.contains("cancelled") || lower.contains("canceled");
        let timebox_elapsed = timeout_seconds
            .map(|limit| elapsed_secs >= limit as f64)
            .unwrap_or(false);
        if timed_out_hint || (!cancellation_requested && cancelled_hint && timebox_elapsed) {
            format!("{action} timed out after {:.1}s: {}", elapsed_secs, message)
        } else if cancellation_requested || cancelled_hint {
            format!("{action} cancelled after {:.1}s: {}", elapsed_secs, message)
        } else {
            format!("{action} failed after {:.1}s: {}", elapsed_secs, message)
        }
    }

    fn prepare_plan_for_mode(
        &self,
        genome_id: &str,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
        mode: GenomePrepareLaunchMode,
    ) -> Result<PrepareGenomePlan, EngineError> {
        match mode {
            GenomePrepareLaunchMode::Prepare => {
                GentleEngine::prepare_reference_genome_plan(genome_id, catalog_path, cache_dir)
            }
            GenomePrepareLaunchMode::ReindexCachedFiles => {
                GentleEngine::reindex_reference_genome_plan(genome_id, catalog_path, cache_dir)
            }
            GenomePrepareLaunchMode::RefreshFromSources => {
                GentleEngine::reinstall_reference_genome_plan(genome_id, catalog_path, cache_dir)
            }
        }
    }

    fn reset_prepare_step_state_from_plan(&mut self, plan: Option<PrepareGenomePlan>) {
        self.genome_prepare_steps = plan
            .map(|plan| {
                plan.steps
                    .into_iter()
                    .map(PrepareGenomeUiStepState::from_plan_step)
                    .collect()
            })
            .unwrap_or_default();
    }

    fn reset_prepare_dialog_preview_state(&mut self) {
        if self.genome_prepare_task.is_none() {
            self.genome_prepare_steps.clear();
            self.genome_prepare_progress = None;
            self.genome_prepare_eta_baseline = None;
            self.genome_prepare_failure_recovery = None;
            self.genome_prepare_status.clear();
        }
    }

    fn clear_prepare_step_state(&mut self) {
        if self.genome_prepare_task.is_none() {
            self.genome_prepare_steps.clear();
            self.genome_prepare_progress = None;
            self.genome_prepare_eta_baseline = None;
        }
    }

    fn selected_prepared_genome_inspection(
        &self,
    ) -> Result<Option<PreparedGenomeInspection>, String> {
        let genome_id = self.genome_id.trim();
        if genome_id.is_empty() {
            return Ok(None);
        }
        let catalog_path = self.genome_catalog_path_resolved();
        let catalog = GenomeCatalog::from_json_file(&catalog_path)?;
        let cache_dir = self.genome_cache_dir_opt();
        catalog.inspect_prepared_genome(genome_id, cache_dir.as_deref())
    }

    fn apply_prepared_genome_inspection_to_steps(
        &mut self,
        inspection: Option<&PreparedGenomeInspection>,
    ) {
        for step in &mut self.genome_prepare_steps {
            step.status = PrepareGenomeUiStepStatus::Pending;
            step.progress_fraction = None;
            step.raw_phase = None;
            step.bytes_done = None;
            step.bytes_total = None;
            step.eta_remaining = None;

            let (completed, detail) = match (step.step_id, inspection) {
                (PrepareGenomeStepId::ResetIndexes, Some(_)) => (
                    false,
                    Some("Will clear derived indexes before rebuilding".to_string()),
                ),
                (PrepareGenomeStepId::ResetIndexes, None) => (false, None),
                (PrepareGenomeStepId::Sequence, Some(inspection)) => (
                    inspection.sequence_present,
                    Some(inspection.sequence_path.clone()),
                ),
                (PrepareGenomeStepId::Annotation, Some(inspection)) => (
                    inspection.annotation_present,
                    Some(inspection.annotation_path.clone()),
                ),
                (PrepareGenomeStepId::FastaIndex, Some(inspection)) => (
                    inspection.fasta_index_ready,
                    Some(inspection.fasta_index_path.clone()),
                ),
                (PrepareGenomeStepId::GeneIndex, Some(inspection)) => {
                    let transcript_required = inspection.transcript_index_path.is_some();
                    let completed = inspection.gene_index_ready
                        && (!transcript_required || inspection.transcript_index_ready);
                    let detail = inspection
                        .transcript_index_path
                        .clone()
                        .unwrap_or_else(|| inspection.gene_index_path.clone());
                    (completed, Some(detail))
                }
                (PrepareGenomeStepId::BlastIndex, Some(inspection)) => (
                    inspection.blast_index_ready,
                    inspection.blast_db_prefix.clone(),
                ),
                (_, None) => (false, None),
            };

            if completed {
                step.status = PrepareGenomeUiStepStatus::Completed;
                step.progress_fraction = Some(1.0);
            }
            step.detail = detail.unwrap_or_default();
        }
    }

    fn ensure_prepare_step_preview_for_current_selection(
        &mut self,
        primary_action: PrepareGenomeDialogPrimaryAction,
    ) {
        if self.genome_prepare_task.is_some()
            || !self.genome_prepare_steps.is_empty()
            || self.genome_prepare_progress.is_some()
        {
            return;
        }
        let genome_id = self.genome_id.trim();
        if genome_id.is_empty() || matches!(primary_action, PrepareGenomeDialogPrimaryAction::None)
        {
            return;
        }
        let mode = match primary_action {
            PrepareGenomeDialogPrimaryAction::Prepare => GenomePrepareLaunchMode::Prepare,
            PrepareGenomeDialogPrimaryAction::Reindex => {
                GenomePrepareLaunchMode::ReindexCachedFiles
            }
            PrepareGenomeDialogPrimaryAction::None => return,
        };
        let catalog_path = self.genome_catalog_path_opt();
        let cache_dir = self.genome_cache_dir_opt();
        let plan = self
            .prepare_plan_for_mode(
                genome_id,
                catalog_path.as_deref(),
                cache_dir.as_deref(),
                mode,
            )
            .ok();
        self.reset_prepare_step_state_from_plan(plan);
        let inspection = self.selected_prepared_genome_inspection().ok().flatten();
        self.apply_prepared_genome_inspection_to_steps(inspection.as_ref());
    }

    fn format_duration_compact(duration: Duration) -> String {
        let total_secs = duration.as_secs();
        if total_secs >= 3600 {
            let hours = total_secs / 3600;
            let minutes = (total_secs % 3600) / 60;
            if minutes > 0 {
                format!("{hours}h {minutes}m")
            } else {
                format!("{hours}h")
            }
        } else if total_secs >= 60 {
            let minutes = total_secs / 60;
            let seconds = total_secs % 60;
            if seconds > 0 {
                format!("{minutes}m {seconds}s")
            } else {
                format!("{minutes}m")
            }
        } else if total_secs > 0 {
            format!("{total_secs}s")
        } else {
            "<1s".to_string()
        }
    }

    fn format_prepare_byte_progress(
        bytes_done: u64,
        bytes_total: u64,
        eta_remaining: Option<Duration>,
    ) -> String {
        let mut text = format!("bytes: {} / {}", bytes_done, bytes_total);
        if let Some(eta) = eta_remaining {
            text.push_str(&format!(" • ETA {}", Self::format_duration_compact(eta)));
        }
        text
    }

    fn update_prepare_eta_from_progress_at(
        &mut self,
        progress: &PrepareGenomeProgress,
        observed_at: Instant,
    ) -> Option<Duration> {
        let Some(step_id) = progress.step_id else {
            self.genome_prepare_eta_baseline = None;
            return None;
        };
        let Some(bytes_total) = progress.bytes_total.filter(|total| *total > 0) else {
            self.genome_prepare_eta_baseline = None;
            return None;
        };
        if progress.phase == "ready" {
            self.genome_prepare_eta_baseline = None;
            return None;
        }

        let bytes_done = progress.bytes_done.min(bytes_total);
        let should_reset = self
            .genome_prepare_eta_baseline
            .as_ref()
            .map(|baseline| {
                baseline.step_id != step_id
                    || baseline.bytes_total != bytes_total
                    || bytes_done < baseline.bytes_done
            })
            .unwrap_or(true);
        if should_reset {
            self.genome_prepare_eta_baseline = Some(PrepareGenomeEtaBaseline {
                step_id,
                observed_at,
                bytes_done,
                bytes_total,
            });
            return None;
        }

        let Some(baseline) = self.genome_prepare_eta_baseline.as_ref() else {
            return None;
        };
        if bytes_done <= baseline.bytes_done {
            return None;
        }
        let elapsed = observed_at.saturating_duration_since(baseline.observed_at);
        if elapsed.is_zero() {
            return None;
        }
        let advanced = bytes_done.saturating_sub(baseline.bytes_done);
        if advanced == 0 {
            return None;
        }
        let remaining = bytes_total.saturating_sub(bytes_done);
        if remaining == 0 {
            return None;
        }
        let throughput = advanced as f64 / elapsed.as_secs_f64();
        if !throughput.is_finite() || throughput <= 0.0 {
            return None;
        }
        let eta_secs = remaining as f64 / throughput;
        if !eta_secs.is_finite() || eta_secs <= 0.0 {
            return None;
        }
        Some(Duration::from_secs_f64(eta_secs))
    }

    fn is_prepare_reinstall_recommended_failure(
        message: &str,
        mode: GenomePrepareLaunchMode,
    ) -> bool {
        if mode != GenomePrepareLaunchMode::ReindexCachedFiles {
            return false;
        }
        let lower = message.to_ascii_lowercase();
        lower.contains("prepared genome")
            && lower.contains("is inconsistent")
            && lower.contains("references contigs missing from prepared sequence")
    }

    fn queue_prepare_failure_reinstall(&mut self, dialog_host: PreparedGenomeReinstallDialogHost) {
        let Some(recovery) = self.genome_prepare_failure_recovery.clone() else {
            return;
        };
        self.pending_prepared_genome_reinstall = Some(PreparedGenomeReinstallRequest {
            genome_id: recovery.genome_id,
            scope: recovery.scope,
            catalog_path: recovery.catalog_path,
            cache_dir: recovery.cache_dir,
            dialog_host,
        });
    }

    fn prepare_progress_fraction(progress: &PrepareGenomeProgress) -> Option<f32> {
        progress
            .percent
            .map(|p| (p as f32 / 100.0).clamp(0.0, 1.0))
            .or_else(|| {
                progress.bytes_total.and_then(|total| {
                    if total == 0 {
                        None
                    } else {
                        Some((progress.bytes_done as f32 / total as f32).clamp(0.0, 1.0))
                    }
                })
            })
    }

    fn apply_prepare_progress_to_steps(&mut self, progress: &PrepareGenomeProgress) {
        self.apply_prepare_progress_to_steps_at(progress, Instant::now());
    }

    fn apply_prepare_progress_to_steps_at(
        &mut self,
        progress: &PrepareGenomeProgress,
        observed_at: Instant,
    ) {
        let Some(step_id) = progress.step_id else {
            self.genome_prepare_eta_baseline = None;
            return;
        };
        let Some(current_index) = self
            .genome_prepare_steps
            .iter()
            .position(|step| step.step_id == step_id)
        else {
            return;
        };

        for step in self.genome_prepare_steps.iter_mut().take(current_index) {
            if matches!(
                step.status,
                PrepareGenomeUiStepStatus::Pending | PrepareGenomeUiStepStatus::Running
            ) {
                step.status = PrepareGenomeUiStepStatus::Completed;
                step.progress_fraction = Some(1.0);
                step.eta_remaining = None;
                if step.detail.is_empty() {
                    step.detail = step.operation_summary.clone();
                }
            }
        }

        let eta_remaining = self.update_prepare_eta_from_progress_at(progress, observed_at);
        let step = &mut self.genome_prepare_steps[current_index];
        step.raw_phase = Some(progress.phase.clone());
        step.detail = progress.item.clone();
        step.progress_fraction = Self::prepare_progress_fraction(progress);
        step.bytes_done = Some(progress.bytes_done);
        step.bytes_total = progress.bytes_total;
        step.eta_remaining = eta_remaining;

        let completed = progress.phase == "ready"
            || Self::prepare_progress_fraction(progress)
                .map(|fraction| (fraction - 1.0).abs() < f32::EPSILON)
                .unwrap_or(false);
        if completed {
            step.status = PrepareGenomeUiStepStatus::Completed;
            step.progress_fraction = Some(1.0);
            step.eta_remaining = None;
        } else {
            step.status = PrepareGenomeUiStepStatus::Running;
        }
    }

    fn finalize_prepare_steps_success(&mut self) {
        for step in &mut self.genome_prepare_steps {
            step.status = PrepareGenomeUiStepStatus::Completed;
            step.progress_fraction = Some(1.0);
            step.eta_remaining = None;
            if step.detail.is_empty() {
                step.detail = step.operation_summary.clone();
            }
        }
        self.genome_prepare_eta_baseline = None;
    }

    fn finalize_prepare_steps_failure(&mut self, cancelled: bool) {
        if let Some(step) = self
            .genome_prepare_steps
            .iter_mut()
            .rev()
            .find(|step| step.status == PrepareGenomeUiStepStatus::Running)
        {
            step.status = if cancelled {
                PrepareGenomeUiStepStatus::Cancelled
            } else {
                PrepareGenomeUiStepStatus::Failed
            };
            step.eta_remaining = None;
            if step.progress_fraction.is_none() {
                step.progress_fraction = Some(0.0);
            }
        }
        self.genome_prepare_eta_baseline = None;
    }

    fn start_prepare_reference_genome(&mut self) {
        self.start_prepare_reference_genome_with_mode(GenomePrepareLaunchMode::Prepare);
    }

    fn start_prepare_reference_genome_for_current_selection(&mut self) {
        let preparable_genomes = match self.unprepared_genomes_for_prepare_dialog() {
            Ok(names) => names,
            Err(e) => {
                self.genome_prepare_status = format!("Prepared-state check error: {e}");
                return;
            }
        };
        let preparable_set: HashSet<String> = preparable_genomes.into_iter().collect();
        match Self::prepare_dialog_primary_action(
            &self.genome_id,
            &self.genome_catalog_genomes,
            &preparable_set,
        ) {
            PrepareGenomeDialogPrimaryAction::Prepare => {
                self.start_prepare_reference_genome_with_mode(GenomePrepareLaunchMode::Prepare);
            }
            PrepareGenomeDialogPrimaryAction::Reindex => {
                self.start_prepare_reference_genome_with_mode(
                    GenomePrepareLaunchMode::ReindexCachedFiles,
                );
            }
            PrepareGenomeDialogPrimaryAction::None => {
                self.genome_prepare_status = "Select a genome first".to_string();
            }
        }
    }

    fn start_prepare_reference_genome_with_mode(&mut self, mode: GenomePrepareLaunchMode) {
        if self.genome_prepare_task.is_some() {
            self.genome_prepare_status = "Genome preparation is already running".to_string();
            return;
        }
        let scope = self.genome_dialog_scope;
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
        let binary_preflight = self
            .engine
            .read()
            .unwrap()
            .blast_external_binary_preflight_report();
        let makeblastdb_preflight =
            Self::summarize_binary_probe(&binary_preflight.makeblastdb, true);
        let catalog_path = self.genome_catalog_path_opt();
        let cache_dir = self.genome_cache_dir_opt();
        let prepare_plan = self
            .prepare_plan_for_mode(
                &genome_id,
                catalog_path.as_deref(),
                cache_dir.as_deref(),
                mode,
            )
            .ok();
        let cancel_requested = Arc::new(AtomicBool::new(false));
        let job_id = self.alloc_background_job_id();
        let (tx, rx) = mpsc::channel::<GenomePrepareTaskMessage>();
        let action_title = mode.progress_label();
        let action_summary = mode.job_summary();
        self.reset_prepare_step_state_from_plan(prepare_plan);
        self.genome_prepare_eta_baseline = None;
        self.genome_prepare_failure_recovery = None;
        self.genome_prepare_progress = Some(PrepareGenomeProgress {
            genome_id: genome_id.clone(),
            phase: "queued".to_string(),
            item: "waiting for worker".to_string(),
            bytes_done: 0,
            bytes_total: None,
            percent: Some(0.0),
            step_id: None,
            step_label: None,
        });
        self.genome_prepare_status = if let Some(timeout) = timeout_seconds {
            format!(
                "{action_title} genome '{genome_id}' in background (timeout: {} s). You can keep using the UI.\npreflight makeblastdb: {}",
                timeout, makeblastdb_preflight
            )
        } else {
            format!(
                "{action_title} genome '{genome_id}' in background. You can keep using the UI.\npreflight makeblastdb: {}",
                makeblastdb_preflight
            )
        };
        self.push_job_event(
            BackgroundJobKind::PrepareGenome,
            BackgroundJobEventPhase::Started,
            Some(job_id),
            format!(
                "{action_summary}: {} (makeblastdb preflight: {})",
                genome_id, makeblastdb_preflight
            ),
        );
        self.genome_prepare_task = Some(GenomePrepareTask {
            job_id,
            started: Instant::now(),
            cancel_requested: cancel_requested.clone(),
            timeout_seconds,
            mode,
            genome_id: genome_id.clone(),
            scope,
            catalog_path: self.genome_catalog_path.clone(),
            cache_dir: self.genome_cache_dir.clone(),
            receiver: rx,
        });
        std::thread::spawn(move || {
            let tx_progress = tx.clone();
            let cancel_flag = cancel_requested.clone();
            let mut last_phase = String::new();
            let mut last_percent_tenths: Option<i64> = None;
            let mut last_bytes_bucket: u64 = 0;
            let mut progress_forwarder = move |p: PrepareGenomeProgress| -> bool {
                if cancel_flag.load(Ordering::Relaxed) {
                    return false;
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
            let outcome = (match mode {
                GenomePrepareLaunchMode::Prepare => GentleEngine::prepare_reference_genome_once(
                    &genome_id,
                    catalog_path.as_deref(),
                    cache_dir.as_deref(),
                    timeout_seconds,
                    &mut progress_forwarder,
                ),
                GenomePrepareLaunchMode::ReindexCachedFiles => {
                    GentleEngine::reindex_reference_genome_once(
                        &genome_id,
                        catalog_path.as_deref(),
                        cache_dir.as_deref(),
                        timeout_seconds,
                        &mut progress_forwarder,
                    )
                }
                GenomePrepareLaunchMode::RefreshFromSources => {
                    GentleEngine::reinstall_reference_genome_once(
                        &genome_id,
                        catalog_path.as_deref(),
                        cache_dir.as_deref(),
                        timeout_seconds,
                        &mut progress_forwarder,
                    )
                }
            })
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
                protocol_cartoon_preview: None,
                genome_annotation_projection: None,
                sequence_alignment: None,
                protein_derivation_report: None,
                reverse_translation_report: None,
                protease_digest_report: None,
                protein_residue_genomic_coordinates: None,
                exon_skip_selection_plan: None,
                exon_skip_materialization: None,
                cdna_assay_test_report: None,
                cdna_assay_product_materialization: None,
                primer_specificity_report: None,
                transcript_qpcr_panel: None,
                construct_reasoning_graph: None,
                sequencing_confirmation_report: None,
                sequencing_primer_overlay_report: None,
                sequencing_trace_import_report: None,
                sequencing_trace_record: None,
                sequencing_trace_summaries: None,
                cutrun_dataset_list: None,
                cutrun_dataset_status: None,
                cutrun_read_report: None,
                cutrun_read_report_summaries: None,
                cutrun_read_coverage_export: None,
                cutrun_regulatory_support: None,
                read_acquisition_report: None,
                cutrun_dataset_projection: None,
                microarray_projection: None,
                genome_coordinate_projection: None,
                rna_read_gene_support_summary: None,
                rna_read_gene_support_audit: None,
                rna_read_target_quality_export: None,
                rna_read_batch_map_report: None,
                rna_read_isoform_preflight: None,
                tfbs_region_summary: None,
                tfbs_score_tracks: None,
                tfbs_track_similarity: None,
                multi_gene_promoter_tfbs: None,
                repeat_annotation_query: None,
                sequence_repeat_overlaps: None,
                repeat_feature_materialization: None,
                repeat_environment_cohort: None,
                window_cohort_tfbs: None,
                tfbs_hit_scan: None,
                restriction_site_scan: None,
                jaspar_remote_metadata_snapshot: None,
                jaspar_catalog_report: None,
                tf_query_resolution_report: None,
                jaspar_entry_expert_view: None,
                jaspar_registry_benchmark: None,
                jaspar_entry_presentation: None,
                sequence_context_view: None,
                sequence_context_bundle: None,
                alternative_promoter_comparison: None,
                variant_promoter_context: None,
                promoter_evidence_matrix: None,
                isoform_promoter_comparison: None,
                promoter_expression_evidence: None,
                promoter_artifact_manifest: None,
                promoter_reporter_candidates: None,
                reporter_catalog: None,
                reporter_recommendation: None,
                reporter_corpus_export: None,
                reporter_construct_handoff: None,
                uniprot_projection_audit: None,
                uniprot_projection_audit_parity: None,
                lab_assistant_instructions: None,
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
        let mut progress_messages: Vec<(u64, PrepareGenomeProgress, bool, &'static str)> = vec![];
        if let Some(task) = &self.genome_prepare_task {
            let active_job_id = task.job_id;
            let action_title = task.mode.progress_label();
            let cancel_requested = task.cancel_requested.clone();
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
                        progress_messages.push((
                            job_id,
                            progress,
                            cancel_requested.load(Ordering::Relaxed),
                            action_title,
                        ));
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

                                cause_chain: vec![],
                            }),
                        ));
                        break;
                    }
                }
            }
        }
        for (_job_id, progress, canceling, action_title) in progress_messages {
            self.genome_prepare_progress = Some(progress.clone());
            self.apply_prepare_progress_to_steps(&progress);
            self.genome_prepare_status = format!(
                "{action_title} genome '{}': {} ({}){}",
                progress.genome_id,
                progress.phase,
                progress.item,
                if canceling {
                    " (cancellation requested)"
                } else {
                    ""
                }
            );
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
            let (elapsed, cancellation_requested, timeout_seconds, mode, failure_recovery) = self
                .genome_prepare_task
                .as_ref()
                .map(|task| {
                    (
                        task.started.elapsed().as_secs_f64(),
                        task.cancel_requested.load(Ordering::Relaxed),
                        task.timeout_seconds,
                        task.mode,
                        PrepareGenomeFailureRecovery {
                            genome_id: task.genome_id.clone(),
                            scope: task.scope,
                            catalog_path: task.catalog_path.clone(),
                            cache_dir: task.cache_dir.clone(),
                        },
                    )
                })
                .unwrap_or((
                    0.0,
                    false,
                    None,
                    GenomePrepareLaunchMode::Prepare,
                    PrepareGenomeFailureRecovery {
                        genome_id: self.genome_id.clone(),
                        scope: self.genome_dialog_scope,
                        catalog_path: self.genome_catalog_path.clone(),
                        cache_dir: self.genome_cache_dir.clone(),
                    },
                ));
            self.genome_prepare_task = None;
            self.genome_prepare_eta_baseline = None;
            match outcome {
                Ok(result) => {
                    self.finalize_prepare_steps_success();
                    self.genome_prepare_failure_recovery = None;
                    let prefix = mode.result_prefix(cancellation_requested);
                    self.genome_prepare_status = format!(
                        "{}\nelapsed: {:.1}s",
                        Self::format_op_result_status(
                            prefix,
                            &result.created_seq_ids,
                            &result.warnings,
                            &result.messages
                        ),
                        elapsed
                    );
                    self.invalidate_genome_genes();
                    self.push_job_event(
                        BackgroundJobKind::PrepareGenome,
                        BackgroundJobEventPhase::Completed,
                        Some(job_id),
                        format!("{prefix} in {:.1}s", elapsed),
                    );
                }
                Err(e) => {
                    self.finalize_prepare_steps_failure(cancellation_requested);
                    self.genome_prepare_failure_recovery =
                        if Self::is_prepare_reinstall_recommended_failure(&e.message, mode) {
                            Some(failure_recovery)
                        } else {
                            None
                        };
                    self.genome_prepare_status = Self::format_prepare_failure_status(
                        &e.message,
                        elapsed,
                        cancellation_requested,
                        timeout_seconds,
                        mode,
                    );
                    self.push_job_event(
                        BackgroundJobKind::PrepareGenome,
                        BackgroundJobEventPhase::Failed,
                        Some(job_id),
                        format!(
                            "{} genome ended in {:.1}s: {}",
                            match mode {
                                GenomePrepareLaunchMode::Prepare => "Prepare",
                                GenomePrepareLaunchMode::ReindexCachedFiles => "Reindex",
                                GenomePrepareLaunchMode::RefreshFromSources => "Refresh",
                            },
                            elapsed,
                            e.message
                        ),
                    );
                }
            }
        }
    }

    fn poll_tutorial_project_task(&mut self, ctx: &egui::Context) {
        if self.tutorial_project_task.is_none() {
            return;
        }
        ctx.request_repaint_after(Duration::from_millis(100));
        let mut done: Option<(u64, Result<TutorialProjectOpenOutcome, String>)> = None;
        let mut stale_job_ids: Vec<u64> = vec![];
        let mut progress_messages: Vec<(u64, TutorialProjectTaskProgress, bool)> = vec![];
        if let Some(task) = &self.tutorial_project_task {
            let active_job_id = task.job_id;
            let cancel_requested = task.cancel_requested.clone();
            const MAX_MESSAGES_PER_TICK: usize = 128;
            for _ in 0..MAX_MESSAGES_PER_TICK {
                match task.receiver.try_recv() {
                    Ok(TutorialProjectTaskMessage::Progress { job_id, progress }) => {
                        if job_id != active_job_id {
                            if !stale_job_ids.contains(&job_id) {
                                stale_job_ids.push(job_id);
                            }
                            continue;
                        }
                        progress_messages.push((
                            job_id,
                            progress,
                            cancel_requested.load(Ordering::Relaxed),
                        ));
                    }
                    Ok(TutorialProjectTaskMessage::Done { job_id, result }) => {
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
                            Err("Tutorial project worker disconnected".to_string()),
                        ));
                        break;
                    }
                }
            }
        }

        for (_job_id, progress, canceling) in progress_messages {
            self.tutorial_project_progress_fraction = progress.percent;
            self.tutorial_project_progress_label = progress.item.clone();
            self.tutorial_project_status = Self::tutorial_project_status_from_progress(&progress);
            if canceling {
                self.tutorial_project_status
                    .push_str(" (cancellation requested)");
            }
            self.app_status = self.tutorial_project_status.clone();
            if let Some(task) = self.tutorial_project_task.as_mut() {
                task.chapter_id = progress.chapter_id.clone();
                task.chapter_title = progress.chapter_title.clone();
            }
        }
        for stale_job_id in stale_job_ids {
            self.push_job_event(
                BackgroundJobKind::OpenTutorialProject,
                BackgroundJobEventPhase::IgnoredStale,
                Some(stale_job_id),
                "Ignored stale tutorial-project worker message",
            );
        }

        if let Some((job_id, outcome)) = done {
            let (elapsed, cancellation_requested, chapter_id, chapter_title) = self
                .tutorial_project_task
                .as_ref()
                .map(|task| {
                    (
                        task.started.elapsed().as_secs_f64(),
                        task.cancel_requested.load(Ordering::Relaxed),
                        task.chapter_id.clone(),
                        task.chapter_title.clone(),
                    )
                })
                .unwrap_or((0.0, false, String::new(), String::new()));
            self.tutorial_project_task = None;
            match outcome {
                Ok(result) => {
                    if cancellation_requested {
                        self.tutorial_project_progress_fraction = None;
                        self.tutorial_project_progress_label.clear();
                        self.tutorial_project_status = format!(
                            "Tutorial opening cancelled after {:.1}s. Temporary project kept at '{}'.",
                            elapsed, result.project_path
                        );
                        self.app_status = self.tutorial_project_status.clone();
                        self.push_job_event(
                            BackgroundJobKind::OpenTutorialProject,
                            BackgroundJobEventPhase::Completed,
                            Some(job_id),
                            format!(
                                "Tutorial '{}' cancelled after {:.1}s",
                                result.chapter_id, elapsed
                            ),
                        );
                        return;
                    }
                    if let Err(err) =
                        self.load_project_from_file_with_recent(&result.project_path, false)
                    {
                        self.tutorial_project_progress_fraction = None;
                        self.tutorial_project_progress_label.clear();
                        self.tutorial_project_status = format!(
                            "Could not open tutorial project '{}': {err}",
                            result.project_path
                        );
                        self.app_status = self.tutorial_project_status.clone();
                        self.push_job_event(
                            BackgroundJobKind::OpenTutorialProject,
                            BackgroundJobEventPhase::Failed,
                            Some(job_id),
                            format!(
                                "Tutorial '{}' build finished in {:.1}s but project open failed: {}",
                                result.chapter_id, elapsed, err
                            ),
                        );
                        return;
                    }
                    self.tutorial_project_progress_fraction = None;
                    self.tutorial_project_progress_label.clear();
                    self.tutorial_project_status = format!(
                        "Opened tutorial project: {} ({})",
                        result.chapter_title, result.chapter_id
                    );
                    if let Err(err) = self.open_help_tutorial_path(
                        &result.guide_path,
                        &result.chapter_title,
                        &result.guide_summary,
                    ) {
                        self.tutorial_project_status
                            .push_str(&format!(" | guide unavailable: {err}"));
                    }
                    self.app_status = self.tutorial_project_status.clone();
                    self.push_job_event(
                        BackgroundJobKind::OpenTutorialProject,
                        BackgroundJobEventPhase::Completed,
                        Some(job_id),
                        format!("Opened tutorial '{}' in {:.1}s", result.chapter_id, elapsed),
                    );
                }
                Err(err) => {
                    self.tutorial_project_progress_fraction = None;
                    self.tutorial_project_progress_label.clear();
                    let lower = err.to_ascii_lowercase();
                    let cancelled = cancellation_requested
                        || lower.contains("cancelled")
                        || lower.contains("canceled");
                    self.tutorial_project_status = if cancelled {
                        format!("Tutorial opening cancelled after {:.1}s: {}", elapsed, err)
                    } else if chapter_title.trim().is_empty() {
                        format!("Could not build tutorial project: {err}")
                    } else {
                        format!(
                            "Could not build tutorial project '{}' after {:.1}s: {}",
                            chapter_title, elapsed, err
                        )
                    };
                    self.app_status = self.tutorial_project_status.clone();
                    self.push_job_event(
                        BackgroundJobKind::OpenTutorialProject,
                        BackgroundJobEventPhase::Failed,
                        Some(job_id),
                        if cancelled {
                            format!(
                                "Tutorial '{}' cancelled after {:.1}s",
                                if chapter_id.trim().is_empty() {
                                    "-"
                                } else {
                                    chapter_id.as_str()
                                },
                                elapsed
                            )
                        } else {
                            format!(
                                "Tutorial '{}' failed after {:.1}s: {}",
                                if chapter_id.trim().is_empty() {
                                    "-"
                                } else {
                                    chapter_id.as_str()
                                },
                                elapsed,
                                err
                            )
                        },
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

                                cause_chain: vec![],
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
                    let refreshed_windows = if result.changed_seq_ids.is_empty() {
                        self.refresh_sequence_windows_from_engine_state()
                    } else {
                        self.refresh_sequence_windows_for_seq_ids(&result.changed_seq_ids)
                    };
                    if refreshed_windows > 0 {
                        self.genome_track_status.push_str(&format!(
                            "\nrefreshed {} open sequence window(s)",
                            refreshed_windows
                        ));
                        ctx.request_repaint();
                    }
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

    fn merge_json_objects_recursive(
        base: &mut serde_json::Map<String, serde_json::Value>,
        overlay: &serde_json::Map<String, serde_json::Value>,
    ) {
        for (key, value) in overlay {
            match (base.get_mut(key), value) {
                (
                    Some(serde_json::Value::Object(existing)),
                    serde_json::Value::Object(incoming),
                ) => {
                    Self::merge_json_objects_recursive(existing, incoming);
                }
                _ => {
                    base.insert(key.clone(), value.clone());
                }
            }
        }
    }

    fn build_genome_blast_request_override_json(
        &self,
    ) -> Result<Option<serde_json::Value>, String> {
        let mut merged = serde_json::Map::<String, serde_json::Value>::new();
        if let Some(preset) = self.genome_blast_options_preset.as_request_override_json() {
            let Some(obj) = preset.as_object() else {
                return Err(
                    "Internal preset error: BLAST preset payload is not a JSON object".to_string(),
                );
            };
            Self::merge_json_objects_recursive(&mut merged, obj);
        }
        if let Some(structured_thresholds) =
            self.build_structured_blast_thresholds_override_json()?
        {
            let Some(obj) = structured_thresholds.as_object() else {
                return Err("Internal threshold error: structured BLAST thresholds payload is not a JSON object".to_string());
            };
            Self::merge_json_objects_recursive(&mut merged, obj);
        }
        let raw_advanced = self.genome_blast_options_json.trim();
        if !raw_advanced.is_empty() {
            let parsed: serde_json::Value = serde_json::from_str(raw_advanced)
                .map_err(|e| format!("Invalid advanced BLAST options JSON: {e}"))?;
            let Some(obj) = parsed.as_object() else {
                return Err("Advanced BLAST options must decode to a JSON object".to_string());
            };
            Self::merge_json_objects_recursive(&mut merged, obj);
        }
        if merged.is_empty() {
            Ok(None)
        } else {
            Ok(Some(serde_json::Value::Object(merged)))
        }
    }

    fn build_structured_blast_thresholds_override_json(
        &self,
    ) -> Result<Option<serde_json::Value>, String> {
        let mut thresholds = serde_json::Map::<String, serde_json::Value>::new();

        if self.genome_blast_threshold_use_max_evalue {
            let parsed = self
                .genome_blast_threshold_max_evalue
                .trim()
                .parse::<f64>()
                .map_err(|e| format!("Invalid max_evalue value: {e}"))?;
            if !parsed.is_finite() || parsed < 0.0 {
                return Err("max_evalue must be a finite number >= 0.0".to_string());
            }
            thresholds.insert("max_evalue".to_string(), serde_json::Value::from(parsed));
        }

        if self.genome_blast_threshold_use_min_identity_percent {
            let parsed = self
                .genome_blast_threshold_min_identity_percent
                .trim()
                .parse::<f64>()
                .map_err(|e| format!("Invalid min_identity_percent value: {e}"))?;
            if !parsed.is_finite() || !(0.0..=100.0).contains(&parsed) {
                return Err("min_identity_percent must be within [0, 100]".to_string());
            }
            thresholds.insert(
                "min_identity_percent".to_string(),
                serde_json::Value::from(parsed),
            );
        }

        if self.genome_blast_threshold_use_min_query_coverage_percent {
            let parsed = self
                .genome_blast_threshold_min_query_coverage_percent
                .trim()
                .parse::<f64>()
                .map_err(|e| format!("Invalid min_query_coverage_percent value: {e}"))?;
            if !parsed.is_finite() || !(0.0..=100.0).contains(&parsed) {
                return Err("min_query_coverage_percent must be within [0, 100]".to_string());
            }
            thresholds.insert(
                "min_query_coverage_percent".to_string(),
                serde_json::Value::from(parsed),
            );
        }

        if self.genome_blast_threshold_use_min_alignment_length_bp {
            let parsed = self
                .genome_blast_threshold_min_alignment_length_bp
                .trim()
                .parse::<usize>()
                .map_err(|e| format!("Invalid min_alignment_length_bp value: {e}"))?;
            if parsed == 0 {
                return Err("min_alignment_length_bp must be >= 1".to_string());
            }
            thresholds.insert(
                "min_alignment_length_bp".to_string(),
                serde_json::Value::from(parsed),
            );
        }

        if self.genome_blast_threshold_use_min_bit_score {
            let parsed = self
                .genome_blast_threshold_min_bit_score
                .trim()
                .parse::<f64>()
                .map_err(|e| format!("Invalid min_bit_score value: {e}"))?;
            if !parsed.is_finite() || parsed < 0.0 {
                return Err("min_bit_score must be a finite number >= 0.0".to_string());
            }
            thresholds.insert("min_bit_score".to_string(), serde_json::Value::from(parsed));
        }

        if self.genome_blast_threshold_unique_best_hit {
            thresholds.insert("unique_best_hit".to_string(), serde_json::Value::from(true));
        }

        if thresholds.is_empty() {
            Ok(None)
        } else {
            Ok(Some(serde_json::json!({ "thresholds": thresholds })))
        }
    }

    fn resolve_genome_blast_options_preview(
        &self,
        request_override_json: Option<&serde_json::Value>,
    ) -> Result<crate::engine::BlastResolvedOptions, String> {
        let task_name = self.genome_blast_task_name.trim();
        let legacy_task = if task_name.is_empty() {
            None
        } else {
            Some(task_name)
        };
        let legacy_max_hits = Some(self.genome_blast_max_hits.max(1));
        let engine = self.engine.read().unwrap();
        engine
            .resolve_blast_options_for_request(request_override_json, legacy_task, legacy_max_hits)
            .map_err(|e| e.message)
    }

    fn format_blast_thresholds_summary(
        thresholds: &crate::engine::BlastThresholdOptions,
    ) -> String {
        let mut parts: Vec<String> = vec![];
        if let Some(v) = thresholds.max_evalue {
            parts.push(format!("max_evalue={v:.3e}"));
        }
        if let Some(v) = thresholds.min_identity_percent {
            parts.push(format!("min_identity={v:.2}%"));
        }
        if let Some(v) = thresholds.min_query_coverage_percent {
            parts.push(format!("min_qcov={v:.2}%"));
        }
        if let Some(v) = thresholds.min_alignment_length_bp {
            parts.push(format!("min_aln_bp={v}"));
        }
        if let Some(v) = thresholds.min_bit_score {
            parts.push(format!("min_bit_score={v:.2}"));
        }
        if let Some(v) = thresholds.unique_best_hit {
            parts.push(format!("unique_best_hit={v}"));
        }
        if parts.is_empty() {
            "none".to_string()
        } else {
            parts.join(", ")
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
        let legacy_max_hits = self.genome_blast_max_hits.max(1);
        let blast_task_name = self.genome_blast_task_name.trim().to_string();
        let legacy_task_arg = if blast_task_name.is_empty() {
            None
        } else {
            Some(blast_task_name)
        };
        let request_override_json = match self.build_genome_blast_request_override_json() {
            Ok(v) => v,
            Err(e) => {
                self.genome_blast_status = e;
                return;
            }
        };
        let resolved_preview =
            match self.resolve_genome_blast_options_preview(request_override_json.as_ref()) {
                Ok(v) => v,
                Err(e) => {
                    self.genome_blast_status = format!("BLAST options preflight failed: {e}");
                    return;
                }
            };
        let binary_preflight = self
            .engine
            .read()
            .unwrap()
            .blast_external_binary_preflight_report();
        let blastn_preflight = Self::summarize_binary_probe(&binary_preflight.blastn, true);
        let makeblastdb_preflight =
            Self::summarize_binary_probe(&binary_preflight.makeblastdb, false);
        let blastn_bin_for_status = binary_preflight.blastn.executable.clone();
        let project_state_snapshot = self.engine.read().unwrap().state().clone();
        let job_id = self.alloc_background_job_id();
        let (tx, rx) = mpsc::channel::<GenomeBlastTaskMessage>();
        let cancel_requested = Arc::new(AtomicBool::new(false));
        self.genome_blast_results.clear();
        self.genome_blast_selected_result = 0;
        self.genome_blast_progress_fraction = Some(0.0);
        self.genome_blast_progress_label = format!("0 / {total_queries}");
        self.genome_blast_status = format!(
            "Running BLAST for {} quer{} in background\npreflight blastn: {}\npreflight makeblastdb: {}",
            total_queries,
            if total_queries == 1 { "y" } else { "ies" },
            blastn_preflight,
            makeblastdb_preflight
        );
        self.push_job_event(
            BackgroundJobKind::BlastGenome,
            BackgroundJobEventPhase::Started,
            Some(job_id),
            format!(
                "BLAST started: genome='{}', queries={}, task='{}', max_hits={}, thresholds={}, request_override={}, blastn_preflight={}, makeblastdb_preflight={}",
                genome_id,
                total_queries,
                resolved_preview.task,
                resolved_preview.max_hits,
                Self::format_blast_thresholds_summary(&resolved_preview.thresholds),
                if request_override_json.is_some() {
                    "yes"
                } else {
                    "no"
                },
                blastn_preflight,
                makeblastdb_preflight
            ),
        );
        self.genome_blast_task = Some(GenomeBlastTask {
            job_id,
            started: Instant::now(),
            cancel_requested: cancel_requested.clone(),
            receiver: rx,
        });
        std::thread::spawn(move || {
            let mut reports: Vec<GenomeBlastQueryResult> = vec![];
            let mut failed_queries: Vec<String> = vec![];
            let mut cancelled = false;
            let cancel_flag = cancel_requested.clone();
            let task_name_for_status = resolved_preview.task.clone();
            let max_hits_for_status = resolved_preview.max_hits;

            for (idx, (label, query)) in blast_queries.into_iter().enumerate() {
                if cancel_flag.load(Ordering::Relaxed) {
                    let _ = tx.send(GenomeBlastTaskMessage::Status {
                        job_id,
                        status: format!("BLAST cancellation acknowledged before query '{}'", label),
                    });
                    cancelled = true;
                    break;
                }
                let _ = tx.send(GenomeBlastTaskMessage::Progress {
                    job_id,
                    done_queries: idx,
                    total_queries,
                    current_query_label: label.clone(),
                });

                let query_length = query.len();
                let invocation_template = format!(
                    "{} -db <prepared_db_prefix> -query <temp_query.fa> -task {} -max_target_seqs {}",
                    blastn_bin_for_status, task_name_for_status, max_hits_for_status
                );
                let _ = tx.send(GenomeBlastTaskMessage::Status {
                    job_id,
                    status: format!(
                        "BLAST query '{}' started ({} bp)\ninvocation template: {}",
                        label, query_length, invocation_template
                    ),
                });

                let query_started = Instant::now();
                let (query_tx, query_rx) =
                    mpsc::channel::<Result<GenomeBlastReport, EngineError>>();
                let catalog_for_query = catalog_path.clone();
                let cache_for_query = cache_dir.clone();
                let genome_for_query = genome_id.clone();
                let task_for_query = legacy_task_arg.clone();
                let request_override_for_query = request_override_json.clone();
                let project_state_for_query = project_state_snapshot.clone();
                let cancel_for_query = cancel_flag.clone();
                std::thread::spawn(move || {
                    let blast_engine = GentleEngine::from_state(project_state_for_query);
                    let mut should_cancel = || cancel_for_query.load(Ordering::Relaxed);
                    let outcome = blast_engine
                        .blast_reference_genome_with_project_and_request_options_and_cancel(
                            catalog_for_query.as_deref(),
                            &genome_for_query,
                            &query,
                            request_override_for_query.as_ref(),
                            task_for_query.as_deref(),
                            Some(legacy_max_hits),
                            cache_for_query.as_deref(),
                            &mut should_cancel,
                        );
                    let _ = query_tx.send(outcome);
                });

                let report_result = loop {
                    match query_rx.recv_timeout(Duration::from_millis(250)) {
                        Ok(outcome) => break outcome,
                        Err(mpsc::RecvTimeoutError::Timeout) => {
                            let cancel_suffix = if cancel_flag.load(Ordering::Relaxed) {
                                " (cancellation requested)"
                            } else {
                                ""
                            };
                            let _ = tx.send(GenomeBlastTaskMessage::Status {
                                job_id,
                                status: format!(
                                    "BLAST query '{}' running ({:.1}s){}\ninvocation template: {}",
                                    label,
                                    query_started.elapsed().as_secs_f32(),
                                    cancel_suffix,
                                    invocation_template
                                ),
                            });
                        }
                        Err(mpsc::RecvTimeoutError::Disconnected) => {
                            break Err(EngineError {
                                code: ErrorCode::Io,
                                message: "BLAST worker disconnected".to_string(),

                                cause_chain: vec![],
                            });
                        }
                    }
                };

                match report_result {
                    Ok(report) => {
                        let invocation = Self::format_blast_command_line(
                            &report.blastn_executable,
                            &report.command,
                        );
                        let _ = tx.send(GenomeBlastTaskMessage::Status {
                            job_id,
                            status: format!(
                                "BLAST query '{}' completed in {:.1}s (hits={})\ninvocation: {}",
                                label,
                                query_started.elapsed().as_secs_f32(),
                                report.hit_count,
                                invocation
                            ),
                        });
                        reports.push(GenomeBlastQueryResult {
                            query_label: label.clone(),
                            query_length,
                            report,
                        });
                    }
                    Err(e) => {
                        if crate::genomes::is_blast_cancelled_error(&e.message)
                            || cancel_flag.load(Ordering::Relaxed)
                        {
                            let _ = tx.send(GenomeBlastTaskMessage::Status {
                                job_id,
                                status: format!(
                                    "BLAST cancelled while running query '{}' after {:.1}s",
                                    label,
                                    query_started.elapsed().as_secs_f32()
                                ),
                            });
                            cancelled = true;
                            break;
                        }
                        let _ = tx.send(GenomeBlastTaskMessage::Status {
                            job_id,
                            status: format!(
                                "BLAST query '{}' failed after {:.1}s: {}",
                                label,
                                query_started.elapsed().as_secs_f32(),
                                e.message
                            ),
                        });
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

            let completed_queries = reports.len() + failed_queries.len();
            let done_result = if cancelled {
                Err(format!(
                    "BLAST cancelled by user after {} / {} quer{}",
                    completed_queries,
                    total_queries,
                    if total_queries == 1 { "y" } else { "ies" }
                ))
            } else if reports.is_empty() && !failed_queries.is_empty() {
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
                        let canceling = task.cancel_requested.load(Ordering::Relaxed);
                        self.genome_blast_progress_fraction = Some(fraction);
                        self.genome_blast_progress_label = format!(
                            "{done_queries} / {total_queries} ({current_query_label}){}",
                            if canceling { " [cancel requested]" } else { "" }
                        );
                    }
                    Ok(GenomeBlastTaskMessage::Status { job_id, status }) => {
                        if job_id != active_job_id {
                            if !stale_job_ids.contains(&job_id) {
                                stale_job_ids.push(job_id);
                            }
                            continue;
                        }
                        self.genome_blast_status = status;
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
            let cancellation_requested = self
                .genome_blast_task
                .as_ref()
                .map(|task| task.cancel_requested.load(Ordering::Relaxed))
                .unwrap_or(false);
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
                    if cancellation_requested {
                        self.genome_blast_status = format!(
                            "BLAST finished after cancellation request in {:.1}s: {} query result(s), {} total hit(s)",
                            elapsed,
                            self.genome_blast_results.len(),
                            hit_total
                        );
                        self.push_job_event(
                            BackgroundJobKind::BlastGenome,
                            BackgroundJobEventPhase::Completed,
                            Some(job_id),
                            format!(
                                "BLAST finished after cancellation request in {:.1}s: {} result(s), {} hit(s)",
                                elapsed,
                                self.genome_blast_results.len(),
                                hit_total
                            ),
                        );
                    } else if batch.failed_queries.is_empty() {
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
                    if cancellation_requested
                        || crate::genomes::is_blast_cancelled_error(&e)
                        || e.to_ascii_lowercase().contains("cancel")
                    {
                        self.genome_blast_status =
                            format!("BLAST cancelled after {:.1}s: {}", elapsed, e);
                    } else {
                        self.genome_blast_status =
                            format!("BLAST failed after {:.1}s: {}", elapsed, e);
                    }
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
            blast_provenance: Some(BlastInvocationProvenance {
                genome_id: result.report.genome_id.clone(),
                query_label: result.query_label.clone(),
                query_length: result.query_length,
                max_hits: result.report.max_hits,
                task: result.report.task.clone(),
                blastn_executable: result.report.blastn_executable.clone(),
                blast_db_prefix: result.report.blast_db_prefix.clone(),
                command: result.report.command.clone(),
                command_line: Self::format_blast_command_line(
                    &result.report.blastn_executable,
                    &result.report.command,
                ),
                catalog_path: self.genome_catalog_path_opt(),
                cache_dir: self.genome_cache_dir_opt(),
                options_override_json: result.report.options_override_json.clone(),
                effective_options_json: result.report.effective_options_json.clone(),
            }),
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
        self.genome_retrieve_contig_suggestions.clear();
        if self.genome_output_id.trim().is_empty() || self.genome_output_id_autofilled {
            self.genome_output_id = Self::default_retrieve_genome_output_id(
                &self.genome_id,
                gene,
                self.genome_gene_extract_mode,
                self.current_gene_promoter_upstream_bp(),
            );
            self.genome_output_id_autofilled = true;
        }
    }

    fn normalize_output_id_token(raw: &str) -> String {
        let mut out = String::new();
        let mut previous_underscore = false;
        for ch in raw.chars() {
            if ch.is_ascii_alphanumeric() {
                out.push(ch.to_ascii_lowercase());
                previous_underscore = false;
            } else if !previous_underscore {
                out.push('_');
                previous_underscore = true;
            }
        }
        let trimmed = out.trim_matches('_');
        if trimmed.is_empty() {
            "region".to_string()
        } else {
            trimmed.to_string()
        }
    }

    fn default_retrieve_genome_output_id(
        genome_id: &str,
        gene: &GenomeGeneRecord,
        extract_mode: GenomeGeneExtractMode,
        promoter_upstream_bp: usize,
    ) -> String {
        let label = gene
            .gene_name
            .as_deref()
            .or(gene.gene_id.as_deref())
            .unwrap_or("gene");
        let genome_token = Self::normalize_output_id_token(genome_id);
        let label_token = Self::normalize_output_id_token(label);
        match extract_mode {
            GenomeGeneExtractMode::Gene => format!(
                "{}_{}_{}_{}",
                genome_token, label_token, gene.start_1based, gene.end_1based
            ),
            GenomeGeneExtractMode::CodingWithPromoter => {
                if promoter_upstream_bp == 0 {
                    format!("{genome_token}_{label_token}_coding")
                } else {
                    format!(
                        "{genome_token}_{label_token}_coding_promoter_{}bp",
                        promoter_upstream_bp
                    )
                }
            }
        }
    }

    fn current_gene_promoter_upstream_bp(&self) -> usize {
        self.genome_gene_promoter_upstream_bp
            .trim()
            .parse::<usize>()
            .unwrap_or(0)
    }

    fn refresh_selected_gene_output_id_if_autofilled(&mut self) {
        if !self.genome_output_id.trim().is_empty() && !self.genome_output_id_autofilled {
            return;
        }
        let Some(selected_idx) = self.genome_selected_gene else {
            return;
        };
        let Some(gene) = self.genome_genes.get(selected_idx) else {
            return;
        };
        self.genome_output_id = Self::default_retrieve_genome_output_id(
            &self.genome_id,
            gene,
            self.genome_gene_extract_mode,
            self.current_gene_promoter_upstream_bp(),
        );
        self.genome_output_id_autofilled = true;
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
        if let (Some(gene_name), Some(gene_id)) = (&gene.gene_name, &gene.gene_id)
            && gene_name != gene_id
        {
            name = format!("{gene_name} ({gene_id})");
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

    fn sync_genome_annotation_scope_controls(&mut self) {
        if !self.genome_include_genomic_annotation {
            self.genome_annotation_scope = GenomeAnnotationScope::None;
        } else if matches!(self.genome_annotation_scope, GenomeAnnotationScope::None) {
            self.genome_annotation_scope = GenomeAnnotationScope::Core;
        }
    }

    fn extract_reference_genome_gene(&mut self) {
        let genome_id = self.genome_id.trim().to_string();
        if genome_id.is_empty() {
            self.genome_retrieve_status = "Select a genome first".to_string();
            self.genome_retrieve_contig_suggestions.clear();
            return;
        }
        let Some((gene_query, occurrence)) = self.selected_gene_query_and_occurrence() else {
            self.genome_retrieve_status =
                "Select a gene with a gene_id or gene_name first".to_string();
            self.genome_retrieve_contig_suggestions.clear();
            return;
        };
        let output_id = if self.genome_output_id.trim().is_empty() {
            None
        } else {
            Some(self.genome_output_id.trim().to_string())
        };
        let extract_mode = self.genome_gene_extract_mode;
        let promoter_upstream_bp =
            if matches!(extract_mode, GenomeGeneExtractMode::CodingWithPromoter) {
                if self.genome_gene_promoter_upstream_bp.trim().is_empty() {
                    0
                } else {
                    match self
                        .genome_gene_promoter_upstream_bp
                        .trim()
                        .parse::<usize>()
                    {
                        Ok(parsed) => parsed,
                        Err(_) => {
                            self.genome_retrieve_status =
                                "promoter bp before CDS must be a non-negative integer".to_string();
                            self.genome_retrieve_contig_suggestions.clear();
                            return;
                        }
                    }
                }
            } else {
                0
            };
        self.sync_genome_annotation_scope_controls();
        let annotation_scope = self.genome_annotation_scope;
        let max_annotation_features = if self.genome_max_annotation_features.trim().is_empty() {
            None
        } else {
            match self.genome_max_annotation_features.trim().parse::<usize>() {
                Ok(parsed) => Some(parsed),
                Err(_) => {
                    self.genome_retrieve_status = "max annotation features must be a non-negative integer (empty = unlimited)".to_string();
                    self.genome_retrieve_contig_suggestions.clear();
                    return;
                }
            }
        };
        let op = Operation::ExtractGenomeGene {
            genome_id,
            gene_query,
            occurrence: Some(occurrence),
            output_id,
            extract_mode: Some(extract_mode),
            promoter_upstream_bp: matches!(extract_mode, GenomeGeneExtractMode::CodingWithPromoter)
                .then_some(promoter_upstream_bp),
            annotation_scope: Some(annotation_scope),
            max_annotation_features,
            include_genomic_annotation: Some(self.genome_include_genomic_annotation),
            catalog_path: self.genome_catalog_path_opt(),
            cache_dir: self.genome_cache_dir_opt(),
        };
        let chromosome_for_suggestions = self.genome_chromosome.trim().to_string();
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
                self.genome_retrieve_contig_suggestions.clear();
            }
            Err(e) => {
                self.genome_retrieve_status = format!("Extract gene failed: {}", e.message);
                self.update_retrieve_contig_suggestions_from_error(
                    &chromosome_for_suggestions,
                    &e.message,
                );
            }
        }
    }

    fn extract_reference_genome_region(&mut self) {
        let genome_id = self.genome_id.trim().to_string();
        if genome_id.is_empty() {
            self.genome_retrieve_status = "Select a genome first".to_string();
            self.genome_retrieve_contig_suggestions.clear();
            return;
        }
        let chromosome = self.genome_chromosome.trim().to_string();
        if chromosome.is_empty() {
            self.genome_retrieve_status = "Chromosome cannot be empty".to_string();
            self.genome_retrieve_contig_suggestions.clear();
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
                self.genome_retrieve_contig_suggestions.clear();
                return;
            }
        };
        let output_id = if self.genome_output_id.trim().is_empty() {
            None
        } else {
            Some(self.genome_output_id.trim().to_string())
        };
        self.sync_genome_annotation_scope_controls();
        let annotation_scope = self.genome_annotation_scope;
        let max_annotation_features = if self.genome_max_annotation_features.trim().is_empty() {
            None
        } else {
            match self.genome_max_annotation_features.trim().parse::<usize>() {
                Ok(parsed) => Some(parsed),
                Err(_) => {
                    self.genome_retrieve_status = "max annotation features must be a non-negative integer (empty = unlimited)".to_string();
                    self.genome_retrieve_contig_suggestions.clear();
                    return;
                }
            }
        };
        let catalog_path = self.genome_catalog_path_opt();
        let extract_op = Operation::ExtractGenomeRegion {
            genome_id: genome_id.clone(),
            chromosome: chromosome.clone(),
            start_1based,
            end_1based,
            output_id,
            annotation_scope: Some(annotation_scope),
            max_annotation_features,
            include_genomic_annotation: Some(self.genome_include_genomic_annotation),
            catalog_path,
            cache_dir: self.genome_cache_dir_opt(),
        };
        let result = { self.engine.write().unwrap().apply(extract_op) };
        match result {
            Ok(r) => {
                for seq_id in &r.created_seq_ids {
                    self.open_sequence_window(seq_id);
                }
                self.genome_retrieve_status = Self::format_extract_region_status(&r);
                self.genome_retrieve_contig_suggestions.clear();
            }
            Err(e) => {
                self.genome_retrieve_status = format!("Extract region failed: {}", e.message);
                self.update_retrieve_contig_suggestions_from_error(&chromosome, &e.message);
            }
        }
    }

    fn reopen_pcr_designer_from_operation(
        &mut self,
        op_id: &str,
    ) -> std::result::Result<bool, String> {
        let Some(seq_id) = self.lineage_reopenable_pcr_op_seq_ids.get(op_id).cloned() else {
            return Ok(false);
        };
        self.open_pcr_design_dialog_for_seq_id(&seq_id)?;
        self.app_status = format!(
            "Opened PCR Designer on template '{}' from PCR-related operation '{}'",
            seq_id, op_id
        );
        Ok(true)
    }

    fn render_arrangement_gel_preview_contents(&mut self, ui: &mut Ui) -> bool {
        let mut close_requested = false;
        let close_hover = Self::specialist_window_close_hover_text("Arrangement Gel");
        if self.render_specialist_window_nav_with_close(ui, Some(("Close", close_hover.as_str()))) {
            close_requested = true;
        }
        ui.label("Preview one serial arrangement inside GENtle, tune left/right DNA ladders, then save or export the chosen setup.");
        ui.horizontal(|ui| {
            ui.label("Arrangement");
            ui.monospace(self.arrangement_gel_preview.arrangement_title.clone());
        });
        ui.horizontal(|ui| {
            ui.label("Saved ladders");
            ui.monospace(Self::describe_arrangement_ladders(
                &self.arrangement_gel_preview.saved_ladders,
            ));
            let preview_ladders = self.arrangement_gel_preview_effective_ladders();
            if preview_ladders != self.arrangement_gel_preview.saved_ladders {
                ui.small("Preview differs from the saved arrangement until you click 'Save to Arrangement'.");
            }
        });
        let ladder_names = GentleEngine::inspect_dna_ladders(None)
            .ladders
            .into_iter()
            .map(|ladder| ladder.name)
            .collect::<Vec<_>>();
        let mut selection_changed = false;
        ui.horizontal(|ui| {
            ui.label("Left ladder");
            egui::ComboBox::from_id_salt("arrangement_gel_left_ladder")
                .selected_text(
                    if self
                        .arrangement_gel_preview
                        .left_ladder_name
                        .trim()
                        .is_empty()
                    {
                        "Auto".to_string()
                    } else {
                        self.arrangement_gel_preview.left_ladder_name.clone()
                    },
                )
                .show_ui(ui, |ui| {
                    selection_changed |= ui
                        .selectable_value(
                            &mut self.arrangement_gel_preview.left_ladder_name,
                            String::new(),
                            "Auto",
                        )
                        .changed();
                    for ladder_name in &ladder_names {
                        selection_changed |= ui
                            .selectable_value(
                                &mut self.arrangement_gel_preview.left_ladder_name,
                                ladder_name.clone(),
                                ladder_name,
                            )
                            .changed();
                    }
                });
            ui.label("Right ladder");
            egui::ComboBox::from_id_salt("arrangement_gel_right_ladder")
                .selected_text(
                    if self
                        .arrangement_gel_preview
                        .right_ladder_name
                        .trim()
                        .is_empty()
                    {
                        "Auto".to_string()
                    } else {
                        self.arrangement_gel_preview.right_ladder_name.clone()
                    },
                )
                .show_ui(ui, |ui| {
                    selection_changed |= ui
                        .selectable_value(
                            &mut self.arrangement_gel_preview.right_ladder_name,
                            String::new(),
                            "Auto",
                        )
                        .changed();
                    for ladder_name in &ladder_names {
                        selection_changed |= ui
                            .selectable_value(
                                &mut self.arrangement_gel_preview.right_ladder_name,
                                ladder_name.clone(),
                                ladder_name,
                            )
                            .changed();
                    }
                });
            if ui
                .button("Auto")
                .on_hover_text(
                    "Reset both flanking ladders to automatic selection for this preview",
                )
                .clicked()
            {
                self.arrangement_gel_preview.left_ladder_name.clear();
                self.arrangement_gel_preview.right_ladder_name.clear();
                selection_changed = true;
            }
        });
        if selection_changed {
            self.coerce_arrangement_gel_preview_pair();
            self.refresh_arrangement_gel_preview_svg();
        }

        let mut condition_changed = false;
        ui.horizontal_wrapped(|ui| {
            ui.label("Conditions");
            ui.label("Agarose %");
            condition_changed |= ui
                .add(
                    egui::TextEdit::singleline(
                        &mut self.arrangement_gel_preview.agarose_percent_text,
                    )
                    .desired_width(52.0),
                )
                .changed();
            ui.separator();
            ui.label("Buffer");
            egui::ComboBox::from_id_salt("arrangement_gel_buffer_model")
                .selected_text(
                    self.arrangement_gel_preview
                        .buffer_model
                        .as_str()
                        .to_ascii_uppercase(),
                )
                .show_ui(ui, |ui| {
                    condition_changed |= ui
                        .selectable_value(
                            &mut self.arrangement_gel_preview.buffer_model,
                            crate::engine::GelBufferModel::Tae,
                            "TAE",
                        )
                        .changed();
                    condition_changed |= ui
                        .selectable_value(
                            &mut self.arrangement_gel_preview.buffer_model,
                            crate::engine::GelBufferModel::Tbe,
                            "TBE",
                        )
                        .changed();
                });
            ui.separator();
            condition_changed |= ui
                .checkbox(
                    &mut self.arrangement_gel_preview.topology_aware,
                    "Topology-aware circular migration",
                )
                .changed();
            if ui
                .button("Reset conditions")
                .on_hover_text("Restore the deterministic default gel conditions")
                .clicked()
            {
                let defaults = crate::engine::GelRunConditions::default();
                self.arrangement_gel_preview.agarose_percent_text =
                    format!("{:.1}", defaults.agarose_percent);
                self.arrangement_gel_preview.buffer_model = defaults.buffer_model;
                self.arrangement_gel_preview.topology_aware = defaults.topology_aware;
                condition_changed = true;
            }
        });
        if condition_changed {
            self.refresh_arrangement_gel_preview_svg();
        }

        ui.horizontal(|ui| {
            if ui
                .button("Preview")
                .on_hover_text("Recompute the in-window gel preview using the current ladder choices")
                .clicked()
            {
                self.refresh_arrangement_gel_preview_svg();
            }
            if ui
                .button("Save to Arrangement")
                .on_hover_text("Persist the current ladder choices on this arrangement for future export/reopen")
                .clicked()
            {
                self.save_arrangement_gel_preview_ladders();
            }
            if ui
                .button("Export SVG...")
                .on_hover_text("Export one serial gel SVG using the current preview ladder choices")
                .clicked()
            {
                let stem = Self::sanitize_file_stem(
                    &self.arrangement_gel_preview.arrangement_title,
                    "arrangement_gel",
                );
                self.prompt_export_serial_gel_svg(
                    &stem,
                    None,
                    Some(self.arrangement_gel_preview.arrangement_id.clone()),
                    Some(self.arrangement_gel_preview_effective_ladders()),
                    self.arrangement_gel_preview_effective_conditions().ok(),
                );
            }
        });

        if !self.arrangement_gel_preview.status.trim().is_empty() {
            ui.small(&self.arrangement_gel_preview.status);
        }
        ui.separator();
        ui.label("Gel preview");
        if !self.arrangement_gel_preview.svg_uri.trim().is_empty() {
            ui.add(
                egui::Image::from_uri(self.arrangement_gel_preview.svg_uri.clone())
                    .max_width(ui.available_width())
                    .max_height(520.0)
                    .shrink_to_fit(),
            );
        } else {
            ui.small("No arrangement gel preview is loaded yet.");
        }
        close_requested
    }

    fn render_arrangement_gel_preview_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_arrangement_gel_preview_dialog {
            return;
        }
        let mut open = self.show_arrangement_gel_preview_dialog;
        let title = self.arrangement_gel_preview_title();
        let viewport_id = Self::arrangement_gel_preview_viewport_id();
        let spec = self.hosted_window_spec_for_viewport(
            title.clone(),
            egui::Id::new(("hosted_arrangement_gel_preview_window", viewport_id)),
            viewport_id,
            Vec2::new(980.0, 760.0),
            Vec2::new(760.0, 560.0),
        );
        if ctx.embed_viewports() {
            let mut close_requested = false;
            crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
                egui::ScrollArea::vertical()
                    .id_salt("arrangement_gel_embedded_scroll")
                    .auto_shrink([false, false])
                    .show(ui, |ui| {
                        close_requested = self.render_arrangement_gel_preview_contents(ui);
                    });
            });
            self.clear_viewport_foreground_request_after_render(viewport_id);
            if close_requested || ctx.input(|i| i.key_pressed(Key::Escape)) {
                open = false;
            }
            self.show_arrangement_gel_preview_dialog = open;
            return;
        }
        let builder = crate::egui_compat::viewport_builder_for_hosted_window(&spec);
        ctx.show_viewport_immediate(viewport_id, builder, |ctx, class| {
            self.note_viewport_focus_if_active(ctx, viewport_id);
            if class == egui::ViewportClass::EmbeddedWindow {
                let mut close_requested = false;
                crate::egui_compat::show_hosted_window(&mut *ctx, &spec, &mut open, |ui| {
                    egui::ScrollArea::vertical()
                        .id_salt("arrangement_gel_embedded_scroll")
                        .auto_shrink([false, false])
                        .show(ui, |ui| {
                            close_requested = self.render_arrangement_gel_preview_contents(ui);
                        });
                });
                if close_requested {
                    open = false;
                }
            } else {
                let mut close_requested = false;
                crate::egui_compat::show_central_panel(
                    &mut *ctx,
                    egui::CentralPanel::default(),
                    |ui| {
                        egui::ScrollArea::vertical()
                            .id_salt("arrangement_gel_viewport_scroll")
                            .auto_shrink([false, false])
                            .show(ui, |ui| {
                                close_requested = self.render_arrangement_gel_preview_contents(ui);
                            });
                    },
                );
                if close_requested || Self::viewport_close_requested_or_shortcut(ctx) {
                    open = false;
                }
            }
        });
        if ctx.input(|i| i.key_pressed(Key::Escape)) {
            open = false;
        }
        self.show_arrangement_gel_preview_dialog = open;
    }

    fn render_rack_labels_preview_contents(&mut self, ui: &mut Ui) -> bool {
        let mut close_requested = false;
        ui.horizontal(|ui| {
            if ui
                .button("Close")
                .on_hover_text("Close the rack labels preview")
                .clicked()
            {
                close_requested = true;
            }
            if ui
                .button("Export SVG...")
                .on_hover_text("Export the same rack labels SVG that is shown in this preview")
                .clicked()
            {
                let stem = self
                    .rack_labels_preview
                    .arrangement_title
                    .as_deref()
                    .filter(|value| !value.trim().is_empty())
                    .map(|value| Self::sanitize_file_stem(value, "rack_labels"))
                    .unwrap_or_else(|| {
                        Self::sanitize_file_stem(&self.rack_labels_preview.rack_id, "rack_labels")
                    });
                let arrangement_id = self
                    .rack_labels_preview
                    .arrangement_id
                    .as_deref()
                    .map(str::to_string);
                let rack_id = self.rack_labels_preview.rack_id.clone();
                self.prompt_export_rack_labels_svg_scoped(
                    &rack_id,
                    arrangement_id.as_deref(),
                    &stem,
                );
            }
        });
        ui.separator();
        ui.horizontal(|ui| {
            ui.label("Rack");
            ui.monospace(self.rack_labels_preview.rack_title.clone());
        });
        if let Some(arrangement_title) = self
            .rack_labels_preview
            .arrangement_title
            .as_deref()
            .filter(|value| !value.trim().is_empty())
        {
            ui.horizontal(|ui| {
                ui.label("Arrangement");
                ui.monospace(arrangement_title.to_string());
            });
        }
        let preset_before = self.rack_label_sheet_preset;
        ui.horizontal(|ui| {
            ui.label("Label preset");
            egui::ComboBox::from_id_salt("rack_labels_preview_preset_combo")
                .selected_text(Self::rack_label_sheet_preset_label(
                    self.rack_label_sheet_preset,
                ))
                .show_ui(ui, |ui| {
                    for preset in [
                        RackLabelSheetPreset::CompactCards,
                        RackLabelSheetPreset::PrintA4,
                        RackLabelSheetPreset::WideCards,
                    ] {
                        ui.selectable_value(
                            &mut self.rack_label_sheet_preset,
                            preset,
                            Self::rack_label_sheet_preset_label(preset),
                        );
                    }
                });
        });
        if self.rack_label_sheet_preset != preset_before {
            self.refresh_rack_labels_preview_svg();
        }
        ui.small(
            "Preview the arrangement-scoped or rack-wide label sheet here first, then export the same SVG when it looks ready to print.",
        );
        if !self.rack_labels_preview.status.trim().is_empty() {
            ui.small(self.rack_labels_preview.status.clone());
        }
        ui.separator();
        ui.label("Label preview");
        if !self.rack_labels_preview.svg_uri.trim().is_empty() {
            ui.add(
                egui::Image::from_uri(self.rack_labels_preview.svg_uri.clone())
                    .max_width(ui.available_width())
                    .max_height(620.0)
                    .shrink_to_fit(),
            );
        } else {
            ui.small("No rack label preview is loaded yet.");
        }
        close_requested
    }

    fn render_rack_labels_preview_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_rack_labels_preview_dialog {
            return;
        }
        let mut open = self.show_rack_labels_preview_dialog;
        let title = self.rack_labels_preview_title();
        let viewport_id = Self::rack_labels_preview_viewport_id();
        let spec = self.hosted_window_spec_for_viewport(
            title.clone(),
            egui::Id::new(("hosted_rack_labels_preview_window", viewport_id)),
            viewport_id,
            Vec2::new(980.0, 760.0),
            Vec2::new(760.0, 560.0),
        );
        if ctx.embed_viewports() {
            let mut close_requested = false;
            crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
                egui::ScrollArea::vertical()
                    .auto_shrink([false, false])
                    .show(ui, |ui| {
                        close_requested = self.render_rack_labels_preview_contents(ui);
                    });
            });
            self.clear_viewport_foreground_request_after_render(viewport_id);
            if close_requested || ctx.input(|i| i.key_pressed(Key::Escape)) {
                open = false;
            }
            self.show_rack_labels_preview_dialog = open;
            return;
        }
        let builder = crate::egui_compat::viewport_builder_for_hosted_window(&spec);
        ctx.show_viewport_immediate(viewport_id, builder, |ctx, class| {
            self.note_viewport_focus_if_active(ctx, viewport_id);
            if class == egui::ViewportClass::EmbeddedWindow {
                let mut close_requested = false;
                crate::egui_compat::show_hosted_window(&mut *ctx, &spec, &mut open, |ui| {
                    egui::ScrollArea::vertical()
                        .auto_shrink([false, false])
                        .show(ui, |ui| {
                            close_requested = self.render_rack_labels_preview_contents(ui);
                        });
                });
                if close_requested {
                    open = false;
                }
            } else {
                let mut close_requested = false;
                crate::egui_compat::show_central_panel(
                    &mut *ctx,
                    egui::CentralPanel::default(),
                    |ui| {
                        egui::ScrollArea::vertical()
                            .auto_shrink([false, false])
                            .show(ui, |ui| {
                                close_requested = self.render_rack_labels_preview_contents(ui);
                            });
                    },
                );
                if close_requested || Self::viewport_close_requested_or_shortcut(ctx) {
                    open = false;
                }
            }
        });
        if ctx.input(|i| i.key_pressed(Key::Escape)) {
            open = false;
        }
        self.show_rack_labels_preview_dialog = open;
    }

    fn render_pcr_design_contents(&mut self, ui: &mut Ui, ctx: &egui::Context) -> bool {
        let mut close_requested = false;
        let close_hover = Self::specialist_window_close_hover_text("PCR Designer");
        if self.render_specialist_window_nav_with_close(ui, Some(("Close", close_hover.as_str()))) {
            close_requested = true;
        }
        ui.label("Selection-first pair-PCR specialist. Paint ROI/windows on the map and run deterministic primer-pair design queue operations.");
        if self.pcr_design_seq_id.trim().is_empty() {
            let target = self
                .active_dna_window_context()
                .map(|(seq_id, _)| seq_id)
                .or_else(|| self.project_sequence_ids_for_blast().first().cloned());
            if let Some(seq_id) = target {
                self.pcr_design_seq_id = seq_id;
            }
        }
        if self.pcr_design_seq_id.trim().is_empty() {
            ui.colored_label(
                egui::Color32::from_rgb(190, 70, 70),
                "No sequence is available. Load/open one sequence window first.",
            );
            return close_requested;
        }
        let seq_id = self.pcr_design_seq_id.trim().to_string();
        if self.find_open_sequence_viewport_id(&seq_id).is_none() {
            if ui
                .button("Open target sequence window")
                .on_hover_text("Open the sequence window used as PCR Designer context")
                .clicked()
            {
                self.open_sequence_window(&seq_id);
            }
            ui.small(format!(
                "Target sequence window '{}' is not open yet.",
                seq_id
            ));
            return close_requested;
        }
        let mut rendered = false;
        for window in self.windows.values() {
            let Ok(mut guard) = window.write() else {
                continue;
            };
            if guard.sequence_id().as_deref() == Some(seq_id.as_str()) {
                guard.render_pcr_designer_specialist(ui, ctx);
                rendered = true;
                break;
            }
        }
        if !rendered {
            ui.colored_label(
                egui::Color32::from_rgb(190, 70, 70),
                format!(
                    "Could not lock sequence window '{}' for PCR designer rendering.",
                    seq_id
                ),
            );
        }
        close_requested
    }

    fn render_pcr_design_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_pcr_design_dialog {
            return;
        }
        let mut open = self.show_pcr_design_dialog;
        let title = if self.pcr_design_seq_id.trim().is_empty() {
            "PCR Designer".to_string()
        } else {
            format!("PCR Designer — {}", self.pcr_design_seq_id.trim())
        };
        let viewport_id = Self::pcr_design_viewport_id();
        let spec = self.hosted_window_spec_for_viewport(
            title.clone(),
            egui::Id::new(("hosted_pcr_design_window", viewport_id)),
            viewport_id,
            Vec2::new(1200.0, 840.0),
            Vec2::new(920.0, 620.0),
        );
        if ctx.embed_viewports() {
            let mut close_requested = false;
            crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
                egui::ScrollArea::vertical()
                    .id_salt("pcr_design_embedded_scroll")
                    .auto_shrink([false, false])
                    .show(ui, |ui| {
                        close_requested = self.render_pcr_design_contents(ui, ctx)
                    });
            });
            self.clear_viewport_foreground_request_after_render(viewport_id);
            if close_requested || ctx.input(|i| i.key_pressed(Key::Escape)) {
                open = false;
            }
            self.show_pcr_design_dialog = open;
            return;
        }
        let builder = crate::egui_compat::viewport_builder_for_hosted_window(&spec);
        ctx.show_viewport_immediate(viewport_id, builder, |ctx, class| {
            self.note_viewport_focus_if_active(ctx, viewport_id);
            let viewport_ctx = ctx.ctx().clone();
            if class == egui::ViewportClass::EmbeddedWindow {
                let mut close_requested = false;
                crate::egui_compat::show_hosted_window(&mut *ctx, &spec, &mut open, |ui| {
                    egui::ScrollArea::vertical()
                        .id_salt("pcr_design_embedded_scroll")
                        .auto_shrink([false, false])
                        .show(ui, |ui| {
                            close_requested =
                                self.render_pcr_design_contents(ui, &viewport_ctx)
                        });
                });
                if close_requested {
                    open = false;
                }
            } else {
                let mut close_requested = false;
                crate::egui_compat::show_central_panel(
                    &mut *ctx,
                    egui::CentralPanel::default(),
                    |ui| {
                        egui::ScrollArea::vertical()
                            .id_salt("pcr_design_viewport_scroll")
                            .auto_shrink([false, false])
                            .show(ui, |ui| {
                                close_requested =
                                    self.render_pcr_design_contents(ui, &viewport_ctx)
                            });
                    },
                );
                if close_requested || Self::viewport_close_requested_or_shortcut(ctx) {
                    open = false;
                }
            }
        });
        if ctx.input(|i| i.key_pressed(Key::Escape)) {
            open = false;
        }
        self.show_pcr_design_dialog = open;
    }

    fn render_sequencing_confirmation_contents(
        &mut self,
        ui: &mut Ui,
        ctx: &egui::Context,
    ) -> bool {
        let mut close_requested = false;
        let close_hover = Self::specialist_window_close_hover_text("Sequencing Confirmation");
        if self.render_specialist_window_nav_with_close(ui, Some(("Close", close_hover.as_str()))) {
            close_requested = true;
        }
        ui.label("Construct-confirmation specialist. Review persisted sequencing-confirmation reports or run a new confirmation pass against the current expected construct sequence using called reads and/or imported sequencing traces.");
        if self.sequencing_confirmation_seq_id.trim().is_empty() {
            let target = self
                .active_dna_window_context()
                .map(|(seq_id, _)| seq_id)
                .or_else(|| self.project_sequence_ids_for_blast().first().cloned());
            if let Some(seq_id) = target {
                self.sequencing_confirmation_seq_id = seq_id;
            }
        }
        if self.sequencing_confirmation_seq_id.trim().is_empty() {
            ui.colored_label(
                egui::Color32::from_rgb(190, 70, 70),
                "No sequence is available. Load/open one sequence window first.",
            );
            return close_requested;
        }
        let seq_id = self.sequencing_confirmation_seq_id.trim().to_string();
        if self.find_open_sequence_viewport_id(&seq_id).is_none() {
            if ui
                .button("Open expected construct sequence window")
                .on_hover_text("Open the sequence window used as sequencing-confirmation context")
                .clicked()
            {
                self.open_sequence_window(&seq_id);
            }
            ui.small(format!(
                "Expected construct window '{}' is not open yet.",
                seq_id
            ));
            return close_requested;
        }
        let mut rendered = false;
        for window in self.windows.values() {
            let Ok(mut guard) = window.write() else {
                continue;
            };
            if guard.sequence_id().as_deref() == Some(seq_id.as_str()) {
                guard.render_sequencing_confirmation_specialist(ui, ctx);
                rendered = true;
                break;
            }
        }
        if !rendered {
            ui.colored_label(
                egui::Color32::from_rgb(190, 70, 70),
                format!(
                    "Could not lock sequence window '{}' for sequencing-confirmation rendering.",
                    seq_id
                ),
            );
        }
        close_requested
    }

    fn render_sequencing_confirmation_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_sequencing_confirmation_dialog {
            return;
        }
        let mut open = self.show_sequencing_confirmation_dialog;
        let title = if self.sequencing_confirmation_seq_id.trim().is_empty() {
            "Sequencing Confirmation".to_string()
        } else {
            format!(
                "Sequencing Confirmation — {}",
                self.sequencing_confirmation_seq_id.trim()
            )
        };
        let viewport_id = Self::sequencing_confirmation_viewport_id();
        let spec = self.hosted_window_spec_for_viewport(
            title.clone(),
            egui::Id::new(("hosted_sequencing_confirmation_window", viewport_id)),
            viewport_id,
            Vec2::new(1180.0, 820.0),
            Vec2::new(920.0, 620.0),
        );
        if ctx.embed_viewports() {
            let mut close_requested = false;
            crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
                egui::ScrollArea::vertical()
                    .id_salt("sequencing_confirmation_embedded_scroll")
                    .auto_shrink([false, false])
                    .show(ui, |ui| {
                        close_requested = self.render_sequencing_confirmation_contents(ui, ctx)
                    });
            });
            self.clear_viewport_foreground_request_after_render(viewport_id);
            if close_requested || ctx.input(|i| i.key_pressed(Key::Escape)) {
                open = false;
            }
            self.show_sequencing_confirmation_dialog = open;
            return;
        }
        let builder = crate::egui_compat::viewport_builder_for_hosted_window(&spec);
        ctx.show_viewport_immediate(viewport_id, builder, |ctx, class| {
            self.note_viewport_focus_if_active(ctx, viewport_id);
            let viewport_ctx = ctx.ctx().clone();
            if class == egui::ViewportClass::EmbeddedWindow {
                let mut close_requested = false;
                crate::egui_compat::show_hosted_window(&mut *ctx, &spec, &mut open, |ui| {
                    egui::ScrollArea::vertical()
                        .id_salt("sequencing_confirmation_embedded_scroll")
                        .auto_shrink([false, false])
                        .show(ui, |ui| {
                            close_requested = self
                                .render_sequencing_confirmation_contents(ui, &viewport_ctx)
                        });
                });
                if close_requested {
                    open = false;
                }
            } else {
                let mut close_requested = false;
                crate::egui_compat::show_central_panel(
                    &mut *ctx,
                    egui::CentralPanel::default(),
                    |ui| {
                        egui::ScrollArea::vertical()
                            .id_salt("sequencing_confirmation_viewport_scroll")
                            .auto_shrink([false, false])
                            .show(ui, |ui| {
                                close_requested = self
                                    .render_sequencing_confirmation_contents(ui, &viewport_ctx)
                            });
                    },
                );
                if close_requested || Self::viewport_close_requested_or_shortcut(ctx) {
                    open = false;
                }
            }
        });
        if ctx.input(|i| i.key_pressed(Key::Escape)) {
            open = false;
        }
        self.show_sequencing_confirmation_dialog = open;
    }

    fn render_planning_contents(&mut self, ui: &mut Ui) -> bool {
        let mut close_requested = false;
        let close_hover = Self::specialist_window_close_hover_text("Planning");
        if self.render_specialist_window_nav_with_close(ui, Some(("Close", close_hover.as_str()))) {
            close_requested = true;
        }
        ui.label(
            "Planning meta-layer controls for time/cost/local-fit routine ranking and advisory sync suggestions.",
        );
        ui.small(
            "Effective merge precedence: global -> confirmed_agent_overlay -> project_override",
        );
        ui.horizontal(|ui| {
            if ui
                .button("Reload From Engine")
                .on_hover_text("Refresh all planning editor buffers from current project state")
                .clicked()
            {
                self.refresh_planning_editor_buffers_from_engine();
                self.planning_status = "Reloaded planning buffers from engine state".to_string();
            }
        });
        ui.separator();

        ui.heading("Profiles");
        ui.small("Empty editor text means clear on Apply.");

        ui.label("Global profile");
        ui.add(
            egui::TextEdit::multiline(&mut self.planning_profile_global_json)
                .desired_rows(7)
                .desired_width(f32::INFINITY),
        );
        ui.horizontal(|ui| {
            if ui
                .button("Apply Global")
                .on_hover_text("Set planning profile scope 'global' from this JSON")
                .clicked()
            {
                let payload = self.planning_profile_global_json.clone();
                self.apply_planning_profile_json(PlanningProfileScope::Global, &payload);
            }
            if ui
                .button("Clear Global")
                .on_hover_text("Clear planning profile scope 'global'")
                .clicked()
            {
                self.clear_planning_profile_scope(PlanningProfileScope::Global);
            }
        });

        ui.separator();
        ui.label("Confirmed agent overlay profile");
        ui.add(
            egui::TextEdit::multiline(&mut self.planning_profile_overlay_json)
                .desired_rows(7)
                .desired_width(f32::INFINITY),
        );
        ui.horizontal(|ui| {
            if ui
                .button("Apply Agent Overlay")
                .on_hover_text(
                    "Set planning profile scope 'confirmed_agent_overlay' from this JSON",
                )
                .clicked()
            {
                let payload = self.planning_profile_overlay_json.clone();
                self.apply_planning_profile_json(
                    PlanningProfileScope::ConfirmedAgentOverlay,
                    &payload,
                );
            }
            if ui
                .button("Clear Agent Overlay")
                .on_hover_text("Clear planning profile scope 'confirmed_agent_overlay'")
                .clicked()
            {
                self.clear_planning_profile_scope(PlanningProfileScope::ConfirmedAgentOverlay);
            }
        });

        ui.separator();
        ui.label("Project override profile");
        ui.add(
            egui::TextEdit::multiline(&mut self.planning_profile_project_json)
                .desired_rows(7)
                .desired_width(f32::INFINITY),
        );
        ui.horizontal(|ui| {
            if ui
                .button("Apply Project Override")
                .on_hover_text("Set planning profile scope 'project_override' from this JSON")
                .clicked()
            {
                let payload = self.planning_profile_project_json.clone();
                self.apply_planning_profile_json(PlanningProfileScope::ProjectOverride, &payload);
            }
            if ui
                .button("Clear Project Override")
                .on_hover_text("Clear planning profile scope 'project_override'")
                .clicked()
            {
                self.clear_planning_profile_scope(PlanningProfileScope::ProjectOverride);
            }
        });

        ui.separator();
        ui.label("Effective profile (read-only)");
        let mut effective_profile_json = {
            let engine = self.engine.read().unwrap();
            Self::pretty_json_or_fallback(&engine.planning_effective_profile())
        };
        ui.add_enabled(
            false,
            egui::TextEdit::multiline(&mut effective_profile_json)
                .desired_rows(7)
                .desired_width(f32::INFINITY),
        );

        ui.separator();
        ui.heading("Objective");
        ui.add(
            egui::TextEdit::multiline(&mut self.planning_objective_json)
                .desired_rows(7)
                .desired_width(f32::INFINITY),
        );
        ui.horizontal(|ui| {
            if ui
                .button("Apply Objective")
                .on_hover_text("Set planning objective from this JSON")
                .clicked()
            {
                let payload = self.planning_objective_json.clone();
                self.apply_planning_objective_json(&payload);
            }
            if ui
                .button("Clear Objective")
                .on_hover_text("Clear planning objective and use engine defaults")
                .clicked()
            {
                self.clear_planning_objective();
            }
        });

        ui.separator();
        ui.heading("Host Profile Browser");
        ui.small(format!(
            "Browse starter host/strain records used by construct reasoning: {}",
            DEFAULT_HOST_PROFILE_CATALOG_PATH
        ));
        ui.horizontal(|ui| {
            ui.label("filter");
            ui.add(
                egui::TextEdit::singleline(&mut self.planning_host_filter)
                    .desired_width(f32::INFINITY),
            );
            if ui
                .button("Clear")
                .on_hover_text("Clear the host-profile browser filter")
                .clicked()
            {
                self.planning_host_filter.clear();
            }
        });
        match self.host_profile_catalog_entries(None, Some(&self.planning_host_filter)) {
            Ok(entries) => {
                if entries.is_empty() {
                    ui.small("No host profiles matched the current filter.");
                } else {
                    if !entries.iter().any(|profile| {
                        Self::host_profile_matches_selection(
                            profile,
                            &self.planning_host_selected_id,
                        )
                    }) {
                        self.planning_host_selected_id = entries[0].profile_id.clone();
                    }
                    ui.columns(2, |columns| {
                        columns[0].label("Catalog entries");
                        egui::ScrollArea::vertical()
                            .id_salt("planning_host_profile_entries")
                            .max_height(220.0)
                            .show(&mut columns[0], |ui| {
                                scroll_input_policy::apply_scrollarea_keyboard_navigation(
                                    ui,
                                    scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                                );
                                for profile in &entries {
                                    let selected = Self::host_profile_matches_selection(
                                        profile,
                                        &self.planning_host_selected_id,
                                    );
                                    let label =
                                        format!("{} — {}", profile.profile_id, profile.strain);
                                    if ui.selectable_label(selected, label).clicked() {
                                        self.planning_host_selected_id = profile.profile_id.clone();
                                    }
                                }
                            });
                        columns[1].label("Selected host profile");
                        if let Some(profile) = entries.iter().find(|profile| {
                            Self::host_profile_matches_selection(
                                profile,
                                &self.planning_host_selected_id,
                            )
                        }) {
                            Self::render_host_profile_record_panel(
                                &mut columns[1],
                                profile,
                                "Selected Host Profile",
                            );
                        }
                    });
                }
            }
            Err(e) => {
                ui.colored_label(
                    egui::Color32::from_rgb(190, 70, 70),
                    format!("Host-profile browser error: {e}"),
                );
            }
        }

        ui.separator();
        ui.heading("Helper Construct Browser");
        let (helper_catalog_path, _) = self.scope_genome_paths_resolved(GenomeDialogScope::Helper);
        ui.small(format!(
            "Browse helper semantics from the current helper catalog: {helper_catalog_path}"
        ));
        match self.helper_vector_doctor_issues_for_catalog_path(&helper_catalog_path) {
            Ok(issues) => {
                if issues.is_empty() {
                    ui.small("catalog doctor: no deterministic vector-catalog issues");
                } else {
                    ui.collapsing(format!("Catalog Doctor Issues ({})", issues.len()), |ui| {
                        for issue in &issues {
                            ui.monospace(Self::format_helper_vector_doctor_issue(issue));
                        }
                    });
                }
            }
            Err(e) => {
                ui.colored_label(
                    egui::Color32::from_rgb(190, 70, 70),
                    format!("Helper catalog doctor error: {e}"),
                );
            }
        }
        ui.horizontal(|ui| {
            ui.label("filter");
            ui.add(
                egui::TextEdit::singleline(&mut self.planning_helper_filter)
                    .desired_width(f32::INFINITY),
            );
            if ui
                .button("Clear")
                .on_hover_text("Clear the helper-construct browser filter")
                .clicked()
            {
                self.planning_helper_filter.clear();
            }
        });
        let helper_cards = self.helper_vector_cards_for_catalog_path(
            &helper_catalog_path,
            Some(&self.planning_helper_filter),
        );
        if let Err(e) = helper_cards.as_ref() {
            ui.colored_label(
                egui::Color32::from_rgb(190, 70, 70),
                format!("Helper vector-card error: {e}"),
            );
        }
        match self.helper_catalog_entries_for_catalog_path(
            &helper_catalog_path,
            Some(&self.planning_helper_filter),
        ) {
            Ok(entries) => {
                if entries.is_empty() {
                    ui.small("No helper catalog entries matched the current filter.");
                } else {
                    if !entries.iter().any(|entry| {
                        Self::helper_catalog_entry_matches_selection(
                            entry,
                            &self.planning_helper_selected_id,
                        )
                    }) {
                        self.planning_helper_selected_id = entries[0].genome_id.clone();
                    }
                    ui.columns(2, |columns| {
                        columns[0].label("Catalog entries");
                        egui::ScrollArea::vertical()
                            .id_salt("planning_helper_catalog_entries")
                            .max_height(240.0)
                            .show(&mut columns[0], |ui| {
                                scroll_input_policy::apply_scrollarea_keyboard_navigation(
                                    ui,
                                    scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                                );
                                for entry in &entries {
                                    let selected = Self::helper_catalog_entry_matches_selection(
                                        entry,
                                        &self.planning_helper_selected_id,
                                    );
                                    let label = entry
                                        .summary
                                        .as_deref()
                                        .map(str::trim)
                                        .filter(|value| !value.is_empty())
                                        .map(|summary| format!("{} — {}", entry.genome_id, summary))
                                        .unwrap_or_else(|| entry.genome_id.clone());
                                    if ui.selectable_label(selected, label).clicked() {
                                        self.planning_helper_selected_id = entry.genome_id.clone();
                                    }
                                }
                            });
                        columns[1].label("Selected helper");
                        if let Some(entry) = entries.iter().find(|entry| {
                            Self::helper_catalog_entry_matches_selection(
                                entry,
                                &self.planning_helper_selected_id,
                            )
                        }) {
                            if let Ok(cards) = helper_cards.as_ref()
                                && let Some(card) =
                                    cards.iter().find(|card| card.helper_id == entry.genome_id)
                            {
                                Self::render_helper_vector_card_panel(&mut columns[1], card);
                            }
                            Self::render_helper_catalog_entry_panel(
                                &mut columns[1],
                                entry,
                                "Selected Helper Construct",
                            );
                        }
                    });
                }
            }
            Err(e) => {
                ui.colored_label(
                    egui::Color32::from_rgb(190, 70, 70),
                    format!("Helper browser error: {e}"),
                );
            }
        }

        ui.separator();
        ui.heading("Sync Suggestions");
        ui.small(
            "Register advisory planning patches as pending pull/push suggestions before explicit accept/reject.",
        );
        ui.horizontal(|ui| {
            ui.label("source");
            ui.add(egui::TextEdit::singleline(&mut self.planning_sync_source).desired_width(180.0));
            ui.label("confidence");
            ui.add(
                egui::TextEdit::singleline(&mut self.planning_sync_confidence).desired_width(80.0),
            );
            ui.label("snapshot");
            ui.add(
                egui::TextEdit::singleline(&mut self.planning_sync_snapshot_id)
                    .desired_width(180.0),
            );
        });
        ui.horizontal(|ui| {
            ui.label("message");
            ui.add(
                egui::TextEdit::singleline(&mut self.planning_sync_message)
                    .desired_width(f32::INFINITY),
            );
        });
        ui.label("sync payload JSON ({ profile_patch?, objective_patch? })");
        ui.add(
            egui::TextEdit::multiline(&mut self.planning_sync_payload_json)
                .desired_rows(6)
                .desired_width(f32::INFINITY),
        );
        ui.horizontal(|ui| {
            if ui
                .button("Register Pull Suggestion")
                .on_hover_text("Create pending suggestion with direction=pull")
                .clicked()
            {
                self.register_planning_sync_suggestion("pull");
            }
            if ui
                .button("Register Push Suggestion")
                .on_hover_text("Create pending suggestion with direction=push")
                .clicked()
            {
                self.register_planning_sync_suggestion("push");
            }
        });

        ui.separator();
        ui.horizontal(|ui| {
            ui.label("suggestion filter");
            egui::ComboBox::from_id_salt("planning_suggestions_filter")
                .selected_text(Self::planning_filter_label(
                    self.planning_suggestions_filter,
                ))
                .show_ui(ui, |ui| {
                    ui.selectable_value(&mut self.planning_suggestions_filter, None, "all");
                    ui.selectable_value(
                        &mut self.planning_suggestions_filter,
                        Some(PlanningSuggestionStatus::Pending),
                        "pending",
                    );
                    ui.selectable_value(
                        &mut self.planning_suggestions_filter,
                        Some(PlanningSuggestionStatus::Accepted),
                        "accepted",
                    );
                    ui.selectable_value(
                        &mut self.planning_suggestions_filter,
                        Some(PlanningSuggestionStatus::Rejected),
                        "rejected",
                    );
                });
            ui.label("reject reason");
            ui.add(
                egui::TextEdit::singleline(&mut self.planning_rejection_reason)
                    .desired_width(260.0),
            );
        });
        let (sync_status, suggestions) = {
            let engine = self.engine.read().unwrap();
            (
                engine.planning_sync_status(),
                engine.list_planning_suggestions(self.planning_suggestions_filter),
            )
        };
        ui.small(format!(
            "sync status: pending={} | last_pull={:?} | last_push={:?} | source={} | snapshot={} | last_error={}",
            sync_status.pending_suggestion_count,
            sync_status.last_pull_at_unix_ms,
            sync_status.last_push_at_unix_ms,
            sync_status
                .last_source
                .as_deref()
                .unwrap_or("-"),
            sync_status
                .last_snapshot_id
                .as_deref()
                .unwrap_or("-"),
            sync_status.last_error.as_deref().unwrap_or("-")
        ));

        if suggestions.is_empty() {
            ui.small("No planning suggestions for this filter.");
        } else {
            let mut accept_id: Option<String> = None;
            let mut reject_id: Option<String> = None;
            egui::ScrollArea::vertical()
                .id_salt("planning_suggestions_scroll")
                .max_height(320.0)
                .show(ui, |ui| {
                    scroll_input_policy::apply_scrollarea_keyboard_navigation(
                        ui,
                        scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                    );
                    for suggestion in suggestions.iter().rev() {
                        ui.group(|ui| {
                            ui.horizontal_wrapped(|ui| {
                                ui.monospace(&suggestion.suggestion_id);
                                ui.label(format!("status={}", suggestion.status.as_str()));
                                ui.label(format!("direction={}", suggestion.direction));
                                ui.label(format!("source={}", suggestion.source));
                                if let Some(confidence) = suggestion.confidence {
                                    ui.label(format!("confidence={confidence:.2}"));
                                }
                                ui.label(format!("created={}", suggestion.created_at_unix_ms));
                                if let Some(resolved_at) = suggestion.resolved_at_unix_ms {
                                    ui.label(format!("resolved={resolved_at}"));
                                }
                            });
                            if let Some(message) = suggestion.message.as_deref()
                                && !message.trim().is_empty() {
                                    ui.small(format!("message: {message}"));
                                }
                            if let Some(reason) = suggestion.rejection_reason.as_deref()
                                && !reason.trim().is_empty() {
                                    ui.small(format!("rejection_reason: {reason}"));
                                }
                            let mut diff_json = Self::pretty_json_or_fallback(&suggestion.diff);
                            ui.add_enabled(
                                false,
                                egui::TextEdit::multiline(&mut diff_json)
                                    .desired_rows(4)
                                    .desired_width(f32::INFINITY),
                            );
                            if suggestion.status == PlanningSuggestionStatus::Pending {
                                ui.horizontal(|ui| {
                                    if ui
                                        .button("Accept")
                                        .on_hover_text(
                                            "Accept this suggestion and merge patch into active planning state",
                                        )
                                        .clicked()
                                    {
                                        accept_id = Some(suggestion.suggestion_id.clone());
                                    }
                                    if ui
                                        .button("Reject")
                                        .on_hover_text(
                                            "Reject this suggestion with optional reason text",
                                        )
                                        .clicked()
                                    {
                                        reject_id = Some(suggestion.suggestion_id.clone());
                                    }
                                });
                            }
                        });
                    }
                });
            if let Some(suggestion_id) = accept_id {
                self.accept_planning_suggestion(&suggestion_id);
            }
            if let Some(suggestion_id) = reject_id {
                self.reject_planning_suggestion(&suggestion_id);
            }
        }

        if !self.planning_status.trim().is_empty() {
            ui.separator();
            ui.monospace(self.planning_status.clone());
        }
        close_requested
    }

    fn render_planning_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_planning_dialog {
            return;
        }
        let mut open = self.show_planning_dialog;
        let viewport_id = Self::planning_viewport_id();
        let spec = self.hosted_window_spec_for_viewport(
            "Planning",
            egui::Id::new(("hosted_planning_window", viewport_id)),
            viewport_id,
            Vec2::new(980.0, 760.0),
            Vec2::new(700.0, 520.0),
        );
        if ctx.embed_viewports() {
            let mut close_requested = false;
            crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
                close_requested = self.render_planning_contents(ui)
            });
            self.clear_viewport_foreground_request_after_render(viewport_id);
            if close_requested || ctx.input(|i| i.key_pressed(Key::Escape)) {
                open = false;
            }
            self.show_planning_dialog = open;
            return;
        }
        let builder = crate::egui_compat::viewport_builder_for_hosted_window(&spec);
        ctx.show_viewport_immediate(viewport_id, builder, |ctx, class| {
            self.note_viewport_focus_if_active(ctx, viewport_id);
            if class == egui::ViewportClass::EmbeddedWindow {
                let mut close_requested = false;
                crate::egui_compat::show_hosted_window(&mut *ctx, &spec, &mut open, |ui| {
                    close_requested = self.render_planning_contents(ui)
                });
                if close_requested {
                    open = false;
                }
            } else {
                let mut close_requested = false;
                crate::egui_compat::show_central_panel(
                    &mut *ctx,
                    egui::CentralPanel::default(),
                    |ui| {
                        close_requested = self.render_planning_contents(ui);
                    },
                );
                if close_requested || Self::viewport_close_requested_or_shortcut(ctx) {
                    open = false;
                }
            }
        });
        if ctx.input(|i| i.key_pressed(Key::Escape)) {
            open = false;
        }
        self.show_planning_dialog = open;
    }

    fn refresh_project_restriction_enzymes(&mut self, resource_path: &str) -> Result<usize> {
        let enzymes = enzymes::load_restriction_enzymes_from_path(resource_path)?;
        let mut engine = self.engine.write().unwrap();
        for dna in engine.state_mut().sequences.values_mut() {
            enzymes.clone_into(dna.restriction_enzymes_mut());
            dna.set_max_restriction_enzyme_sites(None);
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

    fn import_workflow_macro_templates_from_path(&mut self, path: &str) {
        let command = ShellCommand::MacrosTemplateImport {
            path: path.to_string(),
        };
        let options = ShellExecutionOptions::from_env();
        let run = {
            let mut engine = self.engine.write().expect("Engine lock poisoned");
            execute_shell_command_with_options(&mut engine, &command, &options)
        };
        match run {
            Ok(out) => {
                let imported_count = out
                    .output
                    .get("imported_count")
                    .and_then(|v| v.as_u64())
                    .unwrap_or(0);
                let source_file_count = out
                    .output
                    .get("source_files")
                    .and_then(|v| v.as_array())
                    .map(|rows| rows.len())
                    .unwrap_or(0);
                self.app_status = format!(
                    "Imported {imported_count} workflow macro template(s) from '{path}' ({source_file_count} source file(s))"
                );
            }
            Err(err) => {
                self.app_status = format!("Macro template import failed for '{path}': {err}");
            }
        }
    }

    fn prompt_import_workflow_macro_templates_file(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("JSON", &["json"])
            .pick_file()
        {
            self.import_workflow_macro_templates_from_path(&path.display().to_string());
        }
    }

    fn prompt_import_workflow_macro_templates_directory(&mut self) {
        if let Some(path) = rfd::FileDialog::new().pick_folder() {
            self.import_workflow_macro_templates_from_path(&path.display().to_string());
        }
    }

    fn humanize_catalog_label(raw: &str) -> String {
        let normalized = raw.replace(['_', '-'], " ");
        let mut out_words = vec![];
        for word in normalized.split_whitespace() {
            let mut chars = word.chars();
            let Some(first) = chars.next() else {
                continue;
            };
            let mut out = first.to_uppercase().collect::<String>();
            out.push_str(chars.as_str());
            out_words.push(out);
        }
        if out_words.is_empty() {
            raw.to_string()
        } else {
            out_words.join(" ")
        }
    }

    fn collect_cloning_pattern_catalog_entries(
        root: &Path,
    ) -> Result<Vec<CloningPatternCatalogEntry>, String> {
        if !root.exists() {
            return Err(format!("Catalog path '{}' does not exist", root.display()));
        }
        if !root.is_dir() {
            return Err(format!(
                "Catalog path '{}' is not a directory",
                root.display()
            ));
        }
        Self::collect_cloning_pattern_catalog_entries_from_dir(root)
    }

    fn collect_cloning_pattern_catalog_entries_from_dir(
        dir: &Path,
    ) -> Result<Vec<CloningPatternCatalogEntry>, String> {
        let mut out = vec![];
        let mut entries = fs::read_dir(dir)
            .map_err(|e| format!("Could not read catalog directory '{}': {e}", dir.display()))?
            .filter_map(|entry| entry.ok())
            .collect::<Vec<_>>();
        entries.sort_by_key(|entry| entry.file_name().to_string_lossy().to_string());
        for entry in entries {
            let path = entry.path();
            if path.is_dir() {
                let children = Self::collect_cloning_pattern_catalog_entries_from_dir(&path)?;
                if children.is_empty() {
                    continue;
                }
                let label = Self::humanize_catalog_label(&entry.file_name().to_string_lossy());
                out.push(CloningPatternCatalogEntry {
                    label,
                    path: path.display().to_string(),
                    is_file: false,
                    children,
                });
                continue;
            }
            if !path
                .extension()
                .and_then(|ext| ext.to_str())
                .is_some_and(|ext| ext.eq_ignore_ascii_case("json"))
            {
                continue;
            }
            let stem = path
                .file_stem()
                .and_then(|s| s.to_str())
                .map(|s| s.to_string())
                .unwrap_or_else(|| entry.file_name().to_string_lossy().to_string());
            out.push(CloningPatternCatalogEntry {
                label: Self::humanize_catalog_label(&stem),
                path: path.display().to_string(),
                is_file: true,
                children: vec![],
            });
        }
        Ok(out)
    }

    fn render_cloning_pattern_catalog_menu_entries(
        ui: &mut Ui,
        entries: &[CloningPatternCatalogEntry],
        selected_path: &mut Option<String>,
    ) {
        for entry in entries {
            if entry.is_file {
                let response = ui
                    .button(entry.label.clone())
                    .on_hover_text(format!("Import macro template(s) from {}", entry.path));
                if response.clicked() {
                    *selected_path = Some(entry.path.clone());
                }
            } else {
                let label = format!("{}/", entry.label);
                ui.menu_button(label, |ui| {
                    Self::render_cloning_pattern_catalog_menu_entries(
                        ui,
                        &entry.children,
                        selected_path,
                    );
                });
            }
        }
    }

    fn execute_shared_shell_command_json(
        &mut self,
        command: &ShellCommand,
    ) -> std::result::Result<(serde_json::Value, bool), String> {
        let options = ShellExecutionOptions::from_env();
        let run = {
            let mut engine = self.engine.write().expect("Engine lock poisoned");
            execute_shell_command_with_options(&mut engine, command, &options)
        }?;
        if run.state_changed {
            self.lineage_cache_valid = false;
        }
        Ok((run.output, run.state_changed))
    }

    fn render_cloning_routine_menu_entries(
        ui: &mut Ui,
        routines: &[CloningRoutineCatalogRow],
        selected_template_path: &mut Option<String>,
        status_message: &mut Option<String>,
    ) {
        for routine in routines {
            let label = format!("{} [{}]", routine.title, routine.status);
            let tags = if routine.vocabulary_tags.is_empty() {
                "-".to_string()
            } else {
                routine.vocabulary_tags.join(", ")
            };
            let template_path = routine
                .template_path
                .as_deref()
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .unwrap_or("-");
            let details_url = routine
                .details_url
                .as_deref()
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .unwrap_or("-");
            let hover = format!(
                "routine_id: {}\nfamily: {}\nstatus: {}\ntemplate: {}\ntemplate_path: {}\ndetails_url: {}\ninput_ports: {}\noutput_ports: {}\ntags: {}\n{}",
                routine.routine_id,
                routine.family,
                routine.status,
                routine.template_name,
                template_path,
                details_url,
                routine.input_ports.len(),
                routine.output_ports.len(),
                tags,
                routine.summary.as_deref().unwrap_or("No summary")
            );
            let response = ui.button(label).on_hover_text(hover);
            if response.clicked() {
                if let Some(path) = routine
                    .template_path
                    .as_deref()
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                {
                    *selected_template_path = Some(path.to_string());
                    *status_message = Some(format!(
                        "Importing template '{}' for routine '{}'",
                        routine.template_name, routine.routine_id
                    ));
                } else {
                    *status_message = Some(format!(
                        "Routine '{}' uses template '{}' (no template_path specified)",
                        routine.routine_id, routine.template_name
                    ));
                }
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
            .set_file_name(self.default_project_save_file_name())
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
            ProjectAction::OpenTutorialChapter(chapter_id) => {
                self.open_tutorial_project_chapter(&chapter_id)
            }
            ProjectAction::Close => self.close_project(),
            ProjectAction::Quit => self.quit_application(),
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

    fn load_project_from_file_with_recent(&mut self, path: &str, track_recent: bool) -> Result<()> {
        let state = ProjectState::load_from_path(path).map_err(|e| anyhow!(e.to_string()))?;

        self.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
        self.set_current_project_path(path, track_recent);
        self.lineage_cache_valid = false;
        self.lineage_rows.clear();
        self.lineage_edges.clear();
        self.lineage_op_label_by_id.clear();
        self.lineage_containers.clear();
        self.lineage_arrangements.clear();
        self.lineage_racks.clear();
        self.lineage_graph_view = true;
        self.lineage_main_split_fraction = DEFAULT_LINEAGE_MAIN_SPLIT_FRACTION;
        self.lineage_graph_zoom = 1.0;
        self.lineage_graph_area_height = 420.0;
        self.lineage_container_area_height = 220.0;
        self.lineage_container_arrangement_split_fraction =
            DEFAULT_LINEAGE_CONTAINER_ARRANGEMENT_SPLIT_FRACTION;
        self.lineage_graph_scroll_offset = Vec2::ZERO;
        self.lineage_graph_pan_origin = None;
        self.lineage_graph_compact_labels = true;
        self.lineage_main_split_drag_origin = None;
        self.lineage_main_split_drag_start_y = None;
        self.lineage_container_arrangement_split_drag_origin = None;
        self.lineage_container_arrangement_split_drag_start_y = None;
        self.lineage_graph_selected_node_id = None;
        self.lineage_graph_node_offsets.clear();
        self.lineage_graph_drag_origin = None;
        self.lineage_graph_offsets_synced_stamp = 0;
        self.tracked_autosync_last_op_count = None;
        self.tutorial_project_task = None;
        self.tutorial_project_progress_fraction = None;
        self.tutorial_project_progress_label.clear();
        self.tutorial_project_status.clear();
        self.new_windows.clear();
        self.windows.clear();
        self.detached_auxiliary_window_hosts.clear();
        self.windows_to_close.write().unwrap().clear();
        self.pending_focus_viewports.clear();
        self.show_gibson_dialog = false;
        self.show_arrangement_gel_preview_dialog = false;
        self.show_rack_labels_preview_dialog = false;
        self.show_rack_dialog = false;
        self.show_place_arrangement_rack_dialog = false;
        self.show_routine_assistant_dialog = false;
        self.show_jaspar_expert_dialog = false;
        self.gibson_destination_seq_id.clear();
        self.gibson_show_all_unique_cutters = false;
        self.gibson_opening_start_0based.clear();
        self.gibson_opening_end_0based_exclusive.clear();
        self.gibson_insert_seq_id.clear();
        self.gibson_extra_inserts.clear();
        self.gibson_output_id_hint.clear();
        self.gibson_status.clear();
        self.gibson_preview_output = None;
        self.rack_labels_preview = RackLabelsPreviewState::default();
        self.gibson_preview_svg_uri.clear();
        self.arrangement_gel_preview = ArrangementGelPreviewState::default();
        self.rack_view_rack_id.clear();
        self.rack_view_status.clear();
        self.rack_view_scroll_offset = Vec2::ZERO;
        self.rack_view_selected_coordinates.clear();
        self.rack_view_selected_arrangement_ids.clear();
        self.rack_view_drag_state = None;
        self.rack_view_hover_target_coordinate = None;
        self.rack_view_recent_drop_ghost = None;
        self.place_arrangement_source_id.clear();
        self.place_arrangement_target_rack_id.clear();
        self.place_arrangement_status.clear();
        self.routine_assistant_stage = RoutineAssistantStage::GoalAndCandidates;
        self.routine_assistant_goal.clear();
        self.routine_assistant_query.clear();
        self.routine_assistant_candidates.clear();
        self.routine_assistant_selected_routine_id.clear();
        self.routine_assistant_compare_routine_id.clear();
        self.routine_assistant_bindings.clear();
        self.routine_assistant_disambiguation_answers.clear();
        self.routine_assistant_explain_output = None;
        self.routine_assistant_compare_output = None;
        self.routine_assistant_preflight_output = None;
        self.routine_assistant_execute_output = None;
        self.routine_assistant_status.clear();
        self.routine_assistant_decision_trace = None;
        self.routine_assistant_trace_counter = 1;
        self.agent_task = None;
        self.agent_model_discovery_task = None;
        self.agent_status.clear();
        self.agent_preflight_output = None;
        self.agent_last_invocation = None;
        self.agent_execution_log.clear();
        self.agent_discovered_models.clear();
        self.agent_discovered_model_pick.clear();
        self.agent_model_discovery_status.clear();
        self.agent_model_discovery_source_key.clear();
        self.agent_model_discovery_failed_source_key.clear();
        self.load_bed_track_subscriptions_from_state();
        self.load_lineage_graph_workspace_from_state();
        self.load_rack_workspace_from_state();
        self.load_background_job_history_from_state();

        self.mark_clean_snapshot();
        Ok(())
    }

    fn load_project_from_file(&mut self, path: &str) -> Result<()> {
        self.load_project_from_file_with_recent(path, true)
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

    fn current_project_file_stem(&self) -> String {
        let display_name = self.current_project_name();
        let trimmed = display_name.trim();
        let raw_stem = trimmed
            .strip_suffix(".gentle.json")
            .or_else(|| trimmed.strip_suffix(".json"))
            .unwrap_or(trimmed);
        Self::sanitize_file_stem(raw_stem, "project")
    }

    fn default_project_save_file_name(&self) -> String {
        if let Some(path) = self.current_project_path.as_ref()
            && let Some(file_name) = Path::new(path).file_name().and_then(|name| name.to_str())
        {
            let trimmed = file_name.trim();
            if !trimmed.is_empty() {
                return trimmed.to_string();
            }
        }
        format!("{}.gentle.json", self.current_project_file_stem())
    }

    fn default_lineage_svg_file_name(&self) -> String {
        format!("{}.lineage.svg", self.current_project_file_stem())
    }

    fn default_lab_assistant_report_file_name(&self) -> String {
        format!(
            "{}.lab_assistant_report.odt",
            self.current_project_file_stem()
        )
    }

    fn load_lineage_graph_workspace_from_state(&mut self) {
        let (workspace_serialized, legacy_offsets_serialized, legacy_groups_serialized) = {
            let engine = self.engine.read().unwrap();
            let metadata = &engine.state().metadata;
            (
                metadata.get(LINEAGE_GRAPH_WORKSPACE_METADATA_KEY).cloned(),
                metadata.get(LINEAGE_NODE_OFFSETS_METADATA_KEY).cloned(),
                metadata.get(LINEAGE_NODE_GROUPS_METADATA_KEY).cloned(),
            )
        };

        self.lineage_graph_zoom = 1.0;
        self.lineage_main_split_fraction = DEFAULT_LINEAGE_MAIN_SPLIT_FRACTION;
        self.lineage_graph_area_height = 420.0;
        self.lineage_container_area_height = 220.0;
        self.lineage_container_arrangement_split_fraction =
            DEFAULT_LINEAGE_CONTAINER_ARRANGEMENT_SPLIT_FRACTION;
        self.lineage_graph_scroll_offset = Vec2::ZERO;
        self.lineage_graph_pan_origin = None;
        self.lineage_graph_view = true;
        self.lineage_graph_compact_labels = true;
        self.lineage_main_split_drag_origin = None;
        self.lineage_main_split_drag_start_y = None;
        self.lineage_container_arrangement_split_drag_origin = None;
        self.lineage_container_arrangement_split_drag_start_y = None;
        self.lineage_graph_selected_node_id = None;
        self.lineage_graph_node_offsets.clear();
        self.lineage_node_groups.clear();
        if let Some(serialized) = workspace_serialized {
            if let Ok(workspace) =
                serde_json::from_value::<PersistedLineageGraphWorkspace>(serialized)
            {
                self.lineage_graph_view = workspace.graph_view.unwrap_or(true);
                if workspace.zoom.is_finite() {
                    self.lineage_graph_zoom = workspace.zoom.clamp(0.35, 4.0);
                }
                let mut loaded_main_split_fraction = workspace
                    .main_split_fraction
                    .filter(|value| value.is_finite())
                    .map(|value| value.clamp(0.2, 0.9));
                if workspace.graph_area_height.is_finite() {
                    self.lineage_graph_area_height = workspace
                        .graph_area_height
                        .clamp(LINEAGE_MAIN_TOP_PANEL_MIN_HEIGHT, 2400.0);
                }
                if workspace.container_area_height.is_finite() {
                    self.lineage_container_area_height =
                        workspace.container_area_height.clamp(120.0, 1600.0);
                }
                if loaded_main_split_fraction.is_none() {
                    let combined_height =
                        self.lineage_graph_area_height + self.lineage_container_area_height;
                    if combined_height.is_finite() && combined_height > 1.0 {
                        loaded_main_split_fraction = Some(
                            (self.lineage_graph_area_height / combined_height).clamp(0.2, 0.9),
                        );
                    }
                }
                if let Some(main_split_fraction) = loaded_main_split_fraction {
                    self.lineage_main_split_fraction = main_split_fraction;
                }
                if let Some(container_arrangement_split_fraction) = workspace
                    .container_arrangement_split_fraction
                    .filter(|value| value.is_finite())
                    .map(|value| value.clamp(0.2, 0.8))
                {
                    self.lineage_container_arrangement_split_fraction =
                        container_arrangement_split_fraction;
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
                self.lineage_node_groups = workspace.node_groups;
            }
        } else if let Some(serialized) = legacy_offsets_serialized
            && let Ok(raw) = serde_json::from_value::<HashMap<String, [f32; 2]>>(serialized)
        {
            for (node_id, pair) in raw {
                if pair[0].is_finite() && pair[1].is_finite() {
                    self.lineage_graph_node_offsets
                        .insert(node_id, Vec2::new(pair[0], pair[1]));
                }
            }
        }
        if self.lineage_node_groups.is_empty()
            && let Some(serialized) = legacy_groups_serialized
            && let Ok(groups) = serde_json::from_value::<Vec<PersistedLineageNodeGroup>>(serialized)
        {
            self.lineage_node_groups = groups;
        }
        self.lineage_graph_drag_origin = None;
        self.lineage_graph_offsets_synced_stamp = 0;
        self.lineage_group_form_editing_id = None;
        self.lineage_group_form_label.clear();
        self.lineage_group_form_members.clear();
        self.lineage_group_form_representative.clear();
        self.lineage_group_status.clear();
        self.lineage_group_marked_nodes.clear();
        self.lineage_node_rename_target_id = None;
        self.lineage_node_rename_text.clear();
        self.lineage_node_remove_target_id = None;
    }

    fn load_rack_workspace_from_state(&mut self) {
        let workspace_serialized = {
            let engine = self.engine.read().unwrap();
            engine
                .state()
                .metadata
                .get(RACK_WORKSPACE_METADATA_KEY)
                .cloned()
        };

        self.rack_help_strip_collapsed = false;
        self.rack_help_strip_pinned_open = false;
        self.rack_help_strip_successful_move_count = 0;
        self.rack_help_strip_auto_minimized = false;
        if let Some(serialized) = workspace_serialized
            && let Ok(workspace) = serde_json::from_value::<PersistedRackWorkspace>(serialized)
        {
            self.rack_help_strip_collapsed = workspace.help_strip_collapsed;
            self.rack_help_strip_pinned_open = workspace.help_strip_pinned_open;
            self.rack_help_strip_successful_move_count = workspace.help_strip_successful_move_count;
            self.rack_help_strip_auto_minimized = workspace.help_strip_auto_minimized;
        }
    }

    fn persist_lineage_graph_workspace_to_state(&mut self) {
        let mut raw: HashMap<String, [f32; 2]> = HashMap::new();
        for (node_id, offset) in &self.lineage_graph_node_offsets {
            if offset.x.is_finite() && offset.y.is_finite() {
                raw.insert(node_id.clone(), [offset.x, offset.y]);
            }
        }

        let workspace = PersistedLineageGraphWorkspace {
            graph_view: Some(self.lineage_graph_view),
            main_split_fraction: Some(self.lineage_main_split_fraction.clamp(0.2, 0.9)),
            container_arrangement_split_fraction: Some(
                self.lineage_container_arrangement_split_fraction
                    .clamp(0.2, 0.8),
            ),
            zoom: self.lineage_graph_zoom.clamp(0.35, 4.0),
            graph_area_height: self
                .lineage_graph_area_height
                .clamp(LINEAGE_MAIN_TOP_PANEL_MIN_HEIGHT, 2400.0),
            container_area_height: self.lineage_container_area_height.clamp(120.0, 1600.0),
            scroll_offset: [
                self.lineage_graph_scroll_offset.x.max(0.0),
                self.lineage_graph_scroll_offset.y.max(0.0),
            ],
            compact_labels: self.lineage_graph_compact_labels,
            node_offsets: raw.clone(),
            node_groups: self.lineage_node_groups.clone(),
        };

        let workspace_is_default = workspace.graph_view.unwrap_or(true)
            && (workspace.zoom - 1.0).abs() <= 0.0001
            && (workspace
                .main_split_fraction
                .unwrap_or(DEFAULT_LINEAGE_MAIN_SPLIT_FRACTION)
                - DEFAULT_LINEAGE_MAIN_SPLIT_FRACTION)
                .abs()
                <= 0.0001
            && (workspace.graph_area_height - 420.0).abs() <= 0.0001
            && (workspace.container_area_height - 220.0).abs() <= 0.0001
            && (workspace
                .container_arrangement_split_fraction
                .unwrap_or(DEFAULT_LINEAGE_CONTAINER_ARRANGEMENT_SPLIT_FRACTION)
                - DEFAULT_LINEAGE_CONTAINER_ARRANGEMENT_SPLIT_FRACTION)
                .abs()
                <= 0.0001
            && workspace.scroll_offset[0].abs() <= 0.0001
            && workspace.scroll_offset[1].abs() <= 0.0001
            && workspace.compact_labels
            && workspace.node_offsets.is_empty()
            && workspace.node_groups.is_empty();

        let mut engine = self.engine.write().unwrap();
        let state = engine.state_mut();
        if workspace_is_default {
            state.metadata.remove(LINEAGE_GRAPH_WORKSPACE_METADATA_KEY);
        } else if let Ok(value) = serde_json::to_value(&workspace)
            && state.metadata.get(LINEAGE_GRAPH_WORKSPACE_METADATA_KEY) != Some(&value)
        {
            state
                .metadata
                .insert(LINEAGE_GRAPH_WORKSPACE_METADATA_KEY.to_string(), value);
        }

        if raw.is_empty() {
            state.metadata.remove(LINEAGE_NODE_OFFSETS_METADATA_KEY);
        } else if let Ok(value) = serde_json::to_value(raw)
            && state.metadata.get(LINEAGE_NODE_OFFSETS_METADATA_KEY) != Some(&value)
        {
            state
                .metadata
                .insert(LINEAGE_NODE_OFFSETS_METADATA_KEY.to_string(), value);
        }
        if self.lineage_node_groups.is_empty() {
            state.metadata.remove(LINEAGE_NODE_GROUPS_METADATA_KEY);
        } else if let Ok(value) = serde_json::to_value(&self.lineage_node_groups)
            && state.metadata.get(LINEAGE_NODE_GROUPS_METADATA_KEY) != Some(&value)
        {
            state
                .metadata
                .insert(LINEAGE_NODE_GROUPS_METADATA_KEY.to_string(), value);
        }
    }

    fn persist_rack_workspace_to_state(&mut self) {
        let workspace = PersistedRackWorkspace {
            help_strip_collapsed: self.rack_help_strip_collapsed,
            help_strip_pinned_open: self.rack_help_strip_pinned_open,
            help_strip_successful_move_count: self.rack_help_strip_successful_move_count,
            help_strip_auto_minimized: self.rack_help_strip_auto_minimized,
        };
        let mut engine = self.engine.write().unwrap();
        let state = engine.state_mut();
        if !workspace.help_strip_collapsed
            && !workspace.help_strip_pinned_open
            && workspace.help_strip_successful_move_count == 0
            && !workspace.help_strip_auto_minimized
        {
            state.metadata.remove(RACK_WORKSPACE_METADATA_KEY);
        } else if let Ok(value) = serde_json::to_value(&workspace)
            && state.metadata.get(RACK_WORKSPACE_METADATA_KEY) != Some(&value)
        {
            state
                .metadata
                .insert(RACK_WORKSPACE_METADATA_KEY.to_string(), value);
        }
    }

    fn specialist_window_close_hover_text(window_title: &str) -> String {
        format!("Close this {window_title} window (Cmd/Ctrl+W)")
    }

    fn render_specialist_window_nav(&mut self, ui: &mut Ui) {
        let _ = self.render_specialist_window_nav_with_close(ui, None);
    }

    fn render_specialist_window_nav_with_close(
        &mut self,
        ui: &mut Ui,
        close_button: Option<(&str, &str)>,
    ) -> bool {
        let mut close_requested = false;
        ui.horizontal(|ui| {
            if ui
                .button(self.tr("button.help"))
                .on_hover_text("Open GUI help (F1 on Windows/Linux, Cmd+Shift+/ on macOS)")
                .clicked()
            {
                self.open_help_doc(HelpDoc::Gui);
            }
            if ui
                .button(self.tr("button.main"))
                .on_hover_text("Bring the main project window to front")
                .clicked()
            {
                self.queue_focus_viewport(ViewportId::ROOT);
            }
            if let Some((label, hover_text)) = close_button
                && ui.button(label).on_hover_text(hover_text).clicked()
            {
                close_requested = true;
            }
        });
        ui.separator();
        close_requested
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
        if (self.last_native_window_entries != native_window_entries
            || self.last_native_active_window_key != active_window_key)
            && !self.native_window_menu_sync_blocked_by_open_probe()
        {
            let sync_started = Instant::now();
            about::sync_native_open_windows_menu(&native_window_entries, active_window_key);
            let sync_elapsed_ms = sync_started.elapsed().as_millis();
            if sync_elapsed_ms >= WINDOW_OPEN_SLOW_THRESHOLD_MS {
                self.app_status = format!(
                    "Native windows menu sync took {sync_elapsed_ms} ms (entries={})",
                    native_window_entries.len()
                );
            }
            self.last_native_window_entries = native_window_entries;
            self.last_native_active_window_key = active_window_key;
        }
        let history_summary = self
            .engine
            .read()
            .map(|engine| engine.history_summary())
            .unwrap_or_default();
        let undo_count = history_summary.undo_count;
        let redo_count = history_summary.redo_count;
        let history_ops_enabled = !self.has_active_background_jobs();
        egui::MenuBar::new().ui(ui, |ui| {
            ui.menu_button(self.tr("menu.file"), |ui| {
                if ui
                    .button("New Project")
                    .on_hover_text("Create a new empty project state")
                    .clicked()
                {
                    self.request_project_action(ProjectAction::New);
                    ui.close();
                }
                if ui
                    .button("Open Project...")
                    .on_hover_text("Open an existing .gentle.json project")
                    .clicked()
                {
                    self.request_project_action(ProjectAction::Open);
                    ui.close();
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
                        ui.close();
                        return;
                    }

                    if let Some(path) = selected_recent_path {
                        self.request_project_action(ProjectAction::OpenPath(path));
                        ui.close();
                    }
                });
                ui.menu_button("Open Tutorial Project...", |ui| {
                    if let Some(task) = &self.tutorial_project_task {
                        let mut cancel_tutorial_clicked = false;
                        let mut show_jobs_clicked = false;
                        ui.horizontal(|ui| {
                            ui.add(egui::Spinner::new());
                            ui.label(format!(
                                "Opening '{}' ({:.1}s)",
                                task.chapter_title,
                                task.started.elapsed().as_secs_f32()
                            ));
                            if task.cancel_requested.load(Ordering::Relaxed) {
                                ui.small("Cancellation requested...");
                            } else if ui
                                .button("Cancel")
                                .on_hover_text(
                                    "Request cancellation of the running tutorial-project build",
                                )
                                .clicked()
                            {
                                cancel_tutorial_clicked = true;
                            }
                            if ui
                                .button("Jobs")
                                .on_hover_text(
                                    "Open the background-jobs panel to inspect tutorial build progress",
                                )
                                .clicked()
                            {
                                show_jobs_clicked = true;
                            }
                        });
                        if !self.tutorial_project_progress_label.trim().is_empty() {
                            ui.small(self.tutorial_project_progress_label.clone());
                        }
                        if !self.tutorial_project_status.trim().is_empty() {
                            ui.small(self.tutorial_project_status.clone());
                        }
                        if cancel_tutorial_clicked {
                            self.request_tutorial_project_task_cancel("tutorial menu");
                        }
                        if show_jobs_clicked {
                            self.open_background_jobs_panel();
                        }
                        return;
                    }

                    let guided_entries = Self::tutorial_project_guided_walkthrough_entries();
                    let mut selected_guided_tutorial: Option<HelpTutorialDocEntry> = None;
                    if !guided_entries.is_empty() {
                        ui.menu_button(
                            format!("Guided walkthroughs ({})", guided_entries.len()),
                            |ui| {
                                ui.small(
                                    "Documentation tutorials open in Help; executable chapters build project state.",
                                );
                                ui.separator();
                                let mut by_group: BTreeMap<
                                    (usize, String),
                                    Vec<&HelpTutorialDocEntry>,
                                > =
                                    BTreeMap::new();
                                for entry in &guided_entries {
                                    let group = Self::tutorial_audience_group_label(entry);
                                    by_group
                                        .entry((
                                            entry.group_order.unwrap_or(usize::MAX),
                                            group,
                                        ))
                                        .or_default()
                                        .push(entry);
                                }
                                for ((_group_order, group), rows) in by_group {
                                    ui.menu_button(format!("{group} ({})", rows.len()), |ui| {
                                        for entry in rows {
                                            let label = Self::tutorial_display_label(
                                                entry.decimal_id.as_deref(),
                                                None,
                                                &entry.title,
                                            );
                                            if ui
                                                .button(label)
                                                .on_hover_text(format!(
                                                    "Open this tutorial guide in Help\n{}\n{}",
                                                    Self::help_tutorial_review_label(entry),
                                                    entry.summary
                                                ))
                                                .clicked()
                                            {
                                                selected_guided_tutorial = Some(entry.clone());
                                            }
                                        }
                                    });
                                }
                            },
                        );
                        ui.separator();
                    }
                    if let Some(entry) = selected_guided_tutorial {
                        if let Err(err) =
                            self.open_help_tutorial_path(&entry.path, &entry.title, &entry.summary)
                        {
                            self.app_status = err;
                        }
                        ui.close();
                        return;
                    }

                    let entries = match Self::load_tutorial_project_entries() {
                        Ok(entries) => entries,
                        Err(err) => {
                            ui.label("Executable tutorial catalog unavailable");
                            ui.small(err);
                            return;
                        }
                    };
                    if entries.is_empty() {
                        ui.add_enabled(false, egui::Button::new("No tutorial chapters"));
                        return;
                    }

                    let mut selected_chapter_id: Option<String> = None;
                    let mut by_group: BTreeMap<(usize, String), Vec<&TutorialProjectEntry>> =
                        BTreeMap::new();
                    for entry in &entries {
                        let group = entry
                            .group_label
                            .clone()
                            .unwrap_or_else(|| entry.tier.as_str().to_string());
                        by_group
                            .entry((entry.group_order.unwrap_or(usize::MAX), group))
                            .or_default()
                            .push(entry);
                    }
                    for ((_group_order, group), group_entries) in by_group {
                        ui.menu_button(format!("{group} ({})", group_entries.len()), |ui| {
                            for entry in group_entries {
                                let mut label = Self::tutorial_display_label(
                                    entry.decimal_id.as_deref(),
                                    Some(entry.chapter_order),
                                    &entry.chapter_title,
                                );
                                if entry.example.test_mode == ExampleTestMode::Online {
                                    label.push_str(" [online]");
                                }
                                let hover = format!(
                                    "chapter_id: {}\nexample_id: {}\ntier: {}\ngroup: {}\n{}\n{}",
                                    entry.chapter_id,
                                    entry.example.id,
                                    entry.tier.as_str(),
                                    group,
                                    crate::workflow_examples::tutorial_review_badge_label(
                                        entry.review_status.as_deref(),
                                        entry.review_stale,
                                        entry.codex_reviewed_at.as_deref(),
                                        entry.human_reviewed_at.as_deref(),
                                        entry.human_reviewer.as_deref(),
                                    ),
                                    if entry.chapter_summary.trim().is_empty() {
                                        "No summary provided."
                                    } else {
                                        entry.chapter_summary.trim()
                                    }
                                );
                                if ui.button(label).on_hover_text(hover).clicked() {
                                    selected_chapter_id = Some(entry.chapter_id.clone());
                                }
                            }
                        });
                    }

                    if let Some(chapter_id) = selected_chapter_id {
                        self.request_project_action(ProjectAction::OpenTutorialChapter(chapter_id));
                        ui.close();
                    }
                });
                if ui
                    .add_enabled(self.can_close_project(), egui::Button::new("Close Project"))
                    .on_hover_text("Close the current project from the workspace")
                    .clicked()
                {
                    self.request_project_action(ProjectAction::Close);
                    ui.close();
                }
                ui.separator();
                if ui
                    .button("New Sequence...")
                    .on_hover_text("Create a project sequence from typed or pasted IUPAC DNA")
                    .clicked()
                {
                    self.open_new_sequence_dialog();
                    ui.close();
                }
                if ui
                    .button("New Sequence from Clipboard...")
                    .on_hover_text(
                        "Read system clipboard text into the new-sequence dialog for review",
                    )
                    .clicked()
                {
                    self.open_new_sequence_from_clipboard();
                    ui.close();
                }
                ui.separator();
                if ui
                    .button("Open Sequence...")
                    .on_hover_text("Import one or more sequence files into the current project")
                    .clicked()
                {
                    self.prompt_open_sequence();
                    ui.close();
                }
                if ui
                    .button("Protein Evidence...")
                    .on_hover_text(
                        "Fetch/import UniProt or Ensembl protein evidence and compare it against sequences",
                    )
                    .clicked()
                {
                    self.open_uniprot_dialog();
                    ui.close();
                }
                if ui
                    .button("Fetch GenBank / dbSNP...")
                    .on_hover_text(
                        "Fetch one GenBank accession or resolve a dbSNP rsID into a genome region",
                    )
                    .clicked()
                {
                    self.open_genbank_dialog();
                    ui.close();
                }
                ui.separator();
                if ui
                    .button(self.tr("menu.file.configuration"))
                    .on_hover_text("Open global app configuration and graphics defaults")
                    .clicked()
                {
                    self.open_configuration_dialog();
                    ui.close();
                }
                if ui
                    .button(self.tr("menu.file.agent_assistant"))
                    .on_hover_text(
                        "Ask configured agent systems and execute suggested shared-shell commands",
                    )
                    .clicked()
                {
                    self.open_agent_assistant_dialog();
                    ui.close();
                }
                if ui
                    .button("Import REBASE Data...")
                    .on_hover_text("Import restriction-enzyme definitions from REBASE JSON/TXT")
                    .clicked()
                {
                    self.prompt_import_rebase_resource();
                    ui.close();
                }
                if ui
                    .button("Import JASPAR Data...")
                    .on_hover_text("Import TF motif library from JASPAR JSON/TXT")
                    .clicked()
                {
                    self.prompt_import_jaspar_resource();
                    ui.close();
                }
                if ui
                    .button("JASPAR Expert...")
                    .on_hover_text("Inspect local JASPAR motifs, sequence logos, and score distributions")
                    .clicked()
                {
                    self.open_jaspar_expert_dialog();
                    ui.close();
                }
                if ui
                    .button("Save Project...")
                    .on_hover_text("Save current project state to disk")
                    .clicked()
                {
                    self.prompt_save_project();
                    ui.close();
                }
                if ui
                    .button("Export DALG SVG...")
                    .on_hover_text("Export lineage graph as SVG")
                    .clicked()
                {
                    self.prompt_export_lineage_svg();
                    ui.close();
                }
                if ui
                    .button("Export Lab Assistant Report...")
                    .on_hover_text("Export a bench-facing report with lineage graphics as ODT, DOCX, or Markdown")
                    .clicked()
                {
                    self.prompt_export_lab_assistant_report();
                    ui.close();
                }
                ui.separator();
                if ui
                    .button(self.tr("menu.file.quit"))
                    .on_hover_text("Quit GENtle (Cmd/Ctrl+Q)")
                    .clicked()
                {
                    self.request_project_action(ProjectAction::Quit);
                    ui.close();
                }
            });
            ui.menu_button(self.tr("menu.edit"), |ui| {
                let undo_resp = self.track_hover_status(
                    ui.add_enabled(
                        history_ops_enabled && undo_count > 0,
                        egui::Button::new(self.tr("menu.edit.undo")),
                    )
                    .on_hover_text("Undo the most recent operation-level state change"),
                    "Edit > Undo",
                );
                if undo_resp.clicked() {
                    self.undo_last_operation();
                    ui.close();
                }
                let redo_resp = self.track_hover_status(
                    ui.add_enabled(
                        history_ops_enabled && redo_count > 0,
                        egui::Button::new(self.tr("menu.edit.redo")),
                    )
                    .on_hover_text("Redo the most recently undone operation"),
                    "Edit > Redo",
                );
                if redo_resp.clicked() {
                    self.redo_last_operation();
                    ui.close();
                }
                if !history_ops_enabled {
                    ui.small("Undo/redo disabled while background jobs are running.");
                } else {
                    ui.small(format!(
                        "Undo {undo_count} | Redo {redo_count} | limit {}",
                        history_summary.history_limit
                    ));
                }
                ui.separator();
                let palette_resp = self.track_hover_status(
                    ui.button(self.tr("menu.edit.command_palette"))
                        .on_hover_text("Open searchable command palette (Cmd/Ctrl+K)"),
                    "Edit > Command Palette",
                );
                if palette_resp.clicked() {
                    self.open_command_palette_dialog();
                    ui.close();
                }
                let history_resp = self.track_hover_status(
                    ui.button(self.tr("menu.edit.operation_history"))
                        .on_hover_text("Show operation history with undo/redo controls"),
                    "Edit > Operation History",
                );
                if history_resp.clicked() {
                    self.history_ui.show_panel = true;
                    ui.close();
                }
            });
            ui.menu_button(self.tr("menu.settings"), |ui| {
                if ui
                    .button(self.tr("menu.file.configuration"))
                    .on_hover_text("Open global app configuration and graphics defaults")
                    .clicked()
                {
                    self.open_configuration_dialog();
                    ui.close();
                }
            });
            ui.menu_button(self.tr("menu.genome"), |ui| {
                if ui
                    .button("Prepare Reference Genome...")
                    .on_hover_text(
                        "Download/index a reference genome for local extraction and BLAST",
                    )
                    .clicked()
                {
                    self.open_reference_genome_prepare_dialog();
                    ui.close();
                }
                if ui
                    .button("Prepared References...")
                    .on_hover_text("Inspect prepared reference/helper genome installations")
                    .clicked()
                {
                    self.open_reference_genome_inspector_dialog();
                    ui.close();
                }
                if ui
                    .button("Clear Caches...")
                    .on_hover_text(
                        "Inspect and conservatively remove prepared genome/helper cache artifacts",
                    )
                    .clicked()
                {
                    self.open_cache_cleanup_dialog();
                    ui.close();
                }
                if ui
                    .button("Retrieve Genomic Sequence...")
                    .on_hover_text(
                        "Extract anchored region/gene sequence from a prepared reference",
                    )
                    .clicked()
                {
                    self.open_reference_genome_retrieve_dialog();
                    ui.close();
                }
                if ui
                    .button("BLAST Genome Sequence...")
                    .on_hover_text("Run BLAST against prepared reference genome indices")
                    .clicked()
                {
                    self.open_reference_genome_blast_dialog();
                    ui.close();
                }
                if ui
                    .button("Import Genome Track...")
                    .on_hover_text("Import BED/BigWig/VCF tracks onto anchored sequences")
                    .clicked()
                {
                    self.open_genome_bed_track_dialog();
                    ui.close();
                }
                ui.separator();
                if ui
                    .button("Prepare Helper Genome...")
                    .on_hover_text("Prepare helper catalog genomes for extraction and BLAST")
                    .clicked()
                {
                    self.open_helper_genome_prepare_dialog();
                    ui.close();
                }
                if ui
                    .button("Retrieve Helper Sequence...")
                    .on_hover_text("Extract sequence from prepared helper genomes")
                    .clicked()
                {
                    self.open_helper_genome_retrieve_dialog();
                    ui.close();
                }
                if ui
                    .button("BLAST Helper Sequence...")
                    .on_hover_text("Run BLAST against prepared helper genome indices")
                    .clicked()
                {
                    self.open_helper_genome_blast_dialog();
                    ui.close();
                }
            });
            ui.menu_button(self.tr("menu.patterns"), |ui| {
                if ui
                    .button("Import Pattern File...")
                    .on_hover_text("Import workflow macro templates from one JSON pattern file")
                    .clicked()
                {
                    self.prompt_import_workflow_macro_templates_file();
                    ui.close();
                }
                if ui
                    .button("Import Pattern Folder...")
                    .on_hover_text(
                        "Import workflow macro templates from all JSON files in one folder tree",
                    )
                    .clicked()
                {
                    self.prompt_import_workflow_macro_templates_directory();
                    ui.close();
                }
                if ui
                    .button("Import Built-in Legacy Pack")
                    .on_hover_text("Import templates from assets/cloning_patterns.json")
                    .clicked()
                {
                    self.import_workflow_macro_templates_from_path(
                        DEFAULT_CLONING_PATTERN_PACK_PATH,
                    );
                    ui.close();
                }
                ui.separator();

                let catalog_root = Path::new(DEFAULT_CLONING_PATTERN_CATALOG_DIR);
                match Self::collect_cloning_pattern_catalog_entries(catalog_root) {
                    Ok(entries) => {
                        if entries.is_empty() {
                            ui.add_enabled(false, egui::Button::new("Catalog is empty"));
                        } else {
                            if ui
                                .button("Import Full Catalog")
                                .on_hover_text(format!(
                                    "Import all catalog templates from {}",
                                    catalog_root.display()
                                ))
                                .clicked()
                            {
                                self.import_workflow_macro_templates_from_path(
                                    &catalog_root.display().to_string(),
                                );
                                ui.close();
                            }
                            ui.separator();
                            ui.small("Catalog hierarchy");
                            let mut selected_path: Option<String> = None;
                            Self::render_cloning_pattern_catalog_menu_entries(
                                ui,
                                &entries,
                                &mut selected_path,
                            );
                            if let Some(path) = selected_path {
                                self.import_workflow_macro_templates_from_path(&path);
                                ui.close();
                            }
                        }
                    }
                    Err(err) => {
                        ui.add_enabled(false, egui::Button::new("Catalog unavailable"));
                        ui.small(err);
                    }
                }

                ui.separator();
                if ui
                    .button("Gibson...")
                    .on_hover_text(
                        "Open destination-first Gibson specialist window with overlap, primer, and cartoon preview",
                    )
                    .clicked()
                {
                    self.open_gibson_dialog();
                    ui.close();
                }
                if ui
                    .button("PCR Designer...")
                    .on_hover_text(
                        "Open paint-first pair-PCR specialist window (ROI + primer windows + queue)",
                    )
                    .clicked()
                {
                    self.open_pcr_design_dialog();
                    ui.close();
                }
                if ui
                    .button("Sequencing Confirmation...")
                    .on_hover_text(
                        "Open construct-confirmation specialist window for called sequencing reads",
                    )
                    .clicked()
                {
                    self.open_sequencing_confirmation_dialog();
                    ui.close();
                }
                if ui
                    .button("Routine Assistant...")
                    .on_hover_text(
                        "Open staged routine-selection workflow (goal, compare, parameters, preflight, run, export)",
                    )
                    .clicked()
                {
                    self.open_routine_assistant_dialog();
                    ui.close();
                }
                if ui
                    .button("Planning...")
                    .on_hover_text(
                        "Open planning profiles/objective editor and suggestion resolution view",
                    )
                    .clicked()
                {
                    self.open_planning_dialog();
                    ui.close();
                }
                ui.separator();
                ui.small("Routine catalog");
                match self.list_cloning_routines(None, None, None) {
                    Ok(mut routines) => {
                        if routines.is_empty() {
                            ui.add_enabled(false, egui::Button::new("No routines found"));
                        } else {
                            routines.sort_by(|left, right| {
                                left.family
                                    .to_ascii_lowercase()
                                    .cmp(&right.family.to_ascii_lowercase())
                                    .then(
                                        left.title
                                            .to_ascii_lowercase()
                                            .cmp(&right.title.to_ascii_lowercase()),
                                    )
                            });
                            if ui
                                .button("Show Routine Catalog Summary")
                                .on_hover_text(
                                    "Show current routine catalog location and routine count",
                                )
                                .clicked()
                            {
                                self.app_status = format!(
                                    "Loaded {} routine(s) from '{}'",
                                    routines.len(),
                                    DEFAULT_CLONING_ROUTINE_CATALOG_PATH
                                );
                                ui.close();
                            }

                            let mut selected_template_path: Option<String> = None;
                            let mut status_message: Option<String> = None;
                            let mut by_family: BTreeMap<String, Vec<CloningRoutineCatalogRow>> =
                                BTreeMap::new();
                            let mut by_status: BTreeMap<String, Vec<CloningRoutineCatalogRow>> =
                                BTreeMap::new();
                            for routine in routines {
                                by_family
                                    .entry(routine.family.clone())
                                    .or_default()
                                    .push(routine.clone());
                                by_status
                                    .entry(routine.status.clone())
                                    .or_default()
                                    .push(routine);
                            }

                            ui.menu_button("Browse by Family", |ui| {
                                for (family, rows) in &by_family {
                                    ui.menu_button(format!("{family} ({})", rows.len()), |ui| {
                                        Self::render_cloning_routine_menu_entries(
                                            ui,
                                            rows,
                                            &mut selected_template_path,
                                            &mut status_message,
                                        );
                                    });
                                }
                            });
                            ui.menu_button("Browse by Status", |ui| {
                                for (status, rows) in &by_status {
                                    ui.menu_button(format!("{status} ({})", rows.len()), |ui| {
                                        Self::render_cloning_routine_menu_entries(
                                            ui,
                                            rows,
                                            &mut selected_template_path,
                                            &mut status_message,
                                        );
                                    });
                                }
                            });

                            if let Some(path) = selected_template_path {
                                self.import_workflow_macro_templates_from_path(&path);
                                if let Some(message) = status_message {
                                    self.app_status = message;
                                }
                                ui.close();
                            } else if let Some(message) = status_message {
                                self.app_status = message;
                                ui.close();
                            }
                        }
                    }
                    Err(err) => {
                        ui.add_enabled(false, egui::Button::new("Routine catalog unavailable"));
                        ui.small(err);
                    }
                }
            });
            ui.menu_button(self.tr("menu.services"), |ui| {
                if ui
                    .button("External Services...")
                    .on_hover_text(
                        "Inspect provider catalog, validate service requests, and prepare quote handoff bundles",
                    )
                    .clicked()
                {
                    self.open_external_services_dialog();
                    ui.close();
                }
                if ui
                    .button("ClawBio...")
                    .on_hover_text("Send the current sequence/selection context to ClawBio")
                    .clicked()
                {
                    self.open_clawbio_dialog();
                    ui.close();
                }
                ui.small("Shared shell routes: services providers/preflight/quote");
            });
            ui.menu_button(self.tr("menu.windows"), |ui| {
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
                    self.toggle_background_jobs_panel();
                    ui.close();
                }
                let history_panel_resp = self.track_hover_status(
                    ui.button(if self.history_ui.show_panel {
                        "Hide Operation History"
                    } else {
                        "Show Operation History"
                    })
                    .on_hover_text("Toggle operation history panel with undo/redo"),
                    "Window > History Panel",
                );
                if history_panel_resp.clicked() {
                    self.history_ui.show_panel = !self.history_ui.show_panel;
                    ui.close();
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
                        ui.close();
                    }
                }
            });
            ui.menu_button(self.tr("menu.help"), |ui| {
                if ui
                    .button("GUI Manual")
                    .on_hover_text("Open the GUI manual in the built-in help window")
                    .clicked()
                {
                    self.open_help_doc(HelpDoc::Gui);
                    ui.close();
                }
                if ui
                    .button("CLI Manual")
                    .on_hover_text("Open the CLI manual in the built-in help window")
                    .clicked()
                {
                    self.open_help_doc(HelpDoc::Cli);
                    ui.close();
                }
                if ui
                    .button("Agent Interface")
                    .on_hover_text("Open the agent interface guide (CLI, MCP, and Agent Assistant)")
                    .clicked()
                {
                    self.open_help_doc(HelpDoc::AgentInterface);
                    ui.close();
                }
                if ui
                    .button("Reviewer Quickstart")
                    .on_hover_text("Open internal preview guide with known limitations and reviewer walkthrough")
                    .clicked()
                {
                    self.open_help_doc(HelpDoc::ReviewerPreview);
                    ui.close();
                }
                if ui
                    .button("Shell Commands")
                    .on_hover_text("Open shell command reference generated from the glossary")
                    .clicked()
                {
                    self.open_help_doc(HelpDoc::Shell);
                    ui.close();
                }
                ui.menu_button("Tutorials", |ui| {
                    if self.help_tutorial_entries.is_empty() {
                        self.help_tutorial_entries = Self::discover_help_tutorial_entries();
                        self.set_help_tutorial_selected(self.help_tutorial_selected);
                    }
                    if self.help_tutorial_entries.is_empty() {
                        ui.add_enabled(false, egui::Button::new("No tutorial docs found"));
                        return;
                    }
                    let mut selected_tutorial: Option<usize> = None;
                    for (index, entry) in self.help_tutorial_entries.iter().enumerate() {
                        if ui
                            .button(entry.title.as_str())
                            .on_hover_text(format!("Open {}\n{}", entry.title, entry.summary))
                            .clicked()
                        {
                            selected_tutorial = Some(index);
                        }
                    }
                    if let Some(index) = selected_tutorial {
                        self.open_help_tutorial_doc(index);
                        ui.close();
                    }
                });
                ui.separator();
                if ui
                    .button("About GENtle")
                    .on_hover_text("Show version and build information")
                    .clicked()
                {
                    self.show_about_dialog = !about::show_native_about_panel();
                    ui.close();
                }
            });
        });
    }

    fn deferred_window_initial_commands(
        initial_position: Option<Pos2>,
    ) -> Vec<egui::ViewportCommand> {
        let mut commands = Vec::new();
        if let Some(position) = initial_position {
            commands.push(egui::ViewportCommand::OuterPosition(position));
        }
        commands
    }

    fn use_immediate_sequence_viewports() -> bool {
        cfg!(target_os = "macos")
    }

    fn macos_native_child_viewports_override_enabled() -> bool {
        env::var(MACOS_NATIVE_CHILD_VIEWPORTS_ENV)
            .ok()
            .map(|value| {
                matches!(
                    value.trim().to_ascii_lowercase().as_str(),
                    "1" | "true" | "yes" | "on"
                )
            })
            .unwrap_or(false)
    }

    fn should_embed_child_viewports() -> bool {
        cfg!(target_os = "macos") && !Self::macos_native_child_viewports_override_enabled()
    }

    fn configure_platform_viewport_mode(ctx: &egui::Context) {
        if Self::should_embed_child_viewports() {
            if !ctx.embed_viewports() {
                ctx.set_embed_viewports(true);
            }
        } else if ctx.embed_viewports() {
            ctx.set_embed_viewports(false);
        }
    }

    fn sequence_window_accepts_native_close_request() -> bool {
        !cfg!(target_os = "macos")
    }

    fn sequence_viewport_class_has_embedded_shell(class: egui::ViewportClass) -> bool {
        matches!(class, egui::ViewportClass::EmbeddedWindow)
    }

    fn take_window_close_requested(window: &Arc<RwLock<Window>>) -> bool {
        window
            .write()
            .map(|mut guard| guard.take_close_requested())
            .unwrap_or(false)
    }

    fn show_window(
        &mut self,
        ctx: &egui::Context,
        id: ViewportId,
        window: Arc<RwLock<Window>>,
        initial_position: Option<Pos2>,
    ) {
        let windows_to_close = self.windows_to_close.clone();
        let window_title = window
            .read()
            .map(|w| w.name())
            .unwrap_or_else(|_| "GENtle".to_string());
        let render_hosted_sequence_in_foreground = self.viewport_foreground_requested(id);
        let builder = egui::ViewportBuilder::default().with_title(window_title.clone());
        let initial_commands = Self::deferred_window_initial_commands(initial_position);
        if ctx.embed_viewports() {
            let update_result = catch_unwind(AssertUnwindSafe(|| {
                if let Ok(mut w) = window.write() {
                    let mut open = true;
                    let min_size = Vec2::new(820.0, 520.0);
                    let spec = crate::egui_compat::HostedWindowSpec::new(
                        window_title.clone(),
                        // Versioned to drop stale persisted egui Resize state
                        // from pre-fix embedded sequence shells.
                        egui::Id::new(("hosted_sequence_window_v2", id)),
                        Vec2::new(1200.0, 860.0),
                        min_size,
                    )
                    .initial_pos(initial_position)
                    .drag_margin(Vec2::new(
                        EMBEDDED_SEQUENCE_WINDOW_DRAG_MARGIN_X_PX,
                        EMBEDDED_SEQUENCE_WINDOW_DRAG_MARGIN_Y_PX,
                    ))
                    .foreground(render_hosted_sequence_in_foreground)
                    .legacy_layer_id(egui::LayerId::new(egui::Order::Middle, egui::Id::new(id)));
                    crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
                        w.update_embedded(ui);
                    });
                    self.clear_viewport_foreground_request_after_render(id);
                    w.update_auxiliary_windows_only(ctx);
                    if !open {
                        w.take_close_requested();
                        if let Ok(mut to_close) = windows_to_close.write() {
                            to_close.push(id);
                        } else {
                            eprintln!("W GENtleApp: close-queue lock poisoned");
                        }
                    }
                } else {
                    eprintln!("W GENtleApp: window lock poisoned; skipping update");
                }
            }));
            if update_result.is_err() {
                eprintln!("E GENtleApp: recovered from panic while updating window");
            }

            let explicit_close_requested = Self::take_window_close_requested(&window);
            let shortcut_triggered = Self::consume_command_or_ctrl_shortcut(ctx, Key::W);
            if explicit_close_requested || shortcut_triggered {
                if let Ok(mut to_close) = windows_to_close.write() {
                    to_close.push(id);
                } else {
                    eprintln!("W GENtleApp: close-queue lock poisoned");
                }
            }
            return;
        }
        if Self::use_immediate_sequence_viewports() {
            ctx.show_viewport_immediate(id, builder, move |ui, class| {
                let _root_ui_guard = crate::egui_compat::install_current_root_ui(ui);
                if !matches!(
                    class,
                    egui::ViewportClass::Immediate | egui::ViewportClass::EmbeddedWindow
                ) {
                    eprintln!(
                        "W GENtleApp: unexpected viewport class, skipping immediate window update"
                    );
                    return;
                }
                if ui.input(|i| i.viewport().focused.unwrap_or(false)) {
                    report_active_viewport_from_ui(id);
                }

                let update_result = catch_unwind(AssertUnwindSafe(|| {
                    if let Ok(mut w) = window.write() {
                        if Self::sequence_viewport_class_has_embedded_shell(class) {
                            w.update_embedded(ui);
                            w.update_auxiliary_windows_only(ui.ctx());
                        } else {
                            w.update(ui.ctx());
                        }
                    } else {
                        eprintln!("W GENtleApp: window lock poisoned; skipping update");
                    }
                }));
                if update_result.is_err() {
                    eprintln!("E GENtleApp: recovered from panic while updating window");
                }

                let explicit_close_requested = Self::take_window_close_requested(&window);
                let shortcut_triggered = Self::consume_command_or_ctrl_shortcut(ui.ctx(), Key::W);
                let native_close_requested = ui.ctx().input(|i| i.viewport().close_requested());
                if native_close_requested && !Self::sequence_window_accepts_native_close_request() {
                    ui.ctx()
                        .send_viewport_cmd(egui::ViewportCommand::CancelClose);
                }
                if explicit_close_requested
                    || shortcut_triggered
                    || (native_close_requested
                        && Self::sequence_window_accepts_native_close_request())
                {
                    if let Ok(mut to_close) = windows_to_close.write() {
                        to_close.push(id);
                    } else {
                        eprintln!("W GENtleApp: close-queue lock poisoned");
                    }
                }
            });
        } else {
            ctx.show_viewport_deferred(id, builder, move |ctx, class| {
                let _root_ui_guard = crate::egui_compat::install_current_root_ui(ctx);
                if !matches!(
                    class,
                    egui::ViewportClass::Deferred | egui::ViewportClass::EmbeddedWindow
                ) {
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
                let explicit_close_requested = Self::take_window_close_requested(&window);
                let shortcut_triggered = Self::consume_command_or_ctrl_shortcut(ctx, Key::W);
                let native_close_requested = ctx.input(|i| i.viewport().close_requested());
                if native_close_requested && !Self::sequence_window_accepts_native_close_request() {
                    ctx.send_viewport_cmd(egui::ViewportCommand::CancelClose);
                }
                if explicit_close_requested
                    || shortcut_triggered
                    || (native_close_requested
                        && Self::sequence_window_accepts_native_close_request())
                {
                    if let Ok(mut to_close) = windows_to_close.write() {
                        to_close.push(id);
                    } else {
                        eprintln!("W GENtleApp: close-queue lock poisoned");
                    }
                }
            });
        }
        for command in initial_commands {
            ctx.send_viewport_cmd_to(id, command);
        }
    }

    fn refresh_lineage_cache_if_needed(&mut self) {
        let stamp = self.current_lineage_change_stamp();
        if self.lineage_cache_valid && self.lineage_cache_stamp == stamp {
            return;
        }

        let (
            rows,
            lineage_edges,
            op_label_by_id,
            reopenable_gibson_op_ids,
            reopenable_pcr_op_seq_ids,
            containers,
            arrangements,
            racks,
        ) = {
            let engine = self.engine.read().unwrap();
            let state = engine.state();
            let anchor_by_seq: HashMap<String, SequenceGenomeAnchorSummary> = engine
                .list_sequence_genome_anchor_summaries()
                .into_iter()
                .map(|summary| (summary.seq_id.clone(), summary))
                .collect();
            let anchor_provenance_by_seq = Self::latest_genome_anchor_provenance_by_seq(state);
            let mut chromosome_length_cache: HashMap<
                GenomeLengthCacheKey,
                Option<Vec<GenomeChromosomeRecord>>,
            > = HashMap::new();
            let mut op_created_count: HashMap<String, usize> = HashMap::new();
            let mut op_created_ids: HashMap<String, Vec<String>> = HashMap::new();
            let mut op_label_by_id: HashMap<String, String> = HashMap::new();
            let mut individually_rendered_multi_output_ops: HashSet<String> = HashSet::new();
            let mut op_retrieval_by_id: HashMap<String, LineageRetrievalDescriptor> =
                HashMap::new();
            let mut dotplot_op_by_id: HashMap<String, (String, String, Option<String>)> =
                HashMap::new();
            let mut flex_track_op_by_id: HashMap<String, (String, String)> = HashMap::new();
            let mut sequencing_confirmation_op_by_report_id: HashMap<
                String,
                (String, String, Option<String>),
            > = HashMap::new();
            let mut reopenable_gibson_op_ids: HashSet<String> = HashSet::new();
            let mut reopenable_pcr_op_seq_ids: HashMap<String, String> = HashMap::new();
            #[derive(Clone)]
            struct PendingSvgExportRow {
                source_seq_ids: Vec<String>,
                seq_id: String,
                display_name: String,
                origin: String,
                created_by_op: String,
                analysis_kind: Option<LineageAnalysisKind>,
                analysis_artifact_id: Option<String>,
                analysis_reference_seq_id: Option<String>,
                analysis_mode: Option<String>,
            }
            let mut pending_svg_export_rows: Vec<PendingSvgExportRow> = vec![];
            let path_file_name = |path: &str| -> String {
                Path::new(path)
                    .file_name()
                    .and_then(|name| name.to_str())
                    .map(|name| name.to_string())
                    .filter(|name| !name.trim().is_empty())
                    .unwrap_or_else(|| path.to_string())
            };
            for rec in engine.operation_log() {
                op_created_count.insert(rec.result.op_id.clone(), rec.result.created_seq_ids.len());
                op_created_ids.insert(rec.result.op_id.clone(), rec.result.created_seq_ids.clone());
                op_label_by_id.insert(rec.result.op_id.clone(), Self::summarize_operation(&rec.op));
                if matches!(rec.op, Operation::ApplyGibsonAssemblyPlan { .. }) {
                    reopenable_gibson_op_ids.insert(rec.result.op_id.clone());
                    individually_rendered_multi_output_ops.insert(rec.result.op_id.clone());
                }
                match &rec.op {
                    Operation::Pcr { template, .. }
                    | Operation::PcrAdvanced { template, .. }
                    | Operation::PcrMutagenesis { template, .. }
                    | Operation::DesignPrimerPairs { template, .. }
                    | Operation::DesignQpcrAssays { template, .. } => {
                        reopenable_pcr_op_seq_ids
                            .insert(rec.result.op_id.clone(), template.clone());
                    }
                    _ => {}
                }
                if let Some(descriptor) = Self::lineage_retrieval_descriptor_from_operation(&rec.op)
                {
                    op_retrieval_by_id.insert(rec.result.op_id.clone(), descriptor);
                }
                match &rec.op {
                    Operation::ComputeDotplot {
                        seq_id,
                        reference_seq_id,
                        store_as,
                        ..
                    } => {
                        let dotplot_id =
                            Self::lineage_dotplot_id_for_operation(store_as, &rec.result.op_id);
                        dotplot_op_by_id.insert(
                            dotplot_id,
                            (
                                rec.result.op_id.clone(),
                                seq_id.clone(),
                                reference_seq_id.clone(),
                            ),
                        );
                    }
                    Operation::ComputeFlexibilityTrack {
                        seq_id, store_as, ..
                    } => {
                        let track_id =
                            Self::lineage_flex_track_id_for_operation(store_as, &rec.result.op_id);
                        flex_track_op_by_id
                            .insert(track_id, (rec.result.op_id.clone(), seq_id.clone()));
                    }
                    Operation::ConfirmConstructReads {
                        expected_seq_id, ..
                    } => {
                        if let Some(report) = rec.result.sequencing_confirmation_report.as_ref() {
                            sequencing_confirmation_op_by_report_id.insert(
                                report.report_id.clone(),
                                (
                                    rec.result.op_id.clone(),
                                    expected_seq_id.clone(),
                                    report.baseline_seq_id.clone(),
                                ),
                            );
                        }
                    }
                    Operation::RenderSequenceSvg { seq_id, mode, path } => {
                        let mode_label = match mode {
                            RenderSvgMode::Linear => "linear",
                            RenderSvgMode::Circular => "circular",
                        };
                        pending_svg_export_rows.push(PendingSvgExportRow {
                            source_seq_ids: vec![seq_id.clone()],
                            seq_id: seq_id.clone(),
                            display_name: path_file_name(path),
                            origin: "SequenceSvgExport".to_string(),
                            created_by_op: rec.result.op_id.clone(),
                            analysis_kind: None,
                            analysis_artifact_id: None,
                            analysis_reference_seq_id: None,
                            analysis_mode: Some(mode_label.to_string()),
                        });
                    }
                    Operation::RenderDotplotSvg {
                        seq_id,
                        dotplot_id,
                        path,
                        ..
                    } => {
                        let mut source_seq_ids = vec![seq_id.clone()];
                        let reference_seq_id = engine
                            .get_dotplot_view(dotplot_id)
                            .ok()
                            .and_then(|view| view.reference_seq_id);
                        if let Some(reference_seq_id) = reference_seq_id.as_ref()
                            && !reference_seq_id.eq_ignore_ascii_case(seq_id)
                        {
                            source_seq_ids.push(reference_seq_id.clone());
                        }
                        pending_svg_export_rows.push(PendingSvgExportRow {
                            source_seq_ids,
                            seq_id: seq_id.clone(),
                            display_name: path_file_name(path),
                            origin: "DotplotSvgExport".to_string(),
                            created_by_op: rec.result.op_id.clone(),
                            analysis_kind: Some(LineageAnalysisKind::Dotplot),
                            analysis_artifact_id: Some(dotplot_id.clone()),
                            analysis_reference_seq_id: reference_seq_id,
                            analysis_mode: Some("svg_export".to_string()),
                        });
                    }
                    Operation::RenderFeatureExpertSvg {
                        seq_id,
                        target,
                        path,
                    } => {
                        let target_label = target.describe();
                        pending_svg_export_rows.push(PendingSvgExportRow {
                            source_seq_ids: vec![seq_id.clone()],
                            seq_id: seq_id.clone(),
                            display_name: path_file_name(path),
                            origin: "FeatureExpertSvgExport".to_string(),
                            created_by_op: rec.result.op_id.clone(),
                            analysis_kind: None,
                            analysis_artifact_id: None,
                            analysis_reference_seq_id: None,
                            analysis_mode: Some(target_label),
                        });
                    }
                    Operation::RenderIsoformArchitectureSvg {
                        seq_id,
                        panel_id,
                        path,
                        ..
                    } => {
                        pending_svg_export_rows.push(PendingSvgExportRow {
                            source_seq_ids: vec![seq_id.clone()],
                            seq_id: seq_id.clone(),
                            display_name: path_file_name(path),
                            origin: "IsoformArchitectureSvgExport".to_string(),
                            created_by_op: rec.result.op_id.clone(),
                            analysis_kind: None,
                            analysis_artifact_id: None,
                            analysis_reference_seq_id: None,
                            analysis_mode: Some(panel_id.clone()),
                        });
                    }
                    Operation::RenderRnaStructureSvg { seq_id, path } => {
                        pending_svg_export_rows.push(PendingSvgExportRow {
                            source_seq_ids: vec![seq_id.clone()],
                            seq_id: seq_id.clone(),
                            display_name: path_file_name(path),
                            origin: "RnaStructureSvgExport".to_string(),
                            created_by_op: rec.result.op_id.clone(),
                            analysis_kind: None,
                            analysis_artifact_id: None,
                            analysis_reference_seq_id: None,
                            analysis_mode: None,
                        });
                    }
                    Operation::RenderPoolGelSvg {
                        inputs,
                        path,
                        container_ids,
                        arrangement_id,
                        ..
                    } => {
                        let mut source_seq_ids: Vec<String> = inputs.clone();
                        if let Some(container_ids) = container_ids {
                            for container_id in container_ids {
                                if let Some(container) =
                                    state.container_state.containers.get(container_id)
                                {
                                    source_seq_ids.extend(container.members.iter().cloned());
                                }
                            }
                        }
                        if let Some(arrangement_id) = arrangement_id
                            && let Some(arrangement) =
                                state.container_state.arrangements.get(arrangement_id)
                        {
                            for container_id in &arrangement.lane_container_ids {
                                if let Some(container) =
                                    state.container_state.containers.get(container_id)
                                {
                                    source_seq_ids.extend(container.members.iter().cloned());
                                }
                            }
                        }
                        source_seq_ids.retain(|seq_id| !seq_id.trim().is_empty());
                        source_seq_ids.sort();
                        source_seq_ids.dedup();
                        let seq_id = source_seq_ids
                            .first()
                            .cloned()
                            .unwrap_or_else(|| "-".to_string());
                        pending_svg_export_rows.push(PendingSvgExportRow {
                            source_seq_ids,
                            seq_id,
                            display_name: path_file_name(path),
                            origin: "GelSvgExport".to_string(),
                            created_by_op: rec.result.op_id.clone(),
                            analysis_kind: None,
                            analysis_artifact_id: None,
                            analysis_reference_seq_id: None,
                            analysis_mode: Some("serial_gel".to_string()),
                        });
                    }
                    _ => {}
                }
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
                    let created_by_multi_output_op = node
                        .created_by_op
                        .as_ref()
                        .is_some_and(|op| individually_rendered_multi_output_ops.contains(op));
                    let pool_size = if created_by_multi_output_op {
                        1
                    } else {
                        node.created_by_op
                            .as_ref()
                            .and_then(|op| op_created_count.get(op))
                            .cloned()
                            .unwrap_or(1)
                    };
                    let pool_members = if created_by_multi_output_op {
                        vec![node.seq_id.clone()]
                    } else {
                        node.created_by_op
                            .as_ref()
                            .and_then(|op| op_created_ids.get(op))
                            .cloned()
                            .unwrap_or_else(|| vec![node.seq_id.clone()])
                    };
                    let genome_anchor_summary = anchor_by_seq.get(&node.seq_id).cloned();
                    let genome_anchor_display = genome_anchor_summary
                        .as_ref()
                        .map(Self::format_genome_anchor_summary);
                    let is_full_genome_sequence = genome_anchor_summary
                        .as_ref()
                        .map(|anchor| {
                            let provenance = anchor_provenance_by_seq.get(&node.seq_id);
                            Self::sequence_represents_full_genome(
                                anchor,
                                length,
                                provenance,
                                &mut chromosome_length_cache,
                            )
                        })
                        .unwrap_or(false);
                    let retrieval_descriptor = node
                        .created_by_op
                        .as_ref()
                        .and_then(|op_id| op_retrieval_by_id.get(op_id))
                        .cloned();
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
                        genome_anchor_summary,
                        genome_anchor_display,
                        is_full_genome_sequence,
                        retrieval_descriptor,
                        analysis_kind: None,
                        analysis_artifact_id: None,
                        analysis_reference_seq_id: None,
                        analysis_mode: None,
                        analysis_status: None,
                        analysis_point_count: None,
                        analysis_bin_count: None,
                        analysis_read_count: None,
                        analysis_trace_count: None,
                        analysis_target_count: None,
                        analysis_variant_count: None,
                        macro_instance_id: None,
                        macro_routine_id: None,
                        macro_template_name: None,
                        macro_status: None,
                        macro_status_message: None,
                        macro_op_ids: vec![],
                        macro_inputs: vec![],
                        macro_outputs: vec![],
                    }
                })
                .collect();

            let mut lineage_edges: Vec<(String, String, String)> = state
                .lineage
                .edges
                .iter()
                .filter_map(|e| {
                    let from = state.lineage.nodes.get(&e.from_node_id)?.node_id.clone();
                    let to = state.lineage.nodes.get(&e.to_node_id)?.node_id.clone();
                    Some((from, to, e.op_id.clone()))
                })
                .collect();

            for instance in &state.lineage.macro_instances {
                let macro_node_id = format!("macro:{}", instance.macro_instance_id);
                let mut parent_labels: Vec<String> = vec![];
                let summarize_binding = |binding: &LineageMacroPortBinding| -> Option<String> {
                    if binding.values.is_empty() {
                        return None;
                    }
                    Some(format!("{}={}", binding.port_id, binding.values.join(",")))
                };
                let sequence_nodes_for_binding =
                    |binding: &LineageMacroPortBinding| -> Vec<String> {
                        let mut node_ids: Vec<String> = vec![];
                        let mut seen: HashSet<String> = HashSet::new();
                        match binding.kind.trim().to_ascii_lowercase().as_str() {
                            "sequence" => {
                                for seq_id in &binding.values {
                                    if let Some(node_id) = state.lineage.seq_to_node.get(seq_id)
                                        && seen.insert(node_id.clone())
                                    {
                                        node_ids.push(node_id.clone());
                                    }
                                }
                            }
                            "container" => {
                                for container_id in &binding.values {
                                    let Some(container) =
                                        state.container_state.containers.get(container_id)
                                    else {
                                        continue;
                                    };
                                    for seq_id in &container.members {
                                        if let Some(node_id) = state.lineage.seq_to_node.get(seq_id)
                                            && seen.insert(node_id.clone())
                                        {
                                            node_ids.push(node_id.clone());
                                        }
                                    }
                                }
                            }
                            _ => {}
                        }
                        node_ids
                    };

                for binding in &instance.bound_inputs {
                    if let Some(label) = summarize_binding(binding) {
                        parent_labels.push(label);
                    }
                    let from_nodes = sequence_nodes_for_binding(binding);
                    let edge_label_id = format!(
                        "macro:{}:in:{}",
                        instance.macro_instance_id, binding.port_id
                    );
                    op_label_by_id
                        .entry(edge_label_id.clone())
                        .or_insert_with(|| format!("Macro input {}", binding.port_id));
                    for from_node_id in from_nodes {
                        lineage_edges.push((
                            from_node_id,
                            macro_node_id.clone(),
                            edge_label_id.clone(),
                        ));
                    }
                }
                for binding in &instance.bound_outputs {
                    let to_nodes = sequence_nodes_for_binding(binding);
                    let edge_label_id = format!(
                        "macro:{}:out:{}",
                        instance.macro_instance_id, binding.port_id
                    );
                    op_label_by_id
                        .entry(edge_label_id.clone())
                        .or_insert_with(|| format!("Macro output {}", binding.port_id));
                    for to_node_id in to_nodes {
                        lineage_edges.push((
                            macro_node_id.clone(),
                            to_node_id,
                            edge_label_id.clone(),
                        ));
                    }
                }

                if parent_labels.is_empty() {
                    parent_labels.push("-".to_string());
                }

                let status_label = match instance.status {
                    crate::engine::MacroInstanceStatus::Ok => "ok",
                    crate::engine::MacroInstanceStatus::Failed => "failed",
                    crate::engine::MacroInstanceStatus::Cancelled => "cancelled",
                };
                let display_name = instance
                    .routine_title
                    .clone()
                    .or_else(|| instance.routine_id.clone())
                    .or_else(|| instance.template_name.clone())
                    .unwrap_or_else(|| "Macro".to_string());
                out.push(LineageRow {
                    kind: LineageNodeKind::Macro,
                    node_id: macro_node_id,
                    seq_id: instance.macro_instance_id.clone(),
                    display_name,
                    origin: format!("Macro ({status_label})"),
                    created_by_op: instance
                        .expanded_op_ids
                        .first()
                        .cloned()
                        .unwrap_or_else(|| "-".to_string()),
                    created_at: instance.created_at_unix_ms,
                    parents: parent_labels,
                    length: 0,
                    circular: false,
                    pool_size: 0,
                    pool_members: vec![],
                    arrangement_id: None,
                    arrangement_mode: None,
                    lane_container_ids: vec![],
                    ladders: vec![],
                    genome_anchor_summary: None,
                    genome_anchor_display: None,
                    is_full_genome_sequence: false,
                    retrieval_descriptor: None,
                    analysis_kind: None,
                    analysis_artifact_id: None,
                    analysis_reference_seq_id: None,
                    analysis_mode: None,
                    analysis_status: None,
                    analysis_point_count: None,
                    analysis_bin_count: None,
                    analysis_read_count: None,
                    analysis_trace_count: None,
                    analysis_target_count: None,
                    analysis_variant_count: None,
                    macro_instance_id: Some(instance.macro_instance_id.clone()),
                    macro_routine_id: instance.routine_id.clone(),
                    macro_template_name: instance.template_name.clone(),
                    macro_status: Some(status_label.to_string()),
                    macro_status_message: instance.status_message.clone(),
                    macro_op_ids: instance.expanded_op_ids.clone(),
                    macro_inputs: instance.bound_inputs.clone(),
                    macro_outputs: instance.bound_outputs.clone(),
                });
            }

            let dotplot_summaries = engine.list_dotplot_views(None);
            for summary in dotplot_summaries {
                let node_id = format!("analysis:dotplot:{}", summary.dotplot_id);
                let op_binding = dotplot_op_by_id.get(&summary.dotplot_id).cloned();
                let created_by_op = op_binding
                    .as_ref()
                    .map(|(op_id, _, _)| op_id.clone())
                    .unwrap_or_else(|| "-".to_string());
                let edge_op_id = if created_by_op == "-" {
                    format!("analysis:dotplot:{}", summary.dotplot_id)
                } else {
                    created_by_op.clone()
                };
                op_label_by_id
                    .entry(edge_op_id.clone())
                    .or_insert_with(|| format!("Compute dotplot: id={}", summary.dotplot_id));
                let query_seq_id = summary.seq_id.clone();
                let reference_seq_id = summary.reference_seq_id.clone().or_else(|| {
                    op_binding
                        .as_ref()
                        .and_then(|(_, _, reference_seq_id)| reference_seq_id.clone())
                });
                let mut parents = vec![query_seq_id.clone()];
                if let Some(reference_seq_id) = reference_seq_id.as_ref()
                    && !reference_seq_id.eq_ignore_ascii_case(&query_seq_id)
                {
                    parents.push(reference_seq_id.clone());
                }
                out.push(LineageRow {
                    kind: LineageNodeKind::Analysis,
                    node_id: node_id.clone(),
                    seq_id: query_seq_id.clone(),
                    display_name: summary.dotplot_id.clone(),
                    origin: "Dotplot".to_string(),
                    created_by_op,
                    created_at: summary.generated_at_unix_ms,
                    parents,
                    length: 0,
                    circular: false,
                    pool_size: 0,
                    pool_members: vec![],
                    arrangement_id: None,
                    arrangement_mode: None,
                    lane_container_ids: vec![],
                    ladders: vec![],
                    genome_anchor_summary: None,
                    genome_anchor_display: None,
                    is_full_genome_sequence: false,
                    retrieval_descriptor: None,
                    analysis_kind: Some(LineageAnalysisKind::Dotplot),
                    analysis_artifact_id: Some(summary.dotplot_id.clone()),
                    analysis_reference_seq_id: reference_seq_id.clone(),
                    analysis_mode: Some(summary.mode.as_str().to_string()),
                    analysis_status: None,
                    analysis_point_count: Some(summary.point_count),
                    analysis_bin_count: None,
                    analysis_read_count: None,
                    analysis_trace_count: None,
                    analysis_target_count: None,
                    analysis_variant_count: None,
                    macro_instance_id: None,
                    macro_routine_id: None,
                    macro_template_name: None,
                    macro_status: None,
                    macro_status_message: None,
                    macro_op_ids: vec![],
                    macro_inputs: vec![],
                    macro_outputs: vec![],
                });
                let mut seen_sources: HashSet<String> = HashSet::new();
                let source_seq_ids = std::iter::once(query_seq_id)
                    .chain(reference_seq_id)
                    .collect::<Vec<_>>();
                for source_seq_id in source_seq_ids {
                    let Some(source_node_id) = state.lineage.seq_to_node.get(&source_seq_id) else {
                        continue;
                    };
                    if seen_sources.insert(source_node_id.clone()) {
                        lineage_edges.push((
                            source_node_id.clone(),
                            node_id.clone(),
                            edge_op_id.clone(),
                        ));
                    }
                }
            }

            let flexibility_summaries = engine.list_flexibility_tracks(None);
            for summary in flexibility_summaries {
                let node_id = format!("analysis:flex:{}", summary.track_id);
                let op_binding = flex_track_op_by_id.get(&summary.track_id).cloned();
                let created_by_op = op_binding
                    .as_ref()
                    .map(|(op_id, _)| op_id.clone())
                    .unwrap_or_else(|| "-".to_string());
                let edge_op_id = if created_by_op == "-" {
                    format!("analysis:flex:{}", summary.track_id)
                } else {
                    created_by_op.clone()
                };
                op_label_by_id.entry(edge_op_id.clone()).or_insert_with(|| {
                    format!("Compute flexibility track: id={}", summary.track_id)
                });
                let seq_id = summary.seq_id.clone();
                out.push(LineageRow {
                    kind: LineageNodeKind::Analysis,
                    node_id: node_id.clone(),
                    seq_id: seq_id.clone(),
                    display_name: summary.track_id.clone(),
                    origin: "FlexibilityTrack".to_string(),
                    created_by_op,
                    created_at: summary.generated_at_unix_ms,
                    parents: vec![seq_id.clone()],
                    length: 0,
                    circular: false,
                    pool_size: 0,
                    pool_members: vec![],
                    arrangement_id: None,
                    arrangement_mode: None,
                    lane_container_ids: vec![],
                    ladders: vec![],
                    genome_anchor_summary: None,
                    genome_anchor_display: None,
                    is_full_genome_sequence: false,
                    retrieval_descriptor: None,
                    analysis_kind: Some(LineageAnalysisKind::FlexibilityTrack),
                    analysis_artifact_id: Some(summary.track_id.clone()),
                    analysis_reference_seq_id: None,
                    analysis_mode: Some(summary.model.as_str().to_string()),
                    analysis_status: None,
                    analysis_point_count: None,
                    analysis_bin_count: Some(summary.bin_count),
                    analysis_read_count: None,
                    analysis_trace_count: None,
                    analysis_target_count: None,
                    analysis_variant_count: None,
                    macro_instance_id: None,
                    macro_routine_id: None,
                    macro_template_name: None,
                    macro_status: None,
                    macro_status_message: None,
                    macro_op_ids: vec![],
                    macro_inputs: vec![],
                    macro_outputs: vec![],
                });
                if let Some(source_node_id) = state.lineage.seq_to_node.get(&seq_id) {
                    lineage_edges.push((
                        source_node_id.clone(),
                        node_id.clone(),
                        edge_op_id.clone(),
                    ));
                }
            }

            let primer_design_summaries = engine.list_primer_design_reports();
            for summary in primer_design_summaries {
                let node_id = format!("analysis:primer:{}", summary.report_id);
                let created_by_op = summary.op_id.clone().unwrap_or_else(|| "-".to_string());
                let edge_op_id = if created_by_op == "-" {
                    format!("analysis:primer:{}", summary.report_id)
                } else {
                    created_by_op.clone()
                };
                op_label_by_id.entry(edge_op_id.clone()).or_insert_with(|| {
                    format!(
                        "Design primer pairs: template={}, report_id={}",
                        summary.template, summary.report_id
                    )
                });
                let seq_id = summary.template.clone();
                out.push(LineageRow {
                    kind: LineageNodeKind::Analysis,
                    node_id: node_id.clone(),
                    seq_id: seq_id.clone(),
                    display_name: summary.report_id.clone(),
                    origin: "PrimerDesign".to_string(),
                    created_by_op,
                    created_at: summary.generated_at_unix_ms,
                    parents: vec![seq_id.clone()],
                    length: 0,
                    circular: false,
                    pool_size: 0,
                    pool_members: vec![],
                    arrangement_id: None,
                    arrangement_mode: None,
                    lane_container_ids: vec![],
                    ladders: vec![],
                    genome_anchor_summary: None,
                    genome_anchor_display: None,
                    is_full_genome_sequence: false,
                    retrieval_descriptor: None,
                    analysis_kind: Some(LineageAnalysisKind::PrimerDesign),
                    analysis_artifact_id: Some(summary.report_id.clone()),
                    analysis_reference_seq_id: None,
                    analysis_mode: Some(summary.backend_used.clone()),
                    analysis_status: None,
                    analysis_point_count: None,
                    analysis_bin_count: None,
                    analysis_read_count: None,
                    analysis_trace_count: None,
                    analysis_target_count: Some(summary.pair_count),
                    analysis_variant_count: None,
                    macro_instance_id: None,
                    macro_routine_id: None,
                    macro_template_name: None,
                    macro_status: None,
                    macro_status_message: None,
                    macro_op_ids: vec![],
                    macro_inputs: vec![],
                    macro_outputs: vec![],
                });
                if let Some(source_node_id) = state.lineage.seq_to_node.get(&seq_id) {
                    lineage_edges.push((
                        source_node_id.clone(),
                        node_id.clone(),
                        edge_op_id.clone(),
                    ));
                }
            }

            let qpcr_design_summaries = engine.list_qpcr_design_reports();
            for summary in qpcr_design_summaries {
                let node_id = format!("analysis:qpcr:{}", summary.report_id);
                let created_by_op = summary.op_id.clone().unwrap_or_else(|| "-".to_string());
                let edge_op_id = if created_by_op == "-" {
                    format!("analysis:qpcr:{}", summary.report_id)
                } else {
                    created_by_op.clone()
                };
                op_label_by_id.entry(edge_op_id.clone()).or_insert_with(|| {
                    format!(
                        "Design qPCR assays: template={}, report_id={}",
                        summary.template, summary.report_id
                    )
                });
                let seq_id = summary.template.clone();
                out.push(LineageRow {
                    kind: LineageNodeKind::Analysis,
                    node_id: node_id.clone(),
                    seq_id: seq_id.clone(),
                    display_name: summary.report_id.clone(),
                    origin: "QpcrDesign".to_string(),
                    created_by_op,
                    created_at: summary.generated_at_unix_ms,
                    parents: vec![seq_id.clone()],
                    length: 0,
                    circular: false,
                    pool_size: 0,
                    pool_members: vec![],
                    arrangement_id: None,
                    arrangement_mode: None,
                    lane_container_ids: vec![],
                    ladders: vec![],
                    genome_anchor_summary: None,
                    genome_anchor_display: None,
                    is_full_genome_sequence: false,
                    retrieval_descriptor: None,
                    analysis_kind: Some(LineageAnalysisKind::QpcrDesign),
                    analysis_artifact_id: Some(summary.report_id.clone()),
                    analysis_reference_seq_id: None,
                    analysis_mode: Some(summary.backend_used.clone()),
                    analysis_status: None,
                    analysis_point_count: None,
                    analysis_bin_count: None,
                    analysis_read_count: None,
                    analysis_trace_count: None,
                    analysis_target_count: Some(summary.assay_count),
                    analysis_variant_count: None,
                    macro_instance_id: None,
                    macro_routine_id: None,
                    macro_template_name: None,
                    macro_status: None,
                    macro_status_message: None,
                    macro_op_ids: vec![],
                    macro_inputs: vec![],
                    macro_outputs: vec![],
                });
                if let Some(source_node_id) = state.lineage.seq_to_node.get(&seq_id) {
                    lineage_edges.push((
                        source_node_id.clone(),
                        node_id.clone(),
                        edge_op_id.clone(),
                    ));
                }
            }

            let restriction_cloning_summaries = engine.list_restriction_cloning_pcr_handoffs();
            for summary in restriction_cloning_summaries {
                let node_id = format!("analysis:restriction_cloning_pcr:{}", summary.report_id);
                let created_by_op = summary.op_id.clone().unwrap_or_else(|| "-".to_string());
                let edge_op_id = if created_by_op == "-" {
                    format!("analysis:restriction_cloning_pcr:{}", summary.report_id)
                } else {
                    created_by_op.clone()
                };
                op_label_by_id.entry(edge_op_id.clone()).or_insert_with(|| {
                    format!(
                        "Prepare restriction-site cloning PCR handoff: template={}, vector={}, report_id={}",
                        summary.template, summary.destination_vector_seq_id, summary.report_id
                    )
                });
                let seq_id = summary.template.clone();
                out.push(LineageRow {
                    kind: LineageNodeKind::Analysis,
                    node_id: node_id.clone(),
                    seq_id: seq_id.clone(),
                    display_name: summary.report_id.clone(),
                    origin: "RestrictionCloningPcrHandoff".to_string(),
                    created_by_op,
                    created_at: summary.generated_at_unix_ms,
                    parents: vec![seq_id.clone()],
                    length: 0,
                    circular: false,
                    pool_size: 0,
                    pool_members: vec![],
                    arrangement_id: None,
                    arrangement_mode: None,
                    lane_container_ids: vec![],
                    ladders: vec![],
                    genome_anchor_summary: None,
                    genome_anchor_display: None,
                    is_full_genome_sequence: false,
                    retrieval_descriptor: None,
                    analysis_kind: Some(LineageAnalysisKind::RestrictionCloningPcrHandoff),
                    analysis_artifact_id: Some(summary.report_id.clone()),
                    analysis_reference_seq_id: Some(summary.destination_vector_seq_id.clone()),
                    analysis_mode: Some(summary.mode.clone()),
                    analysis_status: Some(summary.compatibility_status.clone()),
                    analysis_point_count: None,
                    analysis_bin_count: None,
                    analysis_read_count: None,
                    analysis_trace_count: None,
                    analysis_target_count: Some(1),
                    analysis_variant_count: None,
                    macro_instance_id: None,
                    macro_routine_id: None,
                    macro_template_name: None,
                    macro_status: None,
                    macro_status_message: None,
                    macro_op_ids: vec![],
                    macro_inputs: vec![],
                    macro_outputs: vec![],
                });
                if let Some(source_node_id) = state.lineage.seq_to_node.get(&seq_id) {
                    lineage_edges.push((
                        source_node_id.clone(),
                        node_id.clone(),
                        edge_op_id.clone(),
                    ));
                }
            }

            let protein_derivation_summaries = engine.list_protein_derivation_reports(None);
            for summary in protein_derivation_summaries {
                let node_id = format!("analysis:protein_derive:{}", summary.report_id);
                let created_by_op = summary.op_id.clone().unwrap_or_else(|| "-".to_string());
                let edge_op_id = if created_by_op == "-" {
                    format!("analysis:protein_derive:{}", summary.report_id)
                } else {
                    created_by_op.clone()
                };
                op_label_by_id.entry(edge_op_id.clone()).or_insert_with(|| {
                    format!(
                        "Derive protein sequences: seq={}, report_id={}",
                        summary.seq_id, summary.report_id
                    )
                });
                let seq_id = summary.seq_id.clone();
                out.push(LineageRow {
                    kind: LineageNodeKind::Analysis,
                    node_id: node_id.clone(),
                    seq_id: seq_id.clone(),
                    display_name: summary.effective_output_prefix.clone(),
                    origin: "ProteinDerivation".to_string(),
                    created_by_op,
                    created_at: summary.generated_at_unix_ms,
                    parents: vec![seq_id.clone()],
                    length: 0,
                    circular: false,
                    pool_size: 0,
                    pool_members: vec![],
                    arrangement_id: None,
                    arrangement_mode: None,
                    lane_container_ids: vec![],
                    ladders: vec![],
                    genome_anchor_summary: None,
                    genome_anchor_display: None,
                    is_full_genome_sequence: false,
                    retrieval_descriptor: None,
                    analysis_kind: Some(LineageAnalysisKind::ProteinDerivation),
                    analysis_artifact_id: Some(summary.report_id.clone()),
                    analysis_reference_seq_id: None,
                    analysis_mode: Some(summary.derivation_mode_summary.clone()),
                    analysis_status: None,
                    analysis_point_count: None,
                    analysis_bin_count: None,
                    analysis_read_count: None,
                    analysis_trace_count: None,
                    analysis_target_count: Some(summary.derived_count),
                    analysis_variant_count: None,
                    macro_instance_id: None,
                    macro_routine_id: None,
                    macro_template_name: None,
                    macro_status: None,
                    macro_status_message: None,
                    macro_op_ids: vec![],
                    macro_inputs: vec![],
                    macro_outputs: vec![],
                });
                if let Some(source_node_id) = state.lineage.seq_to_node.get(&seq_id) {
                    lineage_edges.push((
                        source_node_id.clone(),
                        node_id.clone(),
                        edge_op_id.clone(),
                    ));
                }
            }

            let reverse_translation_summaries = engine.list_reverse_translation_reports(None);
            for summary in reverse_translation_summaries {
                let node_id = format!("analysis:reverse_translate:{}", summary.report_id);
                let created_by_op = summary.op_id.clone().unwrap_or_else(|| "-".to_string());
                let edge_op_id = if created_by_op == "-" {
                    format!("analysis:reverse_translate:{}", summary.report_id)
                } else {
                    created_by_op.clone()
                };
                op_label_by_id.entry(edge_op_id.clone()).or_insert_with(|| {
                    format!(
                        "Reverse translate protein sequence: protein={}, report_id={}",
                        summary.protein_seq_id, summary.report_id
                    )
                });
                let protein_seq_id = summary.protein_seq_id.clone();
                out.push(LineageRow {
                    kind: LineageNodeKind::Analysis,
                    node_id: node_id.clone(),
                    seq_id: protein_seq_id.clone(),
                    display_name: summary.coding_seq_id.clone(),
                    origin: "ReverseTranslation".to_string(),
                    created_by_op,
                    created_at: summary.generated_at_unix_ms,
                    parents: vec![protein_seq_id.clone()],
                    length: 0,
                    circular: false,
                    pool_size: 0,
                    pool_members: vec![],
                    arrangement_id: None,
                    arrangement_mode: None,
                    lane_container_ids: vec![],
                    ladders: vec![],
                    genome_anchor_summary: None,
                    genome_anchor_display: None,
                    is_full_genome_sequence: false,
                    retrieval_descriptor: None,
                    analysis_kind: Some(LineageAnalysisKind::ReverseTranslation),
                    analysis_artifact_id: Some(summary.report_id.clone()),
                    analysis_reference_seq_id: Some(summary.coding_seq_id.clone()),
                    analysis_mode: Some(summary.speed_profile_summary.clone()),
                    analysis_status: Some(summary.diagnostics_summary.clone()),
                    analysis_point_count: Some(summary.translation_table),
                    analysis_bin_count: None,
                    analysis_read_count: None,
                    analysis_trace_count: None,
                    analysis_target_count: Some(summary.coding_length_bp),
                    analysis_variant_count: Some(summary.protein_length_aa),
                    macro_instance_id: None,
                    macro_routine_id: None,
                    macro_template_name: None,
                    macro_status: None,
                    macro_status_message: None,
                    macro_op_ids: vec![],
                    macro_inputs: vec![],
                    macro_outputs: vec![],
                });
                if let Some(source_node_id) = state.lineage.seq_to_node.get(&protein_seq_id) {
                    lineage_edges.push((
                        source_node_id.clone(),
                        node_id.clone(),
                        edge_op_id.clone(),
                    ));
                }
            }

            let construct_reasoning_summaries =
                engine.list_construct_reasoning_graph_summaries(None);
            for summary in construct_reasoning_summaries {
                let node_id = format!("analysis:construct_reasoning:{}", summary.graph_id);
                let created_by_op = summary.op_id.clone().unwrap_or_else(|| "-".to_string());
                let edge_op_id = if created_by_op == "-" {
                    format!("analysis:construct_reasoning:{}", summary.graph_id)
                } else {
                    created_by_op.clone()
                };
                op_label_by_id.entry(edge_op_id.clone()).or_insert_with(|| {
                    format!(
                        "Construct reasoning: seq={}, graph_id={}",
                        summary.seq_id, summary.graph_id
                    )
                });
                let seq_id = summary.seq_id.clone();
                let display_name = if summary.objective_title.trim().is_empty() {
                    if summary.contains_protein_to_dna_handoff {
                        format!("Protein→DNA handoff ({})", summary.graph_id)
                    } else {
                        summary.graph_id.clone()
                    }
                } else {
                    summary.objective_title.clone()
                };
                out.push(LineageRow {
                    kind: LineageNodeKind::Analysis,
                    node_id: node_id.clone(),
                    seq_id: seq_id.clone(),
                    display_name,
                    origin: if summary.contains_protein_to_dna_handoff {
                        "ConstructReasoning / ProteinToDnaHandoff".to_string()
                    } else {
                        "ConstructReasoning".to_string()
                    },
                    created_by_op,
                    created_at: summary.generated_at_unix_ms,
                    parents: vec![seq_id.clone()],
                    length: 0,
                    circular: false,
                    pool_size: 0,
                    pool_members: vec![],
                    arrangement_id: None,
                    arrangement_mode: None,
                    lane_container_ids: vec![],
                    ladders: vec![],
                    genome_anchor_summary: None,
                    genome_anchor_display: None,
                    is_full_genome_sequence: false,
                    retrieval_descriptor: None,
                    analysis_kind: Some(LineageAnalysisKind::ConstructReasoning),
                    analysis_artifact_id: Some(summary.graph_id.clone()),
                    analysis_reference_seq_id: None,
                    analysis_mode: Some(summary.objective_id.clone()),
                    analysis_status: Some(summary.objective_goal.clone()),
                    analysis_point_count: Some(summary.evidence_count),
                    analysis_bin_count: None,
                    analysis_read_count: None,
                    analysis_trace_count: None,
                    analysis_target_count: Some(summary.candidate_count),
                    analysis_variant_count: Some(summary.decision_count),
                    macro_instance_id: None,
                    macro_routine_id: None,
                    macro_template_name: None,
                    macro_status: None,
                    macro_status_message: None,
                    macro_op_ids: vec![],
                    macro_inputs: vec![],
                    macro_outputs: vec![],
                });
                if let Some(source_node_id) = state.lineage.seq_to_node.get(&seq_id) {
                    lineage_edges.push((
                        source_node_id.clone(),
                        node_id.clone(),
                        edge_op_id.clone(),
                    ));
                }
            }

            let rna_read_reports = engine.list_rna_read_reports(None);
            for summary in rna_read_reports {
                let node_id = format!("analysis:rna_reads:{}", summary.report_id);
                let created_by_op = summary.op_id.clone().unwrap_or_else(|| "-".to_string());
                let edge_op_id = if created_by_op == "-" {
                    format!("analysis:rna_reads:{}", summary.report_id)
                } else {
                    created_by_op.clone()
                };
                op_label_by_id.entry(edge_op_id.clone()).or_insert_with(|| {
                    format!(
                        "Interpret RNA reads: seq={}, report_id={}",
                        summary.seq_id, summary.report_id
                    )
                });
                let seq_id = summary.seq_id.clone();
                out.push(LineageRow {
                    kind: LineageNodeKind::Analysis,
                    node_id: node_id.clone(),
                    seq_id: seq_id.clone(),
                    display_name: summary.report_id.clone(),
                    origin: "RnaReadInterpretation".to_string(),
                    created_by_op,
                    created_at: summary.generated_at_unix_ms,
                    parents: vec![seq_id.clone()],
                    length: 0,
                    circular: false,
                    pool_size: 0,
                    pool_members: vec![],
                    arrangement_id: None,
                    arrangement_mode: None,
                    lane_container_ids: vec![],
                    ladders: vec![],
                    genome_anchor_summary: None,
                    genome_anchor_display: None,
                    is_full_genome_sequence: false,
                    retrieval_descriptor: None,
                    analysis_kind: Some(LineageAnalysisKind::RnaReadInterpretation),
                    analysis_artifact_id: Some(summary.report_id.clone()),
                    analysis_reference_seq_id: None,
                    analysis_mode: Some(summary.profile.as_str().to_string()),
                    analysis_status: Some(format!(
                        "{} / {}",
                        summary.report_mode.as_str(),
                        summary.origin_mode.as_str()
                    )),
                    analysis_point_count: Some(summary.read_count_seed_passed),
                    analysis_bin_count: None,
                    analysis_read_count: Some(summary.read_count_total),
                    analysis_trace_count: None,
                    analysis_target_count: Some(summary.target_gene_count),
                    analysis_variant_count: Some(summary.read_count_aligned),
                    macro_instance_id: None,
                    macro_routine_id: None,
                    macro_template_name: None,
                    macro_status: None,
                    macro_status_message: None,
                    macro_op_ids: vec![],
                    macro_inputs: vec![],
                    macro_outputs: vec![],
                });
                if let Some(source_node_id) = state.lineage.seq_to_node.get(&seq_id) {
                    lineage_edges.push((
                        source_node_id.clone(),
                        node_id.clone(),
                        edge_op_id.clone(),
                    ));
                }
            }

            let uniprot_projection_summaries = engine.list_uniprot_genome_projections(None);
            for summary in uniprot_projection_summaries {
                let node_id = format!("analysis:uniprot:{}", summary.projection_id);
                let created_by_op = summary.op_id.clone().unwrap_or_else(|| "-".to_string());
                let edge_op_id = if created_by_op == "-" {
                    format!("analysis:uniprot:{}", summary.projection_id)
                } else {
                    created_by_op.clone()
                };
                op_label_by_id.entry(edge_op_id.clone()).or_insert_with(|| {
                    format!(
                        "Project UniProt to genome: entry={}, seq={}, projection_id={}",
                        summary.entry_id, summary.seq_id, summary.projection_id
                    )
                });
                let seq_id = summary.seq_id.clone();
                out.push(LineageRow {
                    kind: LineageNodeKind::Analysis,
                    node_id: node_id.clone(),
                    seq_id: seq_id.clone(),
                    display_name: summary.projection_id.clone(),
                    origin: "UniprotProjection".to_string(),
                    created_by_op,
                    created_at: summary.created_at_unix_ms,
                    parents: vec![seq_id.clone()],
                    length: 0,
                    circular: false,
                    pool_size: 0,
                    pool_members: vec![],
                    arrangement_id: None,
                    arrangement_mode: None,
                    lane_container_ids: vec![],
                    ladders: vec![],
                    genome_anchor_summary: None,
                    genome_anchor_display: None,
                    is_full_genome_sequence: false,
                    retrieval_descriptor: None,
                    analysis_kind: Some(LineageAnalysisKind::UniprotProjection),
                    analysis_artifact_id: Some(summary.projection_id.clone()),
                    analysis_reference_seq_id: None,
                    analysis_mode: Some(summary.entry_id.clone()),
                    analysis_status: summary
                        .transcript_id_filter
                        .clone()
                        .or_else(|| Some("all_transcripts".to_string())),
                    analysis_point_count: None,
                    analysis_bin_count: None,
                    analysis_read_count: None,
                    analysis_trace_count: None,
                    analysis_target_count: Some(summary.transcript_projection_count),
                    analysis_variant_count: None,
                    macro_instance_id: None,
                    macro_routine_id: None,
                    macro_template_name: None,
                    macro_status: None,
                    macro_status_message: None,
                    macro_op_ids: vec![],
                    macro_inputs: vec![],
                    macro_outputs: vec![],
                });
                if let Some(source_node_id) = state.lineage.seq_to_node.get(&seq_id) {
                    lineage_edges.push((
                        source_node_id.clone(),
                        node_id.clone(),
                        edge_op_id.clone(),
                    ));
                }
            }

            let sequencing_confirmation_reports =
                GentleEngine::sequencing_confirmation_reports_from_state(state);
            for report in sequencing_confirmation_reports {
                let node_id = format!("analysis:seq_confirm:{}", report.report_id);
                let op_binding = sequencing_confirmation_op_by_report_id
                    .get(&report.report_id)
                    .cloned();
                let created_by_op = op_binding
                    .as_ref()
                    .map(|(op_id, _, _)| op_id.clone())
                    .unwrap_or_else(|| "-".to_string());
                let edge_op_id = if created_by_op == "-" {
                    format!("analysis:seq_confirm:{}", report.report_id)
                } else {
                    created_by_op.clone()
                };
                op_label_by_id.entry(edge_op_id.clone()).or_insert_with(|| {
                    format!(
                        "Sequencing confirmation: expected={}, baseline={}, report_id={}",
                        report.expected_seq_id,
                        report.baseline_seq_id.as_deref().unwrap_or("-"),
                        report.report_id
                    )
                });
                let expected_seq_id = report.expected_seq_id.clone();
                let baseline_seq_id = report.baseline_seq_id.clone().or_else(|| {
                    op_binding
                        .as_ref()
                        .and_then(|(_, _, baseline_seq_id)| baseline_seq_id.clone())
                });
                let mut parents = vec![expected_seq_id.clone()];
                if let Some(baseline_seq_id) = baseline_seq_id.as_ref()
                    && !baseline_seq_id.eq_ignore_ascii_case(&expected_seq_id)
                {
                    parents.push(baseline_seq_id.clone());
                }
                out.push(LineageRow {
                    kind: LineageNodeKind::Analysis,
                    node_id: node_id.clone(),
                    seq_id: expected_seq_id.clone(),
                    display_name: report.report_id.clone(),
                    origin: "SequencingConfirmation".to_string(),
                    created_by_op,
                    created_at: report.generated_at_unix_ms,
                    parents,
                    length: 0,
                    circular: false,
                    pool_size: 0,
                    pool_members: vec![],
                    arrangement_id: None,
                    arrangement_mode: None,
                    lane_container_ids: vec![],
                    ladders: vec![],
                    genome_anchor_summary: None,
                    genome_anchor_display: None,
                    is_full_genome_sequence: false,
                    retrieval_descriptor: None,
                    analysis_kind: Some(LineageAnalysisKind::SequencingConfirmation),
                    analysis_artifact_id: Some(report.report_id.clone()),
                    analysis_reference_seq_id: baseline_seq_id.clone(),
                    analysis_mode: None,
                    analysis_status: Some(report.overall_status.as_str().to_string()),
                    analysis_point_count: None,
                    analysis_bin_count: None,
                    analysis_read_count: Some(report.read_seq_ids.len()),
                    analysis_trace_count: Some(report.trace_ids.len()),
                    analysis_target_count: Some(report.targets.len()),
                    analysis_variant_count: Some(report.variants.len()),
                    macro_instance_id: None,
                    macro_routine_id: None,
                    macro_template_name: None,
                    macro_status: None,
                    macro_status_message: None,
                    macro_op_ids: vec![],
                    macro_inputs: vec![],
                    macro_outputs: vec![],
                });
                let mut seen_sources: HashSet<String> = HashSet::new();
                for source_seq_id in std::iter::once(expected_seq_id)
                    .chain(baseline_seq_id)
                    .collect::<Vec<_>>()
                {
                    let Some(source_node_id) = state.lineage.seq_to_node.get(&source_seq_id) else {
                        continue;
                    };
                    if seen_sources.insert(source_node_id.clone()) {
                        lineage_edges.push((
                            source_node_id.clone(),
                            node_id.clone(),
                            edge_op_id.clone(),
                        ));
                    }
                }
            }

            let export_created_at_base = out.iter().map(|row| row.created_at).max().unwrap_or(0);
            for (index, pending) in pending_svg_export_rows.iter().enumerate() {
                let node_id = format!("analysis:export:{}", pending.created_by_op);
                let created_at = export_created_at_base.saturating_add(index as u128 + 1);
                out.push(LineageRow {
                    kind: LineageNodeKind::Analysis,
                    node_id: node_id.clone(),
                    seq_id: pending.seq_id.clone(),
                    display_name: pending.display_name.clone(),
                    origin: pending.origin.clone(),
                    created_by_op: pending.created_by_op.clone(),
                    created_at,
                    parents: pending.source_seq_ids.clone(),
                    length: 0,
                    circular: false,
                    pool_size: 0,
                    pool_members: vec![],
                    arrangement_id: None,
                    arrangement_mode: None,
                    lane_container_ids: vec![],
                    ladders: vec![],
                    genome_anchor_summary: None,
                    genome_anchor_display: None,
                    is_full_genome_sequence: false,
                    retrieval_descriptor: None,
                    analysis_kind: pending.analysis_kind,
                    analysis_artifact_id: pending.analysis_artifact_id.clone(),
                    analysis_reference_seq_id: pending.analysis_reference_seq_id.clone(),
                    analysis_mode: pending.analysis_mode.clone(),
                    analysis_status: None,
                    analysis_point_count: None,
                    analysis_bin_count: None,
                    analysis_read_count: None,
                    analysis_trace_count: None,
                    analysis_target_count: None,
                    analysis_variant_count: None,
                    macro_instance_id: None,
                    macro_routine_id: None,
                    macro_template_name: None,
                    macro_status: None,
                    macro_status_message: None,
                    macro_op_ids: vec![],
                    macro_inputs: vec![],
                    macro_outputs: vec![],
                });
                let mut seen_source_nodes: HashSet<String> = HashSet::new();
                for source_seq_id in &pending.source_seq_ids {
                    let Some(source_node_id) = state.lineage.seq_to_node.get(source_seq_id) else {
                        continue;
                    };
                    if seen_source_nodes.insert(source_node_id.clone()) {
                        lineage_edges.push((
                            source_node_id.clone(),
                            node_id.clone(),
                            pending.created_by_op.clone(),
                        ));
                    }
                }
            }

            out.sort_by(|a, b| {
                a.created_at
                    .cmp(&b.created_at)
                    .then(a.node_id.cmp(&b.node_id))
            });
            let mut containers: Vec<ContainerRow> = state
                .container_state
                .containers
                .iter()
                .map(|(id, c)| ContainerRow {
                    container_id: id.clone(),
                    kind: format!("{:?}", c.kind),
                    declared_contents_exclusive: c.declared_contents_exclusive,
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
                    default_rack_id: arrangement.default_rack_id.clone(),
                })
                .collect();
            arrangements.sort_by(|a, b| a.arrangement_id.cmp(&b.arrangement_id));
            let mut racks: Vec<RackRow> = state
                .container_state
                .racks
                .iter()
                .map(|(id, rack)| RackRow {
                    rack_id: id.clone(),
                    name: rack.name.clone(),
                    profile: rack.profile.kind.as_str().to_string(),
                    occupied_positions: rack.placements.len(),
                    arrangement_ids: {
                        let mut arrangement_ids = rack
                            .placements
                            .iter()
                            .map(|entry| entry.arrangement_id.clone())
                            .collect::<BTreeSet<_>>()
                            .into_iter()
                            .collect::<Vec<_>>();
                        arrangement_ids.sort();
                        arrangement_ids
                    },
                })
                .collect();
            racks.sort_by(|a, b| a.rack_id.cmp(&b.rack_id));
            (
                out,
                lineage_edges,
                op_label_by_id,
                reopenable_gibson_op_ids,
                reopenable_pcr_op_seq_ids,
                containers,
                arrangements,
                racks,
            )
        };

        self.lineage_rows = rows;
        self.lineage_edges = lineage_edges;
        self.lineage_op_label_by_id = op_label_by_id;
        self.lineage_reopenable_gibson_op_ids = reopenable_gibson_op_ids;
        self.lineage_reopenable_pcr_op_seq_ids = reopenable_pcr_op_seq_ids;
        self.lineage_containers = containers;
        self.lineage_arrangements = arrangements;
        self.lineage_racks = racks;
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
                    if let Some(indegree) = indegree_by_node.get_mut(child_id)
                        && *indegree > 0
                    {
                        *indegree -= 1;
                        if *indegree == 0 {
                            ready.push(child_id.clone());
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

    fn normalize_lineage_analysis_id(raw: &str) -> String {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return String::new();
        }
        let mut out = String::with_capacity(trimmed.len());
        for ch in trimmed.chars() {
            if ch.is_ascii_alphanumeric() || matches!(ch, '_' | '-' | '.') {
                out.push(ch);
            } else {
                out.push('_');
            }
        }
        out
    }

    fn lineage_dotplot_id_for_operation(store_as: &Option<String>, op_id: &str) -> String {
        store_as
            .as_deref()
            .map(Self::normalize_lineage_analysis_id)
            .filter(|value| !value.is_empty())
            .unwrap_or_else(|| format!("dotplot_{op_id}"))
    }

    fn lineage_flex_track_id_for_operation(store_as: &Option<String>, op_id: &str) -> String {
        store_as
            .as_deref()
            .map(Self::normalize_lineage_analysis_id)
            .filter(|value| !value.is_empty())
            .unwrap_or_else(|| format!("flex_{op_id}"))
    }

    fn lineage_operation_family(raw: &str) -> &'static str {
        let head = raw.split(':').next().unwrap_or(raw).trim();
        let lower = head.to_ascii_lowercase();
        if lower.contains("reverse complement") {
            return "reverse_complement";
        }
        if lower.contains("digest") {
            return "digest";
        }
        if lower.contains("ligation") {
            return "ligation";
        }
        if lower.contains("gibson") {
            return "gibson";
        }
        if lower.contains("pcr") {
            return "pcr";
        }
        if lower.contains("extract") {
            return "extract";
        }
        if lower.contains("merge") {
            return "merge";
        }
        if lower.contains("filter") {
            return "filter";
        }
        if lower.contains("annotate") {
            return "annotate";
        }
        if lower.contains("dotplot") {
            return "dotplot";
        }
        if lower.contains("flexibility track") {
            return "flexibility";
        }
        if lower.contains("sequencing confirmation") || lower.contains("confirmconstructreads") {
            return "sequencing_confirmation";
        }
        if lower.contains("load") || lower.contains("import") {
            return "import";
        }
        if lower.contains("set parameter") {
            return "set_parameter";
        }
        if lower.contains("branch") {
            return "branch";
        }
        if lower.contains("reverse") {
            return "reverse";
        }
        if lower.contains("complement") {
            return "complement";
        }
        "other"
    }

    fn lineage_operation_symbol_for_family(family: &str) -> &'static str {
        match family {
            "reverse_complement" => "RC",
            "digest" => "D",
            "ligation" => "L",
            "gibson" => "GB",
            "pcr" => "P",
            "extract" => "E",
            "merge" => "M",
            "filter" => "F",
            "annotate" => "A",
            "dotplot" => "DP",
            "flexibility" => "FX",
            "sequencing_confirmation" => "SC",
            "import" => "I",
            "set_parameter" => "S",
            "branch" => "B",
            "reverse" => "R",
            "complement" => "C",
            _ => "O",
        }
    }

    fn lineage_operation_symbol(raw: &str) -> String {
        let family = Self::lineage_operation_family(raw);
        if family != "other" {
            return Self::lineage_operation_symbol_for_family(family).to_string();
        }
        let head = raw.split(':').next().unwrap_or(raw).trim();
        for ch in head.chars() {
            if ch.is_ascii_alphabetic() {
                return ch.to_ascii_uppercase().to_string();
            }
        }
        "O".to_string()
    }

    fn lineage_operation_color_for_family(family: &str) -> egui::Color32 {
        match family {
            "reverse_complement" => egui::Color32::from_rgb(84, 118, 172),
            "digest" => egui::Color32::from_rgb(176, 126, 70),
            "ligation" => egui::Color32::from_rgb(78, 144, 108),
            "gibson" => egui::Color32::from_rgb(62, 142, 124),
            "pcr" => egui::Color32::from_rgb(72, 124, 186),
            "extract" => egui::Color32::from_rgb(153, 121, 74),
            "merge" => egui::Color32::from_rgb(96, 146, 88),
            "filter" => egui::Color32::from_rgb(120, 120, 120),
            "annotate" => egui::Color32::from_rgb(74, 150, 156),
            "dotplot" => egui::Color32::from_rgb(148, 92, 172),
            "flexibility" => egui::Color32::from_rgb(132, 104, 186),
            "sequencing_confirmation" => egui::Color32::from_rgb(74, 126, 176),
            "import" => egui::Color32::from_rgb(111, 131, 170),
            "set_parameter" => egui::Color32::from_rgb(128, 109, 82),
            "branch" => egui::Color32::from_rgb(130, 110, 152),
            "reverse" => egui::Color32::from_rgb(108, 130, 174),
            "complement" => egui::Color32::from_rgb(108, 137, 176),
            _ => egui::Color32::from_rgb(102, 102, 102),
        }
    }

    fn lineage_operation_color(raw: &str) -> egui::Color32 {
        let family = Self::lineage_operation_family(raw);
        Self::lineage_operation_color_for_family(family)
    }

    fn lineage_operation_glyph(op_labels: &[String]) -> (String, egui::Color32) {
        if op_labels.is_empty() {
            return ("O".to_string(), egui::Color32::from_rgb(102, 102, 102));
        }
        let first_symbol = Self::lineage_operation_symbol(&op_labels[0]);
        let first_color = Self::lineage_operation_color(&op_labels[0]);
        let all_same = op_labels
            .iter()
            .all(|label| Self::lineage_operation_symbol(label) == first_symbol);
        if all_same {
            return (first_symbol, first_color);
        }
        let count = op_labels.len().min(99);
        (count.to_string(), egui::Color32::from_rgb(102, 102, 102))
    }

    fn lineage_edge_groups(
        lineage_edges: &[(String, String, String)],
    ) -> Vec<(String, String, Vec<String>)> {
        let mut edge_groups: Vec<(String, String, Vec<String>)> = vec![];
        let mut by_pair: HashMap<(String, String), usize> = HashMap::new();
        for (from_node, to_node, op_id) in lineage_edges {
            let key = (from_node.clone(), to_node.clone());
            if let Some(group_idx) = by_pair.get(&key).copied() {
                let group = &mut edge_groups[group_idx].2;
                if !group.iter().any(|existing| existing == op_id) {
                    group.push(op_id.clone());
                }
                continue;
            }
            by_pair.insert(key, edge_groups.len());
            edge_groups.push((from_node.clone(), to_node.clone(), vec![op_id.clone()]));
        }
        edge_groups
    }

    fn lineage_edge_operation_labels(
        op_ids: &[String],
        op_label_by_id: &HashMap<String, String>,
    ) -> Vec<String> {
        let mut labels: Vec<String> = op_ids
            .iter()
            .map(|op_id| {
                op_label_by_id
                    .get(op_id)
                    .cloned()
                    .unwrap_or_else(|| op_id.clone())
            })
            .collect();
        labels.sort();
        labels.dedup();
        labels
    }

    fn lineage_edge_display_label(op_labels: &[String]) -> String {
        match op_labels.len() {
            0 => "Operation".to_string(),
            1 => op_labels[0].clone(),
            count => format!("{count} operations"),
        }
    }

    fn lineage_edge_accessible_name(op_labels: &[String]) -> String {
        match op_labels.len() {
            0 => "operation".to_string(),
            1 => op_labels[0].clone(),
            count => {
                let preview = op_labels
                    .iter()
                    .take(3)
                    .cloned()
                    .collect::<Vec<_>>()
                    .join(", ");
                if op_labels.len() <= 3 {
                    format!("{count} operations: {preview}")
                } else {
                    format!("{count} operations: {preview}, ...")
                }
            }
        }
    }

    fn lineage_collapsed_group_hidden_op_badges(
        edges: &[(String, String, String)],
        groups: &[PersistedLineageNodeGroup],
        op_label_by_id: &HashMap<String, String>,
    ) -> HashMap<String, LineageGroupHiddenOpBadge> {
        let mut badges: HashMap<String, LineageGroupHiddenOpBadge> = HashMap::new();
        for group in groups {
            if !group.collapsed || group.member_node_ids.is_empty() {
                continue;
            }
            let member_node_ids: HashSet<&str> = group
                .member_node_ids
                .iter()
                .map(|node_id| node_id.as_str())
                .collect();
            let mut seen_op_ids: HashSet<String> = HashSet::new();
            let mut family_counts: HashMap<String, usize> = HashMap::new();
            for (from_node, to_node, op_id) in edges {
                if !member_node_ids.contains(from_node.as_str())
                    && !member_node_ids.contains(to_node.as_str())
                {
                    continue;
                }
                if !seen_op_ids.insert(op_id.clone()) {
                    continue;
                }
                let op_label = op_label_by_id
                    .get(op_id)
                    .cloned()
                    .unwrap_or_else(|| op_id.clone());
                let family = Self::lineage_operation_family(&op_label).to_string();
                *family_counts.entry(family).or_insert(0) += 1;
            }
            if family_counts.is_empty() {
                continue;
            }
            let total_ops = seen_op_ids.len();
            let mut families: Vec<(String, usize)> = family_counts.into_iter().collect();
            families.sort_by(|left, right| right.1.cmp(&left.1).then(left.0.cmp(&right.0)));
            badges.insert(
                group.representative_node_id.clone(),
                LineageGroupHiddenOpBadge {
                    total_ops,
                    families,
                },
            );
        }
        badges
    }

    fn lineage_hidden_op_families_summary(
        badge: &LineageGroupHiddenOpBadge,
        max_families: usize,
    ) -> String {
        if badge.families.is_empty() {
            return "none".to_string();
        }
        let mut entries: Vec<String> = badge
            .families
            .iter()
            .take(max_families)
            .map(|(family, count)| {
                let symbol = Self::lineage_operation_symbol_for_family(family);
                if *count > 1 {
                    format!("{symbol}{count}")
                } else {
                    symbol.to_string()
                }
            })
            .collect();
        if badge.families.len() > max_families {
            entries.push(format!("+{}", badge.families.len() - max_families));
        }
        entries.join(" ")
    }

    fn parse_lineage_group_node_list(raw: &str) -> Vec<String> {
        let mut out = vec![];
        let mut seen: HashSet<String> = HashSet::new();
        for token in raw.split(|c: char| c == ',' || c == ';' || c.is_whitespace()) {
            let trimmed = token.trim();
            if trimmed.is_empty() {
                continue;
            }
            if seen.insert(trimmed.to_string()) {
                out.push(trimmed.to_string());
            }
        }
        out
    }

    fn next_lineage_group_id(groups: &[PersistedLineageNodeGroup]) -> String {
        let existing: HashSet<String> = groups.iter().map(|group| group.group_id.clone()).collect();
        let mut next = 1usize;
        loop {
            let candidate = format!("grp-{next}");
            if !existing.contains(&candidate) {
                return candidate;
            }
            next += 1;
        }
    }

    fn sanitize_lineage_node_groups(
        groups: &[PersistedLineageNodeGroup],
        valid_node_ids: &HashSet<String>,
    ) -> Vec<PersistedLineageNodeGroup> {
        let mut out = vec![];
        let mut used_group_ids: HashSet<String> = HashSet::new();
        let mut used_node_ids: HashSet<String> = HashSet::new();
        for (index, group) in groups.iter().enumerate() {
            let mut group_id = group.group_id.trim().to_string();
            if group_id.is_empty() {
                group_id = format!("grp-{}", index + 1);
            }
            if used_group_ids.contains(&group_id) {
                let base = group_id.clone();
                let mut suffix = 2usize;
                loop {
                    let candidate = format!("{base}-{suffix}");
                    if !used_group_ids.contains(&candidate) {
                        group_id = candidate;
                        break;
                    }
                    suffix += 1;
                }
            }
            used_group_ids.insert(group_id.clone());

            let mut label = group.label.trim().to_string();
            if label.is_empty() {
                label = format!("Group {}", index + 1);
            }

            let mut seen_members: HashSet<String> = HashSet::new();
            let mut members: Vec<String> = vec![];
            for node_id in &group.member_node_ids {
                let trimmed = node_id.trim();
                if trimmed.is_empty() || !valid_node_ids.contains(trimmed) {
                    continue;
                }
                if seen_members.insert(trimmed.to_string()) {
                    members.push(trimmed.to_string());
                }
            }

            let mut representative = group.representative_node_id.trim().to_string();
            if representative.is_empty() || !valid_node_ids.contains(&representative) {
                representative = members.first().cloned().unwrap_or_default();
            }
            if representative.is_empty() {
                continue;
            }
            members.retain(|node_id| node_id != &representative);

            if used_node_ids.contains(&representative) {
                if let Some(index_unused) = members
                    .iter()
                    .position(|node_id| !used_node_ids.contains(node_id))
                {
                    representative = members.remove(index_unused);
                } else {
                    continue;
                }
            }

            members.retain(|node_id| !used_node_ids.contains(node_id));
            if members.is_empty() {
                continue;
            }

            used_node_ids.insert(representative.clone());
            for node_id in &members {
                used_node_ids.insert(node_id.clone());
            }

            out.push(PersistedLineageNodeGroup {
                group_id,
                label,
                representative_node_id: representative,
                member_node_ids: members,
                collapsed: group.collapsed,
            });
        }
        out
    }

    fn lineage_node_group_maps(
        groups: &[PersistedLineageNodeGroup],
    ) -> (
        HashMap<String, PersistedLineageNodeGroup>,
        HashMap<String, String>,
    ) {
        let mut by_group_id: HashMap<String, PersistedLineageNodeGroup> = HashMap::new();
        let mut node_to_group_id: HashMap<String, String> = HashMap::new();
        for group in groups {
            by_group_id.insert(group.group_id.clone(), group.clone());
            node_to_group_id.insert(group.representative_node_id.clone(), group.group_id.clone());
            for node_id in &group.member_node_ids {
                node_to_group_id.insert(node_id.clone(), group.group_id.clone());
            }
        }
        (by_group_id, node_to_group_id)
    }

    fn project_lineage_graph_by_groups(
        rows: &[LineageRow],
        edges: &[(String, String, String)],
        groups: &[PersistedLineageNodeGroup],
    ) -> (Vec<LineageRow>, Vec<(String, String, String)>) {
        let (group_by_id, node_to_group_id) = Self::lineage_node_group_maps(groups);
        let projected_rows: Vec<LineageRow> = rows
            .iter()
            .filter_map(|row| {
                let group = node_to_group_id
                    .get(&row.node_id)
                    .and_then(|group_id| group_by_id.get(group_id));
                if let Some(group) = group
                    && group.collapsed
                    && row.node_id != group.representative_node_id
                {
                    return None;
                }
                Some(row.clone())
            })
            .collect();

        let project_node_id = |node_id: &str| -> String {
            let group = node_to_group_id
                .get(node_id)
                .and_then(|group_id| group_by_id.get(group_id));
            if let Some(group) = group
                && group.collapsed
            {
                return group.representative_node_id.clone();
            }
            node_id.to_string()
        };

        let mut projected_edges: Vec<(String, String, String)> = vec![];
        let mut seen_edges: HashSet<(String, String, String)> = HashSet::new();
        for (from_node, to_node, op_id) in edges {
            let from = project_node_id(from_node);
            let to = project_node_id(to_node);
            if from == to {
                continue;
            }
            let key = (from, to, op_id.clone());
            if seen_edges.insert(key.clone()) {
                projected_edges.push(key);
            }
        }
        (projected_rows, projected_edges)
    }

    fn project_lineage_graph_operation_hubs(
        rows: &[LineageRow],
        edges: &[(String, String, String)],
        op_label_by_id: &HashMap<String, String>,
        hub_op_ids: &HashSet<String>,
    ) -> (
        Vec<LineageRow>,
        Vec<(String, String, String)>,
        HashMap<String, String>,
    ) {
        let valid_node_ids: HashSet<String> = rows.iter().map(|row| row.node_id.clone()).collect();
        let row_by_node: HashMap<String, LineageRow> = rows
            .iter()
            .map(|row| (row.node_id.clone(), row.clone()))
            .collect();
        let mut op_parents: HashMap<String, BTreeSet<String>> = HashMap::new();
        let mut op_children: HashMap<String, BTreeSet<String>> = HashMap::new();
        for (from_node, to_node, op_id) in edges {
            if !hub_op_ids.contains(op_id) {
                continue;
            }
            if !valid_node_ids.contains(from_node) || !valid_node_ids.contains(to_node) {
                continue;
            }
            op_parents
                .entry(op_id.clone())
                .or_default()
                .insert(from_node.clone());
            op_children
                .entry(op_id.clone())
                .or_default()
                .insert(to_node.clone());
        }

        let hubbed_ops: HashSet<String> = hub_op_ids
            .iter()
            .filter(|op_id| {
                op_parents.get(*op_id).map(|rows| rows.len()).unwrap_or(0) >= 2
                    && op_children.get(*op_id).map(|rows| rows.len()).unwrap_or(0) >= 2
            })
            .cloned()
            .collect();
        if hubbed_ops.is_empty() {
            return (rows.to_vec(), edges.to_vec(), op_label_by_id.clone());
        }

        let mut out_rows = rows.to_vec();
        let mut out_edges: Vec<(String, String, String)> = vec![];
        let mut out_op_label_by_id = op_label_by_id.clone();
        let mut seen_edges: HashSet<(String, String, String)> = HashSet::new();

        for op_id in &hubbed_ops {
            let parents = op_parents
                .get(op_id)
                .cloned()
                .unwrap_or_default()
                .into_iter()
                .collect::<Vec<_>>();
            let children = op_children
                .get(op_id)
                .cloned()
                .unwrap_or_default()
                .into_iter()
                .collect::<Vec<_>>();
            let hub_node_id = format!("operation:{op_id}");
            let created_at = parents
                .iter()
                .filter_map(|node_id| row_by_node.get(node_id))
                .map(|row| row.created_at)
                .max()
                .unwrap_or(0);
            out_rows.push(LineageRow {
                kind: LineageNodeKind::Analysis,
                node_id: hub_node_id.clone(),
                seq_id: op_id.clone(),
                display_name: "Gibson cloning".to_string(),
                origin: "OperationHub".to_string(),
                created_by_op: op_id.clone(),
                created_at,
                parents: parents.clone(),
                length: 0,
                circular: false,
                pool_size: 0,
                pool_members: vec![],
                arrangement_id: None,
                arrangement_mode: None,
                lane_container_ids: vec![],
                ladders: vec![],
                genome_anchor_summary: None,
                genome_anchor_display: None,
                is_full_genome_sequence: false,
                retrieval_descriptor: None,
                analysis_kind: None,
                analysis_artifact_id: None,
                analysis_reference_seq_id: None,
                analysis_mode: Some("operation_hub".to_string()),
                analysis_status: None,
                analysis_point_count: None,
                analysis_bin_count: None,
                analysis_read_count: None,
                analysis_trace_count: None,
                analysis_target_count: None,
                analysis_variant_count: None,
                macro_instance_id: None,
                macro_routine_id: None,
                macro_template_name: None,
                macro_status: None,
                macro_status_message: None,
                macro_op_ids: vec![],
                macro_inputs: vec![],
                macro_outputs: vec![],
            });

            let inbound_op_id = format!("{op_id}::hub_in");
            let outbound_op_id = format!("{op_id}::hub_out");
            out_op_label_by_id.insert(inbound_op_id.clone(), "Gibson cloning".to_string());
            out_op_label_by_id.insert(outbound_op_id.clone(), "Gibson cloning".to_string());
            for parent in &parents {
                let key = (parent.clone(), hub_node_id.clone(), inbound_op_id.clone());
                if seen_edges.insert(key.clone()) {
                    out_edges.push(key);
                }
            }
            for child in &children {
                let key = (hub_node_id.clone(), child.clone(), outbound_op_id.clone());
                if seen_edges.insert(key.clone()) {
                    out_edges.push(key);
                }
            }
        }

        for edge in edges {
            if hubbed_ops.contains(&edge.2) {
                continue;
            }
            if seen_edges.insert(edge.clone()) {
                out_edges.push(edge.clone());
            }
        }
        out_rows.sort_by(|a, b| {
            a.created_at
                .cmp(&b.created_at)
                .then(a.node_id.cmp(&b.node_id))
        });
        (out_rows, out_edges, out_op_label_by_id)
    }

    fn current_visible_lineage_graph_for_export(
        &mut self,
    ) -> (
        Vec<LineageRow>,
        Vec<(String, String, String)>,
        HashMap<String, String>,
    ) {
        self.refresh_lineage_cache_if_needed();

        let mut graph_rows = self.lineage_rows.clone();
        let mut graph_edges = self.lineage_edges.clone();
        let mut graph_op_label_by_id = self.lineage_op_label_by_id.clone();
        let seq_node_by_seq_id: HashMap<String, String> = self
            .lineage_rows
            .iter()
            .filter(|row| row.kind == LineageNodeKind::Sequence)
            .map(|row| (row.seq_id.clone(), row.node_id.clone()))
            .collect();
        let container_members_by_id: HashMap<String, Vec<String>> = self
            .lineage_containers
            .iter()
            .map(|row| (row.container_id.clone(), row.members.clone()))
            .collect();
        for arrangement in &self.lineage_arrangements {
            let arrangement_node_id = format!("arr:{}", arrangement.arrangement_id);
            let arrangement_edge_op_id = format!(
                "{}::arrangement:{}",
                arrangement.created_by_op, arrangement.arrangement_id
            );
            let mut source_node_ids: Vec<String> = vec![];
            let mut seen_sources: HashSet<String> = HashSet::new();
            for container_id in &arrangement.lane_container_ids {
                if let Some(members) = container_members_by_id.get(container_id) {
                    for seq_id in members {
                        if let Some(source_node_id) = seq_node_by_seq_id.get(seq_id).cloned()
                            && seen_sources.insert(source_node_id.clone())
                        {
                            source_node_ids.push(source_node_id.clone());
                            graph_edges.push((
                                source_node_id,
                                arrangement_node_id.clone(),
                                arrangement_edge_op_id.clone(),
                            ));
                        }
                    }
                }
            }
            graph_op_label_by_id
                .entry(arrangement_edge_op_id)
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
                genome_anchor_summary: None,
                genome_anchor_display: None,
                is_full_genome_sequence: false,
                retrieval_descriptor: None,
                analysis_kind: None,
                analysis_artifact_id: None,
                analysis_reference_seq_id: None,
                analysis_mode: None,
                analysis_status: None,
                analysis_point_count: None,
                analysis_bin_count: None,
                analysis_read_count: None,
                analysis_trace_count: None,
                analysis_target_count: None,
                analysis_variant_count: None,
                macro_instance_id: None,
                macro_routine_id: None,
                macro_template_name: None,
                macro_status: None,
                macro_status_message: None,
                macro_op_ids: vec![],
                macro_inputs: vec![],
                macro_outputs: vec![],
            });
        }
        graph_rows.sort_by(|a, b| {
            a.created_at
                .cmp(&b.created_at)
                .then(a.node_id.cmp(&b.node_id))
        });

        let valid_lineage_node_ids: HashSet<String> =
            graph_rows.iter().map(|row| row.node_id.clone()).collect();
        let sanitized_groups =
            Self::sanitize_lineage_node_groups(&self.lineage_node_groups, &valid_lineage_node_ids);
        let (projected_graph_rows, projected_graph_edges) =
            Self::project_lineage_graph_by_groups(&graph_rows, &graph_edges, &sanitized_groups);
        Self::project_lineage_graph_operation_hubs(
            &projected_graph_rows,
            &projected_graph_edges,
            &graph_op_label_by_id,
            &self.lineage_reopenable_gibson_op_ids,
        )
    }

    fn lineage_svg_node_kind(row: &LineageRow) -> LineageSvgNodeKind {
        match row.kind {
            LineageNodeKind::Sequence if row.pool_size > 1 => LineageSvgNodeKind::Pool,
            LineageNodeKind::Sequence => LineageSvgNodeKind::Sequence,
            LineageNodeKind::Arrangement => LineageSvgNodeKind::Arrangement,
            LineageNodeKind::Macro => LineageSvgNodeKind::Macro,
            LineageNodeKind::Analysis if Self::is_lineage_operation_hub(row) => {
                LineageSvgNodeKind::OperationHub
            }
            LineageNodeKind::Analysis => LineageSvgNodeKind::Analysis,
        }
    }

    fn lineage_svg_node_title(row: &LineageRow) -> String {
        let display = row.display_name.trim();
        if !display.is_empty() {
            display.to_string()
        } else {
            row.seq_id.clone()
        }
    }

    fn lineage_svg_node_subtitle(row: &LineageRow) -> String {
        match row.kind {
            LineageNodeKind::Sequence if row.pool_size > 1 => {
                format!(
                    "{} | pool n={} | {} bp",
                    row.seq_id, row.pool_size, row.length
                )
            }
            LineageNodeKind::Sequence => {
                let topology = if row.circular { "circular" } else { "linear" };
                format!("{} ({} bp, {topology})", row.seq_id, row.length)
            }
            LineageNodeKind::Arrangement => format!(
                "{} | {} | lanes={}",
                row.arrangement_id.as_deref().unwrap_or(&row.seq_id),
                row.arrangement_mode.as_deref().unwrap_or("-"),
                row.lane_container_ids.len()
            ),
            LineageNodeKind::Macro => format!(
                "{} | ops={}",
                row.macro_instance_id.as_deref().unwrap_or(&row.seq_id),
                row.macro_op_ids.len()
            ),
            LineageNodeKind::Analysis if Self::is_lineage_operation_hub(row) => {
                format!("op={}", row.created_by_op)
            }
            LineageNodeKind::Analysis => {
                if let Some(artifact_id) = row.analysis_artifact_id.as_deref() {
                    format!("{} | {}", row.origin, artifact_id)
                } else {
                    format!("{} | {}", row.origin, row.seq_id)
                }
            }
        }
    }

    fn lineage_node_kind_label(kind: LineageNodeKind) -> &'static str {
        match kind {
            LineageNodeKind::Sequence => "sequence",
            LineageNodeKind::Arrangement => "arrangement",
            LineageNodeKind::Macro => "macro",
            LineageNodeKind::Analysis => "analysis",
        }
    }

    fn lineage_row_primary_id(row: &LineageRow) -> &str {
        match row.kind {
            LineageNodeKind::Sequence => &row.seq_id,
            LineageNodeKind::Arrangement => row.arrangement_id.as_deref().unwrap_or(&row.seq_id),
            LineageNodeKind::Macro => row.macro_instance_id.as_deref().unwrap_or(&row.seq_id),
            LineageNodeKind::Analysis => row.analysis_artifact_id.as_deref().unwrap_or(&row.seq_id),
        }
    }

    fn lineage_row_summary_text(row: &LineageRow) -> String {
        let mut lines = vec![
            format!("node_id: {}", row.node_id),
            format!("kind: {}", Self::lineage_node_kind_label(row.kind)),
            format!("primary_id: {}", Self::lineage_row_primary_id(row)),
            format!("display_name: {}", row.display_name),
            format!("seq_id: {}", row.seq_id),
            format!("origin: {}", row.origin),
            format!("created_by_op: {}", row.created_by_op),
            format!("created_at_unix_ms: {}", row.created_at),
        ];

        if !row.parents.is_empty() {
            lines.push(format!("parents: {}", row.parents.join(" + ")));
        }

        match row.kind {
            LineageNodeKind::Sequence => {
                lines.push(format!("length_bp: {}", row.length));
                lines.push(format!(
                    "topology: {}",
                    if row.circular { "circular" } else { "linear" }
                ));
                if row.pool_size > 1 {
                    lines.push(format!("pool_size: {}", row.pool_size));
                    if !row.pool_members.is_empty() {
                        lines.push(format!("pool_members: {}", row.pool_members.join(", ")));
                    }
                }
                if let Some(anchor) = row.genome_anchor_display.as_deref() {
                    lines.push(format!("genome_anchor: {anchor}"));
                }
                if row.is_full_genome_sequence {
                    lines.push("full_genome_sequence: true".to_string());
                }
                if let Some(descriptor) = row.retrieval_descriptor.as_ref() {
                    lines.push(format!(
                        "retrieval_pattern: {}",
                        Self::lineage_retrieval_pattern_label(descriptor)
                    ));
                }
            }
            LineageNodeKind::Arrangement => {
                if let Some(arrangement_id) = row.arrangement_id.as_deref() {
                    lines.push(format!("arrangement_id: {arrangement_id}"));
                }
                if let Some(mode) = row.arrangement_mode.as_deref() {
                    lines.push(format!("arrangement_mode: {mode}"));
                }
                if !row.lane_container_ids.is_empty() {
                    lines.push(format!(
                        "lane_container_ids: {}",
                        row.lane_container_ids.join(", ")
                    ));
                }
                if !row.ladders.is_empty() {
                    lines.push(format!("ladders: {}", row.ladders.join(", ")));
                }
            }
            LineageNodeKind::Macro => {
                if let Some(instance_id) = row.macro_instance_id.as_deref() {
                    lines.push(format!("macro_instance_id: {instance_id}"));
                }
                if let Some(routine_id) = row.macro_routine_id.as_deref() {
                    lines.push(format!("macro_routine_id: {routine_id}"));
                }
                if let Some(template_name) = row.macro_template_name.as_deref() {
                    lines.push(format!("macro_template_name: {template_name}"));
                }
                if let Some(status) = row.macro_status.as_deref() {
                    lines.push(format!("macro_status: {status}"));
                }
                if let Some(status_message) = row.macro_status_message.as_deref() {
                    lines.push(format!("macro_status_message: {status_message}"));
                }
                if !row.macro_op_ids.is_empty() {
                    lines.push(format!("macro_op_ids: {}", row.macro_op_ids.join(", ")));
                }
            }
            LineageNodeKind::Analysis => {
                if let Some(kind) = row.analysis_kind {
                    lines.push(format!("analysis_kind: {}", kind.as_str()));
                }
                if let Some(artifact_id) = row.analysis_artifact_id.as_deref() {
                    lines.push(format!("analysis_artifact_id: {artifact_id}"));
                }
                if let Some(reference_seq_id) = row.analysis_reference_seq_id.as_deref() {
                    lines.push(format!("analysis_reference_seq_id: {reference_seq_id}"));
                }
                if let Some(mode) = row.analysis_mode.as_deref() {
                    lines.push(format!("analysis_mode: {mode}"));
                }
                if let Some(status) = row.analysis_status.as_deref() {
                    lines.push(format!("analysis_status: {status}"));
                }
                if let Some(count) = row.analysis_target_count {
                    lines.push(format!("analysis_target_count: {count}"));
                }
                if let Some(count) = row.analysis_variant_count {
                    lines.push(format!("analysis_variant_count: {count}"));
                }
            }
        }

        lines.join("\n")
    }

    fn lineage_copy_payload(row: &LineageRow, kind: LineageCopyPayloadKind) -> String {
        match kind {
            LineageCopyPayloadKind::NodeId => row.node_id.clone(),
            LineageCopyPayloadKind::PrimaryId => Self::lineage_row_primary_id(row).to_string(),
            LineageCopyPayloadKind::DisplayLabel => {
                let display_name = row.display_name.trim();
                if display_name.is_empty() {
                    Self::lineage_row_primary_id(row).to_string()
                } else {
                    display_name.to_string()
                }
            }
            LineageCopyPayloadKind::RowSummary => Self::lineage_row_summary_text(row),
        }
    }

    fn copy_lineage_row_payload_to_clipboard(
        &mut self,
        ctx: &egui::Context,
        row: &LineageRow,
        kind: LineageCopyPayloadKind,
    ) {
        ctx.copy_text(Self::lineage_copy_payload(row, kind));
        self.lineage_group_status = format!(
            "Copied {} for lineage node '{}'",
            kind.status_label(),
            row.node_id
        );
    }

    fn copy_lineage_node_id_to_clipboard(&mut self, ctx: &egui::Context, node_id: &str) {
        ctx.copy_text(node_id.to_string());
        self.lineage_group_status = format!("Copied node ID for lineage node '{node_id}'");
    }

    fn render_projected_lineage_svg_text(
        &self,
        rows: &[LineageRow],
        edges: &[(String, String, String)],
        op_label_by_id: &HashMap<String, String>,
    ) -> String {
        let (layout_by_node, _, _) = Self::compute_lineage_dag_layout(rows, edges);
        let nodes: Vec<LineageSvgNode> = rows
            .iter()
            .enumerate()
            .map(|(fallback_rank, row)| {
                let (layer, rank) = layout_by_node
                    .get(&row.node_id)
                    .copied()
                    .unwrap_or((0, fallback_rank));
                let manual_offset = self
                    .lineage_graph_node_offsets
                    .get(&row.node_id)
                    .copied()
                    .unwrap_or(Vec2::ZERO);
                LineageSvgNode {
                    node_id: row.node_id.clone(),
                    title: Self::lineage_svg_node_title(row),
                    subtitle: Self::lineage_svg_node_subtitle(row),
                    kind: Self::lineage_svg_node_kind(row),
                    x: 120.0 + layer as f32 * 220.0 + manual_offset.x,
                    y: 120.0 + rank as f32 * 110.0 + manual_offset.y,
                }
            })
            .collect();
        let svg_edges: Vec<LineageSvgEdge> = edges
            .iter()
            .map(|(from_node_id, to_node_id, op_id)| LineageSvgEdge {
                from_node_id: from_node_id.clone(),
                to_node_id: to_node_id.clone(),
                label: op_label_by_id
                    .get(op_id)
                    .cloned()
                    .unwrap_or_else(|| op_id.clone()),
            })
            .collect();
        let title = format!("GENtle Lineage (DALG) - {}", self.current_project_name());
        export_projected_lineage_svg(&title, &nodes, &svg_edges)
    }

    fn render_current_visible_lineage_svg_text(&mut self) -> String {
        let (rows, edges, op_label_by_id) = self.current_visible_lineage_graph_for_export();
        self.render_projected_lineage_svg_text(&rows, &edges, &op_label_by_id)
    }

    fn export_visible_lineage_svg_to_path(
        &mut self,
        path: &str,
    ) -> std::result::Result<String, String> {
        let svg = self.render_current_visible_lineage_svg_text();
        fs::write(path, svg)
            .map_err(|err| format!("Could not export lineage SVG to '{}': {err}", path))?;
        Ok(format!("Wrote lineage SVG to '{}'", path))
    }

    fn build_lineage_table_entries(
        rows: &[LineageRow],
        groups: &[PersistedLineageNodeGroup],
    ) -> Vec<LineageTableEntry> {
        let (group_by_id, node_to_group_id) = Self::lineage_node_group_maps(groups);
        let row_by_node: HashMap<String, LineageRow> = rows
            .iter()
            .map(|row| (row.node_id.clone(), row.clone()))
            .collect();
        let row_order: HashMap<String, usize> = rows
            .iter()
            .enumerate()
            .map(|(index, row)| (row.node_id.clone(), index))
            .collect();

        let mut emitted: HashSet<String> = HashSet::new();
        let mut entries: Vec<LineageTableEntry> = vec![];

        for row in rows {
            if emitted.contains(&row.node_id) {
                continue;
            }

            let Some(group_id) = node_to_group_id.get(&row.node_id).cloned() else {
                entries.push(LineageTableEntry {
                    row: row.clone(),
                    indent_level: 0,
                    group_id: None,
                    group_label: None,
                    is_group_representative: false,
                    hidden_group_member_count: 0,
                });
                emitted.insert(row.node_id.clone());
                continue;
            };

            let Some(group) = group_by_id.get(&group_id) else {
                entries.push(LineageTableEntry {
                    row: row.clone(),
                    indent_level: 0,
                    group_id: None,
                    group_label: None,
                    is_group_representative: false,
                    hidden_group_member_count: 0,
                });
                emitted.insert(row.node_id.clone());
                continue;
            };

            if row.node_id != group.representative_node_id {
                continue;
            }

            entries.push(LineageTableEntry {
                row: row.clone(),
                indent_level: 0,
                group_id: Some(group.group_id.clone()),
                group_label: Some(group.label.clone()),
                is_group_representative: true,
                hidden_group_member_count: if group.collapsed {
                    group.member_node_ids.len()
                } else {
                    0
                },
            });
            emitted.insert(row.node_id.clone());

            if group.collapsed {
                for member_id in &group.member_node_ids {
                    emitted.insert(member_id.clone());
                }
            } else {
                let mut members = group.member_node_ids.clone();
                members
                    .sort_by_key(|node_id| row_order.get(node_id).copied().unwrap_or(usize::MAX));
                for member_id in members {
                    if emitted.contains(&member_id) {
                        continue;
                    }
                    if let Some(member_row) = row_by_node.get(&member_id) {
                        entries.push(LineageTableEntry {
                            row: member_row.clone(),
                            indent_level: 1,
                            group_id: Some(group.group_id.clone()),
                            group_label: Some(group.label.clone()),
                            is_group_representative: false,
                            hidden_group_member_count: 0,
                        });
                        emitted.insert(member_id);
                    }
                }
            }
        }

        for row in rows {
            if emitted.contains(&row.node_id) {
                continue;
            }
            entries.push(LineageTableEntry {
                row: row.clone(),
                indent_level: 0,
                group_id: None,
                group_label: None,
                is_group_representative: false,
                hidden_group_member_count: 0,
            });
            emitted.insert(row.node_id.clone());
        }

        entries
    }

    fn lineage_group_color(group_id: &str) -> egui::Color32 {
        let palette = [
            egui::Color32::from_rgb(115, 146, 186),
            egui::Color32::from_rgb(173, 132, 84),
            egui::Color32::from_rgb(120, 156, 106),
            egui::Color32::from_rgb(163, 118, 143),
            egui::Color32::from_rgb(119, 155, 160),
            egui::Color32::from_rgb(174, 148, 88),
        ];
        let mut hasher = DefaultHasher::new();
        group_id.hash(&mut hasher);
        let index = (hasher.finish() as usize) % palette.len();
        palette[index]
    }

    fn start_lineage_group_draft_from_marked(
        &mut self,
        representative_node_id: &str,
        valid_node_ids: &HashSet<String>,
    ) -> Result<usize, String> {
        let representative = representative_node_id.trim();
        if representative.is_empty() {
            return Err("Representative node id is required".to_string());
        }
        if !valid_node_ids.contains(representative) {
            return Err(format!(
                "Representative node '{representative}' is not present in current lineage"
            ));
        }
        let mut members: Vec<String> = self
            .lineage_group_marked_nodes
            .iter()
            .filter(|node_id| valid_node_ids.contains(*node_id))
            .cloned()
            .collect();
        members.retain(|node_id| node_id != representative);
        members.sort();
        members.dedup();
        if members.is_empty() {
            return Err(
                "No marked members available. Mark one or more nodes first via context menu."
                    .to_string(),
            );
        }
        self.lineage_group_form_editing_id = None;
        self.lineage_group_form_representative = representative.to_string();
        self.lineage_group_form_members = members.join(", ");
        if self.lineage_group_form_label.trim().is_empty() {
            self.lineage_group_form_label = format!("Group {}", representative);
        }
        Ok(members.len())
    }

    fn apply_lineage_group_form(
        &mut self,
        valid_node_ids: &HashSet<String>,
    ) -> Result<String, String> {
        let representative = self.lineage_group_form_representative.trim().to_string();
        if representative.is_empty() {
            return Err("Representative node id is required".to_string());
        }
        if !valid_node_ids.contains(&representative) {
            return Err(format!(
                "Representative node '{representative}' is not present in current lineage"
            ));
        }

        let mut label = self.lineage_group_form_label.trim().to_string();
        if label.is_empty() {
            label = format!("Group {}", representative);
        }

        let mut members = Self::parse_lineage_group_node_list(&self.lineage_group_form_members);
        members.retain(|node_id| node_id != &representative);
        if members.is_empty() {
            return Err("At least one member node id is required".to_string());
        }

        let mut draft_groups = self.lineage_node_groups.clone();
        if let Some(group_id) = self.lineage_group_form_editing_id.clone() {
            if let Some(group) = draft_groups
                .iter_mut()
                .find(|group| group.group_id == group_id)
            {
                group.label = label.clone();
                group.representative_node_id = representative.clone();
                group.member_node_ids = members.clone();
            } else {
                draft_groups.push(PersistedLineageNodeGroup {
                    group_id,
                    label: label.clone(),
                    representative_node_id: representative.clone(),
                    member_node_ids: members.clone(),
                    collapsed: false,
                });
            }
        } else {
            draft_groups.push(PersistedLineageNodeGroup {
                group_id: Self::next_lineage_group_id(&draft_groups),
                label: label.clone(),
                representative_node_id: representative.clone(),
                member_node_ids: members.clone(),
                collapsed: false,
            });
        }

        let sanitized_draft = Self::sanitize_lineage_node_groups(&draft_groups, valid_node_ids);
        let adjusted = sanitized_draft != draft_groups;
        self.lineage_node_groups = sanitized_draft;
        self.lineage_group_form_editing_id = None;
        self.lineage_group_form_label.clear();
        self.lineage_group_form_members.clear();
        self.lineage_group_form_representative.clear();
        self.lineage_group_marked_nodes.clear();

        if adjusted {
            Ok("Group saved with normalization (overlaps/invalid nodes were removed)".to_string())
        } else {
            Ok(format!("Saved group '{label}'"))
        }
    }

    fn render_lineage_group_context_menu(
        &mut self,
        ui: &mut Ui,
        node_id: &str,
        context_row: Option<&LineageRow>,
        leaf_node_ids: &HashSet<String>,
        valid_node_ids: &HashSet<String>,
        persist_workspace_after_frame: &mut bool,
    ) {
        ui.label(format!("Node: {node_id}"));
        ui.separator();
        if let Some(row) = context_row {
            for kind in [
                LineageCopyPayloadKind::NodeId,
                LineageCopyPayloadKind::PrimaryId,
                LineageCopyPayloadKind::DisplayLabel,
                LineageCopyPayloadKind::RowSummary,
            ] {
                if ui
                    .button(kind.menu_label())
                    .on_hover_text(kind.hover_text())
                    .clicked()
                {
                    self.copy_lineage_row_payload_to_clipboard(ui.ctx(), row, kind);
                    ui.close();
                }
            }
        } else if ui
            .button(LineageCopyPayloadKind::NodeId.menu_label())
            .on_hover_text(LineageCopyPayloadKind::NodeId.hover_text())
            .clicked()
        {
            self.copy_lineage_node_id_to_clipboard(ui.ctx(), node_id);
            ui.close();
        }

        ui.separator();
        let marked = self.lineage_group_marked_nodes.contains(node_id);
        if !marked {
            if ui
                .button("Mark for node-group")
                .on_hover_text("Mark this node as a candidate member for new node groups")
                .clicked()
            {
                self.lineage_group_marked_nodes.insert(node_id.to_string());
                self.lineage_group_status = format!("Marked node '{node_id}' for grouping");
                ui.close();
            }
        } else if ui
            .button("Unmark for node-group")
            .on_hover_text("Remove this node from marked candidate members")
            .clicked()
        {
            self.lineage_group_marked_nodes.remove(node_id);
            self.lineage_group_status = format!("Unmarked node '{node_id}'");
            ui.close();
        }

        if ui
            .button("Use as draft representative")
            .on_hover_text("Fill node-group form from currently marked nodes")
            .clicked()
        {
            match self.start_lineage_group_draft_from_marked(node_id, valid_node_ids) {
                Ok(count) => {
                    self.lineage_group_status = format!(
                        "Draft prepared with representative '{node_id}' and {count} member(s)"
                    );
                }
                Err(err) => {
                    self.lineage_group_status = err;
                }
            }
            ui.close();
        }

        if ui
            .button("Create group now (representative)")
            .on_hover_text("Create a new group immediately from marked nodes")
            .clicked()
        {
            let status = match self.start_lineage_group_draft_from_marked(node_id, valid_node_ids) {
                Ok(_) => match self.apply_lineage_group_form(valid_node_ids) {
                    Ok(status) => {
                        *persist_workspace_after_frame = true;
                        status
                    }
                    Err(err) => err,
                },
                Err(err) => err,
            };
            self.lineage_group_status = status;
            ui.close();
        }

        ui.separator();
        let is_sequence_node = context_row
            .map(|row| matches!(row.kind, LineageNodeKind::Sequence))
            .unwrap_or(false);
        let is_leaf = leaf_node_ids.contains(node_id);
        if is_sequence_node {
            if ui
                .add_enabled(is_leaf, egui::Button::new("Rename (leaf only)"))
                .on_hover_text(if is_leaf {
                    "Rename this leaf node (display name only)"
                } else {
                    "Rename is allowed only for leaf nodes in this first implementation"
                })
                .clicked()
            {
                self.lineage_node_rename_target_id = Some(node_id.to_string());
                self.lineage_node_rename_text = context_row
                    .map(|row| row.display_name.clone())
                    .unwrap_or_default();
                ui.close();
            }
            if ui
                .add_enabled(is_leaf, egui::Button::new("Remove (leaf only)"))
                .on_hover_text(if is_leaf {
                    "Remove this leaf node and its sequence from project state"
                } else {
                    "Remove is allowed only for leaf nodes in this first implementation"
                })
                .clicked()
            {
                self.lineage_group_status = match self.request_remove_leaf_lineage_node(node_id) {
                    Ok(status) => status,
                    Err(err) => err,
                };
                ui.close();
            }
        } else {
            ui.label("Rename/remove is currently available for sequence leaf nodes only.");
        }
    }

    fn lineage_leaf_node_ids(
        rows: &[LineageRow],
        lineage_edges: &[(String, String, String)],
    ) -> HashSet<String> {
        let mut node_ids: HashSet<String> = rows.iter().map(|row| row.node_id.clone()).collect();
        let mut non_leaf_ids: HashSet<String> = HashSet::new();
        for (from_node, to_node, _) in lineage_edges {
            if node_ids.contains(from_node) && node_ids.contains(to_node) {
                non_leaf_ids.insert(from_node.clone());
            }
        }
        node_ids.retain(|node_id| !non_leaf_ids.contains(node_id));
        node_ids
    }

    fn request_remove_leaf_lineage_node(&mut self, node_id: &str) -> Result<String, String> {
        let node_id = node_id.trim();
        if node_id.is_empty() {
            return Err("Node id cannot be empty".to_string());
        }
        let is_leaf = {
            let engine = self.engine.read().unwrap();
            let state = engine.state();
            if !state.lineage.nodes.contains_key(node_id) {
                return Err(format!("Node '{}' not found", node_id));
            }
            !state
                .lineage
                .edges
                .iter()
                .any(|edge| edge.from_node_id == node_id)
        };
        if !is_leaf {
            return Err(format!(
                "Node '{}' is not a leaf; remove downstream nodes first",
                node_id
            ));
        }
        self.lineage_node_remove_target_id = Some(node_id.to_string());
        Ok(format!(
            "Removal requested for leaf node '{}'; confirm in dialog.",
            node_id
        ))
    }

    fn remove_leaf_lineage_node(&mut self, node_id: &str) -> Result<String, String> {
        let node_id = node_id.trim();
        if node_id.is_empty() {
            return Err("Node id cannot be empty".to_string());
        }
        let seq_id = {
            let mut engine = self.engine.write().unwrap();
            let state = engine.state_mut();
            let Some(node) = state.lineage.nodes.get(node_id).cloned() else {
                return Err(format!("Node '{}' not found", node_id));
            };
            if state
                .lineage
                .edges
                .iter()
                .any(|edge| edge.from_node_id == node_id)
            {
                return Err(format!(
                    "Node '{}' is not a leaf; remove downstream nodes first",
                    node_id
                ));
            }
            let seq_id = node.seq_id.clone();
            if state.sequences.remove(&seq_id).is_none() {
                return Err(format!(
                    "Sequence '{}' for node '{}' is missing",
                    seq_id, node_id
                ));
            }
            state.lineage.seq_to_node.remove(&seq_id);
            state.lineage.nodes.remove(node_id);
            state
                .lineage
                .edges
                .retain(|edge| edge.from_node_id != node_id && edge.to_node_id != node_id);

            let mut empty_container_ids: HashSet<String> = HashSet::new();
            for (container_id, container) in state.container_state.containers.iter_mut() {
                container.members.retain(|member| member != &seq_id);
                if container.members.is_empty() {
                    empty_container_ids.insert(container_id.clone());
                }
            }
            for container_id in &empty_container_ids {
                state.container_state.containers.remove(container_id);
            }
            for arrangement in state.container_state.arrangements.values_mut() {
                arrangement
                    .lane_container_ids
                    .retain(|container_id| !empty_container_ids.contains(container_id));
            }
            state
                .container_state
                .arrangements
                .retain(|_, arrangement| !arrangement.lane_container_ids.is_empty());
            let live_sequence_ids: HashSet<String> =
                state.sequences.keys().cloned().collect::<HashSet<_>>();
            let container_members_by_id: HashMap<String, HashSet<String>> = state
                .container_state
                .containers
                .iter()
                .map(|(container_id, container)| {
                    (
                        container_id.clone(),
                        container.members.iter().cloned().collect::<HashSet<_>>(),
                    )
                })
                .collect();
            state
                .container_state
                .seq_to_latest_container
                .retain(|member_seq, container_id| {
                    member_seq != &seq_id
                        && live_sequence_ids.contains(member_seq)
                        && container_members_by_id
                            .get(container_id)
                            .map(|members| members.contains(member_seq))
                            .unwrap_or(false)
                });
            seq_id
        };

        if let Some(viewport_id) = self.find_open_sequence_viewport_id(&seq_id) {
            if let Ok(mut to_close) = self.windows_to_close.write()
                && !to_close.contains(&viewport_id)
            {
                to_close.push(viewport_id);
            }
            self.pending_focus_viewports.retain(|id| *id != viewport_id);
        }
        self.new_windows
            .retain(|window| window.sequence_id().as_deref() != Some(seq_id.as_str()));
        self.lineage_group_marked_nodes.remove(node_id);
        if self
            .lineage_graph_selected_node_id
            .as_deref()
            .is_some_and(|selected| selected == node_id)
        {
            self.lineage_graph_selected_node_id = None;
        }
        if self
            .lineage_node_rename_target_id
            .as_deref()
            .is_some_and(|target| target == node_id)
        {
            self.lineage_node_rename_target_id = None;
            self.lineage_node_rename_text.clear();
        }
        if self
            .lineage_node_remove_target_id
            .as_deref()
            .is_some_and(|target| target == node_id)
        {
            self.lineage_node_remove_target_id = None;
        }
        self.lineage_cache_valid = false;
        Ok(format!(
            "Removed leaf node '{}' (sequence '{}')",
            node_id, seq_id
        ))
    }

    fn rename_leaf_lineage_node(
        &mut self,
        node_id: &str,
        new_name: &str,
    ) -> Result<String, String> {
        let node_id = node_id.trim();
        let new_name = new_name.trim();
        if node_id.is_empty() {
            return Err("Node id cannot be empty".to_string());
        }
        if new_name.is_empty() {
            return Err("New name cannot be empty".to_string());
        }
        let seq_id = {
            let mut engine = self.engine.write().unwrap();
            let state = engine.state_mut();
            let Some(node) = state.lineage.nodes.get(node_id).cloned() else {
                return Err(format!("Node '{}' not found", node_id));
            };
            if state
                .lineage
                .edges
                .iter()
                .any(|edge| edge.from_node_id == node_id)
            {
                return Err(format!(
                    "Node '{}' is not a leaf; rename is leaf-only in this first implementation",
                    node_id
                ));
            }
            let seq_id = node.seq_id.clone();
            let Some(dna) = state.sequences.get_mut(&seq_id) else {
                return Err(format!(
                    "Sequence '{}' for node '{}' is missing",
                    seq_id, node_id
                ));
            };
            dna.set_name(new_name.to_string());
            seq_id
        };
        self.lineage_cache_valid = false;
        Ok(format!(
            "Renamed leaf node '{}' (sequence '{}') to '{}'",
            node_id, seq_id, new_name
        ))
    }

    fn render_main_workspace_host(&mut self, ui: &mut Ui, project_dirty: bool) {
        window_backdrop::paint_window_backdrop(
            ui,
            WindowBackdropKind::Main,
            &self.window_backdrops,
        );
        with_window_content_inset(ui, |ui| {
            let workspace_size = ui.available_size_before_wrap();
            let width = workspace_size.x.max(720.0);
            let height = workspace_size.y.max(420.0);
            ui.allocate_ui_with_layout(
                egui::vec2(width, height),
                egui::Layout::top_down(egui::Align::Min),
                |ui| {
                    let status_height = ui.spacing().interact_size.y + 16.0;
                    let lineage_height = (ui.available_height() - status_height).max(240.0);
                    ui.allocate_ui_with_layout(
                        egui::vec2(width, lineage_height),
                        egui::Layout::top_down(egui::Align::Min),
                        |ui| {
                            self.render_main_lineage(ui);
                        },
                    );
                    ui.separator();
                    if project_dirty {
                        ui.label("Status: unsaved changes");
                    } else {
                        ui.label("Status: saved");
                    }
                },
            );
        });
    }

    fn render_hosted_main_workspace_window(&mut self, ctx: &egui::Context, project_dirty: bool) {
        let min_size = Vec2::new(900.0, 560.0);
        let mut open = true;
        let spec = crate::egui_compat::HostedWindowSpec::new(
            format!("Project — {}", self.current_project_name()),
            Self::main_workspace_hosted_window_id(),
            Vec2::new(1280.0, 860.0),
            min_size,
        )
        .collapsible(false)
        .resizable(true)
        .cleanup_legacy_title_layer(true);
        crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
            self.render_main_workspace_host(ui, project_dirty);
        });
    }

    fn render_root_workspace(&mut self, ctx: &egui::Context, project_dirty: bool) {
        crate::egui_compat::show_central_panel(
            ctx,
            egui::CentralPanel::default().frame(egui::Frame::NONE),
            |ui| {
                let host_rect = ui.max_rect();
                ui.painter()
                    .rect_filled(host_rect, 0.0, egui::Color32::from_gray(218));
                if !Self::should_embed_child_viewports() {
                    self.render_main_workspace_host(ui, project_dirty);
                }
            },
        );
        if Self::should_embed_child_viewports() {
            self.render_hosted_main_workspace_window(ctx, project_dirty);
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
        self.configuration_graphics.show_linear_sequence_panel =
            defaults.show_linear_sequence_panel;
        self.configuration_graphics
            .sequence_panel_max_text_length_bp = defaults.sequence_panel_max_text_length_bp;
        self.configuration_graphics
            .auto_hide_sequence_panel_when_linear_bases_visible =
            defaults.auto_hide_sequence_panel_when_linear_bases_visible;
        self.configuration_graphics.show_map_panel = defaults.show_map_panel;
        self.configuration_graphics.show_features = defaults.show_features;
        self.configuration_graphics.show_cds_features = defaults.show_cds_features;
        self.configuration_graphics.show_gene_features = defaults.show_gene_features;
        self.configuration_graphics.show_mrna_features = defaults.show_mrna_features;
        self.configuration_graphics.show_repeat_features = defaults.show_repeat_features;
        self.configuration_graphics.show_array_features = defaults.show_array_features;
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
        self.configuration_graphics.restriction_enzyme_display_mode =
            defaults.restriction_enzyme_display_mode;
        self.configuration_graphics.preferred_restriction_enzymes =
            defaults.preferred_restriction_enzymes.clone();
        self.configuration_graphics.show_gc_contents = defaults.show_gc_contents;
        self.configuration_graphics.gc_content_bin_size_bp = defaults.gc_content_bin_size_bp;
        self.configuration_graphics.show_open_reading_frames = defaults.show_open_reading_frames;
        self.configuration_graphics.show_methylation_sites = defaults.show_methylation_sites;
        self.configuration_graphics.feature_details_font_size = defaults.feature_details_font_size;
        self.configuration_graphics
            .linear_external_feature_label_font_size =
            defaults.linear_external_feature_label_font_size;
        self.configuration_graphics
            .linear_external_feature_label_background_opacity =
            defaults.linear_external_feature_label_background_opacity;
        self.configuration_graphics.linear_view_start_bp = defaults.linear_view_start_bp;
        self.configuration_graphics.linear_view_span_bp = defaults.linear_view_span_bp;
        self.configuration_graphics.linear_view_vertical_offset_px =
            defaults.linear_view_vertical_offset_px;
        self.configuration_graphics
            .linear_sequence_base_text_max_view_span_bp =
            defaults.linear_sequence_base_text_max_view_span_bp;
        self.configuration_graphics
            .linear_sequence_helical_letters_enabled =
            defaults.linear_sequence_helical_letters_enabled;
        self.configuration_graphics
            .linear_sequence_helical_max_view_span_bp =
            defaults.linear_sequence_helical_max_view_span_bp;
        self.configuration_graphics
            .linear_sequence_condensed_max_view_span_bp =
            defaults.linear_sequence_condensed_max_view_span_bp;
        self.configuration_graphics
            .linear_sequence_letter_layout_mode = defaults.linear_sequence_letter_layout_mode;
        self.configuration_graphics
            .linear_sequence_helical_phase_offset_bp =
            defaults.linear_sequence_helical_phase_offset_bp;
        self.configuration_graphics.linear_show_double_strand_bases =
            defaults.linear_show_double_strand_bases;
        self.configuration_graphics.linear_helical_parallel_strands =
            defaults.linear_helical_parallel_strands;
        self.configuration_graphics
            .linear_hide_backbone_when_sequence_bases_visible =
            defaults.linear_hide_backbone_when_sequence_bases_visible;
        self.configuration_graphics
            .linear_reverse_strand_use_upside_down_letters =
            defaults.linear_reverse_strand_use_upside_down_letters;
        self.configuration_graphics.reverse_strand_visual_opacity =
            defaults.reverse_strand_visual_opacity;
        self.configuration_graphics_dirty = true;
    }

    fn apply_configuration_window_backdrops(&mut self) {
        self.window_backdrops = self.configuration_window_backdrops.clone();
        window_backdrop::set_window_backdrop_settings(self.window_backdrops.clone());
        self.window_backdrop_path_status_cache.clear();
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

    fn apply_configuration_language(&mut self) {
        self.ui_language = self.configuration_ui_language;
        self.i18n.set_language(self.ui_language);
        self.configuration_language_dirty = false;
        self.configuration_status =
            format!("Interface language applied ({})", self.ui_language.label());
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

    fn window_backdrop_row_mut(
        settings: &mut WindowBackdropSettings,
        row_index: usize,
    ) -> (WindowBackdropKind, &'static str, &mut [u8; 3], &mut String) {
        match row_index {
            0 => (
                WindowBackdropKind::Main,
                "Main",
                &mut settings.main_tint_rgb,
                &mut settings.main_image_path,
            ),
            1 => (
                WindowBackdropKind::Sequence,
                "Sequence",
                &mut settings.sequence_tint_rgb,
                &mut settings.sequence_image_path,
            ),
            2 => (
                WindowBackdropKind::Splicing,
                "Splicing",
                &mut settings.splicing_tint_rgb,
                &mut settings.splicing_image_path,
            ),
            3 => (
                WindowBackdropKind::Pool,
                "Pool",
                &mut settings.pool_tint_rgb,
                &mut settings.pool_image_path,
            ),
            4 => (
                WindowBackdropKind::Configuration,
                "Configuration",
                &mut settings.configuration_tint_rgb,
                &mut settings.configuration_image_path,
            ),
            5 => (
                WindowBackdropKind::Help,
                "Help",
                &mut settings.help_tint_rgb,
                &mut settings.help_image_path,
            ),
            6 => (
                WindowBackdropKind::AgentAssistant,
                "Agent Assistant",
                &mut settings.agent_assistant_tint_rgb,
                &mut settings.agent_assistant_image_path,
            ),
            _ => (
                WindowBackdropKind::Main,
                "Unknown",
                &mut settings.main_tint_rgb,
                &mut settings.main_image_path,
            ),
        }
    }

    fn window_backdrop_path_status(
        path: &str,
        status_cache: &mut HashMap<String, (egui::Color32, String)>,
    ) -> (egui::Color32, String) {
        let trimmed = path.trim();
        let cache_key = trimmed.to_string();
        if let Some(cached) = status_cache.get(&cache_key) {
            return cached.clone();
        }

        let status = if trimmed.is_empty() {
            (
                egui::Color32::from_gray(120),
                "No image path configured (text watermark/tint fallback only)".to_string(),
            )
        } else {
            match window_backdrop::validate_window_backdrop_image_path(trimmed) {
                Ok(resolved) => (
                    egui::Color32::from_rgb(20, 140, 45),
                    format!("Resolved image: {resolved}"),
                ),
                Err(message) => (
                    egui::Color32::from_rgb(180, 50, 50),
                    format!("Path check failed: {message}"),
                ),
            }
        };
        status_cache.insert(cache_key, status.clone());
        status
    }

    fn render_window_backdrop_path_table(
        ui: &mut Ui,
        settings: &mut WindowBackdropSettings,
        changed: &mut bool,
        path_status_cache: &mut HashMap<String, (egui::Color32, String)>,
    ) {
        const ROW_COUNT: usize = 7;
        ui.label(crate::i18n::tr("configuration.window_styling.table_title"));
        egui::Grid::new("window_backdrop_path_grid")
            .num_columns(4)
            .striped(true)
            .spacing(egui::vec2(10.0, 6.0))
            .show(ui, |ui| {
                ui.strong(crate::i18n::tr("configuration.window_styling.window"));
                ui.strong(crate::i18n::tr("configuration.window_styling.tint"));
                ui.strong(crate::i18n::tr("configuration.window_styling.path"));
                ui.strong(crate::i18n::tr("configuration.window_styling.actions"));
                ui.end_row();

                for row_index in 0..ROW_COUNT {
                    let (kind, label, tint_rgb, value) =
                        Self::window_backdrop_row_mut(settings, row_index);
                    ui.label(label);

                    ui.horizontal(|ui| {
                        let mut color = egui::Color32::from_rgb(tint_rgb[0], tint_rgb[1], tint_rgb[2]);
                        let response = ui
                            .color_edit_button_srgba(&mut color)
                            .on_hover_text("Pick tint color for this window type");
                        if response.changed() {
                            *tint_rgb = [color.r(), color.g(), color.b()];
                            *changed = true;
                        }
                        ui.monospace(format!("#{:02X}{:02X}{:02X}", color.r(), color.g(), color.b()));
                    });

                    ui.horizontal(|ui| {
                        let field_width =
                            Self::clamp_configuration_backdrop_path_field_width(ui.available_width());
                        if ui
                            .add_sized(
                                [field_width, 0.0],
                                egui::TextEdit::singleline(value)
                                    .hint_text("/absolute/path/to/image.png"),
                            )
                            .on_hover_text(
                                "Optional absolute or working-directory-relative image path for this window type",
                            )
                            .changed()
                        {
                            *changed = true;
                            path_status_cache.clear();
                        }
                        let (status_color, status_message) =
                            Self::window_backdrop_path_status(value, path_status_cache);
                        ui.label(egui::RichText::new("●").color(status_color))
                            .on_hover_text(status_message);
                    });

                    ui.horizontal(|ui| {
                        let default_tint = window_backdrop::default_tint_rgb_for_kind(kind);
                        if ui
                            .small_button(crate::i18n::tr("button.browse"))
                            .on_hover_text("Pick an image file for this window backdrop")
                            .clicked()
                            && let Some(path) = rfd::FileDialog::new()
                                .add_filter("Images", &["png", "jpg", "jpeg", "gif", "bmp", "webp"])
                                .pick_file()
                            {
                                *value = path.display().to_string();
                                *changed = true;
                                path_status_cache.clear();
                            }
                        if ui
                            .add_enabled(
                                !value.trim().is_empty(),
                                egui::Button::new(crate::i18n::tr("button.clear")),
                            )
                            .on_hover_text("Clear custom image path for this window type")
                            .clicked()
                        {
                            value.clear();
                            *changed = true;
                            path_status_cache.clear();
                        }
                        if ui
                            .add_enabled(
                                *tint_rgb != default_tint,
                                egui::Button::new(crate::i18n::tr("button.reset_color")),
                            )
                            .on_hover_text("Reset tint color for this window type")
                            .clicked()
                        {
                            *tint_rgb = default_tint;
                            *changed = true;
                        }
                    });

                    ui.end_row();
                }
            });
        ui.small(crate::i18n::tr(
            "configuration.window_styling.saved_after_apply",
        ));
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

    fn sync_open_windows_if_display_changed(&mut self, ctx: &egui::Context) {
        let stamp = self.current_display_change_stamp();
        if self.last_display_sync_stamp == stamp {
            return;
        }
        self.last_display_sync_stamp = stamp;
        self.refresh_open_sequence_windows(ctx);
    }

    fn render_configuration_external_tab(&mut self, ui: &mut Ui) {
        ui.label(self.tr("configuration.external.description"));
        ui.separator();
        ui.label(self.tr("configuration.external.rnapkin_override"));
        let rnapkin_edit_response = ui.add(
            egui::TextEdit::singleline(&mut self.configuration_rnapkin_executable)
                .hint_text("Leave empty to use PATH lookup for 'rnapkin'"),
        );
        if rnapkin_edit_response.changed() {
            self.clear_rnapkin_validation();
        }
        let active_rnapkin =
            tool_overrides::active_resolution_label("GENTLE_RNAPKIN_BIN", "rnapkin");
        ui.monospace(format!(
            "{}: {active_rnapkin}",
            self.tr("configuration.external.active_resolution")
        ));

        ui.separator();
        ui.label(self.tr("configuration.external.makeblastdb_override"));
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
        ui.monospace(format!(
            "{}: {active_makeblastdb}",
            self.tr("configuration.external.active_makeblastdb")
        ));

        ui.label(self.tr("configuration.external.blastn_override"));
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
        ui.monospace(format!(
            "{}: {active_blastn}",
            self.tr("configuration.external.active_blastn")
        ));

        ui.separator();
        ui.label(self.tr("configuration.external.bigwig_override"));
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
        ui.monospace(format!(
            "{}: {active_bigwig}",
            self.tr("configuration.external.active_bigwig")
        ));

        ui.horizontal_wrapped(|ui| {
            if ui
                .button(self.tr("configuration.external.use_path"))
                .on_hover_text("Clear rnapkin override and use PATH lookup")
                .clicked()
            {
                self.configuration_rnapkin_executable.clear();
                self.clear_rnapkin_validation();
            }
            if ui
                .button(self.tr("configuration.external.use_path_blast"))
                .on_hover_text("Clear makeblastdb/blastn overrides and use PATH lookup")
                .clicked()
            {
                self.configuration_makeblastdb_executable.clear();
                self.configuration_blastn_executable.clear();
                self.clear_blast_validation();
            }
            if ui
                .button(self.tr("configuration.external.use_path_bigwig"))
                .on_hover_text("Clear bigWigToBedGraph override and use PATH lookup")
                .clicked()
            {
                self.configuration_bigwig_to_bedgraph_executable.clear();
            }
            if ui
                .button(self.tr("configuration.external.validate_rnapkin"))
                .on_hover_text("Run rnapkin --version and capture validation status")
                .clicked()
            {
                self.validate_rnapkin_executable();
            }
            if ui
                .button(self.tr("configuration.external.validate_blast"))
                .on_hover_text("Run makeblastdb/blastn --version checks")
                .clicked()
            {
                self.validate_blast_executables();
            }
            if ui
                .button(self.tr("configuration.external.apply"))
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
        let normalized_detail_font = Self::clamp_feature_details_font_size(
            self.configuration_graphics.feature_details_font_size,
        );
        if (self.configuration_graphics.feature_details_font_size - normalized_detail_font).abs()
            > f32::EPSILON
        {
            self.configuration_graphics.feature_details_font_size = normalized_detail_font;
            self.configuration_graphics_dirty = true;
        }
        let normalized_external_label_font = Self::clamp_external_feature_label_font_size(
            self.configuration_graphics
                .linear_external_feature_label_font_size,
        );
        if (self
            .configuration_graphics
            .linear_external_feature_label_font_size
            - normalized_external_label_font)
            .abs()
            > f32::EPSILON
        {
            self.configuration_graphics
                .linear_external_feature_label_font_size = normalized_external_label_font;
            self.configuration_graphics_dirty = true;
        }
        let normalized_external_label_bg = Self::clamp_external_feature_label_background_opacity(
            self.configuration_graphics
                .linear_external_feature_label_background_opacity,
        );
        if (self
            .configuration_graphics
            .linear_external_feature_label_background_opacity
            - normalized_external_label_bg)
            .abs()
            > f32::EPSILON
        {
            self.configuration_graphics
                .linear_external_feature_label_background_opacity = normalized_external_label_bg;
            self.configuration_graphics_dirty = true;
        }
        let normalized_reverse_opacity = Self::clamp_reverse_strand_visual_opacity(
            self.configuration_graphics.reverse_strand_visual_opacity,
        );
        if (self.configuration_graphics.reverse_strand_visual_opacity - normalized_reverse_opacity)
            .abs()
            > f32::EPSILON
        {
            self.configuration_graphics.reverse_strand_visual_opacity = normalized_reverse_opacity;
            self.configuration_graphics_dirty = true;
        }
        let normalized_helical_phase_offset = Self::clamp_linear_helical_phase_offset_bp(
            self.configuration_graphics
                .linear_sequence_helical_phase_offset_bp,
        );
        if self
            .configuration_graphics
            .linear_sequence_helical_phase_offset_bp
            != normalized_helical_phase_offset
        {
            self.configuration_graphics
                .linear_sequence_helical_phase_offset_bp = normalized_helical_phase_offset;
            self.configuration_graphics_dirty = true;
        }
        let normalized_sequence_panel_text_limit = Self::clamp_sequence_panel_max_text_length_bp(
            self.configuration_graphics
                .sequence_panel_max_text_length_bp,
        );
        if self
            .configuration_graphics
            .sequence_panel_max_text_length_bp
            != normalized_sequence_panel_text_limit
        {
            self.configuration_graphics
                .sequence_panel_max_text_length_bp = normalized_sequence_panel_text_limit;
            self.configuration_graphics_dirty = true;
        }

        ui.label(self.tr("configuration.graphics.description"));
        ui.small(format!(
            "{}: {:.2} px",
            self.tr("configuration.graphics.current_feature_font"),
            self.configuration_graphics.feature_details_font_size
        ));
        ui.separator();
        let mut changed = false;
        let mut live_font_changed = false;
        let mut backdrop_changed = false;

        ui.heading(self.tr("configuration.graphics.panels"));
        changed |= ui
            .checkbox(
                &mut self.configuration_graphics.show_sequence_panel,
                crate::i18n::tr("configuration.graphics.show_sequence_panel_circular"),
            )
            .changed();
        changed |= ui
            .checkbox(
                &mut self.configuration_graphics.show_linear_sequence_panel,
                crate::i18n::tr("configuration.graphics.show_sequence_panel_linear"),
            )
            .changed();
        ui.horizontal(|ui| {
            ui.label(crate::i18n::tr(
                "configuration.graphics.sequence_panel_max_text_length",
            ));
            if ui
                .add(
                    egui::DragValue::new(
                        &mut self
                            .configuration_graphics
                            .sequence_panel_max_text_length_bp,
                    )
                    .range(0..=5_000_000)
                    .speed(1000.0)
                    .suffix(" bp"),
                )
                .on_hover_text(
                    "Maximum sequence length shown in the text panel. 0 means unlimited.",
                )
                .changed()
            {
                changed = true;
            }
        });
        changed |= ui
            .checkbox(
                &mut self.configuration_graphics.show_map_panel,
                crate::i18n::tr("configuration.graphics.show_map_panel"),
            )
            .changed();
        changed |= ui
            .checkbox(
                &mut self
                    .configuration_graphics
                    .auto_hide_sequence_panel_when_linear_bases_visible,
                crate::i18n::tr("configuration.graphics.auto_hide_sequence_panel"),
            )
            .changed();

        ui.separator();
        ui.heading(self.tr("configuration.graphics.linear_dna_base_rendering"));
        changed |= ui
            .checkbox(
                &mut self.configuration_graphics.linear_show_double_strand_bases,
                crate::i18n::tr("configuration.graphics.show_reverse_letters"),
            )
            .changed();
        if self.configuration_graphics.linear_show_double_strand_bases {
            changed |= ui
                .checkbox(
                    &mut self.configuration_graphics.linear_helical_parallel_strands,
                    crate::i18n::tr("sequence.keep_helical_parallel"),
                )
                .on_hover_text(
                    "When enabled, forward and reverse strands move in parallel during helical rendering; disable for mirrored/cross-over slant.",
                )
                .changed();
        }
        let hide_backbone_changed = ui
            .checkbox(
                &mut self
                    .configuration_graphics
                    .linear_hide_backbone_when_sequence_bases_visible,
                crate::i18n::tr("sequence.hide_backbone_when_letters_visible"),
            )
            .changed();
        changed |= hide_backbone_changed;
        let helical_letters_enabled_changed = ui
            .checkbox(
                &mut self
                    .configuration_graphics
                    .linear_sequence_helical_letters_enabled,
                crate::i18n::tr("sequence.enable_compressed_letters"),
            )
            .on_hover_text(
                "In Auto mode: disabled means dense views switch letters OFF instead of helical/condensed. Forced modes ignore this toggle.",
            )
            .changed();
        changed |= helical_letters_enabled_changed;
        ui.horizontal(|ui| {
            ui.label(crate::i18n::tr("sequence.dna_letter_routing_mode"));
            let mut selected = self
                .configuration_graphics
                .linear_sequence_letter_layout_mode;
            let mut layout_changed = false;
            egui::ComboBox::from_id_salt("config_linear_helical_letter_layout_mode")
                .selected_text(match selected {
                    LinearSequenceLetterLayoutMode::AutoAdaptive => "Auto adaptive",
                    LinearSequenceLetterLayoutMode::StandardLinear => "Force standard",
                    LinearSequenceLetterLayoutMode::ContinuousHelical => "Force helical",
                    LinearSequenceLetterLayoutMode::Condensed10Row => "Force condensed 10-row",
                })
                .show_ui(ui, |ui| {
                    layout_changed |= ui
                        .selectable_value(
                            &mut selected,
                            LinearSequenceLetterLayoutMode::AutoAdaptive,
                            "Auto adaptive",
                        )
                        .changed();
                    layout_changed |= ui
                        .selectable_value(
                            &mut selected,
                            LinearSequenceLetterLayoutMode::StandardLinear,
                            "Force standard",
                        )
                        .changed();
                    layout_changed |= ui
                        .selectable_value(
                            &mut selected,
                            LinearSequenceLetterLayoutMode::ContinuousHelical,
                            "Force helical",
                        )
                        .changed();
                    layout_changed |= ui
                        .selectable_value(
                            &mut selected,
                            LinearSequenceLetterLayoutMode::Condensed10Row,
                            "Force condensed 10-row",
                        )
                        .changed();
                })
                .response
                .on_hover_text(
                    "Auto mode picks STANDARD/HELICAL/CONDENSED from viewport density. Force modes pin one route regardless of density.",
                );
            if layout_changed {
                self.configuration_graphics
                    .linear_sequence_letter_layout_mode = selected;
                changed = true;
            }
        });
        changed |= ui
            .checkbox(
                &mut self
                    .configuration_graphics
                    .linear_reverse_strand_use_upside_down_letters,
                crate::i18n::tr("sequence.rotate_reverse_letters"),
            )
            .changed();
        ui.horizontal(|ui| {
            ui.label(crate::i18n::tr("sequence.reverse_letter_opacity"));
            if ui
                .add(
                    egui::Slider::new(
                        &mut self.configuration_graphics.reverse_strand_visual_opacity,
                        0.2..=1.0,
                    )
                    .step_by(0.01),
                )
                .on_hover_text(
                    "Shared emphasis control for reverse-strand letters in linear map and sequence panel",
                )
                .changed()
            {
                changed = true;
            }
        });
        ui.horizontal(|ui| {
            ui.label(crate::i18n::tr("sequence.helical_phase_offset"));
            if ui
                .add(
                    egui::DragValue::new(
                        &mut self
                            .configuration_graphics
                            .linear_sequence_helical_phase_offset_bp,
                    )
                    .range(0..=9)
                    .speed(1.0)
                    .suffix(" bp"),
                )
                .on_hover_text(
                    "Row formula is row=(bp+offset)%10. Increasing offset shifts the top-to-bottom seam without changing DNA base order.",
                )
                .changed()
            {
                changed = true;
            }
        });
        ui.small("Adaptive DNA-letter routing now uses viewport density; these settings persist after Apply.");
        ui.horizontal(|ui| {
            ui.label("Feature detail font size");
            if ui
                .add(
                    egui::Slider::new(
                        &mut self.configuration_graphics.feature_details_font_size,
                        8.0..=24.0,
                    )
                    .step_by(0.25)
                    .suffix(" px"),
                )
                .changed()
            {
                changed = true;
                live_font_changed = true;
            }
        });
        ui.horizontal(|ui| {
            ui.label("Linear external feature label font");
            if ui
                .add(
                    egui::Slider::new(
                        &mut self
                            .configuration_graphics
                            .linear_external_feature_label_font_size,
                        8.0..=24.0,
                    )
                    .step_by(0.25)
                    .suffix(" px"),
                )
                .on_hover_text(
                    "Font size for linear-map externalized feature labels (with connector lines)",
                )
                .changed()
            {
                changed = true;
                live_font_changed = true;
            }
        });
        ui.horizontal(|ui| {
            ui.label("Linear external feature label background");
            if ui
                .add(
                    egui::Slider::new(
                        &mut self
                            .configuration_graphics
                            .linear_external_feature_label_background_opacity,
                        0.0..=1.0,
                    )
                    .step_by(0.01),
                )
                .on_hover_text("Background opacity used behind linear external feature labels")
                .changed()
            {
                changed = true;
                live_font_changed = true;
            }
        });

        ui.separator();
        ui.heading(self.tr("configuration.graphics.feature_layers"));
        changed |= ui
            .checkbox(
                &mut self.configuration_graphics.show_features,
                crate::i18n::tr("configuration.graphics.show_features"),
            )
            .changed();
        changed |= ui
            .checkbox(
                &mut self.configuration_graphics.show_cds_features,
                crate::i18n::tr("configuration.graphics.show_cds"),
            )
            .changed();
        changed |= ui
            .checkbox(
                &mut self.configuration_graphics.show_gene_features,
                crate::i18n::tr("configuration.graphics.show_gene"),
            )
            .changed();
        changed |= ui
            .checkbox(
                &mut self.configuration_graphics.show_mrna_features,
                crate::i18n::tr("configuration.graphics.show_mrna"),
            )
            .changed();
        changed |= ui
            .checkbox(
                &mut self.configuration_graphics.show_repeat_features,
                crate::i18n::tr("configuration.graphics.show_repeats"),
            )
            .on_hover_text(
                "Show RepeatMasker/rmsk-derived repeat_region and mobile-element features",
            )
            .changed();
        changed |= ui
            .checkbox(
                &mut self.configuration_graphics.show_array_features,
                crate::i18n::tr("configuration.graphics.show_arrays"),
            )
            .on_hover_text("Show genome-projected microarray contrast features")
            .changed();
        changed |= ui
            .checkbox(
                &mut self.configuration_graphics.show_tfbs,
                crate::i18n::tr("configuration.graphics.show_tfbs"),
            )
            .changed();
        changed |= ui
            .checkbox(
                &mut self.configuration_graphics.regulatory_tracks_near_baseline,
                crate::i18n::tr("configuration.graphics.regulatory_near_baseline"),
            )
            .changed();
        ui.horizontal(|ui| {
            ui.label(crate::i18n::tr(
                "configuration.graphics.regulatory_max_view_span",
            ));
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
        ui.heading(self.tr("configuration.graphics.overlays"));
        changed |= ui
            .checkbox(
                &mut self.configuration_graphics.show_restriction_enzymes,
                crate::i18n::tr("configuration.graphics.show_restriction_enzymes"),
            )
            .changed();
        ui.horizontal(|ui| {
            ui.label(crate::i18n::tr(
                "configuration.graphics.restriction_display_mode",
            ));
            let mut selected = self.configuration_graphics.restriction_enzyme_display_mode;
            let mut mode_changed = false;
            egui::ComboBox::from_id_salt("config_restriction_display_mode")
                .selected_text(selected.label())
                .show_ui(ui, |ui| {
                    mode_changed |= ui
                        .selectable_value(
                            &mut selected,
                            RestrictionEnzymeDisplayMode::PreferredOnly,
                            RestrictionEnzymeDisplayMode::PreferredOnly.label(),
                        )
                        .changed();
                    mode_changed |= ui
                        .selectable_value(
                            &mut selected,
                            RestrictionEnzymeDisplayMode::PreferredAndUnique,
                            RestrictionEnzymeDisplayMode::PreferredAndUnique.label(),
                        )
                        .changed();
                    mode_changed |= ui
                        .selectable_value(
                            &mut selected,
                            RestrictionEnzymeDisplayMode::UniqueOnly,
                            RestrictionEnzymeDisplayMode::UniqueOnly.label(),
                        )
                        .changed();
                    mode_changed |= ui
                        .selectable_value(
                            &mut selected,
                            RestrictionEnzymeDisplayMode::AllInView,
                            RestrictionEnzymeDisplayMode::AllInView.label(),
                        )
                        .changed();
                })
                .response
                .on_hover_text(
                    "Choose whether the DNA window shows preferred cutters, unique cutters, or every cut site in view.",
                );
            if mode_changed {
                self.configuration_graphics.restriction_enzyme_display_mode = selected;
                changed = true;
            }
        });
        ui.label(self.tr("configuration.graphics.preferred_restriction_enzymes"));
        let mut preferred_csv = self
            .configuration_graphics
            .preferred_restriction_enzymes
            .join(", ");
        if ui
            .add(
                egui::TextEdit::singleline(&mut preferred_csv)
                    .desired_width(360.0)
                    .hint_text("EcoRI, BamHI, SmaI"),
            )
            .on_hover_text("Comma-separated REBASE enzyme names used by preferred restriction-site display modes.")
            .changed()
        {
            self.configuration_graphics.preferred_restriction_enzymes =
                crate::dna_display::DnaDisplay::parse_preferred_restriction_enzymes_csv(
                    &preferred_csv,
                );
            changed = true;
        }
        ui.horizontal(|ui| {
            if ui.small_button("Use pUC MCS defaults").clicked() {
                self.configuration_graphics.restriction_enzyme_display_mode =
                    RestrictionEnzymeDisplayMode::PreferredOnly;
                self.configuration_graphics.preferred_restriction_enzymes =
                    crate::enzymes::default_preferred_restriction_enzyme_names();
                changed = true;
            }
            if ui.small_button("Use Golden Gate Type IIS").clicked() {
                self.configuration_graphics.restriction_enzyme_display_mode =
                    RestrictionEnzymeDisplayMode::PreferredOnly;
                self.configuration_graphics.preferred_restriction_enzymes =
                    crate::enzymes::golden_gate_type_iis_preferred_restriction_enzyme_names();
                changed = true;
            }
            if ui.small_button("Clear preferred").clicked() {
                self.configuration_graphics
                    .preferred_restriction_enzymes
                    .clear();
                changed = true;
            }
        });
        changed |= ui
            .checkbox(
                &mut self.configuration_graphics.show_gc_contents,
                crate::i18n::tr("configuration.graphics.show_gc_contents"),
            )
            .changed();
        ui.horizontal(|ui| {
            ui.label(crate::i18n::tr("sequence.gc_bin_size"));
            if ui
                .add(
                    egui::DragValue::new(&mut self.configuration_graphics.gc_content_bin_size_bp)
                        .range(1..=5_000_000)
                        .speed(10.0)
                        .suffix(" bp"),
                )
                .on_hover_text(
                    "Bin size used for GC-content aggregation in linear/circular maps and SVG export",
                )
                .changed()
            {
                changed = true;
            }
        });
        changed |= ui
            .checkbox(
                &mut self.configuration_graphics.show_open_reading_frames,
                crate::i18n::tr("configuration.graphics.show_orfs"),
            )
            .changed();
        changed |= ui
            .checkbox(
                &mut self.configuration_graphics.show_methylation_sites,
                crate::i18n::tr("configuration.graphics.show_methylation"),
            )
            .changed();

        ui.separator();
        ui.heading(self.tr("configuration.graphics.window_styling"));
        backdrop_changed |= ui
            .checkbox(
                &mut self.configuration_window_backdrops.enabled,
                crate::i18n::tr("configuration.window_styling.enable_backdrops"),
            )
            .on_hover_text("Apply subtle per-window-type color/image styling")
            .changed();
        backdrop_changed |= ui
            .checkbox(
                &mut self.configuration_window_backdrops.draw_images,
                crate::i18n::tr("configuration.window_styling.use_images"),
            )
            .on_hover_text("Render optional image watermark per window type")
            .changed();
        backdrop_changed |= ui
            .checkbox(
                &mut self.configuration_window_backdrops.show_text_watermark,
                crate::i18n::tr("configuration.window_styling.show_text_watermark"),
            )
            .on_hover_text("Show low-contrast type text when no image is configured")
            .changed();
        ui.horizontal(|ui| {
            ui.label(crate::i18n::tr("configuration.window_styling.tint_opacity"));
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
            ui.label(crate::i18n::tr(
                "configuration.window_styling.image_opacity",
            ));
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
        let preload_status =
            window_backdrop::window_backdrop_preload_status(&self.window_backdrops);
        if preload_status.active {
            ui.small(format!(
                "Backdrop preload status: {}/{} image(s) queued in texture cache.",
                preload_status.queued_images, preload_status.configured_images
            ));
        } else {
            ui.small("Backdrop preload status: inactive (images are currently disabled).");
        }
        Self::render_window_backdrop_path_table(
            ui,
            &mut self.configuration_window_backdrops,
            &mut backdrop_changed,
            &mut self.window_backdrop_path_status_cache,
        );

        if changed {
            self.configuration_graphics_dirty = true;
        }
        if live_font_changed {
            self.apply_configuration_graphics_to_engine_state();
            self.configuration_graphics_dirty = true;
            let refreshed = self.refresh_open_sequence_windows(ui.ctx());
            self.configuration_status = format!(
                "Applied font settings live; refreshed {} open sequence window(s)",
                refreshed
            );
        }
        if backdrop_changed {
            self.configuration_window_backdrops_dirty = true;
        }

        ui.horizontal_wrapped(|ui| {
            if ui
                .button(self.tr("configuration.graphics.reset_defaults"))
                .on_hover_text("Reset graphics settings to built-in defaults")
                .clicked()
            {
                self.reset_configuration_graphics_to_defaults();
            }
            if ui
                .add_enabled(
                    self.configuration_graphics_dirty,
                    egui::Button::new(self.tr("configuration.graphics.apply_graphics")),
                )
                .on_hover_text("Apply graphics settings to project state")
                .clicked()
            {
                self.apply_configuration_graphics();
            }
            if ui
                .button(self.tr("configuration.graphics.apply_refresh"))
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
                .button(self.tr("configuration.graphics.reset_window_styling"))
                .on_hover_text("Reset themed backdrop settings to defaults")
                .clicked()
            {
                self.reset_configuration_window_backdrops_to_defaults();
            }
            if ui
                .add_enabled(
                    self.configuration_window_backdrops_dirty,
                    egui::Button::new(self.tr("configuration.graphics.apply_window_styling")),
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

    fn render_unsaved_changes_dialog(&mut self, ctx: &egui::Context) {
        if self.pending_project_action.is_none() {
            return;
        }
        let mut open = true;
        let spec = crate::egui_compat::ModalWindowSpec::new(
            "Unsaved Changes",
            Self::unsaved_changes_dialog_id(),
        );
        crate::egui_compat::show_modal_window(ctx, &spec, &mut open, |ui| {
            ui.label("Save changes to the current project before continuing?");
            ui.horizontal(|ui| {
                if ui
                    .button("Save")
                    .on_hover_text("Save current project, then continue")
                    .clicked()
                    && self.save_current_project()
                    && let Some(action) = self.pending_project_action.take()
                {
                    self.execute_project_action(action);
                }
                if ui
                    .button("Don't Save")
                    .on_hover_text("Continue without saving current project changes")
                    .clicked()
                    && let Some(action) = self.pending_project_action.take()
                {
                    self.execute_project_action(action);
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
        if !open {
            self.pending_project_action = None;
        }
    }

    fn render_about_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_about_dialog {
            return;
        }
        let mut open = self.show_about_dialog;
        let spec = crate::egui_compat::ModalWindowSpec::new(
            "About GENtle",
            egui::Id::new("about_gentle_modal"),
        );
        crate::egui_compat::show_modal_window(ctx, &spec, &mut open, |ui| {
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
        self.show_about_dialog = open;
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
        let viewport_id = Self::command_palette_viewport_id();
        let render_command_palette_in_foreground = self.viewport_foreground_requested(viewport_id);
        let builder = egui::ViewportBuilder::default()
            .with_title("Command Palette")
            .with_inner_size([760.0, 520.0])
            .with_min_inner_size([500.0, 320.0]);
        if ctx.embed_viewports() {
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
                if ui.input(|i| i.key_pressed(Key::ArrowDown)) && !entries.is_empty() {
                    self.command_palette_selected =
                        (self.command_palette_selected + 1) % entries.len();
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
                        .id_salt("command_palette_results_scroll")
                        .max_height(400.0)
                        .show(ui, |ui| {
                            scroll_input_policy::apply_scrollarea_keyboard_navigation(
                                ui,
                                scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                            );
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
            let spec = crate::egui_compat::HostedWindowSpec::new(
                "Command Palette",
                egui::Id::new("Command Palette"),
                Vec2::new(760.0, 520.0),
                Vec2::new(500.0, 320.0),
            )
            .foreground(render_command_palette_in_foreground);
            crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| render_contents(ui));
            self.clear_viewport_foreground_request_after_render(viewport_id);

            if let Some(action) = execute_action {
                self.execute_command_palette_action(ctx, action);
                open = false;
            }

            if ctx.input(|i| i.key_pressed(Key::Escape)) {
                open = false;
            }

            self.show_command_palette_dialog = open;
            return;
        }
        ctx.show_viewport_immediate(viewport_id, builder, |ctx, class| {
            self.note_viewport_focus_if_active(ctx, viewport_id);
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
                if ui.input(|i| i.key_pressed(Key::ArrowDown)) && !entries.is_empty() {
                    self.command_palette_selected =
                        (self.command_palette_selected + 1) % entries.len();
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
                        .id_salt("command_palette_results_scroll")
                        .max_height(400.0)
                        .show(ui, |ui| {
                            scroll_input_policy::apply_scrollarea_keyboard_navigation(
                                ui,
                                scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                            );
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

            if class == egui::ViewportClass::EmbeddedWindow {
                let spec = crate::egui_compat::HostedWindowSpec::new(
                    "Command Palette",
                    egui::Id::new("Command Palette"),
                    Vec2::new(760.0, 520.0),
                    Vec2::new(500.0, 320.0),
                )
                .foreground(render_command_palette_in_foreground);
                crate::egui_compat::show_hosted_window(&mut *ctx, &spec, &mut open, |ui| {
                    render_contents(ui)
                });
            } else {
                crate::egui_compat::show_central_panel(
                    &mut *ctx,
                    egui::CentralPanel::default(),
                    |ui| {
                        render_contents(ui);
                    },
                );
                if Self::viewport_close_requested_or_shortcut(ctx) {
                    open = false;
                }
            }
        });

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
        let viewport_id = Self::background_jobs_viewport_id();
        let spec = self.hosted_window_spec_for_viewport(
            "Background Jobs",
            egui::Id::new(("background_jobs_hosted_window", viewport_id)),
            viewport_id,
            Vec2::new(760.0, 480.0),
            Vec2::new(520.0, 320.0),
        );
        crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
            ui.label("Centralized progress, cancellation, and completion summaries");
            ui.separator();

            ui.strong("Prepare Genome");
            let mut cancel_prepare_clicked = false;
            if self.genome_prepare_task.is_some() {
                if let Some(task) = &self.genome_prepare_task {
                    ui.horizontal(|ui| {
                        ui.add(egui::Spinner::new());
                        let genome_label = self
                            .genome_prepare_progress
                            .as_ref()
                            .map(|progress| progress.genome_id.clone())
                            .unwrap_or_else(|| self.genome_id.clone());
                        ui.label(format!(
                            "{} '{}' ({:.1}s)",
                            task.mode.progress_label(),
                            genome_label,
                            task.started.elapsed().as_secs_f32()
                        ));
                        if ui
                            .button("Cancel")
                            .on_hover_text("Request cancellation of genome prepare job")
                            .clicked()
                        {
                            cancel_prepare_clicked = true;
                        }
                    });
                }
                if let Some((current_step, completed_steps, total_steps, overall, eta)) =
                    self.prepare_step_summary()
                {
                    let mut step_summary = format!(
                        "Current step: {} ({}/{})",
                        current_step, completed_steps, total_steps
                    );
                    if let Some(eta) = eta {
                        step_summary
                            .push_str(&format!(" • ETA {}", Self::format_duration_compact(eta)));
                    }
                    ui.small(step_summary);
                    ui.add(
                        egui::ProgressBar::new(overall)
                            .show_percentage()
                            .text(format!("{}/{} steps", completed_steps, total_steps)),
                    );
                } else if let Some(progress) = &self.genome_prepare_progress {
                    ui.small(format!(
                        "Current phase: {} ({})",
                        progress.phase, progress.item
                    ));
                }
            } else {
                ui.horizontal(|ui| {
                        if self.genome_prepare_failure_recovery.is_some() {
                            ui.small("Last reindex needs reinstall");
                            if ui
                                .button("Reinstall...")
                                .on_hover_text(
                                    "Open the reinstall confirmation flow for the last inconsistent prepared-genome reindex.",
                                )
                                .clicked()
                            {
                                self.queue_prepare_failure_reinstall(
                                    PreparedGenomeReinstallDialogHost::Root,
                                );
                            }
                        } else {
                            ui.small("Idle");
                        }
                        if ui
                            .button("Retry")
                            .on_hover_text(
                                "Run the current prepare/reindex action again using the current dialog settings",
                            )
                            .clicked()
                        {
                            let snapshot_id = self.capture_retry_snapshot(
                                BackgroundJobKind::PrepareGenome,
                                "background jobs panel",
                                self.current_prepare_retry_arguments(),
                            );
                            self.push_job_event(
                                BackgroundJobKind::PrepareGenome,
                                BackgroundJobEventPhase::Retried,
                                None,
                                format!(
                                    "Retry requested from background jobs panel (retry snapshot #{snapshot_id})"
                                ),
                            );
                            self.start_prepare_reference_genome_for_current_selection();
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
            let mut cancel_blast_clicked = false;
            if let Some(task) = &self.genome_blast_task {
                ui.horizontal(|ui| {
                    ui.add(egui::Spinner::new());
                    ui.label(format!(
                        "Running ({:.1}s)",
                        task.started.elapsed().as_secs_f32()
                    ));
                    if task.cancel_requested.load(Ordering::Relaxed) {
                        ui.small("Cancellation requested...");
                    } else if ui
                        .button("Cancel")
                        .on_hover_text("Request cancellation of running BLAST job")
                        .clicked()
                    {
                        cancel_blast_clicked = true;
                    }
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
                            let snapshot_id = self.capture_retry_snapshot(
                                BackgroundJobKind::BlastGenome,
                                "background jobs panel",
                                self.current_blast_retry_arguments(),
                            );
                            self.push_job_event(
                                BackgroundJobKind::BlastGenome,
                                BackgroundJobEventPhase::Retried,
                                None,
                                format!(
                                    "Retry requested from background jobs panel (retry snapshot #{snapshot_id})"
                                ),
                            );
                            self.start_reference_genome_blast();
                        }
                    });
            }
            if cancel_blast_clicked {
                self.request_blast_task_cancel("background jobs panel");
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
                            let snapshot_id = self.capture_retry_snapshot(
                                BackgroundJobKind::TrackImport,
                                "background jobs panel",
                                self.current_track_import_retry_arguments(),
                            );
                            self.push_job_event(
                                BackgroundJobKind::TrackImport,
                                BackgroundJobEventPhase::Retried,
                                None,
                                format!(
                                    "Retry requested from background jobs panel (retry snapshot #{snapshot_id})"
                                ),
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
            ui.strong("Tutorial Project");
            let mut cancel_tutorial_clicked = false;
            if let Some(task) = &self.tutorial_project_task {
                ui.horizontal(|ui| {
                    ui.add(egui::Spinner::new());
                    ui.label(format!(
                        "Opening '{}' ({:.1}s)",
                        task.chapter_title,
                        task.started.elapsed().as_secs_f32()
                    ));
                    if task.cancel_requested.load(Ordering::Relaxed) {
                        ui.small("Cancellation requested...");
                    } else if ui
                        .button("Cancel")
                        .on_hover_text("Request cancellation of the running tutorial-project build")
                        .clicked()
                    {
                        cancel_tutorial_clicked = true;
                    }
                });
                if let Some(fraction) = self.tutorial_project_progress_fraction {
                    ui.add(
                        egui::ProgressBar::new(fraction.clamp(0.0, 1.0))
                            .show_percentage()
                            .text(self.tutorial_project_progress_label.clone()),
                    );
                } else if !self.tutorial_project_progress_label.trim().is_empty() {
                    ui.small(self.tutorial_project_progress_label.clone());
                }
            } else {
                ui.small("Idle");
            }
            if cancel_tutorial_clicked {
                self.request_tutorial_project_task_cancel("background jobs panel");
            }
            if !self.tutorial_project_status.trim().is_empty() {
                ui.small(self.tutorial_project_status.clone());
            }

            ui.separator();
            ui.strong("Agent Assistant");
            if let Some(task) = &self.agent_task {
                ui.horizontal(|ui| {
                    ui.add(egui::Spinner::new());
                    ui.label(format!(
                        "Running ({:.1}s)",
                        task.started.elapsed().as_secs_f32()
                    ));
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
                            let snapshot_id = self.capture_retry_snapshot(
                                BackgroundJobKind::AgentAssist,
                                "background jobs panel",
                                self.current_agent_retry_arguments(),
                            );
                            self.push_job_event(
                                BackgroundJobKind::AgentAssist,
                                BackgroundJobEventPhase::Retried,
                                None,
                                format!(
                                    "Retry requested from background jobs panel (retry snapshot #{snapshot_id})"
                                ),
                            );
                            self.start_agent_assistant_request();
                        }
                    });
            }
            if !self.agent_status.trim().is_empty() {
                self.render_agent_status_message(ui, &self.agent_status, false);
            }

            ui.separator();
            ui.strong("Recent job events");
            egui::ScrollArea::vertical()
                .id_salt(BACKGROUND_JOBS_RECENT_JOB_EVENTS_SCROLL_ID)
                .max_height(180.0)
                .show(ui, |ui| {
                    scroll_input_policy::apply_scrollarea_keyboard_navigation(
                        ui,
                        scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                    );
                    for event in self.job_event_log.iter().rev().take(40) {
                        ui.small(event.to_line());
                    }
                });

            ui.separator();
            ui.strong("Recent retry snapshots");
            self.retry_snapshot_retain_count = self.retry_snapshot_retention_limit();
            ui.horizontal(|ui| {
                    ui.label("Retain newest");
                    ui.add(
                        egui::DragValue::new(&mut self.retry_snapshot_retain_count)
                            .range(1..=MAX_BACKGROUND_JOB_RETRY_SNAPSHOTS)
                            .speed(1.0),
                    )
                    .on_hover_text(
                        "Retention cap for future retry captures (up to in-memory maximum)",
                    );
                    ui.small(format!("max {}", MAX_BACKGROUND_JOB_RETRY_SNAPSHOTS));
                    if ui
                        .button("Prune oldest")
                        .on_hover_text("Drop oldest snapshots until retained count is reached")
                        .clicked()
                    {
                        let before = self.retry_argument_snapshots.len();
                        let removed =
                            self.prune_retry_snapshots_to_limit(self.retry_snapshot_retain_count);
                        let after = self.retry_argument_snapshots.len();
                        self.retry_snapshot_pending_cleanup_action = None;
                        self.retry_snapshot_cleanup_confirm_text.clear();
                        self.app_status = if removed == 0 {
                            "Retry snapshot prune: no snapshots removed".to_string()
                        } else {
                            let audit_id = self.record_retry_snapshot_cleanup_audit(
                                "prune_oldest",
                                "Pruned oldest retry snapshots via retention control",
                                self.retry_snapshot_retain_count,
                                removed,
                                before,
                                after,
                                None,
                            );
                            format!(
                                "Retry snapshot prune: removed {removed} oldest snapshot(s) (audit #{audit_id})"
                            )
                        };
                    }
                    if ui
                        .add_enabled(
                            !self.retry_argument_snapshots.is_empty(),
                            egui::Button::new("Clear all"),
                        )
                        .on_hover_text("Delete all retained retry snapshots")
                        .clicked()
                    {
                        let before = self.retry_argument_snapshots.len();
                        let removed = self.clear_retry_snapshots();
                        let after = self.retry_argument_snapshots.len();
                        self.retry_snapshot_pending_cleanup_action = None;
                        self.retry_snapshot_cleanup_confirm_text.clear();
                        self.app_status = if removed == 0 {
                            "Cleared 0 retry snapshot(s)".to_string()
                        } else {
                            let audit_id = self.record_retry_snapshot_cleanup_audit(
                                "clear_all",
                                "Cleared all retained retry snapshots",
                                before,
                                removed,
                                before,
                                after,
                                None,
                            );
                            format!("Cleared {removed} retry snapshot(s) (audit #{audit_id})")
                        };
                    }
                });
            ui.horizontal(|ui| {
                    ui.label("Kind");
                    egui::ComboBox::from_id_salt("retry_snapshot_kind_filter")
                        .selected_text(self.retry_snapshot_kind_filter.label())
                        .show_ui(ui, |ui| {
                            ui.selectable_value(
                                &mut self.retry_snapshot_kind_filter,
                                RetrySnapshotKindFilter::All,
                                RetrySnapshotKindFilter::All.label(),
                            );
                            ui.selectable_value(
                                &mut self.retry_snapshot_kind_filter,
                                RetrySnapshotKindFilter::PrepareGenome,
                                RetrySnapshotKindFilter::PrepareGenome.label(),
                            );
                            ui.selectable_value(
                                &mut self.retry_snapshot_kind_filter,
                                RetrySnapshotKindFilter::BlastGenome,
                                RetrySnapshotKindFilter::BlastGenome.label(),
                            );
                            ui.selectable_value(
                                &mut self.retry_snapshot_kind_filter,
                                RetrySnapshotKindFilter::TrackImport,
                                RetrySnapshotKindFilter::TrackImport.label(),
                            );
                            ui.selectable_value(
                                &mut self.retry_snapshot_kind_filter,
                                RetrySnapshotKindFilter::OpenTutorialProject,
                                RetrySnapshotKindFilter::OpenTutorialProject.label(),
                            );
                            ui.selectable_value(
                                &mut self.retry_snapshot_kind_filter,
                                RetrySnapshotKindFilter::AgentAssist,
                                RetrySnapshotKindFilter::AgentAssist.label(),
                            );
                        });
                    ui.label("Filter");
                    ui.add(
                        egui::TextEdit::singleline(&mut self.retry_snapshot_text_filter)
                            .desired_width(260.0)
                            .hint_text("kind:blast origin:panel args:hg38"),
                    )
                    .on_hover_text(Self::retry_snapshot_filter_help_text());
                    if ui.button("Clear").clicked() {
                        self.retry_snapshot_text_filter.clear();
                    }
                    let filtered_count = self.filtered_retry_snapshots().len();
                    if ui
                        .add_enabled(filtered_count > 0, egui::Button::new("Export filtered..."))
                        .on_hover_text("Export filtered retry snapshots to JSON")
                        .clicked()
                    {
                        self.prompt_export_filtered_retry_snapshots();
                    }
                    if ui
                        .add_enabled(filtered_count > 0, egui::Button::new("Delete filtered"))
                        .on_hover_text("Delete only snapshots that match the current filters")
                        .clicked()
                    {
                        let action = RetrySnapshotPendingCleanupAction::DeleteFiltered;
                        self.retry_snapshot_pending_cleanup_action = Some(action);
                        self.retry_snapshot_cleanup_confirm_text.clear();
                        let confirm_phrase = action.confirm_phrase(filtered_count);
                        self.app_status = format!(
                            "Delete filtered staged: type '{confirm_phrase}' to confirm"
                        );
                    }
                    if ui
                        .add_enabled(
                            filtered_count > 0,
                            egui::Button::new("Archive & delete filtered..."),
                        )
                        .on_hover_text(
                            "Export filtered snapshots to JSON archive and remove them from retained history",
                        )
                        .clicked()
                    {
                        let action = RetrySnapshotPendingCleanupAction::ArchiveAndDeleteFiltered;
                        self.retry_snapshot_pending_cleanup_action = Some(action);
                        self.retry_snapshot_cleanup_confirm_text.clear();
                        let confirm_phrase = action.confirm_phrase(filtered_count);
                        self.app_status = format!(
                            "Archive & delete filtered staged: type '{confirm_phrase}' to confirm"
                        );
                    }
                });
            let dry_run_diff = self.retry_snapshot_dry_run_diff();
            let filtered_snapshots = &dry_run_diff.removed;
            if filtered_snapshots.is_empty() {
                self.retry_snapshot_pending_cleanup_action = None;
                self.retry_snapshot_cleanup_confirm_text.clear();
            }
            if let Some(action) = self.retry_snapshot_pending_cleanup_action {
                let preview_summary =
                    Self::summarize_retry_snapshot_cleanup_targets(filtered_snapshots);
                let confirm_phrase = action.confirm_phrase(filtered_snapshots.len());
                ui.group(|ui| {
                    ui.label(format!("Pending confirmation: {}", action.label()));
                    ui.small(preview_summary.clone());
                    ui.small(dry_run_diff.summary_line());
                    ui.small(format!("Type '{}' to confirm", confirm_phrase));
                    ui.horizontal(|ui| {
                        ui.label("Confirm phrase");
                        ui.add(
                            egui::TextEdit::singleline(
                                &mut self.retry_snapshot_cleanup_confirm_text,
                            )
                            .desired_width(220.0)
                            .hint_text(confirm_phrase.as_str()),
                        );
                    });
                    let confirm_ready = Self::retry_snapshot_cleanup_confirm_input_matches(
                        action,
                        filtered_snapshots.len(),
                        &self.retry_snapshot_cleanup_confirm_text,
                    );
                    if !confirm_ready {
                        ui.small("Confirmation phrase does not match yet");
                    }
                    ui.columns(2, |columns| {
                        columns[0].small("Would be removed");
                        egui::ScrollArea::vertical()
                            .id_salt(BACKGROUND_JOBS_RETRY_SNAPSHOTS_REMOVED_PREVIEW_SCROLL_ID)
                            .max_height(120.0)
                            .show(&mut columns[0], |ui| {
                                if filtered_snapshots.is_empty() {
                                    ui.small("No snapshots match current filters");
                                } else {
                                    for snapshot in filtered_snapshots.iter().rev().take(12) {
                                        ui.small(snapshot.to_line());
                                    }
                                    if filtered_snapshots.len() > 12 {
                                        ui.small(format!(
                                            "... and {} more",
                                            filtered_snapshots.len() - 12
                                        ));
                                    }
                                }
                            });
                        columns[1].small("Would remain");
                        egui::ScrollArea::vertical()
                            .id_salt(BACKGROUND_JOBS_RETRY_SNAPSHOTS_RETAINED_PREVIEW_SCROLL_ID)
                            .max_height(120.0)
                            .show(&mut columns[1], |ui| {
                                if dry_run_diff.retained.is_empty() {
                                    ui.small("No snapshots would remain");
                                } else {
                                    for snapshot in dry_run_diff.retained.iter().rev().take(12) {
                                        ui.small(snapshot.to_line());
                                    }
                                    if dry_run_diff.retained.len() > 12 {
                                        ui.small(format!(
                                            "... and {} more",
                                            dry_run_diff.retained.len() - 12
                                        ));
                                    }
                                }
                            });
                    });
                    ui.horizontal(|ui| {
                        if ui
                            .add_enabled(
                                !filtered_snapshots.is_empty() && confirm_ready,
                                egui::Button::new("Confirm"),
                            )
                            .clicked()
                        {
                            let before = self.retry_argument_snapshots.len();
                            match action {
                                RetrySnapshotPendingCleanupAction::DeleteFiltered => {
                                    let removed = self.remove_filtered_retry_snapshots();
                                    let after = self.retry_argument_snapshots.len();
                                    if removed > 0 {
                                        let audit_id = self.record_retry_snapshot_cleanup_audit(
                                            "delete_filtered",
                                            &preview_summary,
                                            filtered_snapshots.len(),
                                            removed,
                                            before,
                                            after,
                                            None,
                                        );
                                        self.app_status = format!(
                                            "{} (audit #{}); {}",
                                            Self::format_retry_snapshot_cleanup_status(
                                                action, removed, before, after
                                            ),
                                            audit_id,
                                            preview_summary
                                        );
                                    } else {
                                        self.app_status =
                                            "Delete filtered confirmed: no snapshots removed"
                                                .to_string();
                                    }
                                }
                                RetrySnapshotPendingCleanupAction::ArchiveAndDeleteFiltered => {
                                    let archive_result =
                                        self.prompt_archive_and_delete_filtered_retry_snapshots();
                                    if let Some((removed, archive_path)) = archive_result {
                                        let after = self.retry_argument_snapshots.len();
                                        if removed > 0 {
                                            let audit_id = self
                                                .record_retry_snapshot_cleanup_audit(
                                                    "archive_delete_filtered",
                                                    &preview_summary,
                                                    filtered_snapshots.len(),
                                                    removed,
                                                    before,
                                                    after,
                                                    Some(archive_path.as_str()),
                                                );
                                            self.app_status = format!(
                                                "{} (audit #{}); {}",
                                                Self::format_retry_snapshot_cleanup_status(
                                                    action, removed, before, after
                                                ),
                                                audit_id,
                                                preview_summary
                                            );
                                        }
                                    }
                                }
                            }
                            self.retry_snapshot_pending_cleanup_action = None;
                            self.retry_snapshot_cleanup_confirm_text.clear();
                        }
                        if ui.button("Cancel").clicked() {
                            self.retry_snapshot_pending_cleanup_action = None;
                            self.retry_snapshot_cleanup_confirm_text.clear();
                            self.app_status =
                                "Retry snapshot cleanup confirmation canceled".to_string();
                        }
                    });
                });
            }
            ui.small(format!(
                "Showing {} of {} retained snapshots",
                filtered_snapshots.len(),
                self.retry_argument_snapshots.len()
            ));
            egui::ScrollArea::vertical()
                .id_salt(BACKGROUND_JOBS_RETRY_SNAPSHOTS_SCROLL_ID)
                .max_height(160.0)
                .show(ui, |ui| {
                    scroll_input_policy::apply_scrollarea_keyboard_navigation(
                        ui,
                        scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                    );
                    if self.retry_argument_snapshots.is_empty() {
                        ui.small("No retry snapshots captured yet");
                    } else if filtered_snapshots.is_empty() {
                        ui.small("No retry snapshots match current filters");
                    } else {
                        for snapshot in filtered_snapshots.iter().rev().take(40) {
                            ui.small(snapshot.to_line());
                        }
                    }
                });
            ui.separator();
            ui.strong("Retry cleanup audit");
            if self.retry_snapshot_cleanup_audit.is_empty() {
                self.retry_cleanup_audit_pending_clear_all = false;
                self.retry_cleanup_audit_clear_confirm_text.clear();
            }
            self.retry_cleanup_audit_retain_count = self.retry_cleanup_audit_retention_limit();
            ui.horizontal(|ui| {
                ui.label("Retain newest");
                ui.add(
                    egui::DragValue::new(&mut self.retry_cleanup_audit_retain_count)
                        .range(1..=MAX_RETRY_SNAPSHOT_CLEANUP_AUDIT_ENTRIES)
                        .speed(1.0),
                )
                .on_hover_text("Retention cap for cleanup-audit entries");
                ui.small(format!("max {}", MAX_RETRY_SNAPSHOT_CLEANUP_AUDIT_ENTRIES));
                if ui
                    .button("Prune oldest")
                    .on_hover_text(
                        "Drop oldest cleanup-audit entries until retained count is reached",
                    )
                    .clicked()
                {
                    let removed = self
                        .prune_retry_cleanup_audit_to_limit(self.retry_cleanup_audit_retain_count);
                    self.app_status = if removed == 0 {
                        "Retry cleanup audit prune: no entries removed".to_string()
                    } else {
                        format!("Retry cleanup audit prune: removed {removed} oldest entries")
                    };
                }
                if ui
                    .add_enabled(
                        !self.retry_snapshot_cleanup_audit.is_empty(),
                        egui::Button::new("Clear all"),
                    )
                    .on_hover_text("Delete all retained cleanup-audit entries")
                    .clicked()
                {
                    let target_count = self.retry_snapshot_cleanup_audit.len();
                    let confirm_phrase =
                        Self::retry_cleanup_audit_clear_confirm_phrase(target_count);
                    self.retry_cleanup_audit_pending_clear_all = true;
                    self.retry_cleanup_audit_clear_confirm_text.clear();
                    self.app_status = format!(
                        "Retry cleanup audit clear-all staged: type '{confirm_phrase}' to confirm"
                    );
                }
            });
            ui.horizontal(|ui| {
                    ui.label("Action");
                    egui::ComboBox::from_id_salt("retry_cleanup_audit_action_filter")
                        .selected_text(self.retry_cleanup_audit_action_filter.label())
                        .show_ui(ui, |ui| {
                            ui.selectable_value(
                                &mut self.retry_cleanup_audit_action_filter,
                                RetryCleanupAuditActionFilter::All,
                                RetryCleanupAuditActionFilter::All.label(),
                            );
                            ui.selectable_value(
                                &mut self.retry_cleanup_audit_action_filter,
                                RetryCleanupAuditActionFilter::DeleteFiltered,
                                RetryCleanupAuditActionFilter::DeleteFiltered.label(),
                            );
                            ui.selectable_value(
                                &mut self.retry_cleanup_audit_action_filter,
                                RetryCleanupAuditActionFilter::ArchiveDeleteFiltered,
                                RetryCleanupAuditActionFilter::ArchiveDeleteFiltered.label(),
                            );
                            ui.selectable_value(
                                &mut self.retry_cleanup_audit_action_filter,
                                RetryCleanupAuditActionFilter::PruneOldest,
                                RetryCleanupAuditActionFilter::PruneOldest.label(),
                            );
                            ui.selectable_value(
                                &mut self.retry_cleanup_audit_action_filter,
                                RetryCleanupAuditActionFilter::ClearAll,
                                RetryCleanupAuditActionFilter::ClearAll.label(),
                            );
                        });
                    ui.label("Filter");
                    ui.add(
                        egui::TextEdit::singleline(&mut self.retry_cleanup_audit_text_filter)
                            .desired_width(260.0)
                            .hint_text("action:delete kind:blast summary:pruned"),
                    )
                    .on_hover_text(Self::retry_cleanup_audit_filter_help_text());
                    if ui.button("Clear").clicked() {
                        self.retry_cleanup_audit_text_filter.clear();
                    }
                    let filtered_audit_count = self.filtered_retry_cleanup_audit_entries().len();
                    if ui
                        .add_enabled(
                            filtered_audit_count > 0,
                            egui::Button::new("Export filtered report..."),
                        )
                        .on_hover_text(
                            "Export filtered cleanup-audit report JSON (read-only; does not append new audit entries)",
                        )
                        .clicked()
                    {
                        self.prompt_export_retry_cleanup_audit_report();
                    }
                });
            let filtered_audit_entries = self.filtered_retry_cleanup_audit_entries();
            if self.retry_cleanup_audit_pending_clear_all {
                let target_count = self.retry_snapshot_cleanup_audit.len();
                let confirm_phrase = Self::retry_cleanup_audit_clear_confirm_phrase(target_count);
                let confirm_ready = Self::retry_cleanup_audit_clear_confirm_input_matches(
                    target_count,
                    &self.retry_cleanup_audit_clear_confirm_text,
                );
                ui.group(|ui| {
                    ui.label("Pending confirmation: Clear all cleanup-audit entries");
                    ui.small(format!(
                        "This will remove all {} retained cleanup-audit entries",
                        target_count
                    ));
                    ui.small(format!("Type '{}' to confirm", confirm_phrase));
                    ui.horizontal(|ui| {
                        ui.label("Confirm phrase");
                        ui.add(
                            egui::TextEdit::singleline(
                                &mut self.retry_cleanup_audit_clear_confirm_text,
                            )
                            .desired_width(220.0)
                            .hint_text(confirm_phrase.as_str()),
                        );
                    });
                    if !confirm_ready {
                        ui.small("Confirmation phrase does not match yet");
                    }
                    ui.horizontal(|ui| {
                        if ui
                            .add_enabled(
                                target_count > 0 && confirm_ready,
                                egui::Button::new("Confirm clear all"),
                            )
                            .clicked()
                        {
                            let removed = self.clear_retry_cleanup_audit();
                            self.retry_cleanup_audit_pending_clear_all = false;
                            self.retry_cleanup_audit_clear_confirm_text.clear();
                            self.app_status =
                                format!("Cleared {removed} retry cleanup audit entries");
                        }
                        if ui.button("Cancel").clicked() {
                            self.retry_cleanup_audit_pending_clear_all = false;
                            self.retry_cleanup_audit_clear_confirm_text.clear();
                            self.app_status =
                                "Retry cleanup audit clear-all confirmation canceled".to_string();
                        }
                    });
                });
            }
            ui.small(format!(
                "Showing {} of {} retained audit entries",
                filtered_audit_entries.len(),
                self.retry_snapshot_cleanup_audit.len()
            ));
            egui::ScrollArea::vertical()
                .id_salt(BACKGROUND_JOBS_RETRY_CLEANUP_AUDIT_SCROLL_ID)
                .max_height(120.0)
                .show(ui, |ui| {
                    scroll_input_policy::apply_scrollarea_keyboard_navigation(
                        ui,
                        scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                    );
                    if self.retry_snapshot_cleanup_audit.is_empty() {
                        ui.small("No retry cleanup audit entries yet");
                    } else if filtered_audit_entries.is_empty() {
                        ui.small("No retry cleanup audit entries match current filters");
                    } else {
                        for entry in filtered_audit_entries.iter().rev().take(40) {
                            ui.small(entry.to_line());
                        }
                    }
                });
        });
        self.clear_viewport_foreground_request_after_render(viewport_id);
        self.finalize_viewport_open_probe(viewport_id, "Background Jobs");
        self.show_jobs_panel = open;
    }

    fn render_status_bar(&mut self, ctx: &egui::Context) {
        let history_summary = self
            .engine
            .read()
            .map(|engine| engine.history_summary())
            .unwrap_or_default();
        crate::egui_compat::show_bottom_panel(
            ctx,
            egui::Id::new("gentle_status_bar"),
            egui::Panel::bottom("gentle_status_bar"),
            |ui| {
                ui.horizontal_wrapped(|ui| {
                    let hover_text = if self.hover_status_name.trim().is_empty() {
                        "Hover: -".to_string()
                    } else {
                        format!("Hover: {}", self.hover_status_name)
                    };
                    ui.monospace(hover_text);
                    ui.separator();
                    ui.small(format!(
                        "Undo {} / Redo {}",
                        history_summary.undo_count, history_summary.redo_count
                    ));
                    if !self.app_status.trim().is_empty() {
                        ui.separator();
                        ui.small(self.app_status.clone());
                    }
                });
            },
        );
    }

    fn help_content_width_requires_relayout(previous_width: f32, current_width: f32) -> bool {
        if !current_width.is_finite() || current_width <= 0.0 {
            return false;
        }
        if !previous_width.is_finite() || previous_width <= 0.0 {
            return true;
        }
        (previous_width - current_width).abs() > HELP_MARKDOWN_REFLOW_DELTA_PX
    }

    fn help_markdown_max_image_width(available_width: f32) -> usize {
        let safe_width = if available_width.is_finite() {
            available_width.max(120.0)
        } else {
            480.0
        };
        (safe_width * 0.75).round().clamp(220.0, 1600.0) as usize
    }

    fn clamp_help_topic_combo_width(available_width: f32) -> f32 {
        if available_width.is_finite() {
            (available_width * 0.55).clamp(180.0, 420.0)
        } else {
            280.0
        }
    }

    fn reconcile_embedded_window_open_state(content_open: bool, window_open: bool) -> bool {
        content_open && window_open
    }

    fn clamp_configuration_backdrop_path_field_width(available_width: f32) -> f32 {
        if available_width.is_finite() {
            (available_width * 0.45).clamp(180.0, 260.0)
        } else {
            240.0
        }
    }

    fn render_help_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_help_dialog {
            return;
        }
        let viewport_id = Self::help_viewport_id();
        let title = format!("Help - {}", self.active_help_title());
        if ctx.embed_viewports() {
            let render_started = Instant::now();
            let min_size = Vec2::new(420.0, 320.0);
            let mut open = self.show_help_dialog;
            let render_help_in_foreground = self.viewport_foreground_requested(viewport_id);
            let spec = crate::egui_compat::HostedWindowSpec::new(
                "Help",
                Self::hosted_help_window_id(),
                Vec2::new(860.0, 680.0),
                min_size,
            )
            .foreground(render_help_in_foreground)
            .legacy_layer_id(crate::egui_compat::hosted_window_title_layer_id(
                title.as_str(),
            ))
            .legacy_layer_id(egui::LayerId::new(
                egui::Order::Middle,
                egui::Id::new(Self::help_viewport_id()),
            ));
            crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
                self.render_help_contents(ui);
            });
            self.clear_viewport_foreground_request_after_render(viewport_id);
            self.note_slow_open_phase(
                viewport_id,
                "Help first-frame render",
                render_started.elapsed().as_millis(),
            );
            self.show_help_dialog =
                Self::reconcile_embedded_window_open_state(self.show_help_dialog, open);
            self.finalize_viewport_open_probe(viewport_id, "Help");
            return;
        }
        let builder = egui::ViewportBuilder::default()
            .with_title(title.clone())
            .with_inner_size([860.0, 680.0])
            .with_min_inner_size([420.0, 320.0]);
        ctx.show_viewport_immediate(viewport_id, builder, |ctx, class| {
            self.note_viewport_focus_if_active(ctx, viewport_id);
            if class == egui::ViewportClass::EmbeddedWindow {
                let render_started = Instant::now();
                let min_size = Vec2::new(420.0, 320.0);
                let mut open = self.show_help_dialog;
                let spec = crate::egui_compat::HostedWindowSpec::new(
                    "Help",
                    Self::hosted_help_window_id(),
                    Vec2::new(860.0, 680.0),
                    min_size,
                )
                .legacy_layer_id(crate::egui_compat::hosted_window_title_layer_id(
                    title.as_str(),
                ))
                .legacy_layer_id(egui::LayerId::new(
                    egui::Order::Middle,
                    egui::Id::new(Self::help_viewport_id()),
                ));
                crate::egui_compat::show_hosted_window(&mut *ctx, &spec, &mut open, |ui| {
                    self.render_help_contents(ui);
                });
                self.note_slow_open_phase(
                    viewport_id,
                    "Help first-frame render",
                    render_started.elapsed().as_millis(),
                );
                if Self::viewport_close_requested_or_shortcut(ctx) {
                    open = false;
                }
                self.show_help_dialog =
                    Self::reconcile_embedded_window_open_state(self.show_help_dialog, open);
                return;
            }

            let render_started = Instant::now();
            crate::egui_compat::show_central_panel(
                &mut *ctx,
                egui::CentralPanel::default().frame(egui::Frame::NONE),
                |ui| {
                    self.render_help_contents(ui);
                },
            );
            self.note_slow_open_phase(
                viewport_id,
                "Help first-frame render",
                render_started.elapsed().as_millis(),
            );

            if Self::viewport_close_requested_or_shortcut(ctx) {
                self.show_help_dialog = false;
            }
        });
        self.finalize_viewport_open_probe(viewport_id, "Help");
    }

    fn render_help_contents(&mut self, ui: &mut Ui) {
        window_backdrop::paint_window_backdrop(
            ui,
            WindowBackdropKind::Help,
            &self.window_backdrops,
        );
        with_window_content_inset(ui, |ui| {
            let find_shortcut = KeyboardShortcut::new(Modifiers::COMMAND, Key::F);
            if ui.ctx().input_mut(|i| i.consume_shortcut(&find_shortcut)) {
                self.help_focus_search_box = true;
            }

            let mut active_doc_changed = false;
            ui.horizontal_wrapped(|ui| {
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
                    .selectable_label(self.help_doc == HelpDoc::AgentInterface, "Agent Interface")
                    .clicked()
                {
                    self.help_doc = HelpDoc::AgentInterface;
                    active_doc_changed = true;
                }
                if ui
                    .selectable_label(
                        self.help_doc == HelpDoc::ReviewerPreview,
                        "Reviewer Quickstart",
                    )
                    .clicked()
                {
                    self.help_doc = HelpDoc::ReviewerPreview;
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
                if ui
                    .selectable_label(self.help_doc == HelpDoc::Tutorial, "Tutorial")
                    .on_hover_text("Open tutorial markdown docs")
                    .clicked()
                {
                    if self.help_tutorial_entries.is_empty() {
                        self.help_tutorial_entries = Self::discover_help_tutorial_entries();
                    }
                    self.set_help_tutorial_selected(self.help_tutorial_selected);
                    self.help_doc = HelpDoc::Tutorial;
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
                    .button("Copy Page")
                    .on_hover_text("Copy the current help/tutorial markdown source to the clipboard")
                    .clicked()
                {
                    ui.ctx().copy_text(self.active_help_copyable_text());
                }
                ui.separator();
                ui.label("View:");
                if ui
                    .selectable_label(!self.help_selectable_text_mode, "Rendered")
                    .on_hover_text("Show the current document as rendered markdown")
                    .clicked()
                {
                    self.help_selectable_text_mode = false;
                }
                if ui
                    .selectable_label(self.help_selectable_text_mode, "Selectable Text")
                    .on_hover_text(
                        "Show the current document as selectable raw markdown so any label or output id can be copied",
                    )
                    .clicked()
                {
                    self.help_selectable_text_mode = true;
                }
                if ui
                    .button("Close")
                    .on_hover_text("Close help window")
                    .clicked()
                {
                    self.show_help_dialog = false;
                }
            });

            ui.horizontal_wrapped(|ui| {
                if self.help_doc == HelpDoc::Tutorial {
                    if self.help_tutorial_entries.is_empty() {
                        self.help_tutorial_entries = Self::discover_help_tutorial_entries();
                        self.set_help_tutorial_selected(self.help_tutorial_selected);
                    }
                    ui.label("Topic:");
                    let selected_tutorial_title = self
                        .help_tutorial_entries
                        .get(self.help_tutorial_selected)
                        .map(|entry| {
                            Self::tutorial_display_label(
                                entry.decimal_id.as_deref(),
                                None,
                                &entry.title,
                            )
                        })
                        .unwrap_or_else(|| "No tutorial docs found".to_string());
                    let mut selected_tutorial_index = None;
                    egui::ComboBox::from_id_salt("help_tutorial_topic")
                        .width(Self::clamp_help_topic_combo_width(ui.available_width()))
                        .selected_text(selected_tutorial_title)
                        .show_ui(ui, |ui| {
                            let mut last_group: Option<String> = None;
                            for (index, entry) in self.help_tutorial_entries.iter().enumerate() {
                                let group = Self::tutorial_audience_group_label(entry);
                                if last_group.as_deref() != Some(group.as_str()) {
                                    if last_group.is_some() {
                                        ui.separator();
                                    }
                                    ui.label(egui::RichText::new(&group).strong());
                                    last_group = Some(group);
                                }
                                if ui
                                    .selectable_label(
                                        index == self.help_tutorial_selected,
                                        Self::tutorial_display_label(
                                            entry.decimal_id.as_deref(),
                                            None,
                                            &entry.title,
                                        ),
                                    )
                                    .on_hover_text(format!(
                                        "{}\n{}",
                                        Self::help_tutorial_review_label(entry),
                                        entry.summary
                                    ))
                                    .clicked()
                                {
                                    selected_tutorial_index = Some(index);
                                }
                            }
                        });
                    if let Some(index) = selected_tutorial_index {
                        self.set_help_tutorial_selected(index);
                        active_doc_changed = true;
                    }
                    if ui
                        .button("Copy Feedback Context")
                        .on_hover_text(
                            "Copy tutorial id, source, workflow, artifact, version, platform, and current search context for issue reports",
                        )
                        .clicked()
                    {
                        ui.ctx().copy_text(self.selected_tutorial_feedback_context_text());
                    }
                    ui.separator();
                }
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
                self.help_markdown_cache = CommonMarkCache::default();
                self.help_last_content_width = 0.0;
                self.refresh_help_search_matches();
            }
            if let Some(current) = self.help_search_matches.get(self.help_search_selected) {
                ui.small(format!("Line {}: {}", current.line_number, current.snippet));
            }

            ui.separator();
            let help_body_size = ui.available_size_before_wrap();
            ui.allocate_ui_with_layout(
                Vec2::new(help_body_size.x.max(120.0), help_body_size.y.max(180.0)),
                egui::Layout::top_down(egui::Align::Min),
                |ui| {
                    let help_body_height = ui.available_height().max(180.0);
                    egui::ScrollArea::vertical()
                        .id_salt("help_body_scroll")
                        .auto_shrink([false, false])
                        .max_height(help_body_height)
                        .show(ui, |ui| {
                            scroll_input_policy::apply_scrollarea_keyboard_navigation(
                                ui,
                                scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                            );
                            let content_width = ui.available_width().max(120.0);
                            if Self::help_content_width_requires_relayout(
                                self.help_last_content_width,
                                content_width,
                            ) {
                                self.help_markdown_cache = CommonMarkCache::default();
                            }
                            self.help_last_content_width = content_width;
                            let max_image_width =
                                Self::help_markdown_max_image_width(content_width);
                            let active_doc = self.help_doc;
                            let raw_markdown = self
                                .help_selectable_text_mode
                                .then(|| self.help_source_markdown(active_doc).to_string());
                            let rendered_markdown = (!self.help_selectable_text_mode)
                                .then(|| self.rendered_help_markdown_for(active_doc));
                            ui.set_width(content_width);
                            ui.set_max_width(content_width);
                            ui.set_min_height(help_body_height);
                            ui.allocate_ui_with_layout(
                                Vec2::new(content_width, help_body_height.max(320.0)),
                                egui::Layout::top_down(egui::Align::Min),
                                |ui| {
                                    if self.help_selectable_text_mode {
                                        let mut copyable_text =
                                            raw_markdown.clone().unwrap_or_default();
                                        ui.add_sized(
                                            [content_width, ui.available_height().max(320.0)],
                                            egui::TextEdit::multiline(&mut copyable_text)
                                                .font(egui::TextStyle::Monospace)
                                                .desired_width(content_width)
                                                .lock_focus(true),
                                        );
                                    } else {
                                        CommonMarkViewer::new()
                                            .max_image_width(Some(max_image_width))
                                            .show(
                                                ui,
                                                &mut self.help_markdown_cache,
                                                rendered_markdown.as_deref().unwrap_or_default(),
                                            );
                                    }
                                },
                            );
                        });
                },
            );
        });
    }

    fn summarize_operation(op: &Operation) -> String {
        match op {
            Operation::LoadFile { path, as_id } => match as_id {
                Some(id) => format!("Load file: path={path}, as_id={id}"),
                None => format!("Load file: path={path}"),
            },
            Operation::CreateSequenceFromText {
                output_id,
                name,
                circular,
                ..
            } => format!(
                "Create inline {} sequence: output_id={}, name={}",
                if *circular { "circular" } else { "linear" },
                output_id.clone().unwrap_or_else(|| "-".to_string()),
                name.clone().unwrap_or_else(|| "-".to_string())
            ),
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
            Operation::ApplyGibsonAssemblyPlan { plan_json } => {
                let parsed: Result<GibsonAssemblyPlan, _> = serde_json::from_str(plan_json);
                match parsed {
                    Ok(plan) => {
                        let insert_seq_ids = if plan.fragments.is_empty() {
                            "-".to_string()
                        } else {
                            plan.fragments
                                .iter()
                                .map(|fragment| fragment.seq_id.clone())
                                .collect::<Vec<_>>()
                                .join(",")
                        };
                        let opening_mode = if plan.destination.opening.mode.trim().is_empty() {
                            "-".to_string()
                        } else {
                            plan.destination.opening.mode.trim().to_string()
                        };
                        let output = if plan.product.output_id_hint.trim().is_empty() {
                            "-".to_string()
                        } else {
                            plan.product.output_id_hint.trim().to_string()
                        };
                        format!(
                            "Gibson cloning: destination={}, inserts={}, opening={}, output={}",
                            plan.destination.seq_id, insert_seq_ids, opening_mode, output
                        )
                    }
                    Err(_) => "Gibson cloning".to_string(),
                }
            }
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
                    DisplayTarget::RepeatFeatures => "Repeat features",
                    DisplayTarget::ArrayFeatures => "Array features",
                    DisplayTarget::ConstructReasoningOverlay => "Construct reasoning overlay",
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
            Operation::ComputeDotplot {
                seq_id,
                reference_seq_id,
                span_start_0based,
                span_end_0based,
                reference_span_start_0based,
                reference_span_end_0based,
                mode,
                word_size,
                step_bp,
                max_mismatches,
                tile_bp,
                store_as,
            } => format!(
                "Compute dotplot: seq_id={}, reference_seq_id={}, mode={}, query_span={}..{}, reference_span={}..{}, word={}, step={}, mismatches={}, tile_bp={}, store_as={}",
                seq_id,
                reference_seq_id
                    .as_deref()
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .unwrap_or("-"),
                mode.as_str(),
                span_start_0based
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                span_end_0based
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                reference_span_start_0based
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                reference_span_end_0based
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                word_size,
                step_bp,
                max_mismatches,
                tile_bp
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                store_as
                    .as_deref()
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .unwrap_or("-")
            ),
            Operation::ComputeFlexibilityTrack {
                seq_id,
                span_start_0based,
                span_end_0based,
                model,
                bin_bp,
                smoothing_bp,
                store_as,
            } => format!(
                "Compute flexibility track: seq_id={}, model={}, span={}..{}, bin_bp={}, smoothing_bp={}, store_as={}",
                seq_id,
                model.as_str(),
                span_start_0based
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                span_end_0based
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                bin_bp,
                smoothing_bp
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                store_as
                    .as_deref()
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .unwrap_or("-")
            ),
            Operation::DeriveSplicingReferences {
                seq_id,
                span_start_0based,
                span_end_0based,
                seed_feature_id,
                scope,
                output_prefix,
            } => format!(
                "Derive splicing references: seq_id={}, span={}..{}, seed_feature_id={}, scope={}, output_prefix={}",
                seq_id,
                span_start_0based,
                span_end_0based,
                seed_feature_id
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "auto".to_string()),
                scope.as_str(),
                output_prefix
                    .as_deref()
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .unwrap_or("-")
            ),
            Operation::AlignSequences {
                query,
                target,
                query_seq_id,
                target_seq_id,
                query_span_start_0based,
                query_span_end_0based,
                target_span_start_0based,
                target_span_end_0based,
                mode,
                match_score,
                mismatch_score,
                gap_open,
                gap_extend,
            } => {
                let describe_target =
                    |target: &Option<SequenceScanTarget>, legacy_seq_id: &Option<String>| {
                        match target {
                            Some(SequenceScanTarget::SeqId { seq_id, .. }) => seq_id.clone(),
                            Some(SequenceScanTarget::InlineSequence { id_hint, .. }) => id_hint
                                .as_deref()
                                .map(str::trim)
                                .filter(|value| !value.is_empty())
                                .unwrap_or("inline_sequence")
                                .to_string(),
                            None => legacy_seq_id
                                .as_deref()
                                .map(str::trim)
                                .filter(|value| !value.is_empty())
                                .unwrap_or("-")
                                .to_string(),
                        }
                    };
                format!(
                    "Align sequences: query={}, target={}, mode={}, query_span={}..{}, target_span={}..{}, scores(match={}, mismatch={}, gap_open={}, gap_extend={})",
                    describe_target(query, query_seq_id),
                    describe_target(target, target_seq_id),
                    mode.as_str(),
                    query_span_start_0based
                        .map(|value| value.to_string())
                        .unwrap_or_else(|| "-".to_string()),
                    query_span_end_0based
                        .map(|value| value.to_string())
                        .unwrap_or_else(|| "-".to_string()),
                    target_span_start_0based
                        .map(|value| value.to_string())
                        .unwrap_or_else(|| "-".to_string()),
                    target_span_end_0based
                        .map(|value| value.to_string())
                        .unwrap_or_else(|| "-".to_string()),
                    match_score,
                    mismatch_score,
                    gap_open,
                    gap_extend
                )
            }
            Operation::ConfirmConstructReads {
                expected_seq_id,
                baseline_seq_id,
                read_seq_ids,
                trace_ids,
                targets,
                report_id,
                ..
            } => format!(
                "Sequencing confirmation: expected={}, baseline={}, reads={}, traces={}, targets={}, report_id={}",
                expected_seq_id,
                baseline_seq_id.as_deref().unwrap_or("-"),
                read_seq_ids.len(),
                trace_ids.len(),
                targets.len(),
                report_id.as_deref().unwrap_or("-")
            ),
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
            Operation::RenderDotplotSvg {
                seq_id,
                dotplot_id,
                path,
                flex_track_id,
                display_density_threshold,
                display_intensity_gain,
                overlay_x_axis_mode,
                overlay_anchor_exon,
            } => format!(
                "Render dotplot SVG: seq_id={}, dotplot_id={}, path={}, flex_track_id={}, display_threshold={}, intensity_gain={}, overlay_x_axis={}, overlay_anchor_exon={}",
                seq_id,
                dotplot_id,
                path,
                flex_track_id
                    .as_deref()
                    .map(str::trim)
                    .filter(|v| !v.is_empty())
                    .unwrap_or("-"),
                display_density_threshold
                    .map(|v| format!("{v:.3}"))
                    .unwrap_or_else(|| "default(0.000)".to_string()),
                display_intensity_gain
                    .map(|v| format!("{v:.3}"))
                    .unwrap_or_else(|| "default(1.000)".to_string()),
                overlay_x_axis_mode.as_str(),
                overlay_anchor_exon
                    .as_ref()
                    .map(|exon| exon.token())
                    .unwrap_or_else(|| "-".to_string())
            ),
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
                conditions,
            } => format!(
                "Render serial gel SVG: inputs={}, container_ids={}, arrangement_id={}, path={}, ladders={}, conditions={}",
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
                    .unwrap_or_else(|| "auto".to_string()),
                conditions
                    .as_ref()
                    .map(crate::engine::GelRunConditions::describe)
                    .unwrap_or_else(|| crate::engine::GelRunConditions::default().describe())
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
            Operation::SetArrangementLadders {
                arrangement_id,
                ladders,
            } => format!(
                "Set arrangement ladders: arrangement_id={}, ladders={}",
                arrangement_id.trim(),
                ladders
                    .as_ref()
                    .map(|v| v.join(", "))
                    .filter(|s| !s.is_empty())
                    .unwrap_or_else(|| "auto".to_string())
            ),
            Operation::SetContainerDeclaredContentsExclusive {
                container_id,
                exclusive,
            } => format!(
                "Set container contents mode: container_id={}, exclusive={}",
                container_id, exclusive
            ),
            Operation::CreateRackFromArrangement {
                arrangement_id,
                rack_id,
                name,
                profile,
            } => format!(
                "Create rack from arrangement: arrangement_id={}, rack_id={}, name={}, profile={}",
                arrangement_id.trim(),
                rack_id.as_deref().unwrap_or("auto"),
                name.as_deref().unwrap_or("-"),
                profile.map(|kind| kind.as_str()).unwrap_or("auto")
            ),
            Operation::PlaceArrangementOnRack {
                arrangement_id,
                rack_id,
            } => format!(
                "Place arrangement on rack: arrangement_id={}, rack_id={}",
                arrangement_id.trim(),
                rack_id.trim()
            ),
            Operation::MoveRackPlacement {
                rack_id,
                from_coordinate,
                to_coordinate,
                move_block,
            } => format!(
                "Move rack {}: rack_id={}, from={}, to={}",
                if *move_block { "block" } else { "placement" },
                rack_id.trim(),
                from_coordinate.trim(),
                to_coordinate.trim()
            ),
            Operation::MoveRackSamples {
                rack_id,
                from_coordinates,
                to_coordinate,
            } => format!(
                "Move rack samples: rack_id={}, from={}, to={}",
                rack_id.trim(),
                from_coordinates.join(", "),
                to_coordinate.trim()
            ),
            Operation::MoveRackArrangementBlocks {
                rack_id,
                arrangement_ids,
                to_coordinate,
            } => format!(
                "Move rack blocks: rack_id={}, arrangements={}, to={}",
                rack_id.trim(),
                arrangement_ids.join(", "),
                to_coordinate.trim()
            ),
            Operation::SetRackProfile { rack_id, profile } => format!(
                "Set rack profile: rack_id={}, profile={}",
                rack_id.trim(),
                profile.as_str()
            ),
            Operation::ApplyRackTemplate { rack_id, template } => format!(
                "Apply rack template: rack_id={}, template={}",
                rack_id.trim(),
                template.as_str()
            ),
            Operation::SetRackFillDirection {
                rack_id,
                fill_direction,
            } => format!(
                "Set rack fill direction: rack_id={}, fill_direction={}",
                rack_id.trim(),
                fill_direction.as_str()
            ),
            Operation::SetRackProfileCustom {
                rack_id,
                rows,
                columns,
            } => format!(
                "Set custom rack profile: rack_id={}, rows={}, columns={}",
                rack_id.trim(),
                rows,
                columns
            ),
            Operation::SetRackBlockedCoordinates {
                rack_id,
                blocked_coordinates,
            } => format!(
                "Set blocked rack coordinates: rack_id={}, blocked={}",
                rack_id.trim(),
                if blocked_coordinates.is_empty() {
                    "-".to_string()
                } else {
                    blocked_coordinates.join(", ")
                }
            ),
            Operation::ExportRackLabelsSvg {
                rack_id,
                path,
                arrangement_id,
                preset,
            } => format!(
                "Export rack labels SVG: rack_id={}, arrangement_id={}, preset={}, path={}",
                rack_id.trim(),
                arrangement_id.as_deref().unwrap_or("all"),
                preset.as_str(),
                path
            ),
            Operation::ExportRackFabricationSvg {
                rack_id,
                path,
                template,
            } => format!(
                "Export rack fabrication SVG: rack_id={}, template={}, path={}",
                rack_id.trim(),
                template.as_str(),
                path
            ),
            Operation::ExportRackIsometricSvg {
                rack_id,
                path,
                template,
            } => format!(
                "Export rack isometric SVG: rack_id={}, template={}, path={}",
                rack_id.trim(),
                template.as_str(),
                path
            ),
            Operation::ExportRackOpenScad {
                rack_id,
                path,
                template,
            } => format!(
                "Export rack OpenSCAD: rack_id={}, template={}, path={}",
                rack_id.trim(),
                template.as_str(),
                path
            ),
            Operation::ExportRackCarrierLabelsSvg {
                rack_id,
                path,
                arrangement_id,
                template,
                preset,
            } => format!(
                "Export rack carrier labels SVG: rack_id={}, arrangement_id={}, template={}, preset={}, path={}",
                rack_id.trim(),
                arrangement_id.as_deref().unwrap_or("all"),
                template.as_str(),
                preset.as_str(),
                path
            ),
            Operation::ExportRackSimulationJson {
                rack_id,
                path,
                template,
            } => format!(
                "Export rack simulation JSON: rack_id={}, template={}, path={}",
                rack_id.trim(),
                template.as_str(),
                path
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
            Operation::ExportLabAssistantInstructions {
                path,
                run_id,
                title,
                audience,
                format,
            } => format!(
                "Export lab assistant report: path={}, run_id={}, title={}, audience={}, format={}",
                path,
                run_id.as_deref().unwrap_or("all"),
                title.as_deref().unwrap_or("-"),
                audience.as_deref().unwrap_or("-"),
                format
                    .as_ref()
                    .copied()
                    .map(LabAssistantInstructionsFormat::as_str)
                    .unwrap_or("infer")
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
                annotation_scope,
                max_annotation_features,
                include_genomic_annotation,
                catalog_path,
                cache_dir,
            } => format!(
                "Extract genome region: genome_id={}, chromosome={}, start={}, end={}, output_id={}, annotation_scope={}, max_annotation_features={}, include_genomic_annotation={}, catalog_path={}, cache_dir={}",
                genome_id,
                chromosome,
                start_1based,
                end_1based,
                output_id.clone().unwrap_or_else(|| "-".to_string()),
                annotation_scope
                    .map(|scope| scope.as_str().to_string())
                    .unwrap_or_else(|| "-".to_string()),
                max_annotation_features
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                include_genomic_annotation
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                catalog_path.clone().unwrap_or_else(|| "-".to_string()),
                cache_dir.clone().unwrap_or_else(|| "-".to_string())
            ),
            Operation::ExtractGenomeGene {
                genome_id,
                gene_query,
                occurrence,
                output_id,
                extract_mode,
                promoter_upstream_bp,
                annotation_scope,
                max_annotation_features,
                include_genomic_annotation,
                catalog_path,
                cache_dir,
            } => format!(
                "Extract genome gene: genome_id={}, gene_query={}, occurrence={}, output_id={}, extract_mode={}, promoter_upstream_bp={}, annotation_scope={}, max_annotation_features={}, include_genomic_annotation={}, catalog_path={}, cache_dir={}",
                genome_id,
                gene_query,
                occurrence
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                output_id.clone().unwrap_or_else(|| "-".to_string()),
                extract_mode
                    .map(|mode| mode.as_str().to_string())
                    .unwrap_or_else(|| "-".to_string()),
                promoter_upstream_bp
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                annotation_scope
                    .map(|scope| scope.as_str().to_string())
                    .unwrap_or_else(|| "-".to_string()),
                max_annotation_features
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                include_genomic_annotation
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "-".to_string()),
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
                prepared_genome_id,
            } => format!(
                "Extend genome anchor: seq_id={}, side={:?}, length_bp={}, output_id={}, catalog_path={}, cache_dir={}, prepared_genome_id={}",
                seq_id,
                side,
                length_bp,
                output_id.clone().unwrap_or_else(|| "-".to_string()),
                catalog_path.clone().unwrap_or_else(|| "-".to_string()),
                cache_dir.clone().unwrap_or_else(|| "-".to_string()),
                prepared_genome_id
                    .clone()
                    .unwrap_or_else(|| "-".to_string())
            ),
            Operation::VerifyGenomeAnchor {
                seq_id,
                catalog_path,
                cache_dir,
                prepared_genome_id,
            } => format!(
                "Verify genome anchor: seq_id={}, catalog_path={}, cache_dir={}, prepared_genome_id={}",
                seq_id,
                catalog_path.clone().unwrap_or_else(|| "-".to_string()),
                cache_dir.clone().unwrap_or_else(|| "-".to_string()),
                prepared_genome_id
                    .clone()
                    .unwrap_or_else(|| "-".to_string())
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
            Operation::ProjectMicroarrayTrack {
                seq_id,
                manifest_path,
                contrasts,
                level,
                min_abs_logfc,
                max_adj_p,
                max_features,
                clear_existing,
            } => format!(
                "Project microarray track: seq_id={}, manifest={}, contrasts={}, level={}, min_abs_logfc={:?}, max_adj_p={:?}, max_features={:?}, clear_existing={}",
                seq_id,
                manifest_path,
                if contrasts.is_empty() {
                    "manifest-order".to_string()
                } else {
                    contrasts.join(",")
                },
                level.clone().unwrap_or_else(|| "probeset".to_string()),
                min_abs_logfc,
                max_adj_p,
                max_features,
                clear_existing.unwrap_or(false)
            ),
            Operation::ProjectGenomeInterval {
                source_genome_id,
                target_genome_id,
                projection_path,
                chrom,
                start_1based,
                end_1based,
                strand,
            } => format!(
                "Project genome interval: {}:{}-{}{} from {} to {} via {}",
                chrom,
                start_1based,
                end_1based,
                strand
                    .as_ref()
                    .map(|value| format!(" {value}"))
                    .unwrap_or_default(),
                source_genome_id,
                target_genome_id,
                projection_path
            ),
            Operation::ImportBlastHitsTrack {
                seq_id,
                hits,
                track_name,
                clear_existing,
                blast_provenance,
            } => {
                let invocation = blast_provenance
                    .as_ref()
                    .map(|entry| {
                        let raw = if entry.command_line.trim().is_empty() {
                            Self::format_blast_command_line(
                                &entry.blastn_executable,
                                &entry.command,
                            )
                        } else {
                            entry.command_line.clone()
                        };
                        Self::compact_lineage_node_label(&raw, 92)
                    })
                    .unwrap_or_else(|| "-".to_string());
                format!(
                    "Import BLAST hits track: seq_id={}, track_name={}, hits={}, clear_existing={}, invocation={}",
                    seq_id,
                    track_name.clone().unwrap_or_else(|| "-".to_string()),
                    hits.len(),
                    clear_existing.unwrap_or(false),
                    invocation
                )
            }
            other => format!("{other:?}"),
        }
    }
}

impl eframe::App for GENtleApp {
    fn ui(&mut self, ui: &mut egui::Ui, frame: &mut eframe::Frame) {
        crate::egui_compat::with_current_root_ui(ui, |ctx| {
            self.render_root_ui(ctx, frame);
        });
    }
}

impl GENtleApp {
    fn render_root_ui(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        let update_result = catch_unwind(AssertUnwindSafe(|| {
            Self::configure_platform_viewport_mode(ctx);
            if !self.update_has_run_before {
                egui_extras::install_image_loaders(ctx);
                self.update_has_run_before = true;
            }
            window_backdrop::preload_window_backdrop_images(ctx, &self.window_backdrops);
            about::install_native_help_menu_bridge();
            about::install_native_settings_menu_bridge();
            about::install_native_windows_menu_bridge();
            about::install_native_app_windows_menu_bridge();
            self.consume_native_help_request();
            self.consume_native_settings_request();
            self.consume_native_pcr_design_request();
            self.consume_native_jaspar_expert_request();
            self.consume_native_windows_request();
            self.consume_active_viewport_report();
            if ctx.input(|i| i.viewport().focused.unwrap_or(false)) {
                self.set_active_window_viewport(ViewportId::ROOT);
                self.finalize_viewport_focus_probe(ViewportId::ROOT);
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

            self.handle_root_keyboard_shortcuts(ctx);

            self.poll_prepare_reference_genome_task(ctx);
            self.poll_tutorial_project_task(ctx);
            self.poll_reference_genome_blast_task(ctx);
            self.poll_genome_track_import_task(ctx);
            self.poll_dbsnp_fetch_task(ctx);
            self.poll_agent_assistant_task(ctx);
            self.poll_agent_model_discovery_task(ctx);
            self.poll_clawbio_task(ctx);
            self.sync_tracked_bed_tracks_for_new_anchors();
            self.sync_open_windows_if_display_changed(ctx);
            self.reset_root_auxiliary_areas_if_legacy_title_layers_visible(ctx);

            if self.show_help_dialog {
                Self::reset_root_help_areas_if_legacy_layers_visible(
                    ctx,
                    format!("Help - {}", self.active_help_title()).as_str(),
                );
            }

            // Show menu bar
            crate::egui_compat::show_top_panel(
                ctx,
                egui::Id::new("top"),
                egui::Panel::top("top"),
                |ui| {
                    self.render_menu_bar(ui);
                },
            );

            // Show main/root workspace host, or the transient startup splash.
            if self.splash_should_render_at(Instant::now()) {
                self.render_splash_screen(ctx);
            } else {
                self.dismiss_splash_screen();
                self.render_root_workspace(ctx, project_dirty);
                self.render_reference_genome_prepare_dialog(ctx);
                self.render_reference_genome_retrieve_dialog(ctx);
                self.render_new_sequence_dialog(ctx);
                self.render_uniprot_dialog(ctx);
                self.render_genbank_dialog(ctx);
                self.render_reference_genome_blast_dialog(ctx);
                self.render_reference_genome_inspector_dialog(ctx);
                self.render_cache_cleanup_dialog(ctx);
                self.render_prepared_genome_reinstall_confirm_dialog(
                    ctx,
                    PreparedGenomeReinstallDialogHost::Root,
                );
                self.render_pending_ensembl_catalog_update_dialog(ctx);
                self.render_pending_ensembl_installable_genomes_dialog(ctx);
                self.render_pending_ensembl_quick_install_dialog(ctx);
                self.render_pending_prepared_genome_removal_dialog(ctx);
                self.render_pending_catalog_entry_removal_dialog(ctx);
                self.render_genome_bed_track_dialog(ctx);
                self.render_gibson_dialog(ctx);
                self.render_arrangement_gel_preview_dialog(ctx);
                self.render_rack_labels_preview_dialog(ctx);
                self.render_place_arrangement_on_rack_dialog(ctx);
                self.render_rack_dialog(ctx);
                self.render_pcr_design_dialog(ctx);
                self.render_sequencing_confirmation_dialog(ctx);
                self.render_jaspar_expert_dialog(ctx);
                self.render_planning_dialog(ctx);
                self.render_routine_assistant_dialog(ctx);
                self.render_agent_assistant_dialog(ctx);
                self.render_external_services_dialog(ctx);
                self.render_clawbio_dialog(ctx);
                self.render_configuration_dialog(ctx);
                self.render_help_dialog(ctx);
                self.render_about_dialog(ctx);
                self.render_command_palette_dialog(ctx);
                self.render_jobs_panel(ctx);
                self.render_history_panel(ctx);
                self.render_unsaved_changes_dialog(ctx);
            }
            self.render_status_bar(ctx);

            if self.pending_app_quit {
                ctx.send_viewport_cmd(egui::ViewportCommand::Close);
                self.pending_app_quit = false;
            }

            // Open new windows
            let mut new_windows: Vec<Window> = self.new_windows.drain(..).collect();
            for window in new_windows.drain(..) {
                self.register_window(window);
            }

            // Close windows
            self.process_window_close_queue();

            // Show windows
            let windows_to_show = self
                .windows
                .iter()
                .map(|(id, window)| {
                    (
                        *id,
                        window.clone(),
                        self.pending_window_initial_positions.remove(id),
                    )
                })
                .collect::<Vec<_>>();
            for (id, window, initial_position) in windows_to_show {
                self.show_window(ctx, id, window, initial_position);
                self.finalize_viewport_open_probe(id, "Sequence");
            }

            self.render_detached_auxiliary_window_hosts(ctx);

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
#[path = "app/tests.rs"]
mod tests;
