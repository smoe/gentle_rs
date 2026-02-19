use std::{
    collections::HashMap,
    env, fs,
    hash::{Hash, Hasher},
    io::ErrorKind,
    panic::{catch_unwind, AssertUnwindSafe},
    path::{Path, PathBuf},
    process::Command,
    sync::{
        atomic::{AtomicBool, Ordering},
        mpsc, Arc, RwLock,
    },
    time::{Duration, Instant},
};

use crate::{
    about,
    dna_sequence::{self, DNAsequence},
    engine::{
        BlastHitFeatureInput, DisplaySettings, DisplayTarget, Engine, EngineError, ErrorCode,
        GenomeTrackSource, GenomeTrackSubscription, GenomeTrackSyncReport, GentleEngine, OpResult,
        Operation, ProjectState, BIGWIG_TO_BEDGRAPH_ENV_BIN, DEFAULT_BIGWIG_TO_BEDGRAPH_BIN,
    },
    enzymes,
    genomes::{
        GenomeBlastReport, GenomeCatalog, GenomeGeneRecord, GenomeSourcePlan,
        PrepareGenomeProgress, PreparedGenomeInspection, BLASTN_ENV_BIN, DEFAULT_BLASTN_BIN,
        DEFAULT_GENOME_CACHE_DIR, DEFAULT_GENOME_CATALOG_PATH, DEFAULT_HELPER_GENOME_CATALOG_PATH,
        DEFAULT_MAKEBLASTDB_BIN, MAKEBLASTDB_ENV_BIN,
    },
    icons::APP_ICON,
    resource_sync, tf_motifs,
    window::Window,
    TRANSLATIONS,
};
use anyhow::{anyhow, Result};
use eframe::egui::{self, menu, Key, KeyboardShortcut, Modifiers, Pos2, Ui, Vec2, ViewportId};
use egui_commonmark::{CommonMarkCache, CommonMarkViewer};
use pulldown_cmark::{Event, LinkType, Parser, Tag};
use regex::{Regex, RegexBuilder};
use serde::{Deserialize, Serialize};

const GUI_MANUAL_MD: &str = include_str!("../docs/gui.md");
const CLI_MANUAL_MD: &str = include_str!("../docs/cli.md");
const APP_CONFIGURATION_FILE_NAME: &str = ".gentle_gui_settings.json";
static NATIVE_HELP_OPEN_REQUESTED: AtomicBool = AtomicBool::new(false);
static NATIVE_SETTINGS_OPEN_REQUESTED: AtomicBool = AtomicBool::new(false);

pub fn request_open_help_from_native_menu() {
    NATIVE_HELP_OPEN_REQUESTED.store(true, Ordering::SeqCst);
}

pub fn request_open_settings_from_native_menu() {
    NATIVE_SETTINGS_OPEN_REQUESTED.store(true, Ordering::SeqCst);
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
struct PersistedConfiguration {
    rnapkin_executable: String,
    makeblastdb_executable: String,
    blastn_executable: String,
    bigwig_to_bedgraph_executable: String,
    graphics_defaults: DisplaySettings,
}

impl Default for PersistedConfiguration {
    fn default() -> Self {
        Self {
            rnapkin_executable: String::new(),
            makeblastdb_executable: String::new(),
            blastn_executable: String::new(),
            bigwig_to_bedgraph_executable: String::new(),
            graphics_defaults: DisplaySettings::default(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
enum GenomeTrackSourceSelection {
    #[default]
    Auto,
    Bed,
    BigWig,
}

impl GenomeTrackSourceSelection {
    fn label(self) -> &'static str {
        match self {
            Self::Auto => "Auto (from extension)",
            Self::Bed => "BED",
            Self::BigWig => "BigWig",
        }
    }

    fn resolve(self, path: &str) -> GenomeTrackSource {
        match self {
            Self::Auto => GenomeTrackSource::from_path(path),
            Self::Bed => GenomeTrackSource::Bed,
            Self::BigWig => GenomeTrackSource::BigWig,
        }
    }
}

pub struct GENtleApp {
    engine: Arc<RwLock<GentleEngine>>,
    new_windows: Vec<Window>,
    windows: HashMap<ViewportId, Arc<RwLock<Window>>>,
    windows_to_close: Arc<RwLock<Vec<ViewportId>>>,
    viewport_id_counter: usize,
    update_has_run_before: bool,
    show_about_dialog: bool,
    show_help_dialog: bool,
    help_doc: HelpDoc,
    help_markdown_cache: CommonMarkCache,
    help_gui_markdown: String,
    help_cli_markdown: String,
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
    configuration_status: String,
    current_project_path: Option<String>,
    lineage_graph_view: bool,
    lineage_graph_zoom: f32,
    lineage_cache_stamp: u64,
    lineage_cache_valid: bool,
    lineage_rows: Vec<LineageRow>,
    lineage_edges: Vec<(String, String, String)>,
    lineage_op_label_by_id: HashMap<String, String>,
    lineage_containers: Vec<ContainerRow>,
    clean_state_fingerprint: u64,
    dirty_cache_stamp: u64,
    dirty_cache_valid: bool,
    dirty_cache_value: bool,
    last_applied_window_title: String,
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
    genome_track_seq_id: String,
    genome_track_source_selection: GenomeTrackSourceSelection,
    genome_track_path: String,
    genome_track_name: String,
    genome_track_min_score: String,
    genome_track_max_score: String,
    genome_track_clear_existing: bool,
    genome_track_status: String,
    genome_bed_track_subscriptions: Vec<GenomeTrackSubscription>,
    genome_track_autosync_status: String,
    genome_blast_import_track_name: String,
    genome_blast_import_clear_existing: bool,
}

#[derive(Clone, Copy)]
enum ProjectAction {
    New,
    Open,
    Close,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum HelpDoc {
    Gui,
    Cli,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum ConfigurationTab {
    ExternalApplications,
    Graphics,
}

struct GenomePrepareTask {
    started: Instant,
    receiver: mpsc::Receiver<GenomePrepareTaskMessage>,
}

enum GenomePrepareTaskMessage {
    Progress(PrepareGenomeProgress),
    Done(Result<OpResult, EngineError>),
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum GenomeBlastSourceMode {
    Manual,
    ProjectSequence,
    ProjectPool,
}

struct GenomeBlastTask {
    started: Instant,
    receiver: mpsc::Receiver<GenomeBlastTaskMessage>,
}

enum GenomeBlastTaskMessage {
    Progress {
        done_queries: usize,
        total_queries: usize,
        current_query_label: String,
    },
    Done(Result<GenomeBlastBatchResult, String>),
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

#[derive(Clone)]
struct BlastPoolOption {
    container_id: String,
    label: String,
    members: Vec<String>,
}

#[derive(Clone)]
struct LineageRow {
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
}

#[derive(Clone)]
struct ContainerRow {
    container_id: String,
    kind: String,
    member_count: usize,
    representative: String,
    members: Vec<String>,
}

impl Default for GENtleApp {
    fn default() -> Self {
        Self {
            engine: Arc::new(RwLock::new(GentleEngine::new())),
            new_windows: vec![],
            windows: HashMap::new(),
            windows_to_close: Arc::new(RwLock::new(vec![])),
            viewport_id_counter: 0,
            update_has_run_before: false,
            show_about_dialog: false,
            show_help_dialog: false,
            help_doc: HelpDoc::Gui,
            help_markdown_cache: CommonMarkCache::default(),
            help_gui_markdown: GUI_MANUAL_MD.to_string(),
            help_cli_markdown: CLI_MANUAL_MD.to_string(),
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
            configuration_status: String::new(),
            current_project_path: None,
            lineage_graph_view: false,
            lineage_graph_zoom: 1.0,
            lineage_cache_stamp: 0,
            lineage_cache_valid: false,
            lineage_rows: vec![],
            lineage_edges: vec![],
            lineage_op_label_by_id: HashMap::new(),
            lineage_containers: vec![],
            clean_state_fingerprint: 0,
            dirty_cache_stamp: 0,
            dirty_cache_valid: false,
            dirty_cache_value: false,
            last_applied_window_title: String::new(),
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
            genome_track_seq_id: String::new(),
            genome_track_source_selection: GenomeTrackSourceSelection::Auto,
            genome_track_path: String::new(),
            genome_track_name: String::new(),
            genome_track_min_score: String::new(),
            genome_track_max_score: String::new(),
            genome_track_clear_existing: false,
            genome_track_status: String::new(),
            genome_bed_track_subscriptions: vec![],
            genome_track_autosync_status: String::new(),
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

    fn apply_graphics_settings_to_display(source: &DisplaySettings, target: &mut DisplaySettings) {
        target.show_sequence_panel = source.show_sequence_panel;
        target.show_map_panel = source.show_map_panel;
        target.show_cds_features = source.show_cds_features;
        target.show_gene_features = source.show_gene_features;
        target.show_mrna_features = source.show_mrna_features;
        target.show_tfbs = source.show_tfbs;
        target.regulatory_tracks_near_baseline = source.regulatory_tracks_near_baseline;
        target.show_restriction_enzymes = source.show_restriction_enzymes;
        target.show_gc_contents = source.show_gc_contents;
        target.show_open_reading_frames = source.show_open_reading_frames;
        target.show_methylation_sites = source.show_methylation_sites;
    }

    fn apply_configuration_graphics_to_engine_state(&mut self) {
        let mut guard = self.engine.write().expect("Engine lock poisoned");
        let display = &mut guard.state_mut().display;
        Self::apply_graphics_settings_to_display(&self.configuration_graphics, display);
        self.configuration_graphics_dirty = false;
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
        if self.configuration_rnapkin_executable.is_empty() {
            env::remove_var("GENTLE_RNAPKIN_BIN");
        } else {
            env::set_var(
                "GENTLE_RNAPKIN_BIN",
                self.configuration_rnapkin_executable.clone(),
            );
        }
        if self.configuration_makeblastdb_executable.is_empty() {
            env::remove_var(MAKEBLASTDB_ENV_BIN);
        } else {
            env::set_var(
                MAKEBLASTDB_ENV_BIN,
                self.configuration_makeblastdb_executable.clone(),
            );
        }
        if self.configuration_blastn_executable.is_empty() {
            env::remove_var(BLASTN_ENV_BIN);
        } else {
            env::set_var(BLASTN_ENV_BIN, self.configuration_blastn_executable.clone());
        }
        if self.configuration_bigwig_to_bedgraph_executable.is_empty() {
            env::remove_var(BIGWIG_TO_BEDGRAPH_ENV_BIN);
        } else {
            env::set_var(
                BIGWIG_TO_BEDGRAPH_ENV_BIN,
                self.configuration_bigwig_to_bedgraph_executable.clone(),
            );
        }
        self.configuration_graphics = saved.graphics_defaults;
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

    fn refresh_help_docs(&mut self) {
        self.help_gui_markdown = Self::load_help_doc("docs/gui.md", GUI_MANUAL_MD);
        self.help_cli_markdown = Self::load_help_doc("docs/cli.md", CLI_MANUAL_MD);
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
        self.configuration_rnapkin_executable = env::var("GENTLE_RNAPKIN_BIN").unwrap_or_default();
        self.configuration_makeblastdb_executable =
            env::var(MAKEBLASTDB_ENV_BIN).unwrap_or_default();
        self.configuration_blastn_executable = env::var(BLASTN_ENV_BIN).unwrap_or_default();
        self.configuration_bigwig_to_bedgraph_executable =
            env::var(BIGWIG_TO_BEDGRAPH_ENV_BIN).unwrap_or_default();
        if let Ok(engine) = self.engine.read() {
            self.configuration_graphics = engine.state().display.clone();
            self.configuration_graphics_dirty = false;
        }
        self.clear_rnapkin_validation();
        self.clear_blast_validation();
    }

    fn open_configuration_dialog(&mut self) {
        self.sync_configuration_from_runtime();
        self.configuration_tab = ConfigurationTab::ExternalApplications;
        self.show_configuration_dialog = true;
        self.configuration_status.clear();
    }

    fn open_help_doc(&mut self, doc: HelpDoc) {
        self.refresh_help_docs();
        self.help_doc = doc;
        self.show_help_dialog = true;
    }

    fn active_help_title(&self) -> &'static str {
        match self.help_doc {
            HelpDoc::Gui => "GUI Manual",
            HelpDoc::Cli => "CLI Manual",
        }
    }

    fn active_help_markdown(&self) -> &str {
        match self.help_doc {
            HelpDoc::Gui => &self.help_gui_markdown,
            HelpDoc::Cli => &self.help_cli_markdown,
        }
    }

    pub fn new() -> Self {
        let mut app = Self::default();
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

    fn open_sequence_window(&mut self, seq_id: &str) {
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
        if let Ok(dna) = Self::load_dna_from_genbank_file(path) {
            Ok(dna)
        } else if let Ok(dna) = Self::load_dna_from_fasta_file(path) {
            Ok(dna)
        } else {
            Err(anyhow!("Could not load file '{path}'"))
        }
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
        state.metadata.len().hash(&mut hasher);

        let mut metadata_keys: Vec<&String> = state.metadata.keys().collect();
        metadata_keys.sort_unstable();
        for key in metadata_keys {
            key.hash(&mut hasher);
            if let Some(value) = state.metadata.get(key) {
                value.to_string().hash(&mut hasher);
            }
        }

        let display = &state.display;
        display.show_sequence_panel.hash(&mut hasher);
        display.show_map_panel.hash(&mut hasher);
        display.show_cds_features.hash(&mut hasher);
        display.show_gene_features.hash(&mut hasher);
        display.show_mrna_features.hash(&mut hasher);
        display.show_tfbs.hash(&mut hasher);
        display.regulatory_tracks_near_baseline.hash(&mut hasher);
        display.show_restriction_enzymes.hash(&mut hasher);
        display.show_gc_contents.hash(&mut hasher);
        display.show_open_reading_frames.hash(&mut hasher);
        display.show_methylation_sites.hash(&mut hasher);

        hasher.finish()
    }

    fn mark_clean_snapshot(&mut self) {
        self.clean_state_fingerprint = self.current_state_fingerprint();
        self.dirty_cache_stamp = self.current_state_change_stamp();
        self.dirty_cache_valid = true;
        self.dirty_cache_value = false;
    }

    fn is_project_dirty(&mut self) -> bool {
        let stamp = self.current_state_change_stamp();
        if !self.dirty_cache_valid || stamp != self.dirty_cache_stamp {
            self.dirty_cache_value =
                self.current_state_fingerprint() != self.clean_state_fingerprint;
            self.dirty_cache_stamp = stamp;
            self.dirty_cache_valid = true;
        }
        self.dirty_cache_value
    }

    fn reset_to_empty_project(&mut self) {
        self.engine = Arc::new(RwLock::new(GentleEngine::new()));
        self.apply_configuration_graphics_to_engine_state();
        self.current_project_path = None;
        self.last_applied_window_title.clear();
        self.new_windows.clear();
        self.windows.clear();
        self.windows_to_close.write().unwrap().clear();
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
        self.genome_track_status.clear();
        self.genome_track_autosync_status.clear();
        self.load_bed_track_subscriptions_from_state();
        self.mark_clean_snapshot();
    }

    fn prompt_open_project(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("GENtle project", &["json"])
            .pick_file()
        {
            let path = path.display().to_string();
            let _ = self.load_project_from_file(&path);
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
                self.current_project_path = Some(path);
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

    fn open_reference_genome_prepare_dialog(&mut self) {
        self.show_reference_genome_prepare_dialog = true;
    }

    fn open_reference_genome_retrieve_dialog(&mut self) {
        self.show_reference_genome_retrieve_dialog = true;
    }

    fn open_reference_genome_blast_dialog(&mut self) {
        self.show_reference_genome_blast_dialog = true;
    }

    fn open_reference_genome_inspector_dialog(&mut self) {
        self.show_reference_genome_inspector_dialog = true;
    }

    fn open_genome_bed_track_dialog(&mut self) {
        self.load_bed_track_subscriptions_from_state();
        self.show_genome_bed_track_dialog = true;
    }

    fn open_helper_genome_prepare_dialog(&mut self) {
        self.genome_catalog_path = DEFAULT_HELPER_GENOME_CATALOG_PATH.to_string();
        self.genome_cache_dir = "data/helper_genomes".to_string();
        self.invalidate_genome_genes();
        self.show_reference_genome_prepare_dialog = true;
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
        self.show_reference_genome_blast_dialog = true;
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

    fn anchored_sequence_ids_for_tracks(&self) -> Vec<String> {
        self.engine
            .read()
            .unwrap()
            .list_sequences_with_genome_anchor()
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

    fn apply_bed_track_to_sequence(
        &mut self,
        seq_id: &str,
        subscription: &GenomeTrackSubscription,
    ) -> Result<OpResult, EngineError> {
        let op = match subscription.source {
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
        };
        self.engine.write().unwrap().apply(op)
    }

    fn validate_bigwig_converter_available(&self) -> Result<String> {
        let configured = self.configuration_bigwig_to_bedgraph_executable.trim();
        let executable = if configured.is_empty() {
            env::var(BIGWIG_TO_BEDGRAPH_ENV_BIN)
                .ok()
                .filter(|v| !v.trim().is_empty())
                .unwrap_or_else(|| DEFAULT_BIGWIG_TO_BEDGRAPH_BIN.to_string())
        } else {
            configured.to_string()
        };
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
        let result = self.apply_bed_track_to_sequence(&seq_id, &subscription);
        let source_label = subscription.source.label();
        match result {
            Ok(r) => {
                self.genome_track_status = Self::format_op_result_status(
                    &format!("Import {source_label} track: ok"),
                    &r.created_seq_ids,
                    &r.warnings,
                    &r.messages,
                );
            }
            Err(e) => {
                self.genome_track_status =
                    format!("Import {source_label} track failed: {}", e.message);
            }
        }
    }

    fn import_genome_bed_track_for_all_anchored_sequences(&mut self, track_subscription: bool) {
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
                self.load_bed_track_subscriptions_from_state();
            }
            Err(e) => {
                self.genome_track_status = format!("Import track failed: {}", e.message);
            }
        }
    }

    fn apply_tracked_bed_subscription_to_all_anchored(&mut self, index: usize) {
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
            }
            Err(e) => {
                self.genome_track_status =
                    format!("Could not re-apply tracked subscription: {}", e.message);
            }
        }
    }

    fn sync_tracked_bed_tracks_for_new_anchors(&mut self) {
        match self
            .engine
            .write()
            .unwrap()
            .sync_tracked_genome_track_subscriptions(true)
        {
            Ok(report) => {
                if report.applied_imports > 0 || report.failed_imports > 0 {
                    self.genome_track_autosync_status =
                        Self::format_track_sync_status("Auto-sync", &report);
                }
            }
            Err(e) => {
                self.genome_track_autosync_status =
                    format!("Auto-sync failed unexpectedly: {}", e.message);
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
        let catalog_path = self.genome_catalog_path_opt();
        let cache_dir = self.genome_cache_dir_opt();
        let (tx, rx) = mpsc::channel::<GenomePrepareTaskMessage>();
        self.genome_prepare_progress = None;
        self.genome_prepare_status =
            format!("Preparing genome '{genome_id}' in background. You can keep using the UI.");
        self.genome_prepare_task = Some(GenomePrepareTask {
            started: Instant::now(),
            receiver: rx,
        });
        std::thread::spawn(move || {
            let tx_progress = tx.clone();
            let mut last_phase = String::new();
            let mut last_percent_tenths: Option<i64> = None;
            let mut last_bytes_bucket: u64 = 0;
            let mut progress_forwarder = move |p: PrepareGenomeProgress| {
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
                    let _ = tx_progress.send(GenomePrepareTaskMessage::Progress(p));
                }
            };
            let outcome = GentleEngine::prepare_reference_genome_once(
                &genome_id,
                catalog_path.as_deref(),
                cache_dir.as_deref(),
                &mut progress_forwarder,
            )
            .map(|report| OpResult {
                op_id: "background-prepare-genome".to_string(),
                created_seq_ids: vec![],
                changed_seq_ids: vec![],
                warnings: vec![],
                messages: vec![GentleEngine::format_prepare_genome_message(
                    &genome_id,
                    cache_dir.as_deref(),
                    &report,
                )],
            });
            let _ = tx.send(GenomePrepareTaskMessage::Done(outcome));
        });
    }

    fn poll_prepare_reference_genome_task(&mut self, ctx: &egui::Context) {
        if self.genome_prepare_task.is_none() {
            return;
        }
        ctx.request_repaint_after(Duration::from_millis(100));
        let mut done: Option<Result<OpResult, EngineError>> = None;
        if let Some(task) = &self.genome_prepare_task {
            const MAX_MESSAGES_PER_TICK: usize = 128;
            for _ in 0..MAX_MESSAGES_PER_TICK {
                match task.receiver.try_recv() {
                    Ok(GenomePrepareTaskMessage::Progress(progress)) => {
                        self.genome_prepare_progress = Some(progress.clone());
                        self.genome_prepare_status = format!(
                            "Preparing genome '{}': {} ({})",
                            progress.genome_id, progress.phase, progress.item
                        );
                    }
                    Ok(GenomePrepareTaskMessage::Done(result)) => {
                        done = Some(result);
                        break;
                    }
                    Err(mpsc::TryRecvError::Empty) => break,
                    Err(mpsc::TryRecvError::Disconnected) => {
                        done = Some(Err(EngineError {
                            code: ErrorCode::Internal,
                            message: "Genome preparation worker disconnected".to_string(),
                        }));
                        break;
                    }
                }
            }
        }
        if let Some(outcome) = done {
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
                }
                Err(e) => {
                    self.genome_prepare_status =
                        format!("Prepare genome failed after {:.1}s: {}", elapsed, e.message);
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
        self.genome_blast_task = Some(GenomeBlastTask {
            started: Instant::now(),
            receiver: rx,
        });
        std::thread::spawn(move || {
            let mut reports: Vec<GenomeBlastQueryResult> = vec![];
            let mut failed_queries: Vec<String> = vec![];

            for (idx, (label, query)) in blast_queries.into_iter().enumerate() {
                let _ = tx.send(GenomeBlastTaskMessage::Progress {
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
            let _ = tx.send(GenomeBlastTaskMessage::Done(done_result));
        });
    }

    fn poll_reference_genome_blast_task(&mut self, ctx: &egui::Context) {
        if self.genome_blast_task.is_none() {
            return;
        }
        ctx.request_repaint_after(Duration::from_millis(100));
        let mut done: Option<Result<GenomeBlastBatchResult, String>> = None;
        if let Some(task) = &self.genome_blast_task {
            const MAX_MESSAGES_PER_TICK: usize = 128;
            for _ in 0..MAX_MESSAGES_PER_TICK {
                match task.receiver.try_recv() {
                    Ok(GenomeBlastTaskMessage::Progress {
                        done_queries,
                        total_queries,
                        current_query_label,
                    }) => {
                        let fraction = if total_queries == 0 {
                            0.0
                        } else {
                            (done_queries as f32 / total_queries as f32).clamp(0.0, 1.0)
                        };
                        self.genome_blast_progress_fraction = Some(fraction);
                        self.genome_blast_progress_label =
                            format!("{done_queries} / {total_queries} ({current_query_label})");
                    }
                    Ok(GenomeBlastTaskMessage::Done(result)) => {
                        done = Some(result);
                        break;
                    }
                    Err(mpsc::TryRecvError::Empty) => break,
                    Err(mpsc::TryRecvError::Disconnected) => {
                        done = Some(Err("BLAST worker disconnected".to_string()));
                        break;
                    }
                }
            }
        }

        if let Some(outcome) = done {
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
                    } else {
                        self.genome_blast_status = format!(
                            "BLAST finished in {:.1}s: {} query result(s), {} total hit(s), {} failed query(ies): {}",
                            elapsed,
                            self.genome_blast_results.len(),
                            hit_total,
                            batch.failed_queries.len(),
                            batch.failed_queries.join(" | ")
                        );
                    }
                }
                Err(e) => {
                    self.genome_blast_results.clear();
                    self.genome_blast_selected_result = 0;
                    self.genome_blast_status = format!("BLAST failed after {:.1}s: {}", elapsed, e);
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
        ui.label("Download and index a reference genome once.");
        ui.horizontal(|ui| {
            ui.label("catalog");
            ui.text_edit_singleline(&mut self.genome_catalog_path);
            if ui.button("Browse...").clicked() {
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
            if ui.button("Browse...").clicked() {
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
        let selection_changed =
            Self::choose_genome_from_catalog(ui, &mut self.genome_id, &preparable_genomes);
        if selection_changed {
            self.invalidate_genome_genes();
        }
        match self.selected_genome_source_plan() {
            Ok(Some(plan)) => {
                ui.small(format!(
                    "sources: sequence={} | annotation={}",
                    plan.sequence_source_type, plan.annotation_source_type
                ));
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
        ui.horizontal(|ui| {
            if ui
                .add_enabled(
                    !running && selected_preparable,
                    egui::Button::new("Prepare Genome"),
                )
                .clicked()
            {
                self.start_prepare_reference_genome();
            }
            if ui.button("Close").clicked() {
                self.show_reference_genome_prepare_dialog = false;
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
        egui::Window::new("Retrieve Genome Sequence")
            .open(&mut open)
            .collapsible(false)
            .resizable(true)
            .show(ctx, |ui| {
                ui.label("Retrieve region sequences from a prepared genome.");
                ui.horizontal(|ui| {
                    ui.label("catalog");
                    ui.text_edit_singleline(&mut self.genome_catalog_path);
                    if ui.button("Browse...").clicked() {
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
                    if ui.button("Browse...").clicked() {
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
                        ui.small(format!(
                            "sources: sequence={} | annotation={}",
                            plan.sequence_source_type, plan.annotation_source_type
                        ));
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
                        if ui.button("Clear").clicked() {
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
                            if ui.button("All").clicked() {
                                for enabled in self.genome_biotype_filter.values_mut() {
                                    *enabled = true;
                                }
                                self.genome_gene_filter_page = 0;
                            }
                            if ui.button("None").clicked() {
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
                                    .desired_width(90.0)
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
                                    .desired_width(90.0)
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
                            .clicked()
                        {
                            self.extract_reference_genome_gene();
                        }
                        if ui.button("Extract Region").clicked() {
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
            if ui.button("Browse...").clicked() {
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
            if ui.button("Browse...").clicked() {
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
                ui.small(format!(
                    "sources: sequence={} | annotation={}",
                    plan.sequence_source_type, plan.annotation_source_type
                ));
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
                .clicked()
            {
                self.start_reference_genome_blast();
            }
            if ui.button("Clear Results").clicked() {
                self.genome_blast_results.clear();
                self.genome_blast_selected_result = 0;
            }
            if ui.button("Close").clicked() {
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
                let button =
                    ui.add_enabled(can_import, egui::Button::new("Import Hits To Sequence"));
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
                ui.label("Inspect prepared references and installation integrity metadata.");
                ui.horizontal(|ui| {
                    ui.label("catalog");
                    ui.text_edit_singleline(&mut self.genome_catalog_path);
                    if ui.button("Browse...").clicked() {
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
                    if ui.button("Browse...").clicked() {
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
                                                if ui.small_button("Retrieve").clicked() {
                                                    self.genome_id = inspection.genome_id.clone();
                                                    self.invalidate_genome_genes();
                                                    self.open_reference_genome_retrieve_dialog();
                                                }
                                                ui.end_row();
                                            }
                                        });
                                });
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
        let anchored_seq_ids = self.anchored_sequence_ids_for_tracks();
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
            });
            if let Some(anchor) = self.describe_sequence_genome_anchor(&self.genome_track_seq_id) {
                ui.small(format!("selected anchor: {anchor}"));
            }
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
                });
            ui.small(format!("Detected extension: {}", detected_source.label()));
        });

        ui.horizontal(|ui| {
            ui.label("track_path");
            ui.text_edit_singleline(&mut self.genome_track_path);
            if ui.button("Browse...").clicked() {
                if let Some(path) = rfd::FileDialog::new()
                    .add_filter("Signal tracks", &["bed", "gz", "bw", "bigwig"])
                    .pick_file()
                {
                    self.genome_track_path = path.display().to_string();
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
        ui.horizontal(|ui| {
            if ui
                .add_enabled(
                    !anchored_seq_ids.is_empty(),
                    egui::Button::new("Import To Selected"),
                )
                .clicked()
            {
                self.import_genome_bed_track_for_selected_sequence();
            }
            if ui
                .add_enabled(
                    !anchored_seq_ids.is_empty(),
                    egui::Button::new("Import To All Anchored (One-Time)"),
                )
                .clicked()
            {
                self.import_genome_bed_track_for_all_anchored_sequences(false);
            }
            if ui
                .add_enabled(
                    !anchored_seq_ids.is_empty(),
                    egui::Button::new("Import To All Anchored + Track"),
                )
                .clicked()
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
                    for (index, subscription) in
                        self.genome_bed_track_subscriptions.iter().enumerate()
                    {
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
                            if ui.small_button("Apply now").clicked() {
                                apply_now_index = Some(index);
                            }
                            if ui.small_button("Remove").clicked() {
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
            if ui.button("Clear Tracked Files").clicked() {
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
                self.current_project_path = Some(path);
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
            ProjectAction::Close => self.close_project(),
        }
    }

    fn request_project_action(&mut self, action: ProjectAction) {
        if self.is_project_dirty() {
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
        self.current_project_path = Some(path.to_string());
        self.new_windows.clear();
        self.windows.clear();
        self.windows_to_close.write().unwrap().clear();
        self.load_bed_track_subscriptions_from_state();

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

    pub fn render_menu_bar(&mut self, ui: &mut Ui) {
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
                if ui
                    .button("Close Project")
                    .on_hover_text("Close the current project from the workspace")
                    .clicked()
                {
                    self.request_project_action(ProjectAction::Close);
                    ui.close_menu();
                }
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
                    .button("Retrieve Genome Sequence...")
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
                    .on_hover_text("Import BED/BigWig tracks onto anchored sequences")
                    .clicked()
                {
                    self.open_genome_bed_track_dialog();
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
                    .button("Retrieve Genome Sequence...")
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
                    .on_hover_text("Import BED/BigWig tracks onto anchored sequences")
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
        let stamp = self.current_state_change_stamp();
        if self.lineage_cache_valid && self.lineage_cache_stamp == stamp {
            return;
        }

        let (rows, lineage_edges, op_label_by_id, containers) = {
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

            let mut out: Vec<LineageRow> = state
                .lineage
                .nodes
                .values()
                .map(|node| {
                    let parents: Vec<String> = state
                        .lineage
                        .edges
                        .iter()
                        .filter(|e| e.to_node_id == node.node_id)
                        .filter_map(|e| state.lineage.nodes.get(&e.from_node_id))
                        .map(|p| p.seq_id.clone())
                        .collect();
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
            (out, lineage_edges, op_label_by_id, containers)
        };

        self.lineage_rows = rows;
        self.lineage_edges = lineage_edges;
        self.lineage_op_label_by_id = op_label_by_id;
        self.lineage_containers = containers;
        self.lineage_cache_stamp = stamp;
        self.lineage_cache_valid = true;
    }

    fn render_main_lineage(&mut self, ui: &mut Ui) {
        self.refresh_lineage_cache_if_needed();

        ui.heading("Lineage Graph");
        ui.label(format!("Project: {}", self.current_project_name()));
        ui.label("Project-level sequence lineage (branch and merge aware)");
        ui.horizontal(|ui| {
            if ui
                .button("Prepare Reference Genome...")
                .on_hover_text("Download and index reference genomes")
                .clicked()
            {
                self.open_reference_genome_prepare_dialog();
            }
            if ui
                .button("Retrieve Genome Sequence...")
                .on_hover_text("Extract sequence regions from prepared genomes")
                .clicked()
            {
                self.open_reference_genome_retrieve_dialog();
            }
            if ui
                .button("BLAST Genome Sequence...")
                .on_hover_text("Run BLAST against prepared genome indices")
                .clicked()
            {
                self.open_reference_genome_blast_dialog();
            }
        });
        ui.horizontal(|ui| {
            let table = ui
                .selectable_label(!self.lineage_graph_view, "Table")
                .on_hover_text("Show lineage as table");
            if table.clicked() {
                self.lineage_graph_view = false;
            }
            let graph = ui
                .selectable_label(self.lineage_graph_view, "Graph")
                .on_hover_text("Show lineage as node graph");
            if graph.clicked() {
                self.lineage_graph_view = true;
            }
        });
        ui.separator();

        let mut open_seq: Option<String> = None;
        let mut open_pool: Option<(String, Vec<String>)> = None;
        let mut select_candidate_from: Option<String> = None;
        if self.lineage_graph_view {
            let rows = &self.lineage_rows;
            let lineage_edges = &self.lineage_edges;
            let op_label_by_id = &self.lineage_op_label_by_id;
            let mut graph_zoom = self.lineage_graph_zoom.clamp(0.35, 4.0);
            ui.horizontal(|ui| {
                ui.label("Legend:");
                ui.colored_label(egui::Color32::from_rgb(90, 140, 210), "● single sequence");
                ui.colored_label(
                    egui::Color32::from_rgb(180, 120, 70),
                    "◆ pool (n = expected variants)",
                );
                ui.separator();
                if ui.button("−").on_hover_text("Zoom out").clicked() {
                    graph_zoom = (graph_zoom / 1.15).clamp(0.35, 4.0);
                }
                if ui.button("+").on_hover_text("Zoom in").clicked() {
                    graph_zoom = (graph_zoom * 1.15).clamp(0.35, 4.0);
                }
                if ui.button("Reset").on_hover_text("Reset zoom").clicked() {
                    graph_zoom = 1.0;
                }
                ui.add(
                    egui::Slider::new(&mut graph_zoom, 0.35..=4.0)
                        .logarithmic(true)
                        .text("Zoom"),
                );
                ui.label(format!("{:.0}%", graph_zoom * 100.0));
            });
            ui.separator();
            egui::ScrollArea::both()
                .drag_to_scroll(true)
                .show(ui, |ui| {
                    let base_width = (rows.len().max(1) as f32) * 180.0 + 120.0;
                    let base_height = 440.0;
                    let width = base_width * graph_zoom;
                    let height = base_height * graph_zoom;
                    let (resp, painter) =
                        ui.allocate_painter(Vec2::new(width, height), egui::Sense::click());
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
                    let rect = resp.rect;
                    let edge_stroke_width = (1.0 * graph_zoom).clamp(1.0, 2.5);
                    let op_font_size = (10.0 * graph_zoom).clamp(9.0, 16.0);
                    let name_font_size = (12.0 * graph_zoom).clamp(10.0, 18.0);
                    let details_font_size = (10.0 * graph_zoom).clamp(9.0, 15.0);
                    let node_radius = 16.0 * graph_zoom;
                    let mut pos_by_node: HashMap<String, Pos2> = HashMap::new();
                    for (idx, row) in rows.iter().enumerate() {
                        let mut h = std::collections::hash_map::DefaultHasher::new();
                        row.node_id.hash(&mut h);
                        let lane = (h.finish() % 4) as f32;
                        let x = rect.left() + (80.0 + idx as f32 * 170.0) * graph_zoom;
                        let y = rect.top() + (70.0 + lane * 80.0) * graph_zoom;
                        pos_by_node.insert(row.node_id.clone(), Pos2::new(x, y));
                    }
                    let mut used_label_rects: Vec<egui::Rect> = Vec::new();
                    for (from_node, to_node, op_id) in lineage_edges {
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
                        let mid = Pos2::new((from.x + to.x) * 0.5, (from.y + to.y) * 0.5);
                        let op_label = op_label_by_id
                            .get(op_id)
                            .cloned()
                            .unwrap_or_else(|| op_id.clone());
                        let display = op_label;
                        let galley = painter.layout_no_wrap(
                            display.clone(),
                            egui::FontId::proportional(op_font_size),
                            egui::Color32::BLACK,
                        );
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
                                used_label_rects.push(candidate_rect);
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
                        if row.pool_size > 1 {
                            let points = vec![
                                pos + Vec2::new(0.0, -16.0 * graph_zoom),
                                pos + Vec2::new(16.0 * graph_zoom, 0.0),
                                pos + Vec2::new(0.0, 16.0 * graph_zoom),
                                pos + Vec2::new(-16.0 * graph_zoom, 0.0),
                            ];
                            painter.add(egui::Shape::convex_polygon(
                                points,
                                egui::Color32::from_rgb(180, 120, 70),
                                egui::Stroke::new(
                                    edge_stroke_width,
                                    egui::Color32::from_rgb(235, 196, 150),
                                ),
                            ));
                            painter.text(
                                pos + Vec2::new(19.0 * graph_zoom, -14.0 * graph_zoom),
                                egui::Align2::LEFT_TOP,
                                format!("n={}", row.pool_size),
                                egui::FontId::proportional(details_font_size),
                                egui::Color32::YELLOW,
                            );
                        } else {
                            painter.circle_filled(
                                pos,
                                node_radius,
                                egui::Color32::from_rgb(90, 140, 210),
                            );
                        }
                        painter.text(
                            pos + Vec2::new(22.0 * graph_zoom, -4.0 * graph_zoom),
                            egui::Align2::LEFT_BOTTOM,
                            &row.display_name,
                            egui::FontId::proportional(name_font_size),
                            egui::Color32::BLACK,
                        );
                        painter.text(
                            pos + Vec2::new(22.0 * graph_zoom, 10.0 * graph_zoom),
                            egui::Align2::LEFT_TOP,
                            format!("{} ({} bp)", row.seq_id, row.length),
                            egui::FontId::proportional(details_font_size),
                            egui::Color32::BLACK,
                        );
                    }

                    // Interactions: double-click opens sequence, single-click keeps current behavior.
                    if let Some(pointer) = resp.interact_pointer_pos() {
                        let mut hit_row: Option<LineageRow> = None;
                        for row in rows {
                            let Some(pos) = pos_by_node.get(&row.node_id).cloned() else {
                                continue;
                            };
                            let hit = if row.pool_size > 1 {
                                egui::Rect::from_center_size(
                                    pos,
                                    Vec2::new(30.0 * graph_zoom, 24.0 * graph_zoom),
                                )
                                .contains(pointer)
                            } else {
                                pointer.distance(pos) <= 18.0 * graph_zoom
                            };
                            if hit {
                                hit_row = Some(row.clone());
                                break;
                            }
                        }
                        if let Some(row) = hit_row {
                            if resp.double_clicked() {
                                open_seq = Some(row.seq_id.clone());
                            } else if resp.clicked() {
                                if row.pool_size > 1 {
                                    open_pool =
                                        Some((row.seq_id.clone(), row.pool_members.clone()));
                                } else {
                                    open_seq = Some(row.seq_id.clone());
                                }
                            }
                        }
                    }
                });
            self.lineage_graph_zoom = graph_zoom;
        } else {
            egui::ScrollArea::vertical().show(ui, |ui| {
                egui::Grid::new("lineage_grid")
                    .striped(true)
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
                            if ui.button(&row.seq_id).clicked() {
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
                            if ui.button("Select").clicked() {
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
        egui::ScrollArea::vertical()
            .max_height(180.0)
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
                            if c.member_count > 1 {
                                if ui.button("Open Pool").clicked() {
                                    open_pool = Some((c.representative.clone(), c.members.clone()));
                                }
                            } else if !c.representative.is_empty() {
                                if ui.button("Open Seq").clicked() {
                                    open_seq = Some(c.representative.clone());
                                }
                            } else {
                                ui.label("-");
                            }
                            ui.end_row();
                        }
                    });
            });

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

        if self.configuration_rnapkin_executable.is_empty() {
            env::remove_var("GENTLE_RNAPKIN_BIN");
        } else {
            env::set_var(
                "GENTLE_RNAPKIN_BIN",
                self.configuration_rnapkin_executable.clone(),
            );
        }

        if self.configuration_makeblastdb_executable.is_empty() {
            env::remove_var(MAKEBLASTDB_ENV_BIN);
        } else {
            env::set_var(
                MAKEBLASTDB_ENV_BIN,
                self.configuration_makeblastdb_executable.clone(),
            );
        }

        if self.configuration_blastn_executable.is_empty() {
            env::remove_var(BLASTN_ENV_BIN);
        } else {
            env::set_var(BLASTN_ENV_BIN, self.configuration_blastn_executable.clone());
        }

        if self.configuration_bigwig_to_bedgraph_executable.is_empty() {
            env::remove_var(BIGWIG_TO_BEDGRAPH_ENV_BIN);
        } else {
            env::set_var(
                BIGWIG_TO_BEDGRAPH_ENV_BIN,
                self.configuration_bigwig_to_bedgraph_executable.clone(),
            );
        }

        let rnapkin_status = env::var("GENTLE_RNAPKIN_BIN")
            .ok()
            .filter(|v| !v.trim().is_empty())
            .unwrap_or_else(|| "PATH lookup: rnapkin".to_string());
        let makeblastdb_status = env::var(MAKEBLASTDB_ENV_BIN)
            .ok()
            .filter(|v| !v.trim().is_empty())
            .unwrap_or_else(|| format!("PATH lookup: {DEFAULT_MAKEBLASTDB_BIN}"));
        let blastn_status = env::var(BLASTN_ENV_BIN)
            .ok()
            .filter(|v| !v.trim().is_empty())
            .unwrap_or_else(|| format!("PATH lookup: {DEFAULT_BLASTN_BIN}"));
        let bigwig_status = env::var(BIGWIG_TO_BEDGRAPH_ENV_BIN)
            .ok()
            .filter(|v| !v.trim().is_empty())
            .unwrap_or_else(|| format!("PATH lookup: {DEFAULT_BIGWIG_TO_BEDGRAPH_BIN}"));
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
        self.configuration_graphics.show_cds_features = defaults.show_cds_features;
        self.configuration_graphics.show_gene_features = defaults.show_gene_features;
        self.configuration_graphics.show_mrna_features = defaults.show_mrna_features;
        self.configuration_graphics.show_tfbs = defaults.show_tfbs;
        self.configuration_graphics.regulatory_tracks_near_baseline =
            defaults.regulatory_tracks_near_baseline;
        self.configuration_graphics.show_restriction_enzymes = defaults.show_restriction_enzymes;
        self.configuration_graphics.show_gc_contents = defaults.show_gc_contents;
        self.configuration_graphics.show_open_reading_frames = defaults.show_open_reading_frames;
        self.configuration_graphics.show_methylation_sites = defaults.show_methylation_sites;
        self.configuration_graphics_dirty = true;
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
        let active_rnapkin = env::var("GENTLE_RNAPKIN_BIN")
            .ok()
            .filter(|v| !v.trim().is_empty())
            .unwrap_or_else(|| "PATH lookup: rnapkin".to_string());
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
        let active_makeblastdb = env::var(MAKEBLASTDB_ENV_BIN)
            .ok()
            .filter(|v| !v.trim().is_empty())
            .unwrap_or_else(|| format!("PATH lookup: {DEFAULT_MAKEBLASTDB_BIN}"));
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
        let active_blastn = env::var(BLASTN_ENV_BIN)
            .ok()
            .filter(|v| !v.trim().is_empty())
            .unwrap_or_else(|| format!("PATH lookup: {DEFAULT_BLASTN_BIN}"));
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
        let active_bigwig = env::var(BIGWIG_TO_BEDGRAPH_ENV_BIN)
            .ok()
            .filter(|v| !v.trim().is_empty())
            .unwrap_or_else(|| format!("PATH lookup: {DEFAULT_BIGWIG_TO_BEDGRAPH_BIN}"));
        ui.monospace(format!("Active bigWigToBedGraph: {active_bigwig}"));

        ui.horizontal(|ui| {
            if ui.button("Use PATH").clicked() {
                self.configuration_rnapkin_executable.clear();
                self.clear_rnapkin_validation();
            }
            if ui.button("Use PATH (BLAST)").clicked() {
                self.configuration_makeblastdb_executable.clear();
                self.configuration_blastn_executable.clear();
                self.clear_blast_validation();
            }
            if ui.button("Use PATH (BigWig)").clicked() {
                self.configuration_bigwig_to_bedgraph_executable.clear();
            }
            if ui.button("Validate rnapkin").clicked() {
                self.validate_rnapkin_executable();
            }
            if ui.button("Validate BLAST tools").clicked() {
                self.validate_blast_executables();
            }
            if ui.button("Apply External Settings").clicked() {
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

        if changed {
            self.configuration_graphics_dirty = true;
        }

        ui.horizontal(|ui| {
            if ui.button("Reset Defaults").clicked() {
                self.reset_configuration_graphics_to_defaults();
            }
            if ui
                .add_enabled(
                    self.configuration_graphics_dirty,
                    egui::Button::new("Apply Graphics Settings"),
                )
                .clicked()
            {
                self.apply_configuration_graphics();
            }
            if ui.button("Apply + Refresh Open Windows").clicked() {
                self.apply_configuration_graphics();
                let refreshed = self.refresh_open_sequence_windows(ui.ctx());
                self.configuration_status.push_str(&format!(
                    " | refreshed {} open sequence window(s)",
                    refreshed
                ));
            }
            if ui.button("Reload Current").clicked() {
                self.sync_configuration_from_runtime();
                self.configuration_status =
                    "Reloaded graphics settings from current project".to_string();
            }
        });
    }

    fn render_configuration_contents(&mut self, ui: &mut Ui) {
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
            if ui.button("Close").clicked() {
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
                    if ui.button("Save").clicked() {
                        if self.save_current_project() {
                            if let Some(action) = self.pending_project_action.take() {
                                self.execute_project_action(action);
                            }
                        }
                    }
                    if ui.button("Don't Save").clicked() {
                        if let Some(action) = self.pending_project_action.take() {
                            self.execute_project_action(action);
                        }
                    }
                    if ui.button("Cancel").clicked() {
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
        ui.horizontal(|ui| {
            if ui
                .selectable_label(self.help_doc == HelpDoc::Gui, "GUI Manual")
                .clicked()
            {
                self.help_doc = HelpDoc::Gui;
            }
            if ui
                .selectable_label(self.help_doc == HelpDoc::Cli, "CLI Manual")
                .clicked()
            {
                self.help_doc = HelpDoc::Cli;
            }
            ui.separator();
            if ui.button("Reload").clicked() {
                self.refresh_help_docs();
                self.help_markdown_cache = CommonMarkCache::default();
            }
            if ui.button("Close").clicked() {
                self.show_help_dialog = false;
            }
        });
        ui.separator();
        egui::ScrollArea::vertical()
            .auto_shrink([false, false])
            .show(ui, |ui| {
                let markdown = self.active_help_markdown().to_string();
                CommonMarkViewer::new().max_image_width(Some(1200)).show(
                    ui,
                    &mut self.help_markdown_cache,
                    &markdown,
                );
            });
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
            Operation::RenderSequenceSvg { seq_id, mode, path } => format!(
                "Render sequence SVG: seq_id={seq_id}, mode={mode:?}, path={path}"
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
            } => format!(
                "Render pool gel SVG: inputs={}, path={}, ladders={}",
                inputs.join(", "),
                path,
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
            } => format!(
                "Prepare genome: genome_id={}, catalog_path={}, cache_dir={}",
                genome_id,
                catalog_path.clone().unwrap_or_else(|| "-".to_string()),
                cache_dir.clone().unwrap_or_else(|| "-".to_string())
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
            self.consume_native_help_request();
            self.consume_native_settings_request();
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
            let save_project = KeyboardShortcut::new(Modifiers::COMMAND, Key::S);
            let close_project =
                KeyboardShortcut::new(Modifiers::COMMAND | Modifiers::SHIFT, Key::W);
            let open_configuration = KeyboardShortcut::new(Modifiers::COMMAND, Key::Comma);
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
            if ctx.input_mut(|i| i.consume_shortcut(&save_project)) {
                let _ = self.save_current_project();
            }
            if ctx.input_mut(|i| i.consume_shortcut(&close_project)) {
                self.request_project_action(ProjectAction::Close);
            }
            if ctx.input_mut(|i| i.consume_shortcut(&open_configuration)) {
                self.open_configuration_dialog();
            }

            self.poll_prepare_reference_genome_task(ctx);
            self.poll_reference_genome_blast_task(ctx);
            self.sync_tracked_bed_tracks_for_new_anchors();

            // Show menu bar
            egui::TopBottomPanel::top("top").show(ctx, |ui| {
                self.render_menu_bar(ui);
            });

            // Show main window
            egui::CentralPanel::default().show(ctx, |ui| {
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
            self.render_configuration_dialog(ctx);
            self.render_help_dialog(ctx);
            self.render_about_dialog(ctx);
            self.render_unsaved_changes_dialog(ctx);

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
        }));
        if update_result.is_err() {
            eprintln!("E GENtleApp: recovered from panic in app update");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::GENtleApp;
    use std::fs;
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
}
