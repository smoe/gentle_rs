use std::{
    collections::{HashMap, HashSet},
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
        DisplaySettings, DisplayTarget, Engine, EngineError, ErrorCode, GentleEngine, OpResult,
        Operation, ProjectState,
    },
    enzymes,
    genomes::{
        GenomeCatalog, GenomeGeneRecord, GenomeSourcePlan, PrepareGenomeProgress,
        PreparedGenomeInspection, DEFAULT_GENOME_CACHE_DIR, DEFAULT_GENOME_CATALOG_PATH,
        DEFAULT_HELPER_GENOME_CATALOG_PATH,
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
const GENOME_PREPARE_STEP_LABELS: [&str; 5] = [
    "Resolve/download sequence",
    "Resolve/download annotation",
    "Build FASTA index",
    "Build gene index",
    "Finalize installation",
];

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
    graphics_defaults: DisplaySettings,
}

impl Default for PersistedConfiguration {
    fn default() -> Self {
        Self {
            rnapkin_executable: String::new(),
            graphics_defaults: DisplaySettings::default(),
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
    configuration_rnapkin_validation_ok: Option<bool>,
    configuration_rnapkin_validation_message: String,
    configuration_graphics: DisplaySettings,
    configuration_graphics_dirty: bool,
    configuration_status: String,
    current_project_path: Option<String>,
    lineage_graph_view: bool,
    clean_state_fingerprint: u64,
    pending_project_action: Option<ProjectAction>,
    show_reference_genome_prepare_dialog: bool,
    show_reference_genome_retrieve_dialog: bool,
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
    show_genome_bed_track_dialog: bool,
    genome_track_seq_id: String,
    genome_track_path: String,
    genome_track_name: String,
    genome_track_min_score: String,
    genome_track_max_score: String,
    genome_track_clear_existing: bool,
    genome_track_status: String,
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
            configuration_rnapkin_validation_ok: None,
            configuration_rnapkin_validation_message: String::new(),
            configuration_graphics: DisplaySettings::default(),
            configuration_graphics_dirty: false,
            configuration_status: String::new(),
            current_project_path: None,
            lineage_graph_view: false,
            clean_state_fingerprint: 0,
            pending_project_action: None,
            show_reference_genome_prepare_dialog: false,
            show_reference_genome_retrieve_dialog: false,
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
            show_genome_bed_track_dialog: false,
            genome_track_seq_id: String::new(),
            genome_track_path: String::new(),
            genome_track_name: String::new(),
            genome_track_min_score: String::new(),
            genome_track_max_score: String::new(),
            genome_track_clear_existing: false,
            genome_track_status: String::new(),
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

    fn configuration_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("GENtle Configuration Viewport")
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
        if self.configuration_rnapkin_executable.is_empty() {
            env::remove_var("GENTLE_RNAPKIN_BIN");
        } else {
            env::set_var(
                "GENTLE_RNAPKIN_BIN",
                self.configuration_rnapkin_executable.clone(),
            );
        }
        self.configuration_graphics = saved.graphics_defaults;
        self.configuration_graphics_dirty = false;
        self.configuration_rnapkin_validation_ok = None;
        self.configuration_rnapkin_validation_message.clear();

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

    fn clear_rnapkin_validation(&mut self) {
        self.configuration_rnapkin_validation_ok = None;
        self.configuration_rnapkin_validation_message.clear();
    }

    fn validate_rnapkin_executable(&mut self) {
        let executable = self.resolved_rnapkin_executable();
        match Command::new(&executable).arg("--version").output() {
            Ok(output) => {
                let stdout = String::from_utf8_lossy(&output.stdout).to_string();
                let stderr = String::from_utf8_lossy(&output.stderr).to_string();
                let first_non_empty = stdout
                    .lines()
                    .chain(stderr.lines())
                    .find(|line| !line.trim().is_empty())
                    .map(|s| s.trim().to_string())
                    .unwrap_or_else(|| "No version output".to_string());
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

    fn sync_configuration_from_runtime(&mut self) {
        self.configuration_rnapkin_executable = env::var("GENTLE_RNAPKIN_BIN").unwrap_or_default();
        if let Ok(engine) = self.engine.read() {
            self.configuration_graphics = engine.state().display.clone();
            self.configuration_graphics_dirty = false;
        }
        self.clear_rnapkin_validation();
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

    fn mark_clean_snapshot(&mut self) {
        self.clean_state_fingerprint = self.current_state_fingerprint();
    }

    fn is_project_dirty(&self) -> bool {
        self.current_state_fingerprint() != self.clean_state_fingerprint
    }

    fn reset_to_empty_project(&mut self) {
        self.engine = Arc::new(RwLock::new(GentleEngine::new()));
        self.apply_configuration_graphics_to_engine_state();
        self.current_project_path = None;
        self.new_windows.clear();
        self.windows.clear();
        self.windows_to_close.write().unwrap().clear();
        self.viewport_id_counter = 0;
        self.pending_project_action = None;
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

    fn open_reference_genome_inspector_dialog(&mut self) {
        self.show_reference_genome_inspector_dialog = true;
    }

    fn open_genome_bed_track_dialog(&mut self) {
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

    fn import_genome_bed_track(&mut self) {
        let seq_id = self.genome_track_seq_id.trim().to_string();
        if seq_id.is_empty() {
            self.genome_track_status =
                "Select a genome-anchored sequence before importing BED".to_string();
            return;
        }
        let path = self.genome_track_path.trim().to_string();
        if path.is_empty() {
            self.genome_track_status = "BED path cannot be empty".to_string();
            return;
        }
        let min_score =
            match Self::parse_optional_score_field(&self.genome_track_min_score, "min score") {
                Ok(v) => v,
                Err(e) => {
                    self.genome_track_status = format!("Import BED failed: {e}");
                    return;
                }
            };
        let max_score =
            match Self::parse_optional_score_field(&self.genome_track_max_score, "max score") {
                Ok(v) => v,
                Err(e) => {
                    self.genome_track_status = format!("Import BED failed: {e}");
                    return;
                }
            };
        if min_score
            .zip(max_score)
            .map(|(min, max)| min > max)
            .unwrap_or(false)
        {
            self.genome_track_status =
                "Import BED failed: min score must be <= max score".to_string();
            return;
        }
        let track_name = {
            let trimmed = self.genome_track_name.trim();
            if trimmed.is_empty() {
                None
            } else {
                Some(trimmed.to_string())
            }
        };

        let op = Operation::ImportGenomeBedTrack {
            seq_id,
            path,
            track_name,
            min_score,
            max_score,
            clear_existing: Some(self.genome_track_clear_existing),
        };
        let result = { self.engine.write().unwrap().apply(op) };
        match result {
            Ok(r) => {
                self.genome_track_status = Self::format_op_result_status(
                    "Import BED track: ok",
                    &r.created_seq_ids,
                    &r.warnings,
                    &r.messages,
                );
            }
            Err(e) => {
                self.genome_track_status = format!("Import BED track failed: {}", e.message);
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

    fn render_genome_bed_track_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_genome_bed_track_dialog {
            return;
        }
        let anchored_seq_ids = self.anchored_sequence_ids_for_tracks();
        if !anchored_seq_ids.is_empty() && !anchored_seq_ids.contains(&self.genome_track_seq_id) {
            self.genome_track_seq_id = anchored_seq_ids[0].clone();
        }

        let mut open = self.show_genome_bed_track_dialog;
        egui::Window::new("Import BED Track")
            .open(&mut open)
            .collapsible(false)
            .resizable(true)
            .show(ctx, |ui| {
                ui.label("Map BED intervals onto a genome-anchored sequence.");
                ui.small("Supports .bed and .bed.gz input files.");

                if anchored_seq_ids.is_empty() {
                    ui.colored_label(
                        egui::Color32::from_rgb(190, 70, 70),
                        "No genome-anchored sequence is available. Extract a genome region or gene first.",
                    );
                } else {
                    ui.horizontal(|ui| {
                        ui.label("sequence");
                        egui::ComboBox::from_id_salt("genome_track_seq_id_combo")
                            .selected_text(
                                if self.genome_track_seq_id.trim().is_empty() {
                                    "<select sequence>"
                                } else {
                                    self.genome_track_seq_id.as_str()
                                },
                            )
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
                }

                ui.horizontal(|ui| {
                    ui.label("bed_path");
                    ui.text_edit_singleline(&mut self.genome_track_path);
                    if ui.button("Browse...").clicked() {
                        if let Some(path) = rfd::FileDialog::new()
                            .add_filter("BED", &["bed", "gz"])
                            .pick_file()
                        {
                            self.genome_track_path = path.display().to_string();
                        }
                    }
                });
                ui.horizontal(|ui| {
                    ui.label("track_name");
                    ui.text_edit_singleline(&mut self.genome_track_name);
                    ui.small("(optional)");
                });
                ui.horizontal(|ui| {
                    ui.label("min_score");
                    ui.add(
                        egui::TextEdit::singleline(&mut self.genome_track_min_score)
                            .desired_width(90.0),
                    );
                    ui.label("max_score");
                    ui.add(
                        egui::TextEdit::singleline(&mut self.genome_track_max_score)
                            .desired_width(90.0),
                    );
                });
                ui.checkbox(
                    &mut self.genome_track_clear_existing,
                    "Clear existing imported BED track features first",
                );
                ui.horizontal(|ui| {
                    if ui
                        .add_enabled(
                            !anchored_seq_ids.is_empty(),
                            egui::Button::new("Import BED Track"),
                        )
                        .clicked()
                    {
                        self.import_genome_bed_track();
                    }
                });
                if !self.genome_track_status.is_empty() {
                    ui.separator();
                    ui.monospace(&self.genome_track_status);
                }
            });
        self.show_genome_bed_track_dialog = open;
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
                if ui.button("New Project").clicked() {
                    self.request_project_action(ProjectAction::New);
                    ui.close_menu();
                }
                if ui.button("Open Project...").clicked() {
                    self.request_project_action(ProjectAction::Open);
                    ui.close_menu();
                }
                if ui.button("Close Project").clicked() {
                    self.request_project_action(ProjectAction::Close);
                    ui.close_menu();
                }
                if ui.button("Open Sequence...").clicked() {
                    self.prompt_open_sequence();
                    ui.close_menu();
                }
                ui.separator();
                if ui.button("Configuration...").clicked() {
                    self.open_configuration_dialog();
                    ui.close_menu();
                }
                ui.separator();
                if ui.button("Prepare Reference Genome...").clicked() {
                    self.open_reference_genome_prepare_dialog();
                    ui.close_menu();
                }
                if ui.button("Prepared References...").clicked() {
                    self.open_reference_genome_inspector_dialog();
                    ui.close_menu();
                }
                if ui.button("Retrieve Genome Sequence...").clicked() {
                    self.open_reference_genome_retrieve_dialog();
                    ui.close_menu();
                }
                if ui.button("Import BED Track...").clicked() {
                    self.open_genome_bed_track_dialog();
                    ui.close_menu();
                }
                if ui.button("Prepare Helper Genome...").clicked() {
                    self.open_helper_genome_prepare_dialog();
                    ui.close_menu();
                }
                if ui.button("Retrieve Helper Sequence...").clicked() {
                    self.open_helper_genome_retrieve_dialog();
                    ui.close_menu();
                }
                if ui.button("Import REBASE Data...").clicked() {
                    self.prompt_import_rebase_resource();
                    ui.close_menu();
                }
                if ui.button("Import JASPAR Data...").clicked() {
                    self.prompt_import_jaspar_resource();
                    ui.close_menu();
                }
                if ui.button("Save Project...").clicked() {
                    self.prompt_save_project();
                    ui.close_menu();
                }
                if ui.button("Export DALG SVG...").clicked() {
                    self.prompt_export_lineage_svg();
                    ui.close_menu();
                }
            });
            ui.menu_button("Settings", |ui| {
                if ui.button("Configuration...").clicked() {
                    self.open_configuration_dialog();
                    ui.close_menu();
                }
            });
            ui.menu_button("Genome", |ui| {
                if ui.button("Prepare Reference Genome...").clicked() {
                    self.open_reference_genome_prepare_dialog();
                    ui.close_menu();
                }
                if ui.button("Prepared References...").clicked() {
                    self.open_reference_genome_inspector_dialog();
                    ui.close_menu();
                }
                if ui.button("Retrieve Genome Sequence...").clicked() {
                    self.open_reference_genome_retrieve_dialog();
                    ui.close_menu();
                }
                if ui.button("Import BED Track...").clicked() {
                    self.open_genome_bed_track_dialog();
                    ui.close_menu();
                }
                ui.separator();
                if ui.button("Prepare Helper Genome...").clicked() {
                    self.open_helper_genome_prepare_dialog();
                    ui.close_menu();
                }
                if ui.button("Retrieve Helper Sequence...").clicked() {
                    self.open_helper_genome_retrieve_dialog();
                    ui.close_menu();
                }
            });
            ui.menu_button("Help", |ui| {
                if ui.button("GUI Manual").clicked() {
                    self.open_help_doc(HelpDoc::Gui);
                    ui.close_menu();
                }
                if ui.button("CLI Manual").clicked() {
                    self.open_help_doc(HelpDoc::Cli);
                    ui.close_menu();
                }
                ui.separator();
                if ui.button("About GENtle").clicked() {
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

    fn render_main_lineage(&mut self, ui: &mut Ui) {
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
            ui.horizontal(|ui| {
                ui.label("Legend:");
                ui.colored_label(egui::Color32::from_rgb(90, 140, 210), "● single sequence");
                ui.colored_label(
                    egui::Color32::from_rgb(180, 120, 70),
                    "◆ pool (n = expected variants)",
                );
            });
            ui.separator();
            egui::ScrollArea::both().show(ui, |ui| {
                let width = (rows.len().max(1) as f32) * 180.0 + 120.0;
                let height = 440.0;
                let (resp, painter) =
                    ui.allocate_painter(Vec2::new(width, height), egui::Sense::hover());
                let rect = resp.rect;
                let mut pos_by_node: HashMap<String, Pos2> = HashMap::new();
                for (idx, row) in rows.iter().enumerate() {
                    let mut h = std::collections::hash_map::DefaultHasher::new();
                    row.node_id.hash(&mut h);
                    let lane = (h.finish() % 4) as f32;
                    let x = rect.left() + 80.0 + idx as f32 * 170.0;
                    let y = rect.top() + 70.0 + lane * 80.0;
                    pos_by_node.insert(row.node_id.clone(), Pos2::new(x, y));
                }
                let mut used_label_rects: Vec<egui::Rect> = Vec::new();
                for (from_node, to_node, op_id) in &lineage_edges {
                    let Some(from) = pos_by_node.get(from_node) else {
                        continue;
                    };
                    let Some(to) = pos_by_node.get(to_node) else {
                        continue;
                    };
                    painter.line_segment([*from, *to], egui::Stroke::new(1.0, egui::Color32::GRAY));
                    let mid = Pos2::new((from.x + to.x) * 0.5, (from.y + to.y) * 0.5);
                    let op_label = op_label_by_id
                        .get(op_id)
                        .cloned()
                        .unwrap_or_else(|| op_id.clone());
                    let display = op_label;
                    let galley = painter.layout_no_wrap(
                        display.clone(),
                        egui::FontId::proportional(10.0),
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
                        let candidate_center = mid + perp * (12.0 + step * 12.0) * sign;
                        let candidate_rect = egui::Rect::from_center_size(
                            candidate_center,
                            Vec2::new(galley.size().x + 10.0, galley.size().y + 6.0),
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
                        let fallback_center = mid + perp * 12.0;
                        let rect = egui::Rect::from_center_size(
                            fallback_center,
                            Vec2::new(galley.size().x + 10.0, galley.size().y + 6.0),
                        );
                        (fallback_center, rect)
                    });
                    used_label_rects.push(bg_rect);
                    painter.rect_filled(
                        bg_rect,
                        3.0,
                        egui::Color32::from_rgba_premultiplied(245, 245, 245, 235),
                    );
                    painter.text(
                        label_center,
                        egui::Align2::CENTER_CENTER,
                        display,
                        egui::FontId::proportional(10.0),
                        egui::Color32::BLACK,
                    );
                }
                for row in &rows {
                    let Some(pos) = pos_by_node.get(&row.node_id).cloned() else {
                        continue;
                    };
                    if row.pool_size > 1 {
                        let points = vec![
                            pos + Vec2::new(0.0, -16.0),
                            pos + Vec2::new(16.0, 0.0),
                            pos + Vec2::new(0.0, 16.0),
                            pos + Vec2::new(-16.0, 0.0),
                        ];
                        painter.add(egui::Shape::convex_polygon(
                            points,
                            egui::Color32::from_rgb(180, 120, 70),
                            egui::Stroke::new(1.0, egui::Color32::from_rgb(235, 196, 150)),
                        ));
                        painter.text(
                            pos + Vec2::new(19.0, -14.0),
                            egui::Align2::LEFT_TOP,
                            format!("n={}", row.pool_size),
                            egui::FontId::proportional(10.0),
                            egui::Color32::YELLOW,
                        );
                    } else {
                        painter.circle_filled(pos, 16.0, egui::Color32::from_rgb(90, 140, 210));
                    }
                    painter.text(
                        pos + Vec2::new(22.0, -4.0),
                        egui::Align2::LEFT_BOTTOM,
                        &row.display_name,
                        egui::FontId::proportional(12.0),
                        egui::Color32::BLACK,
                    );
                    painter.text(
                        pos + Vec2::new(22.0, 10.0),
                        egui::Align2::LEFT_TOP,
                        format!("{} ({} bp)", row.seq_id, row.length),
                        egui::FontId::proportional(10.0),
                        egui::Color32::BLACK,
                    );
                }

                // Interactions: single-click opens sequence, double-click selects candidate.
                if let Some(pointer) = resp.interact_pointer_pos() {
                    let mut hit_row: Option<LineageRow> = None;
                    for row in &rows {
                        let Some(pos) = pos_by_node.get(&row.node_id).cloned() else {
                            continue;
                        };
                        let hit = if row.pool_size > 1 {
                            egui::Rect::from_center_size(pos, Vec2::new(30.0, 24.0))
                                .contains(pointer)
                        } else {
                            pointer.distance(pos) <= 18.0
                        };
                        if hit {
                            hit_row = Some(row.clone());
                            break;
                        }
                    }
                    if let Some(row) = hit_row {
                        if resp.double_clicked() {
                            select_candidate_from = Some(row.seq_id.clone());
                        } else if resp.clicked() {
                            if row.pool_size > 1 {
                                open_pool = Some((row.seq_id.clone(), row.pool_members.clone()));
                            } else {
                                open_seq = Some(row.seq_id.clone());
                            }
                        }
                    }
                }
            });
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

                        for row in &rows {
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
                        for c in &containers {
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
        let value = self.configuration_rnapkin_executable.trim().to_string();
        if value.is_empty() {
            env::remove_var("GENTLE_RNAPKIN_BIN");
            self.configuration_status =
                "External app settings applied (rnapkin override cleared; PATH lookup active)"
                    .to_string();
        } else {
            env::set_var("GENTLE_RNAPKIN_BIN", &value);
            self.configuration_status = format!(
                "External app settings applied (rnapkin override: {})",
                value
            );
        }
        self.validate_rnapkin_executable();
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
        let edit_response = ui.add(
            egui::TextEdit::singleline(&mut self.configuration_rnapkin_executable)
                .hint_text("Leave empty to use PATH lookup for 'rnapkin'"),
        );
        if edit_response.changed() {
            self.clear_rnapkin_validation();
        }
        let active = env::var("GENTLE_RNAPKIN_BIN")
            .ok()
            .filter(|v| !v.trim().is_empty())
            .unwrap_or_else(|| "PATH lookup: rnapkin".to_string());
        ui.monospace(format!("Active resolution: {active}"));
        ui.horizontal(|ui| {
            if ui.button("Use PATH").clicked() {
                self.configuration_rnapkin_executable.clear();
                self.clear_rnapkin_validation();
            }
            if ui.button("Validate rnapkin").clicked() {
                self.validate_rnapkin_executable();
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
            let dirty_marker = if self.is_project_dirty() { " *" } else { "" };
            ctx.send_viewport_cmd(egui::ViewportCommand::Title(format!(
                "GENtle - {}{} (v{})",
                self.current_project_name(),
                dirty_marker,
                about::GENTLE_DISPLAY_VERSION
            )));

            let open_project = KeyboardShortcut::new(Modifiers::COMMAND, Key::O);
            let new_project = KeyboardShortcut::new(Modifiers::COMMAND, Key::N);
            let open_sequence =
                KeyboardShortcut::new(Modifiers::COMMAND | Modifiers::SHIFT, Key::O);
            let open_retrieve_genome =
                KeyboardShortcut::new(Modifiers::COMMAND | Modifiers::SHIFT, Key::G);
            let open_prepare_genome =
                KeyboardShortcut::new(Modifiers::COMMAND | Modifiers::SHIFT, Key::P);
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

            // Show menu bar
            egui::TopBottomPanel::top("top").show(ctx, |ui| {
                self.render_menu_bar(ui);
            });

            // Show main window
            egui::CentralPanel::default().show(ctx, |ui| {
                self.render_main_lineage(ui);
                if self.is_project_dirty() {
                    ui.label("Status: unsaved changes");
                } else {
                    ui.label("Status: saved");
                }
            });
            self.render_reference_genome_prepare_dialog(ctx);
            self.render_reference_genome_retrieve_dialog(ctx);
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
