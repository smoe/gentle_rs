//! Curated workflow example payloads and templates.

use crate::engine::{Engine, GentleEngine, Operation, OperationProgress, ProjectState, Workflow};
use serde::{Deserialize, Serialize};
use sha1::{Digest, Sha1};
use std::{
    collections::{BTreeMap, HashMap, HashSet},
    env, fs,
    io::Read,
    path::{Path, PathBuf},
};
use tempfile::TempDir;

pub const WORKFLOW_EXAMPLE_SCHEMA: &str = "gentle.workflow_example.v1";
pub const DEFAULT_WORKFLOW_EXAMPLE_DIR: &str = "docs/examples/workflows";
pub const DEFAULT_WORKFLOW_SNIPPET_DIR: &str = "docs/examples/generated";
pub const ONLINE_EXAMPLE_TEST_ENV: &str = "GENTLE_TEST_ONLINE";
pub const SKIP_REMOTE_TESTS_ENV: &str = "GENTLE_SKIP_REMOTE_TESTS";
pub const TUTORIAL_CATALOG_SCHEMA: &str = "gentle.tutorial_catalog.v1";
pub const TUTORIAL_CATALOG_META_SCHEMA: &str = "gentle.tutorial_catalog_meta.v1";
pub const TUTORIAL_SOURCE_SCHEMA: &str = "gentle.tutorial_source.v2";
pub const TUTORIAL_MANIFEST_SCHEMA: &str = "gentle.tutorial_manifest.v1";
pub const TUTORIAL_GENERATION_REPORT_SCHEMA: &str = "gentle.tutorial_generation_report.v1";
pub const DEFAULT_TUTORIAL_CATALOG_PATH: &str = "docs/tutorial/catalog.json";
pub const DEFAULT_TUTORIAL_CATALOG_META_PATH: &str = "docs/tutorial/sources/catalog_meta.json";
pub const DEFAULT_TUTORIAL_SOURCE_DIR: &str = "docs/tutorial/sources";
pub const DEFAULT_TUTORIAL_MANIFEST_PATH: &str = "docs/tutorial/manifest.json";
pub const DEFAULT_TUTORIAL_OUTPUT_DIR: &str = "docs/tutorial/generated";

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TutorialCatalogGeneratedRuntime {
    pub manifest_path: String,
    pub manifest_schema: String,
    pub generated_readme: String,
    pub generation_report: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TutorialCatalogEntry {
    pub id: String,
    pub title: String,
    pub path: String,
    #[serde(rename = "type")]
    pub entry_type: String,
    pub status: String,
    pub source: String,
    #[serde(default)]
    pub audiences: Vec<String>,
    #[serde(default)]
    pub notes: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TutorialCatalog {
    #[serde(default = "default_tutorial_catalog_schema")]
    pub schema: String,
    #[serde(default)]
    pub description: String,
    #[serde(default)]
    pub entry_page: String,
    pub generated_runtime: TutorialCatalogGeneratedRuntime,
    #[serde(default)]
    pub entries: Vec<TutorialCatalogEntry>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TutorialCatalogMeta {
    #[serde(default = "default_tutorial_catalog_meta_schema")]
    pub schema: String,
    #[serde(default)]
    pub description: String,
    pub entry_page: String,
    pub generated_runtime: TutorialCatalogGeneratedRuntime,
    #[serde(default)]
    pub concepts: Vec<TutorialConcept>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TutorialSourceCatalogSection {
    pub order: usize,
    pub path: String,
    #[serde(rename = "type")]
    pub entry_type: String,
    pub status: String,
    pub source: String,
    #[serde(default)]
    pub audiences: Vec<String>,
    #[serde(default)]
    pub notes: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TutorialSourceGeneratedChapterSection {
    pub order: usize,
    pub example_id: String,
    pub tier: TutorialTier,
    #[serde(default)]
    pub guide_path: Option<String>,
    #[serde(default)]
    pub summary: String,
    #[serde(default)]
    pub narrative: String,
    #[serde(default)]
    pub use_cases: Vec<String>,
    #[serde(default)]
    pub gui_steps: Vec<String>,
    #[serde(default)]
    pub learning_objectives: Vec<String>,
    #[serde(default)]
    pub concepts: Vec<String>,
    #[serde(default)]
    pub parameter_notes: Vec<TutorialParameterNote>,
    #[serde(default)]
    pub follow_up_commands: Vec<String>,
    #[serde(default)]
    pub retain_outputs: Vec<String>,
    #[serde(default)]
    pub checkpoints: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TutorialSourceUnit {
    #[serde(default = "default_tutorial_source_schema")]
    pub schema: String,
    pub id: String,
    pub title: String,
    #[serde(default)]
    pub catalog: Option<TutorialSourceCatalogSection>,
    #[serde(default)]
    pub generated_chapter: Option<TutorialSourceGeneratedChapterSection>,
}

impl TutorialSourceCatalogSection {
    fn into_catalog_entry(self, id: String, title: String) -> TutorialCatalogEntry {
        TutorialCatalogEntry {
            id,
            title,
            path: self.path,
            entry_type: self.entry_type,
            status: self.status,
            source: self.source,
            audiences: self.audiences,
            notes: self.notes,
        }
    }
}

impl TutorialSourceGeneratedChapterSection {
    fn into_manifest_chapter(self, id: String, title: String) -> TutorialChapter {
        TutorialChapter {
            id,
            order: self.order,
            title,
            example_id: self.example_id,
            tier: self.tier,
            guide_path: self.guide_path,
            summary: self.summary,
            narrative: self.narrative,
            use_cases: self.use_cases,
            gui_steps: self.gui_steps,
            learning_objectives: self.learning_objectives,
            concepts: self.concepts,
            parameter_notes: self.parameter_notes,
            follow_up_commands: self.follow_up_commands,
            retain_outputs: self.retain_outputs,
            checkpoints: self.checkpoints,
        }
    }
}

impl TutorialSourceUnit {
    fn into_catalog_entry(self) -> Option<(usize, TutorialCatalogEntry)> {
        let TutorialSourceUnit {
            id, title, catalog, ..
        } = self;
        catalog.map(|catalog| {
            let order = catalog.order;
            let entry = catalog.into_catalog_entry(id, title);
            (order, entry)
        })
    }

    fn into_manifest_chapter(self) -> Option<TutorialChapter> {
        let TutorialSourceUnit {
            id,
            title,
            generated_chapter,
            ..
        } = self;
        generated_chapter.map(|generated| generated.into_manifest_chapter(id, title))
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TutorialConcept {
    pub id: String,
    pub name: String,
    #[serde(default)]
    pub description: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TutorialParameterNote {
    pub parameter: String,
    #[serde(default)]
    pub where_used: String,
    #[serde(default)]
    pub why_it_matters: String,
    #[serde(default)]
    pub how_to_derive: String,
    #[serde(default)]
    pub omit_when: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum TutorialTier {
    Core,
    Advanced,
    Online,
}

impl TutorialTier {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Core => "core",
            Self::Advanced => "advanced",
            Self::Online => "online",
        }
    }

    pub fn should_execute(self, online_enabled: bool) -> bool {
        match self {
            Self::Core | Self::Advanced => true,
            Self::Online => online_enabled,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TutorialChapter {
    pub id: String,
    pub order: usize,
    pub title: String,
    pub example_id: String,
    pub tier: TutorialTier,
    #[serde(default)]
    pub guide_path: Option<String>,
    #[serde(default)]
    pub summary: String,
    #[serde(default)]
    pub narrative: String,
    #[serde(default)]
    pub use_cases: Vec<String>,
    #[serde(default)]
    pub gui_steps: Vec<String>,
    #[serde(default)]
    pub learning_objectives: Vec<String>,
    #[serde(default)]
    pub concepts: Vec<String>,
    #[serde(default)]
    pub parameter_notes: Vec<TutorialParameterNote>,
    #[serde(default)]
    pub follow_up_commands: Vec<String>,
    #[serde(default)]
    pub retain_outputs: Vec<String>,
    #[serde(default)]
    pub checkpoints: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TutorialManifest {
    #[serde(default = "default_tutorial_manifest_schema")]
    pub schema: String,
    #[serde(default)]
    pub description: String,
    #[serde(default)]
    pub concepts: Vec<TutorialConcept>,
    pub chapters: Vec<TutorialChapter>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TutorialGenerationChapter {
    pub id: String,
    pub order: usize,
    pub title: String,
    pub tier: String,
    pub example_id: String,
    pub example_source: String,
    pub concepts: Vec<String>,
    pub executed: bool,
    pub retained_artifacts: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TutorialGenerationReport {
    #[serde(default = "default_tutorial_generation_report_schema")]
    pub schema: String,
    pub manifest: String,
    pub source_examples: String,
    pub chapter_count: usize,
    pub online_enabled: bool,
    pub generated_files: Vec<String>,
    pub file_checksums: BTreeMap<String, String>,
    pub chapters: Vec<TutorialGenerationChapter>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ExampleTestMode {
    Always,
    Skip,
    Online,
}

impl Default for ExampleTestMode {
    fn default() -> Self {
        Self::Always
    }
}

impl ExampleTestMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Always => "always",
            Self::Skip => "skip",
            Self::Online => "online",
        }
    }

    pub fn should_run(self, online_enabled: bool) -> bool {
        match self {
            Self::Always => true,
            Self::Skip => false,
            Self::Online => online_enabled,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WorkflowExample {
    #[serde(default = "default_example_schema")]
    pub schema: String,
    pub id: String,
    pub title: String,
    #[serde(default)]
    pub summary: String,
    #[serde(default)]
    pub test_mode: ExampleTestMode,
    #[serde(default)]
    pub required_files: Vec<String>,
    #[serde(default)]
    pub tags: Vec<String>,
    pub workflow: Workflow,
}

fn default_example_schema() -> String {
    WORKFLOW_EXAMPLE_SCHEMA.to_string()
}

fn default_tutorial_catalog_schema() -> String {
    TUTORIAL_CATALOG_SCHEMA.to_string()
}

fn default_tutorial_catalog_meta_schema() -> String {
    TUTORIAL_CATALOG_META_SCHEMA.to_string()
}

fn default_tutorial_source_schema() -> String {
    TUTORIAL_SOURCE_SCHEMA.to_string()
}

fn default_tutorial_manifest_schema() -> String {
    TUTORIAL_MANIFEST_SCHEMA.to_string()
}

fn default_tutorial_generation_report_schema() -> String {
    TUTORIAL_GENERATION_REPORT_SCHEMA.to_string()
}

#[derive(Debug, Clone)]
pub struct LoadedWorkflowExample {
    pub path: PathBuf,
    pub example: WorkflowExample,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WorkflowExampleDocReport {
    pub source_dir: String,
    pub output_dir: String,
    pub example_count: usize,
    pub generated_files: Vec<String>,
}

fn display_path(path: &Path) -> String {
    path.to_string_lossy().replace('\\', "/")
}

fn markdown_file_stem(id: &str) -> String {
    id.chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '-' || c == '_' {
                c
            } else {
                '_'
            }
        })
        .collect()
}

fn parse_example(path: &Path) -> Result<WorkflowExample, String> {
    let raw = fs::read_to_string(path)
        .map_err(|e| format!("Could not read example '{}': {e}", display_path(path)))?;
    serde_json::from_str::<WorkflowExample>(&raw)
        .map_err(|e| format!("Could not parse example '{}': {e}", display_path(path)))
}

fn validate_example(path: &Path, example: &WorkflowExample) -> Result<(), String> {
    if example.schema != WORKFLOW_EXAMPLE_SCHEMA {
        return Err(format!(
            "Example '{}' uses unsupported schema '{}'; expected '{}'",
            display_path(path),
            example.schema,
            WORKFLOW_EXAMPLE_SCHEMA
        ));
    }
    if example.id.trim().is_empty() {
        return Err(format!("Example '{}' has empty id", display_path(path)));
    }
    if example.title.trim().is_empty() {
        return Err(format!("Example '{}' has empty title", display_path(path)));
    }
    if example.workflow.run_id.trim().is_empty() {
        return Err(format!(
            "Example '{}' has empty workflow.run_id",
            display_path(path)
        ));
    }
    if example.workflow.ops.is_empty() {
        return Err(format!(
            "Example '{}' has empty workflow.ops",
            display_path(path)
        ));
    }
    for required in &example.required_files {
        if required.trim().is_empty() {
            return Err(format!(
                "Example '{}' has blank entry in required_files",
                display_path(path)
            ));
        }
    }
    Ok(())
}

pub fn load_workflow_examples(source_dir: &Path) -> Result<Vec<LoadedWorkflowExample>, String> {
    let mut paths = vec![];
    for entry in fs::read_dir(source_dir)
        .map_err(|e| format!("Could not read '{}': {e}", display_path(source_dir)))?
    {
        let entry = entry.map_err(|e| format!("Could not read directory entry: {e}"))?;
        let path = entry.path();
        if !path.is_file() {
            continue;
        }
        if path
            .extension()
            .map(|ext| ext.eq_ignore_ascii_case("json"))
            .unwrap_or(false)
        {
            paths.push(path);
        }
    }
    paths.sort();
    if paths.is_empty() {
        return Err(format!(
            "No example JSON files found in '{}'",
            display_path(source_dir)
        ));
    }
    let mut seen_ids: HashSet<String> = HashSet::new();
    let mut loaded = vec![];
    for path in paths {
        let example = parse_example(&path)?;
        validate_example(&path, &example)?;
        if !seen_ids.insert(example.id.clone()) {
            return Err(format!(
                "Duplicate example id '{}' (found in '{}')",
                example.id,
                display_path(&path)
            ));
        }
        loaded.push(LoadedWorkflowExample { path, example });
    }
    Ok(loaded)
}

pub fn validate_example_required_files(
    example: &WorkflowExample,
    repo_root: &Path,
) -> Result<(), String> {
    for required in &example.required_files {
        let path = Path::new(required);
        let resolved = if path.is_absolute() {
            path.to_path_buf()
        } else {
            repo_root.join(path)
        };
        if !resolved.exists() {
            return Err(format!(
                "Required example file '{}' not found (example id '{}')",
                display_path(&resolved),
                example.id
            ));
        }
    }
    Ok(())
}

fn bool_env_enabled(key: &str) -> bool {
    let raw = env::var(key).unwrap_or_default();
    matches!(
        raw.trim().to_ascii_lowercase().as_str(),
        "1" | "true" | "yes" | "on"
    )
}

pub fn online_example_tests_enabled() -> bool {
    bool_env_enabled(ONLINE_EXAMPLE_TEST_ENV) && !bool_env_enabled(SKIP_REMOTE_TESTS_ENV)
}

pub fn run_example_workflow(example: &WorkflowExample) -> Result<(), String> {
    let temp = TempDir::new().map_err(|e| format!("Could not create temp directory: {e}"))?;
    run_example_workflow_in_dir(example, Path::new("."), temp.path())
}

/// Executes one workflow example with path rewriting and returns the resulting
/// project state.
pub fn run_example_workflow_for_project_state(
    example: &WorkflowExample,
    repo_root: &Path,
    run_dir: &Path,
) -> Result<ProjectState, String> {
    let mut noop = |_p: OperationProgress| true;
    run_example_workflow_for_project_state_with_progress(example, repo_root, run_dir, &mut noop)
}

/// Executes one workflow example with path rewriting, streams shared engine
/// progress callbacks, and returns the resulting project state.
pub fn run_example_workflow_for_project_state_with_progress<F>(
    example: &WorkflowExample,
    repo_root: &Path,
    run_dir: &Path,
    mut on_progress: F,
) -> Result<ProjectState, String>
where
    F: FnMut(OperationProgress) -> bool,
{
    let rewritten = rewrite_example_paths_for_execution(example, repo_root, run_dir)?;
    let mut engine = GentleEngine::from_state(ProjectState::default());
    engine
        .apply_workflow_with_progress(rewritten.workflow, |progress| on_progress(progress))
        .map_err(|e| format!("Workflow example '{}' failed: {e}", example.id))?;
    Ok(engine.state().clone())
}

fn parse_tutorial_manifest(manifest_path: &Path) -> Result<TutorialManifest, String> {
    let raw = fs::read_to_string(manifest_path).map_err(|e| {
        format!(
            "Could not read tutorial manifest '{}': {e}",
            display_path(manifest_path)
        )
    })?;
    serde_json::from_str::<TutorialManifest>(&raw).map_err(|e| {
        format!(
            "Could not parse tutorial manifest '{}': {e}",
            display_path(manifest_path)
        )
    })
}

fn parse_tutorial_catalog(catalog_path: &Path) -> Result<TutorialCatalog, String> {
    let raw = fs::read_to_string(catalog_path).map_err(|e| {
        format!(
            "Could not read tutorial catalog '{}': {e}",
            display_path(catalog_path)
        )
    })?;
    serde_json::from_str::<TutorialCatalog>(&raw).map_err(|e| {
        format!(
            "Could not parse tutorial catalog '{}': {e}",
            display_path(catalog_path)
        )
    })
}

fn parse_tutorial_catalog_meta(meta_path: &Path) -> Result<TutorialCatalogMeta, String> {
    let raw = fs::read_to_string(meta_path).map_err(|e| {
        format!(
            "Could not read tutorial catalog meta '{}': {e}",
            display_path(meta_path)
        )
    })?;
    serde_json::from_str::<TutorialCatalogMeta>(&raw).map_err(|e| {
        format!(
            "Could not parse tutorial catalog meta '{}': {e}",
            display_path(meta_path)
        )
    })
}

fn parse_tutorial_source_unit(path: &Path) -> Result<TutorialSourceUnit, String> {
    let raw = fs::read_to_string(path).map_err(|e| {
        format!(
            "Could not read tutorial source '{}': {e}",
            display_path(path)
        )
    })?;
    serde_json::from_str::<TutorialSourceUnit>(&raw).map_err(|e| {
        format!(
            "Could not parse tutorial source '{}': {e}",
            display_path(path)
        )
    })
}

pub fn load_tutorial_catalog(catalog_path: &Path) -> Result<TutorialCatalog, String> {
    let catalog = parse_tutorial_catalog(catalog_path)?;
    if catalog.schema != TUTORIAL_CATALOG_SCHEMA {
        return Err(format!(
            "Tutorial catalog '{}' uses unsupported schema '{}'; expected '{}'",
            display_path(catalog_path),
            catalog.schema,
            TUTORIAL_CATALOG_SCHEMA
        ));
    }
    if catalog.entry_page.trim().is_empty() {
        return Err(format!(
            "Tutorial catalog '{}' has empty entry_page",
            display_path(catalog_path)
        ));
    }
    if catalog.generated_runtime.manifest_path.trim().is_empty()
        || catalog.generated_runtime.manifest_schema.trim().is_empty()
        || catalog.generated_runtime.generated_readme.trim().is_empty()
        || catalog
            .generated_runtime
            .generation_report
            .trim()
            .is_empty()
    {
        return Err(format!(
            "Tutorial catalog '{}' has incomplete generated_runtime fields",
            display_path(catalog_path)
        ));
    }
    if catalog.entries.is_empty() {
        return Err(format!(
            "Tutorial catalog '{}' has no entries",
            display_path(catalog_path)
        ));
    }
    let mut seen_ids: HashSet<String> = HashSet::new();
    for entry in &catalog.entries {
        if entry.id.trim().is_empty() {
            return Err(format!(
                "Tutorial catalog '{}' contains entry with empty id",
                display_path(catalog_path)
            ));
        }
        if entry.title.trim().is_empty() {
            return Err(format!(
                "Tutorial catalog entry '{}' has empty title",
                entry.id
            ));
        }
        if entry.path.trim().is_empty() {
            return Err(format!(
                "Tutorial catalog entry '{}' has empty path",
                entry.id
            ));
        }
        if entry.entry_type.trim().is_empty() {
            return Err(format!(
                "Tutorial catalog entry '{}' has empty type",
                entry.id
            ));
        }
        if entry.status.trim().is_empty() {
            return Err(format!(
                "Tutorial catalog entry '{}' has empty status",
                entry.id
            ));
        }
        if entry.source.trim().is_empty() {
            return Err(format!(
                "Tutorial catalog entry '{}' has empty source",
                entry.id
            ));
        }
        if !seen_ids.insert(entry.id.clone()) {
            return Err(format!(
                "Tutorial catalog '{}' has duplicate entry id '{}'",
                display_path(catalog_path),
                entry.id
            ));
        }
    }
    Ok(catalog)
}

pub fn load_tutorial_catalog_meta(meta_path: &Path) -> Result<TutorialCatalogMeta, String> {
    let meta = parse_tutorial_catalog_meta(meta_path)?;
    if meta.schema != TUTORIAL_CATALOG_META_SCHEMA {
        return Err(format!(
            "Tutorial catalog meta '{}' uses unsupported schema '{}'; expected '{}'",
            display_path(meta_path),
            meta.schema,
            TUTORIAL_CATALOG_META_SCHEMA
        ));
    }
    if meta.entry_page.trim().is_empty() {
        return Err(format!(
            "Tutorial catalog meta '{}' has empty entry_page",
            display_path(meta_path)
        ));
    }
    if meta.generated_runtime.manifest_path.trim().is_empty()
        || meta.generated_runtime.manifest_schema.trim().is_empty()
        || meta.generated_runtime.generated_readme.trim().is_empty()
        || meta.generated_runtime.generation_report.trim().is_empty()
    {
        return Err(format!(
            "Tutorial catalog meta '{}' has incomplete generated_runtime fields",
            display_path(meta_path)
        ));
    }
    Ok(meta)
}

pub fn load_tutorial_source_units(source_dir: &Path) -> Result<Vec<TutorialSourceUnit>, String> {
    let entries = fs::read_dir(source_dir).map_err(|e| {
        format!(
            "Could not read tutorial source directory '{}': {e}",
            display_path(source_dir)
        )
    })?;
    let mut json_paths = entries
        .flatten()
        .map(|entry| entry.path())
        .filter(|path| {
            path.is_file()
                && path
                    .extension()
                    .and_then(|ext| ext.to_str())
                    .map(|ext| ext.eq_ignore_ascii_case("json"))
                    .unwrap_or(false)
                && path
                    .file_name()
                    .and_then(|name| name.to_str())
                    .map(|name| name != "catalog_meta.json")
                    .unwrap_or(true)
        })
        .collect::<Vec<_>>();
    json_paths.sort();
    if json_paths.is_empty() {
        return Err(format!(
            "Tutorial source directory '{}' has no source JSON files",
            display_path(source_dir)
        ));
    }
    let mut units = Vec::new();
    let mut seen_ids: HashSet<String> = HashSet::new();
    let mut seen_catalog_orders: HashSet<usize> = HashSet::new();
    let mut seen_chapter_orders: HashSet<usize> = HashSet::new();
    for path in json_paths {
        let unit = parse_tutorial_source_unit(&path)?;
        if unit.schema != TUTORIAL_SOURCE_SCHEMA {
            return Err(format!(
                "Tutorial source '{}' uses unsupported schema '{}'; expected '{}'",
                display_path(&path),
                unit.schema,
                TUTORIAL_SOURCE_SCHEMA
            ));
        }
        if unit.id.trim().is_empty() || unit.title.trim().is_empty() {
            return Err(format!(
                "Tutorial source '{}' has one or more empty required fields",
                display_path(&path)
            ));
        }
        if unit.catalog.is_none() && unit.generated_chapter.is_none() {
            return Err(format!(
                "Tutorial source '{}' must define at least one of `catalog` or `generated_chapter`",
                display_path(&path)
            ));
        }
        if !seen_ids.insert(unit.id.clone()) {
            return Err(format!(
                "Tutorial source directory '{}' has duplicate tutorial id '{}'",
                display_path(source_dir),
                unit.id
            ));
        }
        if let Some(catalog) = &unit.catalog {
            if catalog.path.trim().is_empty()
                || catalog.entry_type.trim().is_empty()
                || catalog.status.trim().is_empty()
                || catalog.source.trim().is_empty()
            {
                return Err(format!(
                    "Tutorial source '{}' has incomplete `catalog` fields",
                    display_path(&path)
                ));
            }
            if !seen_catalog_orders.insert(catalog.order) {
                return Err(format!(
                    "Tutorial source directory '{}' has duplicate tutorial catalog order '{}'",
                    display_path(source_dir),
                    catalog.order
                ));
            }
        }
        if let Some(generated) = &unit.generated_chapter {
            if generated.example_id.trim().is_empty() {
                return Err(format!(
                    "Tutorial source '{}' has incomplete `generated_chapter` fields",
                    display_path(&path)
                ));
            }
            if !seen_chapter_orders.insert(generated.order) {
                return Err(format!(
                    "Tutorial source directory '{}' has duplicate generated chapter order '{}'",
                    display_path(source_dir),
                    generated.order
                ));
            }
        }
        units.push(unit);
    }
    units.sort_by(|left, right| left.id.cmp(&right.id));
    Ok(units)
}

pub fn generate_tutorial_catalog_from_sources(
    meta_path: &Path,
    source_dir: &Path,
) -> Result<TutorialCatalog, String> {
    let meta = load_tutorial_catalog_meta(meta_path)?;
    let units = load_tutorial_source_units(source_dir)?;
    let entries = units
        .into_iter()
        .filter_map(TutorialSourceUnit::into_catalog_entry)
        .collect::<Vec<_>>();
    let mut ordered_entries = entries;
    ordered_entries.sort_by(|left, right| {
        left.0
            .cmp(&right.0)
            .then_with(|| left.1.id.cmp(&right.1.id))
    });
    Ok(TutorialCatalog {
        schema: TUTORIAL_CATALOG_SCHEMA.to_string(),
        description: meta.description,
        entry_page: meta.entry_page,
        generated_runtime: meta.generated_runtime,
        entries: ordered_entries
            .into_iter()
            .map(|(_, entry)| entry)
            .collect(),
    })
}

pub fn generate_tutorial_manifest_from_sources(
    meta_path: &Path,
    source_dir: &Path,
) -> Result<TutorialManifest, String> {
    let meta = load_tutorial_catalog_meta(meta_path)?;
    if meta.concepts.is_empty() {
        return Err(format!(
            "Tutorial catalog meta '{}' must define tutorial concepts for manifest generation",
            display_path(meta_path)
        ));
    }
    let units = load_tutorial_source_units(source_dir)?;
    let mut chapters = units
        .into_iter()
        .filter_map(TutorialSourceUnit::into_manifest_chapter)
        .collect::<Vec<_>>();
    chapters.sort_by(|left, right| {
        left.order
            .cmp(&right.order)
            .then_with(|| left.id.cmp(&right.id))
    });
    let manifest = TutorialManifest {
        schema: TUTORIAL_MANIFEST_SCHEMA.to_string(),
        description: "Generated executable tutorial runtime manifest for GENtle. Used by `tutorial-generate` / `tutorial-check` and the GUI `Open Tutorial Project...` flow to map source-backed chapters onto canonical workflow examples.".to_string(),
        concepts: meta.concepts,
        chapters,
    };
    Ok(manifest)
}

pub fn write_tutorial_catalog_from_sources(
    meta_path: &Path,
    source_dir: &Path,
    output_path: &Path,
) -> Result<TutorialCatalog, String> {
    let catalog = generate_tutorial_catalog_from_sources(meta_path, source_dir)?;
    let serialized = serde_json::to_string_pretty(&catalog)
        .map_err(|e| format!("Could not serialize tutorial catalog: {e}"))?;
    fs::write(output_path, format!("{serialized}\n")).map_err(|e| {
        format!(
            "Could not write tutorial catalog '{}': {e}",
            display_path(output_path)
        )
    })?;
    Ok(catalog)
}

pub fn write_tutorial_manifest_from_sources(
    meta_path: &Path,
    source_dir: &Path,
    output_path: &Path,
) -> Result<TutorialManifest, String> {
    let manifest = generate_tutorial_manifest_from_sources(meta_path, source_dir)?;
    let serialized = serde_json::to_string_pretty(&manifest)
        .map_err(|e| format!("Could not serialize tutorial manifest: {e}"))?;
    fs::write(output_path, format!("{serialized}\n")).map_err(|e| {
        format!(
            "Could not write tutorial manifest '{}': {e}",
            display_path(output_path)
        )
    })?;
    Ok(manifest)
}

pub fn check_tutorial_catalog_generated(
    meta_path: &Path,
    source_dir: &Path,
    output_path: &Path,
) -> Result<TutorialCatalog, String> {
    if !output_path.exists() {
        return Err(format!(
            "Expected tutorial catalog '{}' does not exist. Run `tutorial-catalog-generate` first.",
            display_path(output_path)
        ));
    }
    let expected = generate_tutorial_catalog_from_sources(meta_path, source_dir)?;
    let expected_bytes = format!(
        "{}\n",
        serde_json::to_string_pretty(&expected)
            .map_err(|e| format!("Could not serialize expected tutorial catalog: {e}"))?
    )
    .into_bytes();
    let actual_bytes = fs::read(output_path).map_err(|e| {
        format!(
            "Could not read tutorial catalog '{}': {e}",
            display_path(output_path)
        )
    })?;
    if actual_bytes != expected_bytes {
        return Err(format!(
            "Committed tutorial catalog '{}' does not match generated output. Run `tutorial-catalog-generate`.",
            display_path(output_path)
        ));
    }
    Ok(expected)
}

pub fn check_tutorial_manifest_generated(
    meta_path: &Path,
    source_dir: &Path,
    output_path: &Path,
) -> Result<TutorialManifest, String> {
    if !output_path.exists() {
        return Err(format!(
            "Expected tutorial manifest '{}' does not exist. Run `tutorial-manifest-generate` first.",
            display_path(output_path)
        ));
    }
    let expected = generate_tutorial_manifest_from_sources(meta_path, source_dir)?;
    let expected_bytes = format!(
        "{}\n",
        serde_json::to_string_pretty(&expected)
            .map_err(|e| format!("Could not serialize expected tutorial manifest: {e}"))?
    )
    .into_bytes();
    let actual_bytes = fs::read(output_path).map_err(|e| {
        format!(
            "Could not read tutorial manifest '{}': {e}",
            display_path(output_path)
        )
    })?;
    if actual_bytes != expected_bytes {
        return Err(format!(
            "Committed tutorial manifest '{}' does not match generated output. Run `tutorial-manifest-generate`.",
            display_path(output_path)
        ));
    }
    Ok(expected)
}

pub fn load_tutorial_manifest(manifest_path: &Path) -> Result<TutorialManifest, String> {
    let manifest = parse_tutorial_manifest(manifest_path)?;
    if manifest.schema != TUTORIAL_MANIFEST_SCHEMA {
        return Err(format!(
            "Tutorial manifest '{}' uses unsupported schema '{}'; expected '{}'",
            display_path(manifest_path),
            manifest.schema,
            TUTORIAL_MANIFEST_SCHEMA
        ));
    }
    if manifest.description.trim().is_empty() {
        return Err(format!(
            "Tutorial manifest '{}' has empty description",
            display_path(manifest_path)
        ));
    }
    if manifest.chapters.is_empty() {
        return Err(format!(
            "Tutorial manifest '{}' has no chapters",
            display_path(manifest_path)
        ));
    }
    if manifest.concepts.is_empty() {
        return Err(format!(
            "Tutorial manifest '{}' has no concept definitions",
            display_path(manifest_path)
        ));
    }
    let concept_by_id = tutorial_concept_lookup(&manifest.concepts).map_err(|e| {
        format!(
            "Tutorial manifest '{}' concept validation error: {e}",
            display_path(manifest_path)
        )
    })?;
    let mut seen_ids: HashSet<String> = HashSet::new();
    let mut seen_orders: HashSet<usize> = HashSet::new();
    for chapter in &manifest.chapters {
        if chapter.id.trim().is_empty() {
            return Err(format!(
                "Tutorial manifest '{}' contains chapter with empty id",
                display_path(manifest_path)
            ));
        }
        if chapter.title.trim().is_empty() {
            return Err(format!("Tutorial chapter '{}' has empty title", chapter.id));
        }
        if chapter.example_id.trim().is_empty() {
            return Err(format!(
                "Tutorial chapter '{}' has empty example_id",
                chapter.id
            ));
        }
        if chapter.learning_objectives.is_empty() {
            return Err(format!(
                "Tutorial chapter '{}' must define at least one learning objective",
                chapter.id
            ));
        }
        if chapter.use_cases.is_empty() {
            return Err(format!(
                "Tutorial chapter '{}' must define at least one use-case context",
                chapter.id
            ));
        }
        if chapter.tier != TutorialTier::Online && chapter.gui_steps.is_empty() {
            return Err(format!(
                "Tutorial chapter '{}' must define at least one GUI-first step",
                chapter.id
            ));
        }
        if chapter.concepts.is_empty() {
            return Err(format!(
                "Tutorial chapter '{}' must reference at least one concept",
                chapter.id
            ));
        }
        if let Some(guide_path) = &chapter.guide_path {
            if guide_path.trim().is_empty() {
                return Err(format!(
                    "Tutorial chapter '{}' contains blank guide_path",
                    chapter.id
                ));
            }
        }
        if !seen_ids.insert(chapter.id.clone()) {
            return Err(format!(
                "Tutorial manifest '{}' has duplicate chapter id '{}'",
                display_path(manifest_path),
                chapter.id
            ));
        }
        if !seen_orders.insert(chapter.order) {
            return Err(format!(
                "Tutorial manifest '{}' has duplicate chapter order '{}'",
                display_path(manifest_path),
                chapter.order
            ));
        }
        for retain in &chapter.retain_outputs {
            if retain.trim().is_empty() {
                return Err(format!(
                    "Tutorial chapter '{}' contains blank retain_outputs entry",
                    chapter.id
                ));
            }
        }
        for checkpoint in &chapter.checkpoints {
            if checkpoint.trim().is_empty() {
                return Err(format!(
                    "Tutorial chapter '{}' contains blank checkpoints entry",
                    chapter.id
                ));
            }
        }
        for objective in &chapter.learning_objectives {
            if objective.trim().is_empty() {
                return Err(format!(
                    "Tutorial chapter '{}' contains blank learning objective",
                    chapter.id
                ));
            }
        }
        for use_case in &chapter.use_cases {
            if use_case.trim().is_empty() {
                return Err(format!(
                    "Tutorial chapter '{}' contains blank use_case entry",
                    chapter.id
                ));
            }
        }
        for gui_step in &chapter.gui_steps {
            if gui_step.trim().is_empty() {
                return Err(format!(
                    "Tutorial chapter '{}' contains blank gui_steps entry",
                    chapter.id
                ));
            }
        }
        for concept_id in &chapter.concepts {
            if concept_id.trim().is_empty() {
                return Err(format!(
                    "Tutorial chapter '{}' contains blank concept id",
                    chapter.id
                ));
            }
            if !concept_by_id.contains_key(concept_id) {
                return Err(format!(
                    "Tutorial chapter '{}' references unknown concept id '{}'",
                    chapter.id, concept_id
                ));
            }
        }
        for note in &chapter.parameter_notes {
            if note.parameter.trim().is_empty() {
                return Err(format!(
                    "Tutorial chapter '{}' contains parameter note with empty parameter",
                    chapter.id
                ));
            }
            if note.where_used.trim().is_empty()
                && note.why_it_matters.trim().is_empty()
                && note.how_to_derive.trim().is_empty()
                && note.omit_when.trim().is_empty()
            {
                return Err(format!(
                    "Tutorial chapter '{}' parameter note '{}' needs explanatory text",
                    chapter.id, note.parameter
                ));
            }
        }
        for command in &chapter.follow_up_commands {
            if command.trim().is_empty() {
                return Err(format!(
                    "Tutorial chapter '{}' contains blank follow_up_commands entry",
                    chapter.id
                ));
            }
        }
    }
    Ok(manifest)
}

fn example_lookup(examples: &[LoadedWorkflowExample]) -> HashMap<String, LoadedWorkflowExample> {
    let mut by_id = HashMap::new();
    for loaded in examples {
        by_id.insert(loaded.example.id.clone(), loaded.clone());
    }
    by_id
}

fn tutorial_concept_lookup(
    concepts: &[TutorialConcept],
) -> Result<HashMap<String, TutorialConcept>, String> {
    let mut by_id = HashMap::new();
    for concept in concepts {
        if concept.id.trim().is_empty() {
            return Err("Tutorial concept id must not be empty".to_string());
        }
        if concept.name.trim().is_empty() {
            return Err(format!("Tutorial concept '{}' has empty name", concept.id));
        }
        if by_id.insert(concept.id.clone(), concept.clone()).is_some() {
            return Err(format!("Duplicate tutorial concept id '{}'", concept.id));
        }
    }
    Ok(by_id)
}

fn concept_occurrences(sorted_chapters: &[TutorialChapter]) -> HashMap<String, Vec<usize>> {
    let mut occurrences: HashMap<String, Vec<usize>> = HashMap::new();
    for (chapter_index, chapter) in sorted_chapters.iter().enumerate() {
        for concept_id in &chapter.concepts {
            occurrences
                .entry(concept_id.clone())
                .or_default()
                .push(chapter_index);
        }
    }
    occurrences
}

fn chapter_link_label(chapter: &TutorialChapter) -> String {
    format!("Chapter {}: {}", chapter.order, chapter.title)
}

fn chapter_link_target(chapter: &TutorialChapter, from_readme: bool) -> String {
    let filename = chapter_markdown_filename(chapter.order, &chapter.id);
    if from_readme {
        format!("./chapters/{filename}")
    } else {
        format!("./{filename}")
    }
}

pub fn validate_tutorial_manifest_against_examples(
    manifest: &TutorialManifest,
    examples: &[LoadedWorkflowExample],
) -> Result<(), String> {
    let by_id = example_lookup(examples);
    let concept_by_id = tutorial_concept_lookup(&manifest.concepts)?;
    let mut concept_used: HashSet<String> = HashSet::new();
    for chapter in &manifest.chapters {
        let loaded = by_id.get(&chapter.example_id).ok_or_else(|| {
            format!(
                "Tutorial chapter '{}' references unknown example id '{}'",
                chapter.id, chapter.example_id
            )
        })?;
        if chapter.tier == TutorialTier::Core && loaded.example.test_mode != ExampleTestMode::Always
        {
            return Err(format!(
                "Tutorial chapter '{}' is tier 'core' but example '{}' is test_mode '{}'",
                chapter.id,
                chapter.example_id,
                loaded.example.test_mode.as_str()
            ));
        }
        if chapter.tier == TutorialTier::Online
            && loaded.example.test_mode != ExampleTestMode::Online
        {
            return Err(format!(
                "Tutorial chapter '{}' is tier 'online' but example '{}' is test_mode '{}'",
                chapter.id,
                chapter.example_id,
                loaded.example.test_mode.as_str()
            ));
        }
        for concept_id in &chapter.concepts {
            if !concept_by_id.contains_key(concept_id) {
                return Err(format!(
                    "Tutorial chapter '{}' references unknown concept id '{}'",
                    chapter.id, concept_id
                ));
            }
            concept_used.insert(concept_id.clone());
        }
    }
    for concept_id in concept_by_id.keys() {
        if !concept_used.contains(concept_id) {
            return Err(format!(
                "Tutorial concept '{}' is defined but not used in any chapter",
                concept_id
            ));
        }
    }
    Ok(())
}

fn resolve_input_path(path: &str, repo_root: &Path) -> String {
    let raw = Path::new(path);
    if raw.is_absolute() {
        return display_path(raw);
    }
    display_path(&repo_root.join(raw))
}

fn resolve_output_path(path: &str, run_dir: &Path) -> String {
    let raw = Path::new(path);
    if raw.is_absolute() {
        return display_path(raw);
    }
    display_path(&run_dir.join(raw))
}

fn rewrite_optional_input_path(path: &mut Option<String>, repo_root: &Path) {
    if let Some(value) = path {
        *value = resolve_input_path(value, repo_root);
    }
}

fn rewrite_optional_output_path(path: &mut Option<String>, run_dir: &Path) {
    if let Some(value) = path {
        *value = resolve_output_path(value, run_dir);
    }
}

fn ensure_parent_exists(path: &str) -> Result<(), String> {
    let parent = Path::new(path).parent();
    if let Some(parent) = parent {
        if parent.as_os_str().is_empty() {
            return Ok(());
        }
        fs::create_dir_all(parent).map_err(|e| {
            format!(
                "Could not create output parent '{}': {e}",
                display_path(parent)
            )
        })?;
    }
    Ok(())
}

fn ensure_directory_exists(path: &str) -> Result<(), String> {
    let dir = Path::new(path);
    fs::create_dir_all(dir).map_err(|e| {
        format!(
            "Could not create output directory '{}': {e}",
            display_path(dir)
        )
    })?;
    Ok(())
}

fn rewrite_example_paths_for_execution(
    example: &WorkflowExample,
    repo_root: &Path,
    run_dir: &Path,
) -> Result<WorkflowExample, String> {
    let mut rewritten = example.clone();
    for op in &mut rewritten.workflow.ops {
        if let Operation::LoadFile { path, .. } = op {
            *path = resolve_input_path(path, repo_root);
            continue;
        }
        if let Operation::SaveFile { path, .. } = op {
            *path = resolve_output_path(path, run_dir);
            ensure_parent_exists(path)?;
            continue;
        }
        if let Operation::RenderSequenceSvg { path, .. } = op {
            *path = resolve_output_path(path, run_dir);
            ensure_parent_exists(path)?;
            continue;
        }
        if let Operation::RenderFeatureExpertSvg { path, .. } = op {
            *path = resolve_output_path(path, run_dir);
            ensure_parent_exists(path)?;
            continue;
        }
        if let Operation::RenderRnaStructureSvg { path, .. } = op {
            *path = resolve_output_path(path, run_dir);
            ensure_parent_exists(path)?;
            continue;
        }
        if let Operation::RenderLineageSvg { path } = op {
            *path = resolve_output_path(path, run_dir);
            ensure_parent_exists(path)?;
            continue;
        }
        if let Operation::RenderPoolGelSvg { path, .. } = op {
            *path = resolve_output_path(path, run_dir);
            ensure_parent_exists(path)?;
            continue;
        }
        if let Operation::RenderProteinGelSvg { path, .. } = op {
            *path = resolve_output_path(path, run_dir);
            ensure_parent_exists(path)?;
            continue;
        }
        if let Operation::RenderProtein2dGelSvg { path, .. } = op {
            *path = resolve_output_path(path, run_dir);
            ensure_parent_exists(path)?;
            continue;
        }
        if let Operation::ExportDnaLadders { path, .. } = op {
            *path = resolve_output_path(path, run_dir);
            ensure_parent_exists(path)?;
            continue;
        }
        if let Operation::ExportRnaLadders { path, .. } = op {
            *path = resolve_output_path(path, run_dir);
            ensure_parent_exists(path)?;
            continue;
        }
        if let Operation::ExportPool { path, .. } = op {
            *path = resolve_output_path(path, run_dir);
            ensure_parent_exists(path)?;
            continue;
        }
        if let Operation::FindRestrictionSites { path, .. } = op {
            rewrite_optional_output_path(path, run_dir);
            if let Some(path) = path.as_deref() {
                ensure_parent_exists(path)?;
            }
            continue;
        }
        if let Operation::ScanTfbsHits { path, .. } = op {
            rewrite_optional_output_path(path, run_dir);
            if let Some(path) = path.as_deref() {
                ensure_parent_exists(path)?;
            }
            continue;
        }
        if let Operation::SummarizeTfbsScoreTracks { path, .. } = op {
            rewrite_optional_output_path(path, run_dir);
            if let Some(path) = path.as_deref() {
                ensure_parent_exists(path)?;
            }
            continue;
        }
        if let Operation::RenderTfbsScoreTracksSvg { path, .. } = op {
            *path = resolve_output_path(path, run_dir);
            ensure_parent_exists(path)?;
            continue;
        }
        if let Operation::PrepareGenome {
            catalog_path,
            cache_dir,
            ..
        } = op
        {
            rewrite_optional_input_path(catalog_path, repo_root);
            rewrite_optional_output_path(cache_dir, run_dir);
            if let Some(dir) = cache_dir {
                ensure_directory_exists(dir)?;
            }
            continue;
        }
        if let Operation::ExtractGenomeRegion {
            catalog_path,
            cache_dir,
            ..
        } = op
        {
            rewrite_optional_input_path(catalog_path, repo_root);
            rewrite_optional_output_path(cache_dir, run_dir);
            if let Some(dir) = cache_dir {
                ensure_directory_exists(dir)?;
            }
            continue;
        }
        if let Operation::ExtractGenomeGene {
            catalog_path,
            cache_dir,
            ..
        } = op
        {
            rewrite_optional_input_path(catalog_path, repo_root);
            rewrite_optional_output_path(cache_dir, run_dir);
            if let Some(dir) = cache_dir {
                ensure_directory_exists(dir)?;
            }
            continue;
        }
        if let Operation::ExtendGenomeAnchor {
            catalog_path,
            cache_dir,
            ..
        } = op
        {
            rewrite_optional_input_path(catalog_path, repo_root);
            rewrite_optional_output_path(cache_dir, run_dir);
            if let Some(dir) = cache_dir {
                ensure_directory_exists(dir)?;
            }
            continue;
        }
        if let Operation::ImportGenomeBedTrack { path, .. } = op {
            *path = resolve_input_path(path, repo_root);
            continue;
        }
        if let Operation::ImportGenomeBigWigTrack { path, .. } = op {
            *path = resolve_input_path(path, repo_root);
            continue;
        }
        if let Operation::ImportGenomeVcfTrack { path, .. } = op {
            *path = resolve_input_path(path, repo_root);
            continue;
        }
        if let Operation::ImportIsoformPanel { panel_path, .. } = op {
            *panel_path = resolve_input_path(panel_path, repo_root);
            continue;
        }
        if let Operation::ExportGuideOligos { path, .. } = op {
            *path = resolve_output_path(path, run_dir);
            ensure_parent_exists(path)?;
            continue;
        }
        if let Operation::ExportGuideProtocolText { path, .. } = op {
            *path = resolve_output_path(path, run_dir);
            ensure_parent_exists(path)?;
            continue;
        }
        if let Operation::RenderIsoformArchitectureSvg { path, .. } = op {
            *path = resolve_output_path(path, run_dir);
            ensure_parent_exists(path)?;
            continue;
        }
    }
    Ok(rewritten)
}

fn run_example_workflow_in_dir(
    example: &WorkflowExample,
    repo_root: &Path,
    run_dir: &Path,
) -> Result<(), String> {
    let rewritten = rewrite_example_paths_for_execution(example, repo_root, run_dir)?;
    let mut engine = GentleEngine::from_state(ProjectState::default());
    engine
        .apply_workflow(rewritten.workflow)
        .map(|_| ())
        .map_err(|e| format!("Workflow example '{}' failed: {e}", example.id))
}

fn validate_relative_output_path(path: &str) -> Result<PathBuf, String> {
    let as_path = Path::new(path);
    if as_path.is_absolute() {
        return Err(format!("retain_outputs entry '{}' must be relative", path));
    }
    for component in as_path.components() {
        if matches!(component, std::path::Component::ParentDir) {
            return Err(format!(
                "retain_outputs entry '{}' must not contain '..'",
                path
            ));
        }
    }
    Ok(as_path.to_path_buf())
}

fn copy_file(src: &Path, dst: &Path) -> Result<(), String> {
    if let Some(parent) = dst.parent() {
        fs::create_dir_all(parent)
            .map_err(|e| format!("Could not create '{}': {e}", display_path(parent)))?;
    }
    fs::copy(src, dst).map_err(|e| {
        format!(
            "Could not copy '{}' to '{}': {e}",
            display_path(src),
            display_path(dst)
        )
    })?;
    Ok(())
}

fn copy_directory_recursive(src: &Path, dst: &Path) -> Result<(), String> {
    fs::create_dir_all(dst)
        .map_err(|e| format!("Could not create '{}': {e}", display_path(dst)))?;
    let mut entries = vec![];
    for entry in
        fs::read_dir(src).map_err(|e| format!("Could not read '{}': {e}", display_path(src)))?
    {
        let entry = entry.map_err(|e| format!("Could not read directory entry: {e}"))?;
        entries.push(entry.path());
    }
    entries.sort();
    for path in entries {
        let name = path.file_name().ok_or_else(|| {
            format!(
                "Could not determine file name while copying '{}'",
                display_path(&path)
            )
        })?;
        let target = dst.join(name);
        if path.is_dir() {
            copy_directory_recursive(&path, &target)?;
        } else if path.is_file() {
            copy_file(&path, &target)?;
        }
    }
    Ok(())
}

fn copy_retained_output(
    chapter: &TutorialChapter,
    run_dir: &Path,
    output_dir: &Path,
    retain_output: &str,
) -> Result<String, String> {
    let relative = validate_relative_output_path(retain_output)?;
    let source = run_dir.join(&relative);
    if !source.exists() {
        return Err(format!(
            "Tutorial chapter '{}' expected retained output '{}' in run directory, but it was not created",
            chapter.id, retain_output
        ));
    }
    let chapter_slug = markdown_file_stem(&chapter.id);
    let destination_relative = Path::new("artifacts").join(chapter_slug).join(&relative);
    let destination = output_dir.join(&destination_relative);
    if source.is_dir() {
        copy_directory_recursive(&source, &destination)?;
    } else {
        copy_file(&source, &destination)?;
    }
    Ok(display_path(&destination_relative))
}

fn sorted_manifest_chapters(manifest: &TutorialManifest) -> Vec<TutorialChapter> {
    let mut chapters = manifest.chapters.clone();
    chapters.sort_by(|a, b| a.order.cmp(&b.order).then_with(|| a.id.cmp(&b.id)));
    chapters
}

fn chapter_markdown_filename(order: usize, id: &str) -> String {
    format!("{:02}_{}.md", order, markdown_file_stem(id))
}

fn markdown_path_for_chapter(chapter: &TutorialChapter) -> String {
    chapter_markdown_filename(chapter.order, &chapter.id)
}

fn render_tutorial_chapter_markdown(
    chapter: &TutorialChapter,
    loaded: &LoadedWorkflowExample,
    executed: bool,
    retained_artifacts: &[String],
    online_enabled: bool,
    chapter_index: usize,
    sorted_chapters: &[TutorialChapter],
    concept_by_id: &HashMap<String, TutorialConcept>,
    occurrences_by_concept: &HashMap<String, Vec<usize>>,
) -> Result<String, String> {
    let source_path = display_path(&loaded.path);
    let mut out = String::new();
    out.push_str("# ");
    out.push_str(&chapter.title);
    out.push_str("\n\n");
    out.push_str("- Chapter id: `");
    out.push_str(&chapter.id);
    out.push_str("`\n");
    out.push_str("- Tier: `");
    out.push_str(chapter.tier.as_str());
    out.push_str("`\n");
    out.push_str("- Example id: `");
    out.push_str(&chapter.example_id);
    out.push_str("`\n");
    out.push_str("- Source example: `");
    out.push_str(&source_path);
    out.push_str("`\n");
    out.push_str("- Example test_mode: `");
    out.push_str(loaded.example.test_mode.as_str());
    out.push_str("`\n");
    out.push_str("- Executed during generation: `");
    out.push_str(if executed { "yes" } else { "no" });
    out.push_str("`\n");
    if chapter.tier == TutorialTier::Online && !online_enabled {
        out.push_str("- Execution note: set `GENTLE_TEST_ONLINE=1` before `tutorial-generate` to execute this chapter.\n");
    }
    if !chapter.summary.trim().is_empty() {
        out.push('\n');
        out.push_str(chapter.summary.trim());
        out.push_str("\n");
    }
    if !chapter.narrative.trim().is_empty() {
        out.push('\n');
        out.push_str(chapter.narrative.trim());
        out.push_str("\n");
    }
    out.push_str("\n## When This Routine Is Useful\n\n");
    for use_case in &chapter.use_cases {
        out.push_str("- ");
        out.push_str(use_case.trim());
        out.push('\n');
    }
    out.push_str("\n## What You Learn\n\n");
    for objective in &chapter.learning_objectives {
        out.push_str("- ");
        out.push_str(objective.trim());
        out.push('\n');
    }
    out.push_str("\n## Concepts and Recurrence\n\n");
    for concept_id in &chapter.concepts {
        let concept = concept_by_id.get(concept_id).ok_or_else(|| {
            format!(
                "Tutorial chapter '{}' references unknown concept id '{}'",
                chapter.id, concept_id
            )
        })?;
        out.push_str("- **");
        out.push_str(concept.name.trim());
        out.push_str("** (`");
        out.push_str(concept.id.trim());
        out.push_str("`): ");
        if concept.description.trim().is_empty() {
            out.push_str("No description provided.");
        } else {
            out.push_str(concept.description.trim());
        }
        out.push('\n');
        let positions = occurrences_by_concept.get(concept_id).ok_or_else(|| {
            format!(
                "No recurrence mapping found for concept id '{}'",
                concept_id
            )
        })?;
        let mut previous = vec![];
        let mut upcoming = vec![];
        for position in positions {
            if *position < chapter_index {
                previous.push(*position);
            } else if *position > chapter_index {
                upcoming.push(*position);
            }
        }
        if previous.is_empty() {
            out.push_str("  - Status: introduced in this chapter.\n");
        } else {
            out.push_str("  - Status: reinforced from ");
            for (idx, prev_position) in previous.iter().enumerate() {
                if idx > 0 {
                    out.push_str(", ");
                }
                let prev = sorted_chapters.get(*prev_position).ok_or_else(|| {
                    format!(
                        "Invalid chapter index '{}' in recurrence map",
                        prev_position
                    )
                })?;
                out.push('[');
                out.push_str(&chapter_link_label(prev));
                out.push_str("](");
                out.push_str(&chapter_link_target(prev, false));
                out.push(')');
            }
            out.push_str(".\n");
        }
        if upcoming.is_empty() {
            out.push_str("  - Reoccurs in: no later chapter.\n");
        } else {
            out.push_str("  - Reoccurs in: ");
            for (idx, next_position) in upcoming.iter().enumerate() {
                if idx > 0 {
                    out.push_str(", ");
                }
                let next = sorted_chapters.get(*next_position).ok_or_else(|| {
                    format!(
                        "Invalid chapter index '{}' in recurrence map",
                        next_position
                    )
                })?;
                out.push('[');
                out.push_str(&chapter_link_label(next));
                out.push_str("](");
                out.push_str(&chapter_link_target(next, false));
                out.push(')');
            }
            out.push_str(".\n");
        }
    }
    out.push_str("\n## GUI First\n\n");
    if chapter.gui_steps.is_empty() {
        out.push_str("- No GUI-first steps are required for this chapter.\n");
    } else {
        for (idx, step) in chapter.gui_steps.iter().enumerate() {
            out.push_str(&format!("{}. ", idx + 1));
            out.push_str(step.trim());
            out.push('\n');
        }
    }
    out.push_str("\n## Command Equivalent (After GUI)\n\n");
    out.push_str("Run the same routine non-interactively once the GUI flow is clear:\n\n");
    out.push_str("```bash\n");
    out.push_str("cargo run --bin gentle_cli -- workflow @");
    out.push_str(&source_path);
    out.push_str("\n");
    out.push_str("cargo run --bin gentle_cli -- shell 'workflow @");
    out.push_str(&source_path);
    out.push_str("'\n");
    out.push_str("```\n");
    out.push_str("\n## Parameters That Matter\n\n");
    if chapter.parameter_notes.is_empty() {
        out.push_str("- This chapter intentionally avoids additional command options; run the canonical workflow unchanged.\n");
    } else {
        for note in &chapter.parameter_notes {
            out.push_str("- `");
            out.push_str(note.parameter.trim());
            out.push_str("`");
            if !note.where_used.trim().is_empty() {
                out.push_str(" (where used: ");
                out.push_str(note.where_used.trim());
                out.push(')');
            }
            out.push('\n');
            if !note.why_it_matters.trim().is_empty() {
                out.push_str("  - Why it matters: ");
                out.push_str(note.why_it_matters.trim());
                out.push('\n');
            }
            if !note.how_to_derive.trim().is_empty() {
                out.push_str("  - How to derive it: ");
                out.push_str(note.how_to_derive.trim());
                out.push('\n');
            }
            if !note.omit_when.trim().is_empty() {
                out.push_str("  - Omit when: ");
                out.push_str(note.omit_when.trim());
                out.push('\n');
            }
        }
    }
    if !chapter.follow_up_commands.is_empty() {
        out.push_str("\n## Follow-up Commands\n\n");
        out.push_str("```bash\n");
        for command in &chapter.follow_up_commands {
            out.push_str(command.trim());
            out.push('\n');
        }
        out.push_str("```\n");
    }
    if !chapter.checkpoints.is_empty() {
        out.push_str("\n## Checkpoints\n\n");
        for checkpoint in &chapter.checkpoints {
            out.push_str("- ");
            out.push_str(checkpoint.trim());
            out.push('\n');
        }
    }
    out.push_str("\n## Retained Outputs\n\n");
    if retained_artifacts.is_empty() {
        out.push_str("- None for this chapter.\n");
    } else {
        for retained in retained_artifacts {
            out.push_str("- [`");
            out.push_str(retained);
            out.push_str("`](../");
            out.push_str(retained);
            out.push_str(")\n");
        }
    }
    out.push_str("\n## Canonical Source\n\n");
    out.push_str("- Workflow file: `");
    out.push_str(&source_path);
    out.push_str("`\n");
    out.push_str("- Inspect this JSON file directly when you need full option-level detail.\n");
    Ok(out)
}

fn render_tutorial_index(
    manifest: &TutorialManifest,
    chapters: &[TutorialGenerationChapter],
    online_enabled: bool,
    sorted_chapters: &[TutorialChapter],
    occurrences_by_concept: &HashMap<String, Vec<usize>>,
) -> String {
    let mut out = String::new();
    out.push_str("# GENtle Tutorial (Generated)\n\n");
    if !manifest.description.trim().is_empty() {
        out.push_str(manifest.description.trim());
        out.push_str("\n\n");
    }
    out.push_str("This folder is generated from:\n");
    out.push_str("- `docs/tutorial/sources/catalog_meta.json`\n");
    out.push_str("- `docs/tutorial/sources/*.json`\n");
    out.push_str("- `docs/examples/workflows/*.json`\n\n");
    out.push_str("Generated runtime manifest:\n");
    out.push_str("- `docs/tutorial/manifest.json`\n\n");
    out.push_str("Regenerate with:\n\n");
    out.push_str("```bash\n");
    out.push_str("cargo run --bin gentle_examples_docs -- tutorial-generate\n");
    out.push_str("```\n\n");
    out.push_str("Validate committed generated output:\n\n");
    out.push_str("```bash\n");
    out.push_str("cargo run --bin gentle_examples_docs -- tutorial-check\n");
    out.push_str("```\n\n");
    out.push_str("Recommended progression:\n");
    out.push_str(
        "- Start chapters in the GUI to understand biological intent and visual context.\n",
    );
    out.push_str("- Then use the command equivalents for repeatable runs and automation.\n\n");
    out.push_str("Online execution was ");
    out.push_str(if online_enabled {
        "enabled"
    } else {
        "disabled"
    });
    out.push_str(" (`GENTLE_TEST_ONLINE=");
    out.push_str(if online_enabled { "1`" } else { "0`" });
    out.push_str(" during generation).\n\n");
    out.push_str("## Chapters\n\n");
    let mut ordered = chapters.to_vec();
    ordered.sort_by(|a, b| a.order.cmp(&b.order).then_with(|| a.id.cmp(&b.id)));
    for chapter in &ordered {
        let chapter_file = chapter_markdown_filename(chapter.order, &chapter.id);
        out.push_str("- ");
        out.push_str(&chapter.order.to_string());
        out.push_str(". [");
        out.push_str(&chapter.title);
        out.push_str("](./chapters/");
        out.push_str(&chapter_file);
        out.push_str(") - `");
        out.push_str(&chapter.tier);
        out.push_str("` - example `");
        out.push_str(&chapter.example_id);
        out.push_str("` - executed `");
        out.push_str(if chapter.executed { "yes" } else { "no" });
        out.push_str("`\n");
    }
    out.push_str("\n## Concepts and Where They Recur\n\n");
    for concept in &manifest.concepts {
        out.push_str("- **");
        out.push_str(concept.name.trim());
        out.push_str("** (`");
        out.push_str(concept.id.trim());
        out.push_str("`): ");
        if concept.description.trim().is_empty() {
            out.push_str("No description provided.");
        } else {
            out.push_str(concept.description.trim());
        }
        out.push('\n');
        if let Some(indices) = occurrences_by_concept.get(&concept.id) {
            out.push_str("  - Appears in: ");
            for (idx, chapter_index) in indices.iter().enumerate() {
                if idx > 0 {
                    out.push_str(", ");
                }
                if let Some(chapter) = sorted_chapters.get(*chapter_index) {
                    out.push('[');
                    out.push_str(&chapter_link_label(chapter));
                    out.push_str("](");
                    out.push_str(&chapter_link_target(chapter, true));
                    out.push(')');
                }
            }
            out.push_str(".\n");
        }
    }
    out.push_str("\n## Source Summary\n\n");
    out.push_str("- Tutorial schema: `");
    out.push_str(&manifest.schema);
    out.push_str("`\n");
    out.push_str("- Chapter count: `");
    out.push_str(&manifest.chapters.len().to_string());
    out.push_str("`\n");
    out.push_str("- Generation report: [`report.json`](./report.json)\n");
    out
}

fn collect_regular_files(root: &Path) -> Result<Vec<PathBuf>, String> {
    let mut files = vec![];
    if !root.exists() {
        return Ok(files);
    }
    let mut stack = vec![root.to_path_buf()];
    while let Some(dir) = stack.pop() {
        let mut entries = vec![];
        for entry in fs::read_dir(&dir)
            .map_err(|e| format!("Could not read '{}': {e}", display_path(&dir)))?
        {
            let entry = entry.map_err(|e| format!("Could not read directory entry: {e}"))?;
            entries.push(entry.path());
        }
        entries.sort();
        for path in entries {
            if path.is_dir() {
                stack.push(path);
            } else if path.is_file() {
                if path
                    .file_name()
                    .map(|name| name == ".DS_Store")
                    .unwrap_or(false)
                {
                    continue;
                }
                files.push(path);
            }
        }
    }
    files.sort();
    Ok(files)
}

fn sha1_file(path: &Path) -> Result<String, String> {
    let mut file = fs::File::open(path)
        .map_err(|e| format!("Could not open '{}' for checksum: {e}", display_path(path)))?;
    let mut hasher = Sha1::new();
    let mut buffer = [0u8; 8192];
    loop {
        let read = file
            .read(&mut buffer)
            .map_err(|e| format!("Could not read '{}' for checksum: {e}", display_path(path)))?;
        if read == 0 {
            break;
        }
        hasher.update(&buffer[..read]);
    }
    Ok(format!("{:x}", hasher.finalize()))
}

fn checksum_map(root: &Path, exclude_names: &[&str]) -> Result<BTreeMap<String, String>, String> {
    let mut checksums = BTreeMap::new();
    let files = collect_regular_files(root)?;
    for file in files {
        let rel = file
            .strip_prefix(root)
            .map_err(|e| format!("Could not compute relative path: {e}"))?;
        let rel_str = display_path(rel);
        let excluded = exclude_names
            .iter()
            .any(|name| rel.file_name().map(|v| v == *name).unwrap_or(false));
        if excluded {
            continue;
        }
        checksums.insert(rel_str, sha1_file(&file)?);
    }
    Ok(checksums)
}

fn directory_bytes_map(root: &Path) -> Result<BTreeMap<String, Vec<u8>>, String> {
    let mut map = BTreeMap::new();
    for file in collect_regular_files(root)? {
        let rel = file
            .strip_prefix(root)
            .map_err(|e| format!("Could not compute relative path: {e}"))?;
        let rel_str = display_path(rel);
        let bytes = fs::read(&file)
            .map_err(|e| format!("Could not read '{}': {e}", display_path(&file)))?;
        map.insert(rel_str, bytes);
    }
    Ok(map)
}

pub fn generate_tutorial_docs(
    source_dir: &Path,
    manifest_path: &Path,
    output_dir: &Path,
    repo_root: &Path,
) -> Result<TutorialGenerationReport, String> {
    let examples = load_workflow_examples(source_dir)?;
    let example_by_id = example_lookup(&examples);
    let manifest = load_tutorial_manifest(manifest_path)?;
    validate_tutorial_manifest_against_examples(&manifest, &examples)?;
    let concept_by_id = tutorial_concept_lookup(&manifest.concepts)?;
    let sorted_chapters = sorted_manifest_chapters(&manifest);
    let occurrences_by_concept = concept_occurrences(&sorted_chapters);
    if output_dir.exists() {
        fs::remove_dir_all(output_dir)
            .map_err(|e| format!("Could not clean '{}': {e}", display_path(output_dir)))?;
    }
    fs::create_dir_all(output_dir)
        .map_err(|e| format!("Could not create '{}': {e}", display_path(output_dir)))?;
    let chapters_dir = output_dir.join("chapters");
    fs::create_dir_all(&chapters_dir)
        .map_err(|e| format!("Could not create '{}': {e}", display_path(&chapters_dir)))?;
    let artifacts_dir = output_dir.join("artifacts");
    fs::create_dir_all(&artifacts_dir)
        .map_err(|e| format!("Could not create '{}': {e}", display_path(&artifacts_dir)))?;

    let online_enabled = online_example_tests_enabled();
    let mut chapter_reports: Vec<TutorialGenerationChapter> = vec![];

    for (chapter_index, chapter) in sorted_chapters.iter().enumerate() {
        let loaded = example_by_id
            .get(&chapter.example_id)
            .ok_or_else(|| {
                format!(
                    "Tutorial chapter '{}' references unknown example '{}'",
                    chapter.id, chapter.example_id
                )
            })?
            .clone();
        let executed = chapter.tier.should_execute(online_enabled);
        let run_dir =
            TempDir::new().map_err(|e| format!("Could not create temp run directory: {e}"))?;
        if executed {
            validate_example_required_files(&loaded.example, repo_root)?;
            run_example_workflow_in_dir(&loaded.example, repo_root, run_dir.path())?;
        }
        let mut retained_artifacts = vec![];
        if executed {
            for retain_output in &chapter.retain_outputs {
                let retained =
                    copy_retained_output(chapter, run_dir.path(), output_dir, retain_output)?;
                retained_artifacts.push(retained);
            }
        }
        retained_artifacts.sort();
        let chapter_markdown = render_tutorial_chapter_markdown(
            chapter,
            &loaded,
            executed,
            &retained_artifacts,
            online_enabled,
            chapter_index,
            &sorted_chapters,
            &concept_by_id,
            &occurrences_by_concept,
        )?;
        let chapter_file = markdown_path_for_chapter(chapter);
        let chapter_out_path = chapters_dir.join(&chapter_file);
        fs::write(&chapter_out_path, chapter_markdown)
            .map_err(|e| format!("Could not write '{}': {e}", display_path(&chapter_out_path)))?;
        chapter_reports.push(TutorialGenerationChapter {
            id: chapter.id.clone(),
            order: chapter.order,
            title: chapter.title.clone(),
            tier: chapter.tier.as_str().to_string(),
            example_id: chapter.example_id.clone(),
            example_source: display_path(&loaded.path),
            concepts: chapter.concepts.clone(),
            executed,
            retained_artifacts,
        });
    }

    let readme = render_tutorial_index(
        &manifest,
        &chapter_reports,
        online_enabled,
        &sorted_chapters,
        &occurrences_by_concept,
    );
    let readme_path = output_dir.join("README.md");
    fs::write(&readme_path, readme)
        .map_err(|e| format!("Could not write '{}': {e}", display_path(&readme_path)))?;

    let file_checksums = checksum_map(output_dir, &["report.json"])?;
    let mut generated_files: Vec<String> = file_checksums.keys().cloned().collect();
    generated_files.push("report.json".to_string());
    generated_files.sort();
    let report = TutorialGenerationReport {
        schema: TUTORIAL_GENERATION_REPORT_SCHEMA.to_string(),
        manifest: display_path(manifest_path),
        source_examples: display_path(source_dir),
        chapter_count: chapter_reports.len(),
        online_enabled,
        generated_files,
        file_checksums,
        chapters: chapter_reports,
    };
    let report_path = output_dir.join("report.json");
    let report_json = serde_json::to_string_pretty(&report)
        .map_err(|e| format!("Could not serialize tutorial report: {e}"))?;
    fs::write(&report_path, report_json)
        .map_err(|e| format!("Could not write '{}': {e}", display_path(&report_path)))?;
    Ok(report)
}

pub fn check_tutorial_generated(
    source_dir: &Path,
    manifest_path: &Path,
    expected_output_dir: &Path,
    repo_root: &Path,
) -> Result<TutorialGenerationReport, String> {
    if !expected_output_dir.exists() {
        return Err(format!(
            "Expected generated tutorial directory '{}' does not exist. Run `tutorial-generate` first.",
            display_path(expected_output_dir)
        ));
    }
    let temp = TempDir::new().map_err(|e| format!("Could not create temp directory: {e}"))?;
    let generated_dir = temp.path().join("generated");
    let report = generate_tutorial_docs(source_dir, manifest_path, &generated_dir, repo_root)?;
    let expected = directory_bytes_map(expected_output_dir)?;
    let actual = directory_bytes_map(&generated_dir)?;
    if expected.len() != actual.len() {
        return Err(format!(
            "Tutorial generated file count mismatch: expected {}, got {}",
            expected.len(),
            actual.len()
        ));
    }
    for (path, expected_bytes) in &expected {
        let actual_bytes = actual.get(path).ok_or_else(|| {
            format!(
                "Tutorial generated output is missing file '{}' in check run",
                path
            )
        })?;
        if expected_bytes != actual_bytes {
            return Err(format!(
                "Tutorial generated file '{}' differs from committed version",
                path
            ));
        }
    }
    for path in actual.keys() {
        if !expected.contains_key(path) {
            return Err(format!(
                "Tutorial generated output contains unexpected file '{}'",
                path
            ));
        }
    }
    Ok(report)
}

fn test_mode_doc(test_mode: ExampleTestMode) -> &'static str {
    match test_mode {
        ExampleTestMode::Always => "always (included in default test runs)",
        ExampleTestMode::Skip => "skip (validated for syntax only)",
        ExampleTestMode::Online => {
            "online (run only when GENTLE_TEST_ONLINE=1 with working internet)"
        }
    }
}

pub fn render_example_markdown(loaded: &LoadedWorkflowExample) -> Result<String, String> {
    let workflow_json = serde_json::to_string_pretty(&loaded.example.workflow).map_err(|e| {
        format!(
            "Could not serialize workflow example '{}': {e}",
            loaded.example.id
        )
    })?;
    let source_path = display_path(&loaded.path);
    let mut out = String::new();
    out.push_str("# ");
    out.push_str(&loaded.example.title);
    out.push('\n');
    out.push('\n');
    out.push_str("- Example id: `");
    out.push_str(&loaded.example.id);
    out.push_str("`\n");
    out.push_str("- Source file: `");
    out.push_str(&source_path);
    out.push_str("`\n");
    out.push_str("- Test mode: `");
    out.push_str(loaded.example.test_mode.as_str());
    out.push_str("` (");
    out.push_str(test_mode_doc(loaded.example.test_mode));
    out.push_str(")\n");
    if !loaded.example.required_files.is_empty() {
        out.push_str("- Required files:\n");
        for required in &loaded.example.required_files {
            out.push_str("  - `");
            out.push_str(required);
            out.push_str("`\n");
        }
    }
    if !loaded.example.summary.trim().is_empty() {
        out.push('\n');
        out.push_str(&loaded.example.summary);
        out.push('\n');
    }
    out.push('\n');
    out.push_str("## Canonical Workflow JSON\n\n");
    out.push_str("```json\n");
    out.push_str(&workflow_json);
    out.push_str("\n```\n\n");
    out.push_str("## CLI\n\n");
    out.push_str("```bash\n");
    out.push_str("cargo run --bin gentle_cli -- workflow @");
    out.push_str(&source_path);
    out.push_str("\n```\n\n");
    out.push_str("## Shared Shell\n\n");
    out.push_str("```bash\n");
    out.push_str("cargo run --bin gentle_cli -- shell 'workflow @");
    out.push_str(&source_path);
    out.push_str("'\n```\n\n");
    out.push_str("## JavaScript Shell\n\n");
    out.push_str("```javascript\n");
    out.push_str("const state = { sequences: {}, metadata: {}, display: {}, lineage: {} };\n");
    out.push_str("const workflow = ");
    out.push_str(&workflow_json);
    out.push_str(";\n");
    out.push_str("const run = apply_workflow(state, workflow);\n");
    out.push_str("console.log(run.results.map(r => r.op_id));\n");
    out.push_str("```\n\n");
    out.push_str("## Lua Shell\n\n");
    out.push_str("```lua\n");
    out.push_str("local project = { sequences = {}, metadata = {}, display = {}, lineage = {} }\n");
    out.push_str("local workflow_json = [[\n");
    out.push_str(&workflow_json);
    out.push_str("\n]]\n");
    out.push_str("local run = apply_workflow(project, workflow_json)\n");
    out.push_str("for i, result in ipairs(run.results) do print(i, result.op_id) end\n");
    out.push_str("```\n");
    Ok(out)
}

pub fn render_examples_index(examples: &[LoadedWorkflowExample]) -> String {
    let mut out = String::new();
    out.push_str("# GENtle Protocol-First Workflow Examples\n\n");
    out.push_str("This folder is generated from canonical workflow examples in ");
    out.push_str("`docs/examples/workflows`.\n\n");
    out.push_str("Regenerate with:\n\n");
    out.push_str("```bash\n");
    out.push_str("cargo run --bin gentle_examples_docs -- generate\n");
    out.push_str("```\n\n");
    out.push_str("Validation-only mode:\n\n");
    out.push_str("```bash\n");
    out.push_str("cargo run --bin gentle_examples_docs -- --check\n");
    out.push_str("```\n\n");
    out.push_str("Online test opt-in:\n\n");
    out.push_str("```bash\n");
    out.push_str("GENTLE_TEST_ONLINE=1 cargo test workflow_examples -- --test-threads=1\n");
    out.push_str("```\n\n");
    out.push_str("## Examples\n\n");
    for loaded in examples {
        let stem = markdown_file_stem(&loaded.example.id);
        out.push_str("- [");
        out.push_str(&loaded.example.id);
        out.push_str("](./");
        out.push_str(&stem);
        out.push_str(".md)");
        out.push_str(" - ");
        out.push_str(loaded.example.test_mode.as_str());
        if !loaded.example.summary.trim().is_empty() {
            out.push_str(" - ");
            out.push_str(loaded.example.summary.trim());
        }
        out.push('\n');
    }
    out
}

pub fn generate_workflow_example_docs(
    source_dir: &Path,
    output_dir: &Path,
) -> Result<WorkflowExampleDocReport, String> {
    let examples = load_workflow_examples(source_dir)?;
    fs::create_dir_all(output_dir)
        .map_err(|e| format!("Could not create '{}': {e}", display_path(output_dir)))?;
    let mut generated_files = vec![];
    for loaded in &examples {
        let stem = markdown_file_stem(&loaded.example.id);
        let output_path = output_dir.join(format!("{stem}.md"));
        let markdown = render_example_markdown(loaded)?;
        fs::write(&output_path, markdown)
            .map_err(|e| format!("Could not write '{}': {e}", display_path(&output_path)))?;
        generated_files.push(display_path(&output_path));
    }
    let index_path = output_dir.join("README.md");
    fs::write(&index_path, render_examples_index(&examples))
        .map_err(|e| format!("Could not write '{}': {e}", display_path(&index_path)))?;
    generated_files.push(display_path(&index_path));
    Ok(WorkflowExampleDocReport {
        source_dir: display_path(source_dir),
        output_dir: display_path(output_dir),
        example_count: examples.len(),
        generated_files,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;
    use std::sync::Mutex;

    static TEST_ENV_LOCK: Mutex<()> = Mutex::new(());

    struct EnvVarGuard {
        key: &'static str,
        previous: Option<String>,
    }

    impl EnvVarGuard {
        fn set(key: &'static str, value: &str) -> Self {
            let previous = std::env::var(key).ok();
            // SAFETY: Test-only mutation is serialized by TEST_ENV_LOCK.
            unsafe { std::env::set_var(key, value) };
            Self { key, previous }
        }
    }

    impl Drop for EnvVarGuard {
        fn drop(&mut self) {
            if let Some(value) = &self.previous {
                // SAFETY: Test-only mutation is serialized by TEST_ENV_LOCK.
                unsafe { std::env::set_var(self.key, value) };
            } else {
                // SAFETY: Test-only mutation is serialized by TEST_ENV_LOCK.
                unsafe { std::env::remove_var(self.key) };
            }
        }
    }

    fn example_dir() -> PathBuf {
        PathBuf::from(DEFAULT_WORKFLOW_EXAMPLE_DIR)
    }

    fn tutorial_manifest_path() -> PathBuf {
        PathBuf::from(DEFAULT_TUTORIAL_MANIFEST_PATH)
    }

    fn tutorial_catalog_path() -> PathBuf {
        PathBuf::from(DEFAULT_TUTORIAL_CATALOG_PATH)
    }

    fn tutorial_catalog_meta_path() -> PathBuf {
        PathBuf::from(DEFAULT_TUTORIAL_CATALOG_META_PATH)
    }

    fn tutorial_source_dir() -> PathBuf {
        PathBuf::from(DEFAULT_TUTORIAL_SOURCE_DIR)
    }

    fn tutorial_output_dir() -> PathBuf {
        PathBuf::from(DEFAULT_TUTORIAL_OUTPUT_DIR)
    }

    fn markdown_image_targets(markdown: &str) -> Vec<String> {
        let mut targets = Vec::new();
        for line in markdown.lines() {
            let mut remaining = line;
            loop {
                let Some(image_start) = remaining.find("![") else {
                    break;
                };
                remaining = &remaining[image_start + 2..];
                let Some(target_start) = remaining.find("](") else {
                    break;
                };
                let after_target_start = &remaining[target_start + 2..];
                let Some(target_end) = after_target_start.find(')') else {
                    break;
                };
                let target = after_target_start[..target_end].trim();
                if !target.is_empty() {
                    let path_only = target.split_whitespace().next().unwrap_or(target).trim();
                    if !path_only.is_empty() {
                        targets.push(path_only.to_string());
                    }
                }
                remaining = &after_target_start[target_end + 1..];
            }
        }
        targets
    }

    #[test]
    fn workflow_examples_parse_and_render() {
        let examples = load_workflow_examples(&example_dir()).expect("load workflow examples");
        assert!(
            !examples.is_empty(),
            "Expected at least one workflow example"
        );
        let markdown = render_example_markdown(&examples[0]).expect("render markdown");
        assert!(markdown.contains("Canonical Workflow JSON"));
        assert!(markdown.contains("cargo run --bin gentle_cli -- workflow @"));
    }

    #[test]
    fn workflow_examples_always_mode_executes() {
        let examples = load_workflow_examples(&example_dir()).expect("load workflow examples");
        let repo_root = Path::new(".");
        let mut always_count = 0usize;
        for loaded in &examples {
            validate_example_required_files(&loaded.example, repo_root)
                .expect("required files must be present");
            if loaded.example.test_mode == ExampleTestMode::Always {
                always_count += 1;
                run_example_workflow(&loaded.example).expect("always example should execute");
            }
        }
        assert!(always_count > 0, "Expected at least one always example");
    }

    #[test]
    fn workflow_examples_tp73_isoform_protein_gel_writes_svg() {
        let examples = load_workflow_examples(&example_dir()).expect("load workflow examples");
        let loaded = examples
            .iter()
            .find(|loaded| loaded.example.id == "tp73_isoform_protein_gel_offline")
            .expect("tp73 isoform protein gel example should exist");
        let run_dir = TempDir::new().expect("temp run dir");
        run_example_workflow_in_dir(&loaded.example, Path::new("."), run_dir.path())
            .expect("tp73 protein-gel workflow should execute");
        let svg_path = run_dir.path().join("exports/tp73_isoform_protein_gel.svg");
        assert!(svg_path.exists());
        let svg = fs::read_to_string(&svg_path).expect("read tp73 protein gel svg");
        assert!(svg.contains("Protein Gel Preview"));
        assert!(svg.contains("Protein Ladder 10-100 kDa"));
        assert!(svg.contains("Selection notes"));
    }

    #[test]
    fn workflow_examples_tp73_isoform_protein_2d_gel_writes_svg() {
        let examples = load_workflow_examples(&example_dir()).expect("load workflow examples");
        let loaded = examples
            .iter()
            .find(|loaded| loaded.example.id == "tp73_isoform_protein_2d_gel_offline")
            .expect("tp73 isoform protein 2D gel example should exist");
        let run_dir = TempDir::new().expect("temp run dir");
        run_example_workflow_in_dir(&loaded.example, Path::new("."), run_dir.path())
            .expect("tp73 protein 2D gel workflow should execute");
        let svg_path = run_dir.path().join("exports/tp73_isoform_protein_2d_gel.svg");
        assert!(svg_path.exists());
        let svg = fs::read_to_string(&svg_path).expect("read tp73 protein 2D gel svg");
        assert!(svg.contains("Protein 2D Gel Preview"));
        assert!(svg.contains("Protein Ladder 10-100 kDa"));
        assert!(svg.contains("Selection notes"));
        assert!(svg.contains("Spot details"));
    }

    #[test]
    fn workflow_examples_online_mode_opt_in() {
        let examples = load_workflow_examples(&example_dir()).expect("load workflow examples");
        if !online_example_tests_enabled() {
            return;
        }
        let repo_root = Path::new(".");
        for loaded in &examples {
            if loaded.example.test_mode != ExampleTestMode::Online {
                continue;
            }
            validate_example_required_files(&loaded.example, repo_root)
                .expect("required files must be present");
            run_example_workflow(&loaded.example).expect("online example should execute");
        }
    }

    #[test]
    fn online_example_tests_skip_remote_override_wins() {
        let _env_lock = TEST_ENV_LOCK.lock().expect("env lock");
        let _online = EnvVarGuard::set(ONLINE_EXAMPLE_TEST_ENV, "1");
        let _skip = EnvVarGuard::set(SKIP_REMOTE_TESTS_ENV, "1");
        assert!(
            !online_example_tests_enabled(),
            "GENTLE_SKIP_REMOTE_TESTS=1 should disable online tests even when GENTLE_TEST_ONLINE=1"
        );
    }

    #[test]
    fn tutorial_manifest_parses_and_validates() {
        let manifest =
            load_tutorial_manifest(&tutorial_manifest_path()).expect("load tutorial manifest");
        let examples = load_workflow_examples(&example_dir()).expect("load workflow examples");
        validate_tutorial_manifest_against_examples(&manifest, &examples)
            .expect("manifest validates against examples");
        assert!(
            !manifest.chapters.is_empty(),
            "Expected at least one tutorial chapter"
        );
    }

    #[test]
    fn tutorial_catalog_parses_and_paths_exist() {
        let catalog =
            load_tutorial_catalog(&tutorial_catalog_path()).expect("load tutorial catalog");
        assert!(
            !catalog.entries.is_empty(),
            "Expected at least one tutorial catalog entry"
        );
        for entry in &catalog.entries {
            let path = PathBuf::from(&entry.path);
            assert!(
                path.exists(),
                "tutorial catalog entry '{}' points to missing path '{}'",
                entry.id,
                display_path(&path)
            );
        }
        let generated_manifest = PathBuf::from(&catalog.generated_runtime.manifest_path);
        assert!(
            generated_manifest.exists(),
            "catalog manifest path should exist: '{}'",
            display_path(&generated_manifest)
        );
    }

    #[test]
    fn tutorial_source_units_generate_committed_catalog() {
        let generated = generate_tutorial_catalog_from_sources(
            &tutorial_catalog_meta_path(),
            &tutorial_source_dir(),
        )
        .expect("generate tutorial catalog from sources");
        let committed = load_tutorial_catalog(&tutorial_catalog_path())
            .expect("load committed tutorial catalog");
        let generated_json =
            serde_json::to_string_pretty(&generated).expect("serialize generated tutorial catalog");
        let committed_json =
            serde_json::to_string_pretty(&committed).expect("serialize committed tutorial catalog");
        assert_eq!(
            generated_json, committed_json,
            "committed tutorial catalog should match generated source units"
        );
    }

    #[test]
    fn tutorial_catalog_check_passes_on_committed_tree() {
        check_tutorial_catalog_generated(
            &tutorial_catalog_meta_path(),
            &tutorial_source_dir(),
            &tutorial_catalog_path(),
        )
        .expect("committed tutorial catalog should match source units");
    }

    #[test]
    fn tutorial_source_units_generate_committed_manifest() {
        let generated = generate_tutorial_manifest_from_sources(
            &tutorial_catalog_meta_path(),
            &tutorial_source_dir(),
        )
        .expect("generate tutorial manifest from sources");
        let committed = load_tutorial_manifest(&tutorial_manifest_path())
            .expect("load committed tutorial manifest");
        let generated_json = serde_json::to_string_pretty(&generated)
            .expect("serialize generated tutorial manifest");
        let committed_json = serde_json::to_string_pretty(&committed)
            .expect("serialize committed tutorial manifest");
        assert_eq!(
            generated_json, committed_json,
            "committed tutorial manifest should match generated source units"
        );
    }

    #[test]
    fn tutorial_manifest_check_passes_on_committed_tree() {
        check_tutorial_manifest_generated(
            &tutorial_catalog_meta_path(),
            &tutorial_source_dir(),
            &tutorial_manifest_path(),
        )
        .expect("committed tutorial manifest should match source units");
    }

    #[test]
    fn tutorial_core_chapters_map_to_always_examples() {
        let manifest =
            load_tutorial_manifest(&tutorial_manifest_path()).expect("load tutorial manifest");
        let examples = load_workflow_examples(&example_dir()).expect("load workflow examples");
        let by_id = example_lookup(&examples);
        let mut core_count = 0usize;
        for chapter in &manifest.chapters {
            if chapter.tier != TutorialTier::Core {
                continue;
            }
            core_count += 1;
            let loaded = by_id
                .get(&chapter.example_id)
                .expect("core chapter example should exist");
            assert_eq!(
                loaded.example.test_mode,
                ExampleTestMode::Always,
                "core chapter '{}' should map to always-mode example",
                chapter.id
            );
        }
        assert!(
            core_count > 0,
            "Expected at least one core tutorial chapter"
        );
    }

    #[test]
    fn tutorial_generate_is_deterministic() {
        let source = example_dir();
        let manifest = tutorial_manifest_path();
        let repo_root = Path::new(".");
        let first_dir = TempDir::new().expect("create first temp dir");
        let second_dir = TempDir::new().expect("create second temp dir");
        let first_out = first_dir.path().join("generated");
        let second_out = second_dir.path().join("generated");
        generate_tutorial_docs(&source, &manifest, &first_out, repo_root)
            .expect("first generation should pass");
        generate_tutorial_docs(&source, &manifest, &second_out, repo_root)
            .expect("second generation should pass");
        let first = directory_bytes_map(&first_out).expect("collect first outputs");
        let second = directory_bytes_map(&second_out).expect("collect second outputs");
        assert_eq!(
            first, second,
            "tutorial generation output should be deterministic"
        );
    }

    #[test]
    fn tutorial_generated_chapter_includes_narrative_concepts_and_objectives() {
        let source = example_dir();
        let manifest = tutorial_manifest_path();
        let repo_root = Path::new(".");
        let out_dir = TempDir::new().expect("create temp dir");
        let generated = out_dir.path().join("generated");
        generate_tutorial_docs(&source, &manifest, &generated, repo_root)
            .expect("tutorial generation should pass");
        let chapter = generated.join("chapters/01_load_branch_reverse_complement_pgex_fasta.md");
        let markdown = std::fs::read_to_string(&chapter).expect("read generated chapter markdown");
        assert!(markdown.contains("## What You Learn"));
        assert!(markdown.contains("## Concepts and Recurrence"));
        assert!(markdown.contains("## GUI First"));
        assert!(markdown.contains("## Parameters That Matter"));
    }

    #[test]
    fn tutorial_tp73_cdna_genomic_markdown_image_links_exist() {
        let tutorial_path = PathBuf::from("docs/tutorial/two_sequence_dotplot_gui.md");
        let tutorial_text =
            fs::read_to_string(&tutorial_path).expect("read TP73 cDNA/genomic tutorial");
        let targets = markdown_image_targets(&tutorial_text);
        assert!(
            !targets.is_empty(),
            "expected markdown image references in {}",
            display_path(&tutorial_path)
        );
        let parent = tutorial_path
            .parent()
            .expect("tutorial path should have parent");
        for target in targets {
            if target.starts_with("http://")
                || target.starts_with("https://")
                || target.starts_with("data:")
            {
                continue;
            }
            let resolved = if Path::new(&target).is_absolute() {
                PathBuf::from(&target)
            } else {
                parent.join(&target)
            };
            assert!(
                resolved.exists(),
                "markdown image target '{}' in '{}' does not exist (resolved: '{}')",
                target,
                display_path(&tutorial_path),
                display_path(&resolved)
            );
        }
    }

    #[test]
    fn rewrite_example_paths_handles_isoform_panel_io() {
        let example = WorkflowExample {
            schema: WORKFLOW_EXAMPLE_SCHEMA.to_string(),
            id: "isoform_path_rewrite_test".to_string(),
            title: "isoform path rewrite test".to_string(),
            summary: String::new(),
            test_mode: ExampleTestMode::Skip,
            required_files: vec!["assets/panels/tp53_isoforms_v1.json".to_string()],
            tags: vec![],
            workflow: Workflow {
                run_id: "isoform_path_rewrite_test".to_string(),
                ops: vec![
                    Operation::ImportIsoformPanel {
                        seq_id: "seq_a".to_string(),
                        panel_path: "assets/panels/tp53_isoforms_v1.json".to_string(),
                        panel_id: Some("tp53_isoforms_v1".to_string()),
                        strict: false,
                    },
                    Operation::RenderIsoformArchitectureSvg {
                        seq_id: "seq_a".to_string(),
                        panel_id: "tp53_isoforms_v1".to_string(),
                        path: "exports/tp53_isoform.svg".to_string(),
                    },
                ],
            },
        };
        let repo_root = std::env::current_dir().expect("cwd");
        let run_dir = TempDir::new().expect("temp run dir");
        let rewritten =
            rewrite_example_paths_for_execution(&example, repo_root.as_path(), run_dir.path())
                .expect("rewrite should succeed");
        match &rewritten.workflow.ops[0] {
            Operation::ImportIsoformPanel { panel_path, .. } => {
                assert!(
                    Path::new(panel_path).is_absolute(),
                    "panel path should be rewritten to absolute path"
                );
            }
            other => panic!("unexpected operation: {other:?}"),
        }
        match &rewritten.workflow.ops[1] {
            Operation::RenderIsoformArchitectureSvg { path, .. } => {
                assert!(
                    path.starts_with(&display_path(run_dir.path())),
                    "render output path should be rewritten into run dir"
                );
            }
            other => panic!("unexpected operation: {other:?}"),
        }
    }

    #[test]
    fn tutorial_check_passes_on_committed_tree() {
        check_tutorial_generated(
            &example_dir(),
            &tutorial_manifest_path(),
            &tutorial_output_dir(),
            Path::new("."),
        )
        .expect("committed tutorial output should match generator");
    }

    #[test]
    fn tutorial_core_executes_runtime() {
        let manifest =
            load_tutorial_manifest(&tutorial_manifest_path()).expect("load tutorial manifest");
        let examples = load_workflow_examples(&example_dir()).expect("load workflow examples");
        let by_id = example_lookup(&examples);
        let repo_root = Path::new(".");
        let mut executed = 0usize;
        for chapter in &manifest.chapters {
            if chapter.tier != TutorialTier::Core {
                continue;
            }
            let loaded = by_id
                .get(&chapter.example_id)
                .expect("core chapter example should exist");
            validate_example_required_files(&loaded.example, repo_root)
                .expect("required files should exist");
            run_example_workflow(&loaded.example).expect("core example should execute");
            executed += 1;
        }
        assert!(
            executed > 0,
            "Expected at least one core chapter to execute"
        );
    }

    #[test]
    fn run_example_workflow_for_project_state_returns_sequences() {
        let examples = load_workflow_examples(&example_dir()).expect("load workflow examples");
        let loaded = examples
            .iter()
            .find(|loaded| loaded.example.id == "load_branch_reverse_complement_pgex_fasta")
            .expect("example should exist");
        let run_dir = TempDir::new().expect("temp run dir");
        let state =
            run_example_workflow_for_project_state(&loaded.example, Path::new("."), run_dir.path())
                .expect("workflow should produce project state");
        assert!(state.sequences.contains_key("pgex_fasta"));
        assert!(state.sequences.contains_key("pgex_fasta_branch"));
        assert!(state.sequences.contains_key("pgex_fasta_branch_rc"));
    }

    #[test]
    fn gibson_arrangements_baseline_example_materializes_arrangement_ready_state() {
        let examples = load_workflow_examples(&example_dir()).expect("load workflow examples");
        let loaded = examples
            .iter()
            .find(|loaded| loaded.example.id == "gibson_arrangements_baseline")
            .expect("arrangements tutorial example should exist");
        let run_dir = TempDir::new().expect("temp run dir");
        let state =
            run_example_workflow_for_project_state(&loaded.example, Path::new("."), run_dir.path())
                .expect("workflow should produce project state");

        assert!(state.sequences.contains_key("gibson_destination_pgex"));
        assert!(state.sequences.contains_key("gibson_insert_demo"));
        assert!(
            state
                .sequences
                .contains_key("gibson_destination_pgex_with_gibson_insert_demo")
        );

        let arrangement = state
            .container_state
            .arrangements
            .values()
            .next()
            .expect("arrangement baseline should create one arrangement");
        assert_eq!(arrangement.lane_container_ids.len(), 3);

        let lane_members = arrangement
            .lane_container_ids
            .iter()
            .map(|container_id| {
                state
                    .container_state
                    .containers
                    .get(container_id)
                    .map(|container| container.members.clone())
                    .unwrap_or_default()
            })
            .collect::<Vec<_>>();
        assert_eq!(
            lane_members,
            vec![
                vec!["gibson_destination_pgex".to_string()],
                vec!["gibson_insert_demo".to_string()],
                vec!["gibson_destination_pgex_with_gibson_insert_demo".to_string()],
            ]
        );
    }

    #[test]
    fn oe_substitution_example_materializes_staged_pcr_artifacts() {
        let examples = load_workflow_examples(&example_dir()).expect("load workflow examples");
        let loaded = examples
            .iter()
            .find(|loaded| loaded.example.id == "pcr_overlap_extension_substitution_offline")
            .expect("example should exist");
        let run_dir = TempDir::new().expect("temp run dir");
        let state =
            run_example_workflow_for_project_state(&loaded.example, Path::new("."), run_dir.path())
                .expect("workflow should produce project state");
        assert!(state.sequences.contains_key("pgex_linear"));
        let mut has_outer = false;
        let mut has_stage1_left = false;
        let mut has_stage1_right = false;
        let mut has_mutant = false;
        for seq_id in state.sequences.keys() {
            has_outer |= seq_id.contains("pgex_oe_sub_outer_fwd");
            has_stage1_left |= seq_id.contains("pgex_oe_sub_stage1_left");
            has_stage1_right |= seq_id.contains("pgex_oe_sub_stage1_right");
            has_mutant |= seq_id.contains("pgex_oe_sub_mutant");
        }
        assert!(has_outer, "expected outer primer artifact");
        assert!(has_stage1_left, "expected stage-1 left artifact");
        assert!(has_stage1_right, "expected stage-1 right artifact");
        assert!(has_mutant, "expected mutant product artifact");
    }
}
