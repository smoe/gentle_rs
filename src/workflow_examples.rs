//! Curated workflow example payloads and templates.

use crate::engine::{Engine, GentleEngine, Operation, ProjectState, Workflow};
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
pub const TUTORIAL_MANIFEST_SCHEMA: &str = "gentle.tutorial_manifest.v1";
pub const TUTORIAL_GENERATION_REPORT_SCHEMA: &str = "gentle.tutorial_generation_report.v1";
pub const DEFAULT_TUTORIAL_MANIFEST_PATH: &str = "docs/tutorial/manifest.json";
pub const DEFAULT_TUTORIAL_OUTPUT_DIR: &str = "docs/tutorial/generated";

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
    pub summary: String,
    #[serde(default)]
    pub retain_outputs: Vec<String>,
    #[serde(default)]
    pub checkpoints: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TutorialManifest {
    #[serde(default = "default_tutorial_manifest_schema")]
    pub schema: String,
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

pub fn online_example_tests_enabled() -> bool {
    let raw = env::var(ONLINE_EXAMPLE_TEST_ENV).unwrap_or_default();
    matches!(
        raw.trim().to_ascii_lowercase().as_str(),
        "1" | "true" | "yes" | "on"
    )
}

pub fn run_example_workflow(example: &WorkflowExample) -> Result<(), String> {
    let mut engine = GentleEngine::from_state(ProjectState::default());
    engine
        .apply_workflow(example.workflow.clone())
        .map(|_| ())
        .map_err(|e| format!("Workflow example '{}' failed: {e}", example.id))
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
    if manifest.chapters.is_empty() {
        return Err(format!(
            "Tutorial manifest '{}' has no chapters",
            display_path(manifest_path)
        ));
    }
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

pub fn validate_tutorial_manifest_against_examples(
    manifest: &TutorialManifest,
    examples: &[LoadedWorkflowExample],
) -> Result<(), String> {
    let by_id = example_lookup(examples);
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
    }
    Ok(rewritten)
}

fn run_example_workflow_in_dir(
    example: &WorkflowExample,
    repo_root: &Path,
    run_dir: &Path,
) -> Result<(), String> {
    let rewritten = rewrite_example_paths_for_execution(example, repo_root, run_dir)?;
    run_example_workflow(&rewritten)
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
) -> Result<String, String> {
    let workflow_json = serde_json::to_string_pretty(&loaded.example.workflow).map_err(|e| {
        format!(
            "Could not serialize workflow example '{}': {e}",
            loaded.example.id
        )
    })?;
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
    out.push_str("\n## Run Commands\n\n");
    out.push_str("```bash\n");
    out.push_str("cargo run --bin gentle_cli -- workflow @");
    out.push_str(&source_path);
    out.push_str("\n");
    out.push_str("cargo run --bin gentle_cli -- shell 'workflow @");
    out.push_str(&source_path);
    out.push_str("'\n");
    out.push_str("```\n");
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
    out.push_str("\n## Canonical Workflow JSON\n\n");
    out.push_str("```json\n");
    out.push_str(&workflow_json);
    out.push_str("\n```\n");
    Ok(out)
}

fn render_tutorial_index(
    manifest: &TutorialManifest,
    chapters: &[TutorialGenerationChapter],
    online_enabled: bool,
) -> String {
    let mut out = String::new();
    out.push_str("# GENtle Tutorial (Generated)\n\n");
    out.push_str("This folder is generated from:\n");
    out.push_str("- `docs/tutorial/manifest.json`\n");
    out.push_str("- `docs/examples/workflows/*.json`\n\n");
    out.push_str("Regenerate with:\n\n");
    out.push_str("```bash\n");
    out.push_str("cargo run --bin gentle_examples_docs -- tutorial-generate\n");
    out.push_str("```\n\n");
    out.push_str("Validate committed generated output:\n\n");
    out.push_str("```bash\n");
    out.push_str("cargo run --bin gentle_examples_docs -- tutorial-check\n");
    out.push_str("```\n\n");
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

    for chapter in sorted_manifest_chapters(&manifest) {
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
                    copy_retained_output(&chapter, run_dir.path(), output_dir, retain_output)?;
                retained_artifacts.push(retained);
            }
        }
        retained_artifacts.sort();
        let chapter_markdown = render_tutorial_chapter_markdown(
            &chapter,
            &loaded,
            executed,
            &retained_artifacts,
            online_enabled,
        )?;
        let chapter_file = markdown_path_for_chapter(&chapter);
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
            executed,
            retained_artifacts,
        });
    }

    let readme = render_tutorial_index(&manifest, &chapter_reports, online_enabled);
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

    fn example_dir() -> PathBuf {
        PathBuf::from(DEFAULT_WORKFLOW_EXAMPLE_DIR)
    }

    fn tutorial_manifest_path() -> PathBuf {
        PathBuf::from(DEFAULT_TUTORIAL_MANIFEST_PATH)
    }

    fn tutorial_output_dir() -> PathBuf {
        PathBuf::from(DEFAULT_TUTORIAL_OUTPUT_DIR)
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
}
