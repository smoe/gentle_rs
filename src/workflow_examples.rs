//! Curated workflow example payloads and templates.

use crate::engine::{Engine, GentleEngine, Operation, OperationProgress, ProjectState, Workflow};
use serde::{Deserialize, Serialize};
use sha1::{Digest, Sha1};
use std::{
    collections::{BTreeMap, HashMap, HashSet},
    env, fs,
    io::Read,
    path::{Path, PathBuf},
    time::{SystemTime, UNIX_EPOCH},
};
use tempfile::TempDir;

pub const WORKFLOW_EXAMPLE_SCHEMA: &str = "gentle.workflow_example.v1";
pub const DEFAULT_WORKFLOW_EXAMPLE_DIR: &str = "docs/examples/workflows";
pub const DEFAULT_WORKFLOW_SNIPPET_DIR: &str = "docs/examples/generated";
pub const ONLINE_EXAMPLE_TEST_ENV: &str = "GENTLE_TEST_ONLINE";
pub const SKIP_REMOTE_TESTS_ENV: &str = "GENTLE_SKIP_REMOTE_TESTS";
pub const TUTORIAL_CATALOG_SCHEMA: &str = "gentle.tutorial_catalog.v2";
pub const LEGACY_TUTORIAL_CATALOG_SCHEMA_V1: &str = "gentle.tutorial_catalog.v1";
pub const TUTORIAL_CATALOG_META_SCHEMA: &str = "gentle.tutorial_catalog_meta.v2";
pub const LEGACY_TUTORIAL_CATALOG_META_SCHEMA_V1: &str = "gentle.tutorial_catalog_meta.v1";
pub const TUTORIAL_SOURCE_SCHEMA: &str = "gentle.tutorial_source.v4";
pub const LEGACY_TUTORIAL_SOURCE_SCHEMA_V3: &str = "gentle.tutorial_source.v3";
pub const LEGACY_TUTORIAL_SOURCE_SCHEMA_V2: &str = "gentle.tutorial_source.v2";
pub const TUTORIAL_MANIFEST_SCHEMA: &str = "gentle.tutorial_manifest.v2";
pub const LEGACY_TUTORIAL_MANIFEST_SCHEMA_V1: &str = "gentle.tutorial_manifest.v1";
pub const TUTORIAL_GENERATION_REPORT_SCHEMA: &str = "gentle.tutorial_generation_report.v1";
pub const TUTORIAL_REVIEW_MANIFEST_SCHEMA: &str = "gentle.tutorial_review_manifest.v1";
pub const DEFAULT_TUTORIAL_CATALOG_PATH: &str = "docs/tutorial/catalog.json";
pub const DEFAULT_TUTORIAL_CATALOG_META_PATH: &str = "docs/tutorial/sources/catalog_meta.json";
pub const DEFAULT_TUTORIAL_SOURCE_DIR: &str = "docs/tutorial/sources";
pub const DEFAULT_TUTORIAL_MANIFEST_PATH: &str = "docs/tutorial/manifest.json";
pub const DEFAULT_TUTORIAL_OUTPUT_DIR: &str = "docs/tutorial/generated";
pub const DEFAULT_TUTORIAL_REVIEW_MANIFEST_PATH: &str = "docs/tutorial/review_manifest.json";
pub const TUTORIAL_CONFUSION_ISSUE_TEMPLATE_PATH: &str =
    ".github/ISSUE_TEMPLATE/tutorial-confusion.md";
pub const TUTORIAL_ARTIFACT_ISSUE_TEMPLATE_PATH: &str =
    ".github/ISSUE_TEMPLATE/tutorial-artifact-figure.md";
pub const TUTORIAL_EXECUTION_ISSUE_TEMPLATE_PATH: &str =
    ".github/ISSUE_TEMPLATE/tutorial-execution-failure.md";

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
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub group: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub group_label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub group_order: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub group_position: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub decimal_id: Option<String>,
    #[serde(rename = "type")]
    pub entry_type: String,
    pub status: String,
    pub source: String,
    #[serde(default)]
    pub audiences: Vec<String>,
    #[serde(default)]
    pub notes: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub review_status: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub codex_reviewed_at: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub human_reviewed_at: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub human_reviewer: Option<String>,
    #[serde(default, skip_serializing_if = "is_false")]
    pub review_stale: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub review_stale_reason: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub review_issue_template: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub review_issue_template_path: Option<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub graphics: Vec<TutorialGraphic>,
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
    pub groups: Vec<TutorialGroupDefinition>,
    #[serde(default)]
    pub concepts: Vec<TutorialConcept>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TutorialGroupDefinition {
    pub code: String,
    pub label: String,
    pub order: usize,
}

#[derive(Debug, Clone, Default)]
struct TutorialPlacement {
    group: Option<String>,
    group_label: Option<String>,
    group_order: Option<usize>,
    group_position: Option<usize>,
    decimal_id: Option<String>,
}

#[derive(Debug, Clone, Default)]
struct TutorialReviewProjection {
    status: String,
    codex_reviewed_at: Option<String>,
    human_reviewed_at: Option<String>,
    human_reviewer: Option<String>,
    stale: bool,
    stale_reason: Option<String>,
    issue_template: Option<String>,
    issue_template_path: Option<String>,
}

#[derive(Debug, Clone)]
struct TutorialReviewDependency {
    label: String,
    path: PathBuf,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TutorialSourceCatalogSection {
    pub order: usize,
    pub path: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub group: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub group_position: Option<usize>,
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
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub cli_steps: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub step_expectations: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub prerequisites: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub local_execution_note: Option<String>,
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
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub graphics: Vec<TutorialGraphic>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TutorialSourceUnit {
    #[serde(default = "default_tutorial_source_schema")]
    pub schema: String,
    pub id: String,
    pub title: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub group: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub group_position: Option<usize>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub graphics: Vec<TutorialGraphic>,
    #[serde(default)]
    pub catalog: Option<TutorialSourceCatalogSection>,
    #[serde(default)]
    pub generated_chapter: Option<TutorialSourceGeneratedChapterSection>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TutorialGraphic {
    pub kind: String,
    pub path: String,
    pub caption: String,
    pub illustrates_step: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub capture_date: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub regen_command: Option<String>,
}

impl TutorialSourceCatalogSection {
    fn into_catalog_entry(
        self,
        id: String,
        title: String,
        placement: TutorialPlacement,
        review: TutorialReviewProjection,
        graphics: Vec<TutorialGraphic>,
    ) -> TutorialCatalogEntry {
        TutorialCatalogEntry {
            id,
            title,
            path: self.path,
            group: placement.group,
            group_label: placement.group_label,
            group_order: placement.group_order,
            group_position: placement.group_position,
            decimal_id: placement.decimal_id,
            entry_type: self.entry_type,
            status: self.status,
            source: self.source,
            audiences: self.audiences,
            notes: self.notes,
            review_status: Some(review.status),
            codex_reviewed_at: review.codex_reviewed_at,
            human_reviewed_at: review.human_reviewed_at,
            human_reviewer: review.human_reviewer,
            review_stale: review.stale,
            review_stale_reason: review.stale_reason,
            review_issue_template: review.issue_template,
            review_issue_template_path: review.issue_template_path,
            graphics,
        }
    }
}

impl TutorialSourceGeneratedChapterSection {
    fn into_manifest_chapter(
        self,
        id: String,
        title: String,
        placement: TutorialPlacement,
        graphics: Vec<TutorialGraphic>,
    ) -> TutorialChapter {
        TutorialChapter {
            id,
            order: self.order,
            title,
            group: placement.group,
            group_label: placement.group_label,
            group_order: placement.group_order,
            group_position: placement.group_position,
            decimal_id: placement.decimal_id,
            example_id: self.example_id,
            tier: self.tier,
            guide_path: self.guide_path,
            summary: self.summary,
            narrative: self.narrative,
            use_cases: self.use_cases,
            gui_steps: self.gui_steps,
            cli_steps: self.cli_steps,
            step_expectations: self.step_expectations,
            prerequisites: self.prerequisites,
            local_execution_note: self.local_execution_note,
            learning_objectives: self.learning_objectives,
            concepts: self.concepts,
            parameter_notes: self.parameter_notes,
            follow_up_commands: self.follow_up_commands,
            retain_outputs: self.retain_outputs,
            checkpoints: self.checkpoints,
            graphics,
        }
    }
}

impl TutorialSourceUnit {
    fn into_catalog_entry(
        self,
        group_lookup: &HashMap<String, TutorialGroupDefinition>,
        review_by_id: &HashMap<String, TutorialReviewEntry>,
        review_today: Option<(i32, u32, u32)>,
        review_warn_after_months: u32,
        repo_root: &Path,
        source_path_by_id: &HashMap<String, PathBuf>,
    ) -> Result<Option<(usize, TutorialCatalogEntry)>, String> {
        let TutorialSourceUnit {
            id,
            title,
            group,
            group_position,
            graphics,
            catalog,
            generated_chapter,
            ..
        } = self;
        let Some(catalog) = catalog else {
            return Ok(None);
        };
        let effective_group = catalog.group.as_deref().or(group.as_deref());
        let effective_group_position = catalog.group_position.or(group_position);
        let placement = tutorial_source_placement(
            &id,
            effective_group,
            effective_group_position,
            group_lookup,
        )?;
        let review = tutorial_review_projection(
            review_by_id.get(&id),
            review_today,
            review_warn_after_months,
            tutorial_review_dependency_stale_reason(
                review_by_id.get(&id),
                &tutorial_source_dependencies(
                    &id,
                    source_path_by_id.get(&id),
                    catalog.path.as_str(),
                    generated_chapter
                        .as_ref()
                        .map(|chapter| chapter.example_id.as_str()),
                    &graphics,
                    repo_root,
                ),
            ),
        );
        Ok(Some({
            let order = catalog.order;
            let entry = catalog.into_catalog_entry(id, title, placement, review, graphics);
            (order, entry)
        }))
    }

    fn into_manifest_chapter(
        self,
        group_lookup: &HashMap<String, TutorialGroupDefinition>,
    ) -> Result<Option<TutorialChapter>, String> {
        let TutorialSourceUnit {
            id,
            title,
            group,
            group_position,
            graphics,
            catalog,
            generated_chapter,
            ..
        } = self;
        let effective_group = catalog
            .as_ref()
            .and_then(|catalog| catalog.group.as_deref())
            .or(group.as_deref());
        let effective_group_position = catalog
            .as_ref()
            .and_then(|catalog| catalog.group_position)
            .or(group_position);
        let placement = tutorial_source_placement(
            &id,
            effective_group,
            effective_group_position,
            group_lookup,
        )?;
        Ok(generated_chapter
            .map(|generated| generated.into_manifest_chapter(id, title, placement, graphics)))
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
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub group: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub group_label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub group_order: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub group_position: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub decimal_id: Option<String>,
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
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub cli_steps: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub step_expectations: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub prerequisites: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub local_execution_note: Option<String>,
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
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub graphics: Vec<TutorialGraphic>,
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
pub struct TutorialReviewManifest {
    #[serde(default = "default_tutorial_review_manifest_schema")]
    pub schema: String,
    #[serde(default = "default_tutorial_review_warn_after_months")]
    pub warn_after_months: u32,
    #[serde(default)]
    pub entries: Vec<TutorialReviewEntry>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TutorialReviewEntry {
    pub tutorial_id: String,
    pub tutorial_kind: String,
    pub tutorial_status: String,
    #[serde(default)]
    pub replaced_by: Option<String>,
    #[serde(default)]
    pub codex_reviewed_at: Option<String>,
    #[serde(default)]
    pub human_reviewed_at: Option<String>,
    #[serde(default)]
    pub human_reviewer: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TutorialGenerationChapter {
    pub id: String,
    pub order: usize,
    pub title: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub group: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub group_label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub group_order: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub group_position: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub decimal_id: Option<String>,
    pub tier: String,
    pub example_id: String,
    pub example_source: String,
    pub concepts: Vec<String>,
    pub executed: bool,
    #[serde(default)]
    pub review_status: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub codex_reviewed_at: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub human_reviewed_at: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub human_reviewer: Option<String>,
    #[serde(default, skip_serializing_if = "is_false")]
    pub review_stale: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub review_stale_reason: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub review_issue_template: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub review_issue_template_path: Option<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub graphics: Vec<TutorialGraphic>,
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
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub warnings: Vec<String>,
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

fn is_false(value: &bool) -> bool {
    !*value
}

fn default_tutorial_source_schema() -> String {
    TUTORIAL_SOURCE_SCHEMA.to_string()
}

fn is_supported_tutorial_source_schema(schema: &str) -> bool {
    matches!(
        schema,
        TUTORIAL_SOURCE_SCHEMA
            | LEGACY_TUTORIAL_SOURCE_SCHEMA_V3
            | LEGACY_TUTORIAL_SOURCE_SCHEMA_V2
    )
}

fn is_supported_tutorial_catalog_schema(schema: &str) -> bool {
    matches!(
        schema,
        TUTORIAL_CATALOG_SCHEMA | LEGACY_TUTORIAL_CATALOG_SCHEMA_V1
    )
}

fn is_supported_tutorial_catalog_meta_schema(schema: &str) -> bool {
    matches!(
        schema,
        TUTORIAL_CATALOG_META_SCHEMA | LEGACY_TUTORIAL_CATALOG_META_SCHEMA_V1
    )
}

fn is_supported_tutorial_manifest_schema(schema: &str) -> bool {
    matches!(
        schema,
        TUTORIAL_MANIFEST_SCHEMA | LEGACY_TUTORIAL_MANIFEST_SCHEMA_V1
    )
}

fn default_tutorial_manifest_schema() -> String {
    TUTORIAL_MANIFEST_SCHEMA.to_string()
}

fn default_tutorial_generation_report_schema() -> String {
    TUTORIAL_GENERATION_REPORT_SCHEMA.to_string()
}

fn default_tutorial_review_manifest_schema() -> String {
    TUTORIAL_REVIEW_MANIFEST_SCHEMA.to_string()
}

fn default_tutorial_review_warn_after_months() -> u32 {
    12
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

fn parse_tutorial_review_manifest(manifest_path: &Path) -> Result<TutorialReviewManifest, String> {
    let raw = fs::read_to_string(manifest_path).map_err(|e| {
        format!(
            "Could not read tutorial review manifest '{}': {e}",
            display_path(manifest_path)
        )
    })?;
    serde_json::from_str::<TutorialReviewManifest>(&raw).map_err(|e| {
        format!(
            "Could not parse tutorial review manifest '{}': {e}",
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

fn validate_tutorial_groups(groups: &[TutorialGroupDefinition]) -> Result<(), String> {
    if groups.is_empty() {
        return Err("v2 catalog metadata must define at least one tutorial group".to_string());
    }
    let mut seen_codes = HashSet::new();
    let mut seen_orders = HashSet::new();
    for group in groups {
        if group.code.trim().is_empty() || group.label.trim().is_empty() || group.order == 0 {
            return Err("tutorial groups require non-empty code/label and order > 0".to_string());
        }
        if !seen_codes.insert(group.code.clone()) {
            return Err(format!("duplicate tutorial group code '{}'", group.code));
        }
        if !seen_orders.insert(group.order) {
            return Err(format!("duplicate tutorial group order '{}'", group.order));
        }
    }
    Ok(())
}

fn tutorial_group_lookup(
    groups: &[TutorialGroupDefinition],
) -> Result<HashMap<String, TutorialGroupDefinition>, String> {
    validate_tutorial_groups(groups)?;
    Ok(groups
        .iter()
        .cloned()
        .map(|group| (group.code.clone(), group))
        .collect())
}

fn tutorial_decimal_id(group_order: usize, group_position: usize) -> String {
    format!("{group_order:02}.{group_position:02}")
}

fn tutorial_source_placement(
    tutorial_id: &str,
    group_code: Option<&str>,
    group_position: Option<usize>,
    group_lookup: &HashMap<String, TutorialGroupDefinition>,
) -> Result<TutorialPlacement, String> {
    let Some(raw_group_code) = group_code.map(str::trim).filter(|value| !value.is_empty()) else {
        if group_position.is_some() {
            return Err(format!(
                "Tutorial source '{}' defines group_position without group",
                tutorial_id
            ));
        }
        return Ok(TutorialPlacement::default());
    };
    let group = group_lookup.get(raw_group_code).ok_or_else(|| {
        format!(
            "Tutorial source '{}' references unknown tutorial group '{}'",
            tutorial_id, raw_group_code
        )
    })?;
    let decimal_id = group_position.map(|position| tutorial_decimal_id(group.order, position));
    Ok(TutorialPlacement {
        group: Some(group.code.clone()),
        group_label: Some(group.label.clone()),
        group_order: Some(group.order),
        group_position,
        decimal_id,
    })
}

pub fn load_tutorial_catalog(catalog_path: &Path) -> Result<TutorialCatalog, String> {
    let catalog = parse_tutorial_catalog(catalog_path)?;
    if !is_supported_tutorial_catalog_schema(&catalog.schema) {
        return Err(format!(
            "Tutorial catalog '{}' uses unsupported schema '{}'; expected '{}' or '{}'",
            display_path(catalog_path),
            catalog.schema,
            TUTORIAL_CATALOG_SCHEMA,
            LEGACY_TUTORIAL_CATALOG_SCHEMA_V1
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
    if !is_supported_tutorial_catalog_meta_schema(&meta.schema) {
        return Err(format!(
            "Tutorial catalog meta '{}' uses unsupported schema '{}'; expected '{}' or '{}'",
            display_path(meta_path),
            meta.schema,
            TUTORIAL_CATALOG_META_SCHEMA,
            LEGACY_TUTORIAL_CATALOG_META_SCHEMA_V1
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
    if meta.schema == TUTORIAL_CATALOG_META_SCHEMA {
        validate_tutorial_groups(&meta.groups).map_err(|e| {
            format!(
                "Tutorial catalog meta '{}' group validation error: {e}",
                display_path(meta_path)
            )
        })?;
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
        if !is_supported_tutorial_source_schema(&unit.schema) {
            return Err(format!(
                "Tutorial source '{}' uses unsupported schema '{}'; expected '{}', '{}', or '{}'",
                display_path(&path),
                unit.schema,
                TUTORIAL_SOURCE_SCHEMA,
                LEGACY_TUTORIAL_SOURCE_SCHEMA_V3,
                LEGACY_TUTORIAL_SOURCE_SCHEMA_V2
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
        for graphic in &unit.graphics {
            let kind = graphic.kind.trim();
            if !matches!(kind, "generated" | "screenshot") {
                return Err(format!(
                    "Tutorial source '{}' graphic '{}' has unsupported kind '{}'",
                    display_path(&path),
                    graphic.path,
                    graphic.kind
                ));
            }
            if graphic.path.trim().is_empty() || graphic.caption.trim().is_empty() {
                return Err(format!(
                    "Tutorial source '{}' contains a graphic with blank path or caption",
                    display_path(&path)
                ));
            }
            if graphic.illustrates_step == 0 {
                return Err(format!(
                    "Tutorial source '{}' graphic '{}' must reference a 1-based step",
                    display_path(&path),
                    graphic.path
                ));
            }
            if kind == "screenshot"
                && graphic
                    .capture_date
                    .as_deref()
                    .map(str::trim)
                    .unwrap_or_default()
                    .is_empty()
            {
                return Err(format!(
                    "Tutorial source '{}' screenshot graphic '{}' needs capture_date",
                    display_path(&path),
                    graphic.path
                ));
            }
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
            if (unit.schema == TUTORIAL_SOURCE_SCHEMA
                || unit.schema == LEGACY_TUTORIAL_SOURCE_SCHEMA_V3)
                && generated.tier == TutorialTier::Online
                && generated
                    .local_execution_note
                    .as_deref()
                    .map(str::trim)
                    .unwrap_or_default()
                    .is_empty()
            {
                return Err(format!(
                    "Tutorial source '{}' uses v3 online chapter fields and must define `local_execution_note`",
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
    let group_lookup = if meta.groups.is_empty() {
        HashMap::new()
    } else {
        tutorial_group_lookup(&meta.groups)?
    };
    let (review_by_id, review_warn_after_months) =
        tutorial_review_entries_for_source_dir(source_dir)?;
    let review_today = current_utc_review_date();
    let repo_root = tutorial_repo_root_for_source_dir(source_dir);
    let source_path_by_id = tutorial_source_path_lookup(source_dir);
    let mut entries = Vec::new();
    for unit in units {
        if let Some(entry) = unit.into_catalog_entry(
            &group_lookup,
            &review_by_id,
            review_today,
            review_warn_after_months,
            &repo_root,
            &source_path_by_id,
        )? {
            entries.push(entry);
        }
    }
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
    let group_lookup = if meta.groups.is_empty() {
        HashMap::new()
    } else {
        tutorial_group_lookup(&meta.groups)?
    };
    let mut chapters = Vec::new();
    for unit in units {
        if let Some(chapter) = unit.into_manifest_chapter(&group_lookup)? {
            chapters.push(chapter);
        }
    }
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
    if !is_supported_tutorial_manifest_schema(&manifest.schema) {
        return Err(format!(
            "Tutorial manifest '{}' uses unsupported schema '{}'; expected '{}' or '{}'",
            display_path(manifest_path),
            manifest.schema,
            TUTORIAL_MANIFEST_SCHEMA,
            LEGACY_TUTORIAL_MANIFEST_SCHEMA_V1
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
    let chapter_ids = manifest
        .chapters
        .iter()
        .map(|chapter| chapter.id.clone())
        .collect::<HashSet<_>>();
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
        for graphic in &chapter.graphics {
            let kind = graphic.kind.trim();
            if !matches!(kind, "generated" | "screenshot") {
                return Err(format!(
                    "Tutorial chapter '{}' graphic '{}' has unsupported kind '{}'",
                    chapter.id, graphic.path, graphic.kind
                ));
            }
            if graphic.path.trim().is_empty() {
                return Err(format!(
                    "Tutorial chapter '{}' contains graphic with blank path",
                    chapter.id
                ));
            }
            if graphic.caption.trim().is_empty() {
                return Err(format!(
                    "Tutorial chapter '{}' graphic '{}' needs a caption",
                    chapter.id, graphic.path
                ));
            }
            if graphic.illustrates_step == 0 || graphic.illustrates_step > chapter.gui_steps.len() {
                return Err(format!(
                    "Tutorial chapter '{}' graphic '{}' references step {}, but the chapter has {} GUI steps",
                    chapter.id,
                    graphic.path,
                    graphic.illustrates_step,
                    chapter.gui_steps.len()
                ));
            }
            if kind == "screenshot"
                && graphic
                    .capture_date
                    .as_deref()
                    .map(str::trim)
                    .unwrap_or_default()
                    .is_empty()
            {
                return Err(format!(
                    "Tutorial chapter '{}' screenshot graphic '{}' needs capture_date",
                    chapter.id, graphic.path
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
        if chapter.cli_steps.len() > chapter.gui_steps.len() {
            return Err(format!(
                "Tutorial chapter '{}' has more cli_steps than gui_steps",
                chapter.id
            ));
        }
        for cli_step in &chapter.cli_steps {
            if cli_step.trim().is_empty() {
                return Err(format!(
                    "Tutorial chapter '{}' contains blank cli_steps entry",
                    chapter.id
                ));
            }
        }
        if chapter.step_expectations.len() > chapter.gui_steps.len() {
            return Err(format!(
                "Tutorial chapter '{}' has more step_expectations than gui_steps",
                chapter.id
            ));
        }
        for expectation in &chapter.step_expectations {
            if expectation.trim().is_empty() {
                return Err(format!(
                    "Tutorial chapter '{}' contains blank step_expectations entry",
                    chapter.id
                ));
            }
        }
        if let Some(local_execution_note) = &chapter.local_execution_note {
            if local_execution_note.trim().is_empty() {
                return Err(format!(
                    "Tutorial chapter '{}' contains blank local_execution_note",
                    chapter.id
                ));
            }
        }
        for prerequisite in &chapter.prerequisites {
            if prerequisite.trim().is_empty() {
                return Err(format!(
                    "Tutorial chapter '{}' contains blank prerequisites entry",
                    chapter.id
                ));
            }
            if prerequisite == &chapter.id {
                return Err(format!(
                    "Tutorial chapter '{}' cannot list itself as a prerequisite",
                    chapter.id
                ));
            }
            if !chapter_ids.contains(prerequisite) {
                return Err(format!(
                    "Tutorial chapter '{}' references unknown prerequisite '{}'",
                    chapter.id, prerequisite
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

pub fn load_tutorial_review_manifest(
    manifest_path: &Path,
) -> Result<TutorialReviewManifest, String> {
    let manifest = parse_tutorial_review_manifest(manifest_path)?;
    if manifest.schema != TUTORIAL_REVIEW_MANIFEST_SCHEMA {
        return Err(format!(
            "Tutorial review manifest '{}' uses unsupported schema '{}'; expected '{}'",
            display_path(manifest_path),
            manifest.schema,
            TUTORIAL_REVIEW_MANIFEST_SCHEMA
        ));
    }
    if manifest.warn_after_months == 0 {
        return Err(format!(
            "Tutorial review manifest '{}' must set warn_after_months greater than 0",
            display_path(manifest_path)
        ));
    }
    let mut seen_ids: HashSet<String> = HashSet::new();
    for entry in &manifest.entries {
        if entry.tutorial_id.trim().is_empty()
            || entry.tutorial_kind.trim().is_empty()
            || entry.tutorial_status.trim().is_empty()
        {
            return Err(format!(
                "Tutorial review manifest '{}' contains an entry with empty required fields",
                display_path(manifest_path)
            ));
        }
        if !matches!(
            entry.tutorial_kind.as_str(),
            "guided_walkthrough" | "generated_chapter"
        ) {
            return Err(format!(
                "Tutorial review manifest '{}' entry '{}' has unsupported tutorial_kind '{}'",
                display_path(manifest_path),
                entry.tutorial_id,
                entry.tutorial_kind
            ));
        }
        if !matches!(entry.tutorial_status.as_str(), "active" | "deprecated") {
            return Err(format!(
                "Tutorial review manifest '{}' entry '{}' has unsupported tutorial_status '{}'",
                display_path(manifest_path),
                entry.tutorial_id,
                entry.tutorial_status
            ));
        }
        if !seen_ids.insert(entry.tutorial_id.clone()) {
            return Err(format!(
                "Tutorial review manifest '{}' has duplicate tutorial id '{}'",
                display_path(manifest_path),
                entry.tutorial_id
            ));
        }
        validate_optional_review_date(
            manifest_path,
            &entry.tutorial_id,
            "codex_reviewed_at",
            entry.codex_reviewed_at.as_deref(),
        )?;
        validate_optional_review_date(
            manifest_path,
            &entry.tutorial_id,
            "human_reviewed_at",
            entry.human_reviewed_at.as_deref(),
        )?;
    }
    Ok(manifest)
}

fn validate_optional_review_date(
    manifest_path: &Path,
    tutorial_id: &str,
    field: &str,
    value: Option<&str>,
) -> Result<(), String> {
    if let Some(value) = value {
        parse_review_date(value).map_err(|e| {
            format!(
                "Tutorial review manifest '{}' entry '{}' has invalid {} '{}': {e}",
                display_path(manifest_path),
                tutorial_id,
                field,
                value
            )
        })?;
    }
    Ok(())
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

fn parse_review_date(raw: &str) -> Result<(i32, u32, u32), String> {
    let mut parts = raw.split('-');
    let year = parts
        .next()
        .ok_or_else(|| "missing year".to_string())?
        .parse::<i32>()
        .map_err(|e| format!("invalid year: {e}"))?;
    let month = parts
        .next()
        .ok_or_else(|| "missing month".to_string())?
        .parse::<u32>()
        .map_err(|e| format!("invalid month: {e}"))?;
    let day = parts
        .next()
        .ok_or_else(|| "missing day".to_string())?
        .parse::<u32>()
        .map_err(|e| format!("invalid day: {e}"))?;
    if parts.next().is_some() {
        return Err("expected YYYY-MM-DD".to_string());
    }
    if !(1..=12).contains(&month) {
        return Err("month must be in 1..=12".to_string());
    }
    let max_day = match month {
        1 | 3 | 5 | 7 | 8 | 10 | 12 => 31,
        4 | 6 | 9 | 11 => 30,
        2 if is_leap_year(year) => 29,
        2 => 28,
        _ => unreachable!(),
    };
    if day == 0 || day > max_day {
        return Err(format!("day must be in 1..={max_day}"));
    }
    Ok((year, month, day))
}

fn is_leap_year(year: i32) -> bool {
    (year % 4 == 0 && year % 100 != 0) || year % 400 == 0
}

fn review_date_is_stale(reviewed_at: &str, today: (i32, u32, u32), warn_after_months: u32) -> bool {
    let Ok((review_year, review_month, review_day)) = parse_review_date(reviewed_at) else {
        return false;
    };
    let review_total_months = review_year as i64 * 12 + review_month as i64;
    let today_total_months = today.0 as i64 * 12 + today.1 as i64;
    let month_delta = today_total_months - review_total_months;
    month_delta > warn_after_months as i64
        || (month_delta == warn_after_months as i64 && today.2 >= review_day)
}

fn review_date_after(left: (i32, u32, u32), right: (i32, u32, u32)) -> bool {
    left > right
}

fn system_time_review_date(time: SystemTime) -> Option<(i32, u32, u32)> {
    let seconds = time.duration_since(UNIX_EPOCH).ok()?.as_secs();
    Some(civil_from_unix_days((seconds / 86_400) as i64))
}

fn path_modified_review_date(path: &Path) -> Option<(i32, u32, u32)> {
    fs::metadata(path)
        .ok()
        .and_then(|metadata| metadata.modified().ok())
        .and_then(system_time_review_date)
}

fn tutorial_repo_root_for_source_dir(source_dir: &Path) -> PathBuf {
    source_dir
        .parent()
        .and_then(Path::parent)
        .and_then(Path::parent)
        .map(Path::to_path_buf)
        .unwrap_or_else(|| PathBuf::from("."))
}

fn tutorial_issue_template_name(path: &str) -> &'static str {
    match path {
        TUTORIAL_ARTIFACT_ISSUE_TEMPLATE_PATH => "Tutorial artifact/figure problem",
        TUTORIAL_EXECUTION_ISSUE_TEMPLATE_PATH => "Tutorial execution failure",
        _ => "Tutorial confusion",
    }
}

fn current_utc_review_date() -> Option<(i32, u32, u32)> {
    let seconds = SystemTime::now().duration_since(UNIX_EPOCH).ok()?.as_secs();
    Some(civil_from_unix_days((seconds / 86_400) as i64))
}

fn civil_from_unix_days(days_since_epoch: i64) -> (i32, u32, u32) {
    let z = days_since_epoch + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1_460 + doe / 36_524 - doe / 146_096) / 365;
    let year = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = mp + if mp < 10 { 3 } else { -9 };
    let year = year + if month <= 2 { 1 } else { 0 };
    (year as i32, month as u32, day as u32)
}

fn chapter_link_label(chapter: &TutorialChapter) -> String {
    format!("Chapter {}: {}", chapter.order, chapter.title)
}

fn chapter_link_target(chapter: &TutorialChapter, from_readme: bool) -> String {
    let filename = markdown_path_for_chapter(chapter);
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

#[derive(Debug, Clone, Default)]
struct TutorialReviewContext {
    entries_by_id: HashMap<String, TutorialReviewEntry>,
    stale_reasons_by_id: HashMap<String, String>,
    warn_after_months: u32,
    today: Option<(i32, u32, u32)>,
    warnings: Vec<String>,
}

fn tutorial_review_manifest_path_for_manifest(manifest_path: &Path) -> PathBuf {
    manifest_path
        .parent()
        .map(|parent| parent.join("review_manifest.json"))
        .unwrap_or_else(|| PathBuf::from(DEFAULT_TUTORIAL_REVIEW_MANIFEST_PATH))
}

fn tutorial_review_manifest_path_for_source_dir(source_dir: &Path) -> PathBuf {
    source_dir
        .parent()
        .map(|parent| parent.join("review_manifest.json"))
        .unwrap_or_else(|| PathBuf::from(DEFAULT_TUTORIAL_REVIEW_MANIFEST_PATH))
}

fn tutorial_review_entries_for_source_dir(
    source_dir: &Path,
) -> Result<(HashMap<String, TutorialReviewEntry>, u32), String> {
    let review_path = tutorial_review_manifest_path_for_source_dir(source_dir);
    if !review_path.exists() {
        return Ok((HashMap::new(), default_tutorial_review_warn_after_months()));
    }
    let review_manifest = load_tutorial_review_manifest(&review_path)?;
    Ok((
        review_manifest
            .entries
            .into_iter()
            .map(|entry| (entry.tutorial_id.clone(), entry))
            .collect(),
        review_manifest.warn_after_months,
    ))
}

fn tutorial_known_review_ids(
    manifest_path: &Path,
    manifest: &TutorialManifest,
) -> Result<BTreeMap<String, String>, String> {
    let mut known = BTreeMap::new();
    for chapter in &manifest.chapters {
        known.insert(chapter.id.clone(), "generated_chapter".to_string());
    }
    let tutorial_source_dir = tutorial_source_dir_for_manifest(manifest_path);
    if tutorial_source_dir.exists() {
        for unit in load_tutorial_source_units(&tutorial_source_dir)? {
            let kind = if unit.generated_chapter.is_some() {
                "generated_chapter"
            } else {
                "guided_walkthrough"
            };
            known.insert(unit.id, kind.to_string());
        }
    }
    Ok(known)
}

fn tutorial_review_context(
    manifest_path: &Path,
    manifest: &TutorialManifest,
) -> Result<TutorialReviewContext, String> {
    tutorial_review_context_for_date(manifest_path, manifest, current_utc_review_date())
}

fn tutorial_review_context_for_date(
    manifest_path: &Path,
    manifest: &TutorialManifest,
    today: Option<(i32, u32, u32)>,
) -> Result<TutorialReviewContext, String> {
    let review_path = tutorial_review_manifest_path_for_manifest(manifest_path);
    let known = tutorial_known_review_ids(manifest_path, manifest)?;
    if !review_path.exists() {
        return Ok(TutorialReviewContext {
            entries_by_id: HashMap::new(),
            stale_reasons_by_id: HashMap::new(),
            warn_after_months: default_tutorial_review_warn_after_months(),
            today,
            warnings: vec![format!(
                "Tutorial review manifest '{}' is missing; review freshness was not checked",
                display_path(&review_path)
            )],
        });
    }
    let review_manifest = load_tutorial_review_manifest(&review_path)?;
    let mut warnings = Vec::new();
    let mut entries_by_id = HashMap::new();
    for entry in &review_manifest.entries {
        if !known.contains_key(&entry.tutorial_id) {
            warnings.push(format!(
                "Tutorial review manifest '{}' references unknown tutorial id '{}'",
                display_path(&review_path),
                entry.tutorial_id
            ));
        }
        entries_by_id.insert(entry.tutorial_id.clone(), entry.clone());
    }
    for id in known.keys() {
        if !entries_by_id.contains_key(id) {
            warnings.push(format!(
                "Tutorial review manifest '{}' is missing entry for tutorial id '{}'",
                display_path(&review_path),
                id
            ));
        }
    }
    if let Some(today) = today {
        for entry in &review_manifest.entries {
            if let Some(reviewed_at) = entry.human_reviewed_at.as_deref() {
                if review_date_is_stale(reviewed_at, today, review_manifest.warn_after_months) {
                    warnings.push(format!(
                        "Tutorial '{}' human review date {} is older than {} months",
                        entry.tutorial_id, reviewed_at, review_manifest.warn_after_months
                    ));
                }
            }
        }
    }
    let stale_reasons_by_id =
        tutorial_review_dependency_stale_reasons(manifest_path, manifest, &entries_by_id)?;
    let mut stale_reason_rows = stale_reasons_by_id.iter().collect::<Vec<_>>();
    stale_reason_rows.sort_by(|left, right| left.0.cmp(right.0));
    for (id, reason) in stale_reason_rows {
        warnings.push(format!(
            "Tutorial '{}' review is stale because {}",
            id, reason
        ));
    }
    Ok(TutorialReviewContext {
        entries_by_id,
        stale_reasons_by_id,
        warn_after_months: review_manifest.warn_after_months,
        today,
        warnings,
    })
}

fn tutorial_review_entry_exempts_execution_failure(entry: Option<&TutorialReviewEntry>) -> bool {
    entry
        .map(|entry| {
            entry.tutorial_status == "deprecated"
                || entry
                    .replaced_by
                    .as_deref()
                    .map(str::trim)
                    .is_some_and(|value| !value.is_empty())
        })
        .unwrap_or(false)
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
        if let Operation::RenderProteinGelReportsSvg { path, .. } = op {
            *path = resolve_output_path(path, run_dir);
            ensure_parent_exists(path)?;
            continue;
        }
        if let Operation::RenderProteaseDigestGelSvg { path, .. } = op {
            *path = resolve_output_path(path, run_dir);
            ensure_parent_exists(path)?;
            continue;
        }
        if let Operation::RenderProtein2dGelSvg { path, .. } = op {
            *path = resolve_output_path(path, run_dir);
            ensure_parent_exists(path)?;
            continue;
        }
        if let Operation::RenderProtocolCartoonSvg { path, .. } = op {
            *path = resolve_output_path(path, run_dir);
            ensure_parent_exists(path)?;
            continue;
        }
        if let Operation::ExportPrimerDesignReport { path, .. } = op {
            *path = resolve_output_path(path, run_dir);
            ensure_parent_exists(path)?;
            continue;
        }
        if let Operation::TestCdnaPcr {
            path,
            svg_path,
            product_gel_svg_path,
            ..
        }
        | Operation::TestCdnaQpcr {
            path,
            svg_path,
            product_gel_svg_path,
            ..
        } = op
        {
            rewrite_optional_output_path(path, run_dir);
            if let Some(path) = path.as_deref() {
                ensure_parent_exists(path)?;
            }
            rewrite_optional_output_path(svg_path, run_dir);
            if let Some(path) = svg_path.as_deref() {
                ensure_parent_exists(path)?;
            }
            rewrite_optional_output_path(product_gel_svg_path, run_dir);
            if let Some(path) = product_gel_svg_path.as_deref() {
                ensure_parent_exists(path)?;
            }
            continue;
        }
        if let Operation::TestCdnaQpcrFasta { path, svg_path, .. } = op {
            rewrite_optional_output_path(path, run_dir);
            if let Some(path) = path.as_deref() {
                ensure_parent_exists(path)?;
            }
            rewrite_optional_output_path(svg_path, run_dir);
            if let Some(path) = svg_path.as_deref() {
                ensure_parent_exists(path)?;
            }
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
        if let Operation::ExportLabAssistantInstructions { path, .. } = op {
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
        if let Operation::SummarizeTfbsTrackSimilarity { path, .. } = op {
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
        if let Operation::RenderTfbsScoreTrackCorrelationSvg { path, .. } = op {
            *path = resolve_output_path(path, run_dir);
            ensure_parent_exists(path)?;
            continue;
        }
        if let Operation::SummarizeAlternativePromoterComparison { path, .. } = op {
            rewrite_optional_output_path(path, run_dir);
            if let Some(path) = path.as_deref() {
                ensure_parent_exists(path)?;
            }
            continue;
        }
        if let Operation::SummarizePromoterEvidenceMatrix { path, .. } = op {
            rewrite_optional_output_path(path, run_dir);
            if let Some(path) = path.as_deref() {
                ensure_parent_exists(path)?;
            }
            continue;
        }
        if let Operation::SummarizeIsoformPromoterComparison { path, .. }
        | Operation::SummarizePromoterExpressionEvidence { path, .. } = op
        {
            rewrite_optional_output_path(path, run_dir);
            if let Some(path) = path.as_deref() {
                ensure_parent_exists(path)?;
            }
            continue;
        }
        if let Operation::ExportPromoterArtifactManifest {
            artifacts, path, ..
        } = op
        {
            let resolved_manifest_path = resolve_output_path(path, run_dir);
            let manifest_parent = Path::new(&resolved_manifest_path)
                .parent()
                .map(Path::to_path_buf);
            for artifact in artifacts {
                let resolved_artifact_path =
                    PathBuf::from(resolve_output_path(&artifact.path, run_dir));
                artifact.path = manifest_parent
                    .as_ref()
                    .and_then(|parent| resolved_artifact_path.strip_prefix(parent).ok())
                    .map(display_path)
                    .unwrap_or_else(|| display_path(&resolved_artifact_path));
            }
            *path = resolved_manifest_path;
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
        if let Operation::ProjectMicroarrayTrack { manifest_path, .. } = op {
            *manifest_path = resolve_input_path(manifest_path, repo_root);
            continue;
        }
        if let Operation::QueryRepeatOverlaps {
            rmsk_index_path,
            path,
            ..
        }
        | Operation::MaterializeRepeatFeatures {
            rmsk_index_path,
            path,
            ..
        } = op
        {
            *rmsk_index_path = resolve_input_path(rmsk_index_path, repo_root);
            rewrite_optional_output_path(path, run_dir);
            if let Some(path) = path.as_deref() {
                ensure_parent_exists(path)?;
            }
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
    let bytes =
        fs::read(src).map_err(|e| format!("Could not read '{}': {e}", display_path(src)))?;
    if let Ok(text) = std::str::from_utf8(&bytes) {
        let normalized = normalize_retained_tutorial_artifact_text(text);
        fs::write(dst, normalized).map_err(|e| {
            format!(
                "Could not write normalized retained artifact '{}' to '{}': {e}",
                display_path(src),
                display_path(dst)
            )
        })?;
    } else {
        fs::write(dst, bytes).map_err(|e| {
            format!(
                "Could not copy '{}' to '{}': {e}",
                display_path(src),
                display_path(dst)
            )
        })?;
    }
    Ok(())
}

fn normalize_retained_tutorial_artifact_text(text: &str) -> String {
    let mut lines = text
        .lines()
        .map(|line| {
            if let Some(prefix) = line.strip_suffix('\r') {
                format!("{}\r", normalize_retained_tutorial_artifact_line(prefix))
            } else {
                normalize_retained_tutorial_artifact_line(line)
            }
        })
        .collect::<Vec<_>>();
    while lines.last().is_some_and(|line| line.trim().is_empty()) {
        lines.pop();
    }
    lines.join("\n") + if text.ends_with('\n') { "\n" } else { "" }
}

fn normalize_retained_tutorial_artifact_line(line: &str) -> String {
    if line.contains("\"generated_at_unix_ms\":") {
        let mut normalized = String::new();
        let mut replaced = false;
        for (idx, part) in line.split("\"generated_at_unix_ms\":").enumerate() {
            if idx == 0 {
                normalized.push_str(part);
                continue;
            }
            normalized.push_str("\"generated_at_unix_ms\":");
            let trimmed_start = part.len() - part.trim_start().len();
            normalized.push_str(&part[..trimmed_start]);
            normalized.push('0');
            let rest = &part[trimmed_start..];
            let suffix_start = rest
                .find(|ch: char| !ch.is_ascii_digit())
                .unwrap_or(rest.len());
            normalized.push_str(&rest[suffix_start..]);
            replaced = true;
        }
        if replaced {
            return normalized;
        }
    }
    if let Some(prefix) = line.strip_prefix("- Generated (Unix ms): `") {
        if let Some((_, suffix)) = prefix.split_once('`') {
            return format!("- Generated (Unix ms): `0`{suffix}");
        }
    }
    line.to_string()
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

fn tutorial_chapter_lookup(chapters: &[TutorialChapter]) -> HashMap<String, TutorialChapter> {
    chapters
        .iter()
        .map(|chapter| (chapter.id.clone(), chapter.clone()))
        .collect()
}

fn chapter_markdown_filename_from_parts(
    order: usize,
    id: &str,
    decimal_id: Option<&str>,
) -> String {
    if let Some(decimal_id) = decimal_id.map(str::trim).filter(|value| !value.is_empty()) {
        return format!(
            "{}_{}.md",
            decimal_id.replace('.', "-"),
            markdown_file_stem(id)
        );
    }
    format!("{:02}_{}.md", order, markdown_file_stem(id))
}

fn markdown_path_for_chapter(chapter: &TutorialChapter) -> String {
    chapter_markdown_filename_from_parts(chapter.order, &chapter.id, chapter.decimal_id.as_deref())
}

fn yaml_double_quote(value: &str) -> String {
    format!("\"{}\"", value.replace('\\', "\\\\").replace('"', "\\\""))
}

fn normalize_markdown_one_line(value: &str) -> String {
    value.split_whitespace().collect::<Vec<_>>().join(" ")
}

fn truncate_markdown_one_line(value: &str, max_chars: usize) -> String {
    let normalized = normalize_markdown_one_line(value);
    if normalized.chars().count() <= max_chars {
        return normalized;
    }
    let keep = max_chars.saturating_sub(3);
    let mut truncated = normalized.chars().take(keep).collect::<String>();
    truncated.push_str("...");
    truncated
}

fn truncate_at_a_glance_step(value: &str) -> String {
    truncate_markdown_one_line(&value.replace('`', ""), 80)
}

fn markdown_inline_code(value: &str) -> String {
    format!("`{}`", value.replace('`', "'"))
}

fn retained_artifact_extension(path: &str) -> String {
    Path::new(path)
        .extension()
        .and_then(|ext| ext.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase()
}

fn retained_artifact_file_name(path: &str) -> String {
    Path::new(path)
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or(path)
        .to_string()
}

fn retained_artifact_render_rank(path: &str) -> usize {
    match retained_artifact_extension(path).as_str() {
        "svg" => 0,
        "png" => 1,
        "json" => 2,
        "csv" => 3,
        "txt" | "md" => 4,
        _ => 5,
    }
}

fn retained_artifact_preview(path: &Path) -> Option<String> {
    let text = fs::read_to_string(path).ok()?;
    let extension = path
        .extension()
        .and_then(|ext| ext.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase();
    if extension == "json" {
        if let Ok(value) = serde_json::from_str::<serde_json::Value>(&text) {
            if let Some(schema) = value.get("schema").and_then(|schema| schema.as_str()) {
                return Some(format!("schema: {}", markdown_inline_code(schema)));
            }
        }
    }
    let first_line = text.lines().find(|line| !line.trim().is_empty())?.trim();
    if first_line.is_empty() {
        return None;
    }
    Some(markdown_inline_code(&truncate_markdown_one_line(
        first_line, 120,
    )))
}

fn decode_basic_xml_entities(text: &str) -> String {
    text.replace("&lt;", "<")
        .replace("&gt;", ">")
        .replace("&quot;", "\"")
        .replace("&apos;", "'")
        .replace("&amp;", "&")
}

fn retained_svg_text_preview(path: &Path) -> Option<String> {
    let text = fs::read_to_string(path).ok()?;
    let mut labels = Vec::new();
    let mut search_start = 0usize;
    while let Some(relative_text_start) = text[search_start..].find("<text") {
        let text_start = search_start + relative_text_start;
        let Some(tag_close_relative) = text[text_start..].find('>') else {
            break;
        };
        let tag = &text[text_start..text_start + tag_close_relative + 1];
        let value_start = text_start + tag_close_relative + 1;
        let Some(value_end_relative) = text[value_start..].find("</text>") else {
            search_start = value_start;
            continue;
        };
        let value_end = value_start + value_end_relative;
        if tag.contains("tfbs-score-track-logo-letter") {
            search_start = value_end + "</text>".len();
            continue;
        }
        let label = decode_basic_xml_entities(text[value_start..value_end].trim());
        let label = label.split_whitespace().collect::<Vec<_>>().join(" ");
        if !label.is_empty()
            && !matches!(label.as_str(), "A" | "C" | "G" | "T")
            && !labels
                .iter()
                .any(|existing: &String| existing.eq_ignore_ascii_case(&label))
        {
            labels.push(label);
        }
        if labels.len() >= 8 {
            break;
        }
        search_start = value_end + "</text>".len();
    }
    if labels.is_empty() {
        return None;
    }
    Some(markdown_inline_code(&truncate_markdown_one_line(
        &labels.join(" | "),
        180,
    )))
}

fn render_tutorial_front_matter(
    chapter: &TutorialChapter,
    loaded: &LoadedWorkflowExample,
    source_path: &str,
    executed: bool,
    review_entry: Option<&TutorialReviewEntry>,
    review_stale_reason: Option<&str>,
    generated_artifact_dir: &str,
) -> String {
    let mut out = String::new();
    out.push_str("---\n");
    out.push_str("chapter_id: ");
    out.push_str(&yaml_double_quote(&chapter.id));
    out.push('\n');
    out.push_str("title: ");
    out.push_str(&yaml_double_quote(&chapter.title));
    out.push('\n');
    out.push_str("tier: ");
    out.push_str(&yaml_double_quote(chapter.tier.as_str()));
    out.push('\n');
    out.push_str("example_id: ");
    out.push_str(&yaml_double_quote(&chapter.example_id));
    out.push('\n');
    out.push_str("source_example: ");
    out.push_str(&yaml_double_quote(source_path));
    out.push('\n');
    out.push_str("example_test_mode: ");
    out.push_str(&yaml_double_quote(loaded.example.test_mode.as_str()));
    out.push('\n');
    out.push_str("executed_during_generation: ");
    out.push_str(if executed { "true" } else { "false" });
    out.push('\n');
    out.push_str("automated_status: ");
    out.push_str(&yaml_double_quote(&tutorial_automated_status(
        chapter, executed,
    )));
    out.push('\n');
    out.push_str("review_status: ");
    out.push_str(&yaml_double_quote(&tutorial_review_status(review_entry)));
    out.push('\n');
    out.push_str("review_stale: ");
    let review_projection = tutorial_review_projection(
        review_entry,
        current_utc_review_date(),
        default_tutorial_review_warn_after_months(),
        review_stale_reason.map(str::to_string),
    );
    out.push_str(if review_projection.stale {
        "true"
    } else {
        "false"
    });
    out.push('\n');
    push_yaml_optional_string(
        &mut out,
        "codex_reviewed_at",
        review_entry.and_then(|entry| entry.codex_reviewed_at.as_deref()),
    );
    push_yaml_optional_string(
        &mut out,
        "human_reviewed_at",
        review_entry.and_then(|entry| entry.human_reviewed_at.as_deref()),
    );
    push_yaml_optional_string(
        &mut out,
        "human_reviewer",
        review_entry.and_then(|entry| entry.human_reviewer.as_deref()),
    );
    push_yaml_optional_string(
        &mut out,
        "review_stale_reason",
        review_projection.stale_reason.as_deref(),
    );
    push_yaml_optional_string(
        &mut out,
        "review_issue_template",
        review_projection.issue_template.as_deref(),
    );
    push_yaml_optional_string(
        &mut out,
        "review_issue_template_path",
        review_projection.issue_template_path.as_deref(),
    );
    out.push_str("generated_artifact_dir: ");
    out.push_str(&yaml_double_quote(generated_artifact_dir));
    out.push('\n');
    out.push_str("---\n\n");
    out
}

fn push_yaml_optional_string(out: &mut String, key: &str, value: Option<&str>) {
    out.push_str(key);
    out.push_str(": ");
    if let Some(value) = value {
        out.push_str(&yaml_double_quote(value));
    } else {
        out.push_str("null");
    }
    out.push('\n');
}

fn tutorial_automated_status(chapter: &TutorialChapter, executed: bool) -> String {
    if executed {
        "passing".to_string()
    } else if chapter.tier == TutorialTier::Online {
        "skipped_online".to_string()
    } else {
        "skipped".to_string()
    }
}

fn tutorial_review_status(entry: Option<&TutorialReviewEntry>) -> String {
    let Some(entry) = entry else {
        return "missing_review_manifest_entry".to_string();
    };
    if entry.tutorial_status == "deprecated" {
        return "deprecated".to_string();
    }
    if entry
        .replaced_by
        .as_deref()
        .map(str::trim)
        .is_some_and(|value| !value.is_empty())
    {
        return "replaced".to_string();
    }
    if entry
        .human_reviewed_at
        .as_deref()
        .map(str::trim)
        .is_some_and(|value| !value.is_empty())
    {
        return "human_reviewed".to_string();
    }
    if entry
        .codex_reviewed_at
        .as_deref()
        .map(str::trim)
        .is_some_and(|value| !value.is_empty())
    {
        return "codex_reviewed".to_string();
    }
    "unreviewed".to_string()
}

fn tutorial_review_projection(
    entry: Option<&TutorialReviewEntry>,
    today: Option<(i32, u32, u32)>,
    warn_after_months: u32,
    dependency_stale_reason: Option<String>,
) -> TutorialReviewProjection {
    let status = tutorial_review_status(entry);
    let age_stale_reason = entry
        .and_then(|entry| entry.human_reviewed_at.as_deref())
        .zip(today)
        .and_then(|(reviewed_at, today)| {
            review_date_is_stale(reviewed_at, today, warn_after_months).then(|| {
                format!(
                    "human review date {} is older than {} months",
                    reviewed_at, warn_after_months
                )
            })
        });
    let stale_reason = dependency_stale_reason.or(age_stale_reason);
    let stale = stale_reason.is_some();
    let issue_template_path =
        tutorial_review_issue_template_for_projection(&status, stale_reason.as_deref());
    TutorialReviewProjection {
        status,
        codex_reviewed_at: entry.and_then(|entry| entry.codex_reviewed_at.clone()),
        human_reviewed_at: entry.and_then(|entry| entry.human_reviewed_at.clone()),
        human_reviewer: entry.and_then(|entry| entry.human_reviewer.clone()),
        stale,
        stale_reason,
        issue_template: issue_template_path
            .map(tutorial_issue_template_name)
            .map(str::to_string),
        issue_template_path: issue_template_path.map(str::to_string),
    }
}

fn tutorial_review_issue_template_for_projection(
    status: &str,
    stale_reason: Option<&str>,
) -> Option<&'static str> {
    if matches!(status, "deprecated" | "replaced") {
        return None;
    }
    if let Some(reason) = stale_reason {
        let lower = reason.to_ascii_lowercase();
        if lower.contains("declared graphic")
            || lower.contains(".svg")
            || lower.contains(".png")
            || lower.contains(".jpg")
            || lower.contains(".jpeg")
        {
            return Some(TUTORIAL_ARTIFACT_ISSUE_TEMPLATE_PATH);
        }
        return Some(TUTORIAL_CONFUSION_ISSUE_TEMPLATE_PATH);
    }
    if matches!(status, "unreviewed" | "missing_review_manifest_entry") {
        return Some(TUTORIAL_CONFUSION_ISSUE_TEMPLATE_PATH);
    }
    None
}

fn tutorial_review_badge_markdown(
    status: &str,
    stale: bool,
    codex_reviewed_at: Option<&str>,
    human_reviewed_at: Option<&str>,
    human_reviewer: Option<&str>,
    issue_template_path: Option<&str>,
    issue_link_prefix: &str,
) -> String {
    let mut badge = format!("review `{status}`");
    if stale {
        badge.push_str(" `stale`");
    }
    if let Some(human_reviewed_at) = human_reviewed_at
        .map(str::trim)
        .filter(|value| !value.is_empty())
    {
        badge.push_str(" - human ");
        badge.push_str(human_reviewed_at);
        if let Some(reviewer) = human_reviewer
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            badge.push_str(" by ");
            badge.push_str(reviewer);
        }
    } else if let Some(codex_reviewed_at) = codex_reviewed_at
        .map(str::trim)
        .filter(|value| !value.is_empty())
    {
        badge.push_str(" - codex ");
        badge.push_str(codex_reviewed_at);
    }
    if let Some(path) = issue_template_path
        .map(str::trim)
        .filter(|value| !value.is_empty())
    {
        badge.push_str(" - [file feedback](");
        badge.push_str(issue_link_prefix);
        badge.push_str(path);
        badge.push(')');
    }
    badge
}

pub fn tutorial_review_badge_label(
    status: Option<&str>,
    stale: bool,
    codex_reviewed_at: Option<&str>,
    human_reviewed_at: Option<&str>,
    human_reviewer: Option<&str>,
) -> String {
    let status = status
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .unwrap_or("missing_review_manifest_entry");
    let mut badge = format!("review: {status}");
    if stale {
        badge.push_str(", stale");
    }
    if let Some(human_reviewed_at) = human_reviewed_at
        .map(str::trim)
        .filter(|value| !value.is_empty())
    {
        badge.push_str(", human ");
        badge.push_str(human_reviewed_at);
        if let Some(reviewer) = human_reviewer
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            badge.push_str(" by ");
            badge.push_str(reviewer);
        }
    } else if let Some(codex_reviewed_at) = codex_reviewed_at
        .map(str::trim)
        .filter(|value| !value.is_empty())
    {
        badge.push_str(", codex ");
        badge.push_str(codex_reviewed_at);
    }
    badge
}

pub fn build_tutorial_feedback_context_text(
    entry: &TutorialCatalogEntry,
    manifest: Option<&TutorialManifest>,
    review_manifest: Option<&TutorialReviewManifest>,
    tutorial_source_dir: &Path,
    generated_output_dir: &Path,
    workflow_example_dir: &Path,
    current_search_section: Option<&str>,
    gentle_version: &str,
    platform: &str,
) -> String {
    let review_by_id = review_manifest
        .map(|manifest| {
            manifest
                .entries
                .iter()
                .map(|entry| (entry.tutorial_id.as_str(), entry))
                .collect::<HashMap<_, _>>()
        })
        .unwrap_or_default();
    let review_entry = review_by_id.get(entry.id.as_str()).copied();
    let chapter = manifest.and_then(|manifest| {
        manifest
            .chapters
            .iter()
            .find(|chapter| chapter.id == entry.id)
    });
    let source_json = tutorial_source_path_for_id_lossy(tutorial_source_dir, &entry.id);
    let workflow_path = chapter
        .map(|chapter| {
            workflow_example_path_for_id_lossy(workflow_example_dir, &chapter.example_id)
        })
        .unwrap_or_else(|| "not applicable".to_string());
    let generated_artifact_dir = chapter
        .map(|chapter| {
            display_path(
                &generated_output_dir
                    .join("artifacts")
                    .join(markdown_file_stem(&chapter.id)),
            )
        })
        .unwrap_or_else(|| "not applicable".to_string());
    let tutorial_kind = chapter
        .map(|_| "generated_chapter")
        .unwrap_or("hand_written_or_reference");
    let mut out = String::new();
    out.push_str("Tutorial feedback context\n");
    out.push_str("-------------------------\n");
    out.push_str("Tutorial title: ");
    out.push_str(&entry.title);
    out.push('\n');
    out.push_str("Tutorial id: ");
    out.push_str(&entry.id);
    out.push('\n');
    out.push_str("Tutorial kind: ");
    out.push_str(tutorial_kind);
    out.push('\n');
    out.push_str("Catalog path: ");
    out.push_str(&entry.path);
    out.push('\n');
    out.push_str("Catalog source: ");
    out.push_str(&entry.source);
    out.push('\n');
    out.push_str("Current search section: ");
    out.push_str(
        current_search_section
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .unwrap_or("none"),
    );
    out.push('\n');
    out.push_str("GENtle version: ");
    out.push_str(gentle_version);
    out.push('\n');
    out.push_str("Platform: ");
    out.push_str(platform);
    out.push('\n');
    out.push_str("Source JSON path: ");
    out.push_str(&source_json);
    out.push('\n');
    out.push_str("Workflow path: ");
    out.push_str(&workflow_path);
    out.push('\n');
    out.push_str("Generated artifact dir: ");
    out.push_str(&generated_artifact_dir);
    out.push('\n');
    out.push_str("Review status: ");
    out.push_str(&tutorial_review_status(review_entry));
    out.push('\n');
    out.push_str("Review stale: ");
    out.push_str(if entry.review_stale { "yes" } else { "no" });
    out.push('\n');
    out.push_str("Review stale reason: ");
    out.push_str(
        entry
            .review_stale_reason
            .as_deref()
            .unwrap_or("not recorded"),
    );
    out.push('\n');
    out.push_str("Suggested issue template: ");
    out.push_str(
        entry
            .review_issue_template_path
            .as_deref()
            .unwrap_or("not recorded"),
    );
    out.push('\n');
    out.push_str("Codex reviewed at: ");
    out.push_str(
        review_entry
            .and_then(|entry| entry.codex_reviewed_at.as_deref())
            .unwrap_or("not recorded"),
    );
    out.push('\n');
    out.push_str("Human reviewed at: ");
    out.push_str(
        review_entry
            .and_then(|entry| entry.human_reviewed_at.as_deref())
            .unwrap_or("not recorded"),
    );
    out.push('\n');
    out.push_str("Interface used: GUI / CLI / Agent Assistant / ClawBio\n");
    out.push_str("Step reached:\n");
    out.push_str("Expected vs. actual:\n");
    out
}

fn render_tutorial_concepts_compact(
    chapter: &TutorialChapter,
    concept_by_id: &HashMap<String, TutorialConcept>,
) -> Result<String, String> {
    let mut out = String::new();
    out.push_str("\n## Applied Concepts\n\n");
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
    }
    Ok(out)
}

fn render_tutorial_parameters_that_matter(chapter: &TutorialChapter) -> String {
    let mut out = String::new();
    out.push_str("\n## Parameters That Matter\n\n");
    if chapter.parameter_notes.is_empty() {
        out.push_str("- This chapter intentionally avoids additional command options; run the canonical workflow unchanged.\n");
        return out;
    }
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
    out
}

fn render_tutorial_at_a_glance(chapter: &TutorialChapter) -> String {
    if chapter.gui_steps.len() < 4 {
        return String::new();
    }
    let mut out = String::new();
    out.push_str("\n## At a Glance\n\n");
    for (idx, step) in chapter.gui_steps.iter().enumerate() {
        out.push_str(&format!(
            "{}. {}\n",
            idx + 1,
            truncate_at_a_glance_step(step)
        ));
    }
    out
}

fn render_tutorial_prerequisites(
    chapter: &TutorialChapter,
    chapter_by_id: &HashMap<String, TutorialChapter>,
) -> Result<String, String> {
    if chapter.prerequisites.is_empty() {
        return Ok(String::new());
    }
    let mut out = String::new();
    out.push_str("\n**Prerequisites:** Read ");
    for (idx, prerequisite_id) in chapter.prerequisites.iter().enumerate() {
        let prerequisite = chapter_by_id.get(prerequisite_id).ok_or_else(|| {
            format!(
                "Tutorial chapter '{}' references unknown prerequisite '{}'",
                chapter.id, prerequisite_id
            )
        })?;
        if idx > 0 {
            out.push_str(", ");
        }
        out.push('[');
        out.push_str(&chapter_link_label(prerequisite));
        out.push_str("](");
        out.push_str(&chapter_link_target(prerequisite, false));
        out.push(')');
    }
    out.push_str(" first.");
    out.push_str("\n");
    Ok(out)
}

fn render_tutorial_local_execution_note(chapter: &TutorialChapter) -> String {
    let Some(note) = &chapter.local_execution_note else {
        return String::new();
    };
    let note = note.trim();
    if note.is_empty() {
        return String::new();
    }
    let mut out = String::new();
    out.push_str("\n> **How to Run This Locally**\n");
    for line in note.lines() {
        out.push_str("> ");
        out.push_str(line.trim());
        out.push('\n');
    }
    out
}

fn guided_walkthrough_target_from_generated_chapter(guide_path: &str) -> String {
    let guide_path = guide_path.trim();
    guide_path
        .strip_prefix("docs/tutorial/")
        .map(|relative| format!("../../{relative}"))
        .unwrap_or_else(|| guide_path.to_string())
}

fn render_tutorial_guided_walkthrough_see_also(chapter: &TutorialChapter) -> String {
    let Some(guide_path) = &chapter.guide_path else {
        return String::new();
    };
    let guide_path = guide_path.trim();
    if guide_path.is_empty() {
        return String::new();
    }
    let mut out = String::new();
    out.push_str("\nSee also: guided walkthrough [");
    out.push_str(guide_path);
    out.push_str("](");
    out.push_str(&guided_walkthrough_target_from_generated_chapter(
        guide_path,
    ));
    out.push_str("). Use that page first when you want a human-led path; this chapter is the executable reference.\n");
    out
}

fn tutorial_step_at(steps: &[String], idx: usize) -> Option<&str> {
    steps
        .get(idx)
        .map(|step| step.trim())
        .filter(|step| !step.is_empty())
}

fn tutorial_step_heading(gui_step: &str) -> String {
    let mut heading = gui_step.trim().trim_end_matches('.').replace('`', "");
    if heading.chars().count() > 80 {
        heading = truncate_markdown_one_line(&heading, 80);
    }
    if let Some(first) = heading.get_mut(0..1) {
        first.make_ascii_uppercase();
    }
    heading
}

fn tutorial_chapter_has_complete_cli_steps(chapter: &TutorialChapter) -> bool {
    !chapter.gui_steps.is_empty()
        && chapter
            .gui_steps
            .iter()
            .enumerate()
            .all(|(idx, _)| tutorial_step_at(&chapter.cli_steps, idx).is_some())
}

fn tutorial_graphic_link_target(graphic: &TutorialGraphic) -> String {
    let path = graphic.path.trim();
    if let Some(relative) = path.strip_prefix("docs/tutorial/generated/") {
        return format!("../{relative}");
    }
    if let Some(relative) = path.strip_prefix("docs/tutorial/") {
        return format!("../../{relative}");
    }
    if let Some(relative) = path.strip_prefix("docs/") {
        return format!("../../../{relative}");
    }
    path.to_string()
}

fn tutorial_graphic_output_path(graphic: &TutorialGraphic, output_dir: &Path) -> PathBuf {
    let path = graphic.path.trim();
    if let Some(relative) = path.strip_prefix("docs/tutorial/generated/") {
        return output_dir.join(relative);
    }
    PathBuf::from(path)
}

fn tutorial_generated_graphic_artifact_relative(graphic: &TutorialGraphic) -> Option<String> {
    graphic
        .path
        .trim()
        .strip_prefix("docs/tutorial/generated/")
        .map(ToOwned::to_owned)
}

fn render_tutorial_inline_graphic(graphic: &TutorialGraphic, output_dir: &Path) -> String {
    let mut out = String::new();
    let caption = graphic.caption.trim();
    let target = tutorial_graphic_link_target(graphic);
    let graphic_path = tutorial_graphic_output_path(graphic, output_dir);
    out.push_str("![");
    out.push_str(caption);
    out.push_str("](");
    out.push_str(&target);
    out.push_str(")\n\n");
    out.push_str("*Figure: ");
    out.push_str(caption);
    if let Some(command) = graphic
        .regen_command
        .as_deref()
        .map(str::trim)
        .filter(|v| !v.is_empty())
    {
        out.push_str(" Regenerate with `");
        out.push_str(command);
        out.push_str("`.");
    }
    if let Some(capture_date) = graphic
        .capture_date
        .as_deref()
        .map(str::trim)
        .filter(|v| !v.is_empty())
    {
        out.push_str(" Screenshot captured ");
        out.push_str(capture_date);
        out.push('.');
    }
    out.push_str("*\n\n");
    if retained_artifact_extension(&graphic.path) == "svg" {
        if let Some(preview) = retained_svg_text_preview(&graphic_path) {
            out.push_str("> SVG text labels: ");
            out.push_str(&preview);
            out.push_str(". If the embedded preview omits text in the GUI, open the linked SVG or use these labels as the figure legend.\n\n");
        }
    }
    out
}

fn render_tutorial_step_graphics(
    chapter: &TutorialChapter,
    step_number: usize,
    output_dir: &Path,
) -> String {
    let mut graphics = chapter
        .graphics
        .iter()
        .filter(|graphic| graphic.illustrates_step == step_number)
        .collect::<Vec<_>>();
    graphics.sort_by(|left, right| left.path.cmp(&right.path));
    let mut out = String::new();
    for graphic in graphics {
        out.push_str(&render_tutorial_inline_graphic(graphic, output_dir));
    }
    out
}

fn render_tutorial_gui_steps(chapter: &TutorialChapter, output_dir: &Path) -> String {
    let mut out = String::new();
    out.push_str("\n## GUI First\n\n");
    if chapter.gui_steps.is_empty() {
        out.push_str("- No GUI-first steps are required for this chapter.\n");
        return out;
    }
    if !chapter.cli_steps.is_empty() {
        out.push_str("CLI snippets use GENtle's default `.gentle_state.json` state unless they say otherwise. Add `--state PATH` or `--project PATH` when you want an explicit sandboxed state file for copied commands.\n\n");
    }
    for (idx, step) in chapter.gui_steps.iter().enumerate() {
        let step = step.trim();
        let cli_step = tutorial_step_at(&chapter.cli_steps, idx);
        let expectation = tutorial_step_at(&chapter.step_expectations, idx);
        if cli_step.is_none() && expectation.is_none() {
            out.push_str(&format!("{}. ", idx + 1));
            out.push_str(step);
            out.push('\n');
            out.push_str(&render_tutorial_step_graphics(chapter, idx + 1, output_dir));
            continue;
        }
        out.push_str(&format!(
            "### Step {}: {}\n\n",
            idx + 1,
            tutorial_step_heading(step)
        ));
        out.push_str("GUI: ");
        out.push_str(step);
        out.push_str("\n\n");
        if let Some(cli_step) = cli_step {
            out.push_str("CLI:\n\n");
            out.push_str("```bash\n");
            out.push_str(cli_step);
            out.push_str("\n```\n\n");
        }
        if let Some(expectation) = expectation {
            out.push_str("> Expected: ");
            out.push_str(expectation);
            out.push_str("\n\n");
        }
        out.push_str(&render_tutorial_step_graphics(chapter, idx + 1, output_dir));
    }
    out
}

fn tutorial_inline_artifact_steps(graphics: &[TutorialGraphic]) -> HashMap<String, usize> {
    let mut inline_steps = HashMap::new();
    for graphic in graphics {
        if let Some(relative) = tutorial_generated_graphic_artifact_relative(graphic) {
            inline_steps.insert(relative, graphic.illustrates_step);
        }
    }
    inline_steps
}

fn render_tutorial_produced_artifacts(
    retained_artifacts: &[String],
    output_dir: &Path,
    graphics: &[TutorialGraphic],
) -> String {
    if retained_artifacts.is_empty() {
        return String::new();
    }
    let mut out = String::new();
    out.push_str("\n## What This Chapter Produces\n\n");
    let mut ordered = retained_artifacts.to_vec();
    let inline_steps = tutorial_inline_artifact_steps(graphics);
    ordered.sort_by(|a, b| {
        retained_artifact_render_rank(a)
            .cmp(&retained_artifact_render_rank(b))
            .then_with(|| a.cmp(b))
    });
    for retained in &ordered {
        let artifact_path = output_dir.join(retained);
        let link_target = format!("../{}", retained);
        let file_name = retained_artifact_file_name(retained);
        let extension = retained_artifact_extension(retained);
        let inline_step = inline_steps.get(retained);
        if matches!(extension.as_str(), "svg" | "png") && artifact_path.is_file() {
            out.push_str("- [`");
            out.push_str(retained);
            out.push_str("`](");
            out.push_str(&link_target);
            out.push_str(")\n\n");
            if let Some(step) = inline_step {
                out.push_str("  - Embedded above near Step ");
                out.push_str(&step.to_string());
                out.push_str("; kept here as an audit link.\n\n");
            }
            if extension == "svg" {
                if let Some(preview) = retained_svg_text_preview(&artifact_path) {
                    out.push_str("> SVG text labels: ");
                    out.push_str(&preview);
                    out.push_str(". If this embedded preview omits text in the GUI, open the linked SVG or use these labels as the figure legend.\n\n");
                }
            }
            if inline_step.is_some() {
                continue;
            }
            out.push_str("![");
            out.push_str(&file_name);
            out.push_str("](");
            out.push_str(&link_target);
            out.push_str(")\n\n");
            continue;
        }
        out.push_str("- [`");
        out.push_str(retained);
        out.push_str("`](");
        out.push_str(&link_target);
        out.push(')');
        if let Some(preview) = retained_artifact_preview(&artifact_path) {
            out.push_str(" - ");
            out.push_str(&preview);
        }
        out.push('\n');
    }
    out
}

fn render_tutorial_provenance_section(
    chapter: &TutorialChapter,
    loaded: &LoadedWorkflowExample,
    executed: bool,
    online_enabled: bool,
    tutorial_source_path: &str,
    workflow_path: &str,
    generated_artifact_dir: &str,
    review_entry: Option<&TutorialReviewEntry>,
) -> String {
    let mut out = String::new();
    out.push_str("\n## Tutorial Provenance\n\n");
    out.push_str("- Chapter id: `");
    out.push_str(&chapter.id);
    out.push_str("`\n");
    out.push_str("- Tier: `");
    out.push_str(chapter.tier.as_str());
    out.push_str("`\n");
    out.push_str("- Example id: `");
    out.push_str(&chapter.example_id);
    out.push_str("`\n");
    out.push_str("- Tutorial source JSON: `");
    out.push_str(tutorial_source_path);
    out.push_str("`\n");
    out.push_str("- Workflow file: `");
    out.push_str(workflow_path);
    out.push_str("`\n");
    out.push_str("- Generated artifact dir: `");
    out.push_str(generated_artifact_dir);
    out.push_str("`\n");
    out.push_str("- Example test_mode: `");
    out.push_str(loaded.example.test_mode.as_str());
    out.push_str("`\n");
    out.push_str("- Executed during generation: `");
    out.push_str(if executed { "yes" } else { "no" });
    out.push_str("`\n");
    out.push_str("- Automated status: `");
    out.push_str(&tutorial_automated_status(chapter, executed));
    out.push_str("`\n");
    out.push_str("- Review status: `");
    out.push_str(&tutorial_review_status(review_entry));
    out.push_str("`\n");
    out.push_str("- Codex reviewed at: `");
    out.push_str(
        review_entry
            .and_then(|entry| entry.codex_reviewed_at.as_deref())
            .unwrap_or("not recorded"),
    );
    out.push_str("`\n");
    out.push_str("- Human reviewed at: `");
    out.push_str(
        review_entry
            .and_then(|entry| entry.human_reviewed_at.as_deref())
            .unwrap_or("not recorded"),
    );
    out.push_str("`\n");
    if chapter.tier == TutorialTier::Online && !online_enabled {
        out.push_str(
            "- Execution note: set `GENTLE_TEST_ONLINE=1` before `tutorial-generate` to execute this chapter.\n",
        );
    }
    out.push_str("- Inspect the source JSON when you need full option-level detail.\n");
    out
}

fn render_tutorial_feedback_section(chapter: &TutorialChapter) -> String {
    let mut out = String::new();
    out.push_str("\n## Feedback\n\n");
    out.push_str("If this tutorial is confusing, execution-stale, biologically suspect, or missing a useful figure, please open the matching tutorial issue template and include the context below.\n\n");
    out.push_str("- Tutorial title: `");
    out.push_str(&chapter.title);
    out.push_str("`\n");
    out.push_str("- Tutorial/chapter id: `");
    out.push_str(&chapter.id);
    out.push_str("`\n");
    out.push_str("- Step reached:\n");
    out.push_str("- Expected vs. actual:\n");
    out.push_str("- Interface used: GUI / CLI / Agent Assistant / ClawBio\n\n");
    out.push_str("Paste the Tutorial feedback context here:\n\n");
    out.push_str("```text\n\n```\n");
    out
}

fn render_tutorial_chapter_markdown(
    chapter: &TutorialChapter,
    loaded: &LoadedWorkflowExample,
    executed: bool,
    retained_artifacts: &[String],
    output_dir: &Path,
    provenance_output_dir: &Path,
    tutorial_source_dir: &Path,
    online_enabled: bool,
    concept_by_id: &HashMap<String, TutorialConcept>,
    chapter_by_id: &HashMap<String, TutorialChapter>,
    review_entry: Option<&TutorialReviewEntry>,
    review_stale_reason: Option<&str>,
) -> Result<String, String> {
    let workflow_path = display_path(&loaded.path);
    let tutorial_source_path = tutorial_source_path_for_id_lossy(tutorial_source_dir, &chapter.id);
    let generated_artifact_dir = display_path(
        &provenance_output_dir
            .join("artifacts")
            .join(markdown_file_stem(&chapter.id)),
    );
    let mut out = String::new();
    out.push_str(&render_tutorial_front_matter(
        chapter,
        loaded,
        &workflow_path,
        executed,
        review_entry,
        review_stale_reason,
        &generated_artifact_dir,
    ));
    out.push_str("# ");
    out.push_str(&chapter.title);
    out.push_str("\n\n");
    if !chapter.summary.trim().is_empty() {
        out.push_str(chapter.summary.trim());
        out.push_str("\n");
    }
    if !chapter.narrative.trim().is_empty() {
        out.push('\n');
        out.push_str(chapter.narrative.trim());
        out.push_str("\n");
    }
    out.push_str(&render_tutorial_guided_walkthrough_see_also(chapter));
    out.push_str(&render_tutorial_prerequisites(chapter, chapter_by_id)?);
    out.push_str(&render_tutorial_local_execution_note(chapter));
    out.push_str(&render_tutorial_parameters_that_matter(chapter));
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
    out.push_str(&render_tutorial_concepts_compact(chapter, concept_by_id)?);
    out.push_str(&render_tutorial_at_a_glance(chapter));
    out.push_str(&render_tutorial_gui_steps(chapter, output_dir));
    if !tutorial_chapter_has_complete_cli_steps(chapter) {
        out.push_str("\n## Command Equivalent (After GUI)\n\n");
        out.push_str("Run the same routine non-interactively once the GUI flow is clear:\n\n");
        out.push_str("```bash\n");
        out.push_str("cargo run --bin gentle_cli -- workflow @");
        out.push_str(&workflow_path);
        out.push_str("\n");
        out.push_str("cargo run --bin gentle_cli -- shell 'workflow @");
        out.push_str(&workflow_path);
        out.push_str("'\n");
        out.push_str("```\n");
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
    out.push_str(&render_tutorial_produced_artifacts(
        retained_artifacts,
        output_dir,
        &chapter.graphics,
    ));
    out.push_str(&render_tutorial_provenance_section(
        chapter,
        loaded,
        executed,
        online_enabled,
        &tutorial_source_path,
        &workflow_path,
        &generated_artifact_dir,
        review_entry,
    ));
    out.push_str(&render_tutorial_feedback_section(chapter));
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
    out.push_str("These chapters are generated from executable workflow JSON. They are stable, regenerable, and machine-checked, but they are reference material - not first-time walkthroughs. If you are new to GENtle, start with the guided walkthroughs in [`../README.md`](../README.md).\n\n");
    out.push_str("Each chapter row includes a review badge from `../review_manifest.json`; `review unreviewed` means the workflow is machine-checked but still needs a readability/biology pass from Codex or a human reviewer.\n\n");
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
    ordered.sort_by(|a, b| {
        a.group_order
            .unwrap_or(usize::MAX)
            .cmp(&b.group_order.unwrap_or(usize::MAX))
            .then_with(|| {
                a.group_position
                    .unwrap_or(usize::MAX)
                    .cmp(&b.group_position.unwrap_or(usize::MAX))
            })
            .then_with(|| a.order.cmp(&b.order))
            .then_with(|| a.id.cmp(&b.id))
    });
    let mut last_group: Option<String> = None;
    for chapter in &ordered {
        let group = chapter
            .group_label
            .as_deref()
            .filter(|value| !value.trim().is_empty())
            .unwrap_or("Ungrouped");
        if last_group.as_deref() != Some(group) {
            if last_group.is_some() {
                out.push('\n');
            }
            out.push_str("### ");
            out.push_str(group);
            out.push_str("\n\n");
            last_group = Some(group.to_string());
        }
        let chapter_file = chapter_markdown_filename_from_parts(
            chapter.order,
            &chapter.id,
            chapter.decimal_id.as_deref(),
        );
        out.push_str("- ");
        if let Some(decimal_id) = chapter.decimal_id.as_deref() {
            out.push_str("`");
            out.push_str(decimal_id);
            out.push_str("` ");
        } else {
            out.push_str(&chapter.order.to_string());
            out.push_str(". ");
        }
        out.push_str("[");
        out.push_str(&chapter.title);
        out.push_str("](./chapters/");
        out.push_str(&chapter_file);
        out.push_str(") - `");
        out.push_str(&chapter.tier);
        out.push_str("` - example `");
        out.push_str(&chapter.example_id);
        out.push_str("` - executed `");
        out.push_str(if chapter.executed { "yes" } else { "no" });
        out.push_str("` - ");
        out.push_str(&tutorial_review_badge_markdown(
            &chapter.review_status,
            chapter.review_stale,
            chapter.codex_reviewed_at.as_deref(),
            chapter.human_reviewed_at.as_deref(),
            chapter.human_reviewer.as_deref(),
            chapter.review_issue_template_path.as_deref(),
            "../../../",
        ));
        out.push('\n');
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

fn summarize_path_list(paths: Vec<String>) -> String {
    const MAX_PATHS: usize = 12;
    if paths.is_empty() {
        return "none".to_string();
    }
    let total = paths.len();
    let mut listed: Vec<String> = paths.into_iter().take(MAX_PATHS).collect();
    if total > MAX_PATHS {
        listed.push(format!("... and {} more", total - MAX_PATHS));
    }
    listed.join(", ")
}

fn tutorial_generated_file_count_mismatch_message(
    expected: &BTreeMap<String, Vec<u8>>,
    actual: &BTreeMap<String, Vec<u8>>,
) -> String {
    let generated_only: Vec<String> = actual
        .keys()
        .filter(|path| !expected.contains_key(*path))
        .cloned()
        .collect();
    let committed_only: Vec<String> = expected
        .keys()
        .filter(|path| !actual.contains_key(*path))
        .cloned()
        .collect();
    format!(
        "Tutorial generated file count mismatch: expected {}, got {}; generated-only files: {}; committed-only files: {}",
        expected.len(),
        actual.len(),
        summarize_path_list(generated_only),
        summarize_path_list(committed_only)
    )
}

fn tutorial_source_path_lookup_lossy(source_dir: &Path) -> HashMap<String, String> {
    let mut lookup = HashMap::new();
    let Ok(entries) = fs::read_dir(source_dir) else {
        return lookup;
    };
    for entry in entries.flatten() {
        let path = entry.path();
        if !path
            .extension()
            .and_then(|value| value.to_str())
            .map(|value| value.eq_ignore_ascii_case("json"))
            .unwrap_or(false)
        {
            continue;
        }
        let Ok(raw) = fs::read_to_string(&path) else {
            continue;
        };
        let Ok(unit) = serde_json::from_str::<TutorialSourceUnit>(&raw) else {
            continue;
        };
        lookup.insert(unit.id, display_path(&path));
    }
    lookup
}

fn tutorial_source_path_lookup(source_dir: &Path) -> HashMap<String, PathBuf> {
    let mut lookup = HashMap::new();
    let Ok(entries) = fs::read_dir(source_dir) else {
        return lookup;
    };
    for entry in entries.flatten() {
        let path = entry.path();
        if !path
            .extension()
            .and_then(|value| value.to_str())
            .map(|value| value.eq_ignore_ascii_case("json"))
            .unwrap_or(false)
        {
            continue;
        }
        let Ok(raw) = fs::read_to_string(&path) else {
            continue;
        };
        let Ok(unit) = serde_json::from_str::<TutorialSourceUnit>(&raw) else {
            continue;
        };
        lookup.insert(unit.id, path);
    }
    lookup
}

fn tutorial_source_path_for_id_lossy(source_dir: &Path, tutorial_id: &str) -> String {
    tutorial_source_path_lookup_lossy(source_dir)
        .remove(tutorial_id)
        .unwrap_or_else(|| display_path(&source_dir.join(format!("{tutorial_id}.json"))))
}

fn workflow_example_path_for_id_lossy(example_dir: &Path, example_id: &str) -> String {
    let mut lookup = HashMap::new();
    if let Ok(entries) = fs::read_dir(example_dir) {
        for entry in entries.flatten() {
            let path = entry.path();
            if !path
                .extension()
                .and_then(|value| value.to_str())
                .map(|value| value.eq_ignore_ascii_case("json"))
                .unwrap_or(false)
            {
                continue;
            }
            let Ok(raw) = fs::read_to_string(&path) else {
                continue;
            };
            let Ok(example) = serde_json::from_str::<WorkflowExample>(&raw) else {
                continue;
            };
            lookup.insert(example.id, display_path(&path));
        }
    }
    lookup
        .remove(example_id)
        .unwrap_or_else(|| display_path(&example_dir.join(format!("{example_id}.json"))))
}

fn resolve_tutorial_dependency_path(repo_root: &Path, path: impl AsRef<Path>) -> PathBuf {
    let path = path.as_ref();
    if path.is_absolute() {
        path.to_path_buf()
    } else {
        repo_root.join(path)
    }
}

fn tutorial_source_dependencies(
    tutorial_id: &str,
    source_path: Option<&PathBuf>,
    catalog_path: &str,
    example_id: Option<&str>,
    graphics: &[TutorialGraphic],
    repo_root: &Path,
) -> Vec<TutorialReviewDependency> {
    let mut dependencies = Vec::new();
    if let Some(path) = source_path {
        dependencies.push(TutorialReviewDependency {
            label: "source JSON".to_string(),
            path: path.clone(),
        });
    }
    if !catalog_path.trim().is_empty() {
        dependencies.push(TutorialReviewDependency {
            label: "tutorial Markdown".to_string(),
            path: resolve_tutorial_dependency_path(repo_root, catalog_path.trim()),
        });
    }
    if let Some(example_id) = example_id.map(str::trim).filter(|value| !value.is_empty()) {
        dependencies.push(TutorialReviewDependency {
            label: "workflow JSON".to_string(),
            path: resolve_tutorial_dependency_path(
                repo_root,
                Path::new(DEFAULT_WORKFLOW_EXAMPLE_DIR).join(format!("{example_id}.json")),
            ),
        });
    }
    for graphic in graphics {
        if graphic.path.trim().is_empty() {
            continue;
        }
        dependencies.push(TutorialReviewDependency {
            label: format!("declared graphic for tutorial '{}'", tutorial_id),
            path: resolve_tutorial_dependency_path(repo_root, graphic.path.trim()),
        });
    }
    dependencies
}

fn tutorial_review_dependency_stale_reason(
    entry: Option<&TutorialReviewEntry>,
    dependencies: &[TutorialReviewDependency],
) -> Option<String> {
    let reviewed_at = entry?.human_reviewed_at.as_deref()?.trim();
    if reviewed_at.is_empty() {
        return None;
    }
    let reviewed_date = parse_review_date(reviewed_at).ok()?;
    let mut missing_graphics = dependencies
        .iter()
        .filter(|dependency| dependency.label.starts_with("declared graphic"))
        .filter(|dependency| !dependency.path.exists())
        .map(|dependency| {
            format!(
                "{} '{}' is missing after human review date {}",
                dependency.label,
                display_path(&dependency.path),
                reviewed_at
            )
        })
        .collect::<Vec<_>>();
    missing_graphics.sort();
    if let Some(reason) = missing_graphics.into_iter().next() {
        return Some(reason);
    }
    let mut stale_dependencies = dependencies
        .iter()
        .filter_map(|dependency| {
            let modified = path_modified_review_date(&dependency.path)?;
            review_date_after(modified, reviewed_date).then(|| {
                (
                    modified,
                    format!(
                        "{} '{}' changed after human review date {}",
                        dependency.label,
                        display_path(&dependency.path),
                        reviewed_at
                    ),
                )
            })
        })
        .collect::<Vec<_>>();
    stale_dependencies
        .sort_by(|left, right| right.0.cmp(&left.0).then_with(|| left.1.cmp(&right.1)));
    stale_dependencies
        .into_iter()
        .next()
        .map(|(_, reason)| reason)
}

fn tutorial_review_dependency_stale_reasons(
    manifest_path: &Path,
    manifest: &TutorialManifest,
    review_by_id: &HashMap<String, TutorialReviewEntry>,
) -> Result<HashMap<String, String>, String> {
    let source_dir = tutorial_source_dir_for_manifest(manifest_path);
    let source_path_by_id = tutorial_source_path_lookup(&source_dir);
    let repo_root = tutorial_repo_root_for_source_dir(&source_dir);
    let mut reasons = HashMap::new();
    if source_dir.exists() {
        for unit in load_tutorial_source_units(&source_dir)? {
            let Some(reason) = tutorial_review_dependency_stale_reason(
                review_by_id.get(&unit.id),
                &tutorial_source_dependencies(
                    &unit.id,
                    source_path_by_id.get(&unit.id),
                    unit.catalog
                        .as_ref()
                        .map(|catalog| catalog.path.as_str())
                        .unwrap_or_default(),
                    unit.generated_chapter
                        .as_ref()
                        .map(|chapter| chapter.example_id.as_str()),
                    &unit.graphics,
                    &repo_root,
                ),
            ) else {
                continue;
            };
            reasons.insert(unit.id, reason);
        }
    }
    for chapter in &manifest.chapters {
        if reasons.contains_key(&chapter.id) {
            continue;
        }
        let workflow_dependency = TutorialReviewDependency {
            label: "workflow JSON".to_string(),
            path: resolve_tutorial_dependency_path(
                &repo_root,
                Path::new(DEFAULT_WORKFLOW_EXAMPLE_DIR)
                    .join(format!("{}.json", chapter.example_id)),
            ),
        };
        if let Some(reason) = tutorial_review_dependency_stale_reason(
            review_by_id.get(&chapter.id),
            &[workflow_dependency],
        ) {
            reasons.insert(chapter.id.clone(), reason);
        }
    }
    Ok(reasons)
}

fn tutorial_source_dir_for_manifest(manifest_path: &Path) -> PathBuf {
    manifest_path
        .parent()
        .map(|parent| parent.join("sources"))
        .unwrap_or_else(|| PathBuf::from(DEFAULT_TUTORIAL_SOURCE_DIR))
}

fn tutorial_chapter_for_generated_path<'a>(
    manifest: &'a TutorialManifest,
    rel_path: &str,
) -> Option<&'a TutorialChapter> {
    let rel_path = rel_path.trim();
    for chapter in &manifest.chapters {
        let chapter_file = format!("chapters/{}", markdown_path_for_chapter(chapter));
        if rel_path == chapter_file {
            return Some(chapter);
        }
        let artifact_prefix = format!("artifacts/{}/", markdown_file_stem(&chapter.id));
        if rel_path.starts_with(&artifact_prefix) {
            return Some(chapter);
        }
    }
    None
}

fn tutorial_check_issue_template_for_path(rel_path: Option<&str>) -> &'static str {
    match rel_path {
        Some(path)
            if path.starts_with("artifacts/")
                && (path.ends_with(".svg")
                    || path.ends_with(".png")
                    || path.ends_with(".jpg")
                    || path.ends_with(".jpeg")) =>
        {
            "Tutorial artifact/figure problem"
        }
        Some(path) if path.starts_with("artifacts/") => "Tutorial artifact/figure problem",
        Some(path) if path.starts_with("chapters/") || path.ends_with(".md") => {
            "Tutorial confusion"
        }
        _ => "Tutorial execution failure",
    }
}

fn extract_tutorial_chapter_id_from_error(error: &str) -> Option<String> {
    let marker = "Tutorial chapter '";
    let start = error.find(marker)? + marker.len();
    let rest = &error[start..];
    let end = rest.find('\'')?;
    let id = rest[..end].trim();
    (!id.is_empty()).then(|| id.to_string())
}

fn tutorial_check_feedback_context(
    manifest: &TutorialManifest,
    examples: &[LoadedWorkflowExample],
    tutorial_source_dir: &Path,
    manifest_path: &Path,
    expected_output_dir: &Path,
    rel_path: Option<&str>,
    chapter_id_hint: Option<&str>,
    failing_check: &str,
) -> String {
    let source_paths = tutorial_source_path_lookup_lossy(tutorial_source_dir);
    let example_by_id = example_lookup(examples);
    let chapter = rel_path
        .and_then(|path| tutorial_chapter_for_generated_path(manifest, path))
        .or_else(|| {
            chapter_id_hint.and_then(|id| manifest.chapters.iter().find(|chapter| chapter.id == id))
        });
    let issue_template = tutorial_check_issue_template_for_path(rel_path);
    let mut out = String::new();
    out.push_str("\n\nTutorial feedback context\n");
    out.push_str("-------------------------\n");
    if let Some(chapter) = chapter {
        let workflow_path = example_by_id
            .get(&chapter.example_id)
            .map(|loaded| display_path(&loaded.path))
            .unwrap_or_else(|| "(workflow example not found)".to_string());
        let source_path = source_paths
            .get(&chapter.id)
            .cloned()
            .unwrap_or_else(|| "(tutorial source JSON not found)".to_string());
        out.push_str(&format!("Chapter: {} {}\n", chapter.order, chapter.title));
        out.push_str(&format!("Tutorial id: {}\n", chapter.id));
        out.push_str(&format!("Tier: {}\n", chapter.tier.as_str()));
        out.push_str(&format!("Source JSON: {source_path}\n"));
        out.push_str(&format!("Workflow: {workflow_path}\n"));
        out.push_str(&format!(
            "Generated chapter: {}\n",
            display_path(
                &expected_output_dir
                    .join("chapters")
                    .join(markdown_path_for_chapter(chapter))
            )
        ));
        out.push_str(&format!(
            "Artifact dir: {}\n",
            display_path(
                &expected_output_dir
                    .join("artifacts")
                    .join(markdown_file_stem(&chapter.id))
            )
        ));
    } else {
        out.push_str("Chapter: (not identified)\n");
        out.push_str(&format!("Manifest: {}\n", display_path(manifest_path)));
        out.push_str(&format!(
            "Tutorial source dir: {}\n",
            display_path(tutorial_source_dir)
        ));
        out.push_str(&format!(
            "Generated tutorial dir: {}\n",
            display_path(expected_output_dir)
        ));
    }
    if let Some(path) = rel_path {
        out.push_str(&format!("Affected generated path: {path}\n"));
    }
    out.push_str(&format!("Failing check: {failing_check}\n"));
    out.push_str(&format!("Suggested issue template: {issue_template}\n"));
    out.push_str("When reporting, include what you expected, what happened, and whether you used GUI, CLI, Agent Assistant, or ClawBio.\n");
    out
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
    let review_context = tutorial_review_context(manifest_path, &manifest)?;
    let mut warnings = review_context.warnings.clone();
    let tutorial_source_dir = tutorial_source_dir_for_manifest(manifest_path);
    let provenance_output_dir = manifest_path
        .parent()
        .map(|parent| parent.join("generated"))
        .unwrap_or_else(|| PathBuf::from(DEFAULT_TUTORIAL_OUTPUT_DIR));
    let concept_by_id = tutorial_concept_lookup(&manifest.concepts)?;
    let sorted_chapters = sorted_manifest_chapters(&manifest);
    let chapter_by_id = tutorial_chapter_lookup(&sorted_chapters);
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

    for chapter in &sorted_chapters {
        let loaded = example_by_id
            .get(&chapter.example_id)
            .ok_or_else(|| {
                format!(
                    "Tutorial chapter '{}' references unknown example '{}'",
                    chapter.id, chapter.example_id
                )
            })?
            .clone();
        let mut executed = chapter.tier.should_execute(online_enabled);
        if chapter.tier == TutorialTier::Online
            && !executed
            && chapter
                .local_execution_note
                .as_deref()
                .map(str::trim)
                .unwrap_or_default()
                .is_empty()
        {
            return Err(format!(
                "Tutorial chapter '{}' is online and skipped during generation; add `local_execution_note` so readers know how to run it locally",
                chapter.id
            ));
        }
        let run_dir =
            TempDir::new().map_err(|e| format!("Could not create temp run directory: {e}"))?;
        if executed {
            let execution_result = validate_example_required_files(&loaded.example, repo_root)
                .and_then(|_| {
                    run_example_workflow_in_dir(&loaded.example, repo_root, run_dir.path())
                });
            if let Err(error) = execution_result {
                if tutorial_review_entry_exempts_execution_failure(
                    review_context.entries_by_id.get(&chapter.id),
                ) {
                    warnings.push(format!(
                        "Tutorial '{}' is deprecated or replaced and did not block generation despite execution failure: {}",
                        chapter.id, error
                    ));
                    executed = false;
                } else {
                    return Err(error);
                }
            }
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
            output_dir,
            &provenance_output_dir,
            &tutorial_source_dir,
            online_enabled,
            &concept_by_id,
            &chapter_by_id,
            review_context.entries_by_id.get(&chapter.id),
            review_context
                .stale_reasons_by_id
                .get(&chapter.id)
                .map(String::as_str),
        )?;
        let chapter_file = markdown_path_for_chapter(chapter);
        let chapter_out_path = chapters_dir.join(&chapter_file);
        fs::write(&chapter_out_path, chapter_markdown)
            .map_err(|e| format!("Could not write '{}': {e}", display_path(&chapter_out_path)))?;
        let review_projection = tutorial_review_projection(
            review_context.entries_by_id.get(&chapter.id),
            review_context.today,
            review_context.warn_after_months,
            review_context.stale_reasons_by_id.get(&chapter.id).cloned(),
        );
        chapter_reports.push(TutorialGenerationChapter {
            id: chapter.id.clone(),
            order: chapter.order,
            title: chapter.title.clone(),
            group: chapter.group.clone(),
            group_label: chapter.group_label.clone(),
            group_order: chapter.group_order,
            group_position: chapter.group_position,
            decimal_id: chapter.decimal_id.clone(),
            tier: chapter.tier.as_str().to_string(),
            example_id: chapter.example_id.clone(),
            example_source: display_path(&loaded.path),
            concepts: chapter.concepts.clone(),
            executed,
            review_status: review_projection.status,
            codex_reviewed_at: review_projection.codex_reviewed_at,
            human_reviewed_at: review_projection.human_reviewed_at,
            human_reviewer: review_projection.human_reviewer,
            review_stale: review_projection.stale,
            review_stale_reason: review_projection.stale_reason,
            review_issue_template: review_projection.issue_template,
            review_issue_template_path: review_projection.issue_template_path,
            graphics: chapter.graphics.clone(),
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
        warnings,
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
    let examples = load_workflow_examples(source_dir)?;
    let manifest = load_tutorial_manifest(manifest_path)?;
    let tutorial_source_dir = tutorial_source_dir_for_manifest(manifest_path);
    if !expected_output_dir.exists() {
        return Err(format!(
            "Expected generated tutorial directory '{}' does not exist. Run `tutorial-generate` first.",
            display_path(expected_output_dir)
        ));
    }
    let temp = TempDir::new().map_err(|e| format!("Could not create temp directory: {e}"))?;
    let generated_dir = temp.path().join("generated");
    let report = generate_tutorial_docs(source_dir, manifest_path, &generated_dir, repo_root)
        .map_err(|error| {
            let chapter_hint = extract_tutorial_chapter_id_from_error(&error);
            format!(
                "{}{}",
                error,
                tutorial_check_feedback_context(
                    &manifest,
                    &examples,
                    &tutorial_source_dir,
                    manifest_path,
                    expected_output_dir,
                    None,
                    chapter_hint.as_deref(),
                    "tutorial generation failed"
                )
            )
        })?;
    let expected = directory_bytes_map(expected_output_dir)?;
    let actual = directory_bytes_map(&generated_dir)?;
    if expected.len() != actual.len() {
        let affected_path = actual
            .keys()
            .find(|path| !expected.contains_key(*path))
            .or_else(|| expected.keys().find(|path| !actual.contains_key(*path)))
            .map(String::as_str);
        return Err(format!(
            "{}{}",
            tutorial_generated_file_count_mismatch_message(&expected, &actual),
            tutorial_check_feedback_context(
                &manifest,
                &examples,
                &tutorial_source_dir,
                manifest_path,
                expected_output_dir,
                affected_path,
                None,
                "generated tutorial file set differs from committed output"
            )
        ));
    }
    for (path, expected_bytes) in &expected {
        let actual_bytes = actual.get(path).ok_or_else(|| {
            format!(
                "{}{}",
                format!(
                    "Tutorial generated output is missing file '{}' in check run",
                    path
                ),
                tutorial_check_feedback_context(
                    &manifest,
                    &examples,
                    &tutorial_source_dir,
                    manifest_path,
                    expected_output_dir,
                    Some(path),
                    None,
                    "generated tutorial output is missing a committed file"
                )
            )
        })?;
        if expected_bytes != actual_bytes {
            return Err(format!(
                "{}{}",
                format!(
                    "Tutorial generated file '{}' differs from committed version",
                    path
                ),
                tutorial_check_feedback_context(
                    &manifest,
                    &examples,
                    &tutorial_source_dir,
                    manifest_path,
                    expected_output_dir,
                    Some(path),
                    None,
                    "generated tutorial file differs from committed version"
                )
            ));
        }
    }
    for path in actual.keys() {
        if !expected.contains_key(path) {
            return Err(format!(
                "{}{}",
                format!(
                    "Tutorial generated output contains unexpected file '{}'",
                    path
                ),
                tutorial_check_feedback_context(
                    &manifest,
                    &examples,
                    &tutorial_source_dir,
                    manifest_path,
                    expected_output_dir,
                    Some(path),
                    None,
                    "generated tutorial output contains an unexpected file"
                )
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

    fn minimal_tutorial_chapter(id: &str) -> TutorialChapter {
        TutorialChapter {
            id: id.to_string(),
            order: 1,
            title: "Minimal Tutorial".to_string(),
            group: None,
            group_label: None,
            group_order: None,
            group_position: None,
            decimal_id: None,
            example_id: "minimal_example".to_string(),
            tier: TutorialTier::Core,
            guide_path: None,
            summary: "Synthetic tutorial chapter.".to_string(),
            narrative: "Synthetic tutorial narrative.".to_string(),
            use_cases: vec!["Exercise tutorial review metadata.".to_string()],
            gui_steps: vec!["Open the tutorial.".to_string()],
            cli_steps: vec![],
            step_expectations: vec![],
            prerequisites: vec![],
            local_execution_note: None,
            learning_objectives: vec!["Understand tutorial review metadata.".to_string()],
            concepts: vec!["concept".to_string()],
            parameter_notes: vec![],
            follow_up_commands: vec![],
            retain_outputs: vec![],
            checkpoints: vec![],
            graphics: vec![],
        }
    }

    fn lock_jaspar_registry_for_test() -> std::sync::MutexGuard<'static, ()> {
        crate::tf_motifs::test_registry_lock()
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
    }

    fn first_qualifier(feature: &gb_io::seq::Feature, key: &str) -> Option<String> {
        feature.qualifier_values(key).next().map(str::to_string)
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
        let _serial = lock_jaspar_registry_for_test();
        crate::tf_motifs::reload();
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
    fn workflow_examples_isoform_protein_gel_demo_writes_svg() {
        let examples = load_workflow_examples(&example_dir()).expect("load workflow examples");
        let loaded = examples
            .iter()
            .find(|loaded| loaded.example.id == "tp73_isoform_protein_gel_offline")
            .expect("isoform protein gel demo should exist");
        let run_dir = TempDir::new().expect("temp run dir");
        run_example_workflow_in_dir(&loaded.example, Path::new("."), run_dir.path())
            .expect("tp73 protein-gel workflow should execute");
        let svg_path = run_dir.path().join("exports/tp73_isoform_protein_gel.svg");
        assert!(svg_path.exists());
        let svg = fs::read_to_string(&svg_path).expect("read protein gel demo svg");
        assert!(svg.contains("Protein Gel Preview"));
        assert!(svg.contains("Protein Ladder 10-100 kDa"));
        assert!(svg.contains("Selection notes"));
    }

    #[test]
    fn workflow_examples_isoform_protein_2d_gel_demo_writes_svg() {
        let examples = load_workflow_examples(&example_dir()).expect("load workflow examples");
        let loaded = examples
            .iter()
            .find(|loaded| loaded.example.id == "tp73_isoform_protein_2d_gel_offline")
            .expect("isoform protein 2D gel demo should exist");
        let run_dir = TempDir::new().expect("temp run dir");
        run_example_workflow_in_dir(&loaded.example, Path::new("."), run_dir.path())
            .expect("tp73 protein 2D gel workflow should execute");
        let svg_path = run_dir
            .path()
            .join("exports/tp73_isoform_protein_2d_gel.svg");
        assert!(svg_path.exists());
        let svg = fs::read_to_string(&svg_path).expect("read tp73 protein 2D gel svg");
        assert!(svg.contains("Protein 2D Gel Preview"));
        assert!(svg.contains("Protein Ladder 10-100 kDa"));
        assert!(svg.contains("Selection notes"));
        assert!(svg.contains("Spot details"));
    }

    #[test]
    fn workflow_examples_trypsin_digest_gel_demo_writes_svg() {
        let examples = load_workflow_examples(&example_dir()).expect("load workflow examples");
        let loaded = examples
            .iter()
            .find(|loaded| loaded.example.id == "tp73_variant1_trypsin_digest_gel_offline")
            .expect("tp73 trypsin digest gel example should exist");
        let run_dir = TempDir::new().expect("temp run dir");
        run_example_workflow_in_dir(&loaded.example, Path::new("."), run_dir.path())
            .expect("tp73 trypsin digest workflow should execute");
        let svg_path = run_dir
            .path()
            .join("exports/tp73_variant1_trypsin_digest_gel.svg");
        assert!(svg_path.exists());
        let svg = fs::read_to_string(&svg_path).expect("read tp73 trypsin digest gel svg");
        assert!(svg.contains("Protein Gel Preview"));
        assert!(svg.contains("Proteases: Trypsin"));
        assert!(svg.contains("Peptides shown"));
    }

    #[test]
    fn workflow_examples_simple_pcr_primer_design_exports_figure_and_report() {
        let examples = load_workflow_examples(&example_dir()).expect("load workflow examples");
        let loaded = examples
            .iter()
            .find(|loaded| loaded.example.id == "simple_pcr_primer_design_offline")
            .expect("simple PCR primer-design example should exist");
        let run_dir = TempDir::new().expect("temp run dir");
        run_example_workflow_in_dir(&loaded.example, Path::new("."), run_dir.path())
            .expect("simple PCR primer-design workflow should execute");
        let report_path = run_dir
            .path()
            .join("artifacts/simple_pcr_demo_primers.report.json");
        let figure_path = run_dir
            .path()
            .join("artifacts/simple_pcr_demo_primers.protocol.svg");
        assert!(report_path.exists());
        assert!(figure_path.exists());
        let figure = fs::read_to_string(&figure_path).expect("read simple PCR figure");
        assert!(figure.contains("data-protocol-id=\"pcr.assay.pair\""));
        assert!(figure.contains("Primer Constraints"));
        assert!(figure.contains("Product"));
        let report: serde_json::Value = serde_json::from_str(
            &fs::read_to_string(&report_path).expect("read simple PCR report"),
        )
        .expect("parse simple PCR report");
        assert_eq!(
            report["schema"].as_str(),
            Some("gentle.primer_design_report.v1")
        );
        assert_eq!(report["template"].as_str(), Some("simple_pcr_template"));
        assert_eq!(report["pair_count"].as_u64(), Some(5));
    }

    #[test]
    fn workflow_examples_tp73_evidence_viewer_release_proof_writes_artifacts_and_features() {
        let _serial = lock_jaspar_registry_for_test();
        crate::tf_motifs::reload_builtin_for_test();
        let examples = load_workflow_examples(&example_dir()).expect("load workflow examples");
        let loaded = examples
            .iter()
            .find(|loaded| loaded.example.id == "tp73_genome_evidence_viewer_release_proof")
            .expect("TP73 evidence-viewer proof example should exist");
        let run_dir = TempDir::new().expect("temp run dir");
        let state =
            run_example_workflow_for_project_state(&loaded.example, Path::new("."), run_dir.path())
                .expect("TP73 evidence-viewer proof workflow should execute");

        for relative in [
            "artifacts/tp73_evidence_viewer/tp73_evidence_viewer.linear.svg",
            "artifacts/tp73_evidence_viewer/tp73_evidence_viewer.splicing.expert.svg",
            "artifacts/tp73_evidence_viewer/tp73_evidence_viewer.tfbs_score_tracks.svg",
            "artifacts/tp73_evidence_viewer/tp73_evidence_viewer.repeat_materialization.json",
            "artifacts/tp73_evidence_viewer/tp73_evidence_viewer.tfbs_score_tracks.json",
        ] {
            let path = run_dir.path().join(relative);
            assert!(path.exists(), "expected workflow artifact {relative}");
            assert!(
                path.metadata().expect("artifact metadata").len() > 0,
                "expected non-empty workflow artifact {relative}"
            );
        }

        let repeat_report: serde_json::Value = serde_json::from_str(
            &fs::read_to_string(run_dir.path().join(
                "artifacts/tp73_evidence_viewer/tp73_evidence_viewer.repeat_materialization.json",
            ))
            .expect("read repeat materialization report"),
        )
        .expect("parse repeat materialization report");
        assert_eq!(
            repeat_report["schema"].as_str(),
            Some("gentle.repeat_feature_materialization.v1")
        );
        assert_eq!(repeat_report["added_feature_count"].as_u64(), Some(3));

        let engine = GentleEngine::from_state(state.clone());
        let anchor = engine
            .sequence_genome_anchor_summary("tp73_evidence_viewer")
            .expect("TP73 proof sequence should keep its genome anchor");
        assert_eq!(anchor.genome_id, "GRCh38.p14");
        assert_eq!(anchor.chromosome, "1");
        assert_eq!(anchor.start_1based, 3652516);
        assert_eq!(anchor.end_1based, 3736201);

        let dna = state
            .sequences
            .get("tp73_evidence_viewer")
            .expect("TP73 proof sequence should be present");
        let repeat_count = dna
            .features()
            .iter()
            .filter(|feature| {
                first_qualifier(feature, "gentle_generated").as_deref() == Some("ucsc_rmsk")
            })
            .count();
        assert_eq!(repeat_count, 3);

        assert!(dna.features().iter().any(|feature| {
            first_qualifier(feature, "gentle_track_source").as_deref() == Some("Array")
                && first_qualifier(feature, "gentle_array_dataset").as_deref()
                    == Some("E-MTAB-14704")
                && first_qualifier(feature, "feature_id").as_deref() == Some("PSR_TP73_0001")
                && first_qualifier(feature, "logFC").is_some()
        }));
        assert!(dna.features().iter().any(|feature| {
            first_qualifier(feature, "gentle_track_source").as_deref() == Some("BED")
                && first_qualifier(feature, "gentle_track_name").as_deref()
                    == Some("TP73 CUT&RUN proof BED")
                && first_qualifier(feature, "score").as_deref() == Some("650.000000")
        }));
        assert!(dna.features().iter().any(|feature| {
            feature.kind.to_string().eq_ignore_ascii_case("TFBS")
                || first_qualifier(feature, "gentle_generated").as_deref() == Some("tfbs_scan")
        }));
        assert!(state.display.show_repeat_features);
        assert!(state.display.show_array_features);
        assert!(state.display.show_tfbs);
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
            assert!(
                entry.review_status.as_deref().is_some(),
                "tutorial catalog entry '{}' should surface review status",
                entry.id
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
    fn tutorial_sources_are_v4_with_catalog_group_assignments() {
        let sources =
            load_tutorial_source_units(&tutorial_source_dir()).expect("load tutorial sources");
        assert_eq!(sources.len(), 40);
        let unnumbered_reference_units = HashSet::from([
            "generated_hub".to_string(),
            "tutorial_landscape_overview".to_string(),
        ]);
        for source in &sources {
            assert_eq!(
                source.schema, TUTORIAL_SOURCE_SCHEMA,
                "source '{}' should be migrated to the current tutorial source schema",
                source.id
            );
            let catalog = source
                .catalog
                .as_ref()
                .unwrap_or_else(|| panic!("source '{}' should define a catalog entry", source.id));
            assert!(
                catalog.group.is_some(),
                "source '{}' should define catalog.group",
                source.id
            );
            if unnumbered_reference_units.contains(&source.id) {
                assert!(
                    catalog.group_position.is_none(),
                    "reference source '{}' should be grouped but unnumbered",
                    source.id
                );
            } else {
                assert!(
                    catalog.group_position.is_some(),
                    "tutorial source '{}' should define catalog.group_position",
                    source.id
                );
            }
        }
    }

    #[test]
    fn tutorial_catalog_meta_defines_decimal_group_taxonomy() {
        let meta =
            load_tutorial_catalog_meta(&tutorial_catalog_meta_path()).expect("load catalog meta");
        assert_eq!(meta.schema, TUTORIAL_CATALOG_META_SCHEMA);
        assert_eq!(meta.groups.len(), 10);
        let group_codes = meta
            .groups
            .iter()
            .map(|group| group.code.as_str())
            .collect::<HashSet<_>>();
        assert!(group_codes.contains("01"));
        assert!(group_codes.contains("10"));
    }

    #[test]
    fn tutorial_catalog_entries_have_group_and_decimal_placement() {
        let catalog =
            load_tutorial_catalog(&tutorial_catalog_path()).expect("load tutorial catalog");
        assert_eq!(catalog.entries.len(), 40);
        let unnumbered_reference_units = HashSet::from([
            "generated_hub".to_string(),
            "tutorial_landscape_overview".to_string(),
        ]);
        for entry in &catalog.entries {
            assert!(
                entry.group.is_some(),
                "catalog entry '{}' should include a tutorial group",
                entry.id
            );
            if unnumbered_reference_units.contains(&entry.id) {
                assert!(
                    entry.decimal_id.is_none(),
                    "reference catalog entry '{}' should not have a decimal id",
                    entry.id
                );
            } else {
                assert!(
                    entry.group_position.is_some(),
                    "catalog entry '{}' should include a position within its group",
                    entry.id
                );
                assert!(
                    entry.decimal_id.is_some(),
                    "catalog entry '{}' should include a derived decimal id",
                    entry.id
                );
            }
        }
    }

    #[test]
    fn tutorial_dual_surface_units_are_one_source_unit() {
        let sources =
            load_tutorial_source_units(&tutorial_source_dir()).expect("load tutorial sources");
        for id in [
            "simple_pcr_selection_gui",
            "tp73_uniprot_projection_audit_cli",
        ] {
            let unit = sources
                .iter()
                .find(|unit| unit.id == id)
                .unwrap_or_else(|| panic!("missing source unit '{id}'"));
            assert!(
                unit.catalog.is_some(),
                "dual-surface source unit '{id}' should define a catalog entry"
            );
            assert!(
                unit.generated_chapter.is_some(),
                "dual-surface source unit '{id}' should define a generated chapter"
            );
        }
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
    fn tutorial_review_manifest_stale_date_warns() {
        let dir = TempDir::new().expect("temp tutorial dir");
        let manifest_path = dir.path().join("manifest.json");
        let review_path = dir.path().join("review_manifest.json");
        let manifest = TutorialManifest {
            schema: TUTORIAL_MANIFEST_SCHEMA.to_string(),
            description: "test manifest".to_string(),
            concepts: vec![TutorialConcept {
                id: "concept".to_string(),
                name: "Concept".to_string(),
                description: String::new(),
            }],
            chapters: vec![minimal_tutorial_chapter("stale_review")],
        };
        fs::write(
            &review_path,
            r#"{
  "schema": "gentle.tutorial_review_manifest.v1",
  "warn_after_months": 12,
  "entries": [
    {
      "tutorial_id": "stale_review",
      "tutorial_kind": "generated_chapter",
      "tutorial_status": "active",
      "replaced_by": null,
      "codex_reviewed_at": null,
      "human_reviewed_at": "2024-01-01",
      "human_reviewer": "smoe"
    }
  ]
}
"#,
        )
        .expect("write review manifest");

        let context =
            tutorial_review_context_for_date(&manifest_path, &manifest, Some((2026, 5, 21)))
                .expect("review context");

        assert!(
            context
                .warnings
                .iter()
                .any(|warning| warning.contains("older than 12 months")),
            "expected stale review warning, got {:?}",
            context.warnings
        );
    }

    #[test]
    fn tutorial_review_dependency_staleness_marks_source_changes() {
        let dir = TempDir::new().expect("temp tutorial dir");
        let repo_root = dir.path();
        let tutorial_dir = repo_root.join("docs/tutorial");
        let source_dir = tutorial_dir.join("sources");
        let manifest_path = tutorial_dir.join("manifest.json");
        let review_path = tutorial_dir.join("review_manifest.json");
        fs::create_dir_all(&source_dir).expect("create source dir");
        fs::create_dir_all(repo_root.join("docs/examples/workflows")).expect("create workflow dir");
        fs::write(
            source_dir.join("01_stale_dependency.json"),
            r#"{
  "schema": "gentle.tutorial_source.v4",
  "id": "stale_dependency",
  "title": "Stale Dependency",
  "group": "01",
  "group_position": 1,
  "catalog": {
    "order": 1,
    "path": "docs/tutorial/stale_dependency.md",
    "group": "01",
    "group_position": 1,
    "type": "guided_walkthrough",
    "status": "manual",
    "source": "synthetic",
    "audiences": ["test"],
    "notes": "Synthetic source for dependency staleness tests."
  }
}
"#,
        )
        .expect("write source");
        fs::write(
            tutorial_dir.join("stale_dependency.md"),
            "# Stale Dependency\n",
        )
        .expect("write markdown");
        fs::write(
            &review_path,
            r#"{
  "schema": "gentle.tutorial_review_manifest.v1",
  "warn_after_months": 240,
  "entries": [
    {
      "tutorial_id": "stale_dependency",
      "tutorial_kind": "guided_walkthrough",
      "tutorial_status": "active",
      "replaced_by": null,
      "codex_reviewed_at": null,
      "human_reviewed_at": "2024-01-01",
      "human_reviewer": "smoe"
    }
  ]
}
"#,
        )
        .expect("write review manifest");
        let manifest = TutorialManifest {
            schema: TUTORIAL_MANIFEST_SCHEMA.to_string(),
            description: "test manifest".to_string(),
            concepts: vec![TutorialConcept {
                id: "concept".to_string(),
                name: "Concept".to_string(),
                description: String::new(),
            }],
            chapters: vec![],
        };

        let context =
            tutorial_review_context_for_date(&manifest_path, &manifest, Some((2024, 1, 2)))
                .expect("review context");
        let reason = context
            .stale_reasons_by_id
            .get("stale_dependency")
            .expect("stale dependency reason");

        assert!(
            reason.contains("source JSON") && reason.contains("changed after human review date"),
            "expected source-json stale reason, got {reason}"
        );
        assert!(
            context
                .warnings
                .iter()
                .any(|warning| warning.contains("review is stale because source JSON")),
            "expected dependency stale warning, got {:?}",
            context.warnings
        );
    }

    #[test]
    fn tutorial_catalog_projects_stale_reason_and_feedback_template() {
        let dir = TempDir::new().expect("temp tutorial dir");
        let repo_root = dir.path();
        let tutorial_dir = repo_root.join("docs/tutorial");
        let source_dir = tutorial_dir.join("sources");
        let meta_path = source_dir.join("catalog_meta.json");
        let review_path = tutorial_dir.join("review_manifest.json");
        fs::create_dir_all(&source_dir).expect("create source dir");
        fs::write(
            &meta_path,
            r#"{
  "schema": "gentle.tutorial_catalog_meta.v2",
  "description": "Synthetic tutorial catalog metadata.",
  "entry_page": "docs/tutorial/README.md",
  "generated_runtime": {
    "manifest_path": "docs/tutorial/manifest.json",
    "manifest_schema": "gentle.tutorial_manifest.v2",
    "generated_readme": "docs/tutorial/generated/README.md",
    "generation_report": "docs/tutorial/generated/report.json"
  },
  "groups": [
    { "code": "01", "label": "Getting Started", "order": 1 }
  ],
  "concepts": []
}
"#,
        )
        .expect("write catalog meta");
        fs::write(
            source_dir.join("01_stale_catalog.json"),
            r#"{
  "schema": "gentle.tutorial_source.v4",
  "id": "stale_catalog",
  "title": "Stale Catalog",
  "group": "01",
  "group_position": 1,
  "catalog": {
    "order": 1,
    "path": "docs/tutorial/stale_catalog.md",
    "group": "01",
    "group_position": 1,
    "type": "guided_walkthrough",
    "status": "manual",
    "source": "synthetic",
    "audiences": ["test"],
    "notes": "Synthetic source for stale catalog projection tests."
  }
}
"#,
        )
        .expect("write source");
        fs::write(tutorial_dir.join("stale_catalog.md"), "# Stale Catalog\n")
            .expect("write markdown");
        fs::write(
            &review_path,
            r#"{
  "schema": "gentle.tutorial_review_manifest.v1",
  "warn_after_months": 240,
  "entries": [
    {
      "tutorial_id": "stale_catalog",
      "tutorial_kind": "guided_walkthrough",
      "tutorial_status": "active",
      "replaced_by": null,
      "codex_reviewed_at": null,
      "human_reviewed_at": "2024-01-01",
      "human_reviewer": "smoe"
    }
  ]
}
"#,
        )
        .expect("write review manifest");

        let catalog = generate_tutorial_catalog_from_sources(&meta_path, &source_dir)
            .expect("generate tutorial catalog");
        let entry = catalog
            .entries
            .iter()
            .find(|entry| entry.id == "stale_catalog")
            .expect("catalog entry");

        assert!(entry.review_stale, "expected stale review projection");
        assert!(
            entry
                .review_stale_reason
                .as_deref()
                .is_some_and(|reason| reason.contains("source JSON")),
            "expected source JSON stale reason, got {:?}",
            entry.review_stale_reason
        );
        assert_eq!(
            entry.review_issue_template_path.as_deref(),
            Some(TUTORIAL_CONFUSION_ISSUE_TEMPLATE_PATH)
        );
    }

    #[test]
    fn tutorial_catalog_projects_graphic_metadata_and_missing_graphic_staleness() {
        let dir = TempDir::new().expect("temp tutorial dir");
        let repo_root = dir.path();
        let tutorial_dir = repo_root.join("docs/tutorial");
        let source_dir = tutorial_dir.join("sources");
        let meta_path = source_dir.join("catalog_meta.json");
        let review_path = tutorial_dir.join("review_manifest.json");
        fs::create_dir_all(&source_dir).expect("create source dir");
        fs::write(
            &meta_path,
            r#"{
  "schema": "gentle.tutorial_catalog_meta.v2",
  "description": "Synthetic tutorial catalog metadata.",
  "entry_page": "docs/tutorial/README.md",
  "generated_runtime": {
    "manifest_path": "docs/tutorial/manifest.json",
    "manifest_schema": "gentle.tutorial_manifest.v2",
    "generated_readme": "docs/tutorial/generated/README.md",
    "generation_report": "docs/tutorial/generated/report.json"
  },
  "groups": [
    { "code": "01", "label": "Getting Started", "order": 1 }
  ],
  "concepts": []
}
"#,
        )
        .expect("write catalog meta");
        fs::write(
            source_dir.join("01_missing_graphic.json"),
            r#"{
  "schema": "gentle.tutorial_source.v4",
  "id": "missing_graphic",
  "title": "Missing Graphic",
  "group": "01",
  "group_position": 1,
  "graphics": [
    {
      "kind": "screenshot",
      "path": "docs/screenshots/missing.png",
      "caption": "Missing tutorial screenshot.",
      "illustrates_step": 1,
      "capture_date": "2026-03-20"
    }
  ],
  "catalog": {
    "order": 1,
    "path": "docs/tutorial/missing_graphic.md",
    "group": "01",
    "group_position": 1,
    "type": "guided_walkthrough",
    "status": "manual",
    "source": "synthetic",
    "audiences": ["test"],
    "notes": "Synthetic source for graphic metadata tests."
  }
}
"#,
        )
        .expect("write source");
        fs::write(
            tutorial_dir.join("missing_graphic.md"),
            "# Missing Graphic\n",
        )
        .expect("write markdown");
        fs::write(
            &review_path,
            r#"{
  "schema": "gentle.tutorial_review_manifest.v1",
  "warn_after_months": 240,
  "entries": [
    {
      "tutorial_id": "missing_graphic",
      "tutorial_kind": "guided_walkthrough",
      "tutorial_status": "active",
      "replaced_by": null,
      "codex_reviewed_at": null,
      "human_reviewed_at": "2026-03-21",
      "human_reviewer": "smoe"
    }
  ]
}
"#,
        )
        .expect("write review manifest");

        let catalog = generate_tutorial_catalog_from_sources(&meta_path, &source_dir)
            .expect("generate tutorial catalog");
        let entry = catalog
            .entries
            .iter()
            .find(|entry| entry.id == "missing_graphic")
            .expect("catalog entry");

        assert_eq!(entry.graphics.len(), 1);
        assert!(
            entry.review_stale,
            "missing graphic should stale the review"
        );
        assert!(
            entry
                .review_stale_reason
                .as_deref()
                .is_some_and(|reason| reason.contains("declared graphic")),
            "expected declared-graphic stale reason, got {:?}",
            entry.review_stale_reason
        );
        assert_eq!(
            entry.review_issue_template_path.as_deref(),
            Some(TUTORIAL_ARTIFACT_ISSUE_TEMPLATE_PATH)
        );
    }

    #[test]
    fn tutorial_review_manifest_missing_entry_warns() {
        let dir = TempDir::new().expect("temp tutorial dir");
        let manifest_path = dir.path().join("manifest.json");
        let review_path = dir.path().join("review_manifest.json");
        let manifest = TutorialManifest {
            schema: TUTORIAL_MANIFEST_SCHEMA.to_string(),
            description: "test manifest".to_string(),
            concepts: vec![TutorialConcept {
                id: "concept".to_string(),
                name: "Concept".to_string(),
                description: String::new(),
            }],
            chapters: vec![minimal_tutorial_chapter("missing_review")],
        };
        fs::write(
            &review_path,
            r#"{
  "schema": "gentle.tutorial_review_manifest.v1",
  "warn_after_months": 12,
  "entries": []
}
"#,
        )
        .expect("write review manifest");

        let context =
            tutorial_review_context_for_date(&manifest_path, &manifest, Some((2026, 5, 21)))
                .expect("review context");

        assert!(
            context
                .warnings
                .iter()
                .any(|warning| warning.contains("missing entry for tutorial id 'missing_review'")),
            "expected missing-entry warning, got {:?}",
            context.warnings
        );
    }

    #[test]
    fn deprecated_tutorial_execution_failure_is_nonfatal() {
        let dir = TempDir::new().expect("temp tutorial dir");
        let examples_dir = dir.path().join("examples");
        let tutorial_dir = dir.path().join("tutorial");
        let output_dir = tutorial_dir.join("generated");
        fs::create_dir_all(&examples_dir).expect("create examples dir");
        fs::create_dir_all(&tutorial_dir).expect("create tutorial dir");
        fs::write(
            examples_dir.join("deprecated_example.json"),
            r#"{
  "schema": "gentle.workflow_example.v1",
  "id": "deprecated_example",
  "title": "Deprecated Example",
  "summary": "Synthetic deprecated tutorial example.",
  "test_mode": "always",
  "required_files": [],
  "tags": [],
  "workflow": {
    "run_id": "deprecated_example",
    "ops": [
      {
        "LoadFile": {
          "path": "missing-input.fa",
          "as_id": "missing_input"
        }
      }
    ]
  }
}
"#,
        )
        .expect("write workflow example");
        fs::write(
            tutorial_dir.join("manifest.json"),
            r#"{
  "schema": "gentle.tutorial_manifest.v1",
  "description": "Synthetic tutorial manifest.",
  "concepts": [
    {
      "id": "concept",
      "name": "Concept",
      "description": "Synthetic concept."
    }
  ],
  "chapters": [
    {
      "id": "deprecated_tutorial",
      "order": 1,
      "title": "Deprecated Tutorial",
      "example_id": "deprecated_example",
      "tier": "core",
      "guide_path": null,
      "summary": "Synthetic deprecated tutorial.",
      "narrative": "Synthetic narrative.",
      "use_cases": ["Test deprecated review handling."],
      "gui_steps": ["Open the synthetic tutorial."],
      "cli_steps": [],
      "step_expectations": [],
      "prerequisites": [],
      "local_execution_note": null,
      "learning_objectives": ["Understand deprecated tutorial handling."],
      "concepts": ["concept"],
      "parameter_notes": [],
      "follow_up_commands": [],
      "retain_outputs": [],
      "checkpoints": []
    }
  ]
}
"#,
        )
        .expect("write tutorial manifest");
        fs::write(
            tutorial_dir.join("review_manifest.json"),
            r#"{
  "schema": "gentle.tutorial_review_manifest.v1",
  "warn_after_months": 12,
  "entries": [
    {
      "tutorial_id": "deprecated_tutorial",
      "tutorial_kind": "generated_chapter",
      "tutorial_status": "deprecated",
      "replaced_by": null,
      "codex_reviewed_at": null,
      "human_reviewed_at": null,
      "human_reviewer": null
    }
  ]
}
"#,
        )
        .expect("write review manifest");

        let report = generate_tutorial_docs(
            &examples_dir,
            &tutorial_dir.join("manifest.json"),
            &output_dir,
            Path::new("."),
        )
        .expect("deprecated tutorial execution failure should warn, not fail");

        assert!(
            report
                .warnings
                .iter()
                .any(|warning| warning
                    .contains("did not block generation despite execution failure")),
            "expected deprecated execution warning, got {:?}",
            report.warnings
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
        let _serial = lock_jaspar_registry_for_test();
        crate::tf_motifs::reload_builtin_for_test();
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
        let _serial = lock_jaspar_registry_for_test();
        crate::tf_motifs::reload_builtin_for_test();
        let source = example_dir();
        let manifest = tutorial_manifest_path();
        let repo_root = Path::new(".");
        let out_dir = TempDir::new().expect("create temp dir");
        let generated = out_dir.path().join("generated");
        generate_tutorial_docs(&source, &manifest, &generated, repo_root)
            .expect("tutorial generation should pass");
        let readme = generated.join("README.md");
        let readme_markdown = std::fs::read_to_string(&readme).expect("read generated README");
        assert!(readme_markdown.contains("reference material - not first-time walkthroughs"));
        assert!(readme_markdown.contains("guided walkthroughs in [`../README.md`](../README.md)"));
        assert!(readme_markdown.contains("### Sequence Basics & Lineage"));
        assert!(readme_markdown.contains("`02.01` [Load FASTA, branch, and reverse-complement]"));
        assert!(readme_markdown.contains("review `unreviewed`"));

        let chapter = generated.join("chapters/02-01_load_branch_reverse_complement_pgex_fasta.md");
        let markdown = std::fs::read_to_string(&chapter).expect("read generated chapter markdown");
        assert!(markdown.starts_with("---\nchapter_id: "));
        assert!(markdown.contains("## What You Learn"));
        assert!(markdown.contains("## Applied Concepts"));
        assert!(!markdown.contains("Reoccurs in:"));
        assert!(markdown.contains("## GUI First"));
        assert!(markdown.contains("## Parameters That Matter"));
        assert!(
            markdown.find("## Parameters That Matter").unwrap()
                < markdown.find("## GUI First").unwrap()
        );
        assert!(markdown.contains("## Tutorial Provenance"));
        assert!(markdown.contains("review_status: "));
        assert!(markdown.contains("review_stale: "));
        assert!(markdown.contains("review_issue_template_path: "));
        assert!(markdown.contains("## Feedback"));

        let online_chapter = generated.join("chapters/05-02_prepare_reference_genome_online.md");
        let online_markdown =
            std::fs::read_to_string(&online_chapter).expect("read online chapter markdown");
        assert!(online_markdown.contains("> **How to Run This Locally**"));
        assert!(online_markdown.contains("### Step 1:"));
        assert!(online_markdown.contains("CLI:\n\n```bash\nGENTLE_TEST_ONLINE=1 cargo run"));
        assert!(online_markdown.contains("> Expected:"));
        assert!(!online_markdown.contains("## Command Equivalent (After GUI)"));

        let simple_pcr_chapter = generated.join("chapters/04-01_simple_pcr_selection_gui.md");
        let simple_pcr_markdown =
            std::fs::read_to_string(&simple_pcr_chapter).expect("read simple PCR chapter markdown");
        assert!(simple_pcr_markdown.contains("See also: guided walkthrough"));
        assert!(simple_pcr_markdown.contains("../../04-01_simple_pcr_selection_gui.md"));

        let promoter_chapter =
            generated.join("chapters/08-03_promoter_design_artifact_slice_offline.md");
        let promoter_markdown =
            std::fs::read_to_string(&promoter_chapter).expect("read promoter chapter markdown");
        assert!(promoter_markdown.contains("**Prerequisites:** Read [Chapter 1:"));
        assert!(promoter_markdown.contains("### Step 8:"));
        assert!(promoter_markdown.contains("CLI:\n\n```bash\ncargo run --bin gentle_cli"));
        assert!(promoter_markdown
            .contains("CLI snippets use GENtle's default `.gentle_state.json` state"));
        assert!(promoter_markdown.contains("> Expected: `tfbs_score_tracks.svg`"));
        assert!(!promoter_markdown.contains("## Command Equivalent (After GUI)"));
        assert!(promoter_markdown.contains("## What This Chapter Produces"));
        assert!(promoter_markdown.contains("> SVG text labels:"));
        assert!(promoter_markdown.contains(
            "![TFBS score tracks across the synthetic TP73 promoter slice.](../artifacts/promoter_design_artifact_slice_offline/artifacts/tp73_promoter_artifact_demo.tfbs_score_tracks.svg)"
        ));
        assert!(promoter_markdown.contains("Embedded above near Step 8"));
        let step_8 = promoter_markdown.find("### Step 8:").unwrap();
        let figure = promoter_markdown
            .find("![TFBS score tracks across")
            .unwrap();
        let step_9 = promoter_markdown.find("### Step 9:").unwrap();
        assert!(
            step_8 < figure && figure < step_9,
            "declared graphic should render inline at the step it illustrates"
        );
        assert!(promoter_markdown.contains("schema: `gentle.promoter_artifact_manifest.v1`"));
    }

    #[test]
    fn tutorial_cdna_genomic_demo_markdown_image_links_exist() {
        let tutorial_path = PathBuf::from("docs/tutorial/02-03_tp73_cdna_genomic_dotplot_gui.md");
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
                        expression_tsv_path: None,
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
    fn rewrite_example_paths_handles_promoter_and_handoff_outputs() {
        let example = WorkflowExample {
            schema: WORKFLOW_EXAMPLE_SCHEMA.to_string(),
            id: "promoter_output_rewrite_test".to_string(),
            title: "promoter output rewrite test".to_string(),
            summary: String::new(),
            test_mode: ExampleTestMode::Skip,
            required_files: vec![],
            tags: vec![],
            workflow: Workflow {
                run_id: "promoter_output_rewrite_test".to_string(),
                ops: vec![
                    Operation::SummarizeAlternativePromoterComparison {
                        input: "seq_a".to_string(),
                        gene_label: Some("TP73".to_string()),
                        transcript_id: None,
                        promoter_upstream_bp: 40,
                        promoter_downstream_bp: 15,
                        path: Some("artifacts/alt_promoters.json".to_string()),
                    },
                    Operation::SummarizePromoterEvidenceMatrix {
                        input: "seq_a".to_string(),
                        gene_label: Some("TP73".to_string()),
                        transcript_id: None,
                        promoter_upstream_bp: 40,
                        promoter_downstream_bp: 15,
                        include_feature_overlaps: true,
                        path: Some("artifacts/evidence_matrix.json".to_string()),
                    },
                    Operation::SummarizeIsoformPromoterComparison {
                        input: "seq_a".to_string(),
                        gene_label: Some("TP73".to_string()),
                        transcript_id: None,
                        promoter_upstream_bp: 40,
                        promoter_downstream_bp: 15,
                        include_feature_overlaps: true,
                        path: Some("artifacts/isoform_promoters.json".to_string()),
                    },
                    Operation::SummarizePromoterExpressionEvidence {
                        input: "seq_a".to_string(),
                        gene_label: Some("TP73".to_string()),
                        transcript_id: None,
                        promoter_upstream_bp: 40,
                        promoter_downstream_bp: 15,
                        expression_rows: vec![],
                        expression_source_label: Some("demo_expression".to_string()),
                        path: Some("artifacts/expression_evidence.json".to_string()),
                    },
                    Operation::ExportPromoterArtifactManifest {
                        input: "seq_a".to_string(),
                        gene_label: Some("TP73".to_string()),
                        artifacts: vec![crate::engine::PromoterArtifactManifestEntry {
                            artifact_id: "evidence_matrix".to_string(),
                            artifact_kind: "promoter_evidence_matrix".to_string(),
                            path: "artifacts/evidence_matrix.json".to_string(),
                            schema_hint: Some("gentle.promoter_evidence_matrix.v1".to_string()),
                            label: None,
                            recommended_use: None,
                            required: false,
                        }],
                        path: "artifacts/promoter_manifest.json".to_string(),
                    },
                    Operation::ExportLabAssistantInstructions {
                        path: "artifacts/handoff.md".to_string(),
                        run_id: None,
                        title: None,
                        audience: None,
                        format: None,
                    },
                ],
            },
        };
        let repo_root = std::env::current_dir().expect("cwd");
        let run_dir = TempDir::new().expect("temp run dir");
        let rewritten =
            rewrite_example_paths_for_execution(&example, repo_root.as_path(), run_dir.path())
                .expect("rewrite should succeed");
        for op in &rewritten.workflow.ops {
            let path = match op {
                Operation::SummarizeAlternativePromoterComparison { path, .. }
                | Operation::SummarizePromoterEvidenceMatrix { path, .. }
                | Operation::SummarizeIsoformPromoterComparison { path, .. }
                | Operation::SummarizePromoterExpressionEvidence { path, .. } => path
                    .as_deref()
                    .expect("optional path should remain present"),
                Operation::ExportPromoterArtifactManifest { path, .. } => path.as_str(),
                Operation::ExportLabAssistantInstructions { path, .. } => path.as_str(),
                other => panic!("unexpected operation: {other:?}"),
            };
            assert!(
                path.starts_with(&display_path(run_dir.path())),
                "output path should be rewritten into run dir: {path}"
            );
            if let Operation::ExportPromoterArtifactManifest { artifacts, .. } = op {
                assert!(
                    artifacts
                        .iter()
                        .all(|artifact| !Path::new(&artifact.path).is_absolute()
                            && !artifact.path.starts_with(&display_path(run_dir.path()))),
                    "manifest artifact paths should stay portable relative paths: {artifacts:?}"
                );
            }
        }
    }

    #[test]
    fn retained_tutorial_artifact_normalization_removes_wall_clock_fields() {
        let raw = "- Generated (Unix ms): `1777756509715`\n{\n  \"generated_at_unix_ms\": 1777756511961,\n}\n\n";
        let normalized = normalize_retained_tutorial_artifact_text(raw);
        assert!(normalized.contains("- Generated (Unix ms): `0`"));
        assert!(normalized.contains("\"generated_at_unix_ms\": 0,"));
        assert!(!normalized.contains("1777756509715"));
        assert!(!normalized.contains("1777756511961"));
        assert!(!normalized.ends_with("\n\n"));
    }

    #[test]
    fn tutorial_file_count_mismatch_message_lists_path_delta() {
        let mut expected = BTreeMap::new();
        expected.insert("README.md".to_string(), Vec::new());
        let mut actual = expected.clone();
        actual.insert(
            "artifacts/new_demo/artifacts/new_retained_output.md".to_string(),
            Vec::new(),
        );

        let message = tutorial_generated_file_count_mismatch_message(&expected, &actual);

        assert!(message.contains("expected 1, got 2"));
        assert!(message
            .contains("generated-only files: artifacts/new_demo/artifacts/new_retained_output.md"));
        assert!(message.contains("committed-only files: none"));
    }

    #[test]
    fn tutorial_feedback_context_identifies_generated_chapter_paths() {
        let manifest =
            load_tutorial_manifest(&tutorial_manifest_path()).expect("load tutorial manifest");
        let examples = load_workflow_examples(&example_dir()).expect("load workflow examples");
        let chapter = manifest.chapters.first().expect("at least one chapter");
        let rel_path = format!("chapters/{}", markdown_path_for_chapter(chapter));

        let context = tutorial_check_feedback_context(
            &manifest,
            &examples,
            &tutorial_source_dir(),
            &tutorial_manifest_path(),
            &tutorial_output_dir(),
            Some(&rel_path),
            None,
            "generated tutorial file differs from committed version",
        );

        assert!(context.contains("Tutorial feedback context"));
        assert!(context.contains(&format!("Tutorial id: {}", chapter.id)));
        assert!(context.contains("Source JSON: docs/tutorial/sources/"));
        assert!(context.contains("Workflow: docs/examples/workflows/"));
        assert!(context.contains("Artifact dir: docs/tutorial/generated/artifacts/"));
        assert!(context.contains("Suggested issue template: Tutorial confusion"));
    }

    #[test]
    fn tutorial_feedback_context_suggests_artifact_issue_for_svg_artifacts() {
        let manifest =
            load_tutorial_manifest(&tutorial_manifest_path()).expect("load tutorial manifest");
        let examples = load_workflow_examples(&example_dir()).expect("load workflow examples");
        let chapter = manifest.chapters.first().expect("at least one chapter");
        let rel_path = format!(
            "artifacts/{}/artifacts/example.svg",
            markdown_file_stem(&chapter.id)
        );

        let context = tutorial_check_feedback_context(
            &manifest,
            &examples,
            &tutorial_source_dir(),
            &tutorial_manifest_path(),
            &tutorial_output_dir(),
            Some(&rel_path),
            None,
            "generated tutorial artifact differs from committed version",
        );

        assert!(context.contains(&format!("Tutorial id: {}", chapter.id)));
        assert!(context.contains("Affected generated path: artifacts/"));
        assert!(context.contains("Suggested issue template: Tutorial artifact/figure problem"));
    }

    #[test]
    fn tutorial_check_passes_on_committed_tree() {
        let _serial = lock_jaspar_registry_for_test();
        crate::tf_motifs::reload_builtin_for_test();
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
        let _serial = lock_jaspar_registry_for_test();
        crate::tf_motifs::reload_builtin_for_test();
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
        assert!(state
            .sequences
            .contains_key("gibson_destination_pgex_with_gibson_insert_demo"));

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
