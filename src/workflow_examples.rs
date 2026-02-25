//! Curated workflow example payloads and templates.

use crate::engine::{Engine, GentleEngine, ProjectState, Workflow};
use serde::{Deserialize, Serialize};
use std::{
    collections::HashSet,
    env, fs,
    path::{Path, PathBuf},
};

pub const WORKFLOW_EXAMPLE_SCHEMA: &str = "gentle.workflow_example.v1";
pub const DEFAULT_WORKFLOW_EXAMPLE_DIR: &str = "docs/examples/workflows";
pub const DEFAULT_WORKFLOW_SNIPPET_DIR: &str = "docs/examples/generated";
pub const ONLINE_EXAMPLE_TEST_ENV: &str = "GENTLE_TEST_ONLINE";

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
}
