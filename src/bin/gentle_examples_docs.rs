//! Documentation helper binary that renders workflow examples into generated snippets.

use gentle::workflow_examples::{
    DEFAULT_TUTORIAL_CATALOG_META_PATH, DEFAULT_TUTORIAL_CATALOG_PATH,
    DEFAULT_TUTORIAL_MANIFEST_PATH, DEFAULT_TUTORIAL_OUTPUT_DIR, DEFAULT_TUTORIAL_SOURCE_DIR,
    DEFAULT_WORKFLOW_EXAMPLE_DIR, DEFAULT_WORKFLOW_SNIPPET_DIR, check_tutorial_catalog_generated,
    check_tutorial_generated, generate_tutorial_docs, generate_workflow_example_docs,
    load_workflow_examples, validate_example_required_files, write_tutorial_catalog_from_sources,
};
use serde_json::json;
use std::{
    env,
    path::{Path, PathBuf},
    process,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Mode {
    ExampleGenerate,
    ExampleCheck,
    TutorialGenerate,
    TutorialCheck,
    TutorialCatalogGenerate,
    TutorialCatalogCheck,
}

#[derive(Debug)]
struct CliArgs {
    show_help: bool,
    mode: Mode,
    source_dir: String,
    example_output_dir: String,
    tutorial_catalog: String,
    tutorial_catalog_meta: String,
    tutorial_source_dir: String,
    tutorial_manifest: String,
    tutorial_output_dir: String,
    repo_root: String,
}

impl Default for CliArgs {
    fn default() -> Self {
        Self {
            show_help: false,
            mode: Mode::ExampleGenerate,
            source_dir: DEFAULT_WORKFLOW_EXAMPLE_DIR.to_string(),
            example_output_dir: DEFAULT_WORKFLOW_SNIPPET_DIR.to_string(),
            tutorial_catalog: DEFAULT_TUTORIAL_CATALOG_PATH.to_string(),
            tutorial_catalog_meta: DEFAULT_TUTORIAL_CATALOG_META_PATH.to_string(),
            tutorial_source_dir: DEFAULT_TUTORIAL_SOURCE_DIR.to_string(),
            tutorial_manifest: DEFAULT_TUTORIAL_MANIFEST_PATH.to_string(),
            tutorial_output_dir: DEFAULT_TUTORIAL_OUTPUT_DIR.to_string(),
            repo_root: ".".to_string(),
        }
    }
}

fn usage() {
    eprintln!(
        "Usage:\n  \
gentle_examples_docs [generate] [--source DIR] [--output DIR]\n  \
gentle_examples_docs --check [--source DIR]\n  \
gentle_examples_docs tutorial-generate [--source DIR] [--manifest FILE] [--tutorial-output DIR] [--repo-root DIR]\n  \
gentle_examples_docs tutorial-check [--source DIR] [--manifest FILE] [--tutorial-output DIR] [--repo-root DIR]\n  \
gentle_examples_docs tutorial-catalog-generate [--catalog-meta FILE] [--tutorial-sources DIR] [--tutorial-catalog FILE]\n  \
gentle_examples_docs tutorial-catalog-check [--catalog-meta FILE] [--tutorial-sources DIR] [--tutorial-catalog FILE]\n\n  \
Defaults:\n  \
  --source {}\n  \
  --output {}\n  \
  --tutorial-catalog {}\n  \
  --catalog-meta {}\n  \
  --tutorial-sources {}\n  \
  --manifest {}\n  \
  --tutorial-output {}\n  \
  --repo-root .",
        DEFAULT_WORKFLOW_EXAMPLE_DIR,
        DEFAULT_WORKFLOW_SNIPPET_DIR,
        DEFAULT_TUTORIAL_CATALOG_PATH,
        DEFAULT_TUTORIAL_CATALOG_META_PATH,
        DEFAULT_TUTORIAL_SOURCE_DIR,
        DEFAULT_TUTORIAL_MANIFEST_PATH,
        DEFAULT_TUTORIAL_OUTPUT_DIR
    );
}

fn parse_args(args: &[String]) -> Result<CliArgs, String> {
    let mut parsed = CliArgs::default();
    let mut idx = 0usize;
    while idx < args.len() {
        match args[idx].as_str() {
            "--help" | "-h" => {
                parsed.show_help = true;
                idx += 1;
            }
            "--check" => {
                parsed.mode = Mode::ExampleCheck;
                idx += 1;
            }
            "--source" => {
                if idx + 1 >= args.len() {
                    return Err("Missing DIR after --source".to_string());
                }
                parsed.source_dir = args[idx + 1].clone();
                idx += 2;
            }
            "--output" => {
                if idx + 1 >= args.len() {
                    return Err("Missing DIR after --output".to_string());
                }
                parsed.example_output_dir = args[idx + 1].clone();
                parsed.tutorial_output_dir = args[idx + 1].clone();
                idx += 2;
            }
            "--tutorial-output" => {
                if idx + 1 >= args.len() {
                    return Err("Missing DIR after --tutorial-output".to_string());
                }
                parsed.tutorial_output_dir = args[idx + 1].clone();
                idx += 2;
            }
            "--tutorial-catalog" => {
                if idx + 1 >= args.len() {
                    return Err("Missing FILE after --tutorial-catalog".to_string());
                }
                parsed.tutorial_catalog = args[idx + 1].clone();
                idx += 2;
            }
            "--catalog-meta" => {
                if idx + 1 >= args.len() {
                    return Err("Missing FILE after --catalog-meta".to_string());
                }
                parsed.tutorial_catalog_meta = args[idx + 1].clone();
                idx += 2;
            }
            "--tutorial-sources" => {
                if idx + 1 >= args.len() {
                    return Err("Missing DIR after --tutorial-sources".to_string());
                }
                parsed.tutorial_source_dir = args[idx + 1].clone();
                idx += 2;
            }
            "--manifest" => {
                if idx + 1 >= args.len() {
                    return Err("Missing FILE after --manifest".to_string());
                }
                parsed.tutorial_manifest = args[idx + 1].clone();
                idx += 2;
            }
            "--repo-root" => {
                if idx + 1 >= args.len() {
                    return Err("Missing DIR after --repo-root".to_string());
                }
                parsed.repo_root = args[idx + 1].clone();
                idx += 2;
            }
            "generate" => {
                parsed.mode = Mode::ExampleGenerate;
                idx += 1;
            }
            "tutorial-generate" => {
                parsed.mode = Mode::TutorialGenerate;
                idx += 1;
            }
            "tutorial-check" => {
                parsed.mode = Mode::TutorialCheck;
                idx += 1;
            }
            "tutorial-catalog-generate" => {
                parsed.mode = Mode::TutorialCatalogGenerate;
                idx += 1;
            }
            "tutorial-catalog-check" => {
                parsed.mode = Mode::TutorialCatalogCheck;
                idx += 1;
            }
            other => {
                return Err(format!("Unknown argument '{other}'"));
            }
        }
    }
    Ok(parsed)
}

fn run_check_mode(source_dir: &Path) -> Result<(), String> {
    let examples = load_workflow_examples(source_dir)?;
    let repo_root = Path::new(".");
    for loaded in &examples {
        validate_example_required_files(&loaded.example, repo_root)?;
    }
    let summary = json!({
        "status": "ok",
        "source_dir": source_dir.to_string_lossy(),
        "example_count": examples.len(),
    });
    let pretty = serde_json::to_string_pretty(&summary)
        .map_err(|e| format!("Could not serialize check report: {e}"))?;
    println!("{pretty}");
    Ok(())
}

fn run_generate_mode(source_dir: &Path, output_dir: &Path) -> Result<(), String> {
    let report = generate_workflow_example_docs(source_dir, output_dir)?;
    let pretty = serde_json::to_string_pretty(&report)
        .map_err(|e| format!("Could not serialize generation report: {e}"))?;
    println!("{pretty}");
    Ok(())
}

fn run_tutorial_generate_mode(
    source_dir: &Path,
    manifest_path: &Path,
    output_dir: &Path,
    repo_root: &Path,
) -> Result<(), String> {
    let report = generate_tutorial_docs(source_dir, manifest_path, output_dir, repo_root)?;
    let pretty = serde_json::to_string_pretty(&report)
        .map_err(|e| format!("Could not serialize tutorial generation report: {e}"))?;
    println!("{pretty}");
    Ok(())
}

fn run_tutorial_check_mode(
    source_dir: &Path,
    manifest_path: &Path,
    output_dir: &Path,
    repo_root: &Path,
) -> Result<(), String> {
    let report = check_tutorial_generated(source_dir, manifest_path, output_dir, repo_root)?;
    let summary = json!({
        "status": "ok",
        "mode": "tutorial-check",
        "chapter_count": report.chapter_count,
        "generated_files": report.generated_files,
    });
    let pretty = serde_json::to_string_pretty(&summary)
        .map_err(|e| format!("Could not serialize tutorial check report: {e}"))?;
    println!("{pretty}");
    Ok(())
}

fn run_tutorial_catalog_generate_mode(
    meta_path: &Path,
    source_dir: &Path,
    output_path: &Path,
) -> Result<(), String> {
    let catalog = write_tutorial_catalog_from_sources(meta_path, source_dir, output_path)?;
    let summary = json!({
        "status": "ok",
        "mode": "tutorial-catalog-generate",
        "entry_count": catalog.entries.len(),
        "output_path": output_path.to_string_lossy(),
    });
    let pretty = serde_json::to_string_pretty(&summary)
        .map_err(|e| format!("Could not serialize tutorial catalog generation report: {e}"))?;
    println!("{pretty}");
    Ok(())
}

fn run_tutorial_catalog_check_mode(
    meta_path: &Path,
    source_dir: &Path,
    output_path: &Path,
) -> Result<(), String> {
    let catalog = check_tutorial_catalog_generated(meta_path, source_dir, output_path)?;
    let summary = json!({
        "status": "ok",
        "mode": "tutorial-catalog-check",
        "entry_count": catalog.entries.len(),
        "output_path": output_path.to_string_lossy(),
    });
    let pretty = serde_json::to_string_pretty(&summary)
        .map_err(|e| format!("Could not serialize tutorial catalog check report: {e}"))?;
    println!("{pretty}");
    Ok(())
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let parsed = match parse_args(&args) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("{e}");
            usage();
            process::exit(1);
        }
    };
    if parsed.show_help {
        usage();
        return;
    }

    let source_dir = PathBuf::from(parsed.source_dir);
    let example_output = PathBuf::from(parsed.example_output_dir);
    let tutorial_catalog = PathBuf::from(parsed.tutorial_catalog);
    let tutorial_catalog_meta = PathBuf::from(parsed.tutorial_catalog_meta);
    let tutorial_source_dir = PathBuf::from(parsed.tutorial_source_dir);
    let tutorial_manifest = PathBuf::from(parsed.tutorial_manifest);
    let tutorial_output = PathBuf::from(parsed.tutorial_output_dir);
    let repo_root = PathBuf::from(parsed.repo_root);
    let result = match parsed.mode {
        Mode::ExampleGenerate => run_generate_mode(&source_dir, &example_output),
        Mode::ExampleCheck => run_check_mode(&source_dir),
        Mode::TutorialGenerate => run_tutorial_generate_mode(
            &source_dir,
            &tutorial_manifest,
            &tutorial_output,
            &repo_root,
        ),
        Mode::TutorialCheck => run_tutorial_check_mode(
            &source_dir,
            &tutorial_manifest,
            &tutorial_output,
            &repo_root,
        ),
        Mode::TutorialCatalogGenerate => run_tutorial_catalog_generate_mode(
            &tutorial_catalog_meta,
            &tutorial_source_dir,
            &tutorial_catalog,
        ),
        Mode::TutorialCatalogCheck => run_tutorial_catalog_check_mode(
            &tutorial_catalog_meta,
            &tutorial_source_dir,
            &tutorial_catalog,
        ),
    };
    if let Err(e) = result {
        eprintln!("{e}");
        process::exit(1);
    }
}
