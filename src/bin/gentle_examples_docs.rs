use gentle::workflow_examples::{
    generate_workflow_example_docs, load_workflow_examples, validate_example_required_files,
    DEFAULT_WORKFLOW_EXAMPLE_DIR, DEFAULT_WORKFLOW_SNIPPET_DIR,
};
use serde_json::json;
use std::{
    env,
    path::{Path, PathBuf},
    process,
};

#[derive(Debug, Default)]
struct CliArgs {
    show_help: bool,
    check_only: bool,
    source_dir: String,
    output_dir: String,
}

fn usage() {
    eprintln!(
        "Usage:\n  \
gentle_examples_docs [generate] [--source DIR] [--output DIR]\n  \
gentle_examples_docs --check [--source DIR]\n\n  \
Defaults:\n  \
  --source {}\n  \
  --output {}",
        DEFAULT_WORKFLOW_EXAMPLE_DIR, DEFAULT_WORKFLOW_SNIPPET_DIR
    );
}

fn parse_args(args: &[String]) -> Result<CliArgs, String> {
    let mut parsed = CliArgs {
        source_dir: DEFAULT_WORKFLOW_EXAMPLE_DIR.to_string(),
        output_dir: DEFAULT_WORKFLOW_SNIPPET_DIR.to_string(),
        ..CliArgs::default()
    };
    let mut idx = 0usize;
    while idx < args.len() {
        match args[idx].as_str() {
            "--help" | "-h" => {
                parsed.show_help = true;
                idx += 1;
            }
            "--check" => {
                parsed.check_only = true;
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
                parsed.output_dir = args[idx + 1].clone();
                idx += 2;
            }
            "generate" => {
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
    let output_dir = PathBuf::from(parsed.output_dir);
    let result = if parsed.check_only {
        run_check_mode(&source_dir)
    } else {
        run_generate_mode(&source_dir, &output_dir)
    };
    if let Err(e) = result {
        eprintln!("{e}");
        process::exit(1);
    }
}
