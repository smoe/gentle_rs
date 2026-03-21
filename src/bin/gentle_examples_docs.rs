//! Documentation helper binary that renders workflow examples into generated snippets.

use gentle::workflow_examples::{
    DEFAULT_TUTORIAL_CATALOG_META_PATH, DEFAULT_TUTORIAL_CATALOG_PATH,
    DEFAULT_TUTORIAL_MANIFEST_PATH, DEFAULT_TUTORIAL_OUTPUT_DIR, DEFAULT_TUTORIAL_SOURCE_DIR,
    DEFAULT_WORKFLOW_EXAMPLE_DIR, DEFAULT_WORKFLOW_SNIPPET_DIR, check_tutorial_catalog_generated,
    check_tutorial_generated, check_tutorial_manifest_generated, generate_tutorial_docs,
    generate_workflow_example_docs, load_workflow_examples, validate_example_required_files,
    write_tutorial_catalog_from_sources, write_tutorial_manifest_from_sources,
};
use regex::Regex;
use resvg::{self, tiny_skia, usvg};
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
    SvgPng,
    TutorialGenerate,
    TutorialCheck,
    TutorialCatalogGenerate,
    TutorialCatalogCheck,
    TutorialManifestGenerate,
    TutorialManifestCheck,
}

#[derive(Debug)]
struct CliArgs {
    show_help: bool,
    mode: Mode,
    source_dir: String,
    example_output_dir: String,
    svg_input: String,
    png_output: String,
    svg_scale: f32,
    drop_dotplot_metadata: bool,
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
            svg_input: String::new(),
            png_output: String::new(),
            svg_scale: 1.0,
            drop_dotplot_metadata: false,
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
gentle_examples_docs svg-png INPUT.svg OUTPUT.png [--scale N] [--drop-dotplot-metadata]\n  \
gentle_examples_docs tutorial-generate [--source DIR] [--manifest FILE] [--tutorial-output DIR] [--repo-root DIR]\n  \
gentle_examples_docs tutorial-check [--source DIR] [--manifest FILE] [--tutorial-output DIR] [--repo-root DIR]\n  \
gentle_examples_docs tutorial-catalog-generate [--catalog-meta FILE] [--tutorial-sources DIR] [--tutorial-catalog FILE]\n  \
gentle_examples_docs tutorial-catalog-check [--catalog-meta FILE] [--tutorial-sources DIR] [--tutorial-catalog FILE]\n  \
gentle_examples_docs tutorial-manifest-generate [--catalog-meta FILE] [--tutorial-sources DIR] [--manifest FILE]\n  \
gentle_examples_docs tutorial-manifest-check [--catalog-meta FILE] [--tutorial-sources DIR] [--manifest FILE]\n\n  \
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
            "--scale" => {
                if idx + 1 >= args.len() {
                    return Err("Missing N after --scale".to_string());
                }
                parsed.svg_scale = args[idx + 1]
                    .parse::<f32>()
                    .map_err(|e| format!("Could not parse --scale '{}': {e}", args[idx + 1]))?;
                idx += 2;
            }
            "--drop-dotplot-metadata" => {
                parsed.drop_dotplot_metadata = true;
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
            "svg-png" => {
                parsed.mode = Mode::SvgPng;
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
            "tutorial-manifest-generate" => {
                parsed.mode = Mode::TutorialManifestGenerate;
                idx += 1;
            }
            "tutorial-manifest-check" => {
                parsed.mode = Mode::TutorialManifestCheck;
                idx += 1;
            }
            other => {
                if parsed.mode == Mode::SvgPng && parsed.svg_input.is_empty() {
                    parsed.svg_input = other.to_string();
                    idx += 1;
                    continue;
                }
                if parsed.mode == Mode::SvgPng && parsed.png_output.is_empty() {
                    parsed.png_output = other.to_string();
                    idx += 1;
                    continue;
                }
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

fn strip_dotplot_metadata_text(svg: &str) -> Result<String, String> {
    let metadata_re = Regex::new(
        r#"<text\b[^>]*>(?:Dotplot workspace export:[^<]*|rendered_cells=[^<]*|x:\s[^<]*|y:\s[^<]*|GENtle dotplot SVG export)</text>"#,
    )
    .map_err(|e| format!("Could not compile dotplot metadata regex: {e}"))?;
    Ok(metadata_re.replace_all(svg, "").into_owned())
}

fn run_svg_png_mode(
    input_path: &Path,
    output_path: &Path,
    scale: f32,
    drop_dotplot_metadata: bool,
) -> Result<(), String> {
    if input_path.as_os_str().is_empty() {
        return Err("svg-png requires INPUT.svg".to_string());
    }
    if output_path.as_os_str().is_empty() {
        return Err("svg-png requires OUTPUT.png".to_string());
    }
    if !(scale.is_finite() && scale > 0.0) {
        return Err(format!(
            "svg-png requires a positive finite --scale value, got {scale}"
        ));
    }

    let mut svg_text = std::fs::read_to_string(input_path)
        .map_err(|e| format!("Could not read SVG '{}': {e}", input_path.display()))?;
    if drop_dotplot_metadata {
        svg_text = strip_dotplot_metadata_text(&svg_text)?;
    }

    let mut opt = usvg::Options {
        resources_dir: std::fs::canonicalize(input_path)
            .ok()
            .and_then(|path| path.parent().map(|parent| parent.to_path_buf())),
        ..usvg::Options::default()
    };
    opt.fontdb_mut().load_system_fonts();

    let tree = usvg::Tree::from_str(&svg_text, &opt)
        .map_err(|e| format!("Could not parse SVG '{}': {e}", input_path.display()))?;
    let pixmap_size = tree
        .size()
        .to_int_size()
        .scale_by(scale)
        .ok_or_else(|| format!("Could not scale SVG size by {scale}"))?;
    let mut pixmap = tiny_skia::Pixmap::new(pixmap_size.width(), pixmap_size.height())
        .ok_or_else(|| {
            format!(
                "Could not allocate PNG canvas {}x{}",
                pixmap_size.width(),
                pixmap_size.height()
            )
        })?;
    let transform = if (scale - 1.0).abs() <= f32::EPSILON {
        tiny_skia::Transform::default()
    } else {
        tiny_skia::Transform::from_scale(scale, scale)
    };
    resvg::render(&tree, transform, &mut pixmap.as_mut());
    pixmap
        .save_png(output_path)
        .map_err(|e| format!("Could not write PNG '{}': {e}", output_path.display()))?;

    let summary = json!({
        "status": "ok",
        "mode": "svg-png",
        "input_path": input_path.to_string_lossy(),
        "output_path": output_path.to_string_lossy(),
        "scale": scale,
        "drop_dotplot_metadata": drop_dotplot_metadata,
        "width": pixmap_size.width(),
        "height": pixmap_size.height(),
    });
    let pretty = serde_json::to_string_pretty(&summary)
        .map_err(|e| format!("Could not serialize svg-png summary: {e}"))?;
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

fn run_tutorial_manifest_generate_mode(
    meta_path: &Path,
    source_dir: &Path,
    output_path: &Path,
) -> Result<(), String> {
    let manifest = write_tutorial_manifest_from_sources(meta_path, source_dir, output_path)?;
    let summary = json!({
        "status": "ok",
        "mode": "tutorial-manifest-generate",
        "chapter_count": manifest.chapters.len(),
        "output_path": output_path.to_string_lossy(),
    });
    let pretty = serde_json::to_string_pretty(&summary)
        .map_err(|e| format!("Could not serialize tutorial manifest generation report: {e}"))?;
    println!("{pretty}");
    Ok(())
}

fn run_tutorial_manifest_check_mode(
    meta_path: &Path,
    source_dir: &Path,
    output_path: &Path,
) -> Result<(), String> {
    let manifest = check_tutorial_manifest_generated(meta_path, source_dir, output_path)?;
    let summary = json!({
        "status": "ok",
        "mode": "tutorial-manifest-check",
        "chapter_count": manifest.chapters.len(),
        "output_path": output_path.to_string_lossy(),
    });
    let pretty = serde_json::to_string_pretty(&summary)
        .map_err(|e| format!("Could not serialize tutorial manifest check report: {e}"))?;
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
    let svg_input = PathBuf::from(parsed.svg_input);
    let png_output = PathBuf::from(parsed.png_output);
    let result = match parsed.mode {
        Mode::ExampleGenerate => run_generate_mode(&source_dir, &example_output),
        Mode::ExampleCheck => run_check_mode(&source_dir),
        Mode::SvgPng => run_svg_png_mode(
            &svg_input,
            &png_output,
            parsed.svg_scale,
            parsed.drop_dotplot_metadata,
        ),
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
        Mode::TutorialManifestGenerate => run_tutorial_manifest_generate_mode(
            &tutorial_catalog_meta,
            &tutorial_source_dir,
            &tutorial_manifest,
        ),
        Mode::TutorialManifestCheck => run_tutorial_manifest_check_mode(
            &tutorial_catalog_meta,
            &tutorial_source_dir,
            &tutorial_manifest,
        ),
    };
    if let Err(e) = result {
        eprintln!("{e}");
        process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::{Mode, parse_args, strip_dotplot_metadata_text};

    #[test]
    fn parse_svg_png_mode_with_cleanup_flag() {
        let args = vec![
            "svg-png".to_string(),
            "input.svg".to_string(),
            "output.png".to_string(),
            "--scale".to_string(),
            "1.5".to_string(),
            "--drop-dotplot-metadata".to_string(),
        ];
        let parsed = parse_args(&args).expect("parse svg-png args");
        assert_eq!(parsed.mode, Mode::SvgPng);
        assert_eq!(parsed.svg_input, "input.svg");
        assert_eq!(parsed.png_output, "output.png");
        assert!((parsed.svg_scale - 1.5).abs() < 1e-6);
        assert!(parsed.drop_dotplot_metadata);
    }

    #[test]
    fn strip_dotplot_metadata_preserves_axis_labels() {
        let svg = concat!(
            "<svg>",
            "<text>Dotplot workspace export: dotplot_primary</text>",
            "<text>rendered_cells=42 sampled_points=77 sample_stride=1</text>",
            "<text>x: tp73_cdna</text>",
            "<text>y: tp73_genomic</text>",
            "<text>1</text>",
            "<text>5026</text>",
            "<text>GENtle dotplot SVG export</text>",
            "</svg>"
        );
        let cleaned = strip_dotplot_metadata_text(svg).expect("strip metadata");
        assert!(!cleaned.contains("Dotplot workspace export:"));
        assert!(!cleaned.contains("rendered_cells="));
        assert!(!cleaned.contains("x: tp73_cdna"));
        assert!(!cleaned.contains("y: tp73_genomic"));
        assert!(!cleaned.contains("GENtle dotplot SVG export"));
        assert!(cleaned.contains("<text>1</text>"));
        assert!(cleaned.contains("<text>5026</text>"));
    }
}
