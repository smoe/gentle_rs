//! Shared support helpers for GENtle command-line binaries.
//!
//! This module intentionally stays below the command-family layer: it contains
//! small argument and output helpers reused by multiple binaries without making
//! `gentle_cli`'s legacy direct-command surface the shared adapter API.

use crate::{engine::ProjectState, svg_pdf::SvgPdfRenderSummary, svg_png::SvgPngRenderSummary};
use serde_json::{Value, json};
use std::{fs, path::Path};

/// Parsed flags for binaries that accept at most one project path.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct SingleProjectCliArgs {
    /// Whether `--help` or `-h` was provided.
    pub show_help: bool,
    /// Whether `--version` or `-V` was provided.
    pub show_version: bool,
    /// Optional project path, either from `--project PATH` or one positional.
    pub project_path: Option<String>,
}

/// Policy knobs for [`parse_single_project_cli_args`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SingleProjectCliOptions {
    /// Whether to reject `--allow-screenshots` with GENtle's security-policy
    /// diagnostic. The GUI binary accepts this historical flag syntactically
    /// only to explain that screenshot capture remains disabled.
    pub reject_allow_screenshots: bool,
}

impl SingleProjectCliOptions {
    /// No extra flags beyond help/version/project/positional project.
    pub const fn standard() -> Self {
        Self {
            reject_allow_screenshots: false,
        }
    }

    /// Include the GUI's disabled screenshot-policy diagnostic.
    pub const fn gui() -> Self {
        Self {
            reject_allow_screenshots: true,
        }
    }
}

/// Parses the common `--help`, `--version`, and single-project-path pattern.
pub fn parse_single_project_cli_args(
    args: &[String],
    options: SingleProjectCliOptions,
) -> Result<SingleProjectCliArgs, String> {
    let mut parsed = SingleProjectCliArgs::default();
    let mut idx = 0usize;
    while idx < args.len() {
        match args[idx].as_str() {
            "--help" | "-h" => {
                parsed.show_help = true;
                idx += 1;
            }
            "--version" | "-V" => {
                parsed.show_version = true;
                idx += 1;
            }
            "--project" => {
                if idx + 1 >= args.len() {
                    return Err("Missing PATH after --project".to_string());
                }
                parsed.project_path = Some(args[idx + 1].clone());
                idx += 2;
            }
            "--allow-screenshots" if options.reject_allow_screenshots => {
                return Err("--allow-screenshots is disabled by security policy".to_string());
            }
            arg if arg.starts_with('-') => {
                return Err(format!("Unknown option '{arg}'"));
            }
            path => {
                if parsed.project_path.is_some() {
                    return Err(format!(
                        "Multiple project paths provided ('{}' and '{}')",
                        parsed.project_path.as_deref().unwrap_or_default(),
                        path
                    ));
                }
                parsed.project_path = Some(path.to_string());
                idx += 1;
            }
        }
    }
    Ok(parsed)
}

/// Parsed flags for stdio/stateful tool binaries such as `gentle_mcp`.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct StatePathCliArgs {
    /// Whether `--help` or `-h` was provided.
    pub show_help: bool,
    /// Whether `--version` or `-V` was provided.
    pub show_version: bool,
    /// Optional explicit state/project path.
    pub state_path: Option<String>,
}

/// Parses `--help`, `--version`, and optional `--state|--project PATH`.
pub fn parse_state_path_cli_args(args: &[String]) -> Result<StatePathCliArgs, String> {
    let mut parsed = StatePathCliArgs::default();
    let mut idx = 0usize;
    while idx < args.len() {
        match args[idx].as_str() {
            "--state" | "--project" => {
                if idx + 1 >= args.len() {
                    return Err(format!("Missing PATH after {}", args[idx]));
                }
                parsed.state_path = Some(args[idx + 1].clone());
                idx += 2;
            }
            "--help" | "-h" => {
                parsed.show_help = true;
                idx += 1;
            }
            "--version" | "-V" => {
                parsed.show_version = true;
                idx += 1;
            }
            other => {
                return Err(format!("Unknown argument '{other}'. Use --help for usage."));
            }
        }
    }
    Ok(parsed)
}

/// Loads a project state from `path`, or returns an empty default state when
/// the path does not exist yet.
pub fn load_state_or_default(path: &str) -> Result<ProjectState, String> {
    if !Path::new(path).exists() {
        return Ok(ProjectState::default());
    }
    let raw =
        fs::read_to_string(path).map_err(|e| format!("Could not read state file '{path}': {e}"))?;
    if raw.trim().is_empty() {
        return Ok(ProjectState::default());
    }
    ProjectState::load_from_path(path).map_err(|e| e.to_string())
}

/// Builds the standard machine-readable `svg-png` JSON summary used by CLI
/// binaries that rasterize an SVG into a PNG.
pub fn svg_png_summary_json(summary: &SvgPngRenderSummary) -> Value {
    json!({
        "status": "ok",
        "mode": "svg-png",
        "input_path": &summary.input_path,
        "output_path": &summary.output_path,
        "scale": &summary.scale,
        "drop_dotplot_metadata": summary.drop_dotplot_metadata,
        "width": summary.width,
        "height": summary.height,
        "font_face_count": summary.font_face_count,
    })
}

/// Builds the standard machine-readable `svg-pdf` JSON summary used by CLI
/// binaries that render an SVG into a one-page PDF.
pub fn svg_pdf_summary_json(summary: &SvgPdfRenderSummary) -> Value {
    json!({
        "status": "ok",
        "mode": "svg-pdf",
        "input_path": &summary.input_path,
        "output_path": &summary.output_path,
        "scale": &summary.scale,
        "drop_dotplot_metadata": summary.drop_dotplot_metadata,
        "width": summary.width,
        "height": summary.height,
        "font_face_count": summary.font_face_count,
        "page_width_pt": &summary.page_width_pt,
        "page_height_pt": &summary.page_height_pt,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    fn argv(values: &[&str]) -> Vec<String> {
        values.iter().map(|value| value.to_string()).collect()
    }

    #[test]
    fn single_project_parser_accepts_project_flag() {
        let parsed = parse_single_project_cli_args(
            &argv(&["--project", "demo.gentle.json"]),
            SingleProjectCliOptions::standard(),
        )
        .expect("parse project args");
        assert_eq!(parsed.project_path.as_deref(), Some("demo.gentle.json"));
        assert!(!parsed.show_help);
        assert!(!parsed.show_version);
    }

    #[test]
    fn single_project_parser_rejects_second_project_path() {
        let err = parse_single_project_cli_args(
            &argv(&["--project", "a.gentle.json", "b.gentle.json"]),
            SingleProjectCliOptions::standard(),
        )
        .expect_err("duplicate project paths should fail");
        assert_eq!(
            err,
            "Multiple project paths provided ('a.gentle.json' and 'b.gentle.json')"
        );
    }

    #[test]
    fn single_project_parser_preserves_gui_screenshot_diagnostic() {
        let err = parse_single_project_cli_args(
            &argv(&["--allow-screenshots"]),
            SingleProjectCliOptions::gui(),
        )
        .expect_err("disabled screenshot flag should fail");
        assert_eq!(err, "--allow-screenshots is disabled by security policy");
    }

    #[test]
    fn state_path_parser_accepts_project_alias() {
        let parsed = parse_state_path_cli_args(&argv(&["--project", "project.gentle.json"]))
            .expect("parse state args");
        assert_eq!(parsed.state_path.as_deref(), Some("project.gentle.json"));
    }

    #[test]
    fn load_state_or_default_accepts_missing_and_empty_paths() {
        let td = tempdir().expect("tempdir");
        let missing = td.path().join("missing.gentle.json");
        let missing_state =
            load_state_or_default(&missing.to_string_lossy()).expect("missing path defaults");
        assert!(missing_state.sequences.is_empty());

        let empty = td.path().join("empty.gentle.json");
        fs::write(&empty, "\n").expect("write empty state");
        let empty_state =
            load_state_or_default(&empty.to_string_lossy()).expect("empty path defaults");
        assert!(empty_state.sequences.is_empty());
    }

    #[test]
    fn svg_png_summary_json_keeps_shared_shape() {
        let summary = SvgPngRenderSummary {
            input_path: "in.svg".to_string(),
            output_path: "out.png".to_string(),
            scale: "2".to_string(),
            drop_dotplot_metadata: true,
            width: 100,
            height: 50,
            font_face_count: 7,
        };
        let value = svg_png_summary_json(&summary);
        assert_eq!(value["status"], "ok");
        assert_eq!(value["mode"], "svg-png");
        assert_eq!(value["input_path"], "in.svg");
        assert_eq!(value["drop_dotplot_metadata"], true);
    }

    #[test]
    fn svg_pdf_summary_json_keeps_shared_shape() {
        let summary = SvgPdfRenderSummary {
            input_path: "in.svg".to_string(),
            output_path: "out.pdf".to_string(),
            scale: "2".to_string(),
            drop_dotplot_metadata: true,
            width: 100,
            height: 50,
            font_face_count: 7,
            page_width_pt: "75.00".to_string(),
            page_height_pt: "37.50".to_string(),
        };
        let value = svg_pdf_summary_json(&summary);
        assert_eq!(value["status"], "ok");
        assert_eq!(value["mode"], "svg-pdf");
        assert_eq!(value["output_path"], "out.pdf");
        assert_eq!(value["page_width_pt"], "75.00");
    }
}
