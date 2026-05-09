//! Direct rendering and SVG/PNG export CLI command handling kept out of `run()`.
//!
//! These routes are the direct CLI's rendering surface: format conversion,
//! sequence/dotplot/RNA/lineage SVG export, and RNA structure inspection. They
//! stay together with the dotplot overlay argument parser because that parser is
//! only meaningful for the direct render-dotplot command.

use gentle::{
    cli_support::svg_png_summary_json,
    engine::{
        DotplotOverlayAnchorExonRef, DotplotOverlayXAxisMode, GentleEngine, Operation,
        RenderSvgMode,
    },
    svg_png::{SvgPngRenderOptions, render_svg_file_to_png},
};
use std::path::Path;

use super::*;

const RENDERING_COMMANDS: &[&str] = &[
    "svg-png",
    "render-svg",
    "render-dotplot-svg",
    "render-rna-svg",
    "rna-info",
    "render-lineage-svg",
];

pub(super) fn is_rendering_command(command: &str) -> bool {
    RENDERING_COMMANDS.contains(&command)
}

pub(super) fn handle_rendering_family(
    command: &str,
    args: &[String],
    cmd_idx: usize,
    state_path: &str,
) -> Result<(), String> {
    match command {
        "svg-png" => handle_svg_png(args, cmd_idx),
        "render-svg" => handle_render_svg(args, cmd_idx, state_path),
        "render-dotplot-svg" => handle_render_dotplot_svg(args, cmd_idx, state_path),
        "render-rna-svg" => handle_render_rna_svg(args, cmd_idx, state_path),
        "rna-info" => handle_rna_info(args, cmd_idx, state_path),
        "render-lineage-svg" => handle_render_lineage_svg(args, cmd_idx, state_path),
        _ => Err(format!("Expected rendering command, got '{command}'")),
    }
}

fn handle_svg_png(args: &[String], cmd_idx: usize) -> Result<(), String> {
    if args.len() <= cmd_idx + 2 {
        usage();
        return Err(
            "svg-png requires: INPUT.svg OUTPUT.png [--scale N] [--drop-dotplot-metadata]"
                .to_string(),
        );
    }
    let input = &args[cmd_idx + 1];
    let output = &args[cmd_idx + 2];
    let mut scale = 1.0f32;
    let mut drop_dotplot_metadata = false;
    let mut idx = cmd_idx + 3;
    while idx < args.len() {
        match args[idx].as_str() {
            "--scale" => {
                if idx + 1 >= args.len() {
                    return Err("Missing N after --scale".to_string());
                }
                scale = args[idx + 1]
                    .parse::<f32>()
                    .map_err(|e| format!("Could not parse --scale '{}': {e}", args[idx + 1]))?;
                idx += 2;
            }
            "--drop-dotplot-metadata" => {
                drop_dotplot_metadata = true;
                idx += 1;
            }
            other => {
                return Err(format!("Unknown option '{}' for svg-png", other));
            }
        }
    }
    ensure_parent_dir(output)?;
    let summary = render_svg_file_to_png(
        Path::new(input),
        Path::new(output),
        SvgPngRenderOptions {
            scale,
            drop_dotplot_metadata,
        },
    )?;
    print_json(&svg_png_summary_json(&summary))
}

fn handle_render_svg(args: &[String], cmd_idx: usize, state_path: &str) -> Result<(), String> {
    if args.len() <= cmd_idx + 3 {
        usage();
        return Err("render-svg requires: SEQ_ID linear|circular OUTPUT.svg".to_string());
    }
    let seq_id = &args[cmd_idx + 1];
    let mode = &args[cmd_idx + 2];
    let output = &args[cmd_idx + 3];
    let mode = match mode.as_str() {
        "linear" => RenderSvgMode::Linear,
        "circular" => RenderSvgMode::Circular,
        _ => {
            return Err(format!(
                "Unknown render mode '{mode}', expected 'linear' or 'circular'"
            ));
        }
    };
    let mut engine = GentleEngine::from_state(load_state(state_path)?);
    let result = engine
        .apply(Operation::RenderSequenceSvg {
            seq_id: seq_id.to_string(),
            mode,
            path: output.to_string(),
        })
        .map_err(|e| e.to_string())?;
    save_state_and_print_first_message(&engine, state_path, &result.messages)
}

fn handle_render_dotplot_svg(
    args: &[String],
    cmd_idx: usize,
    state_path: &str,
) -> Result<(), String> {
    if args.len() <= cmd_idx + 3 {
        usage();
        return Err(
            "render-dotplot-svg requires: SEQ_ID DOTPLOT_ID OUTPUT.svg [--flex-track ID] [--display-threshold N] [--intensity-gain N] [--overlay-x-axis percent_length|left_aligned_bp|right_aligned_bp|shared_exon_anchor|query_anchor_bp] [--overlay-anchor-exon START..END]".to_string(),
        );
    }
    let seq_id = args[cmd_idx + 1].trim();
    let dotplot_id = args[cmd_idx + 2].trim();
    if seq_id.is_empty() {
        return Err("render-dotplot-svg requires non-empty SEQ_ID".to_string());
    }
    if dotplot_id.is_empty() {
        return Err("render-dotplot-svg requires non-empty DOTPLOT_ID".to_string());
    }
    let output = &args[cmd_idx + 3];
    let mut flex_track_id: Option<String> = None;
    let mut display_density_threshold: Option<f32> = None;
    let mut display_intensity_gain: Option<f32> = None;
    let mut overlay_x_axis_mode = DotplotOverlayXAxisMode::PercentLength;
    let mut overlay_anchor_exon: Option<DotplotOverlayAnchorExonRef> = None;
    let mut idx = cmd_idx + 4;
    while idx < args.len() {
        match args[idx].as_str() {
            "--flex-track" => {
                if idx + 1 >= args.len() {
                    return Err("Missing value after --flex-track".to_string());
                }
                let value = args[idx + 1].trim();
                if !value.is_empty() {
                    flex_track_id = Some(value.to_string());
                }
                idx += 2;
            }
            "--display-threshold" => {
                if idx + 1 >= args.len() {
                    return Err("Missing value after --display-threshold".to_string());
                }
                let raw = args[idx + 1].trim();
                let parsed = raw.parse::<f32>().map_err(|_| {
                    format!(
                        "Invalid --display-threshold '{}': expected decimal number",
                        raw
                    )
                })?;
                display_density_threshold = Some(parsed);
                idx += 2;
            }
            "--intensity-gain" => {
                if idx + 1 >= args.len() {
                    return Err("Missing value after --intensity-gain".to_string());
                }
                let raw = args[idx + 1].trim();
                let parsed = raw.parse::<f32>().map_err(|_| {
                    format!(
                        "Invalid --intensity-gain '{}': expected decimal number",
                        raw
                    )
                })?;
                display_intensity_gain = Some(parsed);
                idx += 2;
            }
            "--overlay-x-axis" => {
                if idx + 1 >= args.len() {
                    return Err("Missing value after --overlay-x-axis".to_string());
                }
                overlay_x_axis_mode = parse_dotplot_overlay_x_axis_mode_arg(&args[idx + 1])?;
                idx += 2;
            }
            "--overlay-anchor-exon" => {
                if idx + 1 >= args.len() {
                    return Err("Missing value after --overlay-anchor-exon".to_string());
                }
                overlay_anchor_exon = Some(DotplotOverlayAnchorExonRef::parse(&args[idx + 1])?);
                idx += 2;
            }
            other => {
                return Err(format!(
                    "Unknown argument '{other}' for render-dotplot-svg (expected --flex-track/--display-threshold/--intensity-gain/--overlay-x-axis/--overlay-anchor-exon)"
                ));
            }
        }
    }
    let mut engine = GentleEngine::from_state(load_state(state_path)?);
    let result = engine
        .apply(Operation::RenderDotplotSvg {
            seq_id: seq_id.to_string(),
            dotplot_id: dotplot_id.to_string(),
            path: output.to_string(),
            flex_track_id,
            display_density_threshold,
            display_intensity_gain,
            overlay_x_axis_mode,
            overlay_anchor_exon,
        })
        .map_err(|e| e.to_string())?;
    save_state_and_print_first_message(&engine, state_path, &result.messages)
}

fn handle_render_rna_svg(args: &[String], cmd_idx: usize, state_path: &str) -> Result<(), String> {
    if args.len() <= cmd_idx + 2 {
        usage();
        return Err("render-rna-svg requires: SEQ_ID OUTPUT.svg".to_string());
    }
    let seq_id = &args[cmd_idx + 1];
    let output = &args[cmd_idx + 2];
    let mut engine = GentleEngine::from_state(load_state(state_path)?);
    let result = engine
        .apply(Operation::RenderRnaStructureSvg {
            seq_id: seq_id.to_string(),
            path: output.to_string(),
        })
        .map_err(|e| e.to_string())?;
    save_state_and_print_first_message(&engine, state_path, &result.messages)
}

fn handle_rna_info(args: &[String], cmd_idx: usize, state_path: &str) -> Result<(), String> {
    if args.len() <= cmd_idx + 1 {
        usage();
        return Err("rna-info requires: SEQ_ID".to_string());
    }
    let seq_id = &args[cmd_idx + 1];
    let engine = GentleEngine::from_state(load_state(state_path)?);
    let report = engine
        .inspect_rna_structure(seq_id)
        .map_err(|e| e.to_string())?;
    print_json(&report)
}

fn handle_render_lineage_svg(
    args: &[String],
    cmd_idx: usize,
    state_path: &str,
) -> Result<(), String> {
    if args.len() <= cmd_idx + 1 {
        usage();
        return Err("render-lineage-svg requires: OUTPUT.svg".to_string());
    }
    let output = &args[cmd_idx + 1];
    let mut engine = GentleEngine::from_state(load_state(state_path)?);
    let result = engine
        .apply(Operation::RenderLineageSvg {
            path: output.to_string(),
        })
        .map_err(|e| e.to_string())?;
    save_state_and_print_first_message(&engine, state_path, &result.messages)
}

fn parse_dotplot_overlay_x_axis_mode_arg(raw: &str) -> Result<DotplotOverlayXAxisMode, String> {
    match raw.trim() {
        "percent_length" => Ok(DotplotOverlayXAxisMode::PercentLength),
        "left_aligned_bp" => Ok(DotplotOverlayXAxisMode::LeftAlignedBp),
        "right_aligned_bp" => Ok(DotplotOverlayXAxisMode::RightAlignedBp),
        "shared_exon_anchor" => Ok(DotplotOverlayXAxisMode::SharedExonAnchor),
        "query_anchor_bp" => Ok(DotplotOverlayXAxisMode::QueryAnchorBp),
        other => Err(format!(
            "Invalid --overlay-x-axis '{}': expected percent_length, left_aligned_bp, right_aligned_bp, shared_exon_anchor, or query_anchor_bp",
            other
        )),
    }
}

fn save_state_and_print_first_message(
    engine: &GentleEngine,
    state_path: &str,
    messages: &[String],
) -> Result<(), String> {
    engine
        .state()
        .save_to_path(state_path)
        .map_err(|e| e.to_string())?;
    if let Some(msg) = messages.first() {
        println!("{msg}");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    fn argv(values: &[&str]) -> Vec<String> {
        values.iter().map(|value| value.to_string()).collect()
    }

    #[test]
    fn rendering_command_filter_recognizes_moved_routes() {
        assert!(is_rendering_command("svg-png"));
        assert!(is_rendering_command("render-rna-svg"));
        assert!(is_rendering_command("rna-info"));
        assert!(!is_rendering_command("protocol-cartoon"));
    }

    #[test]
    fn svg_png_missing_args_reports_legacy_message() {
        let args = argv(&["gentle_cli", "svg-png"]);
        let err = handle_rendering_family("svg-png", &args, 1, ".gentle_state.json")
            .expect_err("missing args should fail");
        assert_eq!(
            err,
            "svg-png requires: INPUT.svg OUTPUT.png [--scale N] [--drop-dotplot-metadata]"
        );
    }

    #[test]
    fn render_svg_rejects_unknown_mode() {
        let args = argv(&["gentle_cli", "render-svg", "seq1", "flat", "out.svg"]);
        let err = handle_rendering_family("render-svg", &args, 1, ".gentle_state.json")
            .expect_err("unknown mode should fail");
        assert_eq!(
            err,
            "Unknown render mode 'flat', expected 'linear' or 'circular'"
        );
    }

    #[test]
    fn render_dotplot_svg_rejects_empty_seq_id() {
        let args = argv(&[
            "gentle_cli",
            "render-dotplot-svg",
            "",
            "dotplot-1",
            "out.svg",
        ]);
        let err = handle_rendering_family("render-dotplot-svg", &args, 1, ".gentle_state.json")
            .expect_err("empty seq id should fail");
        assert_eq!(err, "render-dotplot-svg requires non-empty SEQ_ID");
    }

    #[test]
    fn parse_dotplot_overlay_x_axis_mode_arg_accepts_shared_exon_anchor() {
        assert_eq!(
            parse_dotplot_overlay_x_axis_mode_arg("shared_exon_anchor")
                .expect("parse shared exon anchor"),
            DotplotOverlayXAxisMode::SharedExonAnchor
        );
    }

    #[test]
    fn parse_dotplot_overlay_x_axis_mode_arg_accepts_query_anchor_bp() {
        assert_eq!(
            parse_dotplot_overlay_x_axis_mode_arg("query_anchor_bp")
                .expect("parse query anchor mode"),
            DotplotOverlayXAxisMode::QueryAnchorBp
        );
    }

    #[test]
    fn render_lineage_svg_writes_svg_and_state() {
        let td = tempdir().expect("tempdir");
        let state_path = td.path().join("state.json");
        let output_path = td.path().join("lineage.svg");
        let output = output_path.to_string_lossy().to_string();
        let args = argv(&["gentle_cli", "render-lineage-svg", &output]);
        handle_rendering_family(
            "render-lineage-svg",
            &args,
            1,
            &state_path.to_string_lossy(),
        )
        .expect("lineage render should execute");
        assert!(output_path.exists(), "lineage route should write SVG");
        assert!(state_path.exists(), "render route should persist state");
    }
}
