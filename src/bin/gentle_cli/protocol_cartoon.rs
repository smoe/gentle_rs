//! Direct `protocol-cartoon` CLI command handling kept out of `run()`.
//!
//! The shared shell already owns protocol-cartoon parsing and execution for
//! shell adapters. This module preserves the legacy direct CLI route and its
//! exact diagnostics while moving another self-contained command family out of
//! the monolithic binary entry point.

use gentle::protocol_cartoon::{ProtocolCartoonKind, protocol_cartoon_catalog_rows};

use super::*;

pub(super) const USAGE: &str = "\
  gentle_cli [--state PATH|--project PATH] protocol-cartoon list
  gentle_cli [--state PATH|--project PATH] protocol-cartoon render-svg PROTOCOL_ID OUTPUT.svg
  gentle_cli [--state PATH|--project PATH] protocol-cartoon render-template-svg TEMPLATE.json OUTPUT.svg
  gentle_cli [--state PATH|--project PATH] protocol-cartoon template-validate TEMPLATE.json
  gentle_cli [--state PATH|--project PATH] protocol-cartoon render-with-bindings TEMPLATE.json BINDINGS.json OUTPUT.svg
  gentle_cli [--state PATH|--project PATH] protocol-cartoon template-export PROTOCOL_ID OUTPUT.json

";

const SUBCOMMANDS: &str = "list, render-svg, render-template-svg, template-validate, render-with-bindings, template-export";

pub(super) fn handle_protocol_cartoon_family(
    args: &[String],
    cmd_idx: usize,
    state_path: &str,
) -> Result<(), String> {
    if args.len() <= cmd_idx + 1 {
        usage();
        return Err(format!(
            "protocol-cartoon requires a subcommand: {SUBCOMMANDS}"
        ));
    }
    match args[cmd_idx + 1].as_str() {
        "list" => {
            if args.len() != cmd_idx + 2 {
                return Err("protocol-cartoon list takes no additional arguments".to_string());
            }
            print_json(&protocol_cartoon_catalog_rows())
        }
        "render-svg" => {
            if args.len() != cmd_idx + 4 {
                return Err(
                    "protocol-cartoon render-svg requires: PROTOCOL_ID OUTPUT.svg".to_string(),
                );
            }
            let protocol_id = args[cmd_idx + 2].trim();
            if protocol_id.is_empty() {
                return Err(
                    "protocol-cartoon render-svg requires non-empty PROTOCOL_ID".to_string()
                );
            }
            let protocol = parse_protocol_id(protocol_id)?;
            let output = &args[cmd_idx + 3];
            let mut engine = GentleEngine::from_state(load_state(state_path)?);
            let result = engine
                .apply(Operation::RenderProtocolCartoonSvg {
                    protocol,
                    path: output.to_string(),
                })
                .map_err(|e| e.to_string())?;
            save_state_and_print_first_message(&engine, state_path, &result.messages)
        }
        "render-template-svg" => {
            if args.len() != cmd_idx + 4 {
                return Err(
                    "protocol-cartoon render-template-svg requires: TEMPLATE.json OUTPUT.svg"
                        .to_string(),
                );
            }
            let template_path = args[cmd_idx + 2].trim();
            if template_path.is_empty() {
                return Err(
                    "protocol-cartoon render-template-svg requires non-empty TEMPLATE.json"
                        .to_string(),
                );
            }
            let output = &args[cmd_idx + 3];
            let mut engine = GentleEngine::from_state(load_state(state_path)?);
            let result = engine
                .apply(Operation::RenderProtocolCartoonTemplateSvg {
                    template_path: args[cmd_idx + 2].clone(),
                    path: output.to_string(),
                })
                .map_err(|e| e.to_string())?;
            save_state_and_print_first_message(&engine, state_path, &result.messages)
        }
        "template-validate" => {
            if args.len() != cmd_idx + 3 {
                return Err(
                    "protocol-cartoon template-validate requires: TEMPLATE.json".to_string()
                );
            }
            let template_path = args[cmd_idx + 2].trim();
            if template_path.is_empty() {
                return Err(
                    "protocol-cartoon template-validate requires non-empty TEMPLATE.json"
                        .to_string(),
                );
            }
            let mut engine = GentleEngine::from_state(load_state(state_path)?);
            let result = engine
                .apply(Operation::ValidateProtocolCartoonTemplate {
                    template_path: args[cmd_idx + 2].clone(),
                })
                .map_err(|e| e.to_string())?;
            save_state_and_print_first_message(&engine, state_path, &result.messages)
        }
        "render-with-bindings" => {
            if args.len() != cmd_idx + 5 {
                return Err(
                    "protocol-cartoon render-with-bindings requires: TEMPLATE.json BINDINGS.json OUTPUT.svg"
                        .to_string(),
                );
            }
            let template_path = args[cmd_idx + 2].trim();
            if template_path.is_empty() {
                return Err(
                    "protocol-cartoon render-with-bindings requires non-empty TEMPLATE.json"
                        .to_string(),
                );
            }
            let bindings_path = args[cmd_idx + 3].trim();
            if bindings_path.is_empty() {
                return Err(
                    "protocol-cartoon render-with-bindings requires non-empty BINDINGS.json"
                        .to_string(),
                );
            }
            let output = &args[cmd_idx + 4];
            let mut engine = GentleEngine::from_state(load_state(state_path)?);
            let result = engine
                .apply(Operation::RenderProtocolCartoonTemplateWithBindingsSvg {
                    template_path: args[cmd_idx + 2].clone(),
                    bindings_path: args[cmd_idx + 3].clone(),
                    path: output.to_string(),
                })
                .map_err(|e| e.to_string())?;
            save_state_and_print_first_message(&engine, state_path, &result.messages)
        }
        "template-export" => {
            if args.len() != cmd_idx + 4 {
                return Err(
                    "protocol-cartoon template-export requires: PROTOCOL_ID OUTPUT.json"
                        .to_string(),
                );
            }
            let protocol_id = args[cmd_idx + 2].trim();
            if protocol_id.is_empty() {
                return Err(
                    "protocol-cartoon template-export requires non-empty PROTOCOL_ID".to_string(),
                );
            }
            let protocol = parse_protocol_id(protocol_id)?;
            let output = &args[cmd_idx + 3];
            let mut engine = GentleEngine::from_state(load_state(state_path)?);
            let result = engine
                .apply(Operation::ExportProtocolCartoonTemplateJson {
                    protocol,
                    path: output.to_string(),
                })
                .map_err(|e| e.to_string())?;
            save_state_and_print_first_message(&engine, state_path, &result.messages)
        }
        other => Err(format!(
            "Unknown protocol-cartoon subcommand '{other}' (expected {SUBCOMMANDS})"
        )),
    }
}

fn parse_protocol_id(protocol_id: &str) -> Result<ProtocolCartoonKind, String> {
    ProtocolCartoonKind::parse_id(protocol_id).ok_or_else(|| {
        format!("Unknown protocol cartoon '{protocol_id}' (run: protocol-cartoon list)")
    })
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
    fn usage_block_keeps_protocol_cartoon_routes() {
        assert!(USAGE.contains("protocol-cartoon list"));
        assert!(USAGE.contains("protocol-cartoon render-with-bindings"));
        assert!(USAGE.contains("protocol-cartoon template-export"));
    }

    #[test]
    fn protocol_cartoon_missing_subcommand_reports_legacy_message() {
        let args = argv(&["gentle_cli", "protocol-cartoon"]);
        let err = handle_protocol_cartoon_family(&args, 1, ".gentle_state.json")
            .expect_err("missing subcommand should fail");
        assert_eq!(
            err,
            "protocol-cartoon requires a subcommand: list, render-svg, render-template-svg, template-validate, render-with-bindings, template-export"
        );
    }

    #[test]
    fn protocol_cartoon_list_rejects_extra_args() {
        let args = argv(&["gentle_cli", "protocol-cartoon", "list", "extra"]);
        let err = handle_protocol_cartoon_family(&args, 1, ".gentle_state.json")
            .expect_err("extra argument should fail");
        assert_eq!(err, "protocol-cartoon list takes no additional arguments");
    }

    #[test]
    fn protocol_cartoon_list_runs_through_direct_wrapper() {
        let td = tempdir().expect("tempdir");
        let state_path = td.path().join("state.json");
        let args = argv(&["gentle_cli", "protocol-cartoon", "list"]);
        handle_protocol_cartoon_family(&args, 1, &state_path.to_string_lossy())
            .expect("protocol-cartoon list should execute");
        assert!(
            !state_path.exists(),
            "read-only protocol-cartoon list should not persist state"
        );
    }

    #[test]
    fn protocol_cartoon_render_svg_rejects_empty_protocol_id() {
        let args = argv(&[
            "gentle_cli",
            "protocol-cartoon",
            "render-svg",
            "",
            "out.svg",
        ]);
        let err = handle_protocol_cartoon_family(&args, 1, ".gentle_state.json")
            .expect_err("empty protocol id should fail");
        assert_eq!(
            err,
            "protocol-cartoon render-svg requires non-empty PROTOCOL_ID"
        );
    }

    #[test]
    fn protocol_cartoon_unknown_subcommand_reports_legacy_message() {
        let args = argv(&["gentle_cli", "protocol-cartoon", "nope"]);
        let err = handle_protocol_cartoon_family(&args, 1, ".gentle_state.json")
            .expect_err("unknown subcommand should fail");
        assert_eq!(
            err,
            "Unknown protocol-cartoon subcommand 'nope' (expected list, render-svg, render-template-svg, template-validate, render-with-bindings, template-export)"
        );
    }

    #[test]
    fn protocol_cartoon_render_svg_writes_svg_and_state() {
        let td = tempdir().expect("tempdir");
        let state_path = td.path().join("state.json");
        let output_path = td.path().join("gibson.svg");
        let output = output_path.to_string_lossy().to_string();
        let args = argv(&[
            "gentle_cli",
            "protocol-cartoon",
            "render-svg",
            "gibson.two_fragment",
            &output,
        ]);
        handle_protocol_cartoon_family(&args, 1, &state_path.to_string_lossy())
            .expect("protocol-cartoon render-svg should execute");
        assert!(output_path.exists(), "render should write SVG output");
        assert!(
            state_path.exists(),
            "mutating render route should persist state"
        );
    }
}
