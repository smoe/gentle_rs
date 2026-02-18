use gentle::render_export::{export_circular_svg, export_linear_svg};
use gentle::{
    about,
    engine::{Engine, GentleEngine, Operation, ProjectState, Workflow},
};
use serde::Serialize;
use std::{env, fs};

const DEFAULT_STATE_PATH: &str = ".gentle_state.json";

#[derive(Serialize)]
struct SequenceSummary {
    id: String,
    name: Option<String>,
    length: usize,
    circular: bool,
}

#[derive(Serialize)]
struct StateSummary {
    sequence_count: usize,
    sequences: Vec<SequenceSummary>,
    display: gentle::engine::DisplaySettings,
}

fn usage() {
    eprintln!(
        "Usage:\n  \
  gentle_cli --version\n  \
  gentle_cli [--state PATH] capabilities\n  \
  gentle_cli [--state PATH] op '<operation-json>'\n  \
  gentle_cli [--state PATH] workflow '<workflow-json>'\n  \
  gentle_cli [--state PATH] state-summary\n  \
  gentle_cli [--state PATH] export-state PATH\n  \
  gentle_cli [--state PATH] import-state PATH\n  \
  gentle_cli [--state PATH] save-project PATH\n  \
  gentle_cli [--state PATH] load-project PATH\n  \
  gentle_cli [--state PATH] render-svg SEQ_ID linear|circular OUTPUT.svg\n\n  \
  Tip: pass @file.json instead of inline JSON"
    );
}

fn load_json_arg(value: &str) -> Result<String, String> {
    if let Some(path) = value.strip_prefix('@') {
        fs::read_to_string(path).map_err(|e| format!("Could not read JSON file '{path}': {e}"))
    } else {
        Ok(value.to_string())
    }
}

fn load_state(path: &str) -> Result<ProjectState, String> {
    if std::path::Path::new(path).exists() {
        ProjectState::load_from_path(path).map_err(|e| e.to_string())
    } else {
        Ok(ProjectState::default())
    }
}

fn print_json<T: Serialize>(value: &T) -> Result<(), String> {
    let text = serde_json::to_string_pretty(value)
        .map_err(|e| format!("Could not serialize JSON output: {e}"))?;
    println!("{text}");
    Ok(())
}

fn parse_global_state_arg(args: &[String]) -> (String, usize) {
    if args.len() >= 3 && args[1] == "--state" {
        return (args[2].clone(), 3);
    }
    (DEFAULT_STATE_PATH.to_string(), 1)
}

fn summarize_state(engine: &GentleEngine) -> StateSummary {
    let mut sequences: Vec<SequenceSummary> = engine
        .state()
        .sequences
        .iter()
        .map(|(id, dna)| SequenceSummary {
            id: id.to_string(),
            name: dna.name().clone(),
            length: dna.len(),
            circular: dna.is_circular(),
        })
        .collect();
    sequences.sort_by(|a, b| a.id.cmp(&b.id));

    StateSummary {
        sequence_count: sequences.len(),
        sequences,
        display: engine.state().display.clone(),
    }
}

fn main() {
    if let Err(e) = run() {
        eprintln!("{e}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), String> {
    let args: Vec<String> = env::args().collect();
    if args.len() <= 1 {
        usage();
        return Err("Missing command".to_string());
    }
    if args.iter().any(|a| a == "--version" || a == "-V") {
        println!("{}", about::version_cli_text());
        return Ok(());
    }

    let (state_path, cmd_idx) = parse_global_state_arg(&args);
    if args.len() <= cmd_idx {
        usage();
        return Err("Missing command".to_string());
    }

    let command = &args[cmd_idx];

    match command.as_str() {
        "capabilities" => {
            print_json(&GentleEngine::capabilities())?;
            Ok(())
        }
        "import-state" | "load-project" => {
            if args.len() <= cmd_idx + 1 {
                usage();
                return Err(format!("Missing path for {command}"));
            }
            let source = &args[cmd_idx + 1];
            let state = ProjectState::load_from_path(source).map_err(|e| e.to_string())?;
            state.save_to_path(&state_path).map_err(|e| e.to_string())?;
            println!("Loaded project from '{source}' into '{state_path}'");
            Ok(())
        }
        "export-state" | "save-project" => {
            if args.len() <= cmd_idx + 1 {
                usage();
                return Err(format!("Missing path for {command}"));
            }
            let target = &args[cmd_idx + 1];
            let state = load_state(&state_path)?;
            state.save_to_path(target).map_err(|e| e.to_string())?;
            println!("Saved project from '{state_path}' to '{target}'");
            Ok(())
        }
        "state-summary" => {
            let state = load_state(&state_path)?;
            let engine = GentleEngine::from_state(state);
            print_json(&summarize_state(&engine))
        }
        "render-svg" => {
            if args.len() <= cmd_idx + 3 {
                usage();
                return Err("render-svg requires: SEQ_ID linear|circular OUTPUT.svg".to_string());
            }
            let seq_id = &args[cmd_idx + 1];
            let mode = &args[cmd_idx + 2];
            let output = &args[cmd_idx + 3];

            let state = load_state(&state_path)?;
            let dna = state
                .sequences
                .get(seq_id)
                .ok_or_else(|| format!("Sequence '{seq_id}' not found in state '{state_path}'"))?;
            let svg = match mode.as_str() {
                "linear" => export_linear_svg(dna, &state.display),
                "circular" => export_circular_svg(dna, &state.display),
                _ => {
                    return Err(format!(
                        "Unknown render mode '{mode}', expected 'linear' or 'circular'"
                    ))
                }
            };
            fs::write(output, svg)
                .map_err(|e| format!("Could not write SVG output '{output}': {e}"))?;
            println!("Wrote {mode} SVG for '{seq_id}' to '{output}'");
            Ok(())
        }
        "op" => {
            if args.len() <= cmd_idx + 1 {
                usage();
                return Err("Missing operation JSON".to_string());
            }
            let json = load_json_arg(&args[cmd_idx + 1])?;
            let op: Operation =
                serde_json::from_str(&json).map_err(|e| format!("Invalid operation JSON: {e}"))?;

            let mut engine = GentleEngine::from_state(load_state(&state_path)?);
            let result = engine.apply(op).map_err(|e| e.to_string())?;
            engine
                .state()
                .save_to_path(&state_path)
                .map_err(|e| e.to_string())?;
            print_json(&result)
        }
        "workflow" => {
            if args.len() <= cmd_idx + 1 {
                usage();
                return Err("Missing workflow JSON".to_string());
            }
            let json = load_json_arg(&args[cmd_idx + 1])?;
            let workflow: Workflow =
                serde_json::from_str(&json).map_err(|e| format!("Invalid workflow JSON: {e}"))?;

            let mut engine = GentleEngine::from_state(load_state(&state_path)?);
            let results = engine.apply_workflow(workflow).map_err(|e| e.to_string())?;
            engine
                .state()
                .save_to_path(&state_path)
                .map_err(|e| e.to_string())?;
            print_json(&results)
        }
        _ => {
            usage();
            Err(format!("Unknown command '{command}'"))
        }
    }
}
