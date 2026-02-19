use crate::engine::{Engine, GentleEngine, Operation, ProjectState, RenderSvgMode, Workflow};
use serde_json::{json, Value};
use std::fs;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ShellCommand {
    Help,
    Capabilities,
    StateSummary,
    LoadProject { path: String },
    SaveProject { path: String },
    RenderSvg {
        seq_id: String,
        mode: RenderSvgMode,
        output: String,
    },
    RenderLineageSvg { output: String },
    RenderPoolGelSvg {
        inputs: Vec<String>,
        output: String,
        ladders: Option<Vec<String>>,
    },
    ExportPool {
        inputs: Vec<String>,
        output: String,
        human_id: Option<String>,
    },
    Op { payload: String },
    Workflow { payload: String },
}

#[derive(Debug, Clone)]
pub struct ShellRunResult {
    pub state_changed: bool,
    pub output: Value,
}

impl ShellCommand {
    pub fn preview(&self) -> String {
        match self {
            Self::Help => "show shell command help".to_string(),
            Self::Capabilities => "inspect engine capabilities".to_string(),
            Self::StateSummary => "show sequence/container state summary".to_string(),
            Self::LoadProject { path } => format!("load project state from '{path}'"),
            Self::SaveProject { path } => format!("save current project state to '{path}'"),
            Self::RenderSvg {
                seq_id,
                mode,
                output,
            } => format!("render {mode:?} SVG for '{seq_id}' to '{output}'"),
            Self::RenderLineageSvg { output } => format!("render lineage SVG to '{output}'"),
            Self::RenderPoolGelSvg {
                inputs,
                output,
                ladders,
            } => {
                let ladders = ladders
                    .as_ref()
                    .map(|v| v.join(","))
                    .unwrap_or_else(|| "auto".to_string());
                format!(
                    "render pool gel SVG from {} input(s) to '{output}' with ladders {ladders}",
                    inputs.len()
                )
            }
            Self::ExportPool {
                inputs,
                output,
                human_id,
            } => {
                let human = human_id.clone().unwrap_or_else(|| "-".to_string());
                format!(
                    "export pool with {} input(s) to '{output}' (human_id={human})",
                    inputs.len()
                )
            }
            Self::Op { .. } => "apply one engine operation from JSON".to_string(),
            Self::Workflow { .. } => "apply engine workflow from JSON".to_string(),
        }
    }

    pub fn is_state_mutating(&self) -> bool {
        matches!(
            self,
            Self::LoadProject { .. } | Self::Op { .. } | Self::Workflow { .. }
        )
    }
}

pub fn shell_help_text() -> &'static str {
    "GENtle Shell commands:\n\
help\n\
capabilities\n\
state-summary\n\
load-project PATH\n\
save-project PATH\n\
render-svg SEQ_ID linear|circular OUTPUT.svg\n\
render-lineage-svg OUTPUT.svg\n\
render-pool-gel-svg IDS OUTPUT.svg [--ladders NAME[,NAME]]\n\
export-pool IDS OUTPUT.pool.gentle.json [HUMAN_ID]\n\
op <operation-json-or-@file>\n\
workflow <workflow-json-or-@file>\n\
IDS is comma-separated sequence IDs"
}

fn split_ids(input: &str) -> Vec<String> {
    input
        .split(',')
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .map(|s| s.to_string())
        .collect()
}

fn parse_mode(mode: &str) -> Result<RenderSvgMode, String> {
    match mode {
        "linear" => Ok(RenderSvgMode::Linear),
        "circular" => Ok(RenderSvgMode::Circular),
        other => Err(format!(
            "Unknown render mode '{other}', expected 'linear' or 'circular'"
        )),
    }
}

fn parse_json_payload(raw: &str) -> Result<String, String> {
    if let Some(path) = raw.strip_prefix('@') {
        fs::read_to_string(path).map_err(|e| format!("Could not read JSON file '{path}': {e}"))
    } else {
        Ok(raw.to_string())
    }
}

fn token_error(command: &str) -> String {
    format!("Invalid '{command}' usage. Try: help")
}

pub fn parse_shell_tokens(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.is_empty() {
        return Err("Missing shell command".to_string());
    }
    let cmd = tokens[0].as_str();
    match cmd {
        "help" | "-h" | "--help" => Ok(ShellCommand::Help),
        "capabilities" => {
            if tokens.len() == 1 {
                Ok(ShellCommand::Capabilities)
            } else {
                Err(token_error(cmd))
            }
        }
        "state-summary" => {
            if tokens.len() == 1 {
                Ok(ShellCommand::StateSummary)
            } else {
                Err(token_error(cmd))
            }
        }
        "load-project" | "import-state" => {
            if tokens.len() == 2 {
                Ok(ShellCommand::LoadProject {
                    path: tokens[1].clone(),
                })
            } else {
                Err(token_error(cmd))
            }
        }
        "save-project" | "export-state" => {
            if tokens.len() == 2 {
                Ok(ShellCommand::SaveProject {
                    path: tokens[1].clone(),
                })
            } else {
                Err(token_error(cmd))
            }
        }
        "render-svg" => {
            if tokens.len() != 4 {
                return Err(token_error(cmd));
            }
            Ok(ShellCommand::RenderSvg {
                seq_id: tokens[1].clone(),
                mode: parse_mode(&tokens[2])?,
                output: tokens[3].clone(),
            })
        }
        "render-lineage-svg" => {
            if tokens.len() == 2 {
                Ok(ShellCommand::RenderLineageSvg {
                    output: tokens[1].clone(),
                })
            } else {
                Err(token_error(cmd))
            }
        }
        "render-pool-gel-svg" => {
            if tokens.len() < 3 {
                return Err(token_error(cmd));
            }
            let inputs = split_ids(&tokens[1]);
            if inputs.is_empty() {
                return Err("render-pool-gel-svg requires at least one sequence id".to_string());
            }
            let output = tokens[2].clone();
            let mut ladders: Option<Vec<String>> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--ladders" => {
                        if idx + 1 >= tokens.len() {
                            return Err("Missing value after --ladders".to_string());
                        }
                        let parsed = split_ids(&tokens[idx + 1]);
                        if !parsed.is_empty() {
                            ladders = Some(parsed);
                        }
                        idx += 2;
                    }
                    other => {
                        return Err(format!(
                            "Unknown argument '{other}' for render-pool-gel-svg"
                        ));
                    }
                }
            }
            Ok(ShellCommand::RenderPoolGelSvg {
                inputs,
                output,
                ladders,
            })
        }
        "export-pool" => {
            if tokens.len() < 3 {
                return Err(token_error(cmd));
            }
            let inputs = split_ids(&tokens[1]);
            if inputs.is_empty() {
                return Err("export-pool requires at least one sequence id".to_string());
            }
            Ok(ShellCommand::ExportPool {
                inputs,
                output: tokens[2].clone(),
                human_id: tokens.get(3).cloned(),
            })
        }
        "op" => {
            let payload = tokens[1..].join(" ");
            if payload.trim().is_empty() {
                return Err("Missing operation JSON".to_string());
            }
            Ok(ShellCommand::Op { payload })
        }
        "workflow" => {
            let payload = tokens[1..].join(" ");
            if payload.trim().is_empty() {
                return Err("Missing workflow JSON".to_string());
            }
            Ok(ShellCommand::Workflow { payload })
        }
        other => Err(format!("Unknown shell command '{other}'. Try: help")),
    }
}

pub fn parse_shell_line(line: &str) -> Result<ShellCommand, String> {
    let tokens = split_shell_words(line)?;
    parse_shell_tokens(&tokens)
}

pub fn split_shell_words(line: &str) -> Result<Vec<String>, String> {
    #[derive(Clone, Copy, PartialEq, Eq)]
    enum Mode {
        Normal,
        SingleQuoted,
        DoubleQuoted,
    }

    let mut out = Vec::new();
    let mut current = String::new();
    let mut mode = Mode::Normal;
    let mut chars = line.chars().peekable();

    while let Some(ch) = chars.next() {
        match mode {
            Mode::Normal => match ch {
                '\'' => mode = Mode::SingleQuoted,
                '"' => mode = Mode::DoubleQuoted,
                '\\' => {
                    if let Some(next) = chars.next() {
                        current.push(next);
                    }
                }
                c if c.is_whitespace() => {
                    if !current.is_empty() {
                        out.push(current.clone());
                        current.clear();
                    }
                }
                _ => current.push(ch),
            },
            Mode::SingleQuoted => {
                if ch == '\'' {
                    mode = Mode::Normal;
                } else {
                    current.push(ch);
                }
            }
            Mode::DoubleQuoted => {
                if ch == '"' {
                    mode = Mode::Normal;
                } else if ch == '\\' {
                    if let Some(next) = chars.next() {
                        current.push(next);
                    }
                } else {
                    current.push(ch);
                }
            }
        }
    }

    if mode != Mode::Normal {
        return Err("Unterminated quoted string in shell command".to_string());
    }
    if !current.is_empty() {
        out.push(current);
    }
    if out.is_empty() {
        return Err("Empty shell command".to_string());
    }
    Ok(out)
}

pub fn execute_shell_command(
    engine: &mut GentleEngine,
    command: &ShellCommand,
) -> Result<ShellRunResult, String> {
    let result = match command {
        ShellCommand::Help => ShellRunResult {
            state_changed: false,
            output: json!({ "help": shell_help_text() }),
        },
        ShellCommand::Capabilities => ShellRunResult {
            state_changed: false,
            output: serde_json::to_value(GentleEngine::capabilities())
                .map_err(|e| format!("Could not serialize capabilities: {e}"))?,
        },
        ShellCommand::StateSummary => ShellRunResult {
            state_changed: false,
            output: serde_json::to_value(engine.summarize_state())
                .map_err(|e| format!("Could not serialize state summary: {e}"))?,
        },
        ShellCommand::LoadProject { path } => {
            let state = ProjectState::load_from_path(path).map_err(|e| e.to_string())?;
            *engine = GentleEngine::from_state(state);
            ShellRunResult {
                state_changed: true,
                output: json!({
                    "message": format!("Loaded project from '{path}'"),
                    "summary": engine.summarize_state()
                }),
            }
        }
        ShellCommand::SaveProject { path } => {
            engine.state().save_to_path(path).map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({ "message": format!("Saved project to '{path}'") }),
            }
        }
        ShellCommand::RenderSvg {
            seq_id,
            mode,
            output,
        } => {
            let op_result = engine
                .apply(Operation::RenderSequenceSvg {
                    seq_id: seq_id.clone(),
                    mode: mode.clone(),
                    path: output.clone(),
                })
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::RenderLineageSvg { output } => {
            let op_result = engine
                .apply(Operation::RenderLineageSvg {
                    path: output.clone(),
                })
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::RenderPoolGelSvg {
            inputs,
            output,
            ladders,
        } => {
            let op_result = engine
                .apply(Operation::RenderPoolGelSvg {
                    inputs: inputs.clone(),
                    path: output.clone(),
                    ladders: ladders.clone(),
                })
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::ExportPool {
            inputs,
            output,
            human_id,
        } => {
            let op_result = engine
                .apply(Operation::ExportPool {
                    inputs: inputs.clone(),
                    path: output.clone(),
                    pool_id: Some("pool_export".to_string()),
                    human_id: human_id.clone(),
                })
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::Op { payload } => {
            let json_text = parse_json_payload(payload)?;
            let op: Operation = serde_json::from_str(&json_text)
                .map_err(|e| format!("Invalid operation JSON: {e}"))?;
            let op_result = engine.apply(op).map_err(|e| e.to_string())?;
            let state_changed =
                !op_result.created_seq_ids.is_empty() || !op_result.changed_seq_ids.is_empty();
            ShellRunResult {
                state_changed,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::Workflow { payload } => {
            let json_text = parse_json_payload(payload)?;
            let workflow: Workflow = serde_json::from_str(&json_text)
                .map_err(|e| format!("Invalid workflow JSON: {e}"))?;
            let results = engine.apply_workflow(workflow).map_err(|e| e.to_string())?;
            let state_changed = results
                .iter()
                .any(|r| !r.created_seq_ids.is_empty() || !r.changed_seq_ids.is_empty());
            ShellRunResult {
                state_changed,
                output: json!({ "results": results }),
            }
        }
    };
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_workflow_payload_keeps_whitespace() {
        let cmd = parse_shell_line("workflow { \"run_id\": \"x\", \"ops\": [] }")
            .expect("workflow command parse");
        match cmd {
            ShellCommand::Workflow { payload } => {
                assert!(payload.contains("\"run_id\""));
                assert!(payload.contains("\"ops\""));
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_render_pool_gel_with_ladders() {
        let cmd = parse_shell_line("render-pool-gel-svg a,b out.svg --ladders 1kb,100bp")
            .expect("parse command");
        match cmd {
            ShellCommand::RenderPoolGelSvg {
                inputs,
                output,
                ladders,
            } => {
                assert_eq!(inputs, vec!["a".to_string(), "b".to_string()]);
                assert_eq!(output, "out.svg".to_string());
                assert_eq!(
                    ladders,
                    Some(vec!["1kb".to_string(), "100bp".to_string()])
                );
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn execute_state_summary_returns_json() {
        let mut engine = GentleEngine::from_state(ProjectState::default());
        let out = execute_shell_command(&mut engine, &ShellCommand::StateSummary)
            .expect("execute state summary");
        assert!(!out.state_changed);
        assert!(out.output.get("sequence_count").is_some());
    }
}
