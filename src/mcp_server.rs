//! Minimal MCP stdio server adapter.
//!
//! This module implements a small MCP surface on top of existing deterministic
//! engine/shell contracts.

use crate::{
    about,
    engine::{Engine, GentleEngine, Operation, ProjectState, Workflow},
    engine_shell::{ShellExecutionOptions, execute_shell_command_with_options, parse_shell_tokens},
    shell_docs::{
        shell_help_json, shell_help_markdown, shell_help_text, shell_topic_help_json,
        shell_topic_help_markdown, shell_topic_help_text,
    },
};
use serde::Deserialize;
use serde_json::{Map, Value, json};
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

const MCP_PROTOCOL_VERSION: &str = "2025-06-18";
const SERVER_NAME: &str = "gentle_mcp";
const SERVER_TITLE: &str = "GENtle MCP";

pub const DEFAULT_MCP_STATE_PATH: &str = ".gentle_state.json";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DispatchOutcome {
    NoResponse,
    Response,
    Exit,
}

#[derive(Debug, Clone, Deserialize)]
struct ToolCallParams {
    name: String,
    #[serde(default)]
    arguments: Value,
}

pub fn run_stdio_server(state_path: &str) -> Result<(), String> {
    let stdin = io::stdin();
    let stdout = io::stdout();
    let mut reader = BufReader::new(stdin.lock());
    let mut writer = BufWriter::new(stdout.lock());
    run_server_loop(state_path, &mut reader, &mut writer)
}

fn run_server_loop<R: BufRead, W: Write>(
    default_state_path: &str,
    reader: &mut R,
    writer: &mut W,
) -> Result<(), String> {
    loop {
        let Some(message) = read_framed_json(reader)? else {
            return Ok(());
        };
        match handle_message(default_state_path, &message, writer)? {
            DispatchOutcome::NoResponse => {}
            DispatchOutcome::Response => {}
            DispatchOutcome::Exit => return Ok(()),
        }
    }
}

fn load_state(path: &str) -> Result<ProjectState, String> {
    if Path::new(path).exists() {
        ProjectState::load_from_path(path).map_err(|e| e.to_string())
    } else {
        Ok(ProjectState::default())
    }
}

fn read_framed_json<R: BufRead>(reader: &mut R) -> Result<Option<Value>, String> {
    let mut content_length: Option<usize> = None;

    loop {
        let mut line = String::new();
        let bytes_read = reader
            .read_line(&mut line)
            .map_err(|e| format!("Could not read MCP header line: {e}"))?;
        if bytes_read == 0 {
            return if content_length.is_some() {
                Err("Unexpected EOF while reading MCP headers".to_string())
            } else {
                Ok(None)
            };
        }
        let line_trimmed = line.trim_end_matches(['\r', '\n']);
        if line_trimmed.is_empty() {
            if content_length.is_some() {
                break;
            }
            continue;
        }
        if let Some(value) = line_trimmed.strip_prefix("Content-Length:") {
            let len = value
                .trim()
                .parse::<usize>()
                .map_err(|e| format!("Invalid Content-Length header '{line_trimmed}': {e}"))?;
            content_length = Some(len);
        }
    }

    let len = content_length.ok_or_else(|| "Missing Content-Length header".to_string())?;
    let mut body = vec![0u8; len];
    reader
        .read_exact(&mut body)
        .map_err(|e| format!("Could not read MCP JSON payload body: {e}"))?;
    serde_json::from_slice::<Value>(&body)
        .map(Some)
        .map_err(|e| format!("Could not parse MCP JSON payload: {e}"))
}

fn write_framed_json<W: Write>(writer: &mut W, payload: &Value) -> Result<(), String> {
    let body = serde_json::to_vec(payload)
        .map_err(|e| format!("Could not serialize MCP response JSON: {e}"))?;
    writer
        .write_all(format!("Content-Length: {}\r\n\r\n", body.len()).as_bytes())
        .map_err(|e| format!("Could not write MCP response header: {e}"))?;
    writer
        .write_all(&body)
        .map_err(|e| format!("Could not write MCP response body: {e}"))?;
    writer
        .flush()
        .map_err(|e| format!("Could not flush MCP response stream: {e}"))?;
    Ok(())
}

fn tool_list() -> Value {
    json!([
        {
            "name": "capabilities",
            "title": "Capabilities",
            "description": "Return shared GENtle engine capabilities.",
            "inputSchema": {
                "type": "object",
                "properties": {},
                "additionalProperties": false
            }
        },
        {
            "name": "state_summary",
            "title": "State Summary",
            "description": "Return a deterministic summary of the current project state.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "state_path": {
                        "type": "string",
                        "description": "Optional project state path. Defaults to server startup state path."
                    }
                },
                "additionalProperties": false
            }
        },
        {
            "name": "op",
            "title": "Apply Operation",
            "description": "Apply one operation via shared engine contract and persist state (requires confirm=true).",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "confirm": {
                        "type": "boolean",
                        "description": "Must be true for mutating MCP tool execution."
                    },
                    "state_path": {
                        "type": "string",
                        "description": "Optional project state path. Defaults to server startup state path."
                    },
                    "operation": {
                        "type": "object",
                        "description": "Operation payload in engine Operation enum JSON shape."
                    },
                    "op": {
                        "type": "object",
                        "description": "Alias for operation."
                    }
                },
                "required": ["confirm"],
                "additionalProperties": false
            }
        },
        {
            "name": "workflow",
            "title": "Apply Workflow",
            "description": "Apply a workflow via shared engine contract and persist state (requires confirm=true).",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "confirm": {
                        "type": "boolean",
                        "description": "Must be true for mutating MCP tool execution."
                    },
                    "state_path": {
                        "type": "string",
                        "description": "Optional project state path. Defaults to server startup state path."
                    },
                    "workflow": {
                        "type": "object",
                        "description": "Workflow payload in engine Workflow JSON shape."
                    },
                    "wf": {
                        "type": "object",
                        "description": "Alias for workflow."
                    }
                },
                "required": ["confirm"],
                "additionalProperties": false
            }
        },
        {
            "name": "help",
            "title": "Help",
            "description": "Return shell command reference content from docs/glossary.json.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "format": {
                        "type": "string",
                        "enum": ["text", "json", "markdown"]
                    },
                    "interface": {
                        "type": "string",
                        "description": "Optional interface filter (all|cli-direct|cli-shell|gui-shell|js|lua|mcp)."
                    },
                    "topic": {
                        "description": "Optional help topic (string path or array of path tokens).",
                        "oneOf": [
                            { "type": "string" },
                            { "type": "array", "items": { "type": "string" } }
                        ]
                    }
                },
                "additionalProperties": false
            }
        },
        {
            "name": "ui_intents",
            "title": "UI Intents Catalog",
            "description": "Return deterministic UI-intent target/command catalog (shared `ui intents` contract).",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "state_path": {
                        "type": "string",
                        "description": "Optional project state path. Defaults to server startup state path."
                    }
                },
                "additionalProperties": false
            }
        },
        {
            "name": "ui_intent",
            "title": "UI Intent",
            "description": "Resolve/record one UI intent through shared `ui open|focus` parser/executor path.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "state_path": {
                        "type": "string",
                        "description": "Optional project state path. Defaults to server startup state path."
                    },
                    "action": {
                        "type": "string",
                        "enum": ["open", "focus"]
                    },
                    "target": {
                        "type": "string",
                        "enum": [
                            "prepared-references",
                            "prepare-reference-genome",
                            "retrieve-genome-sequence",
                            "blast-genome-sequence",
                            "import-genome-track",
                            "agent-assistant",
                            "prepare-helper-genome",
                            "retrieve-helper-sequence",
                            "blast-helper-sequence"
                        ]
                    },
                    "genome_id": {
                        "type": "string"
                    },
                    "helpers": {
                        "type": "boolean"
                    },
                    "catalog_path": {
                        "type": "string"
                    },
                    "cache_dir": {
                        "type": "string"
                    },
                    "filter": {
                        "type": "string"
                    },
                    "species": {
                        "type": "string"
                    },
                    "latest": {
                        "type": "boolean"
                    }
                },
                "required": ["action", "target"],
                "additionalProperties": false
            }
        },
        {
            "name": "ui_prepared_genomes",
            "title": "UI Prepared Genomes Query",
            "description": "Return deterministic prepared-genome rows via shared `ui prepared-genomes` contract.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "state_path": {
                        "type": "string",
                        "description": "Optional project state path. Defaults to server startup state path."
                    },
                    "helpers": {
                        "type": "boolean"
                    },
                    "catalog_path": {
                        "type": "string"
                    },
                    "cache_dir": {
                        "type": "string"
                    },
                    "filter": {
                        "type": "string"
                    },
                    "species": {
                        "type": "string"
                    },
                    "latest": {
                        "type": "boolean"
                    }
                },
                "additionalProperties": false
            }
        },
        {
            "name": "ui_latest_prepared",
            "title": "UI Latest Prepared",
            "description": "Resolve latest prepared genome for a species via shared `ui latest-prepared` contract.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "state_path": {
                        "type": "string",
                        "description": "Optional project state path. Defaults to server startup state path."
                    },
                    "species": {
                        "type": "string"
                    },
                    "helpers": {
                        "type": "boolean"
                    },
                    "catalog_path": {
                        "type": "string"
                    },
                    "cache_dir": {
                        "type": "string"
                    }
                },
                "required": ["species"],
                "additionalProperties": false
            }
        }
    ])
}

fn jsonrpc_response(id: Value, result: Value) -> Value {
    json!({
        "jsonrpc": "2.0",
        "id": id,
        "result": result
    })
}

fn jsonrpc_error(id: Option<Value>, code: i64, message: &str, data: Option<Value>) -> Value {
    let mut error = json!({
        "code": code,
        "message": message
    });
    if let Some(data) = data {
        error["data"] = data;
    }
    json!({
        "jsonrpc": "2.0",
        "id": id.unwrap_or(Value::Null),
        "error": error
    })
}

fn tool_result_text(text: String, format: &'static str, is_error: bool) -> Value {
    json!({
        "content": [
            {
                "type": "text",
                "text": text
            }
        ],
        "structuredContent": {
            "format": format,
            "text": text
        },
        "isError": is_error
    })
}

fn tool_result_json(value: Value, is_error: bool) -> Value {
    let text = serde_json::to_string_pretty(&value).unwrap_or_else(|_| value.to_string());
    json!({
        "content": [
            {
                "type": "text",
                "text": text
            }
        ],
        "structuredContent": value,
        "isError": is_error
    })
}

fn parse_topic(raw: Option<&Value>) -> Result<Option<Vec<String>>, String> {
    let Some(raw) = raw else {
        return Ok(None);
    };
    match raw {
        Value::Null => Ok(None),
        Value::String(value) => {
            let topic = value
                .split_whitespace()
                .map(str::trim)
                .filter(|v| !v.is_empty())
                .map(ToString::to_string)
                .collect::<Vec<_>>();
            if topic.is_empty() {
                Ok(None)
            } else {
                Ok(Some(topic))
            }
        }
        Value::Array(values) => {
            let mut topic = vec![];
            for value in values {
                let token = value.as_str().ok_or_else(|| {
                    "help.topic array must contain only string values".to_string()
                })?;
                let token = token.trim();
                if !token.is_empty() {
                    topic.push(token.to_string());
                }
            }
            if topic.is_empty() {
                Ok(None)
            } else {
                Ok(Some(topic))
            }
        }
        _ => Err("help.topic must be a string or array of strings".to_string()),
    }
}

fn help_tool_result(arguments: &Value) -> Value {
    let args = arguments.as_object().cloned().unwrap_or_default();
    let format = args
        .get("format")
        .and_then(Value::as_str)
        .unwrap_or("text")
        .trim()
        .to_ascii_lowercase();
    let interface = args.get("interface").and_then(Value::as_str);
    let topic = match parse_topic(args.get("topic")) {
        Ok(topic) => topic,
        Err(err) => return tool_result_text(err, "text", true),
    };

    match (format.as_str(), topic) {
        ("json", Some(topic)) => match shell_topic_help_json(&topic, interface) {
            Ok(value) => tool_result_json(value, false),
            Err(err) => tool_result_text(err, "text", true),
        },
        ("json", None) => match shell_help_json(interface) {
            Ok(value) => tool_result_json(value, false),
            Err(err) => tool_result_text(err, "text", true),
        },
        ("markdown", Some(topic)) => match shell_topic_help_markdown(&topic, interface) {
            Ok(text) => tool_result_text(text, "markdown", false),
            Err(err) => tool_result_text(err, "text", true),
        },
        ("markdown", None) => match shell_help_markdown(interface) {
            Ok(text) => tool_result_text(text, "markdown", false),
            Err(err) => tool_result_text(err, "text", true),
        },
        ("text", Some(topic)) => match shell_topic_help_text(&topic, interface) {
            Ok(text) => tool_result_text(text, "text", false),
            Err(err) => tool_result_text(err, "text", true),
        },
        ("text", None) => match shell_help_text(interface) {
            Ok(text) => tool_result_text(text, "text", false),
            Err(err) => tool_result_text(err, "text", true),
        },
        (other, _) => tool_result_text(
            format!("Unsupported help.format '{other}' (expected text|json|markdown)"),
            "text",
            true,
        ),
    }
}

fn state_path_from_args(default_state_path: &str, args: &Map<String, Value>) -> String {
    args.get("state_path")
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|v| !v.is_empty())
        .unwrap_or(default_state_path)
        .to_string()
}

fn require_confirm_true(args: &Map<String, Value>, tool_name: &str) -> Result<(), String> {
    match args.get("confirm") {
        Some(Value::Bool(true)) => Ok(()),
        Some(Value::Bool(false)) | None => Err(format!(
            "Refusing MCP tool '{tool_name}' without explicit confirm=true"
        )),
        Some(_) => Err(format!(
            "MCP tool '{tool_name}' requires boolean confirm=true"
        )),
    }
}

fn optional_string_arg(args: &Map<String, Value>, key: &str) -> Result<Option<String>, String> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::String(value)) => {
            let value = value.trim();
            if value.is_empty() {
                Ok(None)
            } else {
                Ok(Some(value.to_string()))
            }
        }
        Some(_) => Err(format!("MCP argument '{key}' must be a string")),
    }
}

fn required_string_arg(args: &Map<String, Value>, key: &str) -> Result<String, String> {
    match optional_string_arg(args, key)? {
        Some(value) => Ok(value),
        None => Err(format!(
            "MCP argument '{key}' is required and must be non-empty"
        )),
    }
}

fn optional_bool_arg(args: &Map<String, Value>, key: &str) -> Result<Option<bool>, String> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Bool(value)) => Ok(Some(*value)),
        Some(_) => Err(format!("MCP argument '{key}' must be a boolean")),
    }
}

fn append_string_flag(tokens: &mut Vec<String>, flag: &str, value: Option<String>) {
    if let Some(value) = value {
        tokens.push(flag.to_string());
        tokens.push(value);
    }
}

fn run_non_mutating_shell_tool(
    default_state_path: &str,
    args: &Map<String, Value>,
    tokens: Vec<String>,
    tool_name: &str,
) -> Result<Value, String> {
    let state_path = state_path_from_args(default_state_path, args);
    let state = load_state(&state_path)
        .map_err(|err| format!("Could not load state from '{state_path}': {err}"))?;
    let mut engine = GentleEngine::from_state(state);
    let command = parse_shell_tokens(&tokens).map_err(|err| {
        format!("Could not parse shared shell command for MCP tool '{tool_name}': {err}")
    })?;
    let run = execute_shell_command_with_options(
        &mut engine,
        &command,
        &ShellExecutionOptions::default(),
    )
    .map_err(|err| {
        format!("Could not execute shared shell command for MCP tool '{tool_name}': {err}")
    })?;
    if run.state_changed {
        return Err(format!(
            "MCP tool '{tool_name}' expected non-mutating shared shell command but state_changed=true"
        ));
    }
    Ok(run.output)
}

fn ui_intents_tool_result(default_state_path: &str, arguments: &Value) -> Value {
    let args = arguments.as_object().cloned().unwrap_or_default();
    match run_non_mutating_shell_tool(
        default_state_path,
        &args,
        vec!["ui".to_string(), "intents".to_string()],
        "ui_intents",
    ) {
        Ok(output) => tool_result_json(output, false),
        Err(err) => tool_result_text(err, "text", true),
    }
}

fn ui_prepared_genomes_tool_result(default_state_path: &str, arguments: &Value) -> Value {
    let args = arguments.as_object().cloned().unwrap_or_default();
    let helper_mode = match optional_bool_arg(&args, "helpers") {
        Ok(value) => value.unwrap_or(false),
        Err(err) => return tool_result_text(err, "text", true),
    };
    let latest = match optional_bool_arg(&args, "latest") {
        Ok(value) => value.unwrap_or(false),
        Err(err) => return tool_result_text(err, "text", true),
    };
    let catalog_path = match optional_string_arg(&args, "catalog_path") {
        Ok(value) => value,
        Err(err) => return tool_result_text(err, "text", true),
    };
    let cache_dir = match optional_string_arg(&args, "cache_dir") {
        Ok(value) => value,
        Err(err) => return tool_result_text(err, "text", true),
    };
    let filter = match optional_string_arg(&args, "filter") {
        Ok(value) => value,
        Err(err) => return tool_result_text(err, "text", true),
    };
    let species = match optional_string_arg(&args, "species") {
        Ok(value) => value,
        Err(err) => return tool_result_text(err, "text", true),
    };
    let mut tokens = vec!["ui".to_string(), "prepared-genomes".to_string()];
    if helper_mode {
        tokens.push("--helpers".to_string());
    }
    append_string_flag(&mut tokens, "--catalog", catalog_path);
    append_string_flag(&mut tokens, "--cache-dir", cache_dir);
    append_string_flag(&mut tokens, "--filter", filter);
    append_string_flag(&mut tokens, "--species", species);
    if latest {
        tokens.push("--latest".to_string());
    }
    match run_non_mutating_shell_tool(default_state_path, &args, tokens, "ui_prepared_genomes") {
        Ok(output) => tool_result_json(output, false),
        Err(err) => tool_result_text(err, "text", true),
    }
}

fn ui_latest_prepared_tool_result(default_state_path: &str, arguments: &Value) -> Value {
    let args = arguments.as_object().cloned().unwrap_or_default();
    let species = match required_string_arg(&args, "species") {
        Ok(value) => value,
        Err(err) => return tool_result_text(err, "text", true),
    };
    let helper_mode = match optional_bool_arg(&args, "helpers") {
        Ok(value) => value.unwrap_or(false),
        Err(err) => return tool_result_text(err, "text", true),
    };
    let catalog_path = match optional_string_arg(&args, "catalog_path") {
        Ok(value) => value,
        Err(err) => return tool_result_text(err, "text", true),
    };
    let cache_dir = match optional_string_arg(&args, "cache_dir") {
        Ok(value) => value,
        Err(err) => return tool_result_text(err, "text", true),
    };
    let mut tokens = vec!["ui".to_string(), "latest-prepared".to_string(), species];
    if helper_mode {
        tokens.push("--helpers".to_string());
    }
    append_string_flag(&mut tokens, "--catalog", catalog_path);
    append_string_flag(&mut tokens, "--cache-dir", cache_dir);
    match run_non_mutating_shell_tool(default_state_path, &args, tokens, "ui_latest_prepared") {
        Ok(output) => tool_result_json(output, false),
        Err(err) => tool_result_text(err, "text", true),
    }
}

fn ui_intent_tool_result(default_state_path: &str, arguments: &Value) -> Value {
    let args = arguments.as_object().cloned().unwrap_or_default();
    let action = match required_string_arg(&args, "action") {
        Ok(value) => value.to_ascii_lowercase(),
        Err(err) => return tool_result_text(err, "text", true),
    };
    if action != "open" && action != "focus" {
        return tool_result_text(
            format!("MCP argument 'action' must be 'open' or 'focus' (got '{action}')"),
            "text",
            true,
        );
    }
    let target = match required_string_arg(&args, "target") {
        Ok(value) => value,
        Err(err) => return tool_result_text(err, "text", true),
    };
    let helper_mode = match optional_bool_arg(&args, "helpers") {
        Ok(value) => value.unwrap_or(false),
        Err(err) => return tool_result_text(err, "text", true),
    };
    let latest = match optional_bool_arg(&args, "latest") {
        Ok(value) => value.unwrap_or(false),
        Err(err) => return tool_result_text(err, "text", true),
    };
    let genome_id = match optional_string_arg(&args, "genome_id") {
        Ok(value) => value,
        Err(err) => return tool_result_text(err, "text", true),
    };
    let catalog_path = match optional_string_arg(&args, "catalog_path") {
        Ok(value) => value,
        Err(err) => return tool_result_text(err, "text", true),
    };
    let cache_dir = match optional_string_arg(&args, "cache_dir") {
        Ok(value) => value,
        Err(err) => return tool_result_text(err, "text", true),
    };
    let filter = match optional_string_arg(&args, "filter") {
        Ok(value) => value,
        Err(err) => return tool_result_text(err, "text", true),
    };
    let species = match optional_string_arg(&args, "species") {
        Ok(value) => value,
        Err(err) => return tool_result_text(err, "text", true),
    };

    let mut tokens = vec!["ui".to_string(), action, target];
    append_string_flag(&mut tokens, "--genome-id", genome_id);
    if helper_mode {
        tokens.push("--helpers".to_string());
    }
    append_string_flag(&mut tokens, "--catalog", catalog_path);
    append_string_flag(&mut tokens, "--cache-dir", cache_dir);
    append_string_flag(&mut tokens, "--filter", filter);
    append_string_flag(&mut tokens, "--species", species);
    if latest {
        tokens.push("--latest".to_string());
    }
    match run_non_mutating_shell_tool(default_state_path, &args, tokens, "ui_intent") {
        Ok(output) => tool_result_json(output, false),
        Err(err) => tool_result_text(err, "text", true),
    }
}

fn op_tool_result(default_state_path: &str, arguments: &Value) -> Value {
    let args = arguments.as_object().cloned().unwrap_or_default();
    if let Err(err) = require_confirm_true(&args, "op") {
        return tool_result_text(err, "text", true);
    }
    let state_path = state_path_from_args(default_state_path, &args);
    let Some(raw_operation) = args.get("operation").or_else(|| args.get("op")) else {
        return tool_result_text(
            "op requires an 'operation' argument (or alias 'op')".to_string(),
            "text",
            true,
        );
    };
    let operation = match serde_json::from_value::<Operation>(raw_operation.clone()) {
        Ok(op) => op,
        Err(err) => {
            return tool_result_text(
                format!("Could not parse operation payload: {err}"),
                "text",
                true,
            );
        }
    };
    let state = match load_state(&state_path) {
        Ok(state) => state,
        Err(err) => {
            return tool_result_text(
                format!("Could not load state from '{state_path}': {err}"),
                "text",
                true,
            );
        }
    };
    let mut engine = GentleEngine::from_state(state);
    let result = match engine.apply(operation) {
        Ok(result) => result,
        Err(err) => {
            return tool_result_json(json!({ "state_path": state_path, "error": err }), true);
        }
    };
    if let Err(err) = engine.state().save_to_path(&state_path) {
        return tool_result_json(json!({ "state_path": state_path, "error": err }), true);
    }
    tool_result_json(
        json!({
            "state_path": state_path,
            "result": result,
            "state_summary": engine.summarize_state(),
        }),
        false,
    )
}

fn workflow_tool_result(default_state_path: &str, arguments: &Value) -> Value {
    let args = arguments.as_object().cloned().unwrap_or_default();
    if let Err(err) = require_confirm_true(&args, "workflow") {
        return tool_result_text(err, "text", true);
    }
    let state_path = state_path_from_args(default_state_path, &args);
    let Some(raw_workflow) = args.get("workflow").or_else(|| args.get("wf")) else {
        return tool_result_text(
            "workflow requires a 'workflow' argument (or alias 'wf')".to_string(),
            "text",
            true,
        );
    };
    let workflow = match serde_json::from_value::<Workflow>(raw_workflow.clone()) {
        Ok(wf) => wf,
        Err(err) => {
            return tool_result_text(
                format!("Could not parse workflow payload: {err}"),
                "text",
                true,
            );
        }
    };
    let state = match load_state(&state_path) {
        Ok(state) => state,
        Err(err) => {
            return tool_result_text(
                format!("Could not load state from '{state_path}': {err}"),
                "text",
                true,
            );
        }
    };
    let mut engine = GentleEngine::from_state(state);
    let results = match engine.apply_workflow(workflow) {
        Ok(results) => results,
        Err(err) => {
            return tool_result_json(json!({ "state_path": state_path, "error": err }), true);
        }
    };
    if let Err(err) = engine.state().save_to_path(&state_path) {
        return tool_result_json(json!({ "state_path": state_path, "error": err }), true);
    }
    tool_result_json(
        json!({
            "state_path": state_path,
            "result_count": results.len(),
            "results": results,
            "state_summary": engine.summarize_state(),
        }),
        false,
    )
}

fn tool_call_result(default_state_path: &str, params: ToolCallParams) -> Value {
    match params.name.trim() {
        "capabilities" => tool_result_json(json!(GentleEngine::capabilities()), false),
        "state_summary" => {
            let args = params.arguments.as_object().cloned().unwrap_or_default();
            let state_path = state_path_from_args(default_state_path, &args);
            match load_state(&state_path) {
                Ok(state) => {
                    let engine = GentleEngine::from_state(state);
                    tool_result_json(json!(engine.summarize_state()), false)
                }
                Err(err) => tool_result_text(
                    format!("Could not load state from '{state_path}': {err}"),
                    "text",
                    true,
                ),
            }
        }
        "op" => op_tool_result(default_state_path, &params.arguments),
        "workflow" => workflow_tool_result(default_state_path, &params.arguments),
        "help" => help_tool_result(&params.arguments),
        "ui_intents" => ui_intents_tool_result(default_state_path, &params.arguments),
        "ui_intent" => ui_intent_tool_result(default_state_path, &params.arguments),
        "ui_prepared_genomes" => {
            ui_prepared_genomes_tool_result(default_state_path, &params.arguments)
        }
        "ui_latest_prepared" => {
            ui_latest_prepared_tool_result(default_state_path, &params.arguments)
        }
        other => tool_result_text(format!("Unknown MCP tool '{other}'"), "text", true),
    }
}

fn write_response<W: Write>(writer: &mut W, value: Value) -> Result<DispatchOutcome, String> {
    write_framed_json(writer, &value)?;
    Ok(DispatchOutcome::Response)
}

fn handle_message<W: Write>(
    default_state_path: &str,
    message: &Value,
    writer: &mut W,
) -> Result<DispatchOutcome, String> {
    let Some(obj) = message.as_object() else {
        return write_response(
            writer,
            jsonrpc_error(None, -32600, "Invalid Request: expected JSON object", None),
        );
    };
    let id = obj.get("id").cloned();
    let Some(method) = obj.get("method").and_then(Value::as_str) else {
        return write_response(
            writer,
            jsonrpc_error(
                id,
                -32600,
                "Invalid Request: missing method field",
                Some(message.clone()),
            ),
        );
    };

    match method {
        "initialize" => {
            let Some(id) = id else {
                return write_response(
                    writer,
                    jsonrpc_error(
                        None,
                        -32600,
                        "Invalid Request: initialize requires id",
                        None,
                    ),
                );
            };
            let result = json!({
                "protocolVersion": MCP_PROTOCOL_VERSION,
                "capabilities": {
                    "tools": {
                        "listChanged": false
                    }
                },
                "serverInfo": {
                    "name": SERVER_NAME,
                    "title": SERVER_TITLE,
                    "version": about::GENTLE_DISPLAY_VERSION
                }
            });
            write_response(writer, jsonrpc_response(id, result))
        }
        "notifications/initialized" => Ok(DispatchOutcome::NoResponse),
        "ping" => {
            if let Some(id) = id {
                write_response(writer, jsonrpc_response(id, json!({})))
            } else {
                Ok(DispatchOutcome::NoResponse)
            }
        }
        "tools/list" => {
            let Some(id) = id else {
                return Ok(DispatchOutcome::NoResponse);
            };
            write_response(
                writer,
                jsonrpc_response(id, json!({ "tools": tool_list() })),
            )
        }
        "tools/call" => {
            let Some(id) = id else {
                return Ok(DispatchOutcome::NoResponse);
            };
            let params = obj.get("params").cloned().unwrap_or_else(|| json!({}));
            let call = match serde_json::from_value::<ToolCallParams>(params) {
                Ok(call) => call,
                Err(err) => {
                    return write_response(
                        writer,
                        jsonrpc_error(
                            Some(id),
                            -32602,
                            "Invalid params for tools/call",
                            Some(json!({ "details": err.to_string() })),
                        ),
                    );
                }
            };
            let result = tool_call_result(default_state_path, call);
            write_response(writer, jsonrpc_response(id, result))
        }
        "shutdown" => {
            if let Some(id) = id {
                write_response(writer, jsonrpc_response(id, json!({})))
            } else {
                Ok(DispatchOutcome::NoResponse)
            }
        }
        "exit" => Ok(DispatchOutcome::Exit),
        _ => {
            if id.is_none() {
                return Ok(DispatchOutcome::NoResponse);
            }
            write_response(
                writer,
                jsonrpc_error(id, -32601, &format!("Method '{method}' not found"), None),
            )
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        engine_shell::{execute_shell_command, parse_shell_tokens},
        genomes::GenomeCatalog,
    };
    use std::{fs, io::Cursor, thread, time::Duration};
    use tempfile::tempdir;

    fn frame(value: &Value) -> Vec<u8> {
        let body = serde_json::to_vec(value).expect("serialize test message");
        let mut bytes = format!("Content-Length: {}\r\n\r\n", body.len()).into_bytes();
        bytes.extend(body);
        bytes
    }

    fn read_response_body(buffer: &[u8]) -> Value {
        let text = String::from_utf8(buffer.to_vec()).expect("utf8 response");
        let marker = "\r\n\r\n";
        let idx = text.find(marker).expect("response header separator");
        let body = &text[idx + marker.len()..];
        serde_json::from_str(body).expect("response body json")
    }

    fn run_single(default_state_path: &str, request: Value) -> Value {
        let mut reader = Cursor::new(frame(&request));
        let mut writer = Vec::<u8>::new();
        run_server_loop(default_state_path, &mut reader, &mut writer).expect("server loop");
        read_response_body(&writer)
    }

    fn run_tool(default_state_path: &str, name: &str, arguments: Value) -> Value {
        run_single(
            default_state_path,
            json!({
                "jsonrpc": "2.0",
                "id": 88,
                "method": "tools/call",
                "params": {
                    "name": name,
                    "arguments": arguments
                }
            }),
        )
    }

    fn run_shared_ui_command(tokens: Vec<String>) -> Value {
        let command = parse_shell_tokens(&tokens).expect("parse shared ui command");
        let mut engine = GentleEngine::new();
        let run = execute_shell_command(&mut engine, &command).expect("execute shared ui command");
        assert!(!run.state_changed, "ui commands must be non-mutating");
        run.output
    }

    #[test]
    fn initialize_and_tools_list_roundtrip() {
        let init = json!({
            "jsonrpc": "2.0",
            "id": 1,
            "method": "initialize",
            "params": {
                "protocolVersion": MCP_PROTOCOL_VERSION
            }
        });
        let list = json!({
            "jsonrpc": "2.0",
            "id": 2,
            "method": "tools/list",
            "params": {}
        });
        let mut input = frame(&init);
        input.extend(frame(&list));
        let mut reader = Cursor::new(input);
        let mut writer = Vec::<u8>::new();

        run_server_loop(DEFAULT_MCP_STATE_PATH, &mut reader, &mut writer).expect("server loop");

        let output = String::from_utf8(writer).expect("utf8 output");
        let parts = output
            .split("Content-Length:")
            .filter(|part| !part.trim().is_empty())
            .collect::<Vec<_>>();
        assert_eq!(parts.len(), 2);
    }

    #[test]
    fn tool_call_capabilities_returns_structured_payload() {
        let request = json!({
            "jsonrpc": "2.0",
            "id": 5,
            "method": "tools/call",
            "params": {
                "name": "capabilities",
                "arguments": {}
            }
        });
        let mut reader = Cursor::new(frame(&request));
        let mut writer = Vec::<u8>::new();

        run_server_loop(DEFAULT_MCP_STATE_PATH, &mut reader, &mut writer).expect("server loop");
        let response = read_response_body(&writer);
        let protocol = response
            .pointer("/result/structuredContent/protocol_version")
            .and_then(Value::as_str)
            .unwrap_or_default();
        assert_eq!(protocol, "v1");
    }

    #[test]
    fn tool_call_state_summary_uses_default_state_when_missing_file() {
        let request = json!({
            "jsonrpc": "2.0",
            "id": 6,
            "method": "tools/call",
            "params": {
                "name": "state_summary",
                "arguments": {}
            }
        });
        let mut reader = Cursor::new(frame(&request));
        let mut writer = Vec::<u8>::new();

        run_server_loop(
            "test_files/this_state_file_should_not_exist.gentle.json",
            &mut reader,
            &mut writer,
        )
        .expect("server loop");
        let response = read_response_body(&writer);
        let sequence_count = response
            .pointer("/result/structuredContent/sequence_count")
            .and_then(Value::as_u64)
            .unwrap_or(999);
        assert_eq!(sequence_count, 0);
    }

    #[test]
    fn parse_topic_accepts_string_and_array() {
        let one = parse_topic(Some(&json!("help capabilities"))).expect("topic from string");
        assert_eq!(
            one,
            Some(vec!["help".to_string(), "capabilities".to_string()])
        );
        let two = parse_topic(Some(&json!(["help", "state-summary"]))).expect("topic from array");
        assert_eq!(
            two,
            Some(vec!["help".to_string(), "state-summary".to_string()])
        );
    }

    #[test]
    fn tools_call_unknown_tool_returns_tool_error_payload() {
        let request = json!({
            "jsonrpc": "2.0",
            "id": 7,
            "method": "tools/call",
            "params": {
                "name": "unknown_tool",
                "arguments": {}
            }
        });
        let mut reader = Cursor::new(frame(&request));
        let mut writer = Vec::<u8>::new();
        run_server_loop(DEFAULT_MCP_STATE_PATH, &mut reader, &mut writer).expect("server loop");
        let response = read_response_body(&writer);
        let is_error = response
            .pointer("/result/isError")
            .and_then(Value::as_bool)
            .unwrap_or(false);
        assert!(is_error);
    }

    #[test]
    fn op_tool_requires_explicit_confirm_true() {
        let request = json!({
            "jsonrpc": "2.0",
            "id": 8,
            "method": "tools/call",
            "params": {
                "name": "op",
                "arguments": {
                    "operation": {
                        "SetParameter": {
                            "name": "max_fragments_per_container",
                            "value": 123
                        }
                    }
                }
            }
        });
        let response = run_single(DEFAULT_MCP_STATE_PATH, request);
        let is_error = response
            .pointer("/result/isError")
            .and_then(Value::as_bool)
            .unwrap_or(false);
        assert!(is_error);
        let text = response
            .pointer("/result/content/0/text")
            .and_then(Value::as_str)
            .unwrap_or_default();
        assert!(text.contains("confirm=true"));
    }

    #[test]
    fn op_tool_with_confirm_applies_operation_and_persists_state() {
        let temp = tempdir().expect("tempdir");
        let state_path = temp.path().join("mcp_op_state.gentle.json");
        let state_path_str = state_path.to_string_lossy().to_string();
        let request = json!({
            "jsonrpc": "2.0",
            "id": 9,
            "method": "tools/call",
            "params": {
                "name": "op",
                "arguments": {
                    "confirm": true,
                    "state_path": state_path_str,
                    "operation": {
                        "SetParameter": {
                            "name": "max_fragments_per_container",
                            "value": 123
                        }
                    }
                }
            }
        });
        let response = run_single(DEFAULT_MCP_STATE_PATH, request);
        let is_error = response
            .pointer("/result/isError")
            .and_then(Value::as_bool)
            .unwrap_or(true);
        assert!(!is_error);
        let persisted = ProjectState::load_from_path(&state_path.to_string_lossy())
            .expect("load persisted state");
        assert_eq!(persisted.parameters.max_fragments_per_container, 123);
    }

    #[test]
    fn workflow_tool_with_confirm_applies_workflow_and_persists_state() {
        let temp = tempdir().expect("tempdir");
        let state_path = temp.path().join("mcp_workflow_state.gentle.json");
        let state_path_str = state_path.to_string_lossy().to_string();
        let request = json!({
            "jsonrpc": "2.0",
            "id": 10,
            "method": "tools/call",
            "params": {
                "name": "workflow",
                "arguments": {
                    "confirm": true,
                    "state_path": state_path_str,
                    "workflow": {
                        "run_id": "mcp-test",
                        "ops": [
                            {
                                "SetParameter": {
                                    "name": "max_fragments_per_container",
                                    "value": 111
                                }
                            },
                            {
                                "SetParameter": {
                                    "name": "max_fragments_per_container",
                                    "value": 222
                                }
                            }
                        ]
                    }
                }
            }
        });
        let response = run_single(DEFAULT_MCP_STATE_PATH, request);
        let is_error = response
            .pointer("/result/isError")
            .and_then(Value::as_bool)
            .unwrap_or(true);
        assert!(!is_error);
        let result_count = response
            .pointer("/result/structuredContent/result_count")
            .and_then(Value::as_u64)
            .unwrap_or(0);
        assert_eq!(result_count, 2);
        let persisted = ProjectState::load_from_path(&state_path.to_string_lossy())
            .expect("load persisted state");
        assert_eq!(persisted.parameters.max_fragments_per_container, 222);
    }

    #[test]
    fn mcp_ui_intent_tools_match_shared_shell_outputs() {
        let td = tempdir().expect("tempdir");
        let fasta_113 = td.path().join("h113.fa");
        let ann_113 = td.path().join("h113.gtf");
        let fasta_116 = td.path().join("h116.fa");
        let ann_116 = td.path().join("h116.gtf");
        let cache_dir = td.path().join("cache");
        fs::write(&fasta_113, ">chr1\nACGT\n").expect("write fasta 113");
        fs::write(
            &ann_113,
            "chr1\tsrc\tgene\t1\t4\t.\t+\t.\tgene_id \"G113\"; gene_name \"G113\";\n",
        )
        .expect("write ann 113");
        fs::write(&fasta_116, ">chr1\nACGTACGT\n").expect("write fasta 116");
        fs::write(
            &ann_116,
            "chr1\tsrc\tgene\t1\t8\t.\t+\t.\tgene_id \"G116\"; gene_name \"G116\";\n",
        )
        .expect("write ann 116");

        let catalog_path = td.path().join("catalog.json");
        let catalog_json = format!(
            r#"{{
  "Human GRCh38 Ensembl 113": {{
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }},
  "Human GRCh38 Ensembl 116": {{
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
            fasta_113.display(),
            ann_113.display(),
            cache_dir.display(),
            fasta_116.display(),
            ann_116.display(),
            cache_dir.display()
        );
        fs::write(&catalog_path, catalog_json).expect("write catalog");
        let catalog = GenomeCatalog::from_json_file(catalog_path.to_string_lossy().as_ref())
            .expect("open catalog");
        catalog
            .prepare_genome_once("Human GRCh38 Ensembl 113")
            .expect("prepare 113");
        thread::sleep(Duration::from_millis(2));
        catalog
            .prepare_genome_once("Human GRCh38 Ensembl 116")
            .expect("prepare 116");

        let catalog_path_str = catalog_path.to_string_lossy().to_string();

        let mcp_intents = run_tool(DEFAULT_MCP_STATE_PATH, "ui_intents", json!({}));
        assert_eq!(
            mcp_intents
                .pointer("/result/isError")
                .and_then(Value::as_bool),
            Some(false)
        );
        let expected_intents = run_shared_ui_command(vec!["ui".to_string(), "intents".to_string()]);
        assert_eq!(mcp_intents["result"]["structuredContent"], expected_intents);

        let mcp_prepared = run_tool(
            DEFAULT_MCP_STATE_PATH,
            "ui_prepared_genomes",
            json!({
                "catalog_path": catalog_path_str,
                "species": "human"
            }),
        );
        assert_eq!(
            mcp_prepared
                .pointer("/result/isError")
                .and_then(Value::as_bool),
            Some(false)
        );
        let expected_prepared = run_shared_ui_command(vec![
            "ui".to_string(),
            "prepared-genomes".to_string(),
            "--catalog".to_string(),
            catalog_path.to_string_lossy().to_string(),
            "--species".to_string(),
            "human".to_string(),
        ]);
        assert_eq!(
            mcp_prepared["result"]["structuredContent"],
            expected_prepared
        );

        let mcp_latest = run_tool(
            DEFAULT_MCP_STATE_PATH,
            "ui_latest_prepared",
            json!({
                "catalog_path": catalog_path.to_string_lossy().to_string(),
                "species": "human"
            }),
        );
        assert_eq!(
            mcp_latest
                .pointer("/result/isError")
                .and_then(Value::as_bool),
            Some(false)
        );
        let expected_latest = run_shared_ui_command(vec![
            "ui".to_string(),
            "latest-prepared".to_string(),
            "human".to_string(),
            "--catalog".to_string(),
            catalog_path.to_string_lossy().to_string(),
        ]);
        assert_eq!(mcp_latest["result"]["structuredContent"], expected_latest);

        let mcp_intent = run_tool(
            DEFAULT_MCP_STATE_PATH,
            "ui_intent",
            json!({
                "action": "open",
                "target": "prepared-references",
                "catalog_path": catalog_path.to_string_lossy().to_string(),
                "species": "human",
                "latest": true
            }),
        );
        assert_eq!(
            mcp_intent
                .pointer("/result/isError")
                .and_then(Value::as_bool),
            Some(false)
        );
        let expected_intent = run_shared_ui_command(vec![
            "ui".to_string(),
            "open".to_string(),
            "prepared-references".to_string(),
            "--catalog".to_string(),
            catalog_path.to_string_lossy().to_string(),
            "--species".to_string(),
            "human".to_string(),
            "--latest".to_string(),
        ]);
        assert_eq!(mcp_intent["result"]["structuredContent"], expected_intent);
    }

    #[test]
    fn mcp_ui_intent_invalid_target_options_return_shared_parse_error() {
        let response = run_tool(
            DEFAULT_MCP_STATE_PATH,
            "ui_intent",
            json!({
                "action": "open",
                "target": "agent-assistant",
                "latest": true
            }),
        );
        assert_eq!(
            response.pointer("/result/isError").and_then(Value::as_bool),
            Some(true)
        );
        let text = response
            .pointer("/result/content/0/text")
            .and_then(Value::as_str)
            .unwrap_or_default();
        assert!(
            text.contains("only supports --helpers/--catalog/--cache-dir/--filter/--species/--latest when TARGET is prepared-references"),
            "unexpected error text: {text}"
        );
    }
}
