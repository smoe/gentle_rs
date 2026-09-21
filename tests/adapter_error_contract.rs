//! Real-binary adapter boundary tests with deterministic synthetic requests and
//! an inline BamHI sequence saved only in a temporary project. No network or
//! external fixtures are required; requests exercise stdout framing, errors,
//! state persistence and MCP startup independently of Cargo's stack environment.

use gentle::{dna_sequence::DNAsequence, engine::ProjectState};
use gentle_protocol::{CapabilityAdapter, CapabilitySource, EngineError};
use serde_json::{Value, json};
use std::{
    collections::BTreeSet,
    io::Write,
    process::{Command, Output, Stdio},
};

fn frame(payload: &Value) -> Vec<u8> {
    let body = serde_json::to_vec(payload).expect("request json");
    let mut framed = format!("Content-Length: {}\r\n\r\n", body.len()).into_bytes();
    framed.extend(body);
    framed
}

fn read_framed_response(output: &[u8]) -> Value {
    let text = String::from_utf8_lossy(output);
    let split_at = text
        .find("\r\n\r\n")
        .expect("MCP response has header terminator")
        + 4;
    serde_json::from_str(&text[split_at..]).expect("MCP JSON body")
}

fn run_mcp_input(input: &[u8]) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_gentle_mcp"))
        .env_remove("RUST_MIN_STACK")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn gentle_mcp");
    child
        .stdin
        .as_mut()
        .expect("mcp stdin")
        .write_all(input)
        .expect("write MCP frame");
    child.wait_with_output().expect("wait gentle_mcp")
}

fn run_mcp_once(request: Value) -> Value {
    let output = run_mcp_input(&frame(&request));
    assert!(
        output.status.success(),
        "gentle_mcp failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        output.stdout.starts_with(b"Content-Length: "),
        "MCP stdout must contain framed protocol output only: {}",
        String::from_utf8_lossy(&output.stdout)
    );
    read_framed_response(&output.stdout)
}

#[test]
fn mcp_digest_keeps_stdout_protocol_framed() {
    let temp = tempfile::tempdir().expect("tempdir");
    let state_path = temp.path().join("digest_state.gentle.json");
    let mut state = ProjectState::default();
    state.sequences.insert(
        "digest_input".to_string(),
        DNAsequence::from_sequence("ATGGATCCGCATGGATCCGC").expect("digest input sequence"),
    );
    state
        .save_to_path(&state_path.to_string_lossy())
        .expect("save digest state");

    let response = run_mcp_once(json!({
        "jsonrpc": "2.0",
        "id": 77,
        "method": "tools/call",
        "params": {
            "name": "op",
            "arguments": {
                "confirm": true,
                "state_path": state_path.to_string_lossy(),
                "operation": {
                    "Digest": {
                        "input": "digest_input",
                        "enzymes": ["BamHI"],
                        "output_prefix": "digest_fragment"
                    }
                }
            }
        }
    }));
    assert_eq!(
        response.pointer("/result/isError").and_then(Value::as_bool),
        Some(false)
    );
    assert_eq!(
        response
            .pointer("/result/structuredContent/result/created_seq_ids")
            .and_then(Value::as_array)
            .map(Vec::len),
        Some(3)
    );
    let persisted = ProjectState::load_from_path(&state_path.to_string_lossy())
        .expect("load persisted digest state");
    for seq_id in response["result"]["structuredContent"]["result"]["created_seq_ids"]
        .as_array()
        .expect("created digest fragments")
    {
        assert!(
            persisted
                .sequences
                .contains_key(seq_id.as_str().expect("fragment id"))
        );
    }
}

#[test]
fn mcp_framing_failure_exits_without_stdout_noise() {
    let output = run_mcp_input(b"Content-Length: invalid\r\n\r\n");
    assert!(!output.status.success());
    assert!(output.stdout.is_empty(), "MCP errors must stay off stdout");
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("Invalid Content-Length"),
        "MCP worker error must reach stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn cli_error_boundary_emits_engine_error_payload() {
    let output = Command::new(env!("CARGO_BIN_EXE_gentle_cli"))
        .arg("definitely-not-a-command")
        .output()
        .expect("run failing gentle_cli command");
    assert!(!output.status.success());

    let payload: Value = serde_json::from_slice(&output.stderr).expect("CLI stderr JSON error");
    let _error: EngineError = serde_json::from_value(payload.clone()).expect("EngineError payload");
    assert_eq!(payload["code"].as_str(), Some("InvalidInput"));
    assert_eq!(
        payload["message"].as_str(),
        Some("CLI adapter command failed")
    );
    let cause_chain = payload["cause_chain"]
        .as_array()
        .expect("cause_chain array")
        .iter()
        .map(|value| value.as_str().expect("cause string"))
        .collect::<Vec<_>>();
    assert!(
        cause_chain
            .iter()
            .any(|cause| cause.contains("definitely-not-a-command")),
        "CLI cause chain must preserve the failed command: {cause_chain:#?}"
    );
}

#[test]
fn mcp_tools_list_projects_protocol_registry_names() {
    let response = run_mcp_once(json!({
        "jsonrpc": "2.0",
        "id": 1,
        "method": "tools/list",
        "params": {}
    }));
    let tools = response["result"]["tools"]
        .as_array()
        .expect("tools/list array");
    let actual = tools
        .iter()
        .map(|tool| tool["name"].as_str().expect("tool name").to_string())
        .collect::<BTreeSet<_>>();
    let expected = gentle_protocol::capability_registry_for_adapter(CapabilityAdapter::Mcp)
        .into_iter()
        .filter(|descriptor| descriptor.source == CapabilitySource::McpTool)
        .map(|descriptor| descriptor.name)
        .collect::<BTreeSet<_>>();
    assert_eq!(actual, expected);
}

#[test]
fn mcp_tool_error_boundary_embeds_engine_error_payload() {
    let response = run_mcp_once(json!({
        "jsonrpc": "2.0",
        "id": 2,
        "method": "tools/call",
        "params": {
            "name": "unknown_tool",
            "arguments": {}
        }
    }));
    assert_eq!(response["result"]["isError"].as_bool(), Some(true));
    let payload = &response["result"]["structuredContent"]["error"];
    let _error: EngineError =
        serde_json::from_value(payload.clone()).expect("EngineError-shaped MCP error");
    assert_eq!(payload["code"].as_str(), Some("NotFound"));
    let cause_chain = payload["cause_chain"]
        .as_array()
        .expect("cause_chain array")
        .iter()
        .map(|value| value.as_str().expect("cause string"))
        .collect::<Vec<_>>();
    assert!(
        cause_chain.contains(&"MCP adapter boundary"),
        "MCP cause chain must include adapter boundary: {cause_chain:#?}"
    );
    assert!(
        cause_chain
            .iter()
            .any(|cause| cause.contains("Unknown MCP tool 'unknown_tool'")),
        "MCP cause chain must preserve inner tool failure: {cause_chain:#?}"
    );
}
