use gentle_protocol::{CapabilityAdapter, CapabilitySource, EngineError};
use serde_json::{Value, json};
use std::{
    collections::BTreeSet,
    io::Write,
    process::{Command, Stdio},
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

fn run_mcp_once(request: Value) -> Value {
    let mut child = Command::new(env!("CARGO_BIN_EXE_gentle_mcp"))
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn gentle_mcp");
    child
        .stdin
        .as_mut()
        .expect("mcp stdin")
        .write_all(&frame(&request))
        .expect("write MCP frame");
    let output = child.wait_with_output().expect("wait gentle_mcp");
    assert!(
        output.status.success(),
        "gentle_mcp failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    read_framed_response(&output.stdout)
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
        cause_chain
            .iter()
            .any(|cause| *cause == "MCP adapter boundary"),
        "MCP cause chain must include adapter boundary: {cause_chain:#?}"
    );
    assert!(
        cause_chain
            .iter()
            .any(|cause| cause.contains("Unknown MCP tool 'unknown_tool'")),
        "MCP cause chain must preserve inner tool failure: {cause_chain:#?}"
    );
}
