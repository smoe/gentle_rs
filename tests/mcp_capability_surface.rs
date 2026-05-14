use serde_json::Value;
use std::collections::{BTreeMap, BTreeSet};

const GLOSSARY_JSON: &str = include_str!("../docs/glossary.json");
const AGENT_INTERFACE_MD: &str = include_str!("../docs/agent_interface.md");

fn glossary_commands() -> Vec<Value> {
    let glossary: Value = serde_json::from_str(GLOSSARY_JSON).expect("valid glossary json");
    glossary["commands"]
        .as_array()
        .expect("glossary commands array")
        .clone()
}

fn glossary_summaries_by_path() -> BTreeMap<String, String> {
    let mut summaries = BTreeMap::new();
    for command in glossary_commands() {
        let Some(path) = command["path"].as_str() else {
            continue;
        };
        let Some(summary) = command["summary"].as_str() else {
            continue;
        };
        summaries
            .entry(path.to_string())
            .or_insert_with(|| summary.to_string());
    }
    summaries
}

fn shell_command_paths() -> BTreeSet<String> {
    let shell_interfaces = ["cli-shell", "gui-shell", "cli-direct"]
        .into_iter()
        .collect::<BTreeSet<_>>();
    glossary_commands()
        .into_iter()
        .filter(|command| {
            command["interfaces"].as_array().is_some_and(|interfaces| {
                interfaces
                    .iter()
                    .filter_map(Value::as_str)
                    .any(|interface| shell_interfaces.contains(interface))
            })
        })
        .filter_map(|command| command["path"].as_str().map(ToString::to_string))
        .collect()
}

fn intentionally_mcp_excluded_paths() -> BTreeSet<String> {
    let marker = "#### Intentionally MCP-excluded shell commands";
    let section = AGENT_INTERFACE_MD
        .split(marker)
        .nth(1)
        .expect("intentional MCP exclusion section");
    let section = section.split("\n### ").next().unwrap_or(section);
    section
        .lines()
        .filter_map(|line| line.trim().strip_prefix("- `"))
        .filter_map(|rest| rest.split_once('`').map(|(path, _)| path.to_string()))
        .collect()
}

fn mcp_tools() -> Vec<Value> {
    gentle::mcp_server::mcp_tool_list_for_capability_surface_tests()
        .as_array()
        .expect("MCP tools array")
        .clone()
}

fn tool_command_paths(tools: &[Value]) -> BTreeSet<String> {
    tools
        .iter()
        .flat_map(|tool| {
            tool["commandPaths"]
                .as_array()
                .into_iter()
                .flatten()
                .filter_map(Value::as_str)
                .map(ToString::to_string)
                .collect::<Vec<_>>()
        })
        .collect()
}

fn assert_engine_error_payload(value: &Value) {
    assert!(
        value["code"].as_str().is_some(),
        "EngineError.code must be present: {value:#}"
    );
    assert!(
        value["message"].as_str().is_some(),
        "EngineError.message must be present: {value:#}"
    );
}

#[test]
fn tools_list_has_complete_capability_metadata() {
    let summaries = glossary_summaries_by_path();
    let tools = mcp_tools();
    assert!(!tools.is_empty(), "MCP tools/list must not be empty");

    for tool in &tools {
        let name = tool["name"]
            .as_str()
            .expect("tool descriptor has stable name");
        assert!(
            tool["description"]
                .as_str()
                .is_some_and(|value| !value.is_empty()),
            "tool {name} must have a non-empty description"
        );
        assert!(
            tool["inputSchema"].is_object(),
            "tool {name} must expose an inputSchema"
        );
        assert!(
            tool["outputSchema"].is_object(),
            "tool {name} must expose an outputSchema"
        );
        assert!(
            tool["outputSchema"]
                .pointer("/properties/error/properties/code")
                .is_some(),
            "tool {name} outputSchema must document EngineError.code"
        );
        match &tool["mutating"] {
            Value::Bool(_) => {}
            Value::String(value) if value == "external" => {}
            other => panic!("tool {name} has invalid mutating descriptor: {other:#}"),
        }

        if let Some(command_paths) = tool["commandPaths"].as_array() {
            let first_path = command_paths
                .first()
                .and_then(Value::as_str)
                .expect("commandPaths entries are strings");
            let expected = summaries.get(first_path).unwrap_or_else(|| {
                panic!("tool {name} references unknown glossary path {first_path}")
            });
            assert_eq!(
                tool["description"].as_str(),
                Some(expected.as_str()),
                "tool {name} description must come from docs/glossary.json"
            );
            for path in command_paths {
                let path = path.as_str().expect("command path string");
                assert!(
                    summaries.contains_key(path),
                    "tool {name} references unknown glossary path {path}"
                );
            }
        }
    }
}

#[test]
fn every_glossary_shell_command_is_mcp_tool_or_documented_exclusion() {
    let tools = mcp_tools();
    let covered = tool_command_paths(&tools);
    let excluded = intentionally_mcp_excluded_paths();
    let shell_paths = shell_command_paths();

    let missing = shell_paths
        .iter()
        .filter(|path| !covered.contains(*path) && !excluded.contains(*path))
        .cloned()
        .collect::<Vec<_>>();
    assert!(
        missing.is_empty(),
        "shell commands must be covered by tools/list commandPaths or documented as intentionally MCP-excluded: {missing:#?}"
    );

    let stale_exclusions = excluded
        .iter()
        .filter(|path| !shell_paths.contains(*path))
        .cloned()
        .collect::<Vec<_>>();
    assert!(
        stale_exclusions.is_empty(),
        "intentional MCP exclusions must refer to current glossary shell commands: {stale_exclusions:#?}"
    );
}

#[test]
fn mcp_tool_errors_return_typed_engine_error_payloads() {
    let cases = [
        gentle::mcp_server::mcp_tool_call_for_capability_surface_tests(
            gentle::mcp_server::DEFAULT_MCP_STATE_PATH,
            "unknown_tool",
            serde_json::json!({}),
        ),
        gentle::mcp_server::mcp_tool_call_for_capability_surface_tests(
            gentle::mcp_server::DEFAULT_MCP_STATE_PATH,
            "op",
            serde_json::json!({}),
        ),
        gentle::mcp_server::mcp_tool_call_for_capability_surface_tests(
            gentle::mcp_server::DEFAULT_MCP_STATE_PATH,
            "help",
            serde_json::json!({ "format": "bogus" }),
        ),
        gentle::mcp_server::mcp_tool_call_for_capability_surface_tests(
            gentle::mcp_server::DEFAULT_MCP_STATE_PATH,
            "restriction_site_detail",
            serde_json::json!({}),
        ),
    ];

    for result in cases {
        assert_eq!(result["isError"], Value::Bool(true), "{result:#}");
        assert_engine_error_payload(&result["structuredContent"]["error"]);
    }
}

#[test]
fn mcp_jsonrpc_errors_carry_typed_engine_error_data() {
    let response = gentle::mcp_server::mcp_jsonrpc_error_for_capability_surface_tests(
        -32602,
        "Invalid params for tools/call",
        Some(serde_json::json!({ "details": "bad params" })),
    );
    assert_engine_error_payload(&response["error"]["data"]["engine_error"]);
}
