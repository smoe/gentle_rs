use gentle::engine::GentleEngine;
use gentle_protocol::{
    AdapterSurfacing, CapabilityAdapter, CapabilityDescriptor, CapabilitySource,
    capability_parity_adapters, capability_registry,
};
use serde::Deserialize;
use serde_json::Value;
use std::{collections::BTreeSet, process::Command};

#[derive(Debug, Deserialize)]
struct Glossary {
    commands: Vec<GlossaryCommand>,
}

#[derive(Debug, Deserialize)]
struct GlossaryCommand {
    path: String,
    summary: String,
    interfaces: Vec<String>,
}

#[derive(Debug, Deserialize)]
struct ClawBioIntents {
    routes: Vec<ClawBioRoute>,
}

#[derive(Debug, Deserialize)]
struct ClawBioRoute {
    intent_id: String,
}

fn expected_surfacing(interfaces: &[String], adapter: CapabilityAdapter) -> AdapterSurfacing {
    let interfaces = interfaces
        .iter()
        .map(String::as_str)
        .collect::<BTreeSet<_>>();
    match adapter {
        CapabilityAdapter::Gui => {
            if interfaces.contains("gui-shell") {
                AdapterSurfacing::ShellPassthrough
            } else {
                AdapterSurfacing::NotApplicable
            }
        }
        CapabilityAdapter::Cli => {
            if interfaces.contains("cli-direct") {
                AdapterSurfacing::Prominent
            } else if interfaces.contains("cli-shell") {
                AdapterSurfacing::ShellPassthrough
            } else {
                AdapterSurfacing::NotApplicable
            }
        }
        CapabilityAdapter::Mcp => {
            if interfaces.contains("mcp") {
                AdapterSurfacing::Prominent
            } else {
                AdapterSurfacing::NotApplicable
            }
        }
        CapabilityAdapter::Js => {
            if interfaces.contains("js") {
                AdapterSurfacing::Prominent
            } else {
                AdapterSurfacing::NotApplicable
            }
        }
        CapabilityAdapter::Lua => {
            if interfaces.contains("lua") {
                AdapterSurfacing::Prominent
            } else {
                AdapterSurfacing::NotApplicable
            }
        }
        CapabilityAdapter::Clawbio => AdapterSurfacing::NotApplicable,
    }
}

fn command_descriptor<'a>(
    registry: &'a [CapabilityDescriptor],
    command: &GlossaryCommand,
) -> &'a CapabilityDescriptor {
    registry
        .iter()
        .find(|descriptor| {
            descriptor.source == CapabilitySource::GlossaryCommand
                && descriptor.name == command.path
                && descriptor.description == command.summary
                && capability_parity_adapters().iter().all(|adapter| {
                    descriptor.surfacing_for_adapter(*adapter)
                        == expected_surfacing(&command.interfaces, *adapter)
                })
        })
        .unwrap_or_else(|| {
            panic!(
                "missing registry entry for glossary command `{}`",
                command.path
            )
        })
}

fn mcp_tool_names() -> BTreeSet<String> {
    gentle::mcp_server::mcp_tool_list_for_capability_surface_tests()
        .as_array()
        .expect("MCP tools/list array")
        .iter()
        .map(|tool| tool["name"].as_str().expect("MCP tool name").to_string())
        .collect()
}

fn mcp_tool_command_paths() -> BTreeSet<String> {
    gentle::mcp_server::mcp_tool_list_for_capability_surface_tests()
        .as_array()
        .expect("MCP tools/list array")
        .iter()
        .flat_map(|tool| {
            tool["commandPaths"]
                .as_array()
                .into_iter()
                .flatten()
                .filter_map(|path| path.as_str().map(ToString::to_string))
                .collect::<Vec<_>>()
        })
        .collect()
}

fn clawbio_intent_ids() -> BTreeSet<String> {
    let intents: ClawBioIntents = serde_json::from_str(include_str!(
        "../integrations/clawbio/skills/gentle-cloning/INTENTS.json"
    ))
    .expect("parse ClawBio INTENTS.json");
    intents
        .routes
        .into_iter()
        .map(|route| route.intent_id)
        .collect()
}

#[test]
fn registry_surfacing_covers_all_parity_adapters() {
    for descriptor in capability_registry() {
        for adapter in capability_parity_adapters() {
            let surfacing = descriptor.surfacing_for_adapter(*adapter);
            if surfacing == AdapterSurfacing::NotApplicable {
                let justification = descriptor
                    .surfacing_justifications
                    .get(adapter.as_str())
                    .unwrap_or_else(|| {
                        panic!(
                            "missing not-applicable justification for `{}` on {:?}",
                            descriptor.name, adapter
                        )
                    });
                assert!(
                    !justification.trim().is_empty(),
                    "empty not-applicable justification for `{}` on {:?}",
                    descriptor.name,
                    adapter
                );
            }
        }
    }
}

#[test]
fn registry_for_adapter_includes_every_reachable_projection() {
    for adapter in capability_parity_adapters() {
        let projected = gentle_protocol::capability_registry_for_adapter(*adapter)
            .into_iter()
            .map(|descriptor| (descriptor.source.as_str().to_string(), descriptor.name))
            .collect::<BTreeSet<_>>();
        let expected = capability_registry()
            .iter()
            .filter(|descriptor| {
                descriptor.surfacing_for_adapter(*adapter) != AdapterSurfacing::NotApplicable
            })
            .map(|descriptor| (descriptor.source.as_str().to_string(), descriptor.name.clone()))
            .collect::<BTreeSet<_>>();
        assert_eq!(
            projected, expected,
            "adapter projection mismatch for {:?}",
            adapter
        );
    }
}

#[test]
fn glossary_commands_have_matching_registry_entries() {
    let glossary: Glossary =
        serde_json::from_str(include_str!("../docs/glossary.json")).expect("parse glossary");
    let registry = capability_registry();

    for command in &glossary.commands {
        let descriptor = command_descriptor(registry, command);
        assert_eq!(
            descriptor.description, command.summary,
            "registry description must be sourced from docs/glossary.json for `{}`",
            command.path
        );
        for adapter in capability_parity_adapters() {
            assert_eq!(
                descriptor.surfacing_for_adapter(*adapter),
                expected_surfacing(&command.interfaces, *adapter),
                "surfacing mismatch for `{}` on {:?}",
                command.path,
                adapter
            );
        }
    }
}

#[test]
fn engine_operations_are_reachable_on_shell_route_adapters() {
    for descriptor in capability_registry()
        .iter()
        .filter(|descriptor| descriptor.source == CapabilitySource::EngineOperation)
    {
        for adapter in [
            CapabilityAdapter::Gui,
            CapabilityAdapter::Cli,
            CapabilityAdapter::Mcp,
            CapabilityAdapter::Js,
            CapabilityAdapter::Lua,
        ] {
            assert_eq!(
                descriptor.surfacing_for_adapter(adapter),
                AdapterSurfacing::ShellPassthrough,
                "engine operation `{}` must remain shell-reachable on {:?}",
                descriptor.name,
                adapter
            );
        }
    }
}

#[test]
fn prominent_projections_match_user_facing_listings() {
    let mcp_tools = mcp_tool_names();
    let mcp_command_paths = mcp_tool_command_paths();
    let clawbio_intents = clawbio_intent_ids();

    for descriptor in capability_registry() {
        for adapter in capability_parity_adapters() {
            if descriptor.surfacing_for_adapter(*adapter) != AdapterSurfacing::Prominent {
                continue;
            }
            match adapter {
                CapabilityAdapter::Cli => {
                    assert_eq!(
                        descriptor.source,
                        CapabilitySource::GlossaryCommand,
                        "prominent CLI projection `{}` must come from a CLI glossary command",
                        descriptor.name
                    );
                }
                CapabilityAdapter::Mcp => {
                    if descriptor.source == CapabilitySource::McpTool {
                        assert!(
                            mcp_tools.contains(&descriptor.name),
                            "prominent MCP tool `{}` is absent from tools/list",
                            descriptor.name
                        );
                    } else if descriptor.source == CapabilitySource::GlossaryCommand {
                        assert!(
                            mcp_command_paths.contains(&descriptor.name),
                            "prominent MCP glossary command `{}` is absent from tools/list commandPaths",
                            descriptor.name
                        );
                    }
                }
                CapabilityAdapter::Clawbio => {
                    assert!(
                        clawbio_intents.contains(&descriptor.name),
                        "prominent ClawBio projection `{}` is absent from INTENTS.json",
                        descriptor.name
                    );
                }
                CapabilityAdapter::Gui | CapabilityAdapter::Js | CapabilityAdapter::Lua => {}
            }
        }
    }
}

#[test]
fn mcp_tools_list_has_no_unregistered_prominent_tools() {
    let registered = capability_registry()
        .iter()
        .filter(|descriptor| {
            descriptor.source == CapabilitySource::McpTool
                && descriptor.surfacing_for_adapter(CapabilityAdapter::Mcp)
                    == AdapterSurfacing::Prominent
        })
        .map(|descriptor| descriptor.name.clone())
        .collect::<BTreeSet<_>>();
    assert_eq!(mcp_tool_names(), registered);
}

#[test]
fn engine_capabilities_project_protocol_registry() {
    let capabilities = GentleEngine::capabilities();
    let expected_operations = gentle_protocol::public_engine_operation_names()
        .iter()
        .map(|name| (*name).to_string())
        .collect::<Vec<_>>();
    assert_eq!(capabilities.supported_operations, expected_operations);
    assert_eq!(capabilities.capability_registry, capability_registry());
}

#[test]
fn cli_capabilities_reads_protocol_registry() {
    let output = Command::new(env!("CARGO_BIN_EXE_gentle_cli"))
        .arg("capabilities")
        .output()
        .expect("run gentle_cli capabilities");
    assert!(
        output.status.success(),
        "gentle_cli capabilities failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let payload: Value = serde_json::from_slice(&output.stdout).expect("capabilities JSON");
    assert_eq!(payload["protocol_version"].as_str(), Some("v1"));
    assert_eq!(
        payload["capability_registry"]
            .as_array()
            .expect("capability registry array")
            .len(),
        capability_registry().len()
    );
}
