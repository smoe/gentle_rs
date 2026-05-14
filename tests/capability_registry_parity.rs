use gentle::engine::GentleEngine;
use gentle_protocol::{
    CapabilityAdapter, CapabilityDescriptor, CapabilitySource, capability_registry,
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

fn expected_adapters(interfaces: &[String]) -> Vec<CapabilityAdapter> {
    let mut adapters = BTreeSet::new();
    for interface in interfaces {
        match interface.as_str() {
            "cli-direct" | "cli-shell" => {
                adapters.insert(CapabilityAdapter::Cli);
            }
            "mcp" => {
                adapters.insert(CapabilityAdapter::Mcp);
            }
            "js" => {
                adapters.insert(CapabilityAdapter::Js);
            }
            "lua" => {
                adapters.insert(CapabilityAdapter::Lua);
            }
            _ => {}
        }
    }
    adapters.into_iter().collect()
}

fn command_descriptor<'a>(
    registry: &'a [CapabilityDescriptor],
    command: &GlossaryCommand,
) -> &'a CapabilityDescriptor {
    let expected_adapters = expected_adapters(&command.interfaces);
    registry
        .iter()
        .find(|descriptor| {
            descriptor.source == CapabilitySource::GlossaryCommand
                && descriptor.name == command.path
                && descriptor.description == command.summary
                && descriptor.adapters == expected_adapters
        })
        .unwrap_or_else(|| {
            panic!(
                "missing registry entry for glossary command `{}`",
                command.path
            )
        })
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
        assert_eq!(
            descriptor.adapters,
            expected_adapters(&command.interfaces),
            "adapter projection mismatch for `{}`",
            command.path
        );
    }
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
