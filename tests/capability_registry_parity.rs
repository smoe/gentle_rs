use gentle::{app::gui_prominent_glossary_entries, engine::GentleEngine};
use gentle_protocol::{
    AdapterSurfacing, CapabilityAdapter, CapabilityDescriptor, CapabilitySource,
    capability_parity_adapters, capability_registry,
};
use serde::Deserialize;
use serde_json::Value;
use std::{
    collections::{BTreeMap, BTreeSet},
    process::Command,
};

#[derive(Debug, Deserialize)]
struct Glossary {
    commands: Vec<GlossaryCommand>,
}

#[derive(Debug, Deserialize)]
struct GlossaryCommand {
    path: String,
    summary: String,
    interfaces: Vec<String>,
    #[serde(default)]
    engine_operations: Vec<String>,
}

#[derive(Debug, Deserialize)]
struct ParityMatrixOverrides {
    overrides: Vec<ParityMatrixOverride>,
}

#[derive(Debug, Deserialize)]
struct ParityMatrixOverride {
    source: CapabilitySource,
    name: String,
    adapter: CapabilityAdapter,
    surfacing: AdapterSurfacing,
    reason: String,
}

#[derive(Debug, Deserialize)]
struct ClawBioIntents {
    routes: Vec<ClawBioRoute>,
}

#[derive(Debug, Deserialize)]
struct ClawBioRoute {
    intent_id: String,
}

fn glossary_commands() -> Vec<GlossaryCommand> {
    let glossary: Glossary =
        serde_json::from_str(include_str!("../docs/glossary.json")).expect("parse glossary");
    glossary.commands
}

fn interfaces_by_path(commands: &[GlossaryCommand]) -> BTreeMap<String, BTreeSet<String>> {
    let mut interfaces_by_path = BTreeMap::new();
    for command in commands {
        let interfaces = interfaces_by_path
            .entry(command.path.clone())
            .or_insert_with(BTreeSet::new);
        interfaces.extend(command.interfaces.iter().cloned());
    }
    interfaces_by_path
}

fn parity_matrix_overrides() -> Vec<ParityMatrixOverride> {
    let parsed: ParityMatrixOverrides =
        serde_json::from_str(include_str!("../docs/parity_matrix_overrides.json"))
            .expect("parse parity matrix overrides");
    parsed.overrides
}

fn override_reason<'a>(
    overrides: &'a [ParityMatrixOverride],
    source: CapabilitySource,
    name: &str,
    adapter: CapabilityAdapter,
) -> Option<&'a str> {
    overrides
        .iter()
        .find(|entry| entry.source == source && entry.name == name && entry.adapter == adapter)
        .or_else(|| {
            overrides.iter().find(|entry| {
                entry.source == source && entry.name == "*" && entry.adapter == adapter
            })
        })
        .map(|entry| entry.reason.as_str())
}

fn expected_surfacing_without_overrides(
    command: &GlossaryCommand,
    adapter: CapabilityAdapter,
    interfaces_by_path: &BTreeMap<String, BTreeSet<String>>,
    mcp_command_paths: &BTreeSet<String>,
) -> AdapterSurfacing {
    let interfaces = interfaces_by_path
        .get(&command.path)
        .expect("glossary path interface index");
    match adapter {
        CapabilityAdapter::Gui => {
            if interfaces.contains("gui-menu") {
                AdapterSurfacing::Prominent
            } else if interfaces.contains("gui-shell") {
                AdapterSurfacing::ShellPassthrough
            } else {
                AdapterSurfacing::Gap
            }
        }
        CapabilityAdapter::Cli => {
            if interfaces.contains("cli-direct") {
                AdapterSurfacing::Prominent
            } else if interfaces.contains("cli-shell") || !command.engine_operations.is_empty() {
                AdapterSurfacing::ShellPassthrough
            } else {
                AdapterSurfacing::Gap
            }
        }
        CapabilityAdapter::Mcp => {
            if interfaces.contains("mcp") || mcp_command_paths.contains(&command.path) {
                AdapterSurfacing::Prominent
            } else if !command.engine_operations.is_empty() {
                AdapterSurfacing::ShellPassthrough
            } else {
                AdapterSurfacing::Gap
            }
        }
        CapabilityAdapter::Js => {
            if interfaces.contains("js") {
                AdapterSurfacing::Prominent
            } else if !command.engine_operations.is_empty() {
                AdapterSurfacing::ShellPassthrough
            } else {
                AdapterSurfacing::Gap
            }
        }
        CapabilityAdapter::Lua => {
            if interfaces.contains("lua") {
                AdapterSurfacing::Prominent
            } else if !command.engine_operations.is_empty() {
                AdapterSurfacing::ShellPassthrough
            } else {
                AdapterSurfacing::Gap
            }
        }
        CapabilityAdapter::Clawbio => AdapterSurfacing::Gap,
    }
}

fn expected_surfacing(
    command: &GlossaryCommand,
    adapter: CapabilityAdapter,
    interfaces_by_path: &BTreeMap<String, BTreeSet<String>>,
    mcp_command_paths: &BTreeSet<String>,
    overrides: &[ParityMatrixOverride],
) -> AdapterSurfacing {
    let surfacing = expected_surfacing_without_overrides(
        command,
        adapter,
        interfaces_by_path,
        mcp_command_paths,
    );
    if override_reason(
        overrides,
        CapabilitySource::GlossaryCommand,
        &command.path,
        adapter,
    )
    .is_some()
    {
        assert!(
            !surfacing.is_reachable(),
            "override for `{}` on {:?} must not demote {:?}",
            command.path,
            adapter,
            surfacing
        );
        AdapterSurfacing::NotApplicable
    } else {
        surfacing
    }
}

fn command_descriptor<'a>(
    registry: &'a [CapabilityDescriptor],
    command: &GlossaryCommand,
    interfaces_by_path: &BTreeMap<String, BTreeSet<String>>,
    mcp_command_paths: &BTreeSet<String>,
    overrides: &[ParityMatrixOverride],
) -> &'a CapabilityDescriptor {
    registry
        .iter()
        .find(|descriptor| {
            descriptor.source == CapabilitySource::GlossaryCommand
                && descriptor.name == command.path
                && descriptor.description == command.summary
                && capability_parity_adapters().iter().all(|adapter| {
                    descriptor.surfacing_for_adapter(*adapter)
                        == expected_surfacing(
                            command,
                            *adapter,
                            interfaces_by_path,
                            mcp_command_paths,
                            overrides,
                        )
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

fn open_gap_rows_from_matrix() -> BTreeSet<(String, String, String)> {
    let markdown = gentle_protocol::render_gui_cli_mcp_parity_matrix_markdown();
    let section = markdown
        .split("## Open Gaps")
        .nth(1)
        .expect("open gaps section");
    section
        .lines()
        .filter(|line| line.starts_with('|') && !line.contains("---") && !line.contains("Adapter"))
        .filter_map(|line| {
            let cells = line
                .trim_matches('|')
                .split('|')
                .map(|cell| cell.trim().to_string())
                .collect::<Vec<_>>();
            (cells.len() >= 3 && cells[0] != "_None_")
                .then(|| (cells[0].clone(), cells[1].clone(), cells[2].clone()))
        })
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
                if descriptor.source == CapabilitySource::GlossaryCommand
                    && !descriptor.engine_operations.is_empty()
                {
                    assert!(
                        !justification.contains("does not declare"),
                        "engine-backed n/a for `{}` on {:?} must use a curated reason, not an auto-generated absence note",
                        descriptor.name,
                        adapter
                    );
                }
            }
        }
    }
}

#[test]
fn not_applicable_cells_have_curated_override_entries() {
    let overrides = parity_matrix_overrides();
    for descriptor in capability_registry() {
        for adapter in capability_parity_adapters() {
            if descriptor.surfacing_for_adapter(*adapter) != AdapterSurfacing::NotApplicable {
                continue;
            }
            let override_reason =
                override_reason(&overrides, descriptor.source, &descriptor.name, *adapter)
                    .unwrap_or_else(|| {
                        panic!(
                            "not-applicable cell `{}` on {:?} must have a parity override",
                            descriptor.name, adapter
                        )
                    });
            assert!(
                !override_reason.trim().is_empty(),
                "not-applicable override for `{}` on {:?} must have a reason",
                descriptor.name,
                adapter
            );
            let justification = descriptor
                .surfacing_justifications
                .get(adapter.as_str())
                .unwrap_or_else(|| {
                    panic!(
                        "missing not-applicable justification for `{}` on {:?}",
                        descriptor.name, adapter
                    )
                });
            assert_eq!(
                justification, override_reason,
                "not-applicable justification for `{}` on {:?} must come from the override table",
                descriptor.name, adapter
            );
        }
    }
}

#[test]
fn parity_overrides_do_not_demote_prominent_or_shell_routes() {
    let commands = glossary_commands();
    let interfaces_by_path = interfaces_by_path(&commands);
    let mcp_command_paths = mcp_tool_command_paths();
    let overrides = parity_matrix_overrides();

    for entry in &overrides {
        assert_eq!(
            entry.surfacing,
            AdapterSurfacing::NotApplicable,
            "only explicit not-applicable overrides are supported"
        );
        assert!(
            !entry.reason.trim().is_empty(),
            "override for {:?} `{}` on {:?} must have a reason",
            entry.source,
            entry.name,
            entry.adapter
        );
        match entry.source {
            CapabilitySource::GlossaryCommand => {
                let matching = commands
                    .iter()
                    .filter(|command| entry.name == "*" || command.path == entry.name)
                    .collect::<Vec<_>>();
                assert!(
                    !matching.is_empty(),
                    "glossary override references unknown command `{}`",
                    entry.name
                );
                for command in matching {
                    let default = expected_surfacing_without_overrides(
                        command,
                        entry.adapter,
                        &interfaces_by_path,
                        &mcp_command_paths,
                    );
                    assert!(
                        !default.is_reachable(),
                        "override for glossary command `{}` on {:?} would demote {:?}",
                        command.path,
                        entry.adapter,
                        default
                    );
                }
            }
            CapabilitySource::EngineOperation => {
                let matching = gentle_protocol::public_engine_operation_names()
                    .iter()
                    .filter(|name| entry.name == "*" || **name == entry.name.as_str())
                    .collect::<Vec<_>>();
                assert!(
                    !matching.is_empty(),
                    "engine override references unknown operation `{}`",
                    entry.name
                );
                let default = match entry.adapter {
                    CapabilityAdapter::Clawbio => AdapterSurfacing::NotApplicable,
                    CapabilityAdapter::Gui
                    | CapabilityAdapter::Cli
                    | CapabilityAdapter::Mcp
                    | CapabilityAdapter::Js
                    | CapabilityAdapter::Lua => AdapterSurfacing::ShellPassthrough,
                };
                assert!(
                    !default.is_reachable(),
                    "override for engine operation `{}` on {:?} would demote {:?}",
                    entry.name,
                    entry.adapter,
                    default
                );
            }
            CapabilitySource::McpTool => {
                let matching = capability_registry()
                    .iter()
                    .filter(|descriptor| {
                        descriptor.source == CapabilitySource::McpTool
                            && (entry.name == "*" || descriptor.name == entry.name)
                    })
                    .collect::<Vec<_>>();
                assert!(
                    !matching.is_empty(),
                    "MCP tool override references unknown tool `{}`",
                    entry.name
                );
                let default = match entry.adapter {
                    CapabilityAdapter::Mcp => AdapterSurfacing::Prominent,
                    CapabilityAdapter::Gui
                    | CapabilityAdapter::Cli
                    | CapabilityAdapter::Js
                    | CapabilityAdapter::Lua
                    | CapabilityAdapter::Clawbio => AdapterSurfacing::NotApplicable,
                };
                assert!(
                    !default.is_reachable(),
                    "override for MCP tool `{}` on {:?} would demote {:?}",
                    entry.name,
                    entry.adapter,
                    default
                );
            }
        }
    }
}

#[test]
fn gap_rows_are_listed_in_open_gap_section() {
    let open_gap_rows = open_gap_rows_from_matrix();
    let mut gap_count = 0usize;
    for descriptor in capability_registry() {
        for adapter in capability_parity_adapters() {
            if descriptor.surfacing_for_adapter(*adapter) != AdapterSurfacing::Gap {
                continue;
            }
            gap_count += 1;
            assert_eq!(
                descriptor.source,
                CapabilitySource::GlossaryCommand,
                "only glossary adapter routes should appear as open gaps"
            );
            let row = (
                adapter.label().to_string(),
                descriptor.name.clone(),
                descriptor.source.as_str().to_string(),
            );
            assert!(
                open_gap_rows.contains(&row),
                "gap `{}` on {:?} must be listed in ## Open Gaps",
                descriptor.name,
                adapter
            );
            assert!(
                descriptor
                    .surfacing_justifications
                    .get(adapter.as_str())
                    .is_none(),
                "gap `{}` on {:?} must not carry a not-applicable justification",
                descriptor.name,
                adapter
            );
        }
    }
    assert!(
        gap_count > 0,
        "credible parity matrix should surface current adapter gaps"
    );
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
            .filter(|descriptor| descriptor.is_reachable_from_adapter(*adapter))
            .map(|descriptor| {
                (
                    descriptor.source.as_str().to_string(),
                    descriptor.name.clone(),
                )
            })
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
    let commands = glossary_commands();
    let interfaces_by_path = interfaces_by_path(&commands);
    let mcp_command_paths = mcp_tool_command_paths();
    let overrides = parity_matrix_overrides();
    let registry = capability_registry();

    for command in &commands {
        let descriptor = command_descriptor(
            registry,
            command,
            &interfaces_by_path,
            &mcp_command_paths,
            &overrides,
        );
        assert_eq!(
            descriptor.description, command.summary,
            "registry description must be sourced from docs/glossary.json for `{}`",
            command.path
        );
        for adapter in capability_parity_adapters() {
            assert_eq!(
                descriptor.surfacing_for_adapter(*adapter),
                expected_surfacing(
                    command,
                    *adapter,
                    &interfaces_by_path,
                    &mcp_command_paths,
                    &overrides,
                ),
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
    let gui_prominent_paths = gui_prominent_glossary_entries()
        .iter()
        .map(|entry| entry.glossary_path)
        .collect::<BTreeSet<_>>();

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
                CapabilityAdapter::Gui => {
                    assert_eq!(
                        descriptor.source,
                        CapabilitySource::GlossaryCommand,
                        "prominent GUI projection `{}` must come from a GUI glossary command",
                        descriptor.name
                    );
                    assert!(
                        gui_prominent_paths.contains(descriptor.name.as_str()),
                        "prominent GUI glossary command `{}` is absent from the GUI menu/palette oracle",
                        descriptor.name
                    );
                }
                CapabilityAdapter::Js | CapabilityAdapter::Lua => {}
            }
        }
    }

    for entry in gui_prominent_glossary_entries() {
        let descriptor = capability_registry()
            .iter()
            .find(|descriptor| {
                descriptor.source == CapabilitySource::GlossaryCommand
                    && descriptor.name == entry.glossary_path
            })
            .unwrap_or_else(|| {
                panic!(
                    "GUI menu/palette oracle references missing glossary command `{}`",
                    entry.glossary_path
                )
            });
        assert_eq!(
            descriptor.surfacing_for_adapter(CapabilityAdapter::Gui),
            AdapterSurfacing::Prominent,
            "GUI oracle command `{}` must project as GUI prominent",
            entry.glossary_path
        );
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
