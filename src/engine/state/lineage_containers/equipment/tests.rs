//! Synthetic in-memory racks only; no manufacturer artifact or physical acceptance.

use super::*;
use crate::engine_shell::{ShellCommand, execute_shell_command, parse_shell_line, quote_shell_arg};

fn engine(rows: usize, columns: usize) -> GentleEngine {
    let mut state = ProjectState::default();
    state.container_state.racks.insert(
        "rack 'one'".into(),
        Rack {
            rack_id: "rack 'one'".into(),
            name: "PRIVATE SAMPLE NAME".into(),
            profile: RackProfileSnapshot::custom(rows, columns),
            ..Rack::default()
        },
    );
    GentleEngine::from_state(state)
}

const TEMPLATE: RackPhysicalTemplateKind = RackPhysicalTemplateKind::SafeLockTubeRackFloorSeated;

#[test]
fn rack_equipment_frames_match_existing_simulation_in_metres_for_every_template() {
    let mut engine = engine(8, 12);
    let temp = tempfile::tempdir().unwrap();
    let before = serde_json::to_value(engine.state()).unwrap();
    let history = engine.history_summary();
    for (index, template) in [
        RackPhysicalTemplateKind::StoragePcrTubeRack,
        RackPhysicalTemplateKind::PipettingPcrTubeRack,
        RackPhysicalTemplateKind::CellCulturePlate,
        RackPhysicalTemplateKind::SafeLockTubeRackCollarSeated,
        TEMPLATE,
        RackPhysicalTemplateKind::PcrTube96PositionCarrier,
        RackPhysicalTemplateKind::SafeLockStorageBoxDerived,
    ]
    .into_iter()
    .enumerate()
    {
        let legacy_path = temp.path().join(format!("legacy-{index}.json"));
        engine
            .export_rack_simulation_json("rack 'one'", template, legacy_path.to_str().unwrap())
            .unwrap();
        let legacy_bytes = fs::read(&legacy_path).unwrap();
        let legacy: serde_json::Value = serde_json::from_slice(&legacy_bytes).unwrap();
        let description = engine
            .rack_equipment_description("rack 'one'", template)
            .unwrap();
        assert_eq!(description.schema, RACK_EQUIPMENT_SCHEMA);
        assert_eq!(description.slots.len(), 96);
        assert_eq!(description.length_unit, "metre");
        assert_eq!(description.base_frame.parent, None);
        assert_eq!(description.base_frame.rotation_xyzw, [0.0, 0.0, 0.0, 1.0]);
        for (slot, original) in description
            .slots
            .iter()
            .zip(legacy["slots"].as_array().unwrap())
        {
            assert_eq!(slot.coordinate, original["coordinate"].as_str().unwrap());
            assert_eq!(slot.row, original["row"].as_u64().unwrap() as usize);
            assert_eq!(slot.column, original["column"].as_u64().unwrap() as usize);
            assert_eq!(slot.frame.parent.as_deref(), Some("rack_base"));
            // Compare the existing serialized f32 centre through its typed value.
            for (axis, name) in [(0, "x"), (1, "y")] {
                let mm: f32 = serde_json::from_value(original["center_mm"][name].clone()).unwrap();
                assert_eq!(slot.frame.translation_m[axis], f64::from(mm) / 1000.0);
            }
            assert_eq!(
                slot.frame.translation_m[2],
                f64::from(description.source.physical_template.rack_height_mm) / 1000.0
            );
        }
        GentleEngine::verify_rack_equipment_description(&description).unwrap();
        let path = temp.path().join(format!("equipment-{index}.json"));
        let result = engine
            .apply(Operation::ExportRackEquipmentJson {
                rack_id: "rack 'one'".into(),
                path: path.to_str().unwrap().into(),
                template,
            })
            .unwrap();
        let json = serde_json::to_value(result).unwrap();
        let typed_json: serde_json::Value =
            serde_json::from_slice(&serde_json::to_vec(&description).unwrap()).unwrap();
        assert_eq!(typed_json, serde_json::to_value(&description).unwrap());
        assert_eq!(
            json["rack_equipment"],
            serde_json::to_value(&description).unwrap()
        );
        assert_eq!(
            fs::read(&path).unwrap(),
            serde_json::to_vec_pretty(&serde_json::to_value(&description).unwrap()).unwrap()
        );
        let reopened: RackEquipmentDescription =
            serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
        GentleEngine::verify_rack_equipment_description(&reopened).unwrap();
        let second = temp.path().join(format!("legacy-again-{index}.json"));
        engine
            .export_rack_simulation_json("rack 'one'", template, second.to_str().unwrap())
            .unwrap();
        assert_eq!(
            legacy_bytes,
            fs::read(second).unwrap(),
            "legacy bytes changed"
        );
    }
    assert_eq!(serde_json::to_value(engine.state()).unwrap(), before);
    assert_eq!(engine.history_summary().undo_count, history.undo_count);
}

#[test]
fn rack_equipment_hashes_are_deterministic_and_bind_geometry_and_blocking() {
    let mut engine = engine(3, 4);
    let baseline = engine
        .rack_equipment_description("rack 'one'", TEMPLATE)
        .unwrap();
    assert_eq!(
        serde_json::to_string(&baseline).unwrap(),
        serde_json::to_string(
            &engine
                .rack_equipment_description("rack 'one'", TEMPLATE)
                .unwrap()
        )
        .unwrap()
    );
    let rack = engine
        .state
        .container_state
        .racks
        .get_mut("rack 'one'")
        .unwrap();
    rack.profile.blocked_coordinates = vec!["B2".into(), "A1".into()];
    let blocked = engine
        .rack_equipment_description("rack 'one'", TEMPLATE)
        .unwrap();
    assert_eq!(blocked.slots.len(), baseline.slots.len());
    assert_eq!(blocked.slots.iter().filter(|s| s.blocked).count(), 2);
    assert_eq!(
        blocked.asset_revision_sha256,
        baseline.asset_revision_sha256
    );
    assert_ne!(blocked.source_sha256, baseline.source_sha256);
    assert_ne!(blocked.content_sha256, baseline.content_sha256);
    engine
        .state
        .container_state
        .racks
        .get_mut("rack 'one'")
        .unwrap()
        .profile
        .blocked_coordinates
        .reverse();
    assert_eq!(
        blocked.content_sha256,
        engine
            .rack_equipment_description("rack 'one'", TEMPLATE)
            .unwrap()
            .content_sha256
    );
    engine
        .state
        .container_state
        .racks
        .get_mut("rack 'one'")
        .unwrap()
        .carrier_pitch_override_um = Some(14_000);
    let compact = engine
        .rack_equipment_description("rack 'one'", TEMPLATE)
        .unwrap();
    assert_ne!(compact.asset_revision_sha256, blocked.asset_revision_sha256);
    assert_ne!(
        compact.slots[1].frame.translation_m,
        blocked.slots[1].frame.translation_m
    );
    assert_eq!(compact.source.tube_profiles.len(), 2);
    assert_eq!(compact.mass_kg, None);
    assert_eq!(compact.inertia_kg_m2, None);
    assert!(compact.permitted_actions.is_empty());
    assert!(
        !serde_json::to_string(&compact)
            .unwrap()
            .contains("PRIVATE SAMPLE NAME")
    );
}

#[test]
fn rack_equipment_verification_rejects_tampered_frames_sources_and_actions_even_if_rehashed() {
    let description = engine(2, 3)
        .rack_equipment_description("rack 'one'", TEMPLATE)
        .unwrap();
    for change in 0..6 {
        let mut changed = description.clone();
        match change {
            0 => changed.slots[0].frame.translation_m[0] += 0.01,
            1 => changed.source.physical_template.pitch_x_mm += 1.0,
            2 => changed.permitted_actions.push("move_robot".into()),
            3 => changed.asset_revision_sha256 = "sha256:forged".into(),
            4 => changed.mass_kg = Some(0.1),
            _ => changed.length_unit = "millimetre".into(),
        }
        changed.content_sha256.clear();
        changed.content_sha256 = digest(&changed).unwrap();
        assert!(
            GentleEngine::verify_rack_equipment_description(&changed).is_err(),
            "change {change}"
        );
    }
    let roundtrip: RackEquipmentDescription =
        serde_json::from_slice(&serde_json::to_vec(&description).unwrap()).unwrap();
    GentleEngine::verify_rack_equipment_description(&roundtrip).unwrap();
    let mut unknown = serde_json::to_value(description).unwrap();
    unknown["robot_approved"] = serde_json::json!(true);
    assert!(serde_json::from_value::<RackEquipmentDescription>(unknown).is_err());
}

#[test]
fn rack_equipment_refuses_oversize_invalid_geometry_and_unknown_blocked_positions() {
    assert!(
        engine(2, 3)
            .rack_equipment_description(
                "rack 'one'",
                RackPhysicalTemplateKind::PcrTube96PositionCarrier
            )
            .is_err()
    );
    for (rows, columns) in [(0, 1), (1, 1537), (usize::MAX, 2)] {
        assert!(
            engine(rows, columns)
                .rack_equipment_description("rack 'one'", TEMPLATE)
                .is_err()
        );
    }
    let mut engine = engine(2, 3);
    engine
        .state
        .container_state
        .racks
        .get_mut("rack 'one'")
        .unwrap()
        .profile
        .blocked_coordinates
        .push("Z99".into());
    assert!(
        engine
            .rack_equipment_description("rack 'one'", TEMPLATE)
            .is_err()
    );
    assert!(
        engine
            .rack_equipment_description("absent", TEMPLATE)
            .is_err()
    );
    let baseline = super::tests::engine(2, 3)
        .rack_equipment_description("rack 'one'", TEMPLATE)
        .unwrap();
    for value in [0.0, -1.0, f32::NAN, f32::INFINITY] {
        let mut source = baseline.source.clone();
        source.physical_template.pitch_x_mm = value;
        assert!(
            GentleEngine::build_rack_equipment_description("rack", source, "synthetic").is_err()
        );
    }
}

#[test]
fn rack_equipment_export_never_overwrites_or_leaves_a_staging_file() {
    let engine = engine(2, 3);
    let temp = tempfile::tempdir().unwrap();
    let output = temp.path().join("owned by someone else.json");
    fs::write(&output, b"preserve me").unwrap();
    assert!(
        engine
            .export_rack_equipment_json("rack 'one'", TEMPLATE, output.to_str().unwrap())
            .is_err()
    );
    assert_eq!(fs::read(&output).unwrap(), b"preserve me");
    assert_eq!(fs::read_dir(temp.path()).unwrap().count(), 1);
    assert!(
        engine
            .export_rack_equipment_json(
                "absent",
                TEMPLATE,
                temp.path().join("missing.json").to_str().unwrap()
            )
            .is_err()
    );
    assert_eq!(fs::read_dir(temp.path()).unwrap().count(), 1);
}

#[test]
fn rack_equipment_shared_parser_preserves_native_paths_and_returns_canonical_output() {
    let mut engine = engine(2, 3);
    let temp = tempfile::tempdir().unwrap();
    let path = temp.path().join("equipment's description.json");
    for native in [
        path.to_str().unwrap(),
        r"C:\lab's files\rack.json",
        r"\\server\share\rack.json",
        r"\\?\C:\lab files\rack.json",
    ] {
        let command = parse_shell_line(&format!(
            "racks equipment-json {} {} --template {}",
            quote_shell_arg("rack 'one'"),
            quote_shell_arg(native),
            TEMPLATE.as_str()
        ))
        .unwrap();
        let ShellCommand::Op { payload } = command else {
            panic!("expected shared operation")
        };
        let Operation::ExportRackEquipmentJson {
            rack_id,
            path,
            template,
        } = serde_json::from_str(&payload).unwrap()
        else {
            panic!("wrong operation")
        };
        assert_eq!(rack_id, "rack 'one'");
        assert_eq!(path, native);
        assert_eq!(template, TEMPLATE);
    }
    let command = parse_shell_line(&format!(
        "racks equipment-json {} {} --template {}",
        quote_shell_arg("rack 'one'"),
        quote_shell_arg(path.to_str().unwrap()),
        TEMPLATE.as_str()
    ))
    .unwrap();
    let result = execute_shell_command(&mut engine, &command).unwrap();
    assert!(!result.state_changed);
    let text = fs::read_to_string(&path).unwrap();
    assert!(text.contains(RACK_EQUIPMENT_SCHEMA));
    let saved: serde_json::Value = serde_json::from_str(&text).unwrap();
    assert_eq!(result.output["result"]["rack_equipment"], saved);
    for invalid in [
        "racks equipment-json",
        "racks equipment-json rack output --template",
        "racks equipment-json rack output --move robot",
        "racks equipment-json '' output",
        "racks equipment-json rack ''",
        "racks equipment-json rack output --template made_up",
    ] {
        assert!(parse_shell_line(invalid).is_err(), "{invalid}");
    }
}

#[test]
fn rack_equipment_registry_and_introspection_require_external_write_approval() {
    use gentle_protocol::{CapabilityMutation, CapabilitySource, capability_registry};
    for (source, name) in [
        (CapabilitySource::GlossaryCommand, "racks equipment-json"),
        (CapabilitySource::EngineOperation, "ExportRackEquipmentJson"),
    ] {
        let row = capability_registry()
            .iter()
            .find(|row| row.source == source && row.name == name)
            .unwrap();
        assert_eq!(row.mutating, CapabilityMutation::External);
    }
    let mut engine = engine(1, 1);
    let result = execute_shell_command(
        &mut engine,
        &parse_shell_line("introspect capabilities").unwrap(),
    )
    .unwrap();
    for name in ["racks equipment-json", "ExportRackEquipmentJson"] {
        let row = result.output["capabilities"]
            .as_array()
            .unwrap()
            .iter()
            .find(|row| row["id"] == name)
            .unwrap();
        assert_eq!(row["annotation_status"], "fact_annotated");
        assert_eq!(row["mutating"], "external");
        assert_eq!(row["requires_confirmation"], true);
    }
}

#[test]
fn rack_equipment_forged_agent_flags_cannot_authorize_named_or_raw_writes() {
    use crate::agent_execution::execute_agent_plan_candidate_with_feedback;
    use crate::agent_planner::{
        AgentPlanCandidate, AgentPlanCandidateKind, AgentPlanExecutionMode, AgentPlanResult,
    };
    use crate::engine_shell::ShellExecutionOptions;
    let mut engine = engine(1, 1);
    let temp = tempfile::tempdir().unwrap();
    let path = temp.path().join("requires approval.json");
    let operation = serde_json::to_value(Operation::ExportRackEquipmentJson {
        rack_id: "rack 'one'".into(),
        path: path.to_str().unwrap().into(),
        template: TEMPLATE,
    })
    .unwrap();
    let named = format!(
        "racks equipment-json {} {}",
        quote_shell_arg("rack 'one'"),
        quote_shell_arg(path.to_str().unwrap())
    );
    for (kind, shell_command, operation) in [
        (AgentPlanCandidateKind::Shell, Some(named), None),
        (AgentPlanCandidateKind::Op, None, Some(operation)),
    ] {
        let plan = AgentPlanResult {
            schema: "gentle.agent_plan_result.v1".into(),
            assistant_message: "Synthetic safety test".into(),
            questions: vec![],
            candidates: vec![AgentPlanCandidate {
                candidate_id: "equipment".into(),
                title: "Export equipment".into(),
                kind,
                execution_mode: AgentPlanExecutionMode::Auto,
                mutating: false,
                requires_confirmation: false,
                shell_command,
                operation,
                ..Default::default()
            }],
        };
        let denied = execute_agent_plan_candidate_with_feedback(
            &mut engine,
            &plan,
            "equipment",
            false,
            &ShellExecutionOptions::default(),
        );
        assert!(denied.result.unwrap_err().contains("--confirm"));
        assert_eq!(denied.feedback.before, denied.feedback.after);
        assert!(!path.exists());
    }
}
