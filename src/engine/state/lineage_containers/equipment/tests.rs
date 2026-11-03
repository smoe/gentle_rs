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
        let validation = GentleEngine::validate_rack_equipment_bytes(&fs::read(&path).unwrap());
        assert!(
            validation.internally_consistent && validation.matches_builtin_catalog,
            "{validation:?}"
        );
        assert!(validation.diagnostics.is_empty());
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
    for change in 0..8 {
        let mut changed = description.clone();
        match change {
            0 => changed.slots[0].frame.translation_m[0] += 0.01,
            1 => changed.source.physical_template.pitch_x_mm += 1.0,
            2 => changed.permitted_actions.push("move_robot".into()),
            3 => changed.asset_revision_sha256 = "sha256:forged".into(),
            4 => changed.mass_kg = Some(0.1),
            5 => changed.inertia_kg_m2 = Some([0.1; 6]),
            6 => changed.mass_kg = Some(f64::NAN),
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

#[test]
fn rack_equipment_bytes_reject_unknown_nested_fields_and_lossy_f32_values() {
    let engine = engine(8, 12);
    let temp = tempfile::tempdir().unwrap();
    let path = temp.path().join("written.json");
    engine
        .export_rack_equipment_json("rack 'one'", TEMPLATE, path.to_str().unwrap())
        .unwrap();
    let mut raw: serde_json::Value = serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
    raw["source"]["tube_profiles"][0]["fit_approved"] = serde_json::json!(true);
    let bytes = serde_json::to_vec(&raw).unwrap();
    let result = GentleEngine::validate_rack_equipment_bytes(&bytes);
    assert!(!result.internally_consistent && !result.matches_builtin_catalog);
    assert_eq!(result.input_sha256, sha256_prefixed_bytes(&bytes));
    assert_eq!(
        result.diagnostics,
        ["rack_equipment_unknown_or_noncanonical_fields"]
    );
    let description = engine
        .rack_equipment_description(
            "rack 'one'",
            RackPhysicalTemplateKind::PcrTube96PositionCarrier,
        )
        .unwrap();
    let mut raw = serde_json::to_value(&description).unwrap();
    assert_eq!(description.source.physical_template.pitch_x_mm, 9.0);
    raw["source"]["physical_template"]["pitch_x_mm"] = serde_json::json!(9.0000001);
    let result = GentleEngine::validate_rack_equipment_bytes(&serde_json::to_vec(&raw).unwrap());
    assert!(!result.internally_consistent);
    assert_eq!(
        result.diagnostics,
        ["rack_equipment_unknown_or_noncanonical_fields"]
    );
}

#[test]
fn rack_equipment_fully_rebuilt_sources_can_be_consistent_without_catalog_match() {
    let original = engine(2, 3)
        .rack_equipment_description("rack 'one'", TEMPLATE)
        .unwrap();
    for change in 0..3 {
        let mut source = original.source.clone();
        match change {
            0 => source.physical_template.pitch_x_mm += 0.1,
            1 => source.tube_profiles[0].body_outer_diameter.value_mm = Some(9.9),
            _ => {
                source.physical_template.kind = RackPhysicalTemplateKind::StoragePcrTubeRack;
                source.carrier_pitch_override_um = Some(14_000);
            }
        }
        let rebuilt =
            GentleEngine::build_rack_equipment_description("rack 'one'", source, "synthetic")
                .unwrap();
        let result =
            GentleEngine::validate_rack_equipment_bytes(&serde_json::to_vec(&rebuilt).unwrap());
        assert!(result.internally_consistent, "{result:?}");
        assert!(!result.matches_builtin_catalog);
        assert_eq!(result.diagnostics, ["rack_equipment_catalog_mismatch"]);
        GentleEngine::verify_rack_equipment_description(&rebuilt).unwrap();
    }
}

#[test]
fn rack_equipment_pitch_override_is_recorded_and_old_override_sources_fail_catalog_match() {
    let mut engine = engine(2, 3);
    let original = engine
        .rack_equipment_description("rack 'one'", TEMPLATE)
        .unwrap();
    assert!(
        !serde_json::to_value(&original.source)
            .unwrap()
            .as_object()
            .unwrap()
            .contains_key("carrier_pitch_override_um")
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
    assert_eq!(compact.source.carrier_pitch_override_um, Some(14_000));
    let result =
        GentleEngine::validate_rack_equipment_bytes(&serde_json::to_vec(&compact).unwrap());
    assert!(result.internally_consistent && result.matches_builtin_catalog);
    let mut legacy_source = compact.source.clone();
    legacy_source.carrier_pitch_override_um = None;
    let legacy = GentleEngine::build_rack_equipment_description(
        "rack 'one'",
        legacy_source,
        &compact.producer_version,
    )
    .unwrap();
    assert_eq!(legacy.asset_revision_sha256, compact.asset_revision_sha256);
    assert_ne!(legacy.source_sha256, compact.source_sha256);
    assert_ne!(legacy.content_sha256, compact.content_sha256);
    let result = GentleEngine::validate_rack_equipment_bytes(&serde_json::to_vec(&legacy).unwrap());
    assert!(result.internally_consistent && !result.matches_builtin_catalog);
    assert_eq!(result.diagnostics, ["rack_equipment_catalog_mismatch"]);
}

#[test]
fn rack_equipment_validation_findings_use_fixed_codes_without_parser_or_file_content() {
    let invalid = b"PRIVATE SAMPLE NAME: not JSON, rack 'one'";
    let result = GentleEngine::validate_rack_equipment_bytes(invalid);
    assert_eq!(result.diagnostics, ["rack_equipment_json_invalid"]);
    assert_eq!(result.input_sha256, sha256_prefixed_bytes(invalid));
    assert_eq!(result.description_schema, None);
    assert_eq!(result.claimed_content_sha256, None);
    assert!(!serde_json::to_string(&result).unwrap().contains("PRIVATE"));
    let unsupported = GentleEngine::validate_rack_equipment_bytes(br#"{"schema":"unknown"}"#);
    assert_eq!(
        unsupported.diagnostics,
        ["rack_equipment_schema_unsupported"]
    );
    let malformed = GentleEngine::validate_rack_equipment_bytes(
        format!(r#"{{"schema":"{RACK_EQUIPMENT_SCHEMA}"}}"#).as_bytes(),
    );
    assert_eq!(malformed.diagnostics, ["rack_equipment_record_invalid"]);
    let original = engine(2, 3)
        .rack_equipment_description("rack 'one'", TEMPLATE)
        .unwrap();
    let mut raw = serde_json::to_value(original).unwrap();
    raw["content_sha256"] = serde_json::json!("sha256:wrong");
    let mismatch = GentleEngine::validate_rack_equipment_bytes(&serde_json::to_vec(&raw).unwrap());
    assert_eq!(mismatch.diagnostics, ["rack_equipment_binding_mismatch"]);
    raw["source"]["physical_template"]["pitch_x_mm"] = serde_json::json!(0.0);
    let invalid_source =
        GentleEngine::validate_rack_equipment_bytes(&serde_json::to_vec(&raw).unwrap());
    assert_eq!(
        invalid_source.diagnostics,
        ["rack_equipment_source_invalid"]
    );
}

#[test]
fn rack_equipment_validate_shell_is_read_only_bounded_and_content_findings_are_data() {
    let mut engine = engine(2, 3);
    let temp = tempfile::tempdir().unwrap();
    let path = temp.path().join("equipment's report.json");
    let description = engine
        .export_rack_equipment_json("rack 'one'", TEMPLATE, path.to_str().unwrap())
        .unwrap();
    let state = serde_json::to_value(engine.state()).unwrap();
    let history = engine.history_summary();
    let command = parse_shell_line(&format!(
        "racks equipment-validate {}",
        quote_shell_arg(path.to_str().unwrap())
    ))
    .unwrap();
    assert!(!command.is_state_mutating());
    let ShellCommand::Op { payload } = &command else {
        panic!("expected generic operation")
    };
    assert!(matches!(
        serde_json::from_str::<Operation>(payload).unwrap(),
        Operation::ValidateRackEquipmentJson { .. }
    ));
    let output = execute_shell_command(&mut engine, &command).unwrap();
    assert!(!output.state_changed);
    let result: RackEquipmentValidation =
        serde_json::from_value(output.output["result"]["rack_equipment_validation"].clone())
            .unwrap();
    assert!(result.internally_consistent && result.matches_builtin_catalog);
    assert_eq!(
        result.claimed_content_sha256.as_deref(),
        Some(description.content_sha256.as_str())
    );
    let serialized = output.output.to_string();
    assert!(!serialized.contains("rack 'one'") && !serialized.contains("PRIVATE SAMPLE NAME"));
    assert_eq!(serde_json::to_value(engine.state()).unwrap(), state);
    assert_eq!(engine.history_summary().undo_count, history.undo_count);
    assert_eq!(engine.history_summary().redo_count, history.redo_count);
    for invalid in [
        "racks equipment-validate",
        "racks equipment-validate one two",
        "racks equipment-validate file --extra",
        "racks equipment-validate ''",
    ] {
        assert!(parse_shell_line(invalid).is_err(), "{invalid}");
    }
    assert!(
        engine
            .validate_rack_equipment_json(temp.path().to_str().unwrap())
            .unwrap_err()
            .message
            .contains("not_regular_file")
    );
    let oversized = temp.path().join("large.json");
    fs::File::create(&oversized)
        .unwrap()
        .set_len(VALIDATION_INPUT_LIMIT + 1)
        .unwrap();
    assert_eq!(
        engine
            .validate_rack_equipment_json(oversized.to_str().unwrap())
            .unwrap_err()
            .message,
        "rack_equipment_validate_input_too_large"
    );
    fs::write(&path, b"PRIVATE CONTENT invalid JSON").unwrap();
    let output = execute_shell_command(&mut engine, &command).unwrap();
    assert_eq!(
        output.output["result"]["rack_equipment_validation"]["diagnostics"],
        serde_json::json!(["rack_equipment_json_invalid"])
    );
    assert!(!output.output.to_string().contains("PRIVATE CONTENT"));
}

#[test]
fn rack_equipment_validate_agent_runs_without_additional_confirmation() {
    use crate::agent_execution::execute_agent_plan_candidate_with_feedback;
    use crate::agent_planner::{
        AgentPlanCandidate, AgentPlanCandidateKind, AgentPlanExecutionMode, AgentPlanResult,
    };
    let mut engine = engine(2, 3);
    let temp = tempfile::tempdir().unwrap();
    let path = temp.path().join("validation input.json");
    engine
        .export_rack_equipment_json("rack 'one'", TEMPLATE, path.to_str().unwrap())
        .unwrap();
    let before = serde_json::to_value(engine.state()).unwrap();
    let plan = AgentPlanResult {
        schema: "gentle.agent_plan_result.v1".into(),
        assistant_message: String::new(),
        questions: vec![],
        candidates: vec![AgentPlanCandidate {
            candidate_id: "validate".into(),
            title: "Check consistency".into(),
            kind: AgentPlanCandidateKind::Shell,
            execution_mode: AgentPlanExecutionMode::Auto,
            mutating: false,
            requires_confirmation: false,
            shell_command: Some(format!(
                "racks equipment-validate {}",
                quote_shell_arg(path.to_str().unwrap())
            )),
            ..Default::default()
        }],
    };
    let attempt = execute_agent_plan_candidate_with_feedback(
        &mut engine,
        &plan,
        "validate",
        false,
        &Default::default(),
    );
    let result = attempt.result.unwrap();
    assert!(!result.state_changed);
    assert!(
        result
            .output
            .to_string()
            .contains(RACK_EQUIPMENT_VALIDATION_SCHEMA)
    );
    assert_eq!(serde_json::to_value(engine.state()).unwrap(), before);
}

#[test]
fn rack_equipment_export_reuses_project_and_sidecar_destination_protection() {
    let temp = tempfile::tempdir().unwrap();
    let project = temp.path().join("project.json");
    let mut engine = engine(2, 3);
    engine.set_host_project_path(Some(&project)).unwrap();
    let sidecar = ProjectState::candidate_store_sidecar_dir(&project);
    fs::create_dir(&sidecar).unwrap();
    let nested = sidecar.join("records");
    fs::create_dir(&nested).unwrap();
    let sentinel = nested.join("keep.json");
    fs::write(&sentinel, b"keep unchanged").unwrap();
    for path in [
        project.clone(),
        temp.path().join("PROJECT.JSON"),
        sidecar.join("new.json"),
        nested.join("new.json"),
    ] {
        let error = engine
            .export_rack_equipment_json("rack 'one'", TEMPLATE, path.to_str().unwrap())
            .unwrap_err();
        assert_eq!(error.code, ErrorCode::InvalidInput);
        assert!(
            error
                .message
                .starts_with("rack_equipment_export.protected_path:")
        );
        assert!(!path.exists());
    }
    let raw_parent_leaf = format!("{}{}..", temp.path().display(), std::path::MAIN_SEPARATOR);
    let trailing_separator = format!("{}{}", temp.path().display(), std::path::MAIN_SEPARATOR);
    for raw in [
        raw_parent_leaf,
        trailing_separator,
        "control\nname.json".into(),
    ] {
        let error = engine
            .export_rack_equipment_json("rack 'one'", TEMPLATE, &raw)
            .unwrap_err();
        assert_eq!(error.code, ErrorCode::InvalidInput);
        assert!(
            error
                .message
                .starts_with("rack_equipment_export.invalid_path:")
        );
    }
    assert_eq!(fs::read(sentinel).unwrap(), b"keep unchanged");
}

#[test]
fn rack_equipment_pcr96_characterizes_current_row_a_label_edge_convention() {
    let description = engine(8, 12)
        .rack_equipment_description(
            "rack 'one'",
            RackPhysicalTemplateKind::PcrTube96PositionCarrier,
        )
        .unwrap();
    for (coordinate, expected) in [
        ("A1", [0.01438, 0.01124]),
        ("A12", [0.11338, 0.01124]),
        ("H1", [0.01438, 0.07424]),
    ] {
        let frame = &description
            .slots
            .iter()
            .find(|slot| slot.coordinate == coordinate)
            .unwrap()
            .frame;
        for axis in 0..2 {
            assert!(
                (frame.translation_m[axis] - expected[axis]).abs() < 1e-6,
                "{coordinate}, axis {axis}"
            );
        }
        assert_eq!(frame.id, format!("slot_{coordinate}"));
        assert_eq!(frame.parent.as_deref(), Some("rack_base"));
        assert_eq!(frame.rotation_xyzw, [0.0, 0.0, 0.0, 1.0]);
    }
    assert_eq!(
        description.coordinate_convention,
        "right_handed; origin=front_left_base; +x=increasing_columns; +y=increasing_rows; +z=up; slot_origin=carrier_top; quaternion=xyzw"
    );
}

#[test]
fn rack_equipment_validation_registry_and_introspection_are_read_only() {
    use gentle_protocol::{CapabilityMutation, CapabilitySource, capability_registry};
    for (source, name) in [
        (
            CapabilitySource::GlossaryCommand,
            "racks equipment-validate",
        ),
        (
            CapabilitySource::EngineOperation,
            "ValidateRackEquipmentJson",
        ),
    ] {
        let row = capability_registry()
            .iter()
            .find(|row| row.source == source && row.name == name)
            .unwrap();
        assert_eq!(row.mutating, CapabilityMutation::ReadOnly);
    }
    let mut engine = engine(1, 1);
    let output = execute_shell_command(
        &mut engine,
        &parse_shell_line("introspect capabilities").unwrap(),
    )
    .unwrap();
    for name in ["racks equipment-validate", "ValidateRackEquipmentJson"] {
        let row = output.output["capabilities"]
            .as_array()
            .unwrap()
            .iter()
            .find(|row| row["id"] == name)
            .unwrap();
        assert_eq!(row["annotation_status"], "fact_annotated");
        assert_eq!(row["mutating"], "false");
        assert_eq!(row["requires_confirmation"], false);
        assert!(row["effects"].as_array().unwrap().is_empty());
    }
}

#[test]
fn rack_equipment_projection_pinned_hashes_for_all_seven_templates() {
    // Changes require a projection-version bump or a deliberate re-pin recorded
    // in CHANGELOG. These are synthetic code literals, not checkout-bound files.
    let engine = engine(8, 12);
    let templates = [
        RackPhysicalTemplateKind::StoragePcrTubeRack,
        RackPhysicalTemplateKind::PipettingPcrTubeRack,
        RackPhysicalTemplateKind::CellCulturePlate,
        RackPhysicalTemplateKind::SafeLockTubeRackCollarSeated,
        TEMPLATE,
        RackPhysicalTemplateKind::PcrTube96PositionCarrier,
        RackPhysicalTemplateKind::SafeLockStorageBoxDerived,
    ];
    let mut actual = Vec::new();
    for template in templates {
        let source = engine
            .rack_equipment_description("rack 'one'", template)
            .unwrap()
            .source;
        let pinned =
            GentleEngine::build_rack_equipment_description("rack 'one'", source, "pinned").unwrap();
        actual.push([
            pinned.asset_revision_sha256,
            pinned.source_sha256,
            pinned.content_sha256,
        ]);
    }
    let expected = [
        [
            "sha256:50e352108a02934063be06ff7c6da257a7dc9a464f1b6ad497e6976929f30f53",
            "sha256:192402cc56fe854ea775f7fc93258497962b46517ff07b5706c7a66b4c3b3d00",
            "sha256:b5753275fb7eea637b3bb5705c8ef488ddf69a2fa4a1d37a3ec657c0359efa41",
        ],
        [
            "sha256:c16bffdc91173d98ce66a00e1a5cfff2d83e7b6c71fc4e5528e5bcc90c47cb5d",
            "sha256:647ec6942f183ad4e9f60cc26b469496f9a4e784e9301d1895fd20c03df35f68",
            "sha256:4372cfc0227aa6bf39ce236a4359aba534bacfab785ec82d630f4bd1cebf982c",
        ],
        [
            "sha256:062f28ffb1d7ed7ffbedaf2902323465b738ea73ae6f8d292761b746f12dcec6",
            "sha256:73882c1100d4d6e3ef2854a7b7128315e2acf3f195e89ccb0b03ed9ea813c5f0",
            "sha256:a8642ea4594c286c2eb0c6ae7194c8a6b85a404f772885974dfca4530a0eb636",
        ],
        [
            "sha256:fd8ff448cb23dbb309b73c483ed7d9fe7bfe52aeb6d3a2137e775cf12c2d87e8",
            "sha256:70016ad781fe8806c9e9871fbcf45f8f20c9ecf83a6ef8722a084aea6d3acd7e",
            "sha256:0586e36ee868f1c0cbbf534ed5ae8c22d43c3e74fdb2eef766c5c91afac0c659",
        ],
        [
            "sha256:212e7bf803f65dc23fdcc03a0afa1c2bf90399b09a9631de22fa72abee15be25",
            "sha256:5aebee8cbd2c738d21347332400752ed200290a6b9fef7fe523fb830f3f7acde",
            "sha256:de4ecb0be8f3e9f3d23a83a0c327cea7ebb2fed1d78f8c81f19c09ad44e93dc7",
        ],
        [
            "sha256:a4a249cda9664b206e5898d664ac672d54e00a53a267f0ace2ca330f1dc32e09",
            "sha256:c120cd8dd7553a36a62e350d389337d7afcc6ccad47fb57f8693f7e8e4de77ab",
            "sha256:208874fa12e71cd6a67c0647e729ee9a433cb1d2058e7ae32ba9664eb871273e",
        ],
        [
            "sha256:01af0cf6c062429fae8b2ddd88135d2f8fbab49eef464f4624507259fd2d0513",
            "sha256:b9bf856bf802a2b7542084e8c6205873dd427843ee7d6f698152f61d623bde0c",
            "sha256:845a8b89dcb5916f8102bf02e5d3b4147bf2c2d7d38c432dbacb4863ff287f6e",
        ],
    ]
    .map(|row| row.map(String::from));
    assert_eq!(actual.as_slice(), expected.as_slice());
}
