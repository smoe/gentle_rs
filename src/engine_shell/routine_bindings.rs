//! Explicit routine-port bindings shared by GUI setup and macro preflight.

use super::*;

pub(crate) fn routine_parameter_name<'a>(template: &str, port: &'a str) -> &'a str {
    match (template, port) {
        ("grna_anchor_window_scan", "anchor_a") => "anchor_a_pos",
        ("grna_anchor_window_scan", "anchor_b") => "anchor_b_pos",
        ("grna_practical_filter_and_oligos", "output_oligo_set_id") => "oligo_set_id",
        _ => port,
    }
}

pub(crate) fn grna_binding_readiness(
    engine: &GentleEngine,
    template: &str,
    bindings: &HashMap<String, String>,
) -> Result<(), String> {
    let required = |name: &str| {
        bindings
            .get(name)
            .map(String::as_str)
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .ok_or_else(|| format!("Required binding '{name}' is missing"))
    };
    match template {
        "grna_anchor_window_scan" | "grna_candidate_priority_scan" => {
            let seq_id = required("seq_id")?;
            if engine.sequence_kind(seq_id) != Some("dna") {
                return Err(format!(
                    "'{seq_id}' is missing or not DNA; select a DNA target"
                ));
            }
            let positive = |name: &str, default: usize| -> Result<usize, String> {
                let value = bindings
                    .get(name)
                    .map(|v| v.parse::<usize>())
                    .transpose()
                    .map_err(|_| format!("'{name}' must be a positive integer"))?
                    .unwrap_or(default);
                if value == 0 {
                    return Err(format!("'{name}' must be positive"));
                }
                Ok(value)
            };
            let length = positive("length_bp", 20)?;
            positive("step_bp", 1)?;
            positive("limit", 20000)?;
            if template == "grna_anchor_window_scan" {
                let a = required("anchor_a_pos")?
                    .parse::<usize>()
                    .map_err(|_| "Anchor A must be a zero-based boundary integer".to_string())?;
                let b = required("anchor_b_pos")?
                    .parse::<usize>()
                    .map_err(|_| "Anchor B must be a zero-based boundary integer".to_string())?;
                GentleEngine::candidate_anchor_positions(
                    &engine.state().sequences[seq_id],
                    &SequenceAnchor::Position { zero_based: a },
                    &SequenceAnchor::Position { zero_based: b },
                    length,
                )
                .map_err(|e| e.to_string())?;
            }
        }
        "grna_practical_filter_and_oligos" => {
            let id = required("guide_set_id")?;
            if !engine
                .guide_set_input_ids()
                .iter()
                .any(|candidate| candidate == id)
            {
                return Err(format!(
                    "Guide set '{id}' is missing; a generic candidate set is not a guide set"
                ));
            }
        }
        _ => {}
    }
    Ok(())
}

/// Background/import-time validation, never a menu paint operation.
#[cfg(feature = "desktop-gui")]
pub(crate) fn validate_routine_template_parameters(
    path: &str,
    name: &str,
    parameters: &[&str],
) -> Result<(), String> {
    let (templates, _) = load_cloning_pattern_templates_from_path(path)?;
    let template = templates
        .iter()
        .find(|template| template.name == name)
        .ok_or_else(|| format!("Template '{name}' is absent from '{path}'"))?;
    for parameter in parameters {
        if !template.parameters.iter().any(|p| p.name == *parameter) {
            return Err(format!("Template '{name}' has no parameter '{parameter}'"));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn engine() -> GentleEngine {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "dna".into(),
            crate::dna_sequence::DNAsequence::from_sequence(&"ACGT".repeat(12)).unwrap(),
        );
        GentleEngine::from_state(state)
    }

    #[test]
    fn anchor_bindings_use_explicit_mapping_and_shared_boundary_validation() {
        let engine = engine();
        let mut bindings = HashMap::from([("seq_id".into(), "dna".into())]);
        assert!(grna_binding_readiness(&engine, "grna_anchor_window_scan", &bindings).is_err());
        for (a, b, valid) in [
            (0, 48, true),
            (48, 0, true),
            (0, 0, false),
            (0, 19, false),
            (0, 49, false),
        ] {
            bindings.insert(
                routine_parameter_name("grna_anchor_window_scan", "anchor_a").into(),
                a.to_string(),
            );
            bindings.insert(
                routine_parameter_name("grna_anchor_window_scan", "anchor_b").into(),
                b.to_string(),
            );
            assert_eq!(
                grna_binding_readiness(&engine, "grna_anchor_window_scan", &bindings).is_ok(),
                valid,
                "{a}..{b}"
            );
        }
        assert_eq!(routine_parameter_name("unrelated", "anchor_a"), "anchor_a");
    }

    #[test]
    fn canonical_shell_anchor_preflight_executes_and_invalid_input_only_records_failure() {
        let mut engine = engine();
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("assets/cloning_patterns_catalog/crispr/guides/candidate_scans/grna_anchor_window_scan.json");
        execute_shell_command(
            &mut engine,
            &ShellCommand::MacrosTemplateImport {
                path: path.to_string_lossy().into_owned(),
            },
        )
        .unwrap();
        let command = parse_shell_line("macros template-run grna_anchor_window_scan --bind seq_id=dna --bind anchor_a_pos=0 --bind anchor_b_pos=48 --validate-only").unwrap();
        let result = execute_shell_command(&mut engine, &command).unwrap();
        assert_eq!(result.output["can_execute"], true, "{}", result.output);
        assert!(!result.state_changed);
        let invalid = parse_shell_line("macros template-run grna_anchor_window_scan --bind seq_id=dna --bind anchor_a_pos=0 --bind anchor_b_pos=49 --transactional").unwrap();
        let before = serde_json::to_value(engine.state()).unwrap();
        let operations_before = engine.operation_log().len();
        assert!(execute_shell_command(&mut engine, &invalid).is_err());
        let after = serde_json::to_value(engine.state()).unwrap();
        for key in ["sequences", "metadata", "operations", "container_state"] {
            assert_eq!(after[key], before[key], "scientific state changed: {key}");
        }
        let failed = engine.state().lineage.macro_instances.last().unwrap();
        assert_eq!(failed.status, MacroInstanceStatus::Failed);
        assert!(failed.expanded_op_ids.is_empty());
        assert_eq!(operations_before, engine.operation_log().len());
        let valid = parse_shell_line("macros template-run grna_anchor_window_scan --bind seq_id=dna --bind anchor_a_pos=0 --bind anchor_b_pos=48 --transactional").unwrap();
        assert!(
            execute_shell_command(&mut engine, &valid)
                .unwrap()
                .state_changed
        );
    }

    #[test]
    fn descriptor_lookup_is_static_and_matches_published_contracts() {
        let descriptors = annotated_introspection_capability_descriptors();
        for descriptor in &descriptors {
            let id = descriptor["id"].as_str().unwrap();
            let canonical = canonical_introspection_capability_id(id);
            let expected = descriptors
                .iter()
                .find(|row| row["id"].as_str() == Some(canonical));
            if expected.is_none() {
                assert!(capability_descriptor(id).is_none());
                continue;
            }
            let first = capability_descriptor(id).unwrap();
            assert_eq!(first, expected.unwrap());
            assert!(std::ptr::eq(first, capability_descriptor(id).unwrap()));
        }
        assert!(capability_descriptor("not-a-capability").is_none());
    }
}
