//! Synthetic planning-arithmetic regressions; no real protocol or laboratory data.

use super::*;

fn planning_wetlab_test_routine(id: &str) -> CloningRoutineDefinition {
    CloningRoutineDefinition {
        routine_id: id.into(),
        family: "restriction".into(),
        base_time_hours: Some(6.0),
        ..Default::default()
    }
}

fn planning_wetlab_test_estimate(routine: &CloningRoutineDefinition) -> PlanningEstimate {
    estimate_routine_planning(
        &GentleEngine::from_state(ProjectState::default()),
        &PlanningObjective::default(),
        routine,
        &RoutinePreferenceContext::default(),
    )
    .unwrap()
}

fn planning_wetlab_test_catalog() -> serde_json::Value {
    serde_json::from_str(include_str!("../../assets/cloning_wetlab_steps.json")).unwrap()
}

fn planning_wetlab_test_parse(
    value: &Value,
) -> Result<planning_wetlab::Catalog, planning_wetlab::CatalogError> {
    let known = ["restriction.digest_ligate_extract_sticky".to_string()]
        .into_iter()
        .collect();
    planning_wetlab::Catalog::parse(&value.to_string(), &known)
}

#[test]
fn planning_wetlab_incubation_is_not_human_work() {
    let mut value = planning_wetlab_test_catalog();
    let mut step = value["routines"][0]["steps"][1].clone();
    step["depends_on"] = json!([]);
    value["routines"][0]["steps"] = json!([step]);
    let catalog = planning_wetlab_test_parse(&value).unwrap();
    let model = catalog
        .model("restriction.digest_ligate_extract_sticky")
        .unwrap();
    let schedule = model.schedule(0.0).unwrap();
    assert_eq!(schedule.hands_on.typical, 0.0);
    assert_eq!(schedule.unattended.typical, 2.0);
    assert_eq!(schedule.elapsed.typical, 2.0);
    let delayed = model.schedule(10.0).unwrap();
    assert_eq!(delayed.hands_on, schedule.hands_on);
    assert_eq!(delayed.unattended, schedule.unattended);
    assert_eq!(delayed.elapsed.typical, 12.0);
}

#[test]
fn planning_wetlab_bundled_parallel_branch_is_not_summed_into_elapsed() {
    let estimate = planning_wetlab_test_estimate(&planning_wetlab_test_routine(
        "restriction.digest_ligate_extract_sticky",
    ));
    assert_eq!(estimate.schema, "gentle.planning_estimate.v2");
    assert_eq!(
        estimate.step_model,
        crate::engine::PlanningStepModel::Present
    );
    assert_eq!(
        estimate.estimated_time_hours_basis,
        crate::engine::PlanningTimeBasis::StepModelCriticalPath
    );
    let hands = estimate.hands_on_hours.unwrap();
    let unattended = estimate.unattended_hours.unwrap();
    let elapsed = estimate.elapsed_hours.unwrap();
    assert!(elapsed.typical < hands.typical + unattended.typical);
    assert!((elapsed.typical - 3.2).abs() < 1e-9);
    assert_eq!(estimate.estimated_time_hours, elapsed.typical);
    assert!(
        estimate
            .schedule
            .iter()
            .filter(|step| step.depends_on.is_empty())
            .count()
            >= 2
    );
    let report = serde_json::to_value(&estimate).unwrap();
    for row in std::iter::once(&report).chain(report["schedule"].as_array().unwrap().iter()) {
        for field in [
            "hands_on_hours",
            "unattended_hours",
            "elapsed_hours",
            "earliest_start_hours",
            "earliest_finish_hours",
        ] {
            if let Some(range) = row.get(field).filter(|v| !v.is_null()) {
                assert!(range["min"].as_f64().unwrap() <= range["typical"].as_f64().unwrap());
                assert!(range["typical"].as_f64().unwrap() <= range["max"].as_f64().unwrap());
            }
        }
    }
}

#[test]
fn planning_wetlab_range_scenarios_choose_independent_critical_paths() {
    // Hand-crafted graph: the longer branch changes between range scenarios.
    let mut value = planning_wetlab_test_catalog();
    let mut left = value["routines"][0]["steps"][0].clone();
    left["step_id"] = json!("alpha");
    left["hands_on_hours"] = json!({"min":1.0,"typical":5.0,"max":5.0});
    let mut right = left.clone();
    right["step_id"] = json!("beta");
    right["hands_on_hours"] = json!({"min":2.0,"typical":3.0,"max":10.0});
    value["routines"][0]["steps"] = json!([right, left]);
    let first = planning_wetlab_test_parse(&value).unwrap().routines[0]
        .schedule(0.0)
        .unwrap();
    assert_eq!(
        first.elapsed,
        crate::engine::TimeRange {
            min: 2.0,
            typical: 5.0,
            max: 10.0
        }
    );
    assert_eq!(
        first.critical_paths,
        [
            vec!["beta".to_string()],
            vec!["alpha".into()],
            vec!["beta".into()]
        ]
    );
    value["routines"][0]["steps"]
        .as_array_mut()
        .unwrap()
        .reverse();
    let reordered = planning_wetlab_test_parse(&value).unwrap().routines[0]
        .schedule(0.0)
        .unwrap();
    assert_eq!(first.rows, reordered.rows);
    assert_eq!(first.critical_paths, reordered.critical_paths);
}

#[test]
fn planning_wetlab_cost_components_reconcile_without_changing_machine_basis() {
    for id in [
        "restriction.digest_ligate_extract_sticky",
        "synthetic.legacy_routine",
    ] {
        let mut engine = GentleEngine::from_state(ProjectState::default());
        execute_shell_command(
            &mut engine,
            &ShellCommand::PlanningProfileSet {
                scope: PlanningProfileScope::ProjectOverride,
                payload_json: json!({"schema":"gentle.planning_profile.v1", "currency":"EUR",
                    "inventory":{"stock":{"available":true,"unit_cost":2.0},
                                 "missing":{"available":false,"unit_cost":7.0},
                                 "unpriced":{"available":true}},
                    "machine_availability":{"machine":{"available":true,"run_cost_per_hour":5.0}}
                })
                .to_string(),
            },
        )
        .unwrap();
        let mut routine = planning_wetlab_test_routine(id);
        routine.base_cost = Some(3.0);
        routine.required_material_classes = vec![
            "stock".into(),
            "missing".into(),
            "unknown".into(),
            "unpriced".into(),
        ];
        routine.required_machine_classes = vec!["machine".into()];
        let estimate = estimate_routine_planning(
            &engine,
            &PlanningObjective::default(),
            &routine,
            &RoutinePreferenceContext::default(),
        )
        .unwrap();
        assert_eq!(estimate.stock_drawn_material_cost, 2.0);
        assert_eq!(estimate.missing_material_purchase_proxy_cost, 7.0);
        assert_eq!(estimate.machine_run_cost, 30.0);
        assert_eq!(estimate.routine_base_cost, 3.0);
        assert_eq!(estimate.estimated_cost, 42.0);
        assert!(
            (estimate.estimated_cost
                - (estimate.routine_base_cost
                    + estimate.stock_drawn_material_cost
                    + estimate.missing_material_purchase_proxy_cost
                    + estimate.machine_run_cost))
                .abs()
                < 1e-9
        );
        assert_eq!(estimate.explanation["machine_run_cost_basis_hours"], 6.0);
        assert_eq!(
            estimate.explanation["unknown_material_classes"],
            json!(["unknown"])
        );
        assert_eq!(
            estimate.explanation["unpriced_material_classes"],
            json!(["unpriced"])
        );
        assert_eq!(
            estimate.explanation["stock_drawn_material_classes"],
            json!(["stock", "unpriced"])
        );
        assert_eq!(
            estimate.explanation["missing_material_proxy_classes"],
            json!(["missing"])
        );
        assert_eq!(
            estimate.explanation["cost_components"]["estimated_cost"],
            42.0
        );
        assert!(!estimate.cost_assumptions.is_empty());
    }
}

#[test]
fn planning_wetlab_procurement_delay_changes_elapsed_not_hands_on() {
    let mut engine = GentleEngine::from_state(ProjectState::default());
    let mut routine = planning_wetlab_test_routine("restriction.digest_ligate_extract_sticky");
    routine.required_material_classes = vec!["mix".into()];
    let before = estimate_routine_planning(
        &engine,
        &PlanningObjective::default(),
        &routine,
        &RoutinePreferenceContext::default(),
    )
    .unwrap();
    execute_shell_command(
        &mut engine,
        &ShellCommand::PlanningProfileSet {
            scope: PlanningProfileScope::ProjectOverride,
            payload_json: json!({"schema":"gentle.planning_profile.v1",
                "inventory":{"mix":{"available":false,"procurement_business_days":2.0}},
                "machine_availability":{"thermocycler":{"available":true,"queue_business_days":0.5}}
            })
            .to_string(),
        },
    )
    .unwrap();
    let after = estimate_routine_planning(
        &engine,
        &PlanningObjective::default(),
        &routine,
        &RoutinePreferenceContext::default(),
    )
    .unwrap();
    assert_eq!(after.hands_on_hours, before.hands_on_hours);
    assert_eq!(after.unattended_hours, before.unattended_hours);
    let prefix = planning_business_days_to_elapsed_hours(2.5);
    assert!(
        (after.elapsed_hours.unwrap().typical - before.elapsed_hours.unwrap().typical - prefix)
            .abs()
            < 1e-9
    );
    assert!(
        after
            .schedule
            .iter()
            .all(|step| step.earliest_start_hours.min + 1e-9 >= prefix)
    );
}

#[test]
fn planning_wetlab_unmodeled_ranges_remain_unknown_and_legacy_scalar_unchanged() {
    let mut routine = planning_wetlab_test_routine("synthetic.legacy_routine");
    routine.base_time_hours = None;
    let estimate = planning_wetlab_test_estimate(&routine);
    assert_eq!(
        estimate.step_model,
        crate::engine::PlanningStepModel::Absent
    );
    assert_eq!(
        estimate.estimated_time_hours_basis,
        crate::engine::PlanningTimeBasis::LegacyFamilyScalar
    );
    assert_eq!(estimate.estimated_time_hours, 6.0);
    assert!(estimate.hands_on_hours.is_none());
    assert!(estimate.unattended_hours.is_none());
    assert!(estimate.elapsed_hours.is_none());
    assert!(estimate.schedule.is_empty());
    let value = serde_json::to_value(&estimate).unwrap();
    for field in ["hands_on_hours", "unattended_hours", "elapsed_hours"] {
        assert!(value.get(field).is_none_or(Value::is_null));
    }
    let text = planning_estimate_display_lines(&estimate).join("\n");
    assert!(text.contains("Hands-on: unknown"));
    assert!(text.contains("legacy family scalar"));
    let mut zero = estimate.clone();
    zero.hands_on_hours = Some(crate::engine::TimeRange::default());
    let restored: PlanningEstimate =
        serde_json::from_value(serde_json::to_value(zero).unwrap()).unwrap();
    assert_eq!(restored.hands_on_hours.unwrap().typical, 0.0);
    let restored: PlanningEstimate = serde_json::from_value(value).unwrap();
    assert!(restored.hands_on_hours.is_none());
    let mut legacy = serde_json::to_value(&estimate).unwrap();
    legacy["schema"] = json!("gentle.planning_estimate.v1");
    for field in [
        "hands_on_hours",
        "unattended_hours",
        "elapsed_hours",
        "step_model",
        "estimated_time_hours_basis",
        "schedule",
        "stock_drawn_material_cost",
        "missing_material_purchase_proxy_cost",
        "machine_run_cost",
        "routine_base_cost",
        "cost_assumptions",
    ] {
        legacy.as_object_mut().unwrap().remove(field);
    }
    let restored: PlanningEstimate = serde_json::from_value(legacy).unwrap();
    let lines = planning_estimate_display_lines(&restored).join("\n");
    assert!(lines.contains("Hands-on: unknown"));
    assert!(lines.contains("Cost components: unknown"));
}

#[test]
fn planning_wetlab_all_bundled_routines_reconcile_and_only_one_has_a_model() {
    let catalog: CloningRoutineCatalog =
        serde_json::from_str(include_str!("../../assets/cloning_routines.json")).unwrap();
    assert_eq!(catalog.routines.len(), 22);
    let mut profile = json!({"schema":"gentle.planning_profile.v1",
        "inventory":{}, "machine_availability":{}});
    for routine in &catalog.routines {
        for (index, class) in infer_routine_material_classes(routine).iter().enumerate() {
            profile["inventory"][class] = json!({"available":index % 2 == 0,
                "unit_cost":(index + 1) as f64 / 10.0});
        }
        for (index, class) in infer_routine_machine_classes(routine).iter().enumerate() {
            profile["machine_availability"][class] = json!({"available":true,
                "run_cost_per_hour":(index + 1) as f64 / 7.0});
        }
    }
    let mut engine = GentleEngine::from_state(ProjectState::default());
    execute_shell_command(
        &mut engine,
        &ShellCommand::PlanningProfileSet {
            scope: PlanningProfileScope::ProjectOverride,
            payload_json: profile.to_string(),
        },
    )
    .unwrap();
    let profile = engine.planning_effective_profile();
    let mut unmodeled = 0;
    for routine in &catalog.routines {
        let estimate = estimate_routine_planning(
            &engine,
            &PlanningObjective::default(),
            routine,
            &RoutinePreferenceContext::default(),
        )
        .unwrap();
        // Independent replay of the pre-v2 total's addition order, with decimals.
        let mut legacy_cost = default_routine_base_cost(routine);
        let mut legacy_elapsed = default_routine_base_time_hours(routine);
        let mut purchase_days = 0.0;
        for class in infer_routine_material_classes(routine) {
            let item = &profile.inventory[&class];
            legacy_cost += item.unit_cost.unwrap();
            if !item.available {
                purchase_days += profile.procurement_business_days_default;
            }
        }
        for class in infer_routine_machine_classes(routine) {
            let machine = &profile.machine_availability[&class];
            legacy_cost +=
                machine.run_cost_per_hour.unwrap() * default_routine_base_time_hours(routine);
            legacy_elapsed += planning_business_days_to_elapsed_hours(machine.queue_business_days);
        }
        legacy_elapsed += planning_business_days_to_elapsed_hours(purchase_days);
        assert_eq!(estimate.estimated_cost.to_bits(), legacy_cost.to_bits());
        assert!(
            (estimate.estimated_cost
                - (estimate.routine_base_cost
                    + estimate.stock_drawn_material_cost
                    + estimate.missing_material_purchase_proxy_cost
                    + estimate.machine_run_cost))
                .abs()
                < 1e-9
        );
        if routine.routine_id != "restriction.digest_ligate_extract_sticky" {
            unmodeled += 1;
            assert_eq!(
                estimate.step_model,
                crate::engine::PlanningStepModel::Absent
            );
            assert!(estimate.hands_on_hours.is_none());
            assert!(estimate.unattended_hours.is_none());
            assert!(estimate.elapsed_hours.is_none());
            assert!(estimate.schedule.is_empty());
            assert_eq!(
                estimate.estimated_time_hours.to_bits(),
                legacy_elapsed.to_bits()
            );
        }
    }
    assert_eq!(unmodeled, 21);
}

#[test]
fn planning_wetlab_missing_machine_price_remains_unknown_not_free() {
    let estimate = planning_wetlab_test_estimate(&planning_wetlab_test_routine(
        "restriction.digest_ligate_extract_sticky",
    ));
    assert_eq!(estimate.machine_run_cost, 0.0);
    assert!(
        !estimate.explanation["unpriced_machine_classes"]
            .as_array()
            .unwrap()
            .is_empty()
    );
    assert!(
        estimate.explanation["machine_cost_rows"]
            .as_array()
            .unwrap()
            .iter()
            .all(|row| row["cost"].is_null())
    );
    assert_eq!(
        estimate.explanation["cost_coverage"],
        "incomplete_known_unit_subtotals"
    );
}

#[test]
fn planning_wetlab_catalog_rejects_unknown_routines_cycles_and_invalid_ranges() {
    let mut value = planning_wetlab_test_catalog();
    value["routines"][0]["routine_id"] = json!("unknown.routine");
    let err = planning_wetlab_test_parse(&value).unwrap_err();
    assert_eq!(err.code, "unknown_routine_id");
    assert_eq!(
        serde_json::from_str::<Value>(&err.to_string()).unwrap()["code"],
        "unknown_routine_id"
    );
    let mut value = planning_wetlab_test_catalog();
    value["routines"][0]["steps"][0]["depends_on"] = json!(["review_recovery"]);
    assert_eq!(
        planning_wetlab_test_parse(&value).unwrap_err().code,
        "wetlab_dependency_cycle"
    );
    for range in [
        json!({"min":-1.0,"typical":0.0,"max":1.0}),
        json!({"min":2.0,"typical":1.0,"max":3.0}),
        json!({"min":0.0,"typical":3.0,"max":2.0}),
    ] {
        let mut value = planning_wetlab_test_catalog();
        value["routines"][0]["steps"][0]["hands_on_hours"] = range;
        assert_eq!(
            planning_wetlab_test_parse(&value).unwrap_err().code,
            "invalid_wetlab_duration"
        );
    }
    let mut value = planning_wetlab_test_catalog();
    value["routines"][0]["steps"][0]["depends_on"] = json!(["missing"]);
    assert_eq!(
        planning_wetlab_test_parse(&value).unwrap_err().code,
        "unknown_wetlab_dependency"
    );
    let mut value = planning_wetlab_test_catalog();
    value["routines"][0]["steps"][0]["hands_on_hours"]
        .as_object_mut()
        .unwrap()
        .remove("min");
    assert_eq!(
        planning_wetlab_test_parse(&value).unwrap_err().code,
        "invalid_wetlab_catalog"
    );
    let mut value = planning_wetlab_test_catalog();
    let duplicate = value["routines"][0].clone();
    value["routines"].as_array_mut().unwrap().push(duplicate);
    assert_eq!(
        planning_wetlab_test_parse(&value).unwrap_err().code,
        "duplicate_wetlab_routine"
    );
    let mut value = planning_wetlab_test_catalog();
    let duplicate = value["routines"][0]["steps"][0].clone();
    value["routines"][0]["steps"]
        .as_array_mut()
        .unwrap()
        .push(duplicate);
    assert_eq!(
        planning_wetlab_test_parse(&value).unwrap_err().code,
        "duplicate_wetlab_step"
    );
}

#[test]
fn planning_wetlab_shared_consultation_and_comparison_remain_read_only() {
    let mut engine = GentleEngine::from_state(ProjectState::default());
    let before = serde_json::to_value(engine.state()).unwrap();
    let report = execute_shell_command(
        &mut engine,
        &parse_shell_line("planning consult cloning --format text").unwrap(),
    )
    .unwrap();
    assert!(!report.state_changed);
    assert_eq!(
        report.output["strategy_candidates"]
            .as_array()
            .unwrap()
            .len(),
        11
    );
    let text = report.output["text_report"].as_str().unwrap();
    for label in [
        "Hands-on:",
        "Elapsed:",
        "Stock-drawn material",
        "Missing-material purchase proxy",
        "Machine run",
        "legacy family scalar",
    ] {
        assert!(text.contains(label), "Missing {label}");
    }
    let comparison = execute_shell_command(&mut engine, &parse_shell_line("routines compare restriction.digest_ligate_extract_sticky gibson.two_fragment_overlap_preview").unwrap()).unwrap();
    assert!(!comparison.state_changed);
    let axes = comparison.output["comparison"]["difference_matrix"]
        .as_array()
        .unwrap();
    let basis = axes
        .iter()
        .find(|row| row["axis"] == "estimated_time_hours_basis")
        .unwrap();
    assert_eq!(basis["same"], false);
    let hands = axes
        .iter()
        .find(|row| row["axis"] == "hands_on_hours")
        .unwrap();
    assert_eq!(hands["right"], "unknown");
    assert_eq!(hands["comparable"], false);
    assert_eq!(serde_json::to_value(engine.state()).unwrap(), before);
}
