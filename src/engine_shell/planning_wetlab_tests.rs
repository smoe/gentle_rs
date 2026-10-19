//! Synthetic planning-arithmetic regressions; no real protocol or laboratory data.

use super::*;

type PlanningRankingSnapshotRow = (usize, &'static str, &'static str, f64);

// Hand-authored policy snapshots of the bundled catalog/default objective, not
// laboratory data or method recommendations. Changing these literals requires
// an owner decision and a docs/CHANGELOG.md entry, not automatic re-blessing.
// Seven default rows form two exact-score ties: their order is alphabetical,
// not a biological preference. Six distinct scores = two ties + four singletons.
const PLANNING_RANKING_DEFAULT_SNAPSHOT: [PlanningRankingSnapshotRow; 11] = [
    (1, "sequence", "sequence.branch_reverse_complement", -0.15),
    (2, "crispr", "crispr.guides.anchor_window_scan", -3.0),
    (3, "golden_gate", "golden_gate.type_iis_multi_insert", -3.0),
    (4, "pcr", "pcr.site_insertion_then_digest", -3.0),
    (5, "ta_gc", "ta_gc.gc_single_insert", -3.0),
    (6, "topo", "topo.blunt_single_insert", -3.0),
    (7, "infusion", "infusion.multi_fragment_overlap", -3.15),
    (
        8,
        "nebuilder_hifi",
        "nebuilder_hifi.multi_fragment_overlap",
        -3.15,
    ),
    (9, "gibson", "gibson.two_fragment_overlap_preview", -4.15),
    (
        10,
        "restriction",
        "restriction.digest_ligate_extract_sticky",
        -5.3,
    ),
    (11, "gateway", "gateway.bp_single_insert", -7.15),
];

// Independent literals for docs/examples/planning/cloning_waiting.json. Its
// ligase and equipment allowances change both family order and the Golden Gate
// winner. Do not compute this expectation from the production score/sort.
const PLANNING_RANKING_WAITING_SNAPSHOT: [PlanningRankingSnapshotRow; 11] = [
    (1, "sequence", "sequence.branch_reverse_complement", -0.15),
    (2, "crispr", "crispr.guides.anchor_window_scan", -3.0),
    (3, "ta_gc", "ta_gc.gc_single_insert", -3.0),
    (4, "topo", "topo.blunt_single_insert", -3.0),
    (5, "pcr", "pcr.site_insertion_then_digest", -5.0),
    (
        6,
        "golden_gate",
        "golden_gate.type_iis_destination_cleanup",
        -5.15,
    ),
    (7, "infusion", "infusion.multi_fragment_overlap", -39.8),
    (
        8,
        "nebuilder_hifi",
        "nebuilder_hifi.multi_fragment_overlap",
        -39.8,
    ),
    (9, "gibson", "gibson.two_fragment_overlap_preview", -47.8),
    (10, "gateway", "gateway.bp_single_insert", -63.8),
    (
        11,
        "restriction",
        "restriction.digest_ligate_extract_sticky",
        -128.35,
    ),
];

fn planning_ranking_test_engine(waiting: bool) -> GentleEngine {
    let mut engine = GentleEngine::from_state(ProjectState::default());
    if waiting {
        execute_shell_command(
            &mut engine,
            &ShellCommand::PlanningProfileSet {
                scope: PlanningProfileScope::ProjectOverride,
                payload_json: include_str!("../../docs/examples/planning/cloning_waiting.json")
                    .into(),
            },
        )
        .unwrap();
    }
    engine
}

fn planning_ranking_test_consult(
    engine: &mut GentleEngine,
    objective: Option<&PlanningObjective>,
) -> PlanningCloningConsultation {
    let before = serde_json::to_value(engine.state()).unwrap();
    let run = execute_shell_command(
        engine,
        &ShellCommand::PlanningConsultCloning {
            seq_id: None,
            objective_json: objective.map(|value| serde_json::to_string(value).unwrap()),
            profile_scope: PlanningProfileScope::Effective,
            output_format: "json".into(),
        },
    )
    .unwrap();
    assert!(!run.state_changed);
    assert_eq!(serde_json::to_value(engine.state()).unwrap(), before);
    serde_json::from_value(run.output).unwrap()
}

fn planning_ranking_snapshot_check(
    report: &PlanningCloningConsultation,
    expected: &[PlanningRankingSnapshotRow],
) -> Result<(), String> {
    let actual = &report.strategy_candidates;
    if actual.len() == expected.len()
        && actual
            .iter()
            .zip(expected)
            .all(|(row, &(rank, family, routine_id, score))| {
                row.rank == rank
                    && row.family == family
                    && row.routine_id == routine_id
                    && (row.estimate.composite_meta_score - score).abs() <= 1e-9
            })
    {
        return Ok(());
    }
    let mut lines = vec![format!(
        "Planning policy snapshot drift: {} expected rows, {} actual. Owner approval and a changelog entry are required to update literals.\nexpected (rank family routine score) | actual",
        expected.len(),
        actual.len()
    )];
    for index in 0..expected.len().max(actual.len()) {
        let left = expected
            .get(index)
            .map(|&(rank, family, id, score)| format!("{rank} {family} {id} {score:.12}"))
            .unwrap_or_else(|| "<none>".into());
        let right = actual
            .get(index)
            .map(|row| {
                format!(
                    "{} {} {} {:.12}",
                    row.rank, row.family, row.routine_id, row.estimate.composite_meta_score
                )
            })
            .unwrap_or_else(|| "<none>".into());
        lines.push(format!("{left} | {right}"));
    }
    Err(lines.join("\n"))
}

#[test]
fn planning_wetlab_ranking_default_absolute_policy_snapshot() {
    let report = planning_ranking_test_consult(&mut planning_ranking_test_engine(false), None);
    planning_ranking_snapshot_check(&report, &PLANNING_RANKING_DEFAULT_SNAPSHOT)
        .unwrap_or_else(|diff| panic!("{diff}"));
}

#[test]
fn planning_wetlab_ranking_waiting_absolute_policy_snapshot() {
    let report = planning_ranking_test_consult(&mut planning_ranking_test_engine(true), None);
    planning_ranking_snapshot_check(&report, &PLANNING_RANKING_WAITING_SNAPSHOT)
        .unwrap_or_else(|diff| panic!("{diff}"));
}

#[test]
fn planning_wetlab_ranking_ties_are_alphabetical_not_scientific_preferences() {
    let report = planning_ranking_test_consult(&mut planning_ranking_test_engine(false), None);
    planning_ranking_snapshot_check(&report, &PLANNING_RANKING_DEFAULT_SNAPSHOT)
        .unwrap_or_else(|diff| panic!("{diff}"));
    // Bit comparisons group only actual values from this one run; no committed
    // platform-specific f64 bit patterns or raw equality to decimal literals.
    let mut groups: Vec<Vec<&PlanningCloningStrategyCandidate>> = vec![];
    for row in &report.strategy_candidates {
        if let Some(group) = groups.iter_mut().find(|group| {
            group[0].estimate.composite_meta_score.to_bits()
                == row.estimate.composite_meta_score.to_bits()
        }) {
            group.push(row);
        } else {
            groups.push(vec![row]);
        }
    }
    assert_eq!(
        groups.iter().map(Vec::len).collect::<Vec<_>>(),
        [1, 5, 2, 1, 1, 1]
    );
    let tied_families = groups
        .iter()
        .filter(|group| group.len() > 1)
        .map(|group| {
            assert!(group.windows(2).all(|rows| rows[0].family < rows[1].family));
            group
                .iter()
                .map(|row| row.family.as_str())
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    assert_eq!(
        tied_families,
        vec![
            vec!["crispr", "golden_gate", "pcr", "ta_gc", "topo"],
            vec!["infusion", "nebuilder_hifi"]
        ]
    );
    let singletons = groups
        .iter()
        .filter(|group| group.len() == 1)
        .map(|group| group[0].estimate.composite_meta_score)
        .collect::<Vec<_>>();
    assert_eq!(singletons.len(), 4);
    for (index, left) in singletons.iter().enumerate() {
        for right in &singletons[index + 1..] {
            assert!((left - right).abs() > 1e-9);
        }
    }
}

#[test]
fn planning_wetlab_ranking_weight_sweep_is_bit_equal_across_bench_models() {
    let catalog: CloningRoutineCatalog =
        serde_json::from_str(include_str!("../../assets/cloning_routines.json")).unwrap();
    let routine = catalog
        .routines
        .iter()
        .find(|row| row.routine_id == "restriction.digest_ligate_extract_sticky")
        .unwrap();
    let normal = planning_wetlab_test_parse(&planning_wetlab_test_catalog()).unwrap();
    let mut empty = planning_wetlab_test_catalog();
    empty["routines"] = json!([]);
    let absent = planning_wetlab_test_parse(&empty).unwrap();
    let mut altered = planning_wetlab_test_catalog();
    for step in altered["routines"][0]["steps"].as_array_mut().unwrap() {
        step["hands_on_hours"] = json!({"min":10.0,"typical":20.0,"max":40.0});
    }
    let altered = planning_wetlab_test_parse(&altered).unwrap();
    for waiting in [false, true] {
        let engine = planning_ranking_test_engine(waiting);
        let before = serde_json::to_value(engine.state()).unwrap();
        for weight in [0.0, 0.03125, 0.375, 1.0, 2.0, 1e6] {
            let objective = PlanningObjective {
                weight_time: weight,
                ..Default::default()
            };
            let estimate = |bench| {
                estimate_routine_planning_with_bench_catalog(
                    &engine,
                    &objective,
                    routine,
                    &RoutinePreferenceContext::default(),
                    bench,
                )
                .unwrap()
            };
            let legacy = estimate(&absent);
            for bench in [&normal, &altered] {
                let modeled = estimate(bench);
                assert_eq!(
                    modeled.composite_meta_score.to_bits(),
                    legacy.composite_meta_score.to_bits(),
                    "weight_time={weight}, waiting={waiting}"
                );
                assert_eq!(
                    modeled.explanation["time_policy"],
                    legacy.explanation["time_policy"]
                );
                assert!((modeled.estimated_time_hours - legacy.estimated_time_hours).abs() > 1e-9);
            }
        }
        assert_eq!(serde_json::to_value(engine.state()).unwrap(), before);
    }
}

#[test]
fn planning_wetlab_ranking_time_weight_has_the_declared_policy_slope() {
    for waiting in [false, true] {
        let mut engine = planning_ranking_test_engine(waiting);
        let zero = PlanningObjective {
            weight_time: 0.0,
            ..Default::default()
        };
        let intercepts = planning_ranking_test_consult(&mut engine, Some(&zero));
        assert_eq!(intercepts.strategy_candidates.len(), 11);
        for weight in [0.375, 2.0] {
            let objective = PlanningObjective {
                weight_time: weight,
                ..Default::default()
            };
            let report = planning_ranking_test_consult(&mut engine, Some(&objective));
            assert_eq!(report.strategy_candidates.len(), 11);
            for row in &report.strategy_candidates {
                let intercept = intercepts
                    .strategy_candidates
                    .iter()
                    .find(|initial| initial.family == row.family)
                    .unwrap();
                assert_eq!(row.routine_id, intercept.routine_id);
                let policy_hours = row.estimate.explanation["time_policy"]["hours"]
                    .as_f64()
                    .unwrap();
                assert!(
                    (policy_hours
                        - intercept.estimate.explanation["time_policy"]["hours"]
                            .as_f64()
                            .unwrap())
                    .abs()
                        <= 1e-9
                );
                let expected = intercept.estimate.composite_meta_score - weight * policy_hours;
                assert!(
                    (row.estimate.composite_meta_score - expected).abs() <= 1e-9,
                    "family={}, weight_time={weight}, waiting={waiting}: expected={expected:.12}, actual={:.12}",
                    row.family,
                    row.estimate.composite_meta_score
                );
            }
        }
    }
}

#[test]
fn planning_wetlab_ranking_snapshot_rejects_a_real_weight_induced_reorder() {
    let mut engine = planning_ranking_test_engine(false);
    let baseline = planning_ranking_test_consult(&mut engine, None);
    planning_ranking_snapshot_check(&baseline, &PLANNING_RANKING_DEFAULT_SNAPSHOT).unwrap();
    let objective = PlanningObjective {
        weight_time: 0.0,
        ..Default::default()
    };
    let changed = planning_ranking_test_consult(&mut engine, Some(&objective));
    let identities = |report: &PlanningCloningConsultation| {
        report
            .strategy_candidates
            .iter()
            .map(|row| (row.rank, row.family.clone(), row.routine_id.clone()))
            .collect::<Vec<_>>()
    };
    assert_ne!(
        identities(&baseline),
        identities(&changed),
        "negative control must reorder candidates, not merely change score values"
    );
    let diff = planning_ranking_snapshot_check(&changed, &PLANNING_RANKING_DEFAULT_SNAPSHOT)
        .expect_err("guard must reject the changed consultation");
    assert!(diff.contains("expected (rank family routine score) | actual"));
}

#[test]
fn planning_wetlab_ranking_all_rows_declare_the_same_policy_basis() {
    for waiting in [false, true] {
        let report =
            planning_ranking_test_consult(&mut planning_ranking_test_engine(waiting), None);
        assert_eq!(report.strategy_candidates.len(), 11);
        for row in &report.strategy_candidates {
            let policy = &row.estimate.explanation["time_policy"];
            assert_eq!(
                policy["basis"], "legacy_family_scalar",
                "family={}",
                row.family
            );
            let consumers = policy["used_by"].as_array().unwrap();
            for consumer in ["composite_meta_score", "max_time_hours"] {
                assert!(
                    consumers.iter().any(|value| value == consumer),
                    "family={} lacks {consumer}",
                    row.family
                );
            }
            assert_eq!(
                row.estimate.estimated_time_hours_basis,
                if row.family == "restriction" {
                    crate::engine::PlanningTimeBasis::StepModelCriticalPath
                } else {
                    crate::engine::PlanningTimeBasis::LegacyFamilyScalar
                }
            );
        }
    }
}

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
fn planning_wetlab_model_presence_and_duration_do_not_change_ranking_or_guardrails() {
    let normal = planning_wetlab_test_parse(&planning_wetlab_test_catalog()).unwrap();
    let mut empty = planning_wetlab_test_catalog();
    empty["routines"] = json!([]);
    let absent = planning_wetlab_test_parse(&empty).unwrap();
    let mut slow = planning_wetlab_test_catalog();
    for step in slow["routines"][0]["steps"].as_array_mut().unwrap() {
        step["hands_on_hours"] = json!({"min":10.0,"typical":20.0,"max":40.0});
    }
    let slow = planning_wetlab_test_parse(&slow).unwrap();

    // Both sides of the new critical-path value, with and without waiting.
    for base_hours in [2.0, 6.0] {
        for waiting in [false, true] {
            let mut engine = GentleEngine::from_state(ProjectState::default());
            execute_shell_command(&mut engine, &ShellCommand::PlanningProfileSet {
                scope: PlanningProfileScope::ProjectOverride,
                payload_json: json!({"schema":"gentle.planning_profile.v1",
                    "inventory":{"mix":{"available":!waiting,"unit_cost":7.0,
                        "procurement_business_days":2.0}},
                    "machine_availability":{"machine":{"available":true,
                        "queue_business_days":if waiting {0.5} else {0.0},"run_cost_per_hour":5.0}}
                }).to_string(),
            }).unwrap();
            let before = serde_json::to_value(engine.state()).unwrap();
            let mut routine =
                planning_wetlab_test_routine("restriction.digest_ligate_extract_sticky");
            routine.base_time_hours = Some(base_hours);
            routine.base_cost = Some(3.0);
            routine.required_material_classes = vec!["mix".into()];
            routine.required_machine_classes = vec!["machine".into()];
            let delay = if waiting { 84.0 } else { 0.0 };
            let objective = PlanningObjective {
                weight_time: 1.25,
                weight_cost: 0.125,
                weight_local_fit: 2.0,
                max_time_hours: Some(5.0 + delay),
                ..Default::default()
            };
            let estimate = |catalog| {
                estimate_routine_planning_with_bench_catalog(
                    &engine,
                    &objective,
                    &routine,
                    &RoutinePreferenceContext::default(),
                    catalog,
                )
                .unwrap()
            };
            let legacy = estimate(&absent);
            assert!((legacy.estimated_time_hours - (base_hours + delay)).abs() < 1e-9);
            for catalog in [&normal, &slow] {
                let modeled = estimate(catalog);
                assert_ne!(modeled.estimated_time_hours, legacy.estimated_time_hours);
                assert_eq!(
                    modeled.composite_meta_score.to_bits(),
                    legacy.composite_meta_score.to_bits()
                );
                assert_eq!(modeled.passes_guardrails, base_hours < 5.0);
                assert_eq!(modeled.guardrail_failures, legacy.guardrail_failures);
                assert_eq!(
                    modeled.estimated_cost.to_bits(),
                    legacy.estimated_cost.to_bits()
                );
                assert_eq!(
                    modeled.explanation["time_policy"],
                    legacy.explanation["time_policy"]
                );
                assert_eq!(
                    modeled.explanation["time_policy"]["hours"]
                        .as_f64()
                        .unwrap()
                        .to_bits(),
                    legacy.estimated_time_hours.to_bits()
                );
                assert_eq!(
                    modeled.explanation["time_policy"]["basis"],
                    "legacy_family_scalar"
                );
                assert!(
                    planning_estimate_display_lines(&modeled)
                        .iter()
                        .any(|line| line.starts_with("Ranking/time-guardrail proxy:"))
                );
            }
            assert_eq!(serde_json::to_value(engine.state()).unwrap(), before);
        }
    }
}

#[test]
fn planning_wetlab_consultation_order_matches_the_no_model_counterfactual() {
    let catalog: CloningRoutineCatalog =
        serde_json::from_str(include_str!("../../assets/cloning_routines.json")).unwrap();
    let mut empty = planning_wetlab_test_catalog();
    empty["routines"] = json!([]);
    let absent = planning_wetlab_test_parse(&empty).unwrap();
    let mut engine = GentleEngine::from_state(ProjectState::default());
    for waiting in [false, true] {
        if waiting {
            execute_shell_command(
                &mut engine,
                &ShellCommand::PlanningProfileSet {
                    scope: PlanningProfileScope::ProjectOverride,
                    payload_json: include_str!("../../docs/examples/planning/cloning_waiting.json")
                        .into(),
                },
            )
            .unwrap();
        }
        for preferred in [vec![], vec!["restriction".to_string()]] {
            let objective = PlanningObjective {
                preferred_routine_families: preferred,
                ..Default::default()
            };
            let context = GentleEngine::planning_objective_routine_preference_context(&objective);
            let mut legacy = catalog
                .routines
                .iter()
                .filter(|routine| {
                    matches!(
                        routine.family.as_str(),
                        "restriction"
                            | "gibson"
                            | "golden_gate"
                            | "gateway"
                            | "infusion"
                            | "nebuilder_hifi"
                            | "pcr"
                            | "topo"
                            | "ta_gc"
                            | "crispr"
                            | "sequence"
                    )
                })
                .map(|routine| {
                    let estimate = estimate_routine_planning_with_bench_catalog(
                        &engine, &objective, routine, &context, &absent,
                    )
                    .unwrap();
                    (
                        routine.family.clone(),
                        routine.routine_id.clone(),
                        estimate.composite_meta_score,
                    )
                })
                .collect::<Vec<_>>();
            // Independent selection: sort every routine, then keep the first
            // occurrence of each family instead of production's two-stage map.
            legacy.sort_by(|a, b| {
                b.2.total_cmp(&a.2)
                    .then_with(|| a.0.cmp(&b.0))
                    .then_with(|| a.1.cmp(&b.1))
            });
            let mut seen = BTreeSet::new();
            legacy.retain(|row| seen.insert(row.0.clone()));
            assert_eq!(legacy.len(), 11);
            let before = serde_json::to_value(engine.state()).unwrap();
            let report = execute_shell_command(
                &mut engine,
                &ShellCommand::PlanningConsultCloning {
                    seq_id: None,
                    objective_json: Some(serde_json::to_string(&objective).unwrap()),
                    profile_scope: PlanningProfileScope::Effective,
                    output_format: "json".into(),
                },
            )
            .unwrap();
            assert!(!report.state_changed);
            let rows = report.output["strategy_candidates"].as_array().unwrap();
            assert_eq!(rows.len(), legacy.len());
            for (index, (actual, expected)) in rows.iter().zip(&legacy).enumerate() {
                assert_eq!(actual["rank"], index + 1);
                assert_eq!(actual["family"], expected.0);
                assert_eq!(actual["routine_id"], expected.1);
                assert_eq!(
                    actual["estimate"]["composite_meta_score"]
                        .as_f64()
                        .unwrap()
                        .to_bits(),
                    expected.2.to_bits()
                );
            }
            assert_eq!(serde_json::to_value(engine.state()).unwrap(), before);
        }
    }
}

#[test]
fn planning_wetlab_finite_display_time_cannot_hide_policy_time_overflow() {
    let mut engine = GentleEngine::from_state(ProjectState::default());
    execute_shell_command(
        &mut engine,
        &ShellCommand::PlanningProfileSet {
            scope: PlanningProfileScope::ProjectOverride,
            payload_json: json!({"schema":"gentle.planning_profile.v1",
                "machine_availability":{"machine":{"available":true,"queue_business_days":1e306}}
            })
            .to_string(),
        },
    )
    .unwrap();
    let mut routine = planning_wetlab_test_routine("restriction.digest_ligate_extract_sticky");
    routine.base_time_hours = Some(f64::MAX);
    routine.required_machine_classes = vec!["machine".into()];
    let err = estimate_routine_planning(
        &engine,
        &PlanningObjective::default(),
        &routine,
        &RoutinePreferenceContext::default(),
    )
    .unwrap_err();
    assert_eq!(
        serde_json::from_str::<Value>(&err).unwrap()["code"],
        "planning_estimate_overflow"
    );
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
    legacy["explanation"]
        .as_object_mut()
        .unwrap()
        .remove("time_policy");
    let restored: PlanningEstimate = serde_json::from_value(legacy).unwrap();
    let lines = planning_estimate_display_lines(&restored).join("\n");
    assert!(lines.contains("Hands-on: unknown"));
    assert!(lines.contains("Cost components: unknown"));
    assert!(lines.contains("Ranking/time-guardrail proxy: unknown"));
    let axes = planning_wetlab::comparison_axes(&restored, &estimate);
    let policy = axes
        .iter()
        .find(|row| row["axis"] == "time_policy_proxy_hours")
        .unwrap();
    assert!(policy["left_hours"].is_null());
    assert_eq!(policy["comparable"], false);
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
    let objective = PlanningObjective {
        weight_time: 1.25,
        weight_cost: 0.125,
        weight_local_fit: 2.0,
        ..Default::default()
    };
    let mut unmodeled = 0;
    for routine in &catalog.routines {
        let estimate = estimate_routine_planning(
            &engine,
            &objective,
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
        assert_eq!(
            estimate.explanation["time_policy"]["hours"]
                .as_f64()
                .unwrap()
                .to_bits(),
            legacy_elapsed.to_bits()
        );
        let legacy_score = (objective.weight_local_fit * estimate.local_fit_score) + 0.0
            - (objective.weight_time * legacy_elapsed)
            - (objective.weight_cost * legacy_cost);
        assert_eq!(
            estimate.composite_meta_score.to_bits(),
            legacy_score.to_bits()
        );
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
        "Ranking/time-guardrail proxy:",
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
    let policy = axes
        .iter()
        .find(|row| row["axis"] == "time_policy_proxy_hours")
        .unwrap();
    assert_eq!(policy["left_hours"], 6.0);
    assert_eq!(policy["right_hours"], 5.0);
    assert_eq!(policy["comparable"], true);
    let hands = axes
        .iter()
        .find(|row| row["axis"] == "hands_on_hours")
        .unwrap();
    assert_eq!(hands["right"], "unknown");
    assert_eq!(hands["comparable"], false);
    assert_eq!(serde_json::to_value(engine.state()).unwrap(), before);
}
