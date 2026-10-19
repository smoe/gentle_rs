//! Hand-crafted synthetic dependency graphs; no supplier protocol or private DNA.
//! Recreated inline and consumed by the shared projection/JSON/SVG boundary tests.

use super::*;
use crate::engine::PlanningBenchStepSource;
use serde_json::{Value, json};

fn range(hours: f64) -> TimeRange {
    TimeRange {
        min: hours,
        typical: hours,
        max: hours,
    }
}

fn synthetic() -> PlanningEstimate {
    let step =
        |id: &str, dependencies: Vec<String>, hands: f64, waiting: f64, start: f64| ScheduledStep {
            step_id: id.into(),
            label: format!("Synthetic {id}"),
            depends_on: dependencies,
            hands_on_hours: range(hands),
            unattended_hours: range(waiting),
            earliest_start_hours: range(start),
            earliest_finish_hours: range(start + hands + waiting),
            source: PlanningBenchStepSource {
                kind: "synthetic_planning_default".into(),
                note: "Hand-crafted graph test, not laboratory instructions".into(),
            },
        };
    PlanningEstimate {
        step_model: PlanningStepModel::Present,
        schedule: vec![
            step("prepare_b", vec![], 0.0, 2.0, 0.0),
            step(
                "combine",
                vec!["prepare_b".into(), "prepare_a".into()],
                0.5,
                0.0,
                2.0,
            ),
            step("prepare_a", vec![], 1.0, 0.0, 0.0),
        ],
        explanation: json!({"critical_paths":[["prepare_b","combine"],["prepare_b","combine"],["prepare_b","combine"]]}),
        ..Default::default()
    }
}

fn explanation(estimate: &PlanningEstimate) -> Value {
    json!({"routine":{"routine_id":"synthetic"},"planning":{"estimate":estimate,
        "protocol_overview":build("synthetic",estimate).unwrap()}})
}

#[test]
fn planning_protocol_deterministic_parallel_geometry_and_svg_semantics() {
    let estimate = synthetic();
    let before = serde_json::to_value(&estimate).unwrap();
    let report = build("synthetic", &estimate).unwrap();
    assert_eq!(report.status, ProtocolOverviewStatus::Available);
    assert_eq!(report.nodes.len(), 3);
    assert_eq!(report.edges.len(), 2);
    let a = report
        .nodes
        .iter()
        .find(|node| node.step.step_id == "prepare_a")
        .unwrap();
    let b = report
        .nodes
        .iter()
        .find(|node| node.step.step_id == "prepare_b")
        .unwrap();
    let c = report
        .nodes
        .iter()
        .find(|node| node.step.step_id == "combine")
        .unwrap();
    assert_eq!(a.x, b.x);
    assert_ne!(a.y, b.y);
    assert!(c.x > a.x + a.width);
    assert!(!a.typical_critical_path);
    assert!(b.typical_critical_path && c.typical_critical_path);
    assert_eq!(
        report
            .edges
            .iter()
            .filter(|edge| edge.typical_critical_path)
            .count(),
        1
    );
    for node in &report.nodes {
        assert!(node.x + node.width <= report.width);
        assert!(node.y + node.height <= report.height);
        assert!(
            report
                .svg
                .contains(&format!("data-step-id=\"{}\"", node.step.step_id))
        );
        assert!(report.svg.contains(&xml(&node_details(node))));
    }
    let mut reader = quick_xml::Reader::from_str(&report.svg);
    loop {
        match reader.read_event().unwrap() {
            quick_xml::events::Event::Eof => break,
            _ => {}
        }
    }
    assert_eq!(
        serde_json::to_string(&report).unwrap(),
        serde_json::to_string(&build("synthetic", &estimate).unwrap()).unwrap()
    );
    assert_eq!(serde_json::to_value(&estimate).unwrap(), before);
    let mut shuffled = estimate.clone();
    shuffled.schedule.reverse();
    shuffled
        .schedule
        .iter_mut()
        .for_each(|step| step.depends_on.reverse());
    let reordered = build("synthetic", &shuffled).unwrap();
    let positions = |report: &ProtocolOverview| {
        report
            .nodes
            .iter()
            .map(|node| {
                (
                    node.step.step_id.clone(),
                    node.x,
                    node.y,
                    node.width,
                    node.height,
                    node.typical_critical_path,
                )
            })
            .collect::<Vec<_>>()
    };
    assert_eq!(positions(&report), positions(&reordered));
    assert_eq!(
        serde_json::to_value(&report.edges).unwrap(),
        serde_json::to_value(&reordered.edges).unwrap()
    );
    assert_ne!(
        report.source_estimate_sha256,
        reordered.source_estimate_sha256
    );
}

#[test]
fn planning_protocol_skip_column_edges_do_not_cross_step_interiors() {
    let mut estimate = synthetic();
    let mut middle = estimate.schedule[2].clone();
    middle.step_id = "middle".into();
    middle.depends_on = vec!["prepare_a".into()];
    middle.earliest_start_hours = range(1.0);
    middle.earliest_finish_hours = range(2.0);
    estimate.schedule.push(middle);
    estimate.schedule[1].depends_on = vec!["middle".into(), "prepare_b".into()];
    let report = build("synthetic", &estimate).unwrap();
    assert_eq!(
        report
            .edges
            .iter()
            .find(|edge| edge.from == "prepare_b")
            .unwrap()
            .points
            .len(),
        6
    );
    for edge in &report.edges {
        for pair in edge.points.windows(2) {
            let [a, b] = [pair[0], pair[1]];
            assert!(a[0] == b[0] || a[1] == b[1]);
            for node in &report.nodes {
                let crosses = if a[0] == b[0] {
                    a[0] > node.x
                        && a[0] < node.x + node.width
                        && a[1].min(b[1]) < node.y + node.height
                        && a[1].max(b[1]) > node.y
                } else {
                    a[1] > node.y
                        && a[1] < node.y + node.height
                        && a[0].min(b[0]) < node.x + node.width
                        && a[0].max(b[0]) > node.x
                };
                assert!(
                    !crosses,
                    "{} -> {} crosses {}",
                    edge.from, edge.to, node.step.step_id
                );
            }
        }
    }
}

#[test]
fn planning_protocol_legacy_and_unmodeled_are_unavailable_not_empty_protocols() {
    let report = build("legacy", &PlanningEstimate::default()).unwrap();
    assert_eq!(report.status, ProtocolOverviewStatus::Unavailable);
    assert!(report.nodes.is_empty() && report.edges.is_empty());
    assert!(report.typical_critical_path_ids.is_none());
    assert!(report.svg.contains("Bench-step overview unavailable"));
    assert!(from_explanation(&json!({"routine":{"routine_id":"legacy"},"planning":{"estimate":PlanningEstimate::default()}})).is_err());
    let mut estimate = synthetic();
    estimate.explanation = json!({});
    let report = build("no_path", &estimate).unwrap();
    assert!(report.typical_critical_path_ids.is_none());
    assert!(report.svg.contains("Typical critical path unavailable"));
    assert!(report.nodes.iter().all(|node| !node.typical_critical_path));
}

#[test]
fn planning_protocol_rejects_malformed_graphs_and_unbounded_work() {
    for (kind, expected) in [
        ("cycle", "protocol_overview_dependency_cycle"),
        ("unknown", "invalid_protocol_overview_dependency"),
        ("duplicate", "invalid_protocol_overview_step"),
        ("negative", "invalid_protocol_overview_step"),
        ("score", "invalid_protocol_overview_estimate"),
        ("infinite", "invalid_protocol_overview_step"),
        ("empty", "protocol_overview_step_budget"),
        ("budget", "protocol_overview_step_budget"),
        ("source", "invalid_protocol_overview_step"),
        ("path", "invalid_protocol_overview_critical_path"),
        ("absent", "protocol_overview_model_mismatch"),
    ] {
        let mut estimate = synthetic();
        match kind {
            "cycle" => {
                estimate.explanation = json!({});
                estimate.schedule[0].depends_on = vec!["combine".into()];
            }
            "unknown" => estimate.schedule[0].depends_on = vec!["absent".into()],
            "duplicate" => estimate.schedule.push(estimate.schedule[0].clone()),
            "negative" => estimate.schedule[0].hands_on_hours = range(-1.0),
            "score" => estimate.composite_meta_score = f64::INFINITY,
            "infinite" => estimate.schedule[0].hands_on_hours = range(f64::INFINITY),
            "empty" => estimate.schedule.clear(),
            "budget" => estimate.schedule = vec![estimate.schedule[0].clone(); MAX_STEPS + 1],
            "source" => estimate.schedule[0].source.kind = "unverified_protocol".into(),
            "path" => estimate.explanation = json!({"critical_paths":[[],["absent"],[]]}),
            "absent" => estimate.step_model = PlanningStepModel::Absent,
            _ => unreachable!(),
        }
        let error: Value =
            serde_json::from_str(&build("synthetic", &estimate).unwrap_err()).unwrap();
        assert_eq!(error["code"], expected, "{kind}");
    }
}

#[test]
fn planning_protocol_binding_rejects_changed_source_routine_svg_and_redigested_payload() {
    let original = explanation(&synthetic());
    from_explanation(&original).unwrap();
    for target in ["routine", "source", "svg", "geometry"] {
        let mut changed = original.clone();
        match target {
            "routine" => changed["routine"]["routine_id"] = json!("other"),
            "source" => {
                changed["planning"]["estimate"]["schedule"][0]["source"]["note"] = json!("changed")
            }
            "svg" => changed["planning"]["protocol_overview"]["svg"] = json!("<svg/>"),
            "geometry" => changed["planning"]["protocol_overview"]["nodes"][0]["x"] = json!(9000),
            _ => unreachable!(),
        }
        let mut redigested: ProtocolOverview =
            serde_json::from_value(changed["planning"]["protocol_overview"].clone()).unwrap();
        redigested.content_sha256.clear();
        redigested.content_sha256 =
            sha256_prefixed_bytes(canonical_json_string(&redigested).unwrap().as_bytes());
        changed["planning"]["protocol_overview"] = serde_json::to_value(redigested).unwrap();
        let error: Value = serde_json::from_str(&from_explanation(&changed).unwrap_err()).unwrap();
        assert_eq!(
            error["code"], "protocol_overview_binding_mismatch",
            "{target}"
        );
    }
    for newline in ["\n", "\r\n"] {
        let text = serde_json::to_string_pretty(&original)
            .unwrap()
            .replace('\n', newline);
        from_explanation(&serde_json::from_str(&text).unwrap()).unwrap();
    }
}

#[test]
fn planning_protocol_svg_escapes_source_text_and_shares_hover_details() {
    let mut estimate = synthetic();
    estimate.schedule[0].label = "<script>&\"' synthetic long label with additional words".into();
    estimate.schedule[0].source.note = "<unsafe>&\"'\nHand-crafted test".into();
    let report = build("<synthetic>&\"'", &estimate).unwrap();
    assert!(!report.svg.contains("<script>"));
    assert!(!report.svg.contains("<unsafe>"));
    assert!(report.svg.contains("&lt;script&gt;&amp;&quot;&apos;"));
    assert!(report.svg.contains(&xml(legend(&report))));
}

#[test]
fn planning_protocol_shared_shell_is_read_only_and_binds_the_exact_existing_estimate() {
    let mut engine = crate::engine::GentleEngine::default();
    let before = serde_json::to_value(engine.state()).unwrap();
    for (id, count, status) in [
        (
            "restriction.digest_ligate_extract_sticky",
            7,
            ProtocolOverviewStatus::Available,
        ),
        (
            "gibson.two_fragment_overlap_preview",
            0,
            ProtocolOverviewStatus::Unavailable,
        ),
    ] {
        let command =
            crate::engine_shell::parse_shell_line(&format!("routines explain {id}")).unwrap();
        let result = crate::engine_shell::execute_shell_command(&mut engine, &command).unwrap();
        assert!(!result.state_changed);
        let report = from_explanation(&result.output).unwrap();
        assert_eq!(report.status, status);
        assert_eq!(report.nodes.len(), count);
        let estimate: PlanningEstimate =
            serde_json::from_value(result.output["planning"]["estimate"].clone()).unwrap();
        for node in &report.nodes {
            let source = estimate
                .schedule
                .iter()
                .find(|step| step.step_id == node.step.step_id)
                .unwrap();
            assert_eq!(
                serde_json::to_value(&node.step).unwrap(),
                serde_json::to_value(source).unwrap()
            );
        }
        assert_eq!(serde_json::to_value(engine.state()).unwrap(), before);
    }
}
