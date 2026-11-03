//! Inline synthetic expression-review fixtures, never a natural gene or provider run.
//!
//! Recreated by adapter tests from the literal 27-bp tutorial 09.03 starter and
//! the reviewed synthetic requirements in docs/examples/plans/README.md.

use crate::dna_sequence::DNAsequence;
use crate::engine::{GentleEngine, PlanningObjective, PlanningProfileScope, ProjectState};
use crate::engine_shell::{ShellCommand, execute_shell_command};
use serde_json::{Value, json};

pub(crate) const PRODUCT: &str = "expression_review_product";

pub(crate) fn fixture_state() -> ProjectState {
    let mut state = ProjectState::default();
    state.sequences.insert(
        PRODUCT.into(),
        DNAsequence::from_sequence("ATGCACCACCACCACCACCACGCTTAA").unwrap(),
    );
    let mut engine = GentleEngine::from_state(state);
    engine
        .set_planning_objective(Some(PlanningObjective {
            biological_intent: Some("saved_product_review".into()),
            ..PlanningObjective::default()
        }))
        .unwrap();
    engine.state().clone()
}

pub(crate) fn cases() -> Vec<Value> {
    let in_house: Value = serde_json::from_str(include_str!(
        "../../docs/examples/plans/protein_expression_requirements_in_house.json"
    ))
    .unwrap();
    let reviewed: Value = serde_json::from_str(include_str!(
        "../../docs/examples/plans/protein_expression_requirements_review_only.json"
    ))
    .unwrap();
    vec![
        json!({}),
        json!({"seq_id":"unresolved_synthetic_product"}),
        json!({"seq_id":PRODUCT}),
        json!({"seq_id":PRODUCT,"requirements":in_house}),
        json!({"seq_id":PRODUCT,"requirements":in_house,"provider":"geneart"}),
        json!({"seq_id":PRODUCT,"requirements":reviewed,"provider":"geneart"}),
        json!({"seq_id":PRODUCT,"provider":"geneart","requirements":{"schema":"gentle.protein_expression_requirements.v1","outsourcing":{"allowed":true,"approved_providers":["geneart"]}}}),
        json!({"seq_id":PRODUCT,"objective":{}}),
        json!({"seq_id":PRODUCT,"profile_scope":"effective","objective":{"biological_intent":"explicit_product_review","required_capabilities":[]},"requirements":{"schema":"gentle.protein_expression_requirements.v1","folding":{},"outsourcing":{"allowed":false,"approved_providers":[]},"reviewed_at_unix_ms":42}}),
    ]
}

/// Independent operand construction: do not use the new transport admission seam.
pub(crate) fn oracle(state: &ProjectState, options: &Value) -> Value {
    let string = |key| options.get(key).and_then(Value::as_str).map(str::to_string);
    let command = ShellCommand::PlanningProteinExpressionHandoff {
        seq_id: string("seq_id"),
        objective_json: options.get("objective").map(Value::to_string),
        requirements_json: options.get("requirements").map(Value::to_string),
        provider: string("provider"),
        profile_scope: PlanningProfileScope::Effective,
        output_format: string("format").unwrap_or_else(|| "json".into()),
    };
    let mut engine = GentleEngine::from_state(state.clone());
    let before = serde_json::to_value(engine.state()).unwrap();
    let run = execute_shell_command(&mut engine, &command).unwrap();
    assert!(!run.state_changed);
    assert_eq!(serde_json::to_value(engine.state()).unwrap(), before);
    run.output
}

pub(crate) fn canonical(mut report: Value) -> Value {
    assert!(
        report
            .as_object_mut()
            .unwrap()
            .remove("generated_at_unix_ms")
            .unwrap()
            .is_number()
    );
    report
}
