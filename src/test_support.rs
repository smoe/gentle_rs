//! Test-only fixtures shared across unit-test modules.
//!
//! These helpers keep adapter-parity tests anchored to one synthetic fixture so
//! future contract changes only need to be updated in one place.

use crate::dna_sequence::DNAsequence;
use crate::engine::{
    ProjectState, ROUTINE_DECISION_TRACE_SCHEMA, ROUTINE_DECISION_TRACE_STORE_SCHEMA,
    ROUTINE_DECISION_TRACES_METADATA_KEY, RoutineDecisionTrace,
    RoutineDecisionTraceDisambiguationAnswer, RoutineDecisionTraceDisambiguationQuestion,
    RoutineDecisionTracePreflightSnapshot, RoutineDecisionTraceStore,
};

/// Synthetic project state with one routine-decision trace used in parity tests.
pub(crate) fn decision_trace_fixture_state() -> ProjectState {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "s".to_string(),
        DNAsequence::from_sequence("ATGCCA").expect("sequence"),
    );
    state.metadata.insert(
        ROUTINE_DECISION_TRACES_METADATA_KEY.to_string(),
        serde_json::to_value(RoutineDecisionTraceStore {
            schema: ROUTINE_DECISION_TRACE_STORE_SCHEMA.to_string(),
            traces: vec![RoutineDecisionTrace {
                schema: ROUTINE_DECISION_TRACE_SCHEMA.to_string(),
                trace_id: "adapter_trace_1".to_string(),
                source: "gui_routine_assistant".to_string(),
                status: "preflight_failed".to_string(),
                created_at_unix_ms: 10,
                updated_at_unix_ms: 20,
                goal_text: "Assemble insert".to_string(),
                query_text: "golden gate".to_string(),
                disambiguation_questions_presented: vec![
                    RoutineDecisionTraceDisambiguationQuestion {
                        question_id: "question_a".to_string(),
                        question_text: "Question A?".to_string(),
                    },
                ],
                disambiguation_answers: vec![RoutineDecisionTraceDisambiguationAnswer {
                    question_id: "question_a".to_string(),
                    answer_text: "Answer A".to_string(),
                }],
                preflight_history: vec![RoutineDecisionTracePreflightSnapshot {
                    can_execute: false,
                    warnings: vec![],
                    errors: vec!["missing sequence".to_string()],
                    contract_source: Some("routine_catalog".to_string()),
                }],
                preflight_snapshot: None,
                ..RoutineDecisionTrace::default()
            }],
        })
        .expect("trace store"),
    );
    state
}
