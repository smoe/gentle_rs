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
use serde_json::json;
use std::{
    fs,
    path::{Path, PathBuf},
};

const DEMO_REBASE_WITHREFM: &str = "<1>EcoRI\n<2>EcoRI\n<3>GAATTC (1/5)\n<7>N\n//\n";
const DEMO_JASPAR_PFM: &str =
    ">MA0001.1 TEST\nA [ 10 0 0 0 ]\nC [ 0 10 0 0 ]\nG [ 0 0 10 0 ]\nT [ 0 0 0 10 ]\n";

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

/// Writes a minimal deterministic REBASE `.withrefm` fixture and returns its path.
pub(crate) fn write_demo_rebase_withrefm(dir: &Path) -> PathBuf {
    let path = dir.join("rebase.withrefm");
    fs::write(&path, DEMO_REBASE_WITHREFM).expect("write rebase input");
    path
}

/// Writes a minimal deterministic JASPAR PFM fixture and returns its path.
pub(crate) fn write_demo_jaspar_pfm(dir: &Path) -> PathBuf {
    let path = dir.join("motifs.pfm");
    fs::write(&path, DEMO_JASPAR_PFM).expect("write jaspar input");
    path
}

/// Writes a minimal deterministic pool-export fixture and returns its path.
pub(crate) fn write_demo_pool_json(dir: &Path) -> PathBuf {
    let path = dir.join("demo.pool.gentle.json");
    let pool_json = json!({
        "schema": "gentle.pool.v1",
        "pool_id": "demo_pool",
        "human_id": "demo",
        "member_count": 1,
        "members": [
            {
                "seq_id": "member_1",
                "human_id": "member_1",
                "name": "Member One",
                "sequence": "ATGCATGC",
                "length_bp": 8,
                "topology": "linear",
                "ends": {
                    "end_type": "blunt",
                    "forward_5": "",
                    "forward_3": "",
                    "reverse_5": "",
                    "reverse_3": ""
                }
            }
        ]
    });
    fs::write(
        &path,
        serde_json::to_string_pretty(&pool_json).expect("serialize pool json"),
    )
    .expect("write pool json");
    path
}
