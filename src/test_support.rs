//! Test-only fixtures shared across unit-test modules.
//!
//! These helpers keep adapter-parity tests anchored to one synthetic fixture so
//! future contract changes only need to be updated in one place.

use crate::dna_sequence::DNAsequence;
use crate::engine::{
    ConstructObjective, GentleEngine, HostLifecycleRole, HostRouteStep, ProjectState,
    ROUTINE_DECISION_TRACE_SCHEMA, ROUTINE_DECISION_TRACE_STORE_SCHEMA,
    ROUTINE_DECISION_TRACES_METADATA_KEY, RoutineDecisionTrace,
    RoutineDecisionTraceDisambiguationAnswer, RoutineDecisionTraceDisambiguationQuestion,
    RoutineDecisionTracePreflightSnapshot, RoutineDecisionTraceStore,
};
#[cfg(test)]
use crate::enzymes::Enzymes;
#[cfg(test)]
use gb_io::seq::Location;
#[cfg(test)]
use serde::Deserialize;
use serde_json::json;
#[cfg(test)]
use std::collections::BTreeMap;
use std::{
    fs,
    path::{Path, PathBuf},
};

const DEMO_REBASE_WITHREFM: &str = "<1>EcoRI\n<2>EcoRI\n<3>GAATTC (1/5)\n<7>N\n//\n";
const DEMO_JASPAR_PFM: &str =
    ">MA0001.1 TEST\nA [ 10 0 0 0 ]\nC [ 0 10 0 0 ]\nG [ 0 0 10 0 ]\nT [ 0 0 0 10 ]\n";
#[cfg(test)]
const DENSE_PLASMID_VISUAL_BENCHMARK: &str =
    include_str!("../test_files/fixtures/visual_benchmarks/dense_plasmid_map.json");
#[cfg(test)]
const ANTISENSE_REPEAT_DOTPLOT_CONTEXT_VISUAL_BENCHMARK: &str =
    include_str!("../test_files/fixtures/visual_benchmarks/antisense_repeat_dotplot_context.json");

#[cfg(test)]
#[derive(Debug, Deserialize)]
struct VisualBenchmarkFile {
    schema: String,
    fixture_id: String,
    sequences: Vec<VisualBenchmarkSequence>,
}

#[cfg(test)]
#[derive(Debug, Deserialize)]
struct VisualBenchmarkSequence {
    id: String,
    name: Option<String>,
    topology: Option<String>,
    #[serde(default)]
    restriction_enzymes: Vec<String>,
    sequence_segments: Vec<VisualBenchmarkSequenceSegment>,
    #[serde(default)]
    features: Vec<VisualBenchmarkFeature>,
}

#[cfg(test)]
#[derive(Debug, Deserialize)]
struct VisualBenchmarkSequenceSegment {
    literal: Option<String>,
    repeat: Option<String>,
    count: Option<usize>,
}

#[cfg(test)]
#[derive(Debug, Deserialize)]
struct VisualBenchmarkFeature {
    kind: String,
    label: Option<String>,
    strand: Option<String>,
    ranges: Vec<[usize; 2]>,
    #[serde(default)]
    qualifiers: BTreeMap<String, String>,
}

#[cfg(test)]
fn visual_benchmark_sequence_text(segments: &[VisualBenchmarkSequenceSegment]) -> String {
    let mut sequence = String::new();
    for segment in segments {
        match (&segment.literal, &segment.repeat, segment.count) {
            (Some(literal), None, None) => sequence.push_str(literal),
            (None, Some(unit), Some(count)) => {
                for _ in 0..count {
                    sequence.push_str(unit);
                }
            }
            _ => {
                panic!("visual benchmark sequence segment must have either literal or repeat+count")
            }
        }
    }
    sequence
}

#[cfg(test)]
fn visual_benchmark_location(feature: &VisualBenchmarkFeature) -> Location {
    assert!(
        !feature.ranges.is_empty(),
        "visual benchmark feature requires at least one range"
    );
    let mut parts = feature
        .ranges
        .iter()
        .map(|[start, end]| {
            assert!(
                end > start,
                "visual benchmark feature range must be non-empty"
            );
            Location::simple_range(*start as i64, *end as i64)
        })
        .collect::<Vec<_>>();
    let base = if parts.len() == 1 {
        parts.remove(0)
    } else {
        Location::Join(parts)
    };
    if feature.strand.as_deref() == Some("-") {
        Location::Complement(Box::new(base))
    } else {
        base
    }
}

#[cfg(test)]
fn visual_benchmark_dna(sequence_fixture: &VisualBenchmarkSequence) -> DNAsequence {
    let sequence_text = visual_benchmark_sequence_text(&sequence_fixture.sequence_segments);
    let mut dna = DNAsequence::from_sequence(&sequence_text).expect("visual benchmark sequence");
    if let Some(name) = sequence_fixture.name.as_ref() {
        dna.set_name(name);
    }
    if sequence_fixture
        .topology
        .as_deref()
        .is_some_and(|topology| topology.eq_ignore_ascii_case("circular"))
    {
        dna.set_circular(true);
    }
    for feature in &sequence_fixture.features {
        let mut qualifiers = Vec::new();
        if let Some(label) = feature.label.as_ref() {
            qualifiers.push(("label".into(), Some(label.clone())));
        }
        qualifiers.extend(
            feature
                .qualifiers
                .iter()
                .filter(|(key, _)| key.as_str() != "label")
                .map(|(key, value)| (key.clone().into(), Some(value.clone()))),
        );
        dna.features_mut().push(gb_io::seq::Feature {
            kind: feature.kind.clone().into(),
            location: visual_benchmark_location(feature),
            qualifiers,
        });
    }
    if !sequence_fixture.restriction_enzymes.is_empty() {
        let enzyme_names = sequence_fixture
            .restriction_enzymes
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>();
        let enzymes = Enzymes::default().restriction_enzymes_by_name(&enzyme_names);
        assert_eq!(
            enzymes.len(),
            enzyme_names.len(),
            "visual benchmark restriction-enzyme fixture references missing enzymes"
        );
        *dna.restriction_enzymes_mut() = enzymes;
        dna.set_max_restriction_enzyme_sites(None);
    }
    dna.update_computed_features();
    dna
}

#[cfg(test)]
fn visual_benchmark_state_from_json(text: &str, expected_fixture_id: &str) -> ProjectState {
    let fixture: VisualBenchmarkFile =
        serde_json::from_str(text).expect("parse visual benchmark fixture");
    assert_eq!(fixture.schema, "gentle.visual_benchmark_fixture.v1");
    assert_eq!(fixture.fixture_id, expected_fixture_id);
    let mut state = ProjectState::default();
    for sequence in &fixture.sequences {
        state
            .sequences
            .insert(sequence.id.clone(), visual_benchmark_dna(sequence));
    }
    state
}

/// Synthetic dense plasmid-map benchmark state.
#[cfg(test)]
pub fn dense_plasmid_visual_benchmark_state() -> ProjectState {
    visual_benchmark_state_from_json(DENSE_PLASMID_VISUAL_BENCHMARK, "dense_plasmid_map")
}

/// Synthetic genomic dotplot-context benchmark state.
#[cfg(test)]
pub fn antisense_repeat_dotplot_context_visual_benchmark_state() -> ProjectState {
    visual_benchmark_state_from_json(
        ANTISENSE_REPEAT_DOTPLOT_CONTEXT_VISUAL_BENCHMARK,
        "antisense_repeat_dotplot_context",
    )
}

fn push_zip_u16(out: &mut Vec<u8>, value: u16) {
    out.extend_from_slice(&value.to_le_bytes());
}

fn push_zip_u32(out: &mut Vec<u8>, value: u32) {
    out.extend_from_slice(&value.to_le_bytes());
}

fn zip_crc32(bytes: &[u8]) -> u32 {
    let mut crc = 0xffff_ffffu32;
    for &byte in bytes {
        crc ^= byte as u32;
        for _ in 0..8 {
            let mask = (crc & 1).wrapping_neg();
            crc = (crc >> 1) ^ (0xedb8_8320 & mask);
        }
    }
    !crc
}

/// Write a minimal ZIP archive with stored members for tests.
///
/// This avoids depending on a host `zip` binary while still exercising the
/// production unzip-backed resource import paths.
pub fn write_stored_zip_archive(path: &Path, members: &[(&str, &[u8])]) {
    let mut out = Vec::new();
    let mut central = Vec::new();

    for (name, data) in members {
        let name_bytes = name.as_bytes();
        assert!(
            name_bytes.len() <= u16::MAX as usize,
            "synthetic ZIP member name too long"
        );
        assert!(
            data.len() <= u32::MAX as usize,
            "synthetic ZIP member too large"
        );
        assert!(
            out.len() <= u32::MAX as usize,
            "synthetic ZIP local offset too large"
        );

        let local_offset = out.len() as u32;
        let crc = zip_crc32(data);
        let size = data.len() as u32;
        let name_len = name_bytes.len() as u16;

        push_zip_u32(&mut out, 0x0403_4b50);
        push_zip_u16(&mut out, 20);
        push_zip_u16(&mut out, 0);
        push_zip_u16(&mut out, 0);
        push_zip_u16(&mut out, 0);
        push_zip_u16(&mut out, 0);
        push_zip_u32(&mut out, crc);
        push_zip_u32(&mut out, size);
        push_zip_u32(&mut out, size);
        push_zip_u16(&mut out, name_len);
        push_zip_u16(&mut out, 0);
        out.extend_from_slice(name_bytes);
        out.extend_from_slice(data);

        push_zip_u32(&mut central, 0x0201_4b50);
        push_zip_u16(&mut central, 20);
        push_zip_u16(&mut central, 20);
        push_zip_u16(&mut central, 0);
        push_zip_u16(&mut central, 0);
        push_zip_u16(&mut central, 0);
        push_zip_u16(&mut central, 0);
        push_zip_u32(&mut central, crc);
        push_zip_u32(&mut central, size);
        push_zip_u32(&mut central, size);
        push_zip_u16(&mut central, name_len);
        push_zip_u16(&mut central, 0);
        push_zip_u16(&mut central, 0);
        push_zip_u16(&mut central, 0);
        push_zip_u16(&mut central, 0);
        push_zip_u32(&mut central, 0);
        push_zip_u32(&mut central, local_offset);
        central.extend_from_slice(name_bytes);
    }

    assert!(
        members.len() <= u16::MAX as usize,
        "synthetic ZIP has too many members"
    );
    assert!(
        out.len() <= u32::MAX as usize && central.len() <= u32::MAX as usize,
        "synthetic ZIP central directory too large"
    );

    let central_offset = out.len() as u32;
    let central_size = central.len() as u32;
    out.extend_from_slice(&central);
    push_zip_u32(&mut out, 0x0605_4b50);
    push_zip_u16(&mut out, 0);
    push_zip_u16(&mut out, 0);
    push_zip_u16(&mut out, members.len() as u16);
    push_zip_u16(&mut out, members.len() as u16);
    push_zip_u32(&mut out, central_size);
    push_zip_u32(&mut out, central_offset);
    push_zip_u16(&mut out, 0);

    fs::write(path, out).expect("write synthetic ZIP archive");
}

/// Synthetic project state with one routine-decision trace used in parity tests.
pub fn decision_trace_fixture_state() -> ProjectState {
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

/// Synthetic project state with one routine trace plus construct-reasoning metadata.
pub fn decision_trace_with_construct_reasoning_fixture_state() -> ProjectState {
    let state = decision_trace_fixture_state();
    let mut engine = GentleEngine::from_state(state);
    let objective = engine
        .upsert_construct_objective(ConstructObjective {
            title: "Adapter construct reasoning".to_string(),
            goal: "Expose host/helper selection context in run bundles".to_string(),
            propagation_host_profile_id: Some("ecoli_k12".to_string()),
            expression_host_profile_id: Some("ecoli_k12".to_string()),
            helper_profile_id: Some("pUC19".to_string()),
            medium_conditions: vec!["ampicillin".to_string()],
            host_route: vec![HostRouteStep {
                step_id: "storage".to_string(),
                host_profile_id: "ecoli_k12".to_string(),
                role: HostLifecycleRole::Storage,
                rationale: "Retain one documented propagation/storage step.".to_string(),
                notes: vec![],
            }],
            ..ConstructObjective::default()
        })
        .expect("fixture objective");
    engine
        .build_construct_reasoning_graph("s", Some(&objective.objective_id), None)
        .expect("fixture graph");
    engine.state().clone()
}

/// Returns a minimal deterministic REBASE `.withrefm` text fixture.
pub fn demo_rebase_withrefm_text() -> &'static str {
    DEMO_REBASE_WITHREFM
}

/// Returns a minimal deterministic JASPAR PFM text fixture.
pub fn demo_jaspar_pfm_text() -> &'static str {
    DEMO_JASPAR_PFM
}

/// Writes a minimal workflow JSON file with the given `run_id`.
pub fn write_demo_workflow_json(dir: &Path, filename: &str, run_id: &str) -> PathBuf {
    let path = dir.join(filename);
    let payload = format!(r#"{{"run_id":"{run_id}","ops":[]}}"#);
    fs::write(&path, payload).expect("write workflow");
    path
}

/// Writes a shebang-prefixed workflow JSON file with the given `run_id`.
pub fn write_demo_workflow_with_shebang(dir: &Path, filename: &str, run_id: &str) -> PathBuf {
    let path = dir.join(filename);
    let payload =
        format!("#!/usr/bin/env -S gentle_cli workflow\n{{\"run_id\":\"{run_id}\",\"ops\":[]}}\n");
    fs::write(&path, payload).expect("write workflow");
    path
}

/// Writes a minimal deterministic REBASE `.withrefm` fixture and returns its path.
pub fn write_demo_rebase_withrefm(dir: &Path) -> PathBuf {
    let path = dir.join("rebase.withrefm");
    fs::write(&path, demo_rebase_withrefm_text()).expect("write rebase input");
    path
}

/// Writes a minimal deterministic JASPAR PFM fixture and returns its path.
pub fn write_demo_jaspar_pfm(dir: &Path) -> PathBuf {
    let path = dir.join("motifs.pfm");
    fs::write(&path, demo_jaspar_pfm_text()).expect("write jaspar input");
    path
}

/// Writes a minimal deterministic pool-export fixture and returns its path.
pub fn write_demo_pool_json(dir: &Path) -> PathBuf {
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
