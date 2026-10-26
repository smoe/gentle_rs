//! Synthetic lifecycle controls: tutorial 06.07's literal CDS and a six-base invented SNV.
//! Recreated inline/from the existing tutorial input; no biological or runtime acceptance.

use super::*;
use crate::dna_sequence::DNAsequence;
use serde_json::{Value, json};

pub(crate) fn plan(path: Option<String>, loaded: bool) -> Operation {
    let mut request: Value = serde_json::from_str(include_str!(
        "../../docs/tutorial/inputs/synthetic_sequence_design.json"
    ))
    .unwrap();
    if loaded {
        request["target"] = json!({"kind":"loaded_sequence","seq_id":"synthetic_source"});
    }
    serde_json::from_value(json!({"PlanDnaSequenceDesign":{"request":request,"path":path}}))
        .unwrap()
}

pub(crate) fn snv() -> Operation {
    serde_json::from_value(json!({"EvaluateGenomeEditIntent":{"request":{
        "schema":"gentle.genome_edit_intent_request.v1",
        "target":{"kind":"inline_sequence","sequence":"ACGTAC"},
        "reference":{"kind":"synthetic","namespace":"synthetic:slice3","reference_id":"six-base-window"},
        "anchor":{"start_1based":1001,"end_1based":1006,"strand":"plus"},
        "expected_sequence_sha256":"sha256:e76797dcfb9f3a4639fdf013e8a72d0bd943a105bfc397387335c6540c496559",
        "design_haplotype":"reference", "context":null,
        "intent":{"start_0based":1001,"end_0based_exclusive":1002,"deleted":"C","inserted":"T"}
    }}})).unwrap()
}

pub(crate) fn state() -> ProjectState {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "synthetic_source".into(),
        DNAsequence::from_sequence("ATGTTTAAATAA").unwrap(),
    );
    state
        .metadata
        .insert("synthetic_unreconciled".into(), json!({"keep":true}));
    state
}

pub(crate) fn create() -> Operation {
    Operation::CreateSequenceFromText {
        sequence_text: "ATGC".into(),
        output_id: Some("synthetic_temporary".into()),
        name: None,
        circular: false,
    }
}

pub(crate) fn live_with_redo() -> GentleEngine {
    let op = plan(None, true);
    let mut engine = GentleEngine::from_state_for(state(), LifecycleRequest::Operation(&op));
    engine.apply(create()).unwrap();
    engine.undo().unwrap();
    // Ordinary undo repairs its restored checkpoint. Deliberately inject the
    // synthetic legacy state to test inspection with an existing redo stack.
    engine.state = state();
    assert!(
        engine.state.lineage.nodes.is_empty(),
        "fixture must supply genuinely unreconciled input"
    );
    assert!(engine.redo_available() > 0);
    engine
}

fn value(engine: &GentleEngine) -> Value {
    serde_json::to_value(engine.state()).unwrap()
}

fn check_preserved(engine: &GentleEngine, baseline: &GentleEngine) {
    assert_eq!(value(engine), value(baseline));
    assert_eq!(engine.mutation_revision(), baseline.mutation_revision());
    assert_eq!(engine.structural_revision(), baseline.structural_revision());
    for (now, before) in [
        (&engine.undo_stack, &baseline.undo_stack),
        (&engine.redo_stack, &baseline.redo_stack),
    ] {
        assert_eq!(now.len(), before.len());
        assert!(
            now.iter()
                .zip(before)
                .all(|(a, b)| std::sync::Arc::ptr_eq(a, b)),
            "checkpoint contents must not be replaced"
        );
    }
}

#[test]
fn lifecycle_slice3_fact_covers_only_explicit_operations_and_nonempty_all_covered_workflows() {
    for op in [plan(None, false), plan(Some(String::new()), false), snv()] {
        assert!(LifecycleRequest::Operation(&op).preserves_project_content());
    }
    assert!(!LifecycleRequest::Operation(&create()).preserves_project_content());
    let preview_op = plan(None, true);
    let mut planner =
        GentleEngine::from_state_for(state(), LifecycleRequest::Operation(&preview_op));
    let preview = planner
        .apply(preview_op)
        .unwrap()
        .dna_sequence_design
        .unwrap();
    for legacy in [
        Operation::SaveFile {
            seq_id: "synthetic_source".into(),
            path: "unused.gb".into(),
            format: ExportFormat::GenBank,
        },
        Operation::ApplyDnaSequenceDesign {
            approval_digest: preview.approval_digest.clone().unwrap(),
            proposal: preview,
        },
        Operation::SetParameter {
            name: "max_fragments_per_container".into(),
            value: json!(123),
        },
    ] {
        assert!(!LifecycleRequest::Operation(&legacy).preserves_project_content());
    }
    for ops in [vec![], vec![create(), snv()], vec![snv(), create()]] {
        assert!(
            !LifecycleRequest::Workflow(&Workflow {
                run_id: "legacy".into(),
                ops
            })
            .preserves_project_content()
        );
    }
    for ops in [
        vec![snv(), plan(None, false)],
        vec![plan(None, false), snv()],
    ] {
        assert!(
            LifecycleRequest::Workflow(&Workflow {
                run_id: "covered".into(),
                ops
            })
            .preserves_project_content()
        );
    }
}

#[test]
fn lifecycle_slice3_factory_preserves_content_defaults_identity_and_reseeds_ids() {
    let input = state();
    let expected = serde_json::to_value(&input).unwrap();
    assert_ne!(
        serde_json::to_value(GentleEngine::from_state(input.clone()).state()).unwrap(),
        expected,
        "ordinary construction is a real repair control"
    );
    let op = plan(None, true);
    let mut a = GentleEngine::from_state_for(input.clone(), LifecycleRequest::Operation(&op));
    let b = GentleEngine::from_state_for(input, LifecycleRequest::Operation(&op));
    assert_eq!(value(&a), expected);
    assert_ne!(a.instance_id(), b.instance_id());
    assert_eq!(a.history_limit, GentleEngine::new().history_limit);
    assert_eq!(
        (a.undo_available(), a.redo_available(), a.journal_len()),
        (0, 0, 0)
    );
    a.apply(create()).unwrap();
    assert!(
        !a.state.lineage.nodes.is_empty(),
        "factory is not a sticky no-repair mode"
    );
    let mut seeded = a.state.clone();
    for node in seeded.lineage.nodes.values_mut() {
        node.created_by_op = Some("op-900".into());
    }
    let mut seeded_engine = GentleEngine::from_state_for(seeded, LifecycleRequest::Operation(&op));
    let result = seeded_engine.apply(op).unwrap();
    assert!(
        result
            .op_id
            .strip_prefix("op-")
            .unwrap()
            .parse::<u64>()
            .unwrap()
            > 900
    );
}

#[test]
fn lifecycle_slice3_apply_success_refusal_bad_path_and_cancel_preserve_state_and_redo() {
    let td = tempfile::tempdir().unwrap();
    for case in 0..5 {
        let mut engine = live_with_redo();
        let baseline = engine.clone(); // Test oracle includes inherited Arc checkpoint identities.
        let mut op = if case == 4 { snv() } else { plan(None, true) };
        if let Operation::PlanDnaSequenceDesign { request, path } = &mut op {
            if case == 1 {
                request.target =
                    gentle_protocol::sequence_design::DnaDesignTarget::LoadedSequence {
                        seq_id: "missing".into(),
                    };
            }
            if case == 2 {
                *path = Some(td.path().display().to_string());
            }
            if case == 3 {
                *path = Some(td.path().join("cancelled.json").display().to_string());
            }
        }
        let result = engine.apply_with_progress(op, |_| case != 3);
        if case == 1 || case == 2 {
            assert!(result.is_err());
        } else {
            assert!(result.is_ok());
        }
        if case == 3 {
            let report = result.unwrap().dna_sequence_design.unwrap();
            assert_eq!(
                report.status,
                gentle_protocol::sequence_design::DnaDesignStatus::Cancelled
            );
            assert!(
                td.path().join("cancelled.json").is_file(),
                "explicit cancelled report is still written"
            );
        }
        check_preserved(&engine, &baseline);
        assert_eq!(
            engine.journal_len(),
            baseline.journal_len() + usize::from(case != 1 && case != 2)
        );
        engine.redo().unwrap();
        assert!(engine.state().sequences.contains_key("synthetic_temporary"));
    }
}

#[test]
fn lifecycle_slice3_infeasible_and_exhausted_reports_are_inspections_too() {
    use gentle_protocol::sequence_design::{DesignInterval, DnaDesignStatus};
    for frozen in [false, true] {
        let mut engine = live_with_redo();
        let baseline = engine.clone();
        let mut op = plan(None, true);
        if let Operation::PlanDnaSequenceDesign { request, .. } = &mut op {
            if frozen {
                request.protected_intervals = vec![DesignInterval {
                    start_0based: 0,
                    end_0based_exclusive: 12,
                }];
            } else {
                request.max_evaluations = 1;
            }
        }
        let report = engine.apply(op).unwrap().dna_sequence_design.unwrap();
        assert_eq!(
            report.status,
            if frozen {
                DnaDesignStatus::ProvenInfeasible
            } else {
                DnaDesignStatus::SearchExhausted
            }
        );
        check_preserved(&engine, &baseline);
        assert_eq!(engine.journal_len(), baseline.journal_len() + 1);
    }
}

#[test]
fn lifecycle_slice3_workflows_keep_order_run_ids_and_history_in_both_executors() {
    for progress in [false, true] {
        for reverse in [false, true] {
            let mut ops = vec![snv(), plan(None, true)];
            if reverse {
                ops.reverse();
            }
            let wf = Workflow {
                run_id: "synthetic-lifecycle".into(),
                ops: ops.clone(),
            };
            let mut engine = live_with_redo();
            let baseline = engine.clone();
            let results = if progress {
                engine.apply_workflow_with_progress(wf, |_| true)
            } else {
                engine.apply_workflow(wf)
            }
            .unwrap();
            assert_eq!(results.len(), 2);
            check_preserved(&engine, &baseline);
            assert_eq!(engine.journal_len(), baseline.journal_len() + 2);
            let records = &engine.journal[baseline.journal_len()..];
            for (record, op) in records.iter().zip(ops) {
                assert_eq!(record.run_id, "synthetic-lifecycle");
                assert_eq!(
                    serde_json::to_value(&record.op).unwrap(),
                    serde_json::to_value(op).unwrap()
                );
            }
        }
    }
}

#[test]
fn lifecycle_slice3_workflow_failure_and_cancel_preserve_prefix_and_written_artifacts() {
    let td = tempfile::tempdir().unwrap();
    for mode in 0..4 {
        let report = td.path().join(format!("prefix-{mode}.json"));
        let mut invalid = plan(None, true);
        if let Operation::PlanDnaSequenceDesign { request, .. } = &mut invalid {
            request.target = gentle_protocol::sequence_design::DnaDesignTarget::LoadedSequence {
                seq_id: "missing".into(),
            };
        }
        let wf = Workflow {
            run_id: "prefix".into(),
            ops: vec![plan(Some(report.display().to_string()), true), invalid],
        };
        let mut engine = live_with_redo();
        let baseline = engine.clone();
        let outcome = if mode == 0 {
            engine.apply_workflow(wf)
        } else {
            engine.apply_workflow_with_progress(wf, |p| match p {
                OperationProgress::Workflow { completed, .. } => match mode {
                    1 => false,
                    2 => completed == 0,
                    _ => true,
                },
                _ => true,
            })
        };
        assert!(outcome.is_err());
        check_preserved(&engine, &baseline);
        let prefix = usize::from(mode != 1);
        assert_eq!(engine.journal_len(), baseline.journal_len() + prefix);
        assert_eq!(report.exists(), prefix == 1);
    }
    let mut engine = live_with_redo();
    let baseline = engine.clone();
    let wf = Workflow {
        run_id: "final-cancel".into(),
        ops: vec![snv()],
    };
    assert!(
        engine
            .apply_workflow_with_progress(wf, |p| !matches!(
                p,
                OperationProgress::Workflow { completed: 1, .. }
            ))
            .is_err()
    );
    check_preserved(&engine, &baseline);
    assert_eq!(engine.journal_len(), baseline.journal_len() + 1);
}

#[test]
fn lifecycle_slice3_mixed_and_empty_factory_controls_keep_legacy_repair() {
    for ops in [vec![], vec![snv(), create()], vec![create(), snv()]] {
        let input = state();
        let wf = Workflow {
            run_id: "legacy-control".into(),
            ops,
        };
        let mut engine =
            GentleEngine::from_state_for(input.clone(), LifecycleRequest::Workflow(&wf));
        assert_ne!(value(&engine), serde_json::to_value(input).unwrap());
        engine.apply_workflow(wf).unwrap();
        if engine.journal_len() > 0 {
            assert!(engine.state().sequences.contains_key("synthetic_temporary"));
        }
    }
}
