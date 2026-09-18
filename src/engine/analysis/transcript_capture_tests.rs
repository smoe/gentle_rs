//! Hand-crafted synthetic transcript loci; no clinical/private primer data.
//! Recreate by the literal sequences and mRNA/CDS qualifiers in `fixture`.
//! Exercises shared engine execution, geometry, pool scoring and persistence.

use super::*;
use crate::engine_shell::{execute_shell_command, parse_shell_line};

fn fixture(transcripts: &[(&str, &str, Option<usize>)], reverse: bool) -> GentleEngine {
    let mut sequence = String::new();
    let mut features = Vec::new();
    for (id, cdna, utr) in transcripts {
        let start = sequence.len();
        sequence.push_str(&if reverse {
            GentleEngine::reverse_complement(cdna)
        } else {
            cdna.to_string()
        });
        let end = sequence.len();
        let range = gb_io::seq::Location::simple_range(start as i64, end as i64);
        let mut qualifiers = vec![
            ("gene".into(), Some("synthetic_group".into())),
            ("transcript_id".into(), Some(id.to_string())),
        ];
        if let Some(utr) = utr {
            let (a, b) = if reverse {
                (start, end - utr)
            } else {
                (start + utr, end)
            };
            qualifiers.push(("cds_ranges_1based".into(), Some(format!("{}-{}", a + 1, b))));
        }
        features.push(gb_io::seq::Feature {
            kind: "mRNA".into(),
            location: if reverse {
                gb_io::seq::Location::Complement(Box::new(range))
            } else {
                range
            },
            qualifiers,
        });
        sequence.push_str("NNNNNNNNNNNNNNNN");
    }
    let mut dna = DNAsequence::from_sequence(&sequence).unwrap();
    *dna.features_mut() = features;
    let mut state = ProjectState::default();
    state.sequences.insert("locus".into(), dna);
    GentleEngine::from_state(state)
}

fn request() -> TranscriptCapturePoolRequest {
    TranscriptCapturePoolRequest {
        schema: TRANSCRIPT_CAPTURE_REQUEST_SCHEMA.into(),
        report_id: "capture_test".into(),
        targets: vec![TranscriptCaptureTarget {
            target_id: "arbitrary_target".into(),
            sources: vec![TranscriptCaptureSource {
                seq_id: "locus".into(),
                source_feature_id: 0,
                coverage_universe: TranscriptAssayCoverageUniverse::default(),
                annotation_release: Some("synthetic-v1".into()),
            }],
            role: TranscriptCaptureRole::SenseForward,
            window: TranscriptCaptureWindow::FivePrimeUtr,
            max_primers: 1,
            sharing_group: None,
            stage_ids: vec!["pcr".into()],
            tail_5prime: String::new(),
        }],
        fixed_oligos: vec![TranscriptCaptureFixedOligo {
            oligo_id: "oligo_dt".into(),
            full_oligo_5_to_3: "TTTTTTTTTTTTTTTTTTTTVN".into(),
            stage_ids: vec!["rt".into(), "pcr".into()],
            provenance: "Synthetic protocol oligo; not an actual order".into(),
            reorder_same_sequence: true,
        }],
        coverage_policy: TranscriptAssayCoveragePolicy::RequireAll,
        cdna_synthesis: TranscriptAssayCdnaSynthesis::OligoDt,
        search: TranscriptCaptureSearchPolicy {
            min_length_bp: 8,
            max_length_bp: 8,
            max_candidates_per_target: 32,
            ..Default::default()
        },
    }
}

fn run(
    engine: &mut GentleEngine,
    request: TranscriptCapturePoolRequest,
) -> TranscriptCapturePoolReport {
    *engine
        .apply(Operation::DesignTranscriptCapturePool { request })
        .unwrap()
        .transcript_capture_pool
        .unwrap()
}

#[test]
fn transcript_capture_shared_utr_keeps_downstream_differences_and_both_strands() {
    for reverse in [false, true] {
        let mut engine = fixture(
            &[
                ("T1", "AAACGTGCTATGCCC", Some(10)),
                ("T2", "GGACGTGCTATGCCC", Some(10)),
                ("T3", "TTACGTGCTATGCGG", Some(10)),
            ],
            reverse,
        );
        let report = run(&mut engine, request());
        assert!(report.coverage_satisfied);
        assert_eq!(report.proposed_candidate_ids.len(), 1);
        let chosen = report
            .candidates
            .iter()
            .find(|c| report.proposed_candidate_ids.contains(&c.candidate_id))
            .unwrap();
        assert_eq!(chosen.annealing_5_to_3, "ACGTGCTA");
        assert_eq!(report.captured_equivalence_groups.len(), 2);
        assert!(
            report
                .captured_equivalence_groups
                .iter()
                .any(|g| g.binding_instances.len() == 2)
        );
        for binding in &chosen.bindings {
            assert_eq!(binding.upstream_bases_omitted, 2);
            assert_eq!(binding.retained_length_bp, 13);
            assert_eq!(binding.transcript_start_0based, 2);
            let member = report
                .members
                .iter()
                .find(|m| m.member_id == binding.member_id)
                .unwrap();
            assert_eq!(member.five_prime_utr_length_bp, Some(10));
            let feature_start = member.transcript_feature_id * 31;
            assert_eq!(
                binding.source_ranges_0based,
                vec![if reverse {
                    (feature_start + 5, feature_start + 13)
                } else {
                    (feature_start + 2, feature_start + 10)
                }]
            );
        }
        assert_eq!(report.specificity_status, "unassessed");
    }
}

#[test]
fn transcript_capture_unknown_utr_is_not_dropped_or_silently_replaced() {
    let mut engine = fixture(
        &[
            ("known", "ACGTGCTATGCCC", Some(8)),
            ("unknown", "ACGTGCTATGCCC", None),
        ],
        false,
    );
    let report = run(&mut engine, request());
    assert_eq!(report.members.len(), 2);
    assert!(!report.coverage_satisfied);
    assert_eq!(report.uncovered_member_ids.len(), 1);
    let unknown = report
        .members
        .iter()
        .find(|m| m.transcript_id == "unknown")
        .unwrap();
    assert!(unknown.search_range_0based.is_none());
    let mut explicit = request();
    explicit.targets[0].window = TranscriptCaptureWindow::TranscriptRange {
        start_0based: 0,
        end_0based_exclusive: 8,
    };
    assert!(run(&mut engine, explicit).coverage_satisfied);
}

#[test]
fn transcript_capture_budget_preserves_partial_denominator() {
    let mut engine = fixture(
        &[
            ("T1", "ACGTGCTATGCCC", Some(8)),
            ("T2", "GACCTAGCTGCCC", Some(8)),
        ],
        false,
    );
    let first = run(&mut engine, request());
    assert!(!first.coverage_satisfied);
    assert_eq!(first.proposed_candidate_ids.len(), 1);
    let mut two = request();
    two.targets[0].max_primers = 2;
    let second = run(&mut engine, two);
    assert!(second.coverage_satisfied);
    assert_eq!(second.proposed_candidate_ids.len(), 2);
}

#[test]
fn transcript_capture_sharing_is_permission_not_target_count() {
    let mut engine = fixture(&[("T1", "ACGTGCTATGCCC", Some(8))], false);
    let mut req = request();
    let mut second = req.targets[0].clone();
    second.target_id = "other_target".into();
    req.targets.push(second);
    let separate = run(&mut engine, req.clone());
    assert_eq!(separate.proposed_candidate_ids.len(), 2);
    for target in &mut req.targets {
        target.sharing_group = Some("permitted_pair".into());
    }
    let shared = run(&mut engine, req);
    assert_eq!(shared.proposed_candidate_ids.len(), 1);
    assert!(shared.coverage_satisfied);
}

#[test]
fn transcript_capture_shared_retention_keeps_candidates_for_each_actual_target() {
    let mut engine = fixture(
        &[
            ("large_1", "ACGTGCTATGCCC", Some(8)),
            ("large_2", "ACGTGCTATGAAA", Some(8)),
        ],
        false,
    );
    let small = fixture(&[("small_1", "GACCTAGCTGCCC", Some(8))], false);
    engine
        .state
        .sequences
        .insert("small".into(), small.state.sequences["locus"].clone());
    let mut req = request();
    req.targets[0].sharing_group = Some("sharing_permitted".into());
    let mut second = req.targets[0].clone();
    second.target_id = "small_target".into();
    second.sources[0].seq_id = "small".into();
    req.targets.push(second);
    req.search.max_candidates_per_target = 1;

    let report = run(&mut engine, req);
    assert_eq!(report.members.len(), 3);
    assert!(report.coverage_satisfied);
    assert_eq!(report.proposed_candidate_ids.len(), 2);
    assert_eq!(report.candidates.len(), 2);
    assert!(report.candidates.iter().any(|candidate| {
        candidate.annealing_5_to_3 == "GACCTAGC" && candidate.covered_member_ids.len() == 1
    }));
}

#[test]
fn transcript_capture_retains_repeats_outside_windows_and_internal_a_runs() {
    let mut engine = fixture(
        &[("T1", "ACGTGCTACCCACGTGCTAAAAAAAAAAAAAG", Some(8))],
        false,
    );
    let report = run(&mut engine, request());
    let chosen = &report.candidates[0];
    assert_eq!(chosen.bindings.len(), 2);
    assert!(!chosen.bindings[1].within_requested_window);
    assert!(!chosen.bindings[0].internal_a_runs_0based.is_empty());
    assert!(
        chosen.bindings[0]
            .internal_a_runs_0based
            .iter()
            .all(|(a, b)| b - a >= 12)
    );
}

#[test]
fn transcript_capture_reverse_role_and_fixed_iupac_stages() {
    let mut engine = fixture(&[("T1", "ACGTGCTATGCCC", Some(8))], true);
    let mut req = request();
    req.targets[0].role = TranscriptCaptureRole::AntisenseReverse;
    let report = run(&mut engine, req.clone());
    assert_eq!(report.candidates[0].annealing_5_to_3, "TAGCACGT");
    assert_eq!(report.candidates[0].bindings[0].downstream_bases_omitted, 5);
    assert_eq!(report.interactions.len(), 1);
    req.fixed_oligos[0].stage_ids = vec!["rt".into()];
    assert!(run(&mut engine, req).interactions.is_empty());
    let metrics = dimer(b"CCCCAAAAAAAA", b"TTTTTTTTVN");
    assert!(metrics.max_complementary_run_bp >= 8);
    assert!(metrics.max_3prime_complementary_run_bp >= 8);
    assert_eq!(
        metrics.max_3prime_complementary_run_bp,
        dimer(b"TTTTTTTTVN", b"CCCCAAAAAAAA").max_3prime_complementary_run_bp
    );
}

#[test]
fn transcript_capture_strict_schema_limits_tm_and_ambiguity() {
    let mut engine = fixture(&[("T1", "ACGTNCTATGCCC", Some(8))], false);
    let report = run(&mut engine, request());
    assert!(!report.coverage_satisfied);
    assert_eq!(report.ambiguous_windows_skipped, 1);
    let mut req = request();
    req.search.tm_range = Some(TranscriptCaptureTmRange {
        min_c: 110.0,
        max_c: 120.0,
    });
    let mut canonical_engine = fixture(&[("T1", "ACGTGCTATGCCC", Some(8))], false);
    assert!(
        run(&mut canonical_engine, req.clone())
            .candidates
            .is_empty()
    );
    req.search.beam_width = 0;
    assert!(
        canonical_engine
            .apply(Operation::DesignTranscriptCapturePool { request: req })
            .is_err()
    );
    let mut value = serde_json::to_value(request()).unwrap();
    value["invented_fallback"] = json!(true);
    assert!(serde_json::from_value::<TranscriptCapturePoolRequest>(value).is_err());
}

#[test]
fn transcript_capture_shell_persistence_export_and_undo() {
    let mut engine = fixture(&[("T1", "ACGTGCTATGCCC", Some(8))], false);
    let line = format!(
        "primers design-transcript-capture-pool '{}'",
        serde_json::to_string(&request()).unwrap()
    );
    let command = parse_shell_line(&line).unwrap();
    assert!(matches!(
        &command,
        crate::engine_shell::ShellCommand::Op { .. }
    ));
    let result = execute_shell_command(&mut engine, &command).unwrap();
    assert!(result.state_changed);
    let saved = engine
        .get_transcript_capture_pool_report("capture_test")
        .unwrap();
    assert!(saved.coverage_satisfied);
    assert!(engine.project_fact_graph().facts.iter().any(|fact| {
        fact.fact == "report.exists"
            && fact.subject.id == "capture_test"
            && fact.value == Some(json!("transcript_capture_pool"))
            && fact.basis.as_ref().is_some_and(|basis| {
                basis.report_id == saved.report_id && basis.op_id.as_ref() == Some(&saved.op_id)
            })
    }));
    let state_json = serde_json::to_string(engine.state()).unwrap();
    assert!(state_json.contains(PRIMER_DESIGN_REPORTS_SCHEMA));
    let loaded = GentleEngine::from_state(serde_json::from_str(&state_json).unwrap());
    assert_eq!(
        serde_json::to_value(&saved).unwrap(),
        serde_json::to_value(
            loaded
                .get_transcript_capture_pool_report("capture_test")
                .unwrap()
        )
        .unwrap()
    );
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("capture.json");
    engine
        .export_transcript_capture_pool_report("capture_test", path.to_str().unwrap())
        .unwrap();
    assert_eq!(
        serde_json::from_slice::<serde_json::Value>(&std::fs::read(path).unwrap()).unwrap(),
        serde_json::to_value(&saved).unwrap()
    );
    let repeat = run(&mut engine, request());
    assert_eq!(repeat.request_sha256, saved.request_sha256);
    assert_eq!(
        serde_json::to_value(repeat.candidates).unwrap(),
        serde_json::to_value(saved.candidates).unwrap()
    );
    engine.undo_last_operation().unwrap();
    engine.undo_last_operation().unwrap();
    assert!(
        engine
            .get_transcript_capture_pool_report("capture_test")
            .is_err()
    );
    assert!(parse_shell_line("primers design-transcript-capture-pool '{}' ").is_err());
}
