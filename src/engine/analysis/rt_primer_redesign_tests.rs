//! Synthetic cDNA/CDS/FASTA fixtures invented here and recreated by the helpers below.
//! No real gene, private read, physical oligo or reference database is used. Temporary
//! byte-bound FASTA/request files exercise the actual shared intake and shell routes.

use super::*;
use crate::engine_shell::{execute_shell_command, parse_shell_line, quote_shell_arg};

const TAIL: &str = "ACGTCAGTACGA";
const HANDLE: &str = "TTTTTT";

fn word() -> String {
    GentleEngine::reverse_complement(TAIL)
}
fn member(id: &str) -> RtPrimerRedesignMember {
    RtPrimerRedesignMember {
        seq_id: "locus".into(),
        transcript_id: id.into(),
    }
}
fn request(ids: &[&str], start: usize, end: usize) -> RedesignRtPrimersRequest {
    RedesignRtPrimersRequest {
        schema: RT_PRIMER_REDESIGN_REQUEST_SCHEMA.into(),
        group_id: "invented_group".into(),
        mode: RtPrimerRedesignMode::SharedTranscripts,
        sources: vec![RtPrimerComparisonSource {
            source: TranscriptCaptureSource {
                seq_id: "locus".into(),
                source_feature_id: 0,
                coverage_universe: TranscriptAssayCoverageUniverse {
                    kind: TranscriptAssayCoverageUniverseKind::ExplicitTranscripts,
                    required_transcript_ids: ids.iter().map(|id| (*id).into()).collect(),
                    ..Default::default()
                },
                annotation_release: Some("synthetic-v1".into()),
            },
            expected_source_record_sha256: None,
        }],
        required_members: ids.iter().map(|id| member(id)).collect(),
        anchor: RtPrimerRedesignAnchor {
            member: member(ids[0]),
            start_0based: start,
            end_0based_exclusive: end,
        },
        common_handle: HANDLE.into(),
        constraints: RtPrimerRedesignConstraints {
            min_tail_length_bp: 12,
            max_tail_length_bp: 12,
            min_tail_gc_fraction: 0.0,
            max_tail_gc_fraction: 1.0,
            min_tail_tm_c: -100.0,
            max_tail_tm_c: 200.0,
            max_homopolymer_run_bp: 64,
            max_self_complementary_run_bp: 64,
            max_self_3prime_complementary_run_bp: 64,
        },
        references: vec![RtPrimerBackground {
            reference_id: "bg".into(),
            scope: RtPrimerBackgroundScope::CustomRna,
            source: "Explicitly unavailable invented background".into(),
            fasta: None,
            unavailable_reason: Some("No biological background in this fixture".into()),
            circular_record_ids: vec![],
        }],
        seed_limits: vec![],
        budgets: RtPrimerRedesignBudgets::default(),
        assessment_budgets: RtPrimerAssessmentBudgets::default(),
        comparison_budgets: RtPrimerComparisonBudgets::default(),
    }
}

/// Each tuple supplies biological 5′–3′ cDNA and optional annotated local CDS intervals.
fn fixture(
    transcripts: &[(&str, &str, Option<Vec<(usize, usize)>>)],
    reverse: bool,
) -> GentleEngine {
    let mut sequence = String::new();
    let mut features = Vec::new();
    for (id, cdna, cds) in transcripts {
        let start = sequence.len();
        sequence.push_str(&if reverse {
            GentleEngine::reverse_complement(cdna)
        } else {
            cdna.to_string()
        });
        let end = sequence.len();
        let range = gb_io::seq::Location::simple_range(start as i64, end as i64);
        let mut qualifiers = vec![
            ("gene".into(), Some("invented_group".into())),
            ("transcript_id".into(), Some(id.to_string())),
        ];
        if let Some(cds) = cds {
            let value = cds
                .iter()
                .map(|&(a, z)| {
                    let (a, z) = if reverse {
                        (end - z, end - a)
                    } else {
                        (start + a, start + z)
                    };
                    format!("{}-{}", a + 1, z)
                })
                .collect::<Vec<_>>()
                .join(",");
            qualifiers.push(("cds_ranges_1based".into(), Some(value)));
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
fn run(engine: &GentleEngine, req: RedesignRtPrimersRequest) -> RtPrimerRedesignReport {
    engine.redesign_rt_primers(req, &mut |_| true).unwrap()
}
fn background(req: &mut RedesignRtPrimersRequest, root: &std::path::Path, fasta: &str) {
    let path = root.join("invented background's input.fa");
    std::fs::write(&path, fasta.as_bytes()).unwrap();
    req.references[0].fasta = Some(PrimerPayloadFile {
        path: path.to_str().unwrap().into(),
        sha256: sha256_prefixed_bytes(fasta.as_bytes()),
        bytes: fasta.len() as u64,
    });
    req.references[0].unavailable_reason = None;
}
fn limit(max: u64) -> RtPrimerRedesignSeedLimit {
    RtPrimerRedesignSeedLimit {
        reference_id: "bg".into(),
        seed_length: 8,
        max_mismatches: 0,
        max_compatible_windows: max,
    }
}

#[test]
fn rt_primer_redesign_shared_sites_at_different_coordinates_on_both_strands() {
    for reverse in [false, true] {
        let first = format!("AAAA{}TT", word());
        let second = format!("GG{}CCCC", word());
        let engine = fixture(
            &[("first", &first, None), ("second", &second, None)],
            reverse,
        );
        let report = run(&engine, request(&["first", "second"], 4, 16));
        assert_eq!(report.required_member_ids.len(), 2);
        assert_eq!(report.enumerated_windows, 1);
        assert_eq!(report.qualifying_distinct_tails, 1);
        assert_eq!(report.candidates[0].target_tail, TAIL);
        let comparison = report.comparison.unwrap();
        assert_eq!(comparison.candidates[0].members_with_exact_site, 2);
        let starts = comparison
            .isoforms
            .iter()
            .map(|r| r.bindings[0].transcript_start_0based)
            .collect::<std::collections::BTreeSet<_>>();
        assert_eq!(starts, [2, 4].into_iter().collect());
        assert!(
            comparison
                .members
                .iter()
                .all(|m| m.strand == if reverse { "-" } else { "+" })
        );
        let primer = &comparison.request.assessment.primers[0];
        assert_eq!(primer.sequence, format!("{HANDLE}{TAIL}"));
        assert_eq!(
            primer.core_orientation,
            PrimerCoreOrientation::Antisense3prime
        );
    }
}

#[test]
fn rt_primer_redesign_exon_junction_uses_shared_projection_on_both_strands() {
    for reverse in [false, true] {
        let cdna = format!("AAAA{}TTTT", word());
        let genomic = format!("{}NNNNNNNN{}", &cdna[..10], &cdna[10..]);
        let mut dna = DNAsequence::from_sequence(&if reverse {
            GentleEngine::reverse_complement(&genomic)
        } else {
            genomic
        })
        .unwrap();
        let join = gb_io::seq::Location::Join(vec![
            gb_io::seq::Location::simple_range(0, 10),
            gb_io::seq::Location::simple_range(18, 28),
        ]);
        dna.features_mut().push(gb_io::seq::Feature {
            kind: "mRNA".into(),
            location: if reverse {
                gb_io::seq::Location::Complement(Box::new(join))
            } else {
                join
            },
            qualifiers: vec![
                ("gene".into(), Some("invented_group".into())),
                ("transcript_id".into(), Some("split".into())),
            ],
        });
        let mut state = ProjectState::default();
        state.sequences.insert("locus".into(), dna);
        let report = run(&GentleEngine::from_state(state), request(&["split"], 4, 16));
        let comparison = report.comparison.unwrap();
        let binding = &comparison.isoforms[0].bindings[0];
        assert_eq!(
            binding.source_ranges_0based,
            if reverse {
                vec![(18, 24), (4, 10)]
            } else {
                vec![(4, 10), (18, 24)]
            }
        );
        assert_eq!(report.candidates[0].target_tail, TAIL);
    }
}

#[test]
fn rt_primer_redesign_three_prime_utr_requires_each_member_and_all_sites() {
    let a = format!("AAAA{}", word());
    let b = format!("GG{}", word());
    for reverse in [false, true] {
        let engine = fixture(
            &[
                ("first", &a, Some(vec![(0, 4)])),
                ("second", &b, Some(vec![(0, 2)])),
            ],
            reverse,
        );
        let mut req = request(&["first", "second"], 4, 16);
        req.mode = RtPrimerRedesignMode::ThreePrimeUtr;
        let report = run(&engine, req);
        assert_eq!(report.candidates.len(), 1);
        assert!(
            report
                .comparison
                .unwrap()
                .isoforms
                .iter()
                .flat_map(|r| &r.bindings)
                .all(|b| b.region == RtPrimerBindingRegion::ThreePrimeUtr)
        );
    }
    let repeated = format!("{}AAA{}", word(), word());
    let engine = fixture(
        &[
            ("first", &a, Some(vec![(0, 4)])),
            ("second", &repeated, Some(vec![(0, 12)])),
        ],
        false,
    );
    let mut req = request(&["first", "second"], 4, 16);
    assert_eq!(run(&engine, req.clone()).candidates.len(), 1);
    req.mode = RtPrimerRedesignMode::ThreePrimeUtr;
    let report = run(&engine, req);
    assert!(report.candidates.is_empty());
    assert_eq!(report.group_rejected_tails, 1);
    assert!(report.comparison.is_none());
}

#[test]
fn rt_primer_redesign_unknown_or_ambiguity_compatible_non_utr_sites_do_not_qualify() {
    let a = format!("AAAA{}", word());
    let mut possible = word();
    possible.replace_range(0..1, "N");
    let ambiguous_coding = format!("{possible}AAA{}", word());
    for (second, cds) in [(&a, None), (&ambiguous_coding, Some(vec![(0, 12)]))] {
        let engine = fixture(
            &[("first", &a, Some(vec![(0, 4)])), ("second", second, cds)],
            false,
        );
        let mut req = request(&["first", "second"], 4, 16);
        req.mode = RtPrimerRedesignMode::ThreePrimeUtr;
        let report = run(&engine, req);
        assert!(report.candidates.is_empty());
        assert_eq!(report.group_rejected_tails, 1);
    }
}

#[test]
fn rt_primer_redesign_ambiguous_only_group_member_does_not_replace_exact_coverage() {
    let exact = word();
    let mut possible = exact.clone();
    possible.replace_range(3..4, "N");
    let engine = fixture(
        &[("first", &exact, None), ("second", &possible, None)],
        false,
    );
    let report = run(&engine, request(&["first", "second"], 0, 12));
    assert!(report.candidates.is_empty());
    assert_eq!(report.group_rejected_tails, 1);
    let report = run(&engine, request(&["second"], 0, 12));
    assert_eq!(report.skipped_ambiguous_windows, 1);
    assert!(report.candidates.is_empty());
}

#[test]
fn rt_primer_redesign_comparison_only_members_preserve_larger_denominator() {
    let exact = word();
    let engine = fixture(
        &[
            ("first", &exact, None),
            ("second", &exact, None),
            ("outside", "AAAAAAAAAAAA", None),
        ],
        false,
    );
    let mut req = request(&["first", "second", "outside"], 0, 12);
    req.required_members = vec![member("first"), member("second")];
    let report = run(&engine, req);
    assert_eq!(report.required_member_ids.len(), 2);
    let comparison = report.comparison.unwrap();
    assert_eq!(comparison.members.len(), 3);
    assert_eq!(
        (
            comparison.candidates[0].total_members,
            comparison.candidates[0].members_with_exact_site,
            comparison.candidates[0].members_without_site
        ),
        (3, 2, 1)
    );
}

#[test]
fn rt_primer_redesign_repeated_tails_are_deduplicated_and_shortlist_order_is_explicit() {
    let repeated = format!("{}AAA{}", word(), word());
    let engine = fixture(&[("first", &repeated, None)], false);
    let mut req = request(&["first"], 0, repeated.len());
    req.budgets.max_candidates = 2;
    let report = run(&engine, req);
    assert_eq!(report.duplicate_tail_windows, 1);
    assert_eq!(
        report.qualifying_distinct_tails,
        report.enumerated_windows - 1
    );
    assert!(report.shortlist_truncated);
    assert_eq!(report.candidates.len(), 2);
    assert_eq!(report.candidates[0].anchor_start_0based, 0);
    assert_eq!(report.candidates[1].anchor_start_0based, 1);
    let sequence = "ACGTCAGTACGATC";
    let engine = fixture(&[("first", sequence, None)], false);
    let mut req = request(&["first"], 0, sequence.len());
    req.constraints.min_tail_length_bp = 10;
    req.constraints.max_tail_length_bp = 11;
    req.budgets.max_candidates = 2;
    let report = run(&engine, req);
    assert_eq!(
        (
            report.candidates[0].anchor_start_0based,
            report.candidates[0].anchor_end_0based_exclusive
        ),
        (0, 10)
    );
    assert_eq!(
        (
            report.candidates[1].anchor_start_0based,
            report.candidates[1].anchor_end_0based_exclusive
        ),
        (0, 11)
    );
}

#[test]
fn rt_primer_redesign_metric_constraints_include_fixed_handle() {
    let exact = word();
    let engine = fixture(&[("first", &exact, None)], false);
    let mut req = request(&["first"], 0, 12);
    req.constraints.max_homopolymer_run_bp = 5;
    let report = run(&engine, req);
    assert_eq!(report.metric_rejected_tails, 1);
    assert!(report.candidates.is_empty());
    let mut req = request(&["first"], 0, 12);
    req.constraints.min_tail_gc_fraction = 1.0;
    let report = run(&engine, req);
    assert_eq!(report.metric_rejected_tails, 1);
    assert!(report.candidates.is_empty());
}

#[test]
fn rt_primer_redesign_seed_limits_use_all_sites_and_retain_failed_candidates() {
    let td = tempfile::tempdir().unwrap();
    let exact = word();
    let engine = fixture(&[("first", &exact, None)], false);
    let mut req = request(&["first"], 0, 12);
    background(&mut req, td.path(), &format!(">invented\n{exact}\n"));
    req.seed_limits = vec![limit(0)];
    let report = run(&engine, req.clone());
    assert_eq!(report.candidates.len(), 1);
    assert_eq!(
        report.candidates[0].seed_constraint_status,
        RtPrimerRedesignConstraintStatus::Failed
    );
    assert_eq!(
        report.candidates[0].seed_checks[0].compatible_windows,
        Some(1)
    );
    assert!(
        report
            .comparison
            .as_ref()
            .unwrap()
            .assessment
            .request
            .intended_sites
            .is_empty()
    );
    req.seed_limits[0].max_compatible_windows = 1;
    assert_eq!(
        run(&engine, req.clone()).candidates[0].seed_constraint_status,
        RtPrimerRedesignConstraintStatus::Passed
    );
    // The constraint sums the marginal positional-mismatch bins, not terminal runs.
    let mismatching_seed_template = GentleEngine::reverse_complement("TTGTACGA");
    background(
        &mut req,
        td.path(),
        &format!(">two_mismatches\n{mismatching_seed_template}\n"),
    );
    req.seed_limits[0].max_compatible_windows = 0;
    req.seed_limits[0].max_mismatches = 1;
    assert_eq!(
        run(&engine, req.clone()).candidates[0].seed_constraint_status,
        RtPrimerRedesignConstraintStatus::Passed
    );
    req.seed_limits[0].max_mismatches = 2;
    let report = run(&engine, req);
    assert_eq!(
        report.candidates[0].seed_checks[0].compatible_windows,
        Some(1)
    );
    assert_eq!(
        report.candidates[0].seed_constraint_status,
        RtPrimerRedesignConstraintStatus::Failed
    );
}

#[test]
fn rt_primer_redesign_missing_and_ambiguous_backgrounds_are_not_clean_zero() {
    let td = tempfile::tempdir().unwrap();
    let exact = word();
    let engine = fixture(&[("first", &exact, None)], false);
    let mut req = request(&["first"], 0, 12);
    req.seed_limits = vec![limit(0)];
    let report = run(&engine, req.clone());
    assert_eq!(
        report.candidates[0].seed_constraint_status,
        RtPrimerRedesignConstraintStatus::Unresolved
    );
    assert_eq!(report.candidates[0].seed_checks[0].compatible_windows, None);
    background(&mut req, td.path(), ">ambiguous\nNNNNNNNN\n");
    let report = run(&engine, req.clone());
    assert_eq!(
        report.candidates[0].seed_constraint_status,
        RtPrimerRedesignConstraintStatus::Unresolved
    );
    assert_eq!(
        report.candidates[0].seed_checks[0].ambiguous_windows,
        Some(1)
    );
    background(
        &mut req,
        td.path(),
        &format!(">exact\n{exact}\n>ambiguous\nNNNNNNNN\n"),
    );
    assert_eq!(
        run(&engine, req).candidates[0].seed_constraint_status,
        RtPrimerRedesignConstraintStatus::Failed
    );
}

#[test]
fn rt_primer_redesign_failed_shortlist_is_not_refilled_and_short_seeds_are_unresolved() {
    let td = tempfile::tempdir().unwrap();
    let sequence = format!("{}G", word());
    let engine = fixture(&[("first", &sequence, None)], false);
    let mut req = request(&["first"], 0, sequence.len());
    req.budgets.max_candidates = 1;
    background(
        &mut req,
        td.path(),
        &format!(">first_seed\n{}\n", &word()[..8]),
    );
    req.seed_limits = vec![limit(0)];
    let report = run(&engine, req.clone());
    assert_eq!(report.qualifying_distinct_tails, 2);
    assert!(report.shortlist_truncated);
    assert_eq!(report.candidates.len(), 1);
    assert_eq!(
        report.candidates[0].seed_constraint_status,
        RtPrimerRedesignConstraintStatus::Failed
    );
    req.anchor.start_0based = 1;
    assert_eq!(
        run(&engine, req.clone()).candidates[0].seed_constraint_status,
        RtPrimerRedesignConstraintStatus::Passed
    );
    req.constraints.min_tail_length_bp = 10;
    req.constraints.max_tail_length_bp = 10;
    req.seed_limits[0].seed_length = 12;
    let report = run(&engine, req);
    assert_eq!(
        report.candidates[0].seed_constraint_status,
        RtPrimerRedesignConstraintStatus::Unresolved
    );
    assert_eq!(report.candidates[0].seed_checks[0].compatible_windows, None);
    let mut polls = 0;
    let error = engine
        .redesign_rt_primers(request(&["first"], 0, sequence.len()), &mut |_| {
            polls += 1;
            polls < 8
        })
        .unwrap_err();
    assert!(error.message.contains("cancelled"));
    assert_eq!(polls, 8);
}

#[test]
fn rt_primer_redesign_empty_generation_keeps_provenance_without_opening_backgrounds() {
    let exact = word();
    let engine = fixture(
        &[("first", &exact, None), ("second", "AAAAAAAAAAAA", None)],
        false,
    );
    let mut req = request(&["first", "second"], 0, 12);
    req.references[0].unavailable_reason = None;
    req.references[0].fasta = Some(PrimerPayloadFile {
        path: "/not-an-existing-redesign-fixture.fa".into(),
        sha256: format!("sha256:{}", "0".repeat(64)),
        bytes: 100,
    });
    let report = run(&engine, req);
    assert!(report.candidates.is_empty());
    assert!(report.comparison.is_none());
    assert_eq!(report.members.len(), 2);
    assert!(!report.sources[0].source_record_sha256.is_empty());
    assert_eq!(report.required_member_ids.len(), 2);
}

#[test]
fn rt_primer_redesign_explicit_group_anchor_and_source_guard_fail_closed() {
    let exact = word();
    let mut engine = fixture(&[("first", &exact, None)], false);
    let initial = run(&engine, request(&["first"], 0, 12));
    for which in 0..6 {
        let mut req = request(&["first"], 0, 12);
        match which {
            0 => {
                req.sources[0].source.coverage_universe = TranscriptAssayCoverageUniverse::default()
            }
            1 => req.required_members.push(member("missing")),
            2 => req.required_members.push(member("first")),
            3 => req.anchor.member = member("missing"),
            4 => req.anchor.end_0based_exclusive = 13,
            _ => req.required_members.clear(),
        }
        assert!(
            engine.redesign_rt_primers(req, &mut |_| true).is_err(),
            "case {which}"
        );
    }
    let mut req = request(&["first"], 0, 12);
    req.sources[0].expected_source_record_sha256 =
        Some(initial.sources[0].source_record_sha256.clone());
    engine
        .state_mut()
        .sequences
        .get_mut("locus")
        .unwrap()
        .features_mut()[0]
        .qualifiers
        .push(("note".into(), Some("annotation-only edit".into())));
    assert!(engine.redesign_rt_primers(req, &mut |_| true).is_err());
}

#[test]
fn rt_primer_redesign_budgets_schema_and_cancellation_refuse_instead_of_empty_success() {
    let sequence = format!("{}AAA{}", word(), word());
    let engine = fixture(&[("first", &sequence, None)], false);
    for which in 0..7 {
        let mut req = request(&["first"], 0, sequence.len());
        match which {
            0 => req.budgets.max_enumerated_windows = 1,
            1 => req.budgets.max_generation_work = 1,
            2 => req.budgets.max_candidates = 17,
            3 => req.constraints.min_tail_length_bp = 9,
            4 => req.common_handle = "NNNN".into(),
            5 => {
                req.seed_limits = vec![RtPrimerRedesignSeedLimit {
                    seed_length: 7,
                    ..limit(0)
                }]
            }
            _ => req.comparison_budgets.max_transcript_base_comparisons = 1,
        }
        assert!(
            engine.redesign_rt_primers(req, &mut |_| true).is_err(),
            "case {which}"
        );
    }
    assert!(
        engine
            .redesign_rt_primers(request(&["first"], 0, sequence.len()), &mut |_| false)
            .is_err()
    );
    let mut value = serde_json::to_value(request(&["first"], 0, 12)).unwrap();
    value["implicit_utr_fallback"] = serde_json::json!(true);
    assert!(serde_json::from_value::<RedesignRtPrimersRequest>(value).is_err());
}

#[test]
fn rt_primer_redesign_composes_existing_comparison_and_deterministic_hashes() {
    let exact = word();
    let engine = fixture(&[("first", &exact, None)], false);
    let req = request(&["first"], 0, 12);
    let first = run(&engine, req.clone());
    assert_eq!(first.content_sha256, run(&engine, req).content_sha256);
    let comparison = first.comparison.unwrap();
    let direct = engine
        .compare_rt_primer_candidates(comparison.request.clone(), &mut |_| true)
        .unwrap();
    assert_eq!(
        serde_json::to_value(&*comparison).unwrap(),
        serde_json::to_value(direct).unwrap()
    );
    assert_eq!(comparison.assessment.references[0].status, "not_assessed");
}

#[test]
fn rt_primer_redesign_typed_and_shared_shell_at_file_preserve_project_and_history() {
    let exact = word();
    let mut engine = fixture(&[("first", &exact, None)], true);
    let req = request(&["first"], 0, 12);
    let before = serde_json::to_value(engine.snapshot()).unwrap();
    let revision = engine.mutation_revision();
    let undo = engine.undo_available();
    let result = engine
        .apply(Operation::RedesignRtPrimers {
            request: Box::new(req.clone()),
        })
        .unwrap();
    assert_eq!(serde_json::to_value(engine.snapshot()).unwrap(), before);
    assert_eq!(engine.mutation_revision(), revision);
    assert_eq!(engine.undo_available(), undo);
    let expected = serde_json::to_value(result.rt_primer_redesign.unwrap()).unwrap();
    let journal = engine.journal_len();
    let command = parse_shell_line(&format!(
        "primers redesign-rt {}",
        quote_shell_arg(&serde_json::to_string(&req).unwrap())
    ))
    .unwrap();
    assert!(!command.is_state_mutating());
    let actual = execute_shell_command(&mut engine, &command).unwrap();
    assert!(!actual.state_changed);
    assert_eq!(actual.output["result"]["rt_primer_redesign"], expected);
    assert_eq!(engine.journal_len(), journal);
    let td = tempfile::tempdir().unwrap();
    let path = td.path().join("redesign request's input.json");
    std::fs::write(&path, serde_json::to_vec(&req).unwrap()).unwrap();
    let command = parse_shell_line(&format!(
        "primers redesign-rt {}",
        quote_shell_arg(&format!("@{}", path.display()))
    ))
    .unwrap();
    assert_eq!(
        execute_shell_command(&mut engine, &command).unwrap().output["result"]["rt_primer_redesign"],
        expected
    );
    assert_eq!(engine.journal_len(), journal);
    assert_eq!(serde_json::to_value(engine.snapshot()).unwrap(), before);
}
