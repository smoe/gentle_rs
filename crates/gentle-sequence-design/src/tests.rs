//! Synthetic fixtures; no natural gene is represented by these DNA strings.
//! Recreate from literal flank/codon strings or the fixed-seed LCG in the
//! differential test; it uses explicit L/R synonym lists and no external inputs.
//! The complete standard mapping below is independently declared in TCAG order;
//! coupled-edit oracles enumerate their explicit four strings, not solver choices.

use super::*;

fn window_bounds(
    window_bp: usize,
    min_basis_points: u16,
    max_basis_points: u16,
) -> Option<WindowGcBounds> {
    Some(WindowGcBounds {
        window_bp,
        min_basis_points,
        max_basis_points,
    })
}

#[test]
fn window_gc_checks_all_overlaps_non_codon_edges_and_excludes_flanks() {
    // Literal synthetic MFK: global GC alone favours the F edit; local windows
    // require the K edit instead. G/C-rich flanks must not alter these facts.
    let mut request = input("GGGATGTTTAAATAACCC", "MFK", &[]);
    request.cds = Interval { start: 3, end: 15 };
    request.gc_content = gc_bounds(1666, 1667);
    request.gc_window = window_bounds(6, 1666, 1667);
    for strategy in [
        SearchStrategy::FullEnumeration,
        SearchStrategy::ConflictDirected,
    ] {
        let result = solve_with_strategy(&request, strategy, |_| false);
        assert_eq!(result.sequence.as_deref(), Some("GGGATGTTTAAGTAACCC"));
        assert!(result.minimum_edits_proven);
        let facts = assess_window_gc(&request, result.sequence.as_deref())
            .unwrap()
            .unwrap();
        assert_eq!(
            (
                facts.window_bp,
                facts.window_count,
                facts.minimum_gc_bases,
                facts.maximum_gc_bases
            ),
            (6, 7, 1, 1)
        );
        assert_eq!(
            facts.input.iter().map(|r| r.gc_bases).collect::<Vec<_>>(),
            [1, 1, 1, 0, 0, 0, 0]
        );
        for (offset, row) in facts.output.unwrap().iter().enumerate() {
            assert_eq!(
                row.interval,
                Interval {
                    start: 3 + offset,
                    end: 9 + offset
                }
            );
            assert_eq!(row.gc_bases, 1);
            assert!(row.satisfies_bounds);
        }
        assert_eq!(
            validate_output(&request, "GGGATGTTCAAATAACCC").unwrap_err(),
            "gc_window_bounds_not_satisfied"
        );
    }
    request.gc_window = window_bounds(5, 0, 10000);
    let facts = assess_window_gc(&request, None).unwrap().unwrap();
    assert_eq!(facts.window_count, 8);
    assert_eq!(
        facts.input.last().unwrap().interval,
        Interval { start: 10, end: 15 }
    );
    assert!(facts.output.is_none());
    request.gc_window = window_bounds(12, 1666, 1667);
    let result = solve(&request, |_| false);
    assert_eq!(result.sequence.as_deref(), Some("GGGATGTTCAAATAACCC"));
    assert_eq!(
        assess_window_gc(&request, result.sequence.as_deref())
            .unwrap()
            .unwrap()
            .window_count,
        1
    );
}

#[test]
fn window_gc_invalid_empty_frozen_and_cancelled_are_distinct() {
    let mut request = input("ATGTTTAAATAA", "MFK", &[]);
    for bounds in [
        (0, 0, 10000),
        (13, 0, 10000),
        (6, 5001, 5000),
        (6, 0, 10001),
    ] {
        request.gc_window = window_bounds(bounds.0, bounds.1, bounds.2);
        assert_eq!(solve(&request, |_| false).status, Status::Invalid);
        assert!(assess_window_gc(&request, None).is_err());
    }
    for (bounds, reason) in [
        ((6, 1667, 1667), "gc_window_integer_count_bounds_empty"),
        (
            (1, 10000, 10000),
            "gc_window_violation_entirely_in_frozen_bases",
        ),
        ((1, 0, 0), "gc_window_violation_entirely_in_frozen_bases"),
    ] {
        request.gc_window = window_bounds(bounds.0, bounds.1, bounds.2);
        for strategy in [
            SearchStrategy::FullEnumeration,
            SearchStrategy::ConflictDirected,
        ] {
            let result = solve_with_strategy(&request, strategy, |_| false);
            assert_eq!(result.status, Status::ProvenInfeasible);
            assert_eq!(result.evaluated_candidates, 0);
            assert_eq!(result.reason, reason);
            assert!(
                assess_window_gc(&request, None)
                    .unwrap()
                    .unwrap()
                    .output
                    .is_none()
            );
            assert_eq!(
                solve_with_strategy(&request, strategy, |_| true).status,
                Status::Cancelled
            );
        }
    }
    request.gc_window = window_bounds(6, 1666, 1667);
    request.protected = vec![Interval { start: 8, end: 9 }];
    assert_eq!(solve(&request, |_| false).status, Status::ProvenInfeasible);
    request.protected.clear();
    request.max_evaluations = 1;
    for strategy in [
        SearchStrategy::FullEnumeration,
        SearchStrategy::ConflictDirected,
    ] {
        let result = solve_with_strategy(&request, strategy, |_| false);
        assert_eq!(result.status, Status::SearchExhausted);
        assert!(result.sequence.is_none() && !result.minimum_edits_proven);
        assert!(
            assess_window_gc(&request, result.sequence.as_deref())
                .unwrap()
                .unwrap()
                .output
                .is_none()
        );
        request.max_evaluations = 4096;
        let result = solve_with_strategy(&request, strategy, |n| n == 2);
        assert_eq!(result.status, Status::Cancelled);
        assert!(result.sequence.is_none() && !result.minimum_edits_proven);
        request.max_evaluations = 1;
    }
}

#[test]
fn window_gc_partial_codon_direction_is_not_whole_codon_direction() {
    // CGT -> AGG is GC-neutral globally, but lowers its first base and raises its
    // third. Such a synonym must remain reachable when a partial codon conflicts.
    assert!(gc_window::direction_capable(
        6,
        b"CGT",
        &[*b"AGG"],
        Interval { start: 6, end: 7 },
        GcViolation::AboveMaximum
    ));
    assert!(gc_window::direction_capable(
        6,
        b"CGT",
        &[*b"AGG"],
        Interval { start: 8, end: 9 },
        GcViolation::BelowMinimum
    ));
    assert!(!gc::direction_capable(
        b"CGT",
        &[*b"AGG"],
        GcViolation::BelowMinimum
    ));
    assert!(!gc_window::direction_capable(
        6,
        b"CGT",
        &[*b"AGG"],
        Interval { start: 9, end: 12 },
        GcViolation::BelowMinimum
    ));
}

#[test]
fn window_gc_both_solvers_match_independent_36_state_joint_oracle() {
    // Literal standard-code L/R synonyms, not produced by the solver. Count every
    // byte window independently using cross multiplication, not its integer helper.
    let source = "CATGCTACGTTAAG";
    let variants = ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"]
        .into_iter()
        .flat_map(|l| {
            ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"]
                .into_iter()
                .map(move |r| format!("CATG{l}{r}TAAG"))
        })
        .collect::<Vec<_>>();
    for window_bp in [1, 2, 4, 5, 6, 7, 8, 12] {
        for (min, max) in [
            (0, 0),
            (0, 2500),
            (0, 5000),
            (1428, 4286),
            (1666, 1667),
            (2500, 7500),
            (0, 10000),
        ] {
            for global in [None, gc_bounds(2500, 3334)] {
                for patterns in [vec![], vec!["CTA", "CGT"], vec!["CTACGC", "CTGCGT"]] {
                    for protected in [vec![], vec![Interval { start: 6, end: 7 }]] {
                        let mut request = input(source, "MLR", &patterns);
                        request.cds = Interval { start: 1, end: 13 };
                        request.gc_content = global;
                        request.gc_window = window_bounds(window_bp, min, max);
                        request.protected = protected;
                        for motif in &mut request.motifs {
                            motif.strand = Strand::Forward;
                        }
                        let valid = |dna: &str| {
                            let cds = &dna.as_bytes()[1..13];
                            cds.windows(window_bp).all(|window| {
                                let gc = window.iter().filter(|b| b"GC".contains(b)).count();
                                gc * 10000 >= usize::from(min) * window_bp
                                    && gc * 10000 <= usize::from(max) * window_bp
                            }) && global.is_none_or(|g| {
                                let gc = cds.iter().filter(|b| b"GC".contains(b)).count();
                                gc * 10000 >= usize::from(g.min_basis_points) * 12
                                    && gc * 10000 <= usize::from(g.max_basis_points) * 12
                            }) && patterns.iter().all(|p| !dna.contains(p))
                                && request
                                    .protected
                                    .iter()
                                    .all(|p| dna[p.start..p.end] == source[p.start..p.end])
                        };
                        let expected = variants.iter().filter(|dna| valid(dna)).min_by_key(|dna| {
                            (
                                dna.bytes()
                                    .zip(source.bytes())
                                    .filter(|(a, b)| a != b)
                                    .count(),
                                dna.as_str(),
                            )
                        });
                        for strategy in [
                            SearchStrategy::FullEnumeration,
                            SearchStrategy::ConflictDirected,
                        ] {
                            let result = solve_with_strategy(&request, strategy, |_| false);
                            assert_eq!(
                                result.sequence.as_ref(),
                                expected,
                                "{window_bp} {min} {max} {global:?} {patterns:?} {strategy:?}"
                            );
                            assert_eq!(
                                result.status,
                                if expected.is_some() {
                                    Status::Feasible
                                } else {
                                    Status::ProvenInfeasible
                                }
                            );
                            if expected.is_some() {
                                assert!(result.minimum_edits_proven);
                            }
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn window_gc_work_and_report_are_bounded_by_every_requested_window() {
    let dna = format!("ATG{}TAA", "GCT".repeat(1000));
    let mut request = input(&dna, &format!("M{}", "A".repeat(1000)), &[]);
    request.gc_window = window_bounds(100, 0, 10000);
    request.max_evaluations = MAX_EVALUATIONS;
    let length = dna.len() as u64;
    let count = length - 100 + 1;
    let result = solve(&request, |_| false);
    assert_eq!(
        result.effective_evaluation_budget,
        MAX_TOTAL_MOTIF_WORK / (length * 4 + length + 1 + count)
    );
    assert_eq!(
        assess_window_gc(&request, result.sequence.as_deref())
            .unwrap()
            .unwrap()
            .input
            .len(),
        count as usize
    );
}

fn gc_bounds(min: u16, max: u16) -> Option<GcBounds> {
    Some(GcBounds {
        min_basis_points: min,
        max_basis_points: max,
    })
}

#[test]
fn gc_count_bounds_are_integer_inclusive_and_empty_windows_are_infeasible() {
    // Synthetic M + five arginines, 21 bp; literal variants change only Arg.
    let mut request = input("ATGCGTCGTCGTCGTCGTTAA", "MRRRRR", &[]);
    request.gc_content = gc_bounds(4761, 5239);
    let assessment = assess_gc(&request, None).unwrap().unwrap();
    assert_eq!(assessment.denominator_bases, 21);
    assert_eq!(
        (assessment.minimum_gc_bases, assessment.maximum_gc_bases),
        (10, 11)
    );
    for (first, second, expected_count, accepted) in [
        ("AGA", "CGT", 10, true),
        ("CGT", "CGT", 11, true),
        ("AGA", "AGA", 9, false),
        ("CGC", "CGT", 12, false),
    ] {
        let dna = format!("ATG{first}{second}CGTCGTCGTTAA");
        assert_eq!(
            dna.bytes().filter(|b| b"GC".contains(b)).count(),
            expected_count
        );
        assert_eq!(validate_output(&request, &dna).is_ok(), accepted);
    }
    request.gc_content = gc_bounds(4762, 5238);
    let assessment = assess_gc(&request, None).unwrap().unwrap();
    assert_eq!(
        (assessment.minimum_gc_bases, assessment.maximum_gc_bases),
        (11, 10)
    );
    for strategy in [
        SearchStrategy::FullEnumeration,
        SearchStrategy::ConflictDirected,
    ] {
        let result = solve_with_strategy(&request, strategy, |_| false);
        assert_eq!(result.status, Status::ProvenInfeasible);
        assert_eq!(result.evaluated_candidates, 0);
        assert_eq!(
            result.reason,
            "gc_count_window_empty_for_declared_cds_length"
        );
        assert_eq!(
            solve_with_strategy(&request, strategy, |_| true).status,
            Status::Cancelled
        );
    }
    for (min, max) in [(5001, 5000), (0, 10001)] {
        request.gc_content = gc_bounds(min, max);
        for strategy in [
            SearchStrategy::FullEnumeration,
            SearchStrategy::ConflictDirected,
        ] {
            assert_eq!(
                solve_with_strategy(&request, strategy, |_| false).status,
                Status::Invalid
            );
        }
    }
}

#[test]
fn gc_only_design_and_coupled_constraints_are_fully_validated() {
    // Four literal standard-code MFK variants, not a natural gene.
    for (source, motifs, bounds, expected, edit_count) in [
        ("ATGTTTAAATAA", vec![], (1666, 1667), "ATGTTCAAATAA", 1),
        // Raising GC by changing F first creates this motif; changing K works.
        (
            "ATGTTTAAATAA",
            vec!["TTCAAA"],
            (1666, 2500),
            "ATGTTTAAGTAA",
            1,
        ),
        // Removing AAG loses GC; the coupled F edit restores the requested count.
        (
            "ATGTTTAAGTAA",
            vec!["AAGTAA"],
            (1666, 1667),
            "ATGTTCAAATAA",
            2,
        ),
        // GC above maximum with no motif must branch downward, not falsely fail.
        ("ATGTTCAAGTAA", vec![], (833, 834), "ATGTTTAAATAA", 2),
    ] {
        let mut request = input(source, "MFK", &motifs);
        request.gc_content = gc_bounds(bounds.0, bounds.1);
        for strategy in [
            SearchStrategy::FullEnumeration,
            SearchStrategy::ConflictDirected,
        ] {
            let result = solve_with_strategy(&request, strategy, |_| false);
            assert_eq!(result.status, Status::Feasible, "{source} {strategy:?}");
            assert_eq!(result.sequence.as_deref(), Some(expected));
            assert_eq!(result.edits.len(), edit_count);
            assert!(result.optimization_complete && result.minimum_edits_proven);
            assert_eq!(validate_output(&request, expected).unwrap(), result.edits);
            let facts = assess_gc(&request, result.sequence.as_deref())
                .unwrap()
                .unwrap();
            assert!(facts.output.unwrap().satisfies_bounds);
        }
    }
    let empty = input("ATGTTTAAATAA", "MFK", &[]);
    for strategy in [
        SearchStrategy::FullEnumeration,
        SearchStrategy::ConflictDirected,
    ] {
        let result = solve_with_strategy(&empty, strategy, |_| false);
        assert_eq!(result.status, Status::Invalid);
        assert_eq!(
            result.reason,
            "at_least_one_motif_or_explicit_gc_bounds_required"
        );
        assert!(result.sequence.is_none());
        assert_eq!(result.evaluated_candidates, 0);
    }
}

#[test]
fn gc_extrema_include_frozen_bases_and_exclude_flanks_without_proving_feasibility() {
    let mut request = input("ATGTTTAAATAA", "MFK", &[]);
    for bounds in [(0, 0), (10000, 10000)] {
        request.gc_content = gc_bounds(bounds.0, bounds.1);
        for strategy in [
            SearchStrategy::FullEnumeration,
            SearchStrategy::ConflictDirected,
        ] {
            let result = solve_with_strategy(&request, strategy, |_| false);
            assert_eq!(result.status, Status::ProvenInfeasible);
            assert_eq!(result.evaluated_candidates, 0);
            assert!(result.reason.starts_with("gc_bounds_unreachable"));
        }
    }
    request.gc_content = gc_bounds(1666, 1667);
    request.protected = vec![Interval { start: 0, end: 12 }];
    assert_eq!(solve(&request, |_| false).evaluated_candidates, 0);
    assert_eq!(solve(&request, |_| false).status, Status::ProvenInfeasible);
    request.protected.clear();
    request.sequence = "GGGATGTTTAAATAACCC".into();
    request.cds = Interval { start: 3, end: 15 };
    let result = solve_with_strategy(&request, SearchStrategy::ConflictDirected, |_| false);
    assert_eq!(result.sequence.as_deref(), Some("GGGATGTTCAAATAACCC"));
    let facts = assess_gc(&request, result.sequence.as_deref())
        .unwrap()
        .unwrap();
    assert_eq!(facts.denominator_bases, 12);
    assert_eq!(facts.input.gc_bases, 1);
    assert_eq!(facts.output.unwrap().gc_bases, 2);
    assert_eq!(
        (
            facts.reachable_minimum_gc_bases,
            facts.reachable_maximum_gc_bases
        ),
        (1, 3)
    );
    // GC reachability overlaps, but these literal motifs exclude every variant.
    request.motifs = ["TTTAAA", "TTCAAA", "TTTAAG", "TTCAAG"]
        .map(|pattern| Motif {
            pattern: pattern.into(),
            strand: Strand::Forward,
        })
        .to_vec();
    let result = solve(&request, |_| false);
    assert_eq!(result.status, Status::ProvenInfeasible);
    assert!(result.evaluated_candidates > 0);
    assert_eq!(
        result.reason,
        "complete_enumeration_found_no_feasible_variant"
    );
}

#[test]
fn gc_conflict_branching_agrees_with_unrestricted_literal_eight_state_oracle() {
    // Independent unrestricted F/K/F space: no core synonym generation or GC conversion.
    let variants = ["TTT", "TTC"]
        .into_iter()
        .flat_map(|f| {
            ["AAA", "AAG"].into_iter().flat_map(move |k| {
                ["TTT", "TTC"]
                    .into_iter()
                    .map(move |last| format!("CATG{f}{k}{last}TAAG"))
            })
        })
        .collect::<Vec<_>>();
    let original = &variants[0];
    for min in [0, 666, 1333, 2000, 2667, 4000] {
        for max in [0, 667, 1334, 2000, 2667, 10000] {
            if min > max {
                continue;
            }
            for patterns in [
                vec![],
                vec!["TTCAAA"],
                vec!["AAGTAA"],
                vec!["TTTAAA", "TTCAAA", "TTTAAG"],
            ] {
                for protected in [
                    vec![],
                    vec![Interval { start: 6, end: 7 }],
                    vec![Interval { start: 12, end: 13 }],
                ] {
                    let mut request = input(original, "MFKF", &patterns);
                    request.cds = Interval { start: 1, end: 16 };
                    request.protected = protected;
                    request.gc_content = gc_bounds(min, max);
                    for motif in &mut request.motifs {
                        motif.strand = Strand::Forward;
                    }
                    let valid = |dna: &str| {
                        let count = dna[1..16].bytes().filter(|b| b"GC".contains(b)).count();
                        count * 10000 >= usize::from(min) * 15
                            && count * 10000 <= usize::from(max) * 15
                            && patterns.iter().all(|p| !dna.contains(p))
                            && request
                                .protected
                                .iter()
                                .all(|p| dna[p.start..p.end] == original[p.start..p.end])
                    };
                    let expected = variants.iter().filter(|dna| valid(dna)).min_by_key(|dna| {
                        (
                            dna.bytes()
                                .zip(original.bytes())
                                .filter(|(a, b)| a != b)
                                .count(),
                            dna.as_str(),
                        )
                    });
                    for strategy in [
                        SearchStrategy::FullEnumeration,
                        SearchStrategy::ConflictDirected,
                    ] {
                        request.max_evaluations = 4096;
                        let result = solve_with_strategy(&request, strategy, |_| false);
                        assert_eq!(
                            result.sequence.as_ref(),
                            expected,
                            "{min} {max} {patterns:?} {strategy:?}"
                        );
                        assert_eq!(
                            result.status,
                            if expected.is_some() {
                                Status::Feasible
                            } else {
                                Status::ProvenInfeasible
                            }
                        );
                        if expected.is_some() {
                            assert!(result.minimum_edits_proven);
                        }
                        request.max_evaluations = 1;
                        let stopped = solve_with_strategy(&request, strategy, |_| false);
                        if stopped.status == Status::Feasible {
                            assert!(valid(stopped.sequence.as_deref().unwrap()));
                        }
                        if stopped.status == Status::ProvenInfeasible {
                            assert!(expected.is_none());
                        }
                        if stopped.status == Status::SearchExhausted {
                            assert!(!stopped.optimization_complete);
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn gc_conflicts_with_mixed_direction_synonyms_match_unrestricted_oracle() {
    // Synthetic MLR: literal six-way L/R sets include raising, neutral and
    // lowering assignments. The oracle uses neither core choices nor rounding.
    let source = "ATGCTACGTTAA";
    let variants = ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"]
        .into_iter()
        .flat_map(|l| {
            ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"]
                .into_iter()
                .map(move |r| format!("ATG{l}{r}TAA"))
        })
        .collect::<Vec<_>>();
    for (min, max) in [(833, 1667), (2500, 2500), (3333, 3334), (4166, 5000)] {
        for patterns in [vec![], vec!["CTA", "CGT"], vec!["CTACGC", "CTGCGT"]] {
            let mut request = input(source, "MLR", &patterns);
            request.gc_content = gc_bounds(min, max);
            for motif in &mut request.motifs {
                motif.strand = Strand::Forward;
            }
            let valid = |dna: &str| {
                let gc = dna.bytes().filter(|b| b"GC".contains(b)).count();
                gc * 10000 >= usize::from(min) * 12
                    && gc * 10000 <= usize::from(max) * 12
                    && patterns.iter().all(|p| !dna.contains(p))
            };
            let expected = variants.iter().filter(|dna| valid(dna)).min_by_key(|dna| {
                (
                    dna.bytes()
                        .zip(source.bytes())
                        .filter(|(a, b)| a != b)
                        .count(),
                    dna.as_str(),
                )
            });
            for strategy in [
                SearchStrategy::FullEnumeration,
                SearchStrategy::ConflictDirected,
            ] {
                let result = solve_with_strategy(&request, strategy, |_| false);
                assert_eq!(
                    result.sequence.as_ref(),
                    expected,
                    "{min} {max} {patterns:?} {strategy:?}"
                );
                assert!(result.optimization_complete);
                assert_eq!(result.minimum_edits_proven, expected.is_some());
                assert_eq!(
                    result.status,
                    if expected.is_some() {
                        Status::Feasible
                    } else {
                        Status::ProvenInfeasible
                    }
                );
            }
        }
    }
}

#[test]
fn gc_budget_and_cancellation_never_publish_output_measurements_or_proofs() {
    let mut request = input("ATGTTTAAATAA", "MFK", &[]);
    request.gc_content = gc_bounds(1666, 1667);
    request.max_evaluations = 1;
    for strategy in [
        SearchStrategy::FullEnumeration,
        SearchStrategy::ConflictDirected,
    ] {
        let result = solve_with_strategy(&request, strategy, |_| false);
        assert_eq!(result.status, Status::SearchExhausted);
        let facts = assess_gc(&request, result.sequence.as_deref())
            .unwrap()
            .unwrap();
        assert!(!facts.input.satisfies_bounds && facts.output.is_none());
        request.max_evaluations = 4096;
        // Node two is already GC-feasible; cancellation still removes everything.
        let result = solve_with_strategy(&request, strategy, |evaluated| evaluated == 2);
        assert_eq!(result.status, Status::Cancelled);
        assert!(result.sequence.is_none() && result.edits.is_empty());
        assert!(!result.minimum_edits_proven && !result.optimization_complete);
        assert!(
            assess_gc(&request, result.sequence.as_deref())
                .unwrap()
                .unwrap()
                .output
                .is_none()
        );
        request.max_evaluations = 1;
    }
}

#[test]
fn gc_work_budget_accounts_exactly_for_cds_pass_and_motif_work() {
    // Synthetic alanine repeat with fixed flanks; literal lengths independently
    // distinguish whole-DNA work from the added CDS-only GC pass, not runtime.
    let mut request = input(
        &format!(
            "{}ATG{}TAA{}",
            "C".repeat(100),
            "GCT".repeat(1000),
            "G".repeat(100)
        ),
        &format!("M{}", "A".repeat(1000)),
        &["GAATTC"],
    );
    request.cds = Interval {
        start: 100,
        end: 3106,
    };
    request.max_evaluations = MAX_EVALUATIONS;
    assert_eq!(request.sequence.len(), 3206);
    let dna_work = 3206 * 4;
    // GAATTC is palindromic; requesting both strands performs one motif pass.
    let motif_work = (3206 - 6 + 1) * 6;
    let cds_gc_work = 3006;
    for strategy in [
        SearchStrategy::FullEnumeration,
        SearchStrategy::ConflictDirected,
    ] {
        let result = solve_with_strategy(&request, strategy, |_| false);
        assert_eq!(result.status, Status::Feasible);
        assert_eq!(
            result.effective_evaluation_budget,
            MAX_TOTAL_MOTIF_WORK / (dna_work + motif_work)
        );
        let mut gc_request = request.clone();
        gc_request.gc_content = gc_bounds(0, 10000);
        let result = solve_with_strategy(&gc_request, strategy, |_| false);
        assert_eq!(result.status, Status::Feasible);
        assert_eq!(
            result.effective_evaluation_budget,
            MAX_TOTAL_MOTIF_WORK / (dna_work + motif_work + cds_gc_work)
        );
        gc_request.motifs.clear();
        let result = solve_with_strategy(&gc_request, strategy, |_| false);
        assert_eq!(result.status, Status::Feasible);
        assert_eq!(
            result.effective_evaluation_budget,
            MAX_TOTAL_MOTIF_WORK / (dna_work + cds_gc_work)
        );
    }
}

fn code() -> GeneticCode {
    let residues = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    let mut codons = BTreeMap::new();
    let mut i = 0;
    for a in b"TCAG" {
        for b in b"TCAG" {
            for c in b"TCAG" {
                codons.insert([*a, *b, *c], residues[i]);
                i += 1;
            }
        }
    }
    GeneticCode { id: 1, codons }
}

fn input(sequence: &str, protein: &str, patterns: &[&str]) -> Input {
    Input {
        sequence: sequence.into(),
        cds: Interval {
            start: 0,
            end: sequence.len(),
        },
        protein: protein.into(),
        code: code(),
        protected: vec![],
        motifs: patterns
            .iter()
            .map(|p| Motif {
                pattern: (*p).into(),
                strand: Strand::Both,
            })
            .collect(),
        max_evaluations: 4096,
        gc_content: None,
        gc_window: None,
    }
}

#[test]
fn coupled_edit_oracle_defeats_strict_greedy_and_proves_minimum() {
    let request = input("ATGTTTAAATAA", "MFK", &["TTTAAA", "TTCAAA", "TTTAAG"]);
    let variants = [
        "ATGTTTAAATAA",
        "ATGTTCAAATAA",
        "ATGTTTAAGTAA",
        "ATGTTCAAGTAA",
    ];
    let violations = variants
        .iter()
        .map(|s| {
            request
                .motifs
                .iter()
                .filter(|m| s.as_bytes().windows(6).any(|w| w == m.pattern.as_bytes()))
                .count()
        })
        .collect::<Vec<_>>();
    assert_eq!(violations, [1, 1, 1, 0]);
    let result = solve(&request, |_| false);
    assert_eq!(result.status, Status::Feasible);
    assert_eq!(result.sequence.as_deref(), Some(variants[3]));
    assert_eq!(result.edits.len(), 2);
    assert_eq!(result.search_space, Some(4));
    assert_eq!(result.evaluated_candidates, 4);
    assert!(result.minimum_edits_proven && result.optimization_complete);
    assert_eq!(
        validate_output(&request, variants[3]).unwrap(),
        result.edits
    );
}

#[test]
fn bounded_results_agree_with_independent_four_sequence_oracle() {
    // Complete standard-code synonymous space for this MFK insert. Neither
    // variant generation nor this literal matcher calls the core evaluator.
    let variants = [
        "ATGTTTAAATAA",
        "ATGTTCAAATAA",
        "ATGTTTAAGTAA",
        "ATGTTCAAGTAA",
    ];
    let forbidden = |dna: &str, pattern: &str, strand: Strand| {
        let reverse = pattern
            .bytes()
            .rev()
            .map(|b| match b {
                b'A' => 'T',
                b'C' => 'G',
                b'G' => 'C',
                b'T' => 'A',
                _ => panic!("oracle uses literal patterns"),
            })
            .collect::<String>();
        (strand != Strand::Reverse && dna.contains(pattern))
            || (strand != Strand::Forward && dna.contains(&reverse))
    };
    for patterns in [
        vec!["TTTAAA"],
        vec!["TTTAAA", "TTCAAA", "TTTAAG"],
        vec!["TTTAAA", "TTCAAA", "TTTAAG", "TTCAAG"],
    ] {
        for direction in [Strand::Forward, Strand::Reverse, Strand::Both] {
            for protected in [vec![], vec![Interval { start: 5, end: 6 }]] {
                let mut request = input(variants[0], "MFK", &patterns);
                request.protected = protected;
                for motif in &mut request.motifs {
                    motif.strand = direction;
                }
                let expected = variants
                    .iter()
                    .copied()
                    .filter(|dna| {
                        request
                            .protected
                            .iter()
                            .all(|p| dna[p.start..p.end] == variants[0][p.start..p.end])
                            && patterns.iter().all(|p| !forbidden(dna, p, direction))
                    })
                    .min_by_key(|dna| {
                        (
                            dna.bytes()
                                .zip(variants[0].bytes())
                                .filter(|(a, b)| a != b)
                                .count(),
                            *dna,
                        )
                    });
                for budget in 1..=4 {
                    request.max_evaluations = budget;
                    for strategy in [
                        SearchStrategy::FullEnumeration,
                        SearchStrategy::ConflictDirected,
                    ] {
                        let result = solve_with_strategy(&request, strategy, |_| false);
                        match result.status {
                            Status::Feasible => {
                                let dna = result.sequence.as_deref().unwrap();
                                assert!(variants.contains(&dna));
                                assert!(patterns.iter().all(|p| !forbidden(dna, p, direction)));
                                assert!(expected.is_some());
                                if result.minimum_edits_proven {
                                    assert_eq!(Some(dna), expected);
                                }
                            }
                            Status::ProvenInfeasible => assert!(expected.is_none()),
                            Status::SearchExhausted => assert!(!result.optimization_complete),
                            status => panic!("unexpected {status:?}"),
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn conflict_search_matches_an_independent_degenerate_motif_oracle() {
    // Eight explicitly constructed F/K/F variants. These declarations and the
    // base membership tests below do not call the core's mapper or matcher.
    let variants = ["TTT", "TTC"]
        .into_iter()
        .flat_map(|first| {
            ["AAA", "AAG"].into_iter().flat_map(move |middle| {
                ["TTT", "TTC"]
                    .into_iter()
                    .map(move |last| format!("CATG{first}{middle}{last}TAA"))
            })
        })
        .collect::<Vec<_>>();
    let contains = |dna: &str, pattern: &str| {
        dna.as_bytes().windows(pattern.len()).any(|window| {
            window
                .iter()
                .zip(pattern.bytes())
                .all(|(base, symbol)| match symbol {
                    b'R' => b"AG".contains(base),
                    b'Y' => b"CT".contains(base),
                    b'N' => b"ACGT".contains(base),
                    _ => *base == symbol,
                })
        })
    };
    let reverse = |pattern: &str| {
        pattern
            .bytes()
            .rev()
            .map(|symbol| match symbol {
                b'A' => 'T',
                b'T' => 'A',
                b'C' => 'G',
                b'G' => 'C',
                b'R' => 'Y',
                b'Y' => 'R',
                b'N' => 'N',
                _ => panic!("unsupported independent oracle symbol"),
            })
            .collect::<String>()
    };
    for patterns in [
        vec!["TTTAAA"],
        vec!["TTYAAA"],
        vec!["TTTAAA", "TTCAAA", "TTTAAG"],
        vec!["AAA", "AAG"],
        vec!["CATGTTT"],
        vec!["TT", "AAR"],
        vec!["NATG"],
    ] {
        for direction in [Strand::Forward, Strand::Reverse, Strand::Both] {
            for protected in [
                vec![],
                vec![Interval { start: 6, end: 7 }],
                vec![Interval { start: 4, end: 13 }],
            ] {
                let mut request = input(&variants[0], "MFKF", &patterns);
                request.cds = Interval {
                    start: 1,
                    end: variants[0].len(),
                };
                request.protected = protected;
                for motif in &mut request.motifs {
                    motif.strand = direction;
                }
                let expected = variants
                    .iter()
                    .filter(|dna| {
                        request
                            .protected
                            .iter()
                            .all(|p| dna[p.start..p.end] == variants[0][p.start..p.end])
                            && patterns.iter().all(|pattern| {
                                (direction == Strand::Reverse || !contains(dna, pattern))
                                    && (direction == Strand::Forward
                                        || !contains(dna, &reverse(pattern)))
                            })
                    })
                    .min_by_key(|dna| {
                        (
                            dna.bytes()
                                .zip(variants[0].bytes())
                                .filter(|(a, b)| a != b)
                                .count(),
                            *dna,
                        )
                    });
                for strategy in [
                    SearchStrategy::FullEnumeration,
                    SearchStrategy::ConflictDirected,
                ] {
                    let result = solve_with_strategy(&request, strategy, |_| false);
                    assert_eq!(
                        result.sequence.as_ref(),
                        expected,
                        "{strategy:?}: {request:?}"
                    );
                    if expected.is_some() {
                        assert_eq!(result.status, Status::Feasible);
                        assert!(result.minimum_edits_proven && result.optimization_complete);
                    } else {
                        assert_eq!(result.status, Status::ProvenInfeasible);
                    }
                }
            }
        }
    }
}

#[test]
fn conflict_search_can_repair_new_violations_outside_the_original_match() {
    let mut request = input("ATGCTGGAATTCTAA", "MLEF", &["GAATTC", "CTGGAG"]);
    request.protected = vec![Interval { start: 9, end: 12 }];
    let expected = solve(&request, |_| false);
    let actual = solve_with_strategy(&request, SearchStrategy::ConflictDirected, |_| false);
    assert_eq!(actual.status, Status::Feasible);
    assert_eq!(actual.sequence, expected.sequence);
    assert_eq!(actual.sequence.as_deref(), Some("ATGCTAGAGTTCTAA"));
    assert!(actual.minimum_edits_proven && actual.optimization_complete);
    assert_eq!(actual.edits.len(), 2);
    assert!(actual.edits.iter().any(|edit| edit.position < 6));
    assert_eq!(
        validate_output(&request, actual.sequence.as_deref().unwrap()).unwrap(),
        actual.edits
    );
}

#[test]
fn conflict_search_replaces_an_equal_cost_lexicographically_larger_first_solution() {
    let mut request = input("ATGAAAAAATAA", "MKK", &["AAAAAA"]);
    request.motifs[0].strand = Strand::Forward;
    request.max_evaluations = 2;
    let first = solve_with_strategy(&request, SearchStrategy::ConflictDirected, |_| false);
    assert_eq!(first.status, Status::Feasible);
    assert_eq!(first.sequence.as_deref(), Some("ATGAAGAAATAA"));
    assert_eq!(first.edits.len(), 1);
    assert!(!first.minimum_edits_proven && !first.optimization_complete);

    request.max_evaluations = 3;
    let complete = solve_with_strategy(&request, SearchStrategy::ConflictDirected, |_| false);
    assert_eq!(complete.status, Status::Feasible);
    assert_eq!(complete.sequence.as_deref(), Some("ATGAAAAAGTAA"));
    assert!(complete.sequence < first.sequence);
    assert_eq!(complete.edits.len(), first.edits.len());
    assert_eq!(complete.evaluated_candidates, 3);
    assert!(complete.minimum_edits_proven && complete.optimization_complete);
    assert_eq!(
        validate_output(&request, complete.sequence.as_deref().unwrap()).unwrap(),
        complete.edits
    );
}

#[test]
fn conflict_search_prunes_infeasible_equal_cost_nodes_without_visiting_their_children() {
    let mut request = input("ATGAAAAAATAA", "MKK", &["AAAAAA", "AAAAAG"]);
    for motif in &mut request.motifs {
        motif.strand = Strand::Forward;
    }
    request.max_evaluations = 3;
    // The first child is feasible with one edit. The second has one edit but
    // remains infeasible; its two-edit child is valid, yet cannot win.
    assert_eq!(validate_output(&request, "ATGAAGAAATAA").unwrap().len(), 1);
    assert_eq!(edits(request.sequence.as_bytes(), b"ATGAAAAAGTAA").len(), 1);
    assert_eq!(
        validate_output(&request, "ATGAAAAAGTAA").unwrap_err(),
        "forbidden_motifs_remain"
    );
    assert_eq!(validate_output(&request, "ATGAAGAAGTAA").unwrap().len(), 2);
    let result = solve_with_strategy(&request, SearchStrategy::ConflictDirected, |_| false);
    assert_eq!(result.status, Status::Feasible);
    assert_eq!(result.sequence.as_deref(), Some("ATGAAGAAATAA"));
    assert_eq!(result.evaluated_candidates, 3);
    assert!(result.minimum_edits_proven && result.optimization_complete);
}

#[test]
fn edited_candidate_match_limit_is_incomplete_not_invalid_or_infeasible() {
    // Synthetic lysine repeats: AAG occurs once; the 15 poly-A patterns have
    // sum(282 - length + 1, length=3..17) = 4095 hits. AAG -> AAA creates
    // a 285-A run with 4140 hits and removes the single AAG match.
    let dna = format!("ATGAAG{}TAA", "AAA".repeat(94));
    let mut request = input(&dna, &format!("M{}", "K".repeat(95)), &["AAG"]);
    request.motifs[0].strand = Strand::Forward;
    request.motifs.extend((3..=17).map(|length| Motif {
        pattern: "A".repeat(length),
        strand: Strand::Forward,
    }));
    request.max_evaluations = 2;
    let changed = format!("ATG{}TAA", "AAA".repeat(95));
    let literal_count = |sequence: &str| {
        request
            .motifs
            .iter()
            .map(|motif| {
                sequence
                    .as_bytes()
                    .windows(motif.pattern.len())
                    .filter(|window| *window == motif.pattern.as_bytes())
                    .count()
            })
            .sum::<usize>()
    };
    assert_eq!(literal_count(&dna), 4096);
    assert_eq!(literal_count(&changed), 4140);
    assert_eq!(edits(dna.as_bytes(), changed.as_bytes()).len(), 1);
    assert_eq!(
        translate(dna.as_bytes(), request.cds, &request.code),
        translate(changed.as_bytes(), request.cds, &request.code)
    );
    assert_eq!(
        validate_output(&request, &changed).unwrap_err(),
        "match_diagnostic_limit_exceeded"
    );
    for strategy in [
        SearchStrategy::FullEnumeration,
        SearchStrategy::ConflictDirected,
    ] {
        let result = solve_with_strategy(&request, strategy, |_| false);
        assert_eq!(result.status, Status::SearchExhausted, "{strategy:?}");
        assert_eq!(result.reason, "match_diagnostic_limit_exceeded");
        assert_eq!(result.initial_matches.len(), MAX_REPORTED_MATCHES);
        assert_eq!(result.evaluated_candidates, 2);
        assert_eq!(result.effective_evaluation_budget, 2);
        assert!(result.sequence.is_none() && result.edits.is_empty());
        assert!(!result.minimum_edits_proven && !result.optimization_complete);
        let cancelled = solve_with_strategy(&request, strategy, |n| n == 2);
        assert_eq!(cancelled.status, Status::Cancelled);
        assert!(cancelled.sequence.is_none() && !cancelled.minimum_edits_proven);

        let mut inadmissible = request.clone();
        inadmissible.sequence = changed.clone();
        let rejected = solve_with_strategy(&inadmissible, strategy, |_| false);
        assert_eq!(rejected.status, Status::Invalid);
        assert_eq!(rejected.reason, "match_diagnostic_limit_exceeded");
        assert_eq!(rejected.evaluated_candidates, 0);
    }
}

#[test]
fn conflict_search_keeps_a_validated_best_after_a_later_match_limit() {
    // M accepts A/C, not G. These distinct 32-base masks therefore have the
    // same hits on this A/G repeat: 15 * (303 - 31) + one GAAAAAG = 4081.
    let dna = format!("ATGAAAAAG{}TAA", "AAA".repeat(101));
    let mut request = input(&dna, &format!("M{}", "K".repeat(103)), &["GAAAAAG"]);
    request.motifs[0].strand = Strand::Forward;
    request.motifs.push(Motif {
        pattern: "A".repeat(32),
        strand: Strand::Forward,
    });
    for position in 0..14 {
        let mut pattern = "A".repeat(32);
        pattern.replace_range(position..position + 1, "M");
        request.motifs.push(Motif {
            pattern,
            strand: Strand::Forward,
        });
    }
    let mut expected = dna.as_bytes().to_vec();
    expected[5] = b'G';
    for index in 0..101 {
        let start = 9 + index * 3;
        if index % 10 == 8 {
            expected[start + 2] = b'G';
        } else {
            request.protected.push(Interval {
                start,
                end: start + 3,
            });
        }
    }
    let expected = String::from_utf8(expected).unwrap();
    assert_eq!(validate_output(&request, &expected).unwrap().len(), 11);
    request.max_evaluations = 12;
    let before = solve_with_strategy(&request, SearchStrategy::ConflictDirected, |_| false);
    assert_eq!(before.status, Status::Feasible);
    assert_eq!(before.sequence.as_deref(), Some(expected.as_str()));
    assert_eq!(
        before.reason,
        "validated_best_candidate_with_incomplete_optimization"
    );

    // The other root branch edits AAG -> AAA, creating 15 * (309 - 31) = 4170
    // matches. The first branch's independently validated best must survive.
    request.max_evaluations = 13;
    let after = solve_with_strategy(&request, SearchStrategy::ConflictDirected, |_| false);
    assert_eq!(after.status, Status::Feasible);
    assert_eq!(after.sequence, before.sequence);
    assert_eq!(after.edits, before.edits);
    assert_eq!(after.initial_matches.len(), 4081);
    assert_eq!(after.evaluated_candidates, 13);
    assert_eq!(
        after.reason,
        "validated_best_candidate_with_match_diagnostic_limit"
    );
    assert!(!after.minimum_edits_proven && !after.optimization_complete);
    let cancelled = solve_with_strategy(&request, SearchStrategy::ConflictDirected, |n| n == 13);
    assert_eq!(cancelled.status, Status::Cancelled);
    assert!(cancelled.sequence.is_none() && cancelled.edits.is_empty());
    assert!(!cancelled.minimum_edits_proven && !cancelled.optimization_complete);
}

#[test]
fn both_strategies_prove_literal_optima_in_large_leucine_arginine_spaces() {
    for count in 4..=6 {
        let mut dna = String::from("ATG");
        let mut expected = String::from("ATG");
        let mut protein = String::from("M");
        for position in 0..count {
            let (codon, best, residue) = if position % 2 == 0 {
                ("CTT", "CTA", 'L')
            } else {
                ("CGT", "CGA", 'R')
            };
            dna.push_str(codon);
            expected.push_str(best);
            protein.push(residue);
        }
        dna.push_str("TAA");
        expected.push_str("TAA");
        let mut request = input(&dna, &protein, &["CTT", "CGT"]);
        request.max_evaluations = MAX_EVALUATIONS;
        for motif in &mut request.motifs {
            motif.strand = Strand::Forward;
        }
        // Each original codon is forbidden, requiring >= count edits. CTA/CGA
        // are the lexicographically first one-edit synonyms and satisfy all motifs.
        for strategy in [
            SearchStrategy::FullEnumeration,
            SearchStrategy::ConflictDirected,
        ] {
            let result = solve_with_strategy(&request, strategy, |_| false);
            assert_eq!(result.status, Status::Feasible, "{count}: {strategy:?}");
            assert_eq!(result.search_space, Some(6u64.pow(count as u32)));
            assert_eq!(result.sequence.as_deref(), Some(expected.as_str()));
            assert_eq!(result.edits.len(), count);
            assert!(result.minimum_edits_proven && result.optimization_complete);
            assert_eq!(validate_output(&request, &expected).unwrap(), result.edits);
        }
    }
}

#[test]
fn fixed_seed_joint_constraints_agree_between_complete_strategies() {
    fn draw(state: &mut u64, upper: usize) -> usize {
        *state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
        ((*state >> 32) as usize) % upper
    }
    let families = [
        ('L', ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"]),
        ('R', ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"]),
    ];
    let mut seed = 0x2026_1015_5eed_u64;
    let mut feasible = 0;
    let mut infeasible = 0;
    for case in 0..256 {
        let count = 2 + case % 3;
        let mut dna = String::from("ATG");
        let mut protein = String::from("M");
        for _ in 0..count {
            let (residue, synonyms) = families[draw(&mut seed, families.len())];
            dna.push_str(synonyms[draw(&mut seed, synonyms.len())]);
            protein.push(residue);
        }
        dna.push_str("TAA");
        let mut request = input(&dna, &protein, &[]);
        request.max_evaluations = MAX_EVALUATIONS;
        let length = 3 + draw(&mut seed, 3);
        let start = 3 + draw(&mut seed, count * 3 - length + 1);
        let mut pattern = dna[start..start + length].to_string();
        if case % 5 == 0 {
            pattern.replace_range(1..2, "N");
        }
        request.motifs.push(Motif {
            pattern,
            strand: [Strand::Forward, Strand::Reverse, Strand::Both][draw(&mut seed, 3)],
        });
        if case % 4 != 0 {
            let start = 3 + draw(&mut seed, count * 3);
            request.protected.push(Interval {
                start,
                end: start + 1,
            });
        }
        if case % 4 == 1 || case % 4 == 3 {
            let min = 1000 + draw(&mut seed, 5) as u16 * 1000;
            request.gc_content = Some(GcBounds {
                min_basis_points: min,
                max_basis_points: min + 2000 + draw(&mut seed, 4) as u16 * 1000,
            });
        }
        if case % 4 >= 2 {
            let min = draw(&mut seed, 5) as u16 * 1000;
            request.gc_window = window_bounds(
                4 + draw(&mut seed, 4),
                min,
                min + 4000 + draw(&mut seed, 3) as u16 * 1000,
            );
        }
        let enumeration = solve(&request, |_| false);
        let directed = solve_with_strategy(&request, SearchStrategy::ConflictDirected, |_| false);
        assert_eq!(
            directed.status, enumeration.status,
            "case {case}: {request:?}"
        );
        assert_eq!(
            directed.sequence, enumeration.sequence,
            "case {case}: {request:?}"
        );
        assert_eq!(directed.edits, enumeration.edits, "case {case}");
        assert_eq!(
            directed.optimization_complete, enumeration.optimization_complete,
            "case {case}"
        );
        assert_eq!(
            directed.minimum_edits_proven, enumeration.minimum_edits_proven,
            "case {case}"
        );
        match enumeration.status {
            Status::Feasible => {
                feasible += 1;
                assert!(enumeration.minimum_edits_proven && enumeration.optimization_complete);
                assert_eq!(
                    validate_output(&request, enumeration.sequence.as_deref().unwrap()).unwrap(),
                    enumeration.edits
                );
            }
            Status::ProvenInfeasible => infeasible += 1,
            other => panic!("case {case} must finish, got {other:?}: {request:?}"),
        }
    }
    // Independently enumerated with literal L/R lists, regex motif matching
    // and integer GC cross-products, not either solver's choices/predicates.
    assert_eq!((feasible, infeasible), (198, 58));
}

#[test]
fn conflict_search_can_finish_without_enumerating_unrelated_synonyms() {
    let dna = format!("ATG{}GAATTCTAA", "GCT".repeat(100));
    let request = input(&dna, &format!("M{}EF", "A".repeat(100)), &["GAATTC"]);
    let mut bounded = request.clone();
    bounded.max_evaluations = 16;
    let enumeration = solve(&bounded, |_| false);
    let directed = solve_with_strategy(&bounded, SearchStrategy::ConflictDirected, |_| false);
    assert_eq!(enumeration.status, Status::Feasible);
    assert_eq!(directed.sequence, enumeration.sequence);
    assert_eq!(directed.search_space, None);
    assert!(!enumeration.minimum_edits_proven && !enumeration.optimization_complete);
    assert!(directed.minimum_edits_proven && directed.optimization_complete);
    assert!(directed.evaluated_candidates < enumeration.evaluated_candidates);
    assert_eq!(
        directed,
        solve_with_strategy(&bounded, SearchStrategy::ConflictDirected, |_| false)
    );
}

#[test]
fn conflict_search_budget_stop_cancellation_and_proofs_remain_distinct() {
    let mut coupled = input("ATGTTTAAATAA", "MFK", &["TTTAAA", "TTCAAA", "TTTAAG"]);
    coupled.max_evaluations = 2;
    let unresolved = solve_with_strategy(&coupled, SearchStrategy::ConflictDirected, |_| false);
    assert_eq!(unresolved.status, Status::SearchExhausted);
    assert!(unresolved.sequence.is_none() && !unresolved.optimization_complete);
    coupled.max_evaluations = 3;
    let feasible = solve_with_strategy(&coupled, SearchStrategy::ConflictDirected, |_| false);
    assert_eq!(feasible.status, Status::Feasible);
    assert!(!feasible.minimum_edits_proven && !feasible.optimization_complete);
    assert_eq!(feasible.edits.len(), 2);
    let cancelled = solve_with_strategy(&coupled, SearchStrategy::ConflictDirected, |n| n == 3);
    assert_eq!(cancelled.status, Status::Cancelled);
    assert!(cancelled.sequence.is_none() && cancelled.edits.is_empty());
    assert!(!cancelled.minimum_edits_proven && !cancelled.optimization_complete);
    coupled.max_evaluations = 4096;
    let complete = solve_with_strategy(&coupled, SearchStrategy::ConflictDirected, |_| false);
    assert_eq!(complete.sequence, feasible.sequence);
    assert!(complete.minimum_edits_proven && complete.optimization_complete);
    let impossible = input("ATGTTTTAA", "MF", &["TTY"]);
    let proof = solve_with_strategy(&impossible, SearchStrategy::ConflictDirected, |_| false);
    assert_eq!(proof.status, Status::ProvenInfeasible);
    assert_eq!(
        proof.reason,
        "complete_conflict_search_found_no_feasible_variant"
    );
    assert_eq!(proof.evaluated_candidates, 2);
    let frozen = input("ATGAAATAA", "MK", &["ATG"]);
    assert_eq!(
        solve_with_strategy(&frozen, SearchStrategy::ConflictDirected, |_| false).reason,
        "forbidden_match_entirely_in_frozen_bases"
    );
}

#[test]
fn budget_failure_is_not_infeasibility_and_cancellation_discards_valid_candidates() {
    let mut request = input("ATGTTTAAATAA", "MFK", &["TTTAAA", "TTCAAA", "TTTAAG"]);
    request.max_evaluations = 3;
    let exhausted = solve(&request, |_| false);
    assert_eq!(exhausted.status, Status::SearchExhausted);
    assert!(!exhausted.optimization_complete);
    assert!(exhausted.sequence.is_none());
    request.max_evaluations = 4;
    let cancelled = solve(&request, |n| n == 4);
    assert_eq!(cancelled.status, Status::Cancelled);
    assert!(cancelled.sequence.is_none() && cancelled.edits.is_empty());
    let mut request = input("ATGGAATTCTAA", "MEF", &["GAATTC"]);
    request.max_evaluations = 2;
    let feasible = solve(&request, |_| false);
    assert_eq!(feasible.status, Status::Feasible);
    assert!(!feasible.optimization_complete && !feasible.minimum_edits_proven);
}

#[test]
fn iupac_single_symbols_match_literal_membership_on_every_strand() {
    // Probe the production matcher directly; these one-base strings are not CDSs.
    let mut cases = 0;
    for (symbol, forward, reverse) in IUPAC_LITERAL_MEMBERSHIPS {
        let pattern = String::from_utf8(vec![symbol]).unwrap();
        for strand in [Strand::Forward, Strand::Reverse, Strand::Both] {
            let mut request = input("ATGAAATAA", "MK", &[&pattern]);
            request.motifs[0].strand = strand;
            let prepared = prepare(&request).unwrap_or_else(|outcome| panic!("{outcome:?}"));
            for base in b"ACGT" {
                let sequence = [*base];
                assert_eq!(
                    prepared.matches(&sequence).unwrap(),
                    literal_iupac_hits(&sequence, &[forward], &[reverse], strand),
                    "{pattern} / {} / {strand:?}",
                    *base as char
                );
                cases += 1;
            }
        }
    }
    assert_eq!(cases, 180);
}

#[test]
fn iupac_symbol_pairs_preserve_reverse_order_overlaps_and_palindromes() {
    let mut cases = 0;
    for (first, first_forward, first_reverse) in IUPAC_LITERAL_MEMBERSHIPS {
        for (second, second_forward, second_reverse) in IUPAC_LITERAL_MEMBERSHIPS {
            let pattern = String::from_utf8(vec![first, second]).unwrap();
            let forward = [first_forward, second_forward];
            let reverse = [second_reverse, first_reverse];
            for strand in [Strand::Forward, Strand::Reverse, Strand::Both] {
                let mut request = input("ATGAAATAA", "MK", &[&pattern]);
                request.motifs[0].strand = strand;
                let prepared = prepare(&request).unwrap_or_else(|outcome| panic!("{outcome:?}"));
                for first_base in b"ACGT" {
                    for second_base in b"ACGT" {
                        let sequence = [
                            b'C',
                            *first_base,
                            *second_base,
                            *first_base,
                            *second_base,
                            b'A',
                        ];
                        assert_eq!(
                            prepared.matches(&sequence).unwrap(),
                            literal_iupac_hits(&sequence, &forward, &reverse, strand),
                            "{pattern} / {sequence:?} / {strand:?}"
                        );
                        cases += 1;
                    }
                }
            }
        }
    }
    assert_eq!(cases, 10_800);
}

// Hand-declared synthetic membership/complement sets, not derived from mask(),
// complement(), Prepared patterns, or the Python oracle. Recreate by these rows.
const IUPAC_LITERAL_MEMBERSHIPS: [(u8, &[u8], &[u8]); 15] = [
    (b'A', b"A", b"T"),
    (b'C', b"C", b"G"),
    (b'G', b"G", b"C"),
    (b'T', b"T", b"A"),
    (b'R', b"AG", b"CT"),
    (b'Y', b"CT", b"AG"),
    (b'S', b"CG", b"CG"),
    (b'W', b"AT", b"AT"),
    (b'K', b"GT", b"AC"),
    (b'M', b"AC", b"GT"),
    (b'B', b"CGT", b"ACG"),
    (b'D', b"AGT", b"ACT"),
    (b'H', b"ACT", b"AGT"),
    (b'V', b"ACG", b"CGT"),
    (b'N', b"ACGT", b"ACGT"),
];

fn literal_iupac_hits(
    sequence: &[u8],
    forward: &[&[u8]],
    reverse: &[&[u8]],
    strand: Strand,
) -> Vec<Match> {
    let mut hits = vec![];
    for (start, window) in sequence.windows(forward.len()).enumerate() {
        let f = window
            .iter()
            .zip(forward)
            .all(|(base, members)| members.contains(base));
        let r = window
            .iter()
            .zip(reverse)
            .all(|(base, members)| members.contains(base));
        let directions: &[Strand] = match (strand, f, r, forward == reverse) {
            (Strand::Forward, true, _, _) => &[Strand::Forward],
            (Strand::Reverse, _, true, _) => &[Strand::Reverse],
            (Strand::Both, true, true, true) => &[Strand::Both],
            (Strand::Both, true, true, false) => &[Strand::Forward, Strand::Reverse],
            (Strand::Both, true, false, _) => &[Strand::Forward],
            (Strand::Both, false, true, _) => &[Strand::Reverse],
            _ => &[],
        };
        for direction in directions {
            hits.push(Match {
                motif_index: 0,
                interval: Interval {
                    start,
                    end: start + forward.len(),
                },
                strand: *direction,
            });
        }
    }
    hits
}

#[test]
fn palindromes_overlap_reverse_iupac_and_cds_flank_boundary_are_evaluated() {
    let eco = input("ATGGAATTCTAA", "MEF", &["GAATTC"]);
    let result = solve(&eco, |_| false);
    assert_eq!(result.initial_matches.len(), 1);
    assert_eq!(result.initial_matches[0].strand, Strand::Both);
    let mut boundary = input("AATGAAATAA", "MK", &["AATGAAA"]);
    boundary.cds = Interval { start: 1, end: 10 };
    assert_eq!(solve(&boundary, |_| false).status, Status::Feasible);
    assert_eq!(
        solve(&input("ATGGAATTCTAA", "MEF", &["GAATYC"]), |_| false).status,
        Status::Feasible
    );
    let reverse = input("ATGAAATAA", "MK", &["TTT"]);
    let result = solve(&reverse, |_| false);
    assert!(
        result
            .initial_matches
            .iter()
            .any(|h| h.strand == Strand::Reverse)
    );
    assert_eq!(result.status, Status::Feasible);
    let overlap = input("ATGAAAAAATAA", "MKK", &["AAA"]);
    let result = solve(&overlap, |_| false);
    assert_eq!(
        result
            .initial_matches
            .iter()
            .filter(|h| h.strand == Strand::Forward)
            .count(),
        4
    );
    assert_eq!(result.status, Status::Feasible);
}

#[test]
fn frozen_flanks_start_stop_and_protected_bases_are_never_changed() {
    let mut request = input("ATGGAATTCTAA", "MEF", &["GAATTC"]);
    request.protected = vec![Interval { start: 3, end: 6 }];
    let result = solve(&request, |_| false);
    assert_eq!(result.sequence.as_deref(), Some("ATGGAATTTTAA"));
    assert_eq!(result.edits[0].position, 8);
    request.protected = vec![Interval { start: 3, end: 9 }];
    let result = solve(&request, |_| false);
    assert_eq!(result.status, Status::ProvenInfeasible);
    assert_eq!(result.reason, "forbidden_match_entirely_in_frozen_bases");
    let mut flank = input("GAATTCATGAAATAA", "MK", &["GAATTC"]);
    flank.cds = Interval { start: 6, end: 15 };
    assert_eq!(solve(&flank, |_| false).status, Status::ProvenInfeasible);
    for (pattern, expected) in [
        ("ATG", "forbidden_match_entirely_in_frozen_bases"),
        ("TAA", "forbidden_match_entirely_in_frozen_bases"),
    ] {
        let result = solve(&input("ATGAAATAA", "MK", &[pattern]), |_| false);
        assert_eq!(result.reason, expected);
        assert_eq!(result.status, Status::ProvenInfeasible);
    }
    assert!(validate_output(&input("ATGAAATAA", "MK", &["GAATTC"]), "ATGAAATAG").is_err());
}

#[test]
fn unchanged_input_is_feasible_and_sound_complete_enumeration_can_fail() {
    let request = input("ATGAAATAA", "MK", &["GAATTC"]);
    let result = solve(&request, |_| false);
    assert_eq!(result.sequence.as_deref(), Some(request.sequence.as_str()));
    assert!(result.edits.is_empty() && result.minimum_edits_proven);
    let impossible = input("ATGAAATAA", "MK", &["AAA", "AAG"]);
    let result = solve(&impossible, |_| false);
    assert_eq!(result.status, Status::ProvenInfeasible);
    assert_eq!(
        result.reason,
        "complete_enumeration_found_no_feasible_variant"
    );
    assert_eq!(result.evaluated_candidates, 2);
}

#[test]
fn admission_and_work_are_bounded_before_search() {
    for request in [
        input("ATGNNNTAA", "MK", &["AAAA"]),
        input("ATGAAATAA", "ML", &["AAAA"]),
        input("ATGTAATAA", "MK", &["AAAA"]),
    ] {
        assert_eq!(solve(&request, |_| false).status, Status::Invalid);
    }
    let mut request = input("ATGAAATAA", "MK", &["A+"]);
    assert_eq!(solve(&request, |_| false).status, Status::Invalid);
    request.motifs[0].pattern = "GAATTC".into();
    request.code.id = 2;
    assert_eq!(solve(&request, |_| false).status, Status::Unsupported);
    request.code.id = 1;
    request.code.codons.remove(b"AAA");
    assert_eq!(solve(&request, |_| false).status, Status::Invalid);
    let dna = format!("ATG{}TAA", "GCT".repeat(1000));
    let mut request = input(
        &dna,
        &format!("M{}", "A".repeat(1000)),
        &["GGNCCN", "AAAAAA", "TTTTTT"],
    );
    request.max_evaluations = MAX_EVALUATIONS;
    let result = solve(&request, |_| false);
    assert_eq!(result.status, Status::Feasible);
    assert!(result.effective_evaluation_budget < MAX_EVALUATIONS);
    assert_eq!(result.search_space, None);
}
