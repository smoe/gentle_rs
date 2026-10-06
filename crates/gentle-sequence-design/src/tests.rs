//! Hand-crafted synthetic fixtures; no natural gene is represented by these DNA strings.
//! Recreate by concatenating the literal flank/codon strings in each test.
//! The complete standard mapping below is independently declared in TCAG order;
//! coupled-edit oracles enumerate their explicit four strings, not solver choices.

use super::*;

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
        assert_eq!(solve(&request, |_| false).status, Status::Invalid);
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
    assert_eq!(solve(&empty, |_| false).status, Status::Invalid);
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
    let mut request = input(
        &format!("ATG{}TAA", "GCT".repeat(1000)),
        &format!("M{}", "A".repeat(1000)),
        &["GAATTC"],
    );
    request.max_evaluations = MAX_EVALUATIONS;
    let without_gc = solve(&request, |_| false).effective_evaluation_budget;
    request.gc_content = gc_bounds(0, 10000);
    let with_gc = solve(&request, |_| false).effective_evaluation_budget;
    assert!(with_gc < without_gc && with_gc > 0);
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
