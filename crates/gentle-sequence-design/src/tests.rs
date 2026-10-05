//! Hand-crafted synthetic fixtures; no natural gene is represented by these DNA strings.
//! Recreate by concatenating the literal flank/codon strings in each test.
//! The complete standard mapping below is independently declared in TCAG order;
//! coupled-edit oracles enumerate their explicit four strings, not solver choices.

use super::*;

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
