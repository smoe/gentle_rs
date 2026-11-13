//! Hand-crafted synthetic loci, bases and qualities; no cell-line or biological findings.
//! Recreate through the `cdna_comparison_` library-test filter. No fixture files or I/O.

use super::*;

const POS: [usize; 4] = [101, 103, 121, 123];

fn unit(id: &str, pattern: &str) -> CdnaUnit {
    CdnaUnit {
        unit_id: id.into(),
        all_query_views_assessed: true,
        search_declared_complete: false,
        hypotheses: vec![CdnaHypothesis {
            hypothesis_id: format!("mapping_{id}"),
            eligible: true,
            bases: POS
                .iter()
                .zip(pattern.chars())
                .enumerate()
                .map(|(i, (&p, b))| CdnaBaseObservation {
                    reference_position_0based: p,
                    query_position_0based: i,
                    observed_base: b,
                    original_base: b,
                    query_to_reference_strand: '+',
                    phred: Some(40),
                    mapq: Some(60),
                })
                .collect(),
            junctions: vec![(110, 120)],
        }],
    }
}

fn sample(id: &str, patterns: &[&str]) -> CdnaSample {
    CdnaSample {
        sample_id: id.into(),
        label: format!("invented_{id}"),
        source_sha256: sha(id.as_bytes()),
        analysis_context_sha256: sha(format!("context_{id}").as_bytes()),
        units: patterns
            .iter()
            .enumerate()
            .map(|(i, p)| unit(&format!("unit_{i}"), p))
            .collect(),
    }
}

fn request(samples: Vec<CdnaSample>) -> CdnaComparisonRequest {
    let sequence = "A".repeat(40);
    CdnaComparisonRequest {
        schema: CDNA_COMPARISON_REQUEST_SCHEMA.into(),
        gene: "invented_gene".into(),
        reference: CdnaReference {
            assembly: "synthetic-v1".into(),
            contig: "invented_contig".into(),
            start_0based: 100,
            sha256: sha(sequence.as_bytes()),
            sequence,
        },
        transcripts: vec![
            CdnaTranscript {
                transcript_id: "plus".into(),
                strand: '+',
                exons: vec![(100, 110), (120, 130)],
            },
            CdnaTranscript {
                transcript_id: "minus".into(),
                strand: '-',
                exons: vec![(100, 110), (120, 130)],
            },
            CdnaTranscript {
                transcript_id: "short".into(),
                strand: '+',
                exons: vec![(120, 130)],
            },
        ],
        samples,
        supplied_sites: vec![],
        phase_marker_positions: None,
        policy: CdnaComparisonPolicy::default(),
        budgets: CdnaComparisonBudgets::default(),
    }
}

fn run(r: CdnaComparisonRequest) -> CdnaComparisonReport {
    compare_cdna_alleles(r, &mut |_, _| true).unwrap()
}

fn diploid() -> CdnaComparisonRequest {
    request(vec![sample(
        "a",
        &["AAAA", "AAAA", "AAAA", "GGGG", "GGGG", "GGGG"],
    )])
}

/// Invented partial-read cohorts; indices refer to POS, with no reference-filled bases.
fn partial_request(groups: &[(&[usize], &str, usize)]) -> CdnaComparisonRequest {
    let mut s = sample("a", &[]);
    let mut supplied = BTreeSet::new();
    for (group, &(indices, pattern, copies)) in groups.iter().enumerate() {
        assert_eq!(indices.len(), pattern.len());
        for copy in 0..copies {
            let mut u = unit(&format!("partial_{group}_{copy}"), "AAAA");
            u.hypotheses[0].bases = indices
                .iter()
                .zip(pattern.chars())
                .map(|(&index, base)| {
                    supplied.insert(POS[index]);
                    let mut b = u.hypotheses[0].bases[index].clone();
                    b.observed_base = base;
                    b.original_base = base;
                    b
                })
                .collect();
            u.hypotheses[0].junctions.clear();
            s.units.push(u);
        }
    }
    let mut r = request(vec![s]);
    r.supplied_sites = supplied.into_iter().collect();
    r
}

fn vectors(block: &CdnaPhaseBlock) -> BTreeSet<String> {
    block
        .haplotypes
        .iter()
        .map(|h| h.alleles.iter().map(|a| a.observed_base).collect())
        .collect()
}

#[test]
fn cdna_comparison_discovers_snv_candidates_and_preserves_quality_callability() {
    let mut r = diploid();
    let mut missing = unit("quality_missing", "GGGG");
    for b in &mut missing.hypotheses[0].bases {
        b.phred = None;
    }
    r.samples[0].units.push(missing);
    let report = run(r);
    assert_eq!(report.sites.len(), 4);
    assert_eq!(report.discovery_assessed_positions, 4);
    for site in &report.sites {
        assert_eq!(site.candidate_alternate_bases, ['G']);
        let c = &site.samples[0];
        assert_eq!(
            (c.total_units, c.callable_units, c.unassessed_units),
            (7, 6, 1)
        );
        assert_eq!(c.base_counts[&'A'], 3);
        assert_eq!(c.base_counts[&'G'], 3);
    }
    let s = &report.samples[0];
    assert_eq!(s.units_without_complete_search_declaration, 7);
    assert!(
        s.site_calls
            .iter()
            .filter(|c| c.unit_id == "quality_missing")
            .all(|c| c.status == CdnaCallStatus::Unassessed
                && c.base.is_none()
                && c.reasons.contains(&"base_quality_unavailable".into()))
    );
    assert!(
        s.site_observations
            .iter()
            .any(|o| o.unit_id == "quality_missing"
                && o.observation.phred.is_none()
                && o.observation.original_base == 'G')
    );
}

#[test]
fn cdna_comparison_eligible_alternatives_and_incomplete_views_never_become_reference_calls() {
    let mut r = request(vec![sample("a", &["AAAA", "AAAA", "AAAA"])]);
    r.supplied_sites = POS.to_vec();
    let mut conflict = unit("conflict", "AAAA");
    let mut other = unit("other", "GGGG").hypotheses.remove(0);
    other.hypothesis_id = "other".into();
    conflict.hypotheses.push(other);
    r.samples[0].units.push(conflict);
    let mut partial = unit("partial", "AAAA");
    let mut other = partial.hypotheses[0].clone();
    other.hypothesis_id = "partial_other".into();
    other.bases.remove(0);
    partial.hypotheses.push(other);
    r.samples[0].units.push(partial);
    let mut incomplete = unit("incomplete", "AAAA");
    incomplete.all_query_views_assessed = false;
    r.samples[0].units.push(incomplete);
    let report = run(r);
    let calls = &report.samples[0].site_calls;
    assert!(
        calls
            .iter()
            .filter(|c| c.unit_id == "conflict")
            .all(|c| c.status == CdnaCallStatus::Ambiguous && c.possible_bases == ['A', 'G'])
    );
    assert!(calls.iter().any(|c| c.unit_id == "partial"
        && c.position_0based == 101
        && c.status == CdnaCallStatus::Unassessed));
    assert!(
        calls
            .iter()
            .filter(|c| c.unit_id == "incomplete")
            .all(|c| c.status == CdnaCallStatus::Unassessed)
    );
}

#[test]
fn cdna_comparison_nonoptimal_mapping_is_audit_only() {
    let mut r = request(vec![sample("a", &["AAAA"])]);
    r.supplied_sites = POS.to_vec();
    let mut other = unit("other", "GGGG").hypotheses.remove(0);
    other.hypothesis_id = "nonoptimal".into();
    other.eligible = false;
    r.samples[0].units[0].hypotheses.push(other);
    let report = run(r);
    assert!(
        report.samples[0]
            .site_calls
            .iter()
            .all(|c| c.status == CdnaCallStatus::ReferenceSupported)
    );
    assert!(
        report.samples[0]
            .site_observations
            .iter()
            .any(|o| !o.eligible && o.observation.observed_base == 'G')
    );
    assert!(
        report
            .sites
            .iter()
            .all(|s| s.candidate_alternate_bases.is_empty())
    );
}

#[test]
fn cdna_comparison_reverse_query_and_transcript_projection_complement_once() {
    let mut r = request(vec![sample("a", &["GGGG", "GGGG", "GGGG"])]);
    for u in &mut r.samples[0].units {
        for b in &mut u.hypotheses[0].bases {
            b.query_to_reference_strand = '-';
            b.original_base = 'C';
            b.query_position_0based = 3 - b.query_position_0based;
        }
    }
    let report = run(r);
    let s = report
        .sites
        .iter()
        .find(|s| s.position_0based == 101)
        .unwrap();
    let plus = s
        .transcripts
        .iter()
        .find(|t| t.transcript_id == "plus")
        .unwrap();
    let minus = s
        .transcripts
        .iter()
        .find(|t| t.transcript_id == "minus")
        .unwrap();
    let absent = s
        .transcripts
        .iter()
        .find(|t| t.transcript_id == "short")
        .unwrap();
    assert_eq!(
        (
            plus.cdna_position_1based,
            plus.transcript_reference_base,
            plus.transcript_alternate_bases.clone()
        ),
        (Some(2), Some('A'), vec!['G'])
    );
    assert_eq!(
        (
            minus.cdna_position_1based,
            minus.transcript_reference_base,
            minus.transcript_alternate_bases.clone()
        ),
        (Some(19), Some('T'), vec!['C'])
    );
    assert_eq!(absent.status, "not_present");
    assert_eq!(absent.cdna_position_1based, None);
    assert_eq!(
        report.samples[0].phase_blocks[0].candidate_haplotype_count,
        Some(1)
    );
}

#[test]
fn cdna_comparison_three_supported_patterns_are_not_forced_to_two_or_called_copy_number() {
    let r = request(vec![sample(
        "a",
        &[
            "AAAA", "AAAA", "AAAA", "GGGG", "GGGG", "GGGG", "CCCC", "CCCC", "CCCC",
        ],
    )]);
    let report = run(r);
    let block = &report.samples[0].phase_blocks[0];
    assert_eq!(block.candidate_haplotype_count, Some(3));
    assert!(block.excess_expressed_haplotype_diversity);
    assert_eq!(block.copy_number_status, "not_assessed");
    assert_eq!(
        block
            .haplotypes
            .iter()
            .map(|h| h.assigned_units)
            .sum::<usize>(),
        9
    );
    assert!(block.haplotypes.iter().all(|h| h.alleles.len() == 4));
}

#[test]
fn cdna_comparison_quality_weighted_error_assignment_and_partial_ties_are_retained() {
    let mut r = diploid();
    let mut error = unit("error", "AAAT");
    error.hypotheses[0].bases[3].phred = Some(20);
    r.samples[0].units.push(error);
    let report = run(r);
    let a = report.samples[0].phase_blocks[0]
        .assignments
        .iter()
        .find(|a| a.unit_id == "error")
        .unwrap();
    assert_eq!(a.status, "assigned_error_tolerant");
    assert_eq!(a.mismatch_count, Some(1));
    assert_eq!(a.phred_penalty, Some(20));
    let mut r = request(vec![sample(
        "a",
        &["AAGG", "AAGG", "AAGG", "AATT", "AATT", "AATT"],
    )]);
    r.supplied_sites = POS.to_vec();
    let mut partial = unit("partial", "AAAA");
    partial.hypotheses[0].bases.truncate(2);
    r.samples[0].units.push(partial);
    let report = run(r);
    let a = report.samples[0].phase_blocks[0]
        .assignments
        .iter()
        .find(|a| a.unit_id == "partial")
        .unwrap();
    assert_eq!(a.status, "ambiguous");
    assert_eq!(a.compatible_haplotype_ids.len(), 2);
}

#[test]
fn cdna_comparison_disconnected_phase_blocks_are_not_concatenated() {
    let mut r = request(vec![sample(
        "a",
        &[
            "AAAA", "AAAA", "AAAA", "GGGG", "GGGG", "GGGG", "AAAA", "AAAA", "AAAA", "GGGG", "GGGG",
            "GGGG",
        ],
    )]);
    for (i, u) in r.samples[0].units.iter_mut().enumerate() {
        if i < 6 {
            u.hypotheses[0].bases.truncate(2);
        } else {
            u.hypotheses[0].bases.drain(..2);
        }
    }
    let report = run(r);
    assert_eq!(report.samples[0].phase_blocks.len(), 2);
    assert!(
        report.samples[0]
            .phase_blocks
            .iter()
            .all(|b| b.candidate_haplotype_count == Some(2)
                && !b.excess_expressed_haplotype_diversity)
    );
}

#[test]
fn cdna_comparison_connected_partial_vectors_phase_transitively_without_a_full_read() {
    let mut r = request(vec![sample(
        "a",
        &[
            "AAAA", "AAAA", "AAAA", "AAAA", "AAAA", "AAAA", "AAAA", "AAAA", "AAAA",
        ],
    )]);
    r.supplied_sites = POS.to_vec();
    for (i, u) in r.samples[0].units.iter_mut().enumerate() {
        let left = i / 3;
        u.hypotheses[0].bases = u.hypotheses[0].bases[left..left + 2].to_vec();
    }
    let report = run(r);
    let b = &report.samples[0].phase_blocks[0];
    assert_eq!(b.sites_0based.len(), 4);
    assert_eq!(b.candidate_haplotype_count, Some(1));
    assert_eq!(vectors(b), BTreeSet::from(["AAAA".into()]));
    let h = &b.haplotypes[0];
    assert_eq!(h.phase_method, "transitive_supported_local_patterns");
    assert!(h.exact_seed_units.is_empty());
    assert!(!h.support_links.is_empty());
    assert!(h.support_links.iter().all(|link| link.unit_ids.len() >= 3));
}

#[test]
fn cdna_comparison_explicit_phase_markers_preserve_all_snp_evidence_and_scope_templates() {
    let mut r = request(vec![sample(
        "a",
        &[
            "AAAA", "AAAA", "AAAA", "GGGG", "GGGG", "GGGG", "AAAA", "AAAA", "AAAA", "GGGG", "GGGG",
            "GGGG", "AAAA", "AAAA", "AAAA", "GGGG", "GGGG", "GGGG",
        ],
    )]);
    r.supplied_sites = POS.to_vec();
    // Three overlapping two-marker windows phase four sites transitively.
    // An explicit middle-window question changes the phase basis only.
    for (i, u) in r.samples[0].units.iter_mut().enumerate() {
        let left = i / 6;
        u.hypotheses[0].bases = u.hypotheses[0].bases[left..left + 2].to_vec();
    }
    let broad = run(r.clone());
    let transitive = &broad.samples[0].phase_blocks[0];
    assert_eq!(transitive.sites_0based, POS);
    assert_eq!(transitive.candidate_haplotype_count, Some(2));
    assert_eq!(
        vectors(transitive),
        BTreeSet::from(["AAAA".into(), "GGGG".into()])
    );
    assert!(
        transitive
            .haplotypes
            .iter()
            .all(|h| h.exact_seed_units.is_empty()
                && h.phase_method == "transitive_supported_local_patterns")
    );
    assert!(
        serde_json::to_value(&r)
            .unwrap()
            .get("phase_marker_positions")
            .is_none()
    );
    assert!(
        serde_json::to_value(&broad)
            .unwrap()
            .get("phase_marker_positions")
            .is_none()
    );

    r.phase_marker_positions = Some(vec![121, 103]);
    let scoped = run(r.clone());
    assert_eq!(scoped.phase_marker_positions, Some(vec![103, 121]));
    assert_eq!(scoped.samples[0].total_units, 18);
    assert_eq!(scoped.samples[0].total_units, broad.samples[0].total_units);
    assert_eq!(scoped.sites.len(), 4);
    assert_eq!(scoped.samples[0].site_calls.len(), 72);
    assert_eq!(scoped.samples[0].site_observations.len(), 36);
    assert_eq!(
        serde_json::to_value(&scoped.sites).unwrap(),
        serde_json::to_value(&broad.sites).unwrap()
    );
    assert_eq!(
        serde_json::to_value(&scoped.samples[0].site_calls).unwrap(),
        serde_json::to_value(&broad.samples[0].site_calls).unwrap()
    );
    assert_eq!(
        serde_json::to_value(&scoped.samples[0].site_observations).unwrap(),
        serde_json::to_value(&broad.samples[0].site_observations).unwrap()
    );
    assert_eq!(
        serde_json::to_value(&scoped.policy).unwrap(),
        serde_json::to_value(&broad.policy).unwrap()
    );
    assert_ne!(scoped.request_sha256, broad.request_sha256);
    assert_eq!(scoped.samples[0].phase_blocks.len(), 1);
    let block = &scoped.samples[0].phase_blocks[0];
    assert_eq!(block.sites_0based, [103, 121]);
    assert_eq!(block.candidate_haplotype_count, Some(2));
    assert!(block.haplotypes.iter().all(|h| {
        h.exact_seed_units.len() == 3
            && h.assigned_units == 3
            && h.alleles
                .iter()
                .map(|a| a.position_0based)
                .collect::<Vec<_>>()
                == [103, 121]
    }));
    assert_eq!(block.assignments.len(), 18);
    assert_eq!(block.copy_number_status, "not_assessed");
    assert!(!block.excess_expressed_haplotype_diversity);
    r.phase_marker_positions = Some(vec![103, 121]);
    assert_eq!(
        serde_json::to_value(scoped).unwrap(),
        serde_json::to_value(run(r)).unwrap()
    );
}

#[test]
fn cdna_comparison_explicit_phase_markers_require_unique_bounded_supplied_positions() {
    for which in 0..5 {
        let mut r = diploid();
        r.supplied_sites = POS.to_vec();
        r.phase_marker_positions = Some(vec![101, 103]);
        match which {
            0 => r.phase_marker_positions = Some(vec![]),
            1 => r.phase_marker_positions = Some(vec![101, 101]),
            2 => r.phase_marker_positions = Some(vec![99]),
            3 => r.supplied_sites = vec![101],
            _ => r.budgets.max_phase_sites_per_block = 1,
        }
        let error = compare_cdna_alleles(r, &mut |_, _| true).unwrap_err();
        assert!(
            error.contains("cdna_comparison_phase_marker_positions"),
            "case {which}: {error}"
        );
    }
}

#[test]
fn cdna_comparison_transitive_templates_pool_partial_support_without_inventing_seed_reads() {
    let r = partial_request(&[
        (&[0, 1], "AA", 3),
        (&[0, 1], "GG", 3),
        (&[1, 2], "AA", 3),
        (&[1, 2], "GG", 3),
        // Coverage of this long pair is asymmetric; missing GG here is not
        // evidence against the GGG template linked through the middle marker.
        (&[0, 2], "AA", 3),
    ]);
    let report = run(r.clone());
    let block = &report.samples[0].phase_blocks[0];
    assert_eq!(block.sites_0based, POS[..3]);
    assert_eq!(block.candidate_haplotype_count, Some(2));
    assert_eq!(vectors(block), BTreeSet::from(["AAA".into(), "GGG".into()]));
    assert_eq!(report.samples[0].total_units, 15);
    assert_eq!(report.samples[0].site_calls.len(), 45);
    assert_eq!(report.samples[0].site_observations.len(), 30);
    assert!(block.haplotypes.iter().all(|h| {
        h.phase_method == "transitive_supported_local_patterns"
            && h.exact_seed_units.is_empty()
            && !h.support_links.is_empty()
            && h.support_links.iter().all(|link| link.unit_ids.len() >= 3)
    }));
    assert!(report.samples[0].site_observations.iter().all(|b| {
        b.observation.phred == Some(40)
            && b.observation.mapq == Some(60)
            && b.observation.original_base == b.observation.observed_base
    }));
    assert!(report.samples[0].site_calls.iter().any(|c| {
        c.unit_id == "partial_0_0"
            && c.position_0based == POS[2]
            && c.status == CdnaCallStatus::Unassessed
    }));
    let mut reordered = r;
    reordered.samples[0].units.reverse();
    assert_eq!(
        serde_json::to_value(&report).unwrap(),
        serde_json::to_value(run(reordered)).unwrap()
    );
}

#[test]
fn cdna_comparison_transitive_overlap_ambiguity_requires_an_informative_bridge() {
    let ambiguous = run(partial_request(&[
        (&[0, 1], "AA", 3),
        (&[0, 1], "GA", 3),
        (&[1, 2], "AT", 3),
        (&[1, 2], "AC", 3),
    ]));
    assert_eq!(ambiguous.samples[0].phase_blocks.len(), 1);
    let block = &ambiguous.samples[0].phase_blocks[0];
    assert_eq!(block.candidate_haplotype_count, None);
    assert!(block.haplotypes.is_empty());
    assert!(!block.excess_expressed_haplotype_diversity);

    let star = run(partial_request(&[
        (&[0, 1], "AA", 3),
        (&[1, 2], "AT", 3),
        (&[1, 2], "AC", 3),
    ]));
    assert_eq!(
        star.samples[0].phase_blocks[0].candidate_haplotype_count,
        Some(2)
    );
    assert_eq!(
        vectors(&star.samples[0].phase_blocks[0]),
        BTreeSet::from(["AAT".into(), "AAC".into()])
    );

    let bridged = run(partial_request(&[
        (&[0, 1], "AA", 3),
        (&[0, 1], "GA", 3),
        (&[1, 2], "AT", 3),
        (&[1, 2], "AC", 3),
        (&[0, 2], "AT", 3),
        (&[0, 2], "GC", 3),
    ]));
    assert_eq!(
        bridged.samples[0].phase_blocks[0].candidate_haplotype_count,
        Some(2)
    );
    assert_eq!(
        vectors(&bridged.samples[0].phase_blocks[0]),
        BTreeSet::from(["AAT".into(), "GAC".into()])
    );
    assert!(
        bridged.samples[0].phase_blocks[0]
            .haplotypes
            .iter()
            .all(|h| h.exact_seed_units.is_empty())
    );
}

#[test]
fn cdna_comparison_transitive_conflicting_supported_cycles_remain_unresolved() {
    let report = run(partial_request(&[
        (&[0, 1], "AA", 3),
        (&[0, 1], "CC", 3),
        (&[1, 2], "AA", 3),
        (&[1, 2], "CC", 3),
        (&[0, 2], "AC", 3),
        (&[0, 2], "CA", 3),
    ]));
    assert_eq!(report.samples[0].phase_blocks.len(), 1);
    let block = &report.samples[0].phase_blocks[0];
    assert_eq!(block.candidate_haplotype_count, None);
    assert!(block.haplotypes.is_empty());
    assert_eq!(report.samples[0].total_units, 18);
    assert_eq!(report.samples[0].site_observations.len(), 36);
    assert_eq!(report.sites.len(), 3);
}

#[test]
fn cdna_comparison_transitive_partial_reads_allow_three_supported_classes() {
    let report = run(partial_request(&[
        (&[0, 1], "AA", 3),
        (&[0, 1], "CC", 3),
        (&[0, 1], "GG", 3),
        (&[1, 2], "AT", 3),
        (&[1, 2], "CA", 3),
        (&[1, 2], "GC", 3),
    ]));
    let block = &report.samples[0].phase_blocks[0];
    assert_eq!(block.candidate_haplotype_count, Some(3));
    assert_eq!(
        vectors(block),
        BTreeSet::from(["AAT".into(), "CCA".into(), "GGC".into()])
    );
    assert!(block.excess_expressed_haplotype_diversity);
    assert_eq!(block.copy_number_status, "not_assessed");
    assert!(
        block
            .haplotypes
            .iter()
            .all(|h| h.exact_seed_units.is_empty())
    );
}

#[test]
fn cdna_comparison_weak_bridges_and_duplicate_observations_do_not_create_phase_votes() {
    let mut r = partial_request(&[
        (&[0, 1], "AA", 3),
        (&[0, 1], "GG", 3),
        (&[2, 3], "AA", 3),
        (&[2, 3], "GG", 3),
        (&[1, 2], "AG", 1),
    ]);
    let weak = r.samples[0].units.last_mut().unwrap();
    let duplicate_observations = weak.hypotheses[0].bases.clone();
    for _ in 0..4 {
        weak.hypotheses[0]
            .bases
            .extend(duplicate_observations.clone());
    }
    let report = run(r);
    assert_eq!(report.samples[0].total_units, 13);
    assert_eq!(report.samples[0].phase_blocks.len(), 2);
    assert_eq!(report.samples[0].phase_blocks[0].sites_0based, POS[..2]);
    assert_eq!(report.samples[0].phase_blocks[1].sites_0based, POS[2..]);
    assert!(report.samples[0].phase_blocks.iter().all(|b| {
        b.candidate_haplotype_count == Some(2)
            && b.haplotypes.iter().all(|h| {
                h.support_links
                    .iter()
                    .all(|link| !link.unit_ids.contains(&"partial_4_0".into()))
            })
    }));
}

#[test]
fn cdna_comparison_transitive_templates_never_stitch_other_reads_splices_into_a_switch() {
    let mut r = partial_request(&[
        (&[0, 1], "AA", 4),
        (&[0, 1], "GG", 4),
        (&[1, 2, 3], "AAA", 4),
        (&[1, 2, 3], "GGG", 4),
    ]);
    for u in &mut r.samples[0].units {
        if u.unit_id.starts_with("partial_2_") {
            u.hypotheses[0].junctions = vec![(110, 120)];
        } else if u.unit_id.starts_with("partial_3_") {
            u.hypotheses[0].junctions = vec![(104, 120)];
        }
    }
    let mut inspected = unit("switch", "AAGG");
    inspected.hypotheses[0].junctions.clear();
    r.samples[0].units.push(inspected);
    let report = run(r);
    let block = &report.samples[0].phase_blocks[0];
    assert_eq!(
        vectors(block),
        BTreeSet::from(["AAAA".into(), "GGGG".into()])
    );
    assert!(
        block
            .haplotypes
            .iter()
            .all(|h| h.exact_seed_units.is_empty())
    );
    let candidate = report.samples[0]
        .switch_candidates
        .iter()
        .find(|c| c.unit_id == "switch")
        .unwrap();
    assert_eq!(candidate.status, "haplotype_switch_candidate");
    assert!(
        candidate
            .explanations
            .iter()
            .all(|e| e.same_hypothesis_junctions.is_empty())
    );
    let rows = &report.samples[0].haplotype_junction_support;
    assert!(rows.iter().any(|r| r.junctions == [(110, 120)]));
    assert!(rows.iter().any(|r| r.junctions == [(104, 120)]));
    assert!(rows.iter().all(|r| r.junctions.len() <= 1));
}

#[test]
fn cdna_comparison_transitive_switch_leave_one_out_checks_every_support_link() {
    let mut r = partial_request(&[
        (&[0, 1], "AA", 2),
        (&[0, 1], "GG", 4),
        (&[1, 2, 3], "AAA", 4),
        (&[1, 2, 3], "GGG", 4),
    ]);
    r.samples[0].units.push(unit("switch", "AAGG"));
    let report = run(r.clone());
    let block = &report.samples[0].phase_blocks[0];
    let template = block
        .haplotypes
        .iter()
        .find(|h| {
            h.alleles
                .iter()
                .map(|a| a.observed_base)
                .collect::<String>()
                == "AAAA"
        })
        .unwrap();
    assert!(template.exact_seed_units.is_empty());
    assert!(
        template
            .support_links
            .iter()
            .any(|link| link.unit_ids.len() == 3 && link.unit_ids.contains(&"switch".into()))
    );
    assert!(
        template
            .support_links
            .iter()
            .any(|link| link.unit_ids.len() >= 4)
    );
    assert!(
        !report.samples[0]
            .switch_candidates
            .iter()
            .any(|c| c.unit_id == "switch")
    );

    let mut independent = unit("independent_AB_AA", "AAAA");
    independent.hypotheses[0].bases.truncate(2);
    independent.hypotheses[0].junctions.clear();
    r.samples[0].units.push(independent);
    let supported = run(r);
    let candidate = supported.samples[0]
        .switch_candidates
        .iter()
        .find(|c| c.unit_id == "switch")
        .unwrap();
    assert_eq!(candidate.status, "junction_associated_switch_candidate");
    assert!(
        candidate
            .explanations
            .iter()
            .any(|e| e.left_informative_markers >= 2
                && e.right_informative_markers >= 2
                && e.same_hypothesis_junctions.get("mapping_switch") == Some(&vec![(110, 120)]))
    );
}

#[test]
fn cdna_comparison_switch_is_reported_with_same_hypothesis_junction_not_as_mechanism() {
    let mut r = diploid();
    r.samples[0].units.push(unit("switch", "AAGG"));
    let report = run(r);
    let candidate = report.samples[0]
        .switch_candidates
        .iter()
        .find(|c| c.unit_id == "switch")
        .unwrap();
    assert_eq!(candidate.status, "junction_associated_switch_candidate");
    assert_eq!(candidate.mechanism_status, "not_established");
    assert_eq!(candidate.explanations.len(), 1);
    let e = &candidate.explanations[0];
    assert_eq!(
        (e.left_informative_markers, e.right_informative_markers),
        (2, 2)
    );
    assert_eq!(e.marker_bracket_0based, (103, 122));
    assert_eq!(e.same_hypothesis_junctions["mapping_switch"], [(110, 120)]);
    assert_eq!(e.query_marker_positions["mapping_switch"], [(1, 2)]);
    assert_eq!(candidate.whole_template_mismatch_count, Some(2));
    let mut r = diploid();
    let mut u = unit("switch", "AAGG");
    u.hypotheses[0].junctions.clear();
    let mut nonoptimal = u.hypotheses[0].clone();
    nonoptimal.hypothesis_id = "false_junction".into();
    nonoptimal.eligible = false;
    nonoptimal.junctions = vec![(110, 120)];
    u.hypotheses.push(nonoptimal);
    r.samples[0].units.push(u);
    let report = run(r);
    let candidate = report.samples[0]
        .switch_candidates
        .iter()
        .find(|c| c.unit_id == "switch")
        .unwrap();
    assert_eq!(candidate.status, "haplotype_switch_candidate");
    assert!(
        candidate
            .explanations
            .iter()
            .all(|e| e.same_hypothesis_junctions.is_empty())
    );
}

#[test]
fn cdna_comparison_switch_scores_disagreement_at_shared_template_markers() {
    let mut r = diploid();
    r.samples[0].units.push(unit("switch", "AAGG"));
    for unit in &mut r.samples[0].units {
        let base = if unit.unit_id == "switch" { 'C' } else { 'A' };
        unit.hypotheses[0].bases.push(CdnaBaseObservation {
            reference_position_0based: 125,
            query_position_0based: 4,
            observed_base: base,
            original_base: base,
            query_to_reference_strand: '+',
            phred: Some(40),
            mapq: Some(60),
        });
    }
    r.supplied_sites.push(125);
    r.policy.maximum_assignment_mismatches = 0;
    let exact_only = run(r.clone());
    assert!(exact_only.samples[0].switch_candidates.is_empty());

    r.policy.maximum_assignment_mismatches = 1;
    let tolerant = run(r);
    let switch = tolerant.samples[0]
        .switch_candidates
        .iter()
        .find(|s| s.unit_id == "switch")
        .unwrap();
    assert_eq!(switch.explanations.len(), 1);
    assert_eq!(switch.explanations[0].mismatch_count, 1);
    assert_eq!(switch.explanations[0].phred_penalty, 40);
}

#[test]
fn cdna_comparison_stable_third_haplotype_is_a_competing_switch_explanation() {
    let r = request(vec![sample(
        "a",
        &[
            "AAAA", "AAAA", "AAAA", "GGGG", "GGGG", "GGGG", "AAGG", "AAGG", "AAGG", "AAGG",
        ],
    )]);
    let report = run(r);
    let c = report.samples[0]
        .switch_candidates
        .iter()
        .find(|c| c.unit_id == "unit_6")
        .unwrap();
    assert_eq!(c.status, "whole_haplotype_alternative");
    assert!(c.competing_whole_template);
    assert_eq!(
        report.samples[0].phase_blocks[0].candidate_haplotype_count,
        Some(3)
    );
}

#[test]
fn cdna_comparison_sample_differences_use_alleles_not_local_haplotype_labels() {
    let report = run(request(vec![
        sample("a", &["AAAA", "AAAA", "AAAA"]),
        sample("b", &["GGGG", "GGGG", "GGGG"]),
    ]));
    assert_eq!(report.sample_differences.len(), 4);
    assert!(
        report
            .sample_differences
            .iter()
            .all(|d| d.status == "different_observed_candidate_allele_sets")
    );
    let a = &report.samples[0].phase_blocks[0].haplotypes[0];
    let b = &report.samples[1].phase_blocks[0].haplotypes[0];
    assert_eq!(a.haplotype_id, "H1");
    assert_eq!(b.haplotype_id, "H1");
    assert_ne!(a.fingerprint, b.fingerprint);
    let report = run(request(vec![
        sample("a", &["GGGG", "GGGG", "GGGG"]),
        sample("missing", &[]),
    ]));
    assert!(
        report
            .sample_differences
            .iter()
            .all(|d| d.status == "unassessed_comparison")
    );
}

#[test]
fn cdna_comparison_haplotype_junction_sets_count_unit_once_and_keep_alternatives() {
    let mut r = diploid();
    let mut h = r.samples[0].units[0].hypotheses[0].clone();
    h.hypothesis_id = "other_chain".into();
    h.junctions.clear();
    r.samples[0].units[0].hypotheses.push(h);
    let report = run(r);
    let rows = &report.samples[0].haplotype_junction_support;
    assert!(
        rows.iter()
            .any(|r| r.junctions.is_empty() && r.possible_units == 1 && r.exclusive_units == 0)
    );
    assert!(
        rows.iter().any(|r| r.junctions == [(110, 120)]
            && r.possible_units == 3
            && r.exclusive_units == 2)
    );
}

#[test]
fn cdna_comparison_input_permutations_have_stable_reports_and_hashes() {
    let r = diploid();
    let expected = run(r.clone());
    let mut reordered = r;
    reordered.transcripts.reverse();
    reordered.samples[0].units.reverse();
    for u in &mut reordered.samples[0].units {
        u.hypotheses.reverse();
        u.hypotheses[0].bases.reverse();
    }
    let actual = run(reordered);
    assert_eq!(
        serde_json::to_value(expected).unwrap(),
        serde_json::to_value(actual).unwrap()
    );
}

#[test]
fn cdna_comparison_shared_normalization_matches_the_saved_evidence_hash() {
    let mut request = diploid();
    request.samples.push(sample("b", &["AAAA", "GGGG", "GGGG"]));
    request.supplied_sites = POS.to_vec();
    request.phase_marker_positions = Some(POS.to_vec());
    for sample in &mut request.samples {
        for unit in &mut sample.units {
            let mut alternative = unit.hypotheses[0].clone();
            alternative.hypothesis_id = "audit_alternative".into();
            alternative.eligible = false;
            unit.hypotheses.push(alternative);
        }
    }
    let expected = run(request.clone());
    request.transcripts.reverse();
    request.samples.reverse();
    request.supplied_sites.reverse();
    request.phase_marker_positions.as_mut().unwrap().reverse();
    for sample in &mut request.samples {
        sample.units.reverse();
        for unit in &mut sample.units {
            unit.hypotheses.reverse();
            for hypothesis in &mut unit.hypotheses {
                hypothesis.bases.reverse();
            }
        }
    }
    assert_ne!(hash(&request).unwrap(), expected.request_sha256);
    assert_eq!(
        serde_json::to_value(run(request.clone())).unwrap(),
        serde_json::to_value(&expected).unwrap()
    );
    normalize_cdna_comparison_request(&mut request);
    assert_eq!(hash(&request).unwrap(), expected.request_sha256);
    let once = serde_json::to_vec(&request).unwrap();
    normalize_cdna_comparison_request(&mut request);
    assert_eq!(serde_json::to_vec(&request).unwrap(), once);
}

#[test]
fn cdna_comparison_schema_hash_geometry_and_quality_refuse_invalid_input() {
    for which in 0..8 {
        let mut r = diploid();
        match which {
            0 => r.schema = "wrong".into(),
            1 => r.reference.sha256 = sha(b"wrong"),
            2 => r.transcripts[0].strand = '?',
            3 => r.samples[0].units[0].hypotheses[0].bases[0].phred = Some(94),
            4 => r.samples[0].units[0].hypotheses[0].bases[0].reference_position_0based = 115,
            5 => r.samples[0].units[0].hypotheses[0].bases[0].observed_base = 'G',
            6 => {
                let duplicate = r.samples[0].units[0].clone();
                r.samples[0].units.push(duplicate);
            }
            _ => r.policy.minimum_alternate_fraction = f64::NAN,
        }
        assert!(
            compare_cdna_alleles(r, &mut |_, _| true).is_err(),
            "case {which}"
        );
    }
    let mut json = serde_json::to_value(diploid()).unwrap();
    json["force_diploid"] = serde_json::json!(true);
    assert!(serde_json::from_value::<CdnaComparisonRequest>(json).is_err());
}

#[test]
fn cdna_comparison_source_hash_case_cannot_disguise_reused_units() {
    let mut r = diploid();
    let mut duplicate = r.samples[0].clone();
    duplicate.sample_id = "b".into();
    duplicate.label = "invented_b".into();
    r.samples.push(duplicate);
    assert!(compare_cdna_alleles(r.clone(), &mut |_, _| true).is_err());
    r.samples[1].source_sha256 = format!(
        "sha256:{}",
        r.samples[1]
            .source_sha256
            .strip_prefix("sha256:")
            .unwrap()
            .to_ascii_uppercase()
    );
    assert!(compare_cdna_alleles(r, &mut |_, _| true).is_err());
}

#[test]
fn cdna_comparison_budgets_and_cancellation_refuse_instead_of_partial_success() {
    for which in 0..7 {
        let mut r = diploid();
        match which {
            0 => r.budgets.max_units = 1,
            1 => r.budgets.max_sites = 1,
            2 => r.budgets.max_base_observations = 1,
            3 => r.budgets.max_comparisons = 1,
            4 => r.budgets.max_phase_sites_per_block = 2,
            5 => r.policy.maximum_haplotypes = 1,
            _ => r.budgets.max_output_bytes = 1,
        }
        assert!(
            compare_cdna_alleles(r, &mut |_, _| true).is_err(),
            "case {which}"
        );
    }
    assert_eq!(
        compare_cdna_alleles(diploid(), &mut |_, _| false).unwrap_err(),
        "cdna_comparison_cancelled"
    );
    let mut count = 0;
    assert!(
        compare_cdna_alleles(diploid(), &mut |_, _| {
            count += 1;
            count < 5
        })
        .is_err()
    );
    assert_eq!(count, 5);
}
