//! Inline synthetic SAM fields/loci, independently enumerated base mappings and PAF equivalents.
//!
//! No reads, primers or external records. Recreate with the library's
//! `observed_transcript_structures::sam::tests` filter; no fixture files or I/O.

use super::*;
use crate::engine::{PRIMER_PAYLOAD_ANNOTATION_SCHEMA, PrimerPayloadTranscript};

fn sha(c: char) -> String {
    format!("sha256:{}", c.to_string().repeat(64))
}

fn row(ops: Vec<(char, usize)>) -> Row {
    Row {
        ordinal: 1,
        qname: "alias_not_identity".into(),
        flag: 0,
        target: Some("ctg_a".into()),
        position_1based: 101,
        mapq: 0,
        cigar: Some(ops),
        sequence_length: None,
        quality_length: None,
        has_cg: false,
    }
}

fn input() -> SamInput {
    let transcripts = [
        (
            "gene_a",
            "ctg_a",
            '+',
            vec![(100, 110), (120, 130), (140, 150)],
        ),
        ("gene_b", "ctg_b", '-', vec![(200, 210), (225, 235)]),
    ]
    .into_iter()
    .map(|(gene, contig, strand, exons)| PrimerPayloadTranscript {
        transcript_id: format!("{gene}.1"),
        gene_ids: vec![gene.into()],
        target_name: contig.into(),
        target_length: 1000,
        genomic_contig: Some(contig.into()),
        genomic_strand: Some(strand),
        genomic_exons: exons,
    })
    .collect();
    SamInput {
        source_sha256: sha('4'),
        reference: ReferenceBinding {
            assembly: "synthetic".into(),
            fasta_sha256: sha('1'),
        },
        reference_gene: ReferenceGene {
            source_id: "source".into(),
            gene_id: "gene_a".into(),
        },
        comparison_genes: vec![],
        annotation_sources: vec![AnnotationSource {
            source_id: "source".into(),
            provider: Some("synthetic".into()),
            annotation_sha256: sha('2'),
            annotation: PrimerPayloadAnnotation {
                schema: PRIMER_PAYLOAD_ANNOTATION_SCHEMA.into(),
                reference_id: "declared".into(),
                assembly: "synthetic".into(),
                release: "test_snapshot".into(),
                reference_fasta_sha256: sha('1'),
                coordinate_convention: "zero_based_half_open".into(),
                transcripts,
            },
        }],
        dictionary: BTreeMap::from([("ctg_a".into(), 1000), ("ctg_b".into(), 1000)]),
        unpadded_reference: true,
        units: vec![Unit {
            unit_id: "unit_a".into(),
            declared_axis_length: Some(30),
            masks: Some(vec![]),
            ends: [vec![EndConstraint::Unknown], vec![EndConstraint::Unknown]],
            search_declared_complete: true,
            views_declared_complete: true,
            rows: vec![row(vec![
                ('M', 10),
                ('N', 10),
                ('M', 10),
                ('N', 10),
                ('M', 10),
            ])],
        }],
        required_flank_bases: 2,
        budgets: Budgets {
            max_queries: 100,
            max_records: 100,
            max_segments_per_query: 100,
            max_operations: 10000,
            max_skips_per_chain: 100,
            max_chain_families: 100,
            max_alternatives_per_record: 100,
            max_annotation_rows: 100,
            max_model_comparisons_per_family: 100,
            max_total_work: 1_000_000,
        },
        max_rows: 100,
    }
}

fn run(a: &SamInput) -> Analysis {
    analyze(a, &mut || Ok(())).unwrap()
}

/// Walk individual bases including unavailable clips, then reverse query labels.
/// This does not call the normalizer or use its interval/clip arithmetic.
fn basewise(ops: &[(char, usize)], reverse: bool) -> (Vec<SupportedBlock>, Vec<Interval>) {
    let (mut query, mut target) = (0, 100);
    let (mut pairs, mut skips) = (Vec::new(), Vec::new());
    for &(kind, length) in ops {
        let start = target;
        for _ in 0..length {
            match kind {
                'M' | '=' | 'X' => {
                    pairs.push((query, target));
                    query += 1;
                    target += 1;
                }
                'H' | 'S' | 'I' => query += 1,
                'D' | 'N' => target += 1,
                'P' => (),
                _ => panic!("oracle input"),
            }
        }
        if kind == 'N' {
            skips.push((start, target));
        }
    }
    if reverse {
        for (q, _) in &mut pairs {
            *q = query - 1 - *q;
        }
    }
    pairs.sort();
    let mut blocks: Vec<SupportedBlock> = Vec::new();
    for (q, t) in pairs {
        if let Some(b) = blocks.last_mut()
            && b.query.1 == q
            && ((!reverse && b.target.1 == t) || (reverse && t + 1 == b.target.0))
        {
            b.query.1 += 1;
            b.target = (b.target.0.min(t), b.target.1.max(t + 1));
        } else {
            blocks.push(SupportedBlock {
                query: (q, q + 1),
                target: (t, t + 1),
            });
        }
    }
    (blocks, skips)
}

#[test]
fn sam_geometry_matches_independent_basewise_oracle_both_orientations() {
    let cases = [
        vec![
            ('H', 2),
            ('S', 3),
            ('M', 4),
            ('I', 2),
            ('M', 5),
            ('S', 7),
            ('H', 1),
        ],
        vec![
            ('H', 1),
            ('S', 2),
            ('M', 3),
            ('I', 2),
            ('=', 4),
            ('D', 3),
            ('X', 2),
            ('N', 5),
            ('M', 4),
            ('P', 2),
            ('M', 1),
            ('S', 3),
            ('H', 2),
        ],
        vec![('M', 3), ('D', 2), ('M', 4), ('N', 7), ('=', 2)],
    ];
    for ops in cases {
        for reverse in [false, true] {
            let mut r = row(ops.clone());
            r.flag = if reverse { 0x10 } else { 0 };
            let n = normalize(&r).unwrap().unwrap();
            let g = geometry::walk(
                n.operations.iter().map(|&(c, n)| (c as u8, n)),
                n.query,
                n.target,
                if reverse { '-' } else { '+' },
            )
            .unwrap();
            let (blocks, skips) = basewise(&ops, reverse);
            assert_eq!(g.blocks, blocks, "{ops:?} reverse={reverse}");
            assert_eq!(g.skips, skips);
        }
    }
    let mut r = row(vec![
        ('H', 2),
        ('S', 3),
        ('M', 4),
        ('I', 2),
        ('M', 5),
        ('S', 7),
        ('H', 1),
    ]);
    r.sequence_length = Some(21);
    r.quality_length = Some(21);
    let n = normalize(&r).unwrap().unwrap();
    assert_eq!(
        (n.axis_length, n.query, n.clipping.hard, n.clipping.soft),
        (24, (5, 16), [2, 1], [3, 7])
    );
    r.flag = 0x10;
    assert_eq!(normalize(&r).unwrap().unwrap().query, (8, 19));
}

#[test]
fn sam_clip_cartesian_oracle_and_padding_leave_operations_unchanged() {
    for h in 0..3 {
        for s in 0..3 {
            for i in 0..3 {
                for reverse in [false, true] {
                    let mut ops = vec![];
                    if h > 0 {
                        ops.push(('H', h));
                    }
                    if s > 0 {
                        ops.push(('S', s));
                    }
                    ops.push(('M', 2));
                    if i > 0 {
                        ops.push(('I', i));
                    }
                    ops.extend([
                        ('D', 1),
                        ('=', 2),
                        ('P', 1),
                        ('N', 3),
                        ('X', 2),
                        ('S', 2),
                        ('H', 1),
                    ]);
                    let mut r = row(ops.clone());
                    if reverse {
                        r.flag = 0x10;
                    }
                    let n = normalize(&r).unwrap().unwrap();
                    assert_eq!(n.clipping.padding_operations, 1);
                    assert_eq!(
                        n.operations,
                        ops.iter()
                            .filter(|(c, _)| !"HSP".contains(*c))
                            .copied()
                            .collect::<Vec<_>>()
                    );
                    let g = geometry::walk(
                        n.operations.iter().map(|&(c, n)| (c as u8, n)),
                        n.query,
                        n.target,
                        if reverse { '-' } else { '+' },
                    )
                    .unwrap();
                    let (blocks, skips) = basewise(&ops, reverse);
                    assert_eq!((g.blocks, g.skips), (blocks, skips));
                }
            }
        }
    }
}

#[test]
fn sam_and_paf_share_only_physical_chain_identity() {
    let a = input();
    let sam = run(&a);
    let r = &a.units[0].rows[0];
    let n = normalize(r).unwrap().unwrap();
    let g = geometry::walk(
        n.operations.iter().map(|&(c, n)| (c as u8, n)),
        n.query,
        n.target,
        '+',
    )
    .unwrap();
    let paf = Input {
        reference: a.reference.clone(),
        reference_gene: a.reference_gene.clone(),
        comparison_genes: vec![],
        annotation_sources: a.annotation_sources.clone(),
        required_flank_bases: 2,
        budgets: a.budgets.clone(),
        unassessed_records: vec![],
        queries: vec![QueryEvidence {
            record: ReadRecord {
                source_sha256: a.source_sha256.clone(),
                record_index: 0,
                read_id: r.qname.clone(),
            },
            alignment_sha256: a.source_sha256.clone(),
            query_id: "unit_a".into(),
            query_length: 30,
            scope: PrimerPayloadReferenceScope::Genome,
            target_lengths: a.dictionary.clone(),
            masks: Some(vec![]),
            ends: a.units[0].ends.clone(),
            search_complete: true,
            segments: vec![PrimerPayloadSegment {
                query_interval: n.query,
                target_name: r.target.clone().unwrap(),
                target_interval: n.target,
                strand: '+',
                mapq: Some(0),
                reported_role: None,
                reported_identity_per_mille: 1000,
                edits: None,
                alignment_detail_verified: true,
                alignment_operations: n.operations,
                aligned_query_intervals: union(g.blocks.iter().map(|b| b.query).collect()),
                gene_ids: vec![],
                transcript_ids: vec![],
                splice_edges: vec![],
            }],
        }],
    };
    let legacy = summarize(&paf).unwrap();
    assert_eq!(
        legacy.families[0].chain_id,
        sam.summary.families[0].chain_id
    );
    assert_ne!(legacy.context_sha256, sam.summary.context_sha256);
    assert_ne!(
        legacy.observations[0].observation_id,
        sam.summary.observations[0].observation_id
    );
    assert_eq!(sam.coordinate_authority, "producer_declared");
    assert_eq!(sam.identity_scope, "sam_declared_query_unit");
    let encoded = serde_json::to_string(&sam.summary).unwrap();
    assert!(!encoded.contains("payload_query_length"));
    assert!(!encoded.contains("reported_identity"));
    assert!(!encoded.contains("record_index"));
    assert_eq!(sam.summary.observations[0].reported_mapq, Some(0));
    let mut b = a.clone();
    b.units[0].rows[0].mapq = 255;
    assert_eq!(run(&b).summary.observations[0].reported_mapq, None);
    b.units[0].rows[0]
        .cigar
        .as_mut()
        .unwrap()
        .insert(0, ('H', 3));
    b.units[0].declared_axis_length = Some(33);
    let clipped = run(&b);
    assert_eq!(
        clipped.summary.families[0].chain_id,
        sam.summary.families[0].chain_id
    );
    assert_ne!(
        clipped.summary.observations[0].observation_id,
        sam.summary.observations[0].observation_id
    );
}

#[test]
fn sam_gene_selection_is_configurable_not_a_mapping_filter_or_strand_call() {
    let mut a = input();
    let original = run(&a);
    a.reference_gene.gene_id = "gene_b".into();
    let other = run(&a);
    assert_eq!(
        original.summary.families[0].chain_id,
        other.summary.families[0].chain_id
    );
    assert_ne!(
        original.summary.context_sha256,
        other.summary.context_sha256
    );
    assert_eq!(
        other.summary.observations[0].comparisons[0].models[0].relation,
        ModelRelation::OutsideModel
    );
    a.units[0].declared_axis_length = Some(20);
    let r = &mut a.units[0].rows[0];
    r.target = Some("ctg_b".into());
    r.position_1based = 201;
    r.flag = 0x10;
    r.cigar = Some(vec![('M', 10), ('N', 15), ('M', 10)]);
    let b = run(&a);
    let o = &b.summary.observations[0];
    assert_eq!(
        o.comparisons[0].models[0].relation,
        ModelRelation::ExactFullChain
    );
    assert_eq!(o.comparisons[0].models[0].model_strand, '-');
    assert_eq!(o.skips, [(210, 225)]);
    assert_eq!(o.query_to_reference_strand, '-');
}

#[test]
fn sam_qname_aliases_and_row_roles_do_not_multiply_unit_support() {
    let mut a = input();
    for (i, flag) in [0x100, 0x800, 0x200, 0x400].into_iter().enumerate() {
        let mut r = a.units[0].rows[0].clone();
        r.ordinal = i as u64 + 2;
        r.flag = flag;
        a.units[0].rows.push(r);
    }
    let result = run(&a);
    assert_eq!(result.summary.observations.len(), 5);
    assert_eq!(result.summary.families[0].supporting_keys, 1);
    assert_eq!(result.summary.families[0].definite_keys, 1);
    let mut u = a.units[0].clone();
    u.unit_id = "unit_b".into();
    for r in &mut u.rows {
        r.ordinal += 10;
    }
    a.units.push(u);
    let result = run(&a);
    assert_eq!(result.summary.supports.len(), 2);
    assert_eq!(result.summary.families[0].supporting_keys, 2);
    assert_ne!(
        result.summary.supports[0].key,
        result.summary.supports[1].key
    );
    a.units.reverse();
    a.units.iter_mut().for_each(|u| u.rows.reverse());
    assert_eq!(
        serde_json::to_value(run(&a).summary.observations).unwrap(),
        serde_json::to_value(result.summary.observations).unwrap()
    );
}

#[test]
fn sam_unsupported_competitor_clears_entire_unit_not_just_bad_row() {
    for (case, expected) in [
        Unavailable::Paired,
        Unavailable::Unmapped,
        Unavailable::PaddedReference,
        Unavailable::CigarNotProvided,
        Unavailable::NoSupportedBases,
        Unavailable::CgOrPlaceholder,
        Unavailable::QueryAxisUnsupported,
    ]
    .into_iter()
    .enumerate()
    {
        let mut a = input();
        let mut r = a.units[0].rows[0].clone();
        r.ordinal = 2;
        match case {
            0 => r.flag |= 1,
            1 => {
                r.flag = 4;
                r.cigar = None;
                r.target = None;
                r.position_1based = 0;
            }
            2 => a.unpadded_reference = false,
            3 => r.cigar = None,
            4 => {
                r.cigar = Some(vec![('I', 30)]);
            }
            5 => r.has_cg = true,
            _ => {
                // A hard-clipped valid row over the axis cap; keep both rows' axis consistent.
                for good in &mut a.units[0].rows {
                    good.cigar.as_mut().unwrap().insert(0, ('H', 2_000_000));
                }
                r.cigar.as_mut().unwrap().insert(0, ('H', 2_000_000));
                a.units[0].declared_axis_length = Some(2_000_030);
            }
        }
        a.units[0].rows.push(r);
        let result = run(&a);
        assert!(result.units[0].unavailable.contains(&expected), "{case}");
        assert!(result.summary.observations.is_empty(), "{case}");
        assert!(result.summary.families.is_empty(), "{case}");
        assert!(result.summary.supports[0].possible_chains.is_empty());
        assert_eq!(
            result.summary.supports[0].reconciliation,
            Reconciliation::Unassessed
        );
    }
}

#[test]
fn sam_zero_row_and_unmapped_units_keep_unknown_axes() {
    let mut a = input();
    a.units[0].rows.clear();
    a.units[0].declared_axis_length = None;
    a.units[0].masks = None;
    let result = run(&a);
    assert_eq!(
        result.units[0].unavailable,
        BTreeSet::from([Unavailable::NotProvided])
    );
    assert_eq!(result.units[0].axis_length, None);
    let mut r = row(vec![]);
    r.flag = 4;
    r.cigar = None;
    r.target = None;
    r.position_1based = 0;
    a.units[0].rows.push(r);
    let result = run(&a);
    assert_eq!(result.units[0].rows[0].axis_length, None);
    assert!(result.units[0].unavailable.contains(&Unavailable::Unmapped));
    assert_eq!(result.summary.supports.len(), 1);
}

#[test]
fn sam_masks_and_search_completeness_are_independent_axes() {
    let mut a = input();
    a.units[0].masks = None;
    let missing = run(&a);
    assert_eq!(
        missing.summary.supports[0].reconciliation,
        Reconciliation::Agreement
    );
    assert_eq!(missing.summary.families[0].definite_keys, 1);
    assert_eq!(
        missing.summary.observations[0].flank_qualification,
        vec![Qualification::Unassessed; 2]
    );
    a.units[0].masks = Some(vec![]);
    a.units[0].search_declared_complete = false;
    let incomplete = run(&a);
    assert_eq!(incomplete.summary.families[0].definite_keys, 0);
    assert_eq!(incomplete.summary.families[0].supporting_keys, 1);
    assert_eq!(
        incomplete.summary.observations[0].flank_qualification,
        vec![Qualification::Qualified; 2]
    );
    a.units[0].search_declared_complete = true;
    a.units[0].views_declared_complete = false;
    assert_eq!(
        run(&a).summary.supports[0].reconciliation,
        Reconciliation::Unassessed
    );
    a.units[0].views_declared_complete = true;
    a.units[0].masks = Some(vec![(8, 12)]);
    assert_eq!(
        run(&a).summary.observations[0].flank_qualification,
        [Qualification::InsufficientSupport, Qualification::Qualified]
    );
}

#[test]
fn sam_supplementary_rows_do_not_create_a_junction_between_segments() {
    let mut a = input();
    a.units[0].declared_axis_length = Some(20);
    a.units[0].rows[0].cigar = Some(vec![('M', 10), ('S', 10)]);
    let mut r = a.units[0].rows[0].clone();
    r.ordinal = 2;
    r.flag = 0x800;
    r.position_1based = 401;
    r.cigar = Some(vec![('S', 10), ('M', 10)]);
    a.units[0].rows.push(r);
    let result = run(&a);
    assert_eq!(result.summary.observations.len(), 2);
    assert!(
        result
            .summary
            .observations
            .iter()
            .all(|o| o.skips.is_empty())
    );
    assert!(result.summary.families.is_empty());
    assert_eq!(result.summary.supports[0].topologies, [Topology::Segmented]);
}

#[test]
fn sam_longer_mapping_retains_nonoptimal_observation_without_family_support() {
    let mut a = input();
    a.units[0].declared_axis_length = Some(40);
    a.units[0].rows[0].cigar.as_mut().unwrap().push(('S', 10));
    let mut r = row(vec![('M', 40)]);
    r.ordinal = 2;
    r.target = Some("ctg_b".into());
    r.position_1based = 501;
    a.units[0].rows.push(r);
    let result = run(&a);
    assert_eq!(result.summary.observations.len(), 2);
    let spliced = result
        .summary
        .observations
        .iter()
        .find(|o| !o.skips.is_empty())
        .unwrap();
    assert!(spliced.observed_junction_chain.is_some());
    assert!(
        result.summary.supports[0]
            .hypotheses
            .iter()
            .flatten()
            .all(|id| id != &spliced.observation_id)
    );
    assert!(result.summary.families.is_empty());
}

#[test]
fn sam_invalid_clip_span_sequence_and_binding_are_errors_not_partial_results() {
    for case in 0..17 {
        let mut a = input();
        match case {
            0 => a.units[0].rows[0].cigar = Some(vec![('M', 10), ('S', 1), ('M', 19)]),
            1 => a.units[0].rows[0].cigar = Some(vec![('S', 1), ('H', 1), ('M', 28)]),
            2 => a.units[0].rows[0].cigar = Some(vec![('M', 0)]),
            3 => a.units[0].rows[0].cigar = Some(vec![('B', 30)]),
            4 => a.units[0].rows[0].cigar = Some(vec![('M', usize::MAX), ('I', 1)]),
            5 => a.units[0].rows[0].position_1based = 0,
            6 => a.units[0].rows[0].position_1based = 999,
            7 => a.units[0].rows[0].target = Some("ctg_alias".into()),
            8 => a.units[0].rows[0].sequence_length = Some(29),
            9 => a.units[0].rows[0].quality_length = Some(30),
            10 => a.units[0].declared_axis_length = Some(31),
            11 => {
                let r = a.units[0].rows[0].clone();
                a.units[0].rows.push(r);
            }
            12 => {
                let u = a.units[0].clone();
                a.units.push(u);
            }
            13 => a.units[0].masks = Some(vec![(0, 31)]),
            14 => {
                a.dictionary.insert("ctg_a".into(), 999);
            }
            15 => a.units[0].rows[0].flag = 0x1000,
            _ => {
                let mut r = a.units[0].rows[0].clone();
                r.ordinal = 2;
                r.qname = "contradiction".into();
                a.units[0].rows.push(r);
            }
        }
        assert!(analyze(&a, &mut || Ok(())).is_err(), "invalid case {case}");
    }
}

#[test]
fn sam_placeholder_without_cg_is_unavailable_and_has_no_splice() {
    let mut a = input();
    a.units[0].rows[0].cigar = Some(vec![('S', 30), ('N', 50)]);
    let result = run(&a);
    assert!(
        result.units[0]
            .unavailable
            .contains(&Unavailable::CgOrPlaceholder)
    );
    assert!(result.summary.observations.is_empty());
}

#[test]
fn sam_admission_bounds_exact_one_over_and_cooperative_control() {
    let mut a = input();
    a.max_rows = 1;
    a.budgets.max_queries = 1;
    a.budgets.max_records = 1;
    a.budgets.max_segments_per_query = 1;
    a.budgets.max_operations = 5;
    a.budgets.max_annotation_rows = 2;
    let exact_work = run(&a).summary.comparison_work;
    a.budgets.max_total_work = exact_work;
    assert_eq!(run(&a).summary.comparison_work, exact_work);
    a.budgets.max_total_work -= 1;
    assert_eq!(
        analyze(&a, &mut || Ok(())).unwrap_err().message,
        "primer_payload_chain_comparison_budget"
    );
    for case in 0..5 {
        let mut a = input();
        match case {
            0 => {
                a.max_rows = 1;
                let mut r = a.units[0].rows[0].clone();
                r.ordinal = 2;
                a.units[0].rows.push(r);
            }
            1 => a.budgets.max_operations = 4,
            2 => a.budgets.max_annotation_rows = 1,
            3 => {
                a.budgets.max_records = 1;
                let mut u = a.units[0].clone();
                u.unit_id = "other".into();
                u.rows[0].ordinal = 2;
                a.units.push(u);
            }
            _ => {
                a.budgets.max_segments_per_query = 1;
                let mut r = a.units[0].rows[0].clone();
                r.ordinal = 2;
                a.units[0].rows.push(r);
            }
        }
        assert!(analyze(&a, &mut || Ok(())).is_err(), "over cap {case}");
    }
    let mut polls = 0;
    for stop in 1..10 {
        polls = 0;
        let result = analyze(&input(), &mut || {
            polls += 1;
            if polls == stop {
                Err(bad("synthetic_control_error"))
            } else {
                Ok(())
            }
        });
        assert_eq!(result.unwrap_err().message, "synthetic_control_error");
    }
    assert!(polls > 0);
}
