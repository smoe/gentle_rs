//! Hand-crafted synthetic loci, operations and records; no real reads or primer sequences.
//!
//! Recreate by running the observed_transcript_structures library-test filter.
//! Fixtures are inline and consumed only by the pure S0 kernel/helper regressions.

use super::*;
use crate::engine::{PRIMER_PAYLOAD_ANNOTATION_SCHEMA, PrimerPayloadTranscript};

#[test]
fn observed_annotation_preparation_preserves_legacy_touching_exon_admission() {
    let mut a = input();
    a.annotation_sources[0].annotation.transcripts[0].genomic_exons = vec![(100, 110), (110, 120)];
    let context = core::Context {
        reference: &a.reference,
        reference_gene: &a.reference_gene,
        comparison_genes: &a.comparison_genes,
        annotation_sources: &a.annotation_sources,
        required_flank_bases: a.required_flank_bases,
        budgets: &a.budgets,
    };
    // Directly test the pre-existing validator, not a new preparation helper.
    assert!(comparison::sources(&context).is_ok());
    assert_eq!(
        a.annotation_sources[0].annotation.transcripts[0].genomic_exons,
        [(100, 110), (110, 120)]
    );
}

#[path = "preservation_tests.rs"]
mod preservation;

fn sha(c: char) -> String {
    format!("sha256:{}", c.to_string().repeat(64))
}

fn tx(
    id: &str,
    gene: &str,
    contig: &str,
    strand: char,
    exons: Vec<Interval>,
) -> PrimerPayloadTranscript {
    PrimerPayloadTranscript {
        transcript_id: id.into(),
        gene_ids: vec![gene.into()],
        target_name: contig.into(),
        target_length: 1000,
        genomic_contig: Some(contig.into()),
        genomic_strand: Some(strand),
        genomic_exons: exons,
    }
}

fn segment(
    start: usize,
    ops: Vec<(char, usize)>,
    strand: char,
    masks: &[Interval],
) -> PrimerPayloadSegment {
    let qlen = ops
        .iter()
        .filter(|(c, _)| "M=XI".contains(*c))
        .map(|(_, n)| n)
        .sum::<usize>();
    let tlen = ops
        .iter()
        .filter(|(c, _)| "M=XDN".contains(*c))
        .map(|(_, n)| n)
        .sum::<usize>();
    let g = geometry::walk(
        ops.iter().map(|&(c, n)| (c as u8, n)),
        (0, qlen),
        (start, start + tlen),
        strand,
    )
    .unwrap();
    PrimerPayloadSegment {
        query_interval: (0, qlen),
        target_name: "contig_a".into(),
        target_interval: (start, start + tlen),
        strand,
        mapq: Some(60),
        reported_role: None,
        reported_identity_per_mille: 1000,
        edits: None,
        alignment_detail_verified: true,
        alignment_operations: ops,
        aligned_query_intervals: union(
            g.blocks
                .iter()
                .flat_map(|b| unmasked(b.query, masks))
                .collect(),
        ),
        gene_ids: vec!["untrusted_whole_span_label".into()],
        transcript_ids: vec!["not_a_match".into()],
        splice_edges: vec![],
    }
}

fn input() -> Input {
    let annotation = PrimerPayloadAnnotation {
        schema: PRIMER_PAYLOAD_ANNOTATION_SCHEMA.into(),
        reference_id: "coordinate_alias".into(),
        assembly: "toy_assembly".into(),
        release: "toy_release".into(),
        reference_fasta_sha256: sha('1'),
        coordinate_convention: "zero_based_half_open".into(),
        transcripts: vec![
            tx(
                "a.1",
                "gene_a",
                "contig_a",
                '+',
                vec![(100, 110), (120, 130), (140, 150)],
            ),
            tx(
                "b.2",
                "gene_b",
                "contig_b",
                '-',
                vec![(200, 210), (225, 235), (250, 260)],
            ),
            tx(
                "antisense.1",
                "gene_antisense",
                "contig_a",
                '-',
                vec![(100, 110), (120, 130), (140, 150)],
            ),
            tx("nested.1", "gene_nested", "contig_a", '+', vec![(112, 118)]),
        ],
    };
    Input {
        unassessed_records: vec![],
        reference: ReferenceBinding {
            assembly: "toy_assembly".into(),
            fasta_sha256: sha('1'),
        },
        reference_gene: ReferenceGene {
            source_id: "source_a".into(),
            gene_id: "gene_a".into(),
        },
        comparison_genes: vec![],
        annotation_sources: vec![AnnotationSource {
            source_id: "source_a".into(),
            provider: Some("synthetic".into()),
            annotation_sha256: sha('2'),
            annotation,
        }],
        queries: vec![QueryEvidence {
            record: ReadRecord {
                source_sha256: sha('3'),
                record_index: 0,
                read_id: "reused_name".into(),
            },
            alignment_sha256: sha('4'),
            query_id: "query_0".into(),
            query_length: 30,
            scope: PrimerPayloadReferenceScope::Genome,
            target_lengths: BTreeMap::from([("contig_a".into(), 1000), ("contig_b".into(), 1000)]),
            masks: Some(vec![]),
            ends: [
                vec![EndConstraint::PrimerBounded],
                vec![EndConstraint::ReadEndUnverified],
            ],
            search_complete: true,
            segments: vec![segment(
                100,
                vec![('M', 10), ('N', 10), ('M', 10), ('N', 10), ('M', 10)],
                '+',
                &[],
            )],
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
    }
}

#[test]
fn observed_transcripts_missing_anchor_cannot_promote_an_assessed_anchor() {
    let mut a = input();
    let record = a.queries[0].record.clone();
    assert_eq!(summarize(&a).unwrap().records[0].definite_chains.len(), 1);
    a.unassessed_records.push(record.clone());
    let summary = summarize(&a).unwrap();
    assert_eq!(summary.records[0].record, record);
    assert!(!summary.records[0].all_query_views_assessed);
    assert!(summary.records[0].definite_chains.is_empty());
    assert_eq!(summary.families[0].supporting_records, 1);
    assert_eq!(summary.families[0].definite_records, 0);
    assert!(
        !summary.observations[0].comparisons[0]
            .junctions
            .as_ref()
            .unwrap()[0]
            .annotated_in_models
            .is_empty()
    );
}

fn one(result: &Summary) -> &Observation {
    assert_eq!(result.observations.len(), 1);
    &result.observations[0]
}
fn model(observation: &Observation) -> &TranscriptMatch {
    &observation.comparisons[0].models[0]
}

#[test]
fn observed_structures_gene_is_variable_on_two_loci_and_opposite_annotation_strands() {
    let a = summarize(&input()).unwrap();
    assert_eq!(model(one(&a)).relation, ModelRelation::ExactFullChain);
    let mut b = input();
    b.reference_gene.gene_id = "gene_b".into();
    let mut s = segment(
        200,
        vec![('M', 10), ('N', 15), ('M', 10), ('N', 15), ('M', 10)],
        '-',
        &[],
    );
    s.target_name = "contig_b".into();
    b.queries[0].segments = vec![s];
    let b = summarize(&b).unwrap();
    assert_eq!(model(one(&b)).model_strand, '-');
    assert_eq!(model(one(&b)).relation, ModelRelation::ExactFullChain);
    assert_ne!(a.families[0].chain_id, b.families[0].chain_id);
}

#[test]
fn observed_structures_selected_antisense_gene_cannot_change_physical_identity() {
    let a = input();
    let first = summarize(&a).unwrap();
    let mut b = a.clone();
    b.reference_gene.gene_id = "gene_antisense".into();
    let second = summarize(&b).unwrap();
    assert_eq!(first.families[0].chain_id, second.families[0].chain_id);
    assert_eq!(one(&first).observation_id, one(&second).observation_id);
    assert_ne!(first.context_sha256, second.context_sha256);
    assert_eq!(model(one(&second)).model_strand, '-');
    assert_eq!(one(&second).query_to_reference_strand, '+');
    let json = canonical_json_string(&second).unwrap();
    assert!(!json.contains("protein"));
    assert!(!json.contains("isoform"));
}

#[test]
fn observed_structures_missing_ambiguous_and_unmapped_genes_are_typed() {
    let mut a = input();
    a.reference_gene.gene_id = "absent".into();
    assert_eq!(
        summarize(&a).unwrap().gene_resolutions[0].resolution,
        GeneResolution::MissingGene
    );
    a.reference_gene.source_id = "absent_source".into();
    assert!(
        summarize(&a)
            .unwrap()
            .gene_resolutions
            .iter()
            .all(|s| s.resolution != GeneResolution::Resolved)
    );
    let mut a = input();
    a.annotation_sources[0].annotation.transcripts.push(tx(
        "other.1",
        "gene_a",
        "contig_b",
        '+',
        vec![(200, 210)],
    ));
    assert_eq!(
        summarize(&a).unwrap().gene_resolutions[0].resolution,
        GeneResolution::AmbiguousLocus
    );
}

#[test]
fn observed_structures_sources_need_explicit_gene_ids_and_keep_multi_gene_membership() {
    let mut a = input();
    a.annotation_sources[0].annotation.transcripts[0]
        .gene_ids
        .push("gene_shared".into());
    let mut other = a.annotation_sources[0].clone();
    other.source_id = "source_b".into();
    other.annotation.transcripts[0].gene_ids = vec!["different_source_id".into()];
    a.annotation_sources.push(other);
    let r = summarize(&a).unwrap();
    assert_eq!(
        one(&r).comparisons[1].resolution,
        GeneResolution::UnmappedSource
    );
    assert_eq!(model(one(&r)).gene_ids.len(), 2);
    a.comparison_genes.push(ReferenceGene {
        source_id: "source_b".into(),
        gene_id: "different_source_id".into(),
    });
    let r = summarize(&a).unwrap();
    assert_eq!(
        one(&r).comparisons[1].models[0].relation,
        ModelRelation::ExactFullChain
    );
    assert_eq!(r.families[0].supporting_records, 1);
}

#[test]
fn observed_structures_intronic_blocks_ignore_phase2_whole_span_labels() {
    let mut a = input();
    a.queries[0].segments = vec![segment(112, vec![('M', 6)], '+', &[])];
    a.queries[0].query_length = 6;
    let r = summarize(&a).unwrap();
    assert_eq!(model(one(&r)).relation, ModelRelation::Intronic);
    assert!(r.families.is_empty());
    assert!(!one(&r).observation_id.is_empty());
}

#[test]
fn observed_structures_partial_chain_and_primer_bounded_ends_do_not_become_full_rna() {
    let mut a = input();
    let full = a.queries[0].clone();
    a.queries[0].query_id = "partial".into();
    a.queries[0].query_length = 20;
    a.queries[0].segments = vec![segment(
        120,
        vec![('M', 10), ('N', 10), ('M', 10)],
        '+',
        &[],
    )];
    a.queries.push(full);
    let r = summarize(&a).unwrap();
    assert_eq!(r.families.len(), 2);
    assert_eq!(r.records.len(), 1);
    assert_eq!(
        r.records[0].reconciliation,
        Reconciliation::ContainedPartial
    );
    assert!(
        r.observations
            .iter()
            .all(|o| o.ends[0] == vec![EndConstraint::PrimerBounded])
    );
    assert!(
        r.observations
            .iter()
            .any(|o| model(o).relation == ModelRelation::ExactPartialChain)
    );
}

#[test]
fn observed_structures_repeat_anchors_count_once_but_repeated_names_on_distinct_records_do_not() {
    let mut a = input();
    let mut duplicate = a.queries[0].clone();
    duplicate.query_id = "other_anchor".into();
    a.queries.push(duplicate);
    let r = summarize(&a).unwrap();
    assert_eq!(r.families[0].supporting_records, 1);
    assert_eq!(r.records[0].reconciliation, Reconciliation::Agreement);
    a.queries[1].record.record_index = 1;
    let r = summarize(&a).unwrap();
    assert_eq!(r.families[0].supporting_records, 2);
    assert_eq!(r.families[0].definite_records, 2);
}

#[test]
fn observed_structures_conflicting_anchors_preserve_possible_support_without_union() {
    let mut a = input();
    let mut alternate = a.queries[0].clone();
    alternate.query_id = "alternative".into();
    alternate.segments = vec![segment(
        100,
        vec![('M', 11), ('N', 9), ('M', 10), ('N', 10), ('M', 10)],
        '+',
        &[],
    )];
    alternate.query_length = 31;
    a.queries.push(alternate);
    let r = summarize(&a).unwrap();
    assert_eq!(r.families.len(), 2);
    assert_eq!(r.records[0].reconciliation, Reconciliation::Alternatives);
    assert!(
        r.families
            .iter()
            .all(|f| f.supporting_records == 1 && f.definite_records == 0)
    );
}

#[test]
fn observed_structures_primer_masks_and_missing_masks_have_different_qualification() {
    let mut a = input();
    a.queries[0].masks = Some(vec![(8, 12)]);
    a.queries[0].segments = vec![segment(
        100,
        vec![('M', 10), ('N', 10), ('M', 10), ('N', 10), ('M', 10)],
        '+',
        &[(8, 12)],
    )];
    let r = summarize(&a).unwrap();
    assert_eq!(
        one(&r).flank_qualification[0],
        Qualification::InsufficientSupport
    );
    let mut a = input();
    a.queries[0].masks = None;
    let r = summarize(&a).unwrap();
    assert!(
        one(&r)
            .flank_qualification
            .iter()
            .all(|q| *q == Qualification::Unassessed)
    );
    assert_eq!(model(one(&r)).retained_intron_supported_bases, None);
    assert_eq!(model(one(&r)).relation, ModelRelation::ExactFullChain);
}

#[test]
fn observed_structures_deletions_are_not_supported_retained_introns() {
    let mut a = input();
    a.queries[0].segments = vec![segment(
        100,
        vec![('M', 10), ('D', 10), ('M', 10), ('D', 10), ('M', 10)],
        '+',
        &[],
    )];
    let r = summarize(&a).unwrap();
    assert!(r.families.is_empty());
    assert_eq!(model(one(&r)).retained_intron_supported_bases, Some(0));
    a.queries[0].segments = vec![segment(100, vec![('M', 50)], '+', &[])];
    a.queries[0].query_length = 50;
    let r = summarize(&a).unwrap();
    assert!(r.families.is_empty());
    assert_eq!(model(one(&r)).retained_intron_supported_bases, Some(20));
    assert_eq!(model(one(&r)).relation, ModelRelation::ExonIntronBoundary);
}

#[test]
fn observed_structures_coarse_transcriptome_and_incomplete_search_are_not_complete_negatives() {
    let mut a = input();
    a.queries[0].segments[0].alignment_detail_verified = false;
    a.queries[0].segments[0].alignment_operations.clear();
    a.queries[0].segments[0].aligned_query_intervals.clear();
    let r = summarize(&a).unwrap();
    assert_eq!(one(&r).detail, DetailStatus::CoarseUnassessed);
    assert!(r.families.is_empty());
    let mut a = input();
    a.queries[0].scope = PrimerPayloadReferenceScope::Transcriptome;
    let r = summarize(&a).unwrap();
    assert_eq!(one(&r).detail, DetailStatus::TranscriptomeUnassessed);
    assert!(r.families.is_empty());
    let mut a = input();
    a.queries[0].search_complete = false;
    let r = summarize(&a).unwrap();
    assert_eq!(r.families[0].definite_records, 0);
    assert_eq!(r.records[0].reconciliation, Reconciliation::Unassessed);
}

#[test]
fn observed_structures_known_recombination_is_distinct_from_unannotated_junction() {
    let mut a = input();
    a.annotation_sources[0].annotation.transcripts.push(tx(
        "a_alt.1",
        "gene_a",
        "contig_a",
        '+',
        vec![(120, 140), (150, 160)],
    ));
    a.queries[0].segments = vec![segment(
        100,
        vec![('M', 10), ('N', 10), ('M', 20), ('N', 10), ('M', 10)],
        '+',
        &[],
    )];
    a.queries[0].query_length = 40;
    let r = summarize(&a).unwrap();
    assert_eq!(
        one(&r).comparisons[0].junction_relation,
        Some(JunctionRelation::AllJunctionsAnnotatedNoSingleModelContainsRun)
    );
    a.queries[0].segments = vec![segment(
        100,
        vec![('M', 11), ('N', 9), ('M', 20), ('N', 10), ('M', 10)],
        '+',
        &[],
    )];
    a.queries[0].query_length = 41;
    let r = summarize(&a).unwrap();
    assert_eq!(
        one(&r).comparisons[0].junction_relation,
        Some(JunctionRelation::ContainsUnannotatedJunction)
    );
}

#[test]
fn observed_structures_exact_skips_with_intronic_terminal_extension_are_not_exact_model() {
    let mut a = input();
    a.queries[0].segments = vec![segment(
        99,
        vec![('M', 11), ('N', 10), ('M', 10), ('N', 10), ('M', 11)],
        '+',
        &[],
    )];
    a.queries[0].query_length = 32;
    let r = summarize(&a).unwrap();
    assert_eq!(
        model(one(&r)).relation,
        ModelRelation::TerminalCoverageDivergent
    );
}

#[test]
fn observed_structures_resource_bindings_and_stored_geometry_fail_closed() {
    for kind in 0..8 {
        let mut a = input();
        match kind {
            0 => a.annotation_sources[0].annotation.assembly = "different".into(),
            1 => a.queries[0].segments[0].aligned_query_intervals = vec![(0, 1)],
            2 => a.budgets.max_skips_per_chain = 1,
            3 => a.budgets.max_total_work = 1,
            4 => a.budgets.max_operations = 1,
            5 => a.queries[0].segments[0].target_interval.1 = 149,
            6 => a.queries[0].masks = Some(vec![(8, 12), (11, 15)]),
            _ => a.annotation_sources[0].annotation.transcripts[0].target_length = 999,
        }
        assert!(summarize(&a).is_err(), "case {kind}");
    }
}

#[test]
fn observed_structures_partial_run_terminal_blocks_must_belong_to_adjacent_exons() {
    let mut a = input();
    a.queries[0].segments = vec![segment(
        100,
        vec![('M', 5), ('D', 25), ('N', 10), ('M', 10)],
        '+',
        &[],
    )];
    a.queries[0].query_length = 15;
    let r = summarize(&a).unwrap();
    assert_eq!(
        model(one(&r)).relation,
        ModelRelation::TerminalCoverageDivergent
    );
}

#[test]
fn observed_structures_segment_dag_does_not_invent_a_skip_between_records() {
    let mut a = input();
    let left = segment(100, vec![('M', 10)], '+', &[]);
    let mut right = segment(120, vec![('M', 10)], '+', &[]);
    right.query_interval = (10, 20);
    right.aligned_query_intervals = vec![(10, 20)];
    let mut alternative = segment(200, vec![('M', 20)], '-', &[]);
    alternative.target_name = "contig_b".into();
    a.queries[0].query_length = 20;
    a.queries[0].segments = vec![left, right, alternative];
    let r = summarize(&a).unwrap();
    assert_eq!(r.records[0].hypotheses.len(), 2);
    assert!(r.records[0].topologies.contains(&Topology::Segmented));
    assert!(r.families.is_empty());
}

#[test]
fn observed_structures_alignment_receipts_bind_observations_but_not_physical_chains() {
    let a = input();
    let first = summarize(&a).unwrap();
    let mut b = a.clone();
    b.queries[0].alignment_sha256 = sha('5');
    let second = summarize(&b).unwrap();
    assert_ne!(one(&first).observation_id, one(&second).observation_id);
    assert_eq!(first.families[0].chain_id, second.families[0].chain_id);
    assert_eq!(first.annotation_bindings[0].annotation_sha256, sha('2'));
    assert_eq!(one(&first).reported_mapq, Some(60));
}

#[test]
fn observed_structures_alignment_clipping_is_separate_from_payload_end_constraints() {
    let mut a = input();
    a.queries[0].query_length = 32;
    a.queries[0].segments[0].query_interval = (1, 31);
    a.queries[0].segments[0].aligned_query_intervals = vec![(1, 31)];
    let r = summarize(&a).unwrap();
    assert_eq!(one(&r).alignment_clipped, [true, true]);
    assert_eq!(one(&r).aligned_query_interval, (1, 31));
    assert_eq!(one(&r).ends[0], vec![EndConstraint::PrimerBounded]);
    assert_eq!(model(one(&r)).relation, ModelRelation::ExactFullChain);
}

#[test]
fn observed_structures_rejects_unicode_operation_aliases_and_disagreeing_record_ids() {
    let mut a = input();
    a.queries[0].segments[0].alignment_operations[0].0 = '\u{014d}';
    assert!(summarize(&a).is_err());
    let mut a = input();
    let mut q = a.queries[0].clone();
    q.query_id = "other".into();
    q.record.read_id = "contradiction".into();
    a.queries.push(q);
    assert!(summarize(&a).is_err());
}

#[test]
fn observed_structures_input_order_preserves_physical_results_and_full_digests() {
    let mut a = input();
    let mut q = a.queries[0].clone();
    q.query_id = "other".into();
    q.record.record_index = 1;
    a.queries.push(q);
    let first = summarize(&a).unwrap();
    a.queries.reverse();
    a.annotation_sources[0].annotation.transcripts.reverse();
    let second = summarize(&a).unwrap();
    assert_eq!(
        canonical_json_string(&first.families).unwrap(),
        canonical_json_string(&second.families).unwrap()
    );
    assert_eq!(
        canonical_json_string(&first.observations).unwrap(),
        canonical_json_string(&second.observations).unwrap()
    );
    assert_eq!(first.families[0].chain_id.len(), 71);
}

#[test]
fn observed_structures_dag_matches_exhaustive_legal_paths_not_node_union() {
    for pattern in 0..32usize {
        let intervals = (0..6)
            .map(|i| (i / 2, i / 2 + 1 + ((pattern >> i) & 1)))
            .collect::<Vec<_>>();
        let mut intervals = intervals;
        intervals.sort();
        let weights = intervals
            .iter()
            .enumerate()
            .map(|(i, (a, z))| if (pattern >> i) & 1 == 0 { z - a } else { 0 })
            .collect::<Vec<_>>();
        let mut candidates = Vec::new();
        let mut best = 0;
        for mask in 1..64usize {
            let path = (0..6).filter(|&i| (mask >> i) & 1 == 1).collect::<Vec<_>>();
            if path.iter().any(|&i| weights[i] == 0)
                || !path
                    .windows(2)
                    .all(|w| intervals[w[0]].1 <= intervals[w[1]].0)
            {
                continue;
            }
            let score = path.iter().map(|&i| weights[i]).sum::<usize>();
            if score > best {
                best = score;
                candidates.clear();
            }
            if score == best {
                candidates.push(path);
            }
        }
        candidates.sort();
        let mut work = 0;
        let d = dag::best_equivalent(intervals, weights, &mut work, 100000).unwrap();
        assert_eq!(work, 30);
        assert_eq!(d.best, best);
        assert_eq!(d.paths(100, &mut work, 100000).unwrap(), candidates);
    }
}

#[test]
fn observed_structures_dag_alternative_exhaustion_is_not_a_truncated_success() {
    let mut work = 0;
    let d = dag::best_equivalent(
        vec![(0, 10), (0, 10), (10, 20), (10, 20)],
        vec![10; 4],
        &mut work,
        10000,
    )
    .unwrap();
    assert_eq!(d.paths(4, &mut work, 10000).unwrap().len(), 4);
    assert!(d.paths(3, &mut work, 10000).is_err());
}

#[test]
fn observed_structures_geometry_matches_basewise_oracle_in_both_orientations() {
    let alphabet = [('M', 2), ('=', 2), ('X', 1), ('I', 1), ('D', 2), ('N', 4)];
    for a in alphabet {
        for b in alphabet {
            for c in alphabet {
                let ops = [a, b, c];
                let qlen = ops
                    .iter()
                    .filter(|(c, _)| "M=XI".contains(*c))
                    .map(|(_, n)| n)
                    .sum::<usize>();
                let tlen = ops
                    .iter()
                    .filter(|(c, _)| "M=XDN".contains(*c))
                    .map(|(_, n)| n)
                    .sum::<usize>();
                if qlen == 0 || tlen == 0 {
                    continue;
                }
                for strand in ['+', '-'] {
                    let (mut q, mut t) = (0, 20);
                    let mut expected = Vec::new();
                    let mut skips = Vec::new();
                    for &(kind, n) in &ops {
                        if kind == 'N' {
                            skips.push((t, t + n));
                        }
                        for _ in 0..n {
                            if "M=X".contains(kind) {
                                expected.push((
                                    if strand == '+' {
                                        3 + q
                                    } else {
                                        3 + qlen - 1 - q
                                    },
                                    t,
                                ));
                            }
                            if "M=XI".contains(kind) {
                                q += 1;
                            }
                            if "M=XDN".contains(kind) {
                                t += 1;
                            }
                        }
                    }
                    expected.sort();
                    let g = geometry::walk(
                        ops.iter().map(|&(c, n)| (c as u8, n)),
                        (3, 3 + qlen),
                        (20, 20 + tlen),
                        strand,
                    )
                    .unwrap();
                    let mut actual = Vec::new();
                    for block in &g.blocks {
                        for query in block.query.0..block.query.1 {
                            actual.push((
                                query,
                                if strand == '+' {
                                    block.target.0 + query - block.query.0
                                } else {
                                    block.target.1 - 1 - (query - block.query.0)
                                },
                            ));
                        }
                    }
                    actual.sort();
                    assert_eq!(actual, expected);
                    assert_eq!(g.skips, skips);
                }
            }
        }
    }
}
