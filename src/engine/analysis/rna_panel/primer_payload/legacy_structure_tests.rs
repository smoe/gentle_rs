//! Frozen ce4d33b5 analyzer plus inline synthetic PAF rows: exact Phase 2 byte regression.
//! No real reads, aligner or network. Recreate with the primer_payload_s0_legacy filter.
//! The copied oracle is test-only; production callers share the new geometry/DAG helpers.

use super::super::{
    alignment_tests,
    tests::{limits, prepare, source},
};
use super::*;

fn legacy_analyze(
    query: &PrimerPayloadQuery,
    mut alignments: Vec<paf::Alignment>,
    batch: &PrimerPayloadAlignmentBatch,
    reference: &PrimerPayloadReference,
    cat: &reference::Catalog,
    request: &SummarizePrimerPayloadAlignmentsRequest,
    work: &mut usize,
) -> Result<PrimerPayloadReadEvidence, EngineError> {
    let mut seen = BTreeSet::new();
    alignments
        .retain(|a| seen.insert(serde_json::to_string(&a.segment).expect("serializable segment")));
    alignments.sort_by_key(|a| {
        (
            a.segment.query_interval,
            a.segment.target_name.clone(),
            a.segment.target_interval,
            a.segment.strand,
        )
    });
    let policy = &request.policy;
    let expected = request
        .expected_targets
        .iter()
        .filter(|e| {
            e.reference_id == reference.reference_id && query.primer_ids.contains(&e.primer_id)
        })
        .collect::<Vec<_>>();
    let known = query
        .primer_ids
        .iter()
        .all(|id| expected.iter().any(|e| &e.primer_id == id));
    let n = alignments.len();
    let weight = alignments
        .iter()
        .map(|a| span_bases(&a.segment.aligned_query_intervals))
        .collect::<Vec<_>>();
    let mut prefix = weight.clone();
    let mut suffix = weight.clone();
    for i in 0..n {
        for j in 0..i {
            *work = work
                .checked_add(1)
                .filter(|n| *n <= request.budgets.max_chain_comparisons)
                .ok_or_else(|| bad("primer_payload_chain_comparison_budget"))?;
            if alignments[j].segment.query_interval.1 <= alignments[i].segment.query_interval.0 {
                prefix[i] = prefix[i].max(prefix[j] + weight[i]);
            }
        }
    }
    for i in (0..n).rev() {
        for j in i + 1..n {
            *work = work
                .checked_add(1)
                .filter(|n| *n <= request.budgets.max_chain_comparisons)
                .ok_or_else(|| bad("primer_payload_chain_comparison_budget"))?;
            if alignments[i].segment.query_interval.1 <= alignments[j].segment.query_interval.0 {
                suffix[i] = suffix[i].max(weight[i] + suffix[j]);
            }
        }
    }
    let best = prefix.iter().copied().max().unwrap_or(0);
    let optimal = (0..n)
        .filter(|i| weight[*i] > 0 && prefix[*i] + suffix[*i] - weight[*i] == best)
        .collect::<Vec<_>>();
    let mut reasons = Vec::new();
    let mut genes = BTreeSet::new();
    let mut edges = BTreeSet::new();
    let mut spans = Vec::new();
    let mut boundaries = Vec::new();
    let mut boundary_bytes = 0usize;
    let mut cross = false;
    let mut ambiguous = query.primer_ids.len() != 1;
    let mut insufficient = false;
    let mut coarse = false;
    for &i in &optimal {
        let a = &alignments[i];
        let s = &a.segment;
        genes.extend(s.gene_ids.iter().cloned());
        insufficient |= weight[i] < policy.min_segment_bases
            || s.reported_identity_per_mille < policy.min_identity_per_mille;
        ambiguous |= s.mapq.is_none_or(|q| q < policy.min_mapq || q == 0) || s.gene_ids.len() > 1;
        coarse |= !s.alignment_detail_verified;
        spans.extend(reference::expected_spans(
            a,
            &expected,
            cat,
            reference.scope,
            &query.excluded_query_intervals,
        ));
        if s.query_interval.0 == 0
            && !query
                .excluded_query_intervals
                .first()
                .is_some_and(|(a, _)| *a == 0)
            && prefix[i] == weight[i]
            && let Some(e) = edge(a, true)
        {
            edges.insert(e);
        }
        for &j in optimal.iter().filter(|j| **j > i) {
            *work = work
                .checked_add(1)
                .filter(|n| *n <= request.budgets.max_chain_comparisons)
                .ok_or_else(|| bad("primer_payload_chain_comparison_budget"))?;
            let right = &alignments[j];
            if s.query_interval.1 > right.segment.query_interval.0 {
                if s.gene_ids != right.segment.gene_ids
                    || s.target_name != right.segment.target_name
                    || s.strand != right.segment.strand
                    || a.blocks != right.blocks
                    || a.skips != right.skips
                {
                    ambiguous = true;
                }
                continue;
            }
            if prefix[i] + suffix[j] != best {
                continue;
            }
            if let (Some(left), Some(right_edge)) = (edge(a, false), edge(right, true)) {
                let different = !left.gene_ids.is_empty()
                    && !right_edge.gene_ids.is_empty()
                    && left
                        .gene_ids
                        .iter()
                        .all(|g| !right_edge.gene_ids.contains(g));
                cross |= different
                    && weight[i] >= policy.min_cross_locus_segment_bases
                    && weight[j] >= policy.min_cross_locus_segment_bases
                    && [s, &right.segment].iter().all(|seg| {
                        seg.alignment_detail_verified
                            && seg.reported_identity_per_mille >= policy.min_identity_per_mille
                            && seg.mapq.is_some_and(|q| q > 0 && q >= policy.min_mapq)
                    });
                let boundary = PrimerPayloadBoundary {
                    left_query_end: s.query_interval.1,
                    right_query_start: right.segment.query_interval.0,
                    left,
                    right: right_edge,
                };
                boundary_bytes = boundary_bytes
                    .checked_add(encode(&boundary)?.len())
                    .filter(|n| *n as u64 <= request.budgets.max_output_bytes.min(128_000_000))
                    .ok_or_else(|| bad("primer_payload_boundary_budget"))?;
                boundaries.push(boundary);
            }
        }
    }
    // Mapping-only rows are retained, but they never acquire exact block/edge claims.
    if alignments
        .iter()
        .any(|a| !a.segment.alignment_detail_verified)
    {
        coarse = true;
        reasons.push("mapping_only_alignment_detail_unavailable".into());
    }
    let accounting = batch
        .accounting
        .get(&query.query_id)
        .expect("admitted accounting");
    let expected_bases = if known && !coarse {
        Some(span_bases(&union(spans)))
    } else {
        None
    };
    let fraction = target_fraction(expected_bases, query.usable_bases);
    let expected_genes = expected
        .iter()
        .flat_map(|e| e.gene_ids.iter().cloned())
        .collect::<BTreeSet<_>>();
    let has_expected =
        genes.iter().any(|g| expected_genes.contains(g)) || expected_bases.is_some_and(|n| n > 0);
    let has_foreign = genes.iter().any(|g| !expected_genes.contains(g));
    let classification = if *accounting == PrimerPayloadQueryAccounting::Unassessed {
        "unassessed"
    } else if *accounting == PrimerPayloadQueryAccounting::ProducerReportedUnmapped {
        "producer_reported_unmapped"
    } else if coarse || !known || optimal.is_empty() || insufficient || genes.is_empty() {
        "unresolved"
    } else if ambiguous {
        "ambiguous"
    } else if fraction.is_some_and(|n| n >= policy.dominant_fraction_per_mille)
        && !has_foreign
        && best as u128 * 1000
            >= query.usable_bases as u128 * policy.dominant_fraction_per_mille as u128
    {
        "expected_dominant"
    } else if !has_expected
        && best as u128 * 1000
            >= query.usable_bases as u128 * policy.dominant_fraction_per_mille as u128
    {
        "off_target_dominant"
    } else if has_expected && has_foreign {
        "mixed"
    } else {
        "unresolved"
    };
    if !known {
        reasons.push("expected_target_crosswalk_unavailable".into());
    }
    if alignments.iter().any(|a| a.segment.gene_ids.is_empty()) {
        reasons.push("aligned_region_outside_declared_gene_annotation".into());
    }
    if ambiguous {
        reasons.push("primer_or_mapping_alternatives_retained".into());
    }
    if edges.is_empty() {
        reasons.push("primer_adjacent_exact_edge_unavailable".into());
    }
    if !batch.declared_complete {
        reasons.push("producer_search_declared_incomplete".into());
    }
    reasons.push("producer_hit_limits_do_not_establish_reference_uniqueness".into());
    let mut flags = Vec::new();
    if !query.excluded_query_intervals.is_empty() {
        flags.push("observed_additional_primer_intervals_masked".into());
    }
    if query.primer_ids.len() > 1 {
        flags.push("ambiguous_primer_attribution".into());
    }
    if cross {
        reasons.push("ordered_cross_locus_payload_candidate_not_mechanism".into());
    }
    Ok(PrimerPayloadReadEvidence {
        query_id: query.query_id.clone(),
        batch_id: batch.batch_id.clone(),
        reference_id: reference.reference_id.clone(),
        primer_ids: query.primer_ids.clone(),
        file_id: query.file_id.clone(),
        record_index: query.record_index,
        read_id: query.read_id.clone(),
        sample_id: query.sample_id.clone(),
        barcode: query.barcode.clone(),
        priming: query.priming.clone(),
        usable_bases: query.usable_bases,
        expected_target_bases: expected_bases,
        expected_target_fraction_per_mille: fraction,
        classification: classification.into(),
        partner_gene_ids: genes.into_iter().collect(),
        segments: alignments.into_iter().map(|a| a.segment).collect(),
        primer_edge_alternatives: edges.into_iter().collect(),
        internal_query_boundaries: boundaries,
        terminal_seeds: BTreeMap::new(),
        cross_locus_payload_evidence: cross,
        reasons,
        artifact_flags: flags,
    })
}

#[test]
fn primer_payload_s0_legacy_serialized_rows_and_work_are_byte_identical() {
    let td = tempfile::tempdir().unwrap();
    let (engine, trace) = source(td.path(), true, &"A".repeat(50));
    let prepared = engine
        .prepare_primer_payload_alignment_queries(prepare(trace.clone(), td.path()), &mut |_| true)
        .unwrap();
    let req = alignment_tests::request(td.path(), trace, &prepared);
    let reference = &req.references[0];
    let cat = reference::catalog(reference, &limits()).unwrap();
    for case in 0..10 {
        let mut query = prepared.queries[0].clone();
        if case == 7 {
            query.excluded_query_intervals = vec![(8, 12)];
            query.usable_bases = 46;
        }
        if case == 8 {
            query.primer_ids.push("another_attribution".into());
        }
        let mut batch = req.batches[0].clone();
        batch.accounting =
            BTreeMap::from([(query.query_id.clone(), PrimerPayloadQueryAccounting::Mapped)]);
        batch.declared_complete = case != 9;
        let full = |target: &str, strand: &str| {
            format!("q_synthetic\t50\t0\t50\t{strand}\t{target}\t1000\t0\t50\t50\t50\t60\tcg:Z:50M")
        };
        let rows = match case {
            0 => vec![full("toy_expected","+")],
            1 => vec![full("toy_foreign","+")],
            2 => vec![full("toy_expected","+"),full("toy_foreign","+")],
            3 => vec![
                "q_synthetic\t50\t0\t25\t+\ttoy_expected\t1000\t0\t25\t25\t25\t60\tcg:Z:25M".into(),
                "q_synthetic\t50\t25\t50\t+\ttoy_foreign\t1000\t25\t50\t25\t25\t60\tcg:Z:25M".into()],
            4 => vec!["q_synthetic\t50\t2\t17\t-\ttoy_foreign\t1000\t10\t26\t14\t17\t255\tcg:Z:5M1I4M2D5M\tcs:Z::5+a:4-cc:5\tNM:i:3".into()],
            5 => vec!["q_synthetic\t50\t0\t20\t+\ttoy_expected\t1000\t0\t40\t20\t20\t60\tcg:Z:10M20N10M".into()],
            6 => vec!["q_synthetic\t50\t0\t50\t+\ttoy_foreign\t1000\t0\t50\t50\t50\t0".into()],
            _ => vec![full("toy_expected","-")],
        };
        let alignments = rows
            .iter()
            .map(|row| {
                let target = row.split('\t').nth(5).unwrap();
                let row = row.replacen("q_synthetic", &query.query_id, 1).replacen(
                    "\t1000\t",
                    &format!("\t{}\t", cat.lengths[target]),
                    1,
                );
                paf::row(&row, &query, false, &limits()).unwrap().unwrap().0
            })
            .collect::<Vec<_>>();
        let (mut old_work, mut new_work) = (0, 0);
        let old = legacy_analyze(
            &query,
            alignments.clone(),
            &batch,
            reference,
            &cat,
            &req,
            &mut old_work,
        )
        .unwrap();
        let new = analyze(
            &query,
            alignments,
            &batch,
            reference,
            &cat,
            &req,
            &mut new_work,
        )
        .unwrap();
        assert_eq!(
            serde_json::to_vec(&new).unwrap(),
            serde_json::to_vec(&old).unwrap(),
            "case {case}"
        );
        assert_eq!(old_work, new_work, "work accounting case {case}");
        assert_eq!(
            encode(&new).unwrap(),
            encode(&old).unwrap(),
            "canonical bytes case {case}"
        );
    }
}
