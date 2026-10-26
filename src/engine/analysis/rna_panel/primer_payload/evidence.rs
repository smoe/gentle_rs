//! Best-equivalent segment DAGs and distinct read support, not a synthetic mapping winner.

use super::*;

#[cfg(test)]
#[path = "legacy_structure_tests.rs"]
mod legacy_structure_tests;

pub(super) fn target_fraction(expected: Option<usize>, usable: usize) -> Option<usize> {
    match expected {
        Some(n) if usable > 0 && n <= usable => Some((n as u128 * 1000 / usable as u128) as usize),
        _ => None,
    }
}

fn edge(a: &paf::Alignment, head: bool) -> Option<PrimerPayloadEdge> {
    let s = &a.segment;
    let block = if head {
        a.blocks.iter().find(|b| b.query.0 == s.query_interval.0)
    } else {
        a.blocks.iter().find(|b| b.query.1 == s.query_interval.1)
    }?;
    let coordinate = match (s.strand, head) {
        ('+', true) | ('-', false) => block.target.0,
        _ => block.target.1,
    };
    let expected = match (s.strand, head) {
        ('+', true) | ('-', false) => s.target_interval.0,
        _ => s.target_interval.1,
    };
    if coordinate != expected {
        return None;
    }
    Some(PrimerPayloadEdge {
        target_name: s.target_name.clone(),
        coordinate,
        strand: s.strand,
        gene_ids: s.gene_ids.clone(),
        transcript_ids: s.transcript_ids.clone(),
    })
}

pub(super) fn analyze(
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
    let dag = crate::observed_transcript_structures::dag::best_equivalent(
        alignments
            .iter()
            .map(|a| a.segment.query_interval)
            .collect(),
        alignments
            .iter()
            .map(|a| span_bases(&a.segment.aligned_query_intervals))
            .collect(),
        work,
        request.budgets.max_chain_comparisons,
    )?;
    let (weight, prefix, suffix, best, optimal) =
        (dag.weight, dag.prefix, dag.suffix, dag.best, dag.optimal);
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

pub(super) fn context_ranges(
    rows: &[PrimerPayloadReadEvidence],
    _panel: &RnaPrimerPanel,
    lengths: &BTreeMap<String, usize>,
    policy: &PrimerPayloadPolicy,
) -> BTreeSet<reference::ContextKey> {
    rows.iter()
        .flat_map(|row| row.primer_edge_alternatives.iter())
        .flat_map(|e| {
            policy.terminal_seed_lengths.iter().filter_map(move |&n| {
                let length = *lengths.get(&e.target_name)?;
                if e.strand == '+' {
                    e.coordinate
                        .checked_sub(n)
                        .map(|a| (e.target_name.clone(), a, e.coordinate))
                } else {
                    e.coordinate
                        .checked_add(n)
                        .filter(|z| *z <= length)
                        .map(|z| (e.target_name.clone(), e.coordinate, z))
                }
            })
        })
        .collect()
}

pub(super) fn seeds(
    rows: &mut [PrimerPayloadReadEvidence],
    panel: &RnaPrimerPanel,
    contexts: &BTreeMap<reference::ContextKey, Vec<u8>>,
    policy: &PrimerPayloadPolicy,
) -> Result<(), EngineError> {
    for row in rows {
        for edge in &row.primer_edge_alternatives {
            for id in &row.primer_ids {
                let primer = panel
                    .primers
                    .iter()
                    .find(|p| &p.primer_id == id)
                    .ok_or_else(|| bad("primer_payload_seed_primer"))?;
                let core = primer
                    .components
                    .iter()
                    .find(|c| c.role == "target_core")
                    .ok_or_else(|| bad("primer_payload_seed_core"))?;
                let core = &primer.sequence.as_bytes()[core.start..core.end];
                let values = policy
                    .terminal_seed_lengths
                    .iter()
                    .map(|&length| {
                        let key = if edge.strand == '+' {
                            edge.coordinate
                                .checked_sub(length)
                                .map(|a| (edge.target_name.clone(), a, edge.coordinate))
                        } else {
                            edge.coordinate
                                .checked_add(length)
                                .map(|z| (edge.target_name.clone(), edge.coordinate, z))
                        };
                        let context = key.as_ref().and_then(|k| contexts.get(k)).map(|s| {
                            if edge.strand == '+' {
                                s.clone()
                            } else {
                                GentleEngine::reverse_complement_bytes(s)
                            }
                        });
                        let (mismatches, reason) = if core.len() < length {
                            (None, Some("core_shorter_than_requested_seed".into()))
                        } else if let Some(ctx) = context.as_ref().filter(|s| {
                            s.len() >= length
                                && s[s.len() - length..].iter().all(|x| b"ACGT".contains(x))
                        }) {
                            (
                                Some(
                                    core[core.len() - length..]
                                        .iter()
                                        .zip(&ctx[ctx.len() - length..])
                                        .filter(|(a, z)| a != z)
                                        .count(),
                                ),
                                None,
                            )
                        } else {
                            (
                                None,
                                Some("verified_unambiguous_reference_flank_unavailable".into()),
                            )
                        };
                        PrimerPayloadSeed {
                            length,
                            mismatches,
                            reason,
                        }
                    })
                    .collect();
                row.terminal_seeds.insert(
                    serde_json::to_string(&(id, &edge.target_name, edge.coordinate, edge.strand))
                        .map_err(|e| bad(e.to_string()))?,
                    values,
                );
            }
        }
    }
    Ok(())
}

type Base = (String, String, String, String, Vec<String>);
type ReadKey = (String, u64);
fn base(row: &PrimerPayloadReadEvidence) -> Base {
    (
        row.reference_id.clone(),
        row.sample_id.clone(),
        row.barcode.clone(),
        row.priming.clone(),
        row.primer_ids.clone(),
    )
}

pub(super) fn aggregate(
    rows: &[PrimerPayloadReadEvidence],
    prepared: &PrimerPayloadQueries,
    trace: &RnaPrimerTraceReport,
    references: &[String],
    tolerance: usize,
) -> (
    Vec<PrimerPayloadPartnerProfile>,
    Vec<PrimerPayloadBreakpointCluster>,
) {
    let mut groups = BTreeMap::<Base, Vec<&PrimerPayloadReadEvidence>>::new();
    for row in rows {
        groups.entry(base(row)).or_default().push(row);
    }
    let mut excluded = BTreeMap::<Base, BTreeSet<ReadKey>>::new();
    for x in &prepared.unassessed {
        if let Some(unit) = trace
            .files_completed
            .iter()
            .find(|f| f.file_id == x.file_id)
        {
            for reference in references {
                excluded
                    .entry((
                        reference.clone(),
                        unit.input.sample_id.clone(),
                        unit.input.barcode.clone(),
                        unit.input.priming.clone(),
                        x.primer_ids.clone(),
                    ))
                    .or_default()
                    .insert((x.file_id.clone(), x.record_index));
            }
        }
    }
    for key in excluded.keys() {
        groups.entry(key.clone()).or_default();
    }
    let mut profiles = Vec::new();
    let mut clusters = Vec::new();
    for ((reference_id, sample_id, barcode, priming, primer_ids), group) in groups {
        let key = (
            reference_id.clone(),
            sample_id.clone(),
            barcode.clone(),
            priming.clone(),
            primer_ids.clone(),
        );
        let unique = |predicate: &dyn Fn(&PrimerPayloadReadEvidence) -> bool| {
            group
                .iter()
                .filter(|r| predicate(r))
                .map(|r| (r.file_id.clone(), r.record_index))
                .collect::<BTreeSet<_>>()
        };
        let eligible = unique(&|_| true);
        let assessed =
            unique(&|r| !matches!(r.classification.as_str(), "unassessed" | "unresolved"));
        let ambiguous = unique(&|r| r.classification == "ambiguous");
        let unmapped = unique(&|r| r.classification == "producer_reported_unmapped");
        let unassessed =
            unique(&|r| matches!(r.classification.as_str(), "unassessed" | "unresolved"));
        let excluded = excluded.get(&key).map_or(0, BTreeSet::len);
        let mut partners = BTreeMap::<Vec<String>, BTreeSet<ReadKey>>::new();
        if group.is_empty() {
            partners.entry(Vec::new()).or_default();
        }
        for row in &group {
            let support = partners.entry(row.partner_gene_ids.clone()).or_default();
            if !row.partner_gene_ids.is_empty()
                && !matches!(row.classification.as_str(), "unassessed" | "unresolved")
            {
                support.insert((row.file_id.clone(), row.record_index));
            }
        }
        for (partner_gene_ids, supporting) in partners {
            profiles.push(PrimerPayloadPartnerProfile {
                reference_id: reference_id.clone(),
                sample_id: sample_id.clone(),
                barcode: barcode.clone(),
                priming: priming.clone(),
                primer_ids: primer_ids.clone(),
                partner_gene_ids,
                supporting_reads: supporting.len(),
                eligible_reads: eligible.len(),
                assessed_reads: assessed.len(),
                ambiguous_reads: ambiguous.len(),
                unmapped_reads: unmapped.len(),
                unassessed_reads: unassessed.len(),
                excluded_reads: excluded,
                anchor_observations: group.len(),
            });
        }
        let mut edge_groups =
            BTreeMap::<(String, char, String, Vec<String>), Vec<(usize, ReadKey)>>::new();
        for row in group {
            // Alternatives remain per-read evidence; do not turn them into exclusive recurrent breakpoints.
            if !matches!(
                row.classification.as_str(),
                "expected_dominant" | "off_target_dominant" | "mixed"
            ) {
                continue;
            }
            let read = (row.file_id.clone(), row.record_index);
            if row.primer_edge_alternatives.len() == 1 {
                let e = &row.primer_edge_alternatives[0];
                edge_groups
                    .entry((
                        e.target_name.clone(),
                        e.strand,
                        "primer_adjacent".into(),
                        row.partner_gene_ids.clone(),
                    ))
                    .or_default()
                    .push((e.coordinate, read.clone()));
            }
            for boundary in &row.internal_query_boundaries {
                for (ty, e) in [
                    ("internal_left", &boundary.left),
                    ("internal_right", &boundary.right),
                ] {
                    edge_groups
                        .entry((
                            e.target_name.clone(),
                            e.strand,
                            ty.into(),
                            row.partner_gene_ids.clone(),
                        ))
                        .or_default()
                        .push((e.coordinate, read.clone()));
                }
            }
        }
        for ((target_name, strand, edge_type, partner_gene_ids), mut observations) in edge_groups {
            observations.sort();
            observations.dedup();
            let mut start = 0;
            while start < observations.len() {
                let minimum_edge = observations[start].0;
                let mut end = start + 1;
                while end < observations.len() && observations[end].0 - minimum_edge <= tolerance {
                    end += 1;
                }
                let slice = &observations[start..end];
                let maximum_edge = slice.last().unwrap().0;
                let support = slice.iter().map(|(_, r)| r).collect::<BTreeSet<_>>().len();
                clusters.push(PrimerPayloadBreakpointCluster {
                    reference_id: reference_id.clone(),
                    sample_id: sample_id.clone(),
                    barcode: barcode.clone(),
                    priming: priming.clone(),
                    primer_ids: primer_ids.clone(),
                    target_name: target_name.clone(),
                    strand,
                    edge_type: edge_type.clone(),
                    partner_gene_ids: partner_gene_ids.clone(),
                    minimum_edge,
                    maximum_edge,
                    representative_edge: minimum_edge,
                    supporting_reads: support,
                    eligible_reads: eligible.len(),
                    assessed_reads: assessed.len(),
                    ambiguous_reads: ambiguous.len(),
                    unmapped_reads: unmapped.len(),
                    unassessed_reads: unassessed.len(),
                    excluded_reads: excluded,
                });
                start = end;
            }
        }
    }
    profiles.sort_by(|a, z| {
        z.supporting_reads.cmp(&a.supporting_reads).then_with(|| {
            (
                &a.reference_id,
                &a.sample_id,
                &a.barcode,
                &a.priming,
                &a.primer_ids,
                &a.partner_gene_ids,
            )
                .cmp(&(
                    &z.reference_id,
                    &z.sample_id,
                    &z.barcode,
                    &z.priming,
                    &z.primer_ids,
                    &z.partner_gene_ids,
                ))
        })
    });
    (profiles, clusters)
}
