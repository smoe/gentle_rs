//! Legacy admission and byte-preserving projection over the private observation core.

use super::*;

fn admit(input: &Input) -> Result<(), EngineError> {
    let b = &input.budgets;
    b.validate(input.required_flank_bases)?;
    if input.reference.assembly.trim().is_empty()
        || !is_hash(&input.reference.fasta_sha256)
        || input.queries.len() > b.max_queries
        || input.unassessed_records.len() > b.max_records
        || input.annotation_sources.len() > b.max_annotation_rows
        || input.comparison_genes.len() > b.max_annotation_rows
    {
        return Err(bad("observed_structure_admission"));
    }
    let mut query_ids = BTreeSet::new();
    let mut records = BTreeMap::new();
    let mut operations = 0usize;
    let mut lengths = BTreeMap::new();
    for r in &input.unassessed_records {
        if !is_hash(&r.source_sha256)
            || r.read_id.is_empty()
            || records
                .insert((&r.source_sha256, r.record_index), &r.read_id)
                .is_some_and(|id| id != &r.read_id)
            || records.len() > b.max_records
        {
            return Err(bad("observed_structure_unassessed_record"));
        }
    }
    for q in &input.queries {
        let key = (&q.record.source_sha256, q.record.record_index);
        if !is_hash(&q.record.source_sha256)
            || !is_hash(&q.alignment_sha256)
            || q.record.read_id.is_empty()
            || q.query_id.is_empty()
            || !query_ids.insert(&q.query_id)
            || q.query_length == 0
            || q.query_length > 2_000_000
            || q.segments.len() > b.max_segments_per_query
            || records
                .insert(key, &q.record.read_id)
                .is_some_and(|id| id != &q.record.read_id)
            || records.len() > b.max_records
            || q.masks.as_ref().is_some_and(|m| {
                m.len() > b.max_operations
                    || m.iter().any(|&(a, z)| a >= z || z > q.query_length)
                    || !m.windows(2).all(|w| w[0].1 < w[1].0)
            })
        {
            return Err(bad("observed_structure_record_or_query"));
        }
        for (target, &length) in &q.target_lengths {
            if target.is_empty()
                || length == 0
                || lengths
                    .insert(target, length)
                    .is_some_and(|old| old != length)
            {
                return Err(bad("observed_structure_target_length"));
            }
        }
        for s in &q.segments {
            operations = operations
                .checked_add(s.alignment_operations.len())
                .ok_or_else(|| bad("observed_structure_operation_budget"))?;
            if operations > b.max_operations
                || s.query_interval.0 >= s.query_interval.1
                || s.query_interval.1 > q.query_length
                || s.target_interval.0 >= s.target_interval.1
                || q.target_lengths
                    .get(&s.target_name)
                    .is_none_or(|&length| s.target_interval.1 > length)
                || !matches!(s.strand, '+' | '-')
                || s.alignment_detail_verified == s.alignment_operations.is_empty()
                || s.alignment_operations.iter().any(|(c, _)| !c.is_ascii())
            {
                return Err(bad("observed_structure_segment_geometry"));
            }
        }
    }
    for source in &input.annotation_sources {
        for t in &source.annotation.transcripts {
            if lengths
                .get(&t.target_name)
                .is_some_and(|&n| n != t.target_length)
            {
                return Err(bad("observed_structure_reference_length_mismatch"));
            }
        }
    }
    Ok(())
}

struct LegacyIdentity<'a>(&'a Input);

impl core::Identity<ReadRecord> for LegacyIdentity<'_> {
    fn context_id(&self) -> Result<String, EngineError> {
        digest(self.0)
    }

    fn observation_id(
        &self,
        input: &core::Context<'_>,
        query: &core::Query<'_, ReadRecord>,
        segment: &core::Segment<'_>,
        out: &core::CoreObservation<ReadRecord>,
    ) -> Result<String, EngineError> {
        digest(&(
            "gentle.observed_transcript_observation.v1",
            &input.reference,
            &out.key,
            &out.query_id,
            &query.alignment_sha256,
            query.scope,
            &segment.query_interval,
            &out.contig,
            &out.observed_span,
            out.query_to_reference_strand,
            &segment.alignment_operations,
            &out.excluded_query_intervals,
            &out.ends,
        ))
    }
}

/// Interpret supplied in-memory evidence. No references are opened or discovered.
///
/// Hashes bind declared inputs; this is internal consistency, not authenticated
/// provenance, alignment-search replay, natural-end detection or a protein call.
pub fn summarize(input: &Input) -> Result<Summary, EngineError> {
    summarize_with_control(input, &mut || Ok(()))
}

/// Cooperative caller control between bounded query/segment/record units; no I/O here.
pub fn summarize_with_control(
    input: &Input,
    check: &mut dyn FnMut() -> Result<(), EngineError>,
) -> Result<Summary, EngineError> {
    check()?;
    admit(input)?;
    let context = core::Context {
        reference: &input.reference,
        reference_gene: &input.reference_gene,
        comparison_genes: &input.comparison_genes,
        annotation_sources: &input.annotation_sources,
        required_flank_bases: input.required_flank_bases,
        budgets: &input.budgets,
    };
    let queries = input
        .queries
        .iter()
        .map(|q| core::Query {
            key: &q.record,
            alignment_sha256: &q.alignment_sha256,
            query_id: &q.query_id,
            query_length: q.query_length,
            scope: q.scope,
            masks: q.masks.as_deref(),
            ends: &q.ends,
            search_complete: q.search_complete,
            segments: q
                .segments
                .iter()
                .map(|s| core::Segment {
                    query_interval: s.query_interval,
                    target_name: &s.target_name,
                    target_interval: s.target_interval,
                    strand: s.strand,
                    mapq: s.mapq,
                    alignment_detail_verified: s.alignment_detail_verified,
                    alignment_operations: &s.alignment_operations,
                    stored_support: Some(&s.aligned_query_intervals),
                    source_ordinal: None,
                    clipping: None,
                })
                .collect(),
        })
        .collect::<Vec<_>>();
    let result = core::run(
        &context,
        &queries,
        &input.unassessed_records,
        &LegacyIdentity(input),
        0,
        check,
    )?;
    Ok(Summary {
        context_sha256: result.context_sha256,
        reference: result.reference,
        annotation_bindings: result.annotation_bindings,
        reference_gene: result.reference_gene,
        gene_resolutions: result.gene_resolutions,
        observations: result
            .observations
            .into_iter()
            .map(|o| Observation {
                observation_id: o.observation_id,
                record: o.key,
                query_id: o.query_id,
                contig: o.contig,
                observed_span: o.observed_span,
                aligned_query_interval: o.aligned_query_interval,
                payload_query_length: o.query_axis_length,
                alignment_clipped: o.alignment_clipped,
                query_to_reference_strand: o.query_to_reference_strand,
                reported_mapq: o.reported_mapq,
                ends: o.ends,
                detail: o.detail,
                topology: o.topology,
                blocks: o.blocks,
                excluded_query_intervals: o.excluded_query_intervals,
                skips: o.skips,
                flank_qualification: o.flank_qualification,
                observed_junction_chain: o.observed_junction_chain,
                comparisons: o.comparisons,
            })
            .collect(),
        records: result
            .supports
            .into_iter()
            .map(|r| RecordSupport {
                record: r.key,
                hypotheses: r.hypotheses,
                reconciliation: r.reconciliation,
                topologies: r.topologies,
                search_declared_complete: r.search_declared_complete,
                all_query_views_assessed: r.all_query_views_assessed,
                definite_chains: r.definite_chains,
                possible_chains: r.possible_chains,
            })
            .collect(),
        families: result
            .families
            .into_iter()
            .map(|f| ChainFamily {
                chain_id: f.chain_id,
                contig: f.contig,
                skips: f.skips,
                supporting_records: f.supporting_keys,
                definite_records: f.definite_keys,
            })
            .collect(),
        comparison_work: result.comparison_work,
    })
}
