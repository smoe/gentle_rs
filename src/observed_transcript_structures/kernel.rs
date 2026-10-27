//! Bounded, deterministic in-memory observation, record reconciliation and support accounting.

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

fn observe(
    input: &Input,
    query: &QueryEvidence,
    segment: &PrimerPayloadSegment,
    sources: &[comparison::Source<'_>],
    work: &mut usize,
) -> Result<Observation, EngineError> {
    let mut ends = query.ends.clone();
    for constraints in &mut ends {
        constraints.sort();
        constraints.dedup();
        if constraints.is_empty() {
            constraints.push(EndConstraint::Unknown);
        }
    }
    let mut out = Observation {
        observation_id: String::new(),
        record: query.record.clone(),
        query_id: query.query_id.clone(),
        contig: segment.target_name.clone(),
        observed_span: segment.target_interval,
        aligned_query_interval: segment.query_interval,
        payload_query_length: query.query_length,
        alignment_clipped: [
            segment.query_interval.0 != 0,
            segment.query_interval.1 != query.query_length,
        ],
        query_to_reference_strand: segment.strand,
        reported_mapq: segment.mapq,
        ends,
        detail: if query.scope == PrimerPayloadReferenceScope::Transcriptome {
            DetailStatus::TranscriptomeUnassessed
        } else if !segment.alignment_detail_verified {
            DetailStatus::CoarseUnassessed
        } else {
            DetailStatus::Detailed
        },
        blocks: Vec::new(),
        topology: Topology::Unassessed,
        excluded_query_intervals: query.masks.clone(),
        skips: Vec::new(),
        flank_qualification: Vec::new(),
        observed_junction_chain: None,
        comparisons: Vec::new(),
    };
    if segment.alignment_detail_verified {
        let g = geometry::walk(
            segment
                .alignment_operations
                .iter()
                .map(|&(c, n)| (c as u8, n)),
            segment.query_interval,
            segment.target_interval,
            segment.strand,
        )?;
        let intervals = union(
            g.blocks
                .iter()
                .flat_map(|b| unmasked(b.query, query.masks.as_deref().unwrap_or(&[])))
                .collect(),
        );
        if intervals != segment.aligned_query_intervals {
            return Err(bad("observed_structure_stored_support_mismatch"));
        }
        if out.detail == DetailStatus::Detailed {
            out.blocks = g.blocks;
            out.skips = g.skips;
            out.topology = if out.skips.is_empty() {
                Topology::Unspliced
            } else {
                Topology::Spliced
            };
            if out.skips.len() > input.budgets.max_skips_per_chain {
                return Err(bad("observed_structure_skip_budget"));
            }
            for &skip in &out.skips {
                out.flank_qualification.push(comparison::qualify(
                    &out,
                    skip,
                    input.required_flank_bases,
                    work,
                    input.budgets.max_total_work,
                )?);
            }
            if !out.skips.is_empty() {
                out.observed_junction_chain = Some(digest(&(
                    "gentle.observed_junction_chain.v1",
                    &input.reference.assembly,
                    &input.reference.fasta_sha256,
                    &out.contig,
                    &out.skips,
                ))?);
            }
        }
    }
    out.comparisons = comparison::compare(&out, sources, input, work)?;
    out.observation_id = digest(&(
        "gentle.observed_transcript_observation.v1",
        &input.reference,
        &out.record,
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
    ))?;
    Ok(out)
}

fn contains(
    a: &ChainFamily,
    b: &ChainFamily,
    work: &mut usize,
    limit: usize,
) -> Result<bool, EngineError> {
    if a.contig != b.contig || a.skips.len() < b.skips.len() {
        return Ok(false);
    }
    for window in a.skips.windows(b.skips.len()) {
        dag::charge_many(work, limit, b.skips.len())?;
        if window == b.skips {
            return Ok(true);
        }
    }
    Ok(false)
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
    let sources = comparison::sources(input)?;
    let mut work = 0usize;
    let mut observations = BTreeMap::<String, Observation>::new();
    let mut families = BTreeMap::<String, ChainFamily>::new();
    let mut record_paths = BTreeMap::<ReadRecord, (Vec<Vec<String>>, bool, bool)>::new();
    for record in &input.unassessed_records {
        record_paths.insert(record.clone(), (Vec::new(), false, false));
    }
    let mut queries = input.queries.iter().collect::<Vec<_>>();
    queries.sort_by(|a, b| (&a.record, &a.query_id).cmp(&(&b.record, &b.query_id)));
    for query in queries {
        check()?;
        let mut nodes = BTreeMap::new();
        for segment in &query.segments {
            check()?;
            let observation = observe(input, query, segment, &sources, &mut work)?;
            if let Some(chain) = &observation.observed_junction_chain {
                families
                    .entry(chain.clone())
                    .or_insert_with(|| ChainFamily {
                        chain_id: chain.clone(),
                        contig: observation.contig.clone(),
                        skips: observation.skips.clone(),
                        supporting_records: 0,
                        definite_records: 0,
                    });
                if families.len() > input.budgets.max_chain_families {
                    return Err(bad("observed_structure_family_budget"));
                }
            }
            let weight = if observation.detail == DetailStatus::Detailed {
                union(
                    observation
                        .blocks
                        .iter()
                        .flat_map(|b| unmasked(b.query, query.masks.as_deref().unwrap_or(&[])))
                        .collect(),
                )
                .iter()
                .map(|(a, z)| z - a)
                .sum()
            } else {
                0
            };
            nodes.insert(
                (
                    segment.query_interval,
                    segment.target_name.clone(),
                    segment.target_interval,
                    segment.strand,
                    observation.observation_id.clone(),
                ),
                weight,
            );
            observations.insert(observation.observation_id.clone(), observation);
        }
        let nodes = nodes.into_iter().collect::<Vec<_>>();
        let dag = dag::best_equivalent(
            nodes.iter().map(|(k, _)| k.0).collect(),
            nodes.iter().map(|(_, w)| *w).collect(),
            &mut work,
            input.budgets.max_total_work,
        )?;
        let paths = dag.paths(
            input.budgets.max_alternatives_per_record,
            &mut work,
            input.budgets.max_total_work,
        )?;
        let entry = record_paths
            .entry(query.record.clone())
            .or_insert_with(|| (Vec::new(), true, true));
        entry.1 &= !paths.is_empty()
            && query.scope == PrimerPayloadReferenceScope::Genome
            && query.segments.iter().all(|s| s.alignment_detail_verified);
        entry.2 &= query.search_complete;
        for path in paths {
            entry
                .0
                .push(path.iter().map(|&i| nodes[i].0.4.clone()).collect());
        }
        entry.0.sort();
        entry.0.dedup();
        if entry.0.len() > input.budgets.max_alternatives_per_record {
            return Err(bad("observed_structure_alternative_budget"));
        }
    }
    let mut records = Vec::new();
    let mut used_families = BTreeSet::new();
    for (record, (hypotheses, assessed, declared_complete)) in record_paths {
        check()?;
        let complete = assessed && declared_complete;
        let mut possible = BTreeSet::new();
        let mut definite: Option<BTreeSet<String>> = None;
        let mut chain_sets = Vec::new();
        let mut topologies = BTreeSet::new();
        for path in &hypotheses {
            if path.len() > 1 {
                topologies.insert(Topology::Segmented);
            } else if let Some(id) = path.first() {
                topologies.insert(observations[id].topology);
            }
            let chains = path
                .iter()
                .filter_map(|id| observations[id].observed_junction_chain.clone())
                .collect::<BTreeSet<_>>();
            possible.extend(chains.iter().cloned());
            definite = Some(match definite {
                None => chains.clone(),
                Some(old) => old.intersection(&chains).cloned().collect(),
            });
            chain_sets.push(chains);
        }
        let definite = if complete {
            definite.unwrap_or_default()
        } else {
            BTreeSet::new()
        };
        let ids = possible.iter().collect::<Vec<_>>();
        let mut nested = true;
        for (i, a) in ids.iter().enumerate() {
            for b in ids.iter().skip(i + 1) {
                dag::charge(&mut work, input.budgets.max_total_work)?;
                nested &= contains(
                    &families[*a],
                    &families[*b],
                    &mut work,
                    input.budgets.max_total_work,
                )? || contains(
                    &families[*b],
                    &families[*a],
                    &mut work,
                    input.budgets.max_total_work,
                )?;
            }
        }
        let agreement = chain_sets.windows(2).all(|w| w[0] == w[1]);
        let reconciliation = if !complete || hypotheses.is_empty() {
            Reconciliation::Unassessed
        } else if agreement {
            Reconciliation::Agreement
        } else if nested && possible.len() > 1 {
            Reconciliation::ContainedPartial
        } else {
            Reconciliation::Alternatives
        };
        for chain in &possible {
            let family = families.get_mut(chain).expect("observed family");
            family.supporting_records += 1;
            family.definite_records += usize::from(definite.contains(chain));
            used_families.insert(chain.clone());
        }
        records.push(RecordSupport {
            record,
            hypotheses,
            reconciliation,
            topologies: if topologies.is_empty() {
                vec![Topology::Unassessed]
            } else {
                topologies.into_iter().collect()
            },
            search_declared_complete: declared_complete,
            all_query_views_assessed: assessed,
            definite_chains: definite.into_iter().collect(),
            possible_chains: possible.into_iter().collect(),
        });
    }
    Ok(Summary {
        context_sha256: digest(input)?,
        reference: input.reference.clone(),
        annotation_bindings: {
            let mut bindings = input
                .annotation_sources
                .iter()
                .map(|s| AnnotationBinding {
                    source_id: s.source_id.clone(),
                    provider: s.provider.clone(),
                    annotation_sha256: s.annotation_sha256.clone(),
                    reference_id: s.annotation.reference_id.clone(),
                    assembly: s.annotation.assembly.clone(),
                    release: s.annotation.release.clone(),
                    reference_fasta_sha256: s.annotation.reference_fasta_sha256.clone(),
                })
                .collect::<Vec<_>>();
            bindings.sort_by(|a, b| a.source_id.cmp(&b.source_id));
            bindings
        },
        reference_gene: input.reference_gene.clone(),
        gene_resolutions: comparison::resolutions(&sources),
        observations: observations.into_values().collect(),
        records,
        families: families
            .into_iter()
            .filter(|(id, _)| used_families.contains(id))
            .map(|(_, f)| f)
            .collect(),
        comparison_work: work,
    })
}
