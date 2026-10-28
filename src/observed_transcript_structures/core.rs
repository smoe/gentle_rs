//! Private format-neutral observation and identity seam; no I/O or public wire schema.

use super::*;

pub(super) struct Context<'a> {
    pub reference: &'a ReferenceBinding,
    pub reference_gene: &'a ReferenceGene,
    pub comparison_genes: &'a [ReferenceGene],
    pub annotation_sources: &'a [AnnotationSource],
    pub required_flank_bases: usize,
    pub budgets: &'a Budgets,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub(super) struct Clipping {
    pub hard: [usize; 2],
    pub soft: [usize; 2],
    pub padding_operations: usize,
}

pub(super) struct Segment<'a> {
    pub query_interval: Interval,
    pub target_name: &'a str,
    pub target_interval: Interval,
    pub strand: char,
    pub mapq: Option<u8>,
    pub alignment_detail_verified: bool,
    pub alignment_operations: &'a [(char, usize)],
    pub stored_support: Option<&'a [Interval]>,
    pub source_ordinal: Option<u64>,
    pub clipping: Option<Clipping>,
}

pub(super) struct Query<'a, R> {
    pub key: &'a R,
    pub alignment_sha256: &'a str,
    pub query_id: &'a str,
    pub query_length: usize,
    pub scope: PrimerPayloadReferenceScope,
    pub masks: Option<&'a [Interval]>,
    pub ends: &'a [Vec<EndConstraint>; 2],
    pub search_complete: bool,
    pub segments: Vec<Segment<'a>>,
}

pub(super) trait Identity<R> {
    fn context_id(&self) -> Result<String, EngineError>;
    fn observation_id(
        &self,
        context: &Context<'_>,
        query: &Query<'_, R>,
        segment: &Segment<'_>,
        observation: &CoreObservation<R>,
    ) -> Result<String, EngineError>;
}

#[derive(Debug, Clone, Serialize)]
pub(super) struct CoreObservation<R> {
    pub observation_id: String,
    pub key: R,
    pub query_id: String,
    pub contig: String,
    pub observed_span: Interval,
    pub aligned_query_interval: Interval,
    pub query_axis_length: usize,
    pub alignment_clipped: [bool; 2],
    pub query_to_reference_strand: char,
    pub reported_mapq: Option<u8>,
    pub ends: [Vec<EndConstraint>; 2],
    pub detail: DetailStatus,
    pub topology: Topology,
    pub blocks: Vec<SupportedBlock>,
    pub excluded_query_intervals: Option<Vec<Interval>>,
    pub skips: Vec<Interval>,
    pub flank_qualification: Vec<Qualification>,
    pub observed_junction_chain: Option<String>,
    pub comparisons: Vec<SourceComparison>,
    pub source_ordinal: Option<u64>,
    pub clipping: Option<Clipping>,
}

#[derive(Debug, Serialize)]
pub(super) struct CoreSupport<R> {
    pub key: R,
    pub hypotheses: Vec<Vec<String>>,
    pub reconciliation: Reconciliation,
    pub topologies: Vec<Topology>,
    pub search_declared_complete: bool,
    pub all_query_views_assessed: bool,
    pub definite_chains: Vec<String>,
    pub possible_chains: Vec<String>,
}

#[derive(Debug, Serialize)]
pub(super) struct CoreFamily {
    pub chain_id: String,
    pub contig: String,
    pub skips: Vec<Interval>,
    pub supporting_keys: usize,
    pub definite_keys: usize,
}

#[derive(Debug, Serialize)]
pub(super) struct CoreSummary<R> {
    pub context_sha256: String,
    pub reference: ReferenceBinding,
    pub annotation_bindings: Vec<AnnotationBinding>,
    pub reference_gene: ReferenceGene,
    pub gene_resolutions: Vec<SourceComparison>,
    pub observations: Vec<CoreObservation<R>>,
    pub supports: Vec<CoreSupport<R>>,
    pub families: Vec<CoreFamily>,
    pub comparison_work: usize,
}

fn observe<R: Clone>(
    input: &Context<'_>,
    query: &Query<'_, R>,
    segment: &Segment<'_>,
    sources: &[comparison::Source<'_>],
    identity: &impl Identity<R>,
    work: &mut usize,
) -> Result<CoreObservation<R>, EngineError> {
    let mut ends = (*query.ends).clone();
    for constraints in &mut ends {
        constraints.sort();
        constraints.dedup();
        if constraints.is_empty() {
            constraints.push(EndConstraint::Unknown);
        }
    }
    let mut out = CoreObservation {
        observation_id: String::new(),
        key: (*query.key).clone(),
        query_id: query.query_id.to_owned(),
        contig: segment.target_name.to_owned(),
        observed_span: segment.target_interval,
        aligned_query_interval: segment.query_interval,
        query_axis_length: query.query_length,
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
        excluded_query_intervals: query.masks.map(<[_]>::to_vec),
        skips: Vec::new(),
        flank_qualification: Vec::new(),
        observed_junction_chain: None,
        comparisons: Vec::new(),
        source_ordinal: segment.source_ordinal,
        clipping: segment.clipping,
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
                .flat_map(|b| unmasked(b.query, query.masks.unwrap_or(&[])))
                .collect(),
        );
        if segment
            .stored_support
            .is_some_and(|s| intervals.as_slice() != s)
        {
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
    out.observation_id = identity.observation_id(input, query, segment, &out)?;
    Ok(out)
}

fn contains(
    a: &CoreFamily,
    b: &CoreFamily,
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

pub(super) fn run<R: Ord + Clone>(
    input: &Context<'_>,
    queries: &[Query<'_, R>],
    unassessed: &[R],
    identity: &impl Identity<R>,
    initial_work: usize,
    check: &mut dyn FnMut() -> Result<(), EngineError>,
) -> Result<CoreSummary<R>, EngineError> {
    let sources = comparison::sources(input)?;
    let mut work = initial_work;
    let mut observations = BTreeMap::<String, CoreObservation<R>>::new();
    let mut families = BTreeMap::<String, CoreFamily>::new();
    let mut record_paths = BTreeMap::<R, (Vec<Vec<String>>, bool, bool)>::new();
    for record in unassessed {
        record_paths.insert(record.clone(), (Vec::new(), false, false));
    }
    let mut queries = queries.iter().collect::<Vec<_>>();
    queries.sort_by(|a, b| (a.key, a.query_id).cmp(&(b.key, b.query_id)));
    for query in queries {
        check()?;
        let mut nodes = BTreeMap::new();
        for segment in &query.segments {
            check()?;
            let observation = observe(input, query, segment, &sources, identity, &mut work)?;
            if let Some(chain) = &observation.observed_junction_chain {
                families.entry(chain.clone()).or_insert_with(|| CoreFamily {
                    chain_id: chain.clone(),
                    contig: observation.contig.clone(),
                    skips: observation.skips.clone(),
                    supporting_keys: 0,
                    definite_keys: 0,
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
                        .flat_map(|b| unmasked(b.query, query.masks.unwrap_or(&[])))
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
                    segment.target_name.to_owned(),
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
            .entry((*query.key).clone())
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
            family.supporting_keys += 1;
            family.definite_keys += usize::from(definite.contains(chain));
            used_families.insert(chain.clone());
        }
        records.push(CoreSupport {
            key: record,
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
    Ok(CoreSummary {
        context_sha256: identity.context_id()?,
        reference: (*input.reference).clone(),
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
        reference_gene: (*input.reference_gene).clone(),
        gene_resolutions: comparison::resolutions(&sources),
        observations: observations.into_values().collect(),
        supports: records,
        families: families
            .into_iter()
            .filter(|(id, _)| used_families.contains(id))
            .map(|(_, f)| f)
            .collect(),
        comparison_work: work,
    })
}
