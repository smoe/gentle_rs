//! Exact gene-ID resolution and supported-block comparison against existing annotation records.

use super::*;
use crate::engine::{PRIMER_PAYLOAD_ANNOTATION_SCHEMA, PrimerPayloadTranscript};

pub(super) struct Source<'a> {
    id: &'a str,
    gene: Option<&'a str>,
    status: GeneResolution,
    models: Vec<&'a PrimerPayloadTranscript>,
}

pub(super) fn resolutions(sources: &[Source<'_>]) -> Vec<SourceComparison> {
    sources
        .iter()
        .map(|s| SourceComparison {
            source_id: s.id.into(),
            gene_id: s.gene.map(str::to_owned),
            resolution: s.status,
            junction_relation: None,
            junctions: None,
            models: Vec::new(),
        })
        .collect()
}

pub(super) fn sources<'a>(input: &core::Context<'a>) -> Result<Vec<Source<'a>>, EngineError> {
    let mut selectors = BTreeMap::new();
    for gene in std::iter::once(input.reference_gene).chain(input.comparison_genes) {
        if gene.source_id.is_empty()
            || gene.gene_id.is_empty()
            || selectors
                .insert(gene.source_id.as_str(), gene.gene_id.as_str())
                .is_some()
        {
            return Err(bad("observed_structure_gene_selector"));
        }
    }
    let mut sources = BTreeMap::new();
    let mut total_rows = 0usize;
    let mut geometry_rows = 0usize;
    for source in input.annotation_sources {
        let a = &source.annotation;
        total_rows = total_rows
            .checked_add(a.transcripts.len())
            .ok_or_else(|| bad("observed_structure_annotation_budget"))?;
        if source.source_id.is_empty()
            || !is_hash(&source.annotation_sha256)
            || sources.insert(source.source_id.as_str(), source).is_some()
            || a.schema != PRIMER_PAYLOAD_ANNOTATION_SCHEMA
            || a.assembly != input.reference.assembly
            || a.reference_fasta_sha256 != input.reference.fasta_sha256
            || a.coordinate_convention != "zero_based_half_open"
            || a.reference_id.is_empty()
            || a.release.is_empty()
            || total_rows > input.budgets.max_annotation_rows
        {
            return Err(bad("observed_structure_annotation_binding"));
        }
        let mut ids = BTreeSet::new();
        let mut lengths = BTreeMap::new();
        for t in &a.transcripts {
            geometry_rows = geometry_rows
                .checked_add(t.genomic_exons.len())
                .ok_or_else(|| bad("observed_structure_annotation_budget"))?;
            if !ids.insert(&t.transcript_id)
                || t.transcript_id.is_empty()
                || t.gene_ids.is_empty()
                || t.gene_ids.iter().any(String::is_empty)
                || !t.gene_ids.windows(2).all(|w| w[0] < w[1])
                || t.target_name.is_empty()
                || t.target_length == 0
                || lengths
                    .insert(&t.target_name, t.target_length)
                    .is_some_and(|n| n != t.target_length)
                || t.genomic_contig.as_ref() != Some(&t.target_name)
                || !matches!(t.genomic_strand, Some('+' | '-'))
                || t.genomic_exons.is_empty()
                || t.genomic_exons
                    .iter()
                    .any(|&(s, e)| s >= e || e > t.target_length)
                || !t.genomic_exons.windows(2).all(|w| w[0].1 <= w[1].0)
                || geometry_rows > input.budgets.max_operations
            {
                return Err(bad("observed_structure_annotation_geometry"));
            }
        }
    }
    let mut ids = sources.keys().copied().collect::<BTreeSet<_>>();
    ids.extend(selectors.keys().copied());
    let mut result = Vec::new();
    for id in ids {
        let gene = selectors.get(id).copied();
        let models = sources
            .get(id)
            .into_iter()
            .flat_map(|s| &s.annotation.transcripts)
            .filter(|t| gene.is_some_and(|g| t.gene_ids.iter().any(|id| id == g)))
            .collect::<Vec<_>>();
        let loci = models
            .iter()
            .map(|t| (&t.target_name, t.genomic_strand))
            .collect::<BTreeSet<_>>();
        let status = if !sources.contains_key(id) || gene.is_none() {
            GeneResolution::UnmappedSource
        } else if models.is_empty() {
            GeneResolution::MissingGene
        } else if loci.len() != 1 {
            GeneResolution::AmbiguousLocus
        } else {
            GeneResolution::Resolved
        };
        result.push(Source {
            id,
            gene,
            status,
            models,
        });
    }
    Ok(result)
}

fn support<R>(
    observation: &core::CoreObservation<R>,
    range: Interval,
    work: &mut usize,
    limit: usize,
) -> Result<Option<usize>, EngineError> {
    let Some(masks) = &observation.excluded_query_intervals else {
        return Ok(None);
    };
    let mut count = 0;
    for block in &observation.blocks {
        dag::charge(work, limit)?;
        if let Some(query) = geometry::project(block, range, observation.query_to_reference_strand)
        {
            dag::charge_many(work, limit, masks.len())?;
            count += unmasked(query, masks)
                .iter()
                .map(|(a, z)| z - a)
                .sum::<usize>();
        }
    }
    Ok(Some(count))
}

pub(super) fn qualify<R>(
    observation: &core::CoreObservation<R>,
    (a, z): Interval,
    flank: usize,
    work: &mut usize,
    limit: usize,
) -> Result<Qualification, EngineError> {
    if observation.excluded_query_intervals.is_none() {
        return Ok(Qualification::Unassessed);
    }
    let (Some(left), Some(right)) = (a.checked_sub(flank), z.checked_add(flank)) else {
        return Ok(Qualification::InsufficientSupport);
    };
    Ok(
        if support(observation, (left, a), work, limit)? == Some(flank)
            && support(observation, (z, right), work, limit)? == Some(flank)
        {
            Qualification::Qualified
        } else {
            Qualification::InsufficientSupport
        },
    )
}

pub(super) fn compare<R>(
    observation: &core::CoreObservation<R>,
    sources: &[Source<'_>],
    input: &core::Context<'_>,
    work: &mut usize,
) -> Result<Vec<SourceComparison>, EngineError> {
    let mut result = Vec::new();
    let mut model_count = 0;
    for source in sources {
        let mut out = SourceComparison {
            source_id: source.id.into(),
            gene_id: source.gene.map(str::to_owned),
            resolution: source.status,
            junction_relation: None,
            junctions: None,
            models: Vec::new(),
        };
        if source.status != GeneResolution::Resolved || observation.detail != DetailStatus::Detailed
        {
            result.push(out);
            continue;
        }
        let mut known = BTreeMap::<Interval, BTreeSet<String>>::new();
        let mut has_run = false;
        for t in &source.models {
            model_count += 1;
            if model_count > input.budgets.max_model_comparisons_per_family {
                return Err(bad("observed_structure_model_budget"));
            }
            dag::charge(work, input.budgets.max_total_work)?;
            let introns = if t.target_name == observation.contig {
                t.genomic_exons
                    .windows(2)
                    .map(|w| (w[0].1, w[1].0))
                    .collect::<Vec<_>>()
            } else {
                Vec::new()
            };
            for intron in &introns {
                dag::charge(work, input.budgets.max_total_work)?;
                known
                    .entry(*intron)
                    .or_default()
                    .insert(t.transcript_id.clone());
            }
            let mut consecutive = false;
            let mut adjacent_coverage = false;
            if !observation.skips.is_empty() && introns.len() >= observation.skips.len() {
                for (start, window) in introns.windows(observation.skips.len()).enumerate() {
                    dag::charge_many(work, input.budgets.max_total_work, observation.skips.len())?;
                    if window == observation.skips {
                        consecutive = true;
                        dag::charge_many(
                            work,
                            input.budgets.max_total_work,
                            observation.blocks.len(),
                        )?;
                        let (a, z) = (
                            t.genomic_exons[start].0,
                            t.genomic_exons[start + observation.skips.len()].1,
                        );
                        adjacent_coverage |= observation
                            .blocks
                            .iter()
                            .all(|b| a <= b.target.0 && b.target.1 <= z);
                    }
                }
            }
            has_run |= consecutive;
            let mut all_exonic =
                !observation.blocks.is_empty() && t.target_name == observation.contig;
            let (mut exonic, mut intronic) = (0usize, 0usize);
            for b in &observation.blocks {
                let mut inside = false;
                if t.target_name == observation.contig {
                    for &(a, z) in &t.genomic_exons {
                        dag::charge(work, input.budgets.max_total_work)?;
                        inside |= a <= b.target.0 && b.target.1 <= z;
                        exonic += z.min(b.target.1).saturating_sub(a.max(b.target.0));
                    }
                    for &(a, z) in &introns {
                        dag::charge(work, input.budgets.max_total_work)?;
                        intronic += z.min(b.target.1).saturating_sub(a.max(b.target.0));
                    }
                }
                all_exonic &= inside;
            }
            let relation = if t.target_name != observation.contig {
                ModelRelation::OutsideModel
            } else if observation.skips.is_empty() {
                let within_one = t.genomic_exons.iter().any(|&(a, z)| {
                    observation
                        .blocks
                        .iter()
                        .all(|b| a <= b.target.0 && b.target.1 <= z)
                });
                if within_one && !observation.blocks.is_empty() {
                    ModelRelation::WithinExon
                } else if exonic > 0 && intronic > 0 {
                    ModelRelation::ExonIntronBoundary
                } else if intronic > 0 && exonic == 0 {
                    ModelRelation::Intronic
                } else {
                    ModelRelation::IncompatibleCoverage
                }
            } else if consecutive && all_exonic && adjacent_coverage {
                if observation.skips.len() == introns.len() {
                    ModelRelation::ExactFullChain
                } else {
                    ModelRelation::ExactPartialChain
                }
            } else if consecutive {
                ModelRelation::TerminalCoverageDivergent
            } else {
                ModelRelation::IncompatibleCoverage
            };
            let mut retention = observation
                .excluded_query_intervals
                .as_ref()
                .map(|_| 0usize);
            for &intron in &introns {
                if qualify(
                    observation,
                    intron,
                    input.required_flank_bases,
                    work,
                    input.budgets.max_total_work,
                )? == Qualification::Qualified
                {
                    if let (Some(total), Some(bases)) = (
                        &mut retention,
                        support(observation, intron, work, input.budgets.max_total_work)?,
                    ) {
                        *total += bases;
                    }
                }
            }
            out.models.push(TranscriptMatch {
                transcript_id: t.transcript_id.clone(),
                gene_ids: t.gene_ids.clone(),
                model_strand: t.genomic_strand.expect("validated annotation"),
                relation,
                retained_intron_supported_bases: retention,
            });
        }
        if !observation.skips.is_empty()
            && source
                .models
                .iter()
                .any(|t| t.target_name == observation.contig)
        {
            out.junctions = Some(
                observation
                    .skips
                    .iter()
                    .map(|skip| JunctionMatch {
                        skipped_interval: *skip,
                        annotated_in_models: known
                            .get(skip)
                            .map(|ids| ids.iter().cloned().collect())
                            .unwrap_or_default(),
                    })
                    .collect(),
            );
            out.junction_relation = Some(
                if observation.skips.iter().any(|s| !known.contains_key(s)) {
                    JunctionRelation::ContainsUnannotatedJunction
                } else if has_run {
                    JunctionRelation::ConsecutiveModelRun
                } else {
                    JunctionRelation::AllJunctionsAnnotatedNoSingleModelContainsRun
                },
            );
        }
        out.models
            .sort_by(|a, b| a.transcript_id.cmp(&b.transcript_id));
        result.push(out);
    }
    Ok(result)
}
