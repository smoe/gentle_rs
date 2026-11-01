//! Portable SAM projections of the shared kernel, with optimal-path and unit-qualified support.

use super::*;

#[cfg(test)]
#[path = "sam_report_tests.rs"]
mod tests;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ObservedSamAnnotationAvailability {
    Available,
    MissingGene,
    Refused,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedSamSourceComparison {
    pub source_id: String,
    pub availability: ObservedSamAnnotationAvailability,
    pub comparison: Option<SourceComparison>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedSamClipping {
    pub hard: [usize; 2],
    pub soft: [usize; 2],
    pub padding_operations: usize,
}

impl From<core::Clipping> for ObservedSamClipping {
    fn from(value: core::Clipping) -> Self {
        Self {
            hard: value.hard,
            soft: value.soft,
            padding_operations: value.padding_operations,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedSamObservation {
    pub observation_id: String,
    pub unit: ObservedSamQueryUnitKey,
    pub source_ordinal: u64,
    pub on_optimal_path: bool,
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
    pub comparisons: Vec<ObservedSamSourceComparison>,
    pub clipping: Option<ObservedSamClipping>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedSamRowAssessment {
    pub ordinal: u64,
    pub qname: String,
    pub read_group: Option<String>,
    pub flag: u16,
    pub target: Option<String>,
    pub position_1based: usize,
    pub reported_mapq: Option<u8>,
    pub cigar: Option<Vec<(char, usize)>>,
    pub sequence_length: Option<usize>,
    pub quality_length: Option<usize>,
    pub axis_length: Option<usize>,
    pub clipping: Option<ObservedSamClipping>,
    pub unavailable: Option<ObservedSamUnavailableReason>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedSamUnitAssessment {
    pub unit: ObservedSamQueryUnitKey,
    pub unclipped_axis_length: Option<usize>,
    pub unavailable: Vec<ObservedSamUnavailableReason>,
    pub rows: Vec<ObservedSamRowAssessment>,
    pub search_declared_complete: bool,
    pub views_declared_complete: bool,
    pub all_query_views_assessed: bool,
    pub hypotheses: Vec<Vec<String>>,
    pub reconciliation: Reconciliation,
    pub topologies: Vec<Topology>,
    pub definite_chains: Vec<String>,
    pub possible_chains: Vec<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ObservedSamQualificationGroup {
    QualifiedOnly,
    InsufficientOnly,
    UnassessedOnly,
    Mixed,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ObservedSamMappingCertainty {
    Definite,
    PossibleOnly,
}

#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedSamQualifiedUnitCounts {
    pub qualified_only: usize,
    pub insufficient_only: usize,
    pub unassessed_only: usize,
    pub mixed: usize,
}

#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedSamUnitSupportCounts {
    pub definite: ObservedSamQualifiedUnitCounts,
    pub possible_only: ObservedSamQualifiedUnitCounts,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedSamCatalogUnit {
    pub unit: ObservedSamQueryUnitKey,
    pub certainty: ObservedSamMappingCertainty,
    pub qualification: ObservedSamQualificationGroup,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedSamChainFamily {
    pub chain_id: String,
    pub contig: String,
    pub skips: Vec<Interval>,
    pub units: Vec<ObservedSamCatalogUnit>,
    pub counts: ObservedSamUnitSupportCounts,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedSamJunction {
    pub junction_id: String,
    pub contig: String,
    pub skipped_interval: Interval,
    pub units: Vec<ObservedSamCatalogUnit>,
    pub counts: ObservedSamUnitSupportCounts,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedSamAnalysis {
    pub context_sha256: String,
    pub reference: ReferenceBinding,
    pub reference_gene: ReferenceGene,
    pub annotation_bindings: Vec<AnnotationBinding>,
    pub gene_resolutions: Vec<ObservedSamSourceComparison>,
    pub observations: Vec<ObservedSamObservation>,
    pub units: Vec<ObservedSamUnitAssessment>,
    pub chain_families: Vec<ObservedSamChainFamily>,
    pub junctions: Vec<ObservedSamJunction>,
    pub comparison_work: usize,
    pub projection_work: usize,
}

fn comparisons(
    rows: &[SourceComparison],
    states: &BTreeMap<String, ObservedSamAnnotationAvailability>,
) -> Vec<ObservedSamSourceComparison> {
    rows.iter()
        .map(|row| {
            let availability = states
                .get(&row.source_id)
                .copied()
                .unwrap_or(ObservedSamAnnotationAvailability::Available);
            ObservedSamSourceComparison {
                source_id: row.source_id.clone(),
                availability,
                comparison: (availability == ObservedSamAnnotationAvailability::Available)
                    .then(|| row.clone()),
            }
        })
        .collect()
}

fn qualification(values: impl Iterator<Item = Qualification>) -> ObservedSamQualificationGroup {
    let mut flags = [false; 3];
    for value in values {
        flags[match value {
            Qualification::Qualified => 0,
            Qualification::InsufficientSupport => 1,
            Qualification::Unassessed => 2,
        }] = true;
    }
    match flags {
        [true, false, false] => ObservedSamQualificationGroup::QualifiedOnly,
        [false, true, false] => ObservedSamQualificationGroup::InsufficientOnly,
        [false, false, true] | [false, false, false] => {
            ObservedSamQualificationGroup::UnassessedOnly
        }
        _ => ObservedSamQualificationGroup::Mixed,
    }
}

fn counts(units: &[ObservedSamCatalogUnit]) -> ObservedSamUnitSupportCounts {
    let mut out = ObservedSamUnitSupportCounts::default();
    for unit in units {
        let group = match unit.certainty {
            ObservedSamMappingCertainty::Definite => &mut out.definite,
            ObservedSamMappingCertainty::PossibleOnly => &mut out.possible_only,
        };
        match unit.qualification {
            ObservedSamQualificationGroup::QualifiedOnly => group.qualified_only += 1,
            ObservedSamQualificationGroup::InsufficientOnly => group.insufficient_only += 1,
            ObservedSamQualificationGroup::UnassessedOnly => group.unassessed_only += 1,
            ObservedSamQualificationGroup::Mixed => group.mixed += 1,
        }
    }
    out
}

fn core_admission(
    request: &sam_intake::Request,
    prepared: &sam_intake::Prepared,
) -> Result<(), EngineError> {
    // Bound repeated model projections before the kernel allocates comparison rows.
    // Engineering accounting, not a measured RSS or throughput guarantee.
    let selectors = std::iter::once(&request.reference_gene)
        .chain(&request.comparison_genes)
        .map(|g| (&g.source_id, &g.gene_id))
        .collect::<BTreeMap<_, _>>();
    let mut per_row = 8192usize;
    for source in &request.annotation_sources {
        let Some(gene) = selectors.get(&source.source_id) else {
            continue;
        };
        for model in &source.annotation.transcripts {
            if model.gene_ids.contains(gene) {
                per_row = model
                    .genomic_exons
                    .len()
                    .checked_mul(512)
                    .and_then(|n| n.checked_add(32768))
                    .and_then(|n| n.checked_add(per_row))
                    .ok_or_else(|| bad("observed_sam_projection_memory_budget"))?;
            }
        }
    }
    let mut charge = prepared.receipt.retained_peak_bytes;
    for unit in &prepared.input.units {
        for row in &unit.rows {
            let ops = row.cigar.as_deref().unwrap_or(&[]);
            let skips = ops.iter().filter(|&&(c, _)| c == 'N').count();
            let cost = skips
                .checked_mul(4 * unit.unit_id.len() + 2048)
                .and_then(|n| ops.len().checked_mul(256).and_then(|c| n.checked_add(c)))
                .and_then(|n| n.checked_add(16 * unit.unit_id.len()))
                .and_then(|n| {
                    request
                        .kernel_budgets
                        .max_alternatives_per_record
                        .checked_mul(256)
                        .and_then(|c| n.checked_add(c))
                })
                .and_then(|n| n.checked_add(per_row))
                .ok_or_else(|| bad("observed_sam_projection_memory_budget"))?;
            charge = charge
                .checked_add(cost)
                .filter(|n| *n <= request.limits.max_retained_bytes)
                .ok_or_else(|| bad("observed_sam_projection_memory_budget"))?;
        }
    }
    Ok(())
}

pub(crate) fn summarize_bound(
    request: &sam_intake::Request,
    states: &BTreeMap<String, ObservedSamAnnotationAvailability>,
    max_projection_work: usize,
    check: &mut dyn FnMut() -> Result<(), EngineError>,
) -> Result<(ObservedSamAnalysis, sam_intake::Receipt), EngineError> {
    let prepared = sam_intake::prepare(request, check)?;
    core_admission(request, &prepared)?;
    let analyzed = sam::analyze(&prepared.input, check)?;
    let data = project(&prepared, analyzed, states, max_projection_work, check)?;
    Ok((data, prepared.receipt))
}

fn project(
    prepared: &sam_intake::Prepared,
    analysis: sam::Analysis,
    states: &BTreeMap<String, ObservedSamAnnotationAvailability>,
    max_work: usize,
    check: &mut dyn FnMut() -> Result<(), EngineError>,
) -> Result<ObservedSamAnalysis, EngineError> {
    let summary = analysis.summary;
    let mut work = 0usize;
    let by_id = summary
        .observations
        .iter()
        .map(|o| (o.observation_id.as_str(), o))
        .collect::<BTreeMap<_, _>>();
    let supports = summary
        .supports
        .iter()
        .map(|s| (&s.key, s))
        .collect::<BTreeMap<_, _>>();
    let optimal = summary
        .supports
        .iter()
        .flat_map(|s| s.hypotheses.iter().flatten())
        .collect::<BTreeSet<_>>();
    let mut family_units = BTreeMap::<String, Vec<ObservedSamCatalogUnit>>::new();
    let mut junction_units = BTreeMap::<(String, Interval), Vec<ObservedSamCatalogUnit>>::new();
    for support in &summary.supports {
        check()?;
        dag::charge(&mut work, max_work)?;
        let ids = support.hypotheses.iter().flatten().collect::<BTreeSet<_>>();
        let mut occurrences = BTreeMap::<(String, Interval), Vec<Qualification>>::new();
        let mut chains = BTreeMap::<String, Vec<Qualification>>::new();
        for id in &ids {
            check()?;
            let observation = by_id
                .get(id.as_str())
                .ok_or_else(|| bad("observed_sam_optimal_observation_binding"))?;
            dag::charge_many(&mut work, max_work, observation.skips.len() + 1)?;
            if observation.key != support.key
                || observation.skips.len() != observation.flank_qualification.len()
            {
                return Err(bad("observed_sam_optimal_observation_binding"));
            }
            if let Some(chain) = &observation.observed_junction_chain {
                chains
                    .entry(chain.clone())
                    .or_default()
                    .extend(&observation.flank_qualification);
            }
            for (&skip, &value) in observation
                .skips
                .iter()
                .zip(&observation.flank_qualification)
            {
                occurrences
                    .entry((observation.contig.clone(), skip))
                    .or_default()
                    .push(value);
            }
        }
        for (chain, values) in chains {
            if !support.possible_chains.contains(&chain) {
                return Err(bad("observed_sam_family_binding"));
            }
            family_units
                .entry(chain.clone())
                .or_default()
                .push(ObservedSamCatalogUnit {
                    unit: support.key.clone(),
                    certainty: if support.definite_chains.contains(&chain) {
                        ObservedSamMappingCertainty::Definite
                    } else {
                        ObservedSamMappingCertainty::PossibleOnly
                    },
                    qualification: qualification(values.into_iter()),
                });
        }
        // A junction can be present in every optimal path even when the complete chains differ.
        let mut definite: Option<BTreeSet<(String, Interval)>> = None;
        for path in &support.hypotheses {
            check()?;
            dag::charge_many(&mut work, max_work, path.len() + 1)?;
            let mut keys = BTreeSet::new();
            for id in path {
                let observation = by_id
                    .get(id.as_str())
                    .ok_or_else(|| bad("observed_sam_optimal_observation_binding"))?;
                dag::charge_many(&mut work, max_work, observation.skips.len() + 1)?;
                keys.extend(
                    observation
                        .skips
                        .iter()
                        .map(|&skip| (observation.contig.clone(), skip)),
                );
            }
            definite = Some(match definite {
                None => keys,
                Some(old) => {
                    dag::charge_many(&mut work, max_work, old.len() + keys.len())?;
                    old.intersection(&keys).cloned().collect()
                }
            });
        }
        let definite = if support.search_declared_complete && support.all_query_views_assessed {
            definite.unwrap_or_default()
        } else {
            BTreeSet::new()
        };
        for (key, values) in occurrences {
            let certainty = if definite.contains(&key) {
                ObservedSamMappingCertainty::Definite
            } else {
                ObservedSamMappingCertainty::PossibleOnly
            };
            junction_units
                .entry(key)
                .or_default()
                .push(ObservedSamCatalogUnit {
                    unit: support.key.clone(),
                    certainty,
                    qualification: qualification(values.into_iter()),
                });
        }
    }
    let mut families = Vec::new();
    for family in &summary.families {
        check()?;
        let units = family_units.remove(&family.chain_id).unwrap_or_default();
        if units.len() != family.supporting_keys
            || units
                .iter()
                .filter(|u| u.certainty == ObservedSamMappingCertainty::Definite)
                .count()
                != family.definite_keys
        {
            return Err(bad("observed_sam_family_count_binding"));
        }
        families.push(ObservedSamChainFamily {
            chain_id: family.chain_id.clone(),
            contig: family.contig.clone(),
            skips: family.skips.clone(),
            counts: counts(&units),
            units,
        });
    }
    if !family_units.is_empty() {
        return Err(bad("observed_sam_family_binding"));
    }
    let mut junctions = Vec::new();
    for ((contig, skip), units) in junction_units {
        check()?;
        junctions.push(ObservedSamJunction {
            junction_id: digest(&(
                "gentle.observed_splice_junction.v1",
                &summary.reference,
                &contig,
                skip,
            ))?,
            contig,
            skipped_interval: skip,
            counts: counts(&units),
            units,
        });
    }
    let mut observations = Vec::new();
    for o in &summary.observations {
        check()?;
        dag::charge_many(&mut work, max_work, o.skips.len() + o.comparisons.len() + 1)?;
        observations.push(ObservedSamObservation {
            observation_id: o.observation_id.clone(),
            unit: o.key.clone(),
            source_ordinal: o
                .source_ordinal
                .ok_or_else(|| bad("observed_sam_row_binding"))?,
            on_optimal_path: optimal.contains(&o.observation_id),
            contig: o.contig.clone(),
            observed_span: o.observed_span,
            aligned_query_interval: o.aligned_query_interval,
            query_axis_length: o.query_axis_length,
            alignment_clipped: o.alignment_clipped,
            query_to_reference_strand: o.query_to_reference_strand,
            reported_mapq: o.reported_mapq,
            ends: o.ends.clone(),
            detail: o.detail,
            topology: o.topology,
            blocks: o.blocks.clone(),
            excluded_query_intervals: o.excluded_query_intervals.clone(),
            skips: o.skips.clone(),
            flank_qualification: o.flank_qualification.clone(),
            observed_junction_chain: o.observed_junction_chain.clone(),
            comparisons: comparisons(&o.comparisons, states),
            clipping: o.clipping.map(Into::into),
        });
    }
    let witnesses = prepared
        .receipt
        .query_unit_manifest
        .units
        .iter()
        .flat_map(|u| &u.rows)
        .map(|r| (r.ordinal, r))
        .collect::<BTreeMap<_, _>>();
    let mut units = Vec::new();
    for (input, diagnostic) in prepared.input.units.iter().zip(&analysis.units) {
        check()?;
        let support = supports
            .get(&diagnostic.unit)
            .ok_or_else(|| bad("observed_sam_unit_binding"))?;
        let by_ordinal = diagnostic
            .rows
            .iter()
            .map(|r| (r.ordinal, r))
            .collect::<BTreeMap<_, _>>();
        let mut rows = Vec::new();
        for row in &input.rows {
            check()?;
            dag::charge_many(
                &mut work,
                max_work,
                row.cigar.as_ref().map_or(0, Vec::len) + 1,
            )?;
            let d = by_ordinal
                .get(&row.ordinal)
                .ok_or_else(|| bad("observed_sam_row_binding"))?;
            let witness = witnesses
                .get(&row.ordinal)
                .ok_or_else(|| bad("observed_sam_row_binding"))?;
            rows.push(ObservedSamRowAssessment {
                ordinal: row.ordinal,
                qname: row.qname.clone(),
                read_group: witness.read_group.clone(),
                flag: row.flag,
                target: row.target.clone(),
                position_1based: row.position_1based,
                reported_mapq: (row.mapq != 255).then_some(row.mapq),
                cigar: row.cigar.clone(),
                sequence_length: row.sequence_length,
                quality_length: row.quality_length,
                axis_length: d.axis_length,
                clipping: d.clipping.map(Into::into),
                unavailable: d.unavailable,
            });
        }
        rows.sort_by_key(|r| r.ordinal);
        units.push(ObservedSamUnitAssessment {
            unit: diagnostic.unit.clone(),
            unclipped_axis_length: diagnostic.axis_length,
            unavailable: diagnostic.unavailable.iter().copied().collect(),
            rows,
            search_declared_complete: input.search_declared_complete,
            views_declared_complete: input.views_declared_complete,
            all_query_views_assessed: support.all_query_views_assessed,
            hypotheses: support.hypotheses.clone(),
            reconciliation: support.reconciliation,
            topologies: support.topologies.clone(),
            definite_chains: support.definite_chains.clone(),
            possible_chains: support.possible_chains.clone(),
        });
    }
    units.sort_by(|a, b| a.unit.cmp(&b.unit));
    Ok(ObservedSamAnalysis {
        context_sha256: summary.context_sha256,
        reference: summary.reference,
        reference_gene: summary.reference_gene,
        annotation_bindings: summary.annotation_bindings,
        gene_resolutions: comparisons(&summary.gene_resolutions, states),
        observations,
        units,
        chain_families: families,
        junctions,
        comparison_work: summary.comparison_work,
        projection_work: work,
    })
}
