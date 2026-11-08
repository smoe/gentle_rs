//! Pure V2 transcript-model validation and UTR reach over admitted SAM hypotheses.
//!
//! No file access, hash authentication, alignment search, sample selection, coverage
//! or poly(A) inference. Callers must authenticate saved artifacts before using this
//! experimental Rust API. Existing SAM records and standalone Python wires do not change.

use super::*;

#[cfg(test)]
mod tests;

/// Existing standalone model wire admitted by this first Rust boundary.
pub const MODEL_SCHEMA: &str = "gentle.rna_read_end_transcript_models.v2";

/// The existing standalone V2 model, not a new annotation or reference authority.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ModelDocument {
    pub schema: String,
    pub source_id: String,
    pub assembly: String,
    pub reference_fasta_sha256: String,
    pub annotation_sha256: String,
    pub coordinate_convention: String,
    pub model_source: ModelSource,
    pub transcripts: Vec<TranscriptModel>,
}

/// Explicit GTF snapshot provenance; this module never opens or authenticates it.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ModelSource {
    pub format: String,
    pub name: String,
    pub sha256: String,
}

/// Genomic ascending half-open intervals; raw CDS and stops remain separate.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TranscriptModel {
    pub transcript_id: String,
    pub gene_id: String,
    pub contig: String,
    pub strand: char,
    pub genomic_exons_0based_half_open: Vec<Interval>,
    pub genomic_cds_0based_half_open: Vec<Interval>,
    pub genomic_stop_codon_0based_half_open: Vec<Interval>,
}

/// Bounds apply to this projection, not the upstream parser or allocator RSS.
#[derive(Debug, Clone, Copy)]
pub struct Limits {
    pub max_models: usize,
    pub max_intervals: usize,
    pub max_units: usize,
    pub max_observations: usize,
    pub max_work: usize,
}

impl Default for Limits {
    fn default() -> Self {
        Self {
            max_models: 100_000,
            max_intervals: 2_000_000,
            max_units: 2_000_000,
            max_observations: 2_000_000,
            max_work: 10_000_000,
        }
    }
}

/// Reasons have the same meanings as the standalone projector's detailed table.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Reason {
    AxisLengthUnavailable,
    NoOptimalHypothesis,
    InvalidOptimalHypothesis,
    UnboundOptimalHypothesis,
    NoCompatibleCdsModel,
    CandidateModelUnavailable,
    ModelProjectionIncompatible,
    AmbiguousTranscriptEndGeometry,
    CdsBoundaryUnanchored,
    StopCodonUnavailable,
}

/// Missing has no number; a censored lower bound cannot masquerade as an observation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Metric {
    Observed { value_bp: usize },
    RightCensored { lower_bound_bp: usize },
    Missing { reason: Reason },
}

/// How explicit stop intervals meet the terminal CDS in transcript coordinates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StopEvidence {
    Unknown,
    Partial,
    Included,
    Separate,
    Overlapping,
}

/// Derived transcript coordinates only; no raw CDS extension or stop inference.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ModelBoundary {
    pub transcript_id: String,
    pub transcript_length: usize,
    pub cds: Interval,
    pub utr3_start: Option<usize>,
    pub stop_evidence: StopEvidence,
}

/// Identity remains a declared SAM query unit, not a verified original read/molecule.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UnitProjection {
    pub unit: ObservedSamQueryUnitKey,
    pub read_length: Metric,
    pub utr5: Metric,
    pub utr3: Metric,
    pub candidate_transcript_ids: Vec<String>,
    pub gene_ids: Vec<String>,
}

/// Rust-level result, not a registered operation or frozen report wire.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Projection {
    pub source_id: String,
    pub models: Vec<ModelBoundary>,
    pub units: Vec<UnitProjection>,
    pub work: usize,
}

fn fail(reason: &str) -> EngineError {
    EngineError::invalid_input(format!("rna_read_end_{reason}"))
}

fn token(value: &str) -> bool {
    (1..=128).contains(&value.len())
        && value.as_bytes()[0].is_ascii_alphanumeric()
        && value
            .bytes()
            .all(|c| c.is_ascii_alphanumeric() || b"_.-".contains(&c))
}

fn sha(value: &str) -> bool {
    value.strip_prefix("sha256:").is_some_and(|s| {
        s.len() == 64
            && s.bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    })
}

fn identity(value: &str) -> bool {
    !value.is_empty() && value.len() <= 4096 && !value.chars().any(char::is_control)
}

struct Meter {
    limits: Limits,
    work: usize,
}

impl Meter {
    fn charge(&mut self, n: usize) -> Result<(), EngineError> {
        self.work = self
            .work
            .checked_add(n)
            .filter(|n| *n <= self.limits.max_work)
            .ok_or_else(|| fail("work_budget"))?;
        Ok(())
    }

    fn merge(&mut self, mut ranges: Vec<Interval>) -> Result<Vec<Interval>, EngineError> {
        if ranges.len() > self.limits.max_intervals {
            return Err(fail("interval_budget"));
        }
        let sort_work = ranges
            .len()
            .checked_mul((usize::BITS - ranges.len().leading_zeros()) as usize + 1);
        self.charge(sort_work.ok_or_else(|| fail("overflow"))?)?;
        ranges.sort_unstable();
        let mut merged: Vec<Interval> = Vec::new();
        for (start, end) in ranges {
            if let Some(last) = merged.last_mut()
                && start <= last.1
            {
                last.1 = last.1.max(end);
            } else {
                merged.push((start, end));
            }
        }
        Ok(merged)
    }
}

fn intervals(values: &[Interval], required: bool, meter: &mut Meter) -> Result<(), EngineError> {
    if (required && values.is_empty()) || values.len() > meter.limits.max_intervals {
        return Err(fail("model_intervals"));
    }
    let mut previous = None;
    for &(start, end) in values {
        meter.charge(1)?;
        if start >= end || previous.is_some_and(|last| start <= last) {
            return Err(fail("model_intervals"));
        }
        previous = Some(end);
    }
    Ok(())
}

struct Prepared<'a> {
    model: &'a TranscriptModel,
    exons: Vec<SupportedBlock>,
    boundary: ModelBoundary,
}

impl Prepared<'_> {
    fn project(
        &self,
        ranges: &[Interval],
        meter: &mut Meter,
    ) -> Result<Option<Vec<Interval>>, EngineError> {
        let mut projected = Vec::new();
        for &(start, end) in ranges {
            let mut covered = 0usize;
            for exon in &self.exons {
                meter.charge(1)?;
                if let Some(value) = geometry::project(exon, (start, end), self.model.strand) {
                    covered = covered
                        .checked_add(value.1 - value.0)
                        .ok_or_else(|| fail("overflow"))?;
                    if projected.len() >= meter.limits.max_intervals {
                        return Err(fail("interval_budget"));
                    }
                    projected.push(value);
                }
            }
            if covered != end - start {
                return Ok(None);
            }
        }
        meter.merge(projected).map(Some)
    }
}

fn prepare<'a>(model: &'a TranscriptModel, meter: &mut Meter) -> Result<Prepared<'a>, EngineError> {
    if !token(&model.transcript_id)
        || !token(&model.gene_id)
        || !token(&model.contig)
        || !matches!(model.strand, '+' | '-')
    {
        return Err(fail("model_identity"));
    }
    intervals(&model.genomic_exons_0based_half_open, true, meter)?;
    intervals(&model.genomic_cds_0based_half_open, true, meter)?;
    intervals(&model.genomic_stop_codon_0based_half_open, false, meter)?;
    let mut transcript_length = 0usize;
    let mut exons = Vec::new();
    for position in 0..model.genomic_exons_0based_half_open.len() {
        let index = if model.strand == '+' {
            position
        } else {
            model.genomic_exons_0based_half_open.len() - 1 - position
        };
        let (start, end) = model.genomic_exons_0based_half_open[index];
        let next = transcript_length
            .checked_add(end - start)
            .ok_or_else(|| fail("overflow"))?;
        exons.push(SupportedBlock {
            query: (transcript_length, next),
            target: (start, end),
        });
        transcript_length = next;
    }
    let mut value = Prepared {
        model,
        exons,
        boundary: ModelBoundary {
            transcript_id: model.transcript_id.clone(),
            transcript_length,
            cds: (0, 0),
            utr3_start: None,
            stop_evidence: StopEvidence::Unknown,
        },
    };
    let cds = value
        .project(&model.genomic_cds_0based_half_open, meter)?
        .ok_or_else(|| fail("cds_outside_exons"))?;
    value.boundary.cds = (cds[0].0, cds[cds.len() - 1].1);
    let stop = value
        .project(&model.genomic_stop_codon_0based_half_open, meter)?
        .ok_or_else(|| fail("stop_outside_exons"))?;
    let bases = stop.iter().try_fold(0usize, |n, &(a, z)| {
        n.checked_add(z - a).ok_or_else(|| fail("overflow"))
    })?;
    value.boundary.stop_evidence = match bases {
        0 => StopEvidence::Unknown,
        1 | 2 => StopEvidence::Partial,
        3 => {
            if stop.len() != 1
                || !(value.boundary.cds.0 <= stop[0].0
                    && stop[0].0 <= value.boundary.cds.1
                    && value.boundary.cds.1 <= stop[0].1)
            {
                return Err(fail("stop_geometry"));
            }
            value.boundary.utr3_start = Some(stop[0].1);
            if stop[0].1 == value.boundary.cds.1 {
                StopEvidence::Included
            } else if stop[0].0 == value.boundary.cds.1 {
                StopEvidence::Separate
            } else {
                StopEvidence::Overlapping
            }
        }
        _ => return Err(fail("stop_geometry")),
    };
    Ok(value)
}

fn missing(reason: Reason) -> Metric {
    Metric::Missing { reason }
}

fn reach(ranges: &[Interval], boundary: Option<usize>, five: bool, clipped: bool) -> Metric {
    let Some(boundary) = boundary else {
        return missing(Reason::StopCodonUnavailable);
    };
    let value = ranges.iter().find_map(|&(start, end)| {
        if five && start <= boundary && boundary < end {
            Some(boundary - start)
        } else if !five && start < boundary && boundary <= end {
            Some(end - boundary)
        } else {
            None
        }
    });
    match value {
        None => missing(Reason::CdsBoundaryUnanchored),
        Some(value) if clipped => Metric::RightCensored {
            lower_bound_bp: value,
        },
        Some(value_bp) => Metric::Observed { value_bp },
    }
}

fn one_path(
    path: &[&ObservedSamObservation],
    model: &Prepared<'_>,
    meter: &mut Meter,
) -> Result<Option<[Metric; 2]>, EngineError> {
    let mut ranges = Vec::new();
    let mut query_start = usize::MAX;
    let mut query_end = 0usize;
    let mut axis_orientation = None;
    for observation in path {
        meter.charge(1)?;
        if observation.contig != model.model.contig
            || observation.blocks.is_empty()
            || observation.query_axis_length == 0
            || !matches!(observation.query_to_reference_strand, '+' | '-')
        {
            return Ok(None);
        }
        let current = (
            observation.query_axis_length,
            observation.query_to_reference_strand == model.model.strand,
        );
        if axis_orientation.is_some_and(|previous| previous != current) {
            return Ok(None);
        }
        axis_orientation = Some(current);
        for block in &observation.blocks {
            meter.charge(1)?;
            if block.query.0 >= block.query.1
                || block.target.0 >= block.target.1
                || block.query.1 > current.0
            {
                return Err(fail("observation_geometry"));
            }
            if block.query.1 - block.query.0 != block.target.1 - block.target.0 {
                return Ok(None);
            }
            if ranges.len() >= meter.limits.max_intervals {
                return Err(fail("interval_budget"));
            }
            ranges.push(block.target);
            query_start = query_start.min(block.query.0);
            query_end = query_end.max(block.query.1);
        }
    }
    let Some((axis, same)) = axis_orientation else {
        return Ok(None);
    };
    let Some(projected) = model.project(&ranges, meter)? else {
        return Ok(None);
    };
    Ok(Some([
        reach(
            &projected,
            Some(model.boundary.cds.0),
            true,
            if same {
                query_start > 0
            } else {
                query_end < axis
            },
        ),
        reach(
            &projected,
            model.boundary.utr3_start,
            false,
            if same {
                query_end < axis
            } else {
                query_start > 0
            },
        ),
    ]))
}

fn candidates<'a>(
    observation: &'a ObservedSamObservation,
    source_id: &str,
    meter: &mut Meter,
) -> Result<BTreeSet<&'a str>, EngineError> {
    let mut selected = None;
    for row in &observation.comparisons {
        meter.charge(1)?;
        if row.source_id == source_id {
            if selected.is_some() {
                return Ok(BTreeSet::new());
            }
            selected = Some(row);
        }
    }
    let Some(row) = selected else {
        return Ok(BTreeSet::new());
    };
    let Some(comparison) = &row.comparison else {
        return Ok(BTreeSet::new());
    };
    if row.availability != ObservedSamAnnotationAvailability::Available
        || comparison.resolution != GeneResolution::Resolved
    {
        return Ok(BTreeSet::new());
    }
    let mut ids = BTreeSet::new();
    for model in &comparison.models {
        meter.charge(1)?;
        if matches!(
            model.relation,
            ModelRelation::ExactFullChain
                | ModelRelation::ExactPartialChain
                | ModelRelation::WithinExon
        ) {
            if !identity(&model.transcript_id) {
                return Err(fail("candidate_identity"));
            }
            ids.insert(model.transcript_id.as_str());
        }
    }
    Ok(ids)
}

fn unit_projection(
    unit: &ObservedSamUnitAssessment,
    observations: &BTreeMap<&str, &ObservedSamObservation>,
    models: &BTreeMap<&str, Prepared<'_>>,
    source_id: &str,
    meter: &mut Meter,
) -> Result<UnitProjection, EngineError> {
    let read_length = match unit.unclipped_axis_length {
        Some(0) => return Err(fail("unit_axis")),
        Some(value_bp) => Metric::Observed { value_bp },
        None => missing(Reason::AxisLengthUnavailable),
    };
    let mut out = UnitProjection {
        unit: unit.unit.clone(),
        read_length,
        utr5: missing(Reason::NoOptimalHypothesis),
        utr3: missing(Reason::NoOptimalHypothesis),
        candidate_transcript_ids: Vec::new(),
        gene_ids: Vec::new(),
    };
    let mut paths = Vec::new();
    let mut all_candidates = BTreeSet::new();
    let mut unavailable = None;
    for hypothesis in &unit.hypotheses {
        meter.charge(1)?;
        if hypothesis.is_empty() {
            unavailable = Some(Reason::InvalidOptimalHypothesis);
            break;
        }
        let mut path = Vec::new();
        let mut path_ids = BTreeSet::new();
        let mut common: Option<BTreeSet<&str>> = None;
        for id in hypothesis {
            meter.charge(1)?;
            if !path_ids.insert(id) {
                return Err(fail("hypothesis_identity"));
            }
            let Some(&observation) = observations.get(id.as_str()) else {
                unavailable = Some(Reason::UnboundOptimalHypothesis);
                break;
            };
            if observation.unit != unit.unit
                || !observation.on_optimal_path
                || unit
                    .unclipped_axis_length
                    .is_some_and(|axis| axis != observation.query_axis_length)
            {
                return Err(fail("unit_binding"));
            }
            let current = candidates(observation, source_id, meter)?;
            common = Some(match common {
                None => current,
                Some(previous) => {
                    meter.charge(previous.len())?;
                    previous.intersection(&current).copied().collect()
                }
            });
            path.push(observation);
        }
        if unavailable.is_some() {
            break;
        }
        let common = common.unwrap_or_default();
        if common.is_empty() {
            unavailable = Some(Reason::NoCompatibleCdsModel);
            break;
        }
        meter.charge(common.len())?;
        all_candidates.extend(common.iter().copied());
        paths.push((path, common));
    }
    if let Some(reason) = unavailable {
        out.utr5 = missing(reason);
        out.utr3 = missing(reason);
        return Ok(out);
    }
    out.candidate_transcript_ids = all_candidates.iter().map(|id| (*id).to_owned()).collect();
    meter.charge(all_candidates.len())?;
    if all_candidates.iter().any(|id| !models.contains_key(id)) {
        out.utr5 = missing(Reason::CandidateModelUnavailable);
        out.utr3 = missing(Reason::CandidateModelUnavailable);
        return Ok(out);
    }
    out.gene_ids = all_candidates
        .iter()
        .map(|id| models[id].model.gene_id.clone())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let mut consensus: Option<[Metric; 2]> = None;
    for (path, ids) in paths {
        for id in ids {
            let Some(value) = one_path(&path, &models[id], meter)? else {
                out.utr5 = missing(Reason::ModelProjectionIncompatible);
                out.utr3 = missing(Reason::ModelProjectionIncompatible);
                return Ok(out);
            };
            if let Some(previous) = &mut consensus {
                for index in 0..2 {
                    if previous[index] != value[index] {
                        previous[index] = missing(Reason::AmbiguousTranscriptEndGeometry);
                    }
                }
            } else {
                consensus = Some(value);
            }
        }
    }
    if let Some([five, three]) = consensus {
        out.utr5 = five;
        out.utr3 = three;
    }
    Ok(out)
}

/// Bind a V2 model to admitted SAM analysis and project every declared unit.
///
/// Equality of stored hashes is a binding check, not authentication of source bytes,
/// search completeness, RNA ends or original-read identity. No inputs are modified.
/// V1 remains supported conservatively by the unchanged standalone Python reader;
/// this first Rust slice accepts V2 only, with required explicit stop intervals.
pub fn project_read_ends(
    document: &ModelDocument,
    analysis: &ObservedSamAnalysis,
    limits: Limits,
) -> Result<Projection, EngineError> {
    if limits.max_models == 0
        || limits.max_intervals == 0
        || limits.max_units == 0
        || limits.max_observations == 0
        || limits.max_work == 0
    {
        return Err(fail("limits"));
    }
    if document.schema != MODEL_SCHEMA
        || !token(&document.source_id)
        || document.assembly.is_empty()
        || document.coordinate_convention != "zero_based_half_open"
        || !sha(&document.reference_fasta_sha256)
        || !sha(&document.annotation_sha256)
        || document.model_source.format != "gtf"
        || document.model_source.name.is_empty()
        || !sha(&document.model_source.sha256)
    {
        return Err(fail("model_document"));
    }
    if document.transcripts.is_empty() || document.transcripts.len() > limits.max_models {
        return Err(fail("model_budget"));
    }
    if analysis.units.len() > limits.max_units
        || analysis.observations.len() > limits.max_observations
    {
        return Err(fail("analysis_budget"));
    }
    if analysis.reference.assembly != document.assembly
        || analysis.reference.fasta_sha256 != document.reference_fasta_sha256
    {
        return Err(fail("reference_binding"));
    }
    let mut meter = Meter { limits, work: 0 };
    let mut bindings = Vec::new();
    for binding in &analysis.annotation_bindings {
        meter.charge(1)?;
        if binding.source_id == document.source_id {
            bindings.push(binding);
        }
    }
    if bindings.len() != 1 || bindings[0].annotation_sha256 != document.annotation_sha256 {
        return Err(fail("annotation_binding"));
    }
    let mut models = BTreeMap::new();
    let mut interval_count = 0usize;
    for model in &document.transcripts {
        meter.charge(1)?;
        for values in [
            &model.genomic_exons_0based_half_open,
            &model.genomic_cds_0based_half_open,
            &model.genomic_stop_codon_0based_half_open,
        ] {
            interval_count = interval_count
                .checked_add(values.len())
                .filter(|n| *n <= limits.max_intervals)
                .ok_or_else(|| fail("interval_budget"))?;
        }
        if models
            .insert(model.transcript_id.as_str(), prepare(model, &mut meter)?)
            .is_some()
        {
            return Err(fail("duplicate_model"));
        }
    }
    let mut observations = BTreeMap::new();
    for observation in &analysis.observations {
        meter.charge(1)?;
        if !identity(&observation.observation_id)
            || observations
                .insert(observation.observation_id.as_str(), observation)
                .is_some()
        {
            return Err(fail("observation_identity"));
        }
    }
    let mut units = Vec::new();
    let mut keys = BTreeSet::new();
    for unit in &analysis.units {
        meter.charge(1)?;
        if !identity(&unit.unit.unit_id)
            || !sha(&unit.unit.source_sha256)
            || !keys.insert(&unit.unit)
        {
            return Err(fail("unit_identity"));
        }
        units.push(unit_projection(
            unit,
            &observations,
            &models,
            &document.source_id,
            &mut meter,
        )?);
    }
    Ok(Projection {
        source_id: document.source_id.clone(),
        models: models.values().map(|m| m.boundary.clone()).collect(),
        units,
        work: meter.work,
    })
}
