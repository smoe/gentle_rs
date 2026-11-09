//! Pure V2 transcript-model validation and UTR reach over admitted SAM hypotheses.
//!
//! No file access, hash authentication, alignment search or poly(A) inference.
//! Child modules own explicit saved-input coverage, evidence and pileup views.
//! Callers must authenticate saved artifacts before using this
//! experimental Rust API. Existing SAM records and standalone Python wires do not change.

use super::*;

/// Saved-projection denominators and conditional coverage, without source rereads.
pub mod coverage;
/// Separate unavailable basecall channel and explicitly configured tail imports.
pub mod evidence;
/// Explicit original-query basecalls and literal terminal runs; no implicit FASTQ join.
pub mod oriented;
/// Saved per-unit evidence, explicit view selection and renderer projection.
pub mod pileup;
mod projection;

#[cfg(test)]
mod tests;

/// Existing standalone model wire admitted by this first Rust boundary.
pub const MODEL_SCHEMA: &str = "gentle.rna_read_end_transcript_models.v2";

/// One shared eligibility rule for UTR, coverage and pileup consumers.
pub const ACCEPTED_RELATIONS: [ModelRelation; 3] = [
    ModelRelation::ExactFullChain,
    ModelRelation::ExactPartialChain,
    ModelRelation::WithinExon,
];

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
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
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
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
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
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "status", rename_all = "snake_case", deny_unknown_fields)]
pub enum Metric {
    Observed { value_bp: usize },
    RightCensored { lower_bound_bp: usize },
    Missing { reason: Reason },
}

/// How explicit stop intervals meet the terminal CDS in transcript coordinates.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum StopEvidence {
    Unknown,
    Partial,
    Included,
    Separate,
    Overlapping,
}

/// Derived transcript coordinates only; no raw CDS extension or stop inference.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ModelBoundary {
    pub transcript_id: String,
    pub transcript_length: usize,
    pub cds: Interval,
    pub utr3_start: Option<usize>,
    pub stop_evidence: StopEvidence,
}

/// Identity remains a declared SAM query unit, not a verified original read/molecule.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct UnitProjection {
    pub unit: ObservedSamQueryUnitKey,
    pub read_length: Metric,
    pub utr5: Metric,
    pub utr3: Metric,
    pub candidate_transcript_ids: Vec<String>,
    /// Complete-model identities, empty if any retained candidate lacks a model.
    pub model_gene_ids: Vec<String>,
    /// Saved annotation identities, including candidates without a CDS model.
    pub annotation_gene_ids: Vec<String>,
    pub hypotheses: Vec<HypothesisProjection>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ClipKind {
    None,
    Soft,
    Hard,
    Mixed,
    Unknown,
}

/// Additive evidence only; neither length nor kind reinterprets censoring.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EndClip {
    pub unaligned_query_bases: Option<usize>,
    pub kind: ClipKind,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CandidateRelation {
    pub observation_id: String,
    pub relation: ModelRelation,
    pub gene_ids: Vec<String>,
    pub model_strand: char,
}

/// One retained hypothesis x candidate. No consumer needs to redo its geometry.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CandidateProjection {
    pub transcript_id: String,
    pub annotation_gene_ids: Vec<String>,
    pub model_gene_ids: Vec<String>,
    pub relations: Vec<CandidateRelation>,
    pub transcript_intervals: Vec<Interval>,
    pub utr5: Metric,
    pub utr3: Metric,
    pub clip5: EndClip,
    pub clip3: EndClip,
    pub sense_to_annotation: Option<bool>,
    pub model_available: bool,
    pub projection_available: bool,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HypothesisProjection {
    pub hypothesis_index: usize,
    pub observation_ids: Vec<String>,
    pub candidates: Vec<CandidateProjection>,
    pub unavailable: Option<Reason>,
}

/// Rust-level result, not a registered operation or frozen report wire.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
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

struct Meter<'a> {
    limits: Limits,
    work: usize,
    control: Option<&'a mut dyn FnMut() -> Result<(), EngineError>>,
}

impl Meter<'_> {
    fn charge(&mut self, n: usize) -> Result<(), EngineError> {
        let previous = self.work;
        self.work = self
            .work
            .checked_add(n)
            .filter(|n| *n <= self.limits.max_work)
            .ok_or_else(|| fail("work_budget"))?;
        if (previous == 0 || previous / 1024 != self.work / 1024)
            && let Some(check) = &mut self.control
        {
            check()?;
        }
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

struct OnePath {
    intervals: Vec<Interval>,
    metrics: [Metric; 2],
    clips: [EndClip; 2],
    sense: bool,
    available: bool,
}

fn one_path(
    path: &[&ObservedSamObservation],
    model: &Prepared<'_>,
    meter: &mut Meter,
) -> Result<Option<OnePath>, EngineError> {
    let Some((ends, same)) = projection::query_ends(path, model.model.strand, meter)? else {
        return Ok(None);
    };
    let mut ranges = Vec::new();
    for observation in path {
        meter.charge(1)?;
        if observation.contig != model.model.contig {
            return Ok(None);
        }
        for block in &observation.blocks {
            meter.charge(1)?;
            if block.target.0 >= block.target.1 {
                return Err(fail("observation_geometry"));
            }
            if block.query.1 - block.query.0 != block.target.1 - block.target.0 {
                return Ok(None);
            }
            if ranges.len() >= meter.limits.max_intervals {
                return Err(fail("interval_budget"));
            }
            ranges.push(block.target);
        }
    }
    let projected = model.project(&ranges, meter)?;
    let available = projected.is_some();
    let projected = projected.unwrap_or_default();
    let metrics = if available {
        [
            reach(
                &projected,
                Some(model.boundary.cds.0),
                true,
                ends[0].unaligned_query_bases.is_some_and(|n| n > 0),
            ),
            reach(
                &projected,
                model.boundary.utr3_start,
                false,
                ends[1].unaligned_query_bases.is_some_and(|n| n > 0),
            ),
        ]
    } else {
        [
            missing(Reason::ModelProjectionIncompatible),
            missing(Reason::ModelProjectionIncompatible),
        ]
    };
    Ok(Some(OnePath {
        intervals: projected,
        metrics,
        clips: ends,
        sense: same,
        available,
    }))
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
        if ACCEPTED_RELATIONS.contains(&model.relation) {
            if !identity(&model.transcript_id) {
                return Err(fail("candidate_identity"));
            }
            ids.insert(model.transcript_id.as_str());
        }
    }
    Ok(ids)
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
    annotation: &crate::engine::PrimerPayloadAnnotation,
    limits: Limits,
) -> Result<Projection, EngineError> {
    project_read_ends_with_control(document, analysis, annotation, limits, &mut || Ok(()))
}

/// Same pure calculation with a caller-owned cooperative observation hook.
pub fn project_read_ends_with_control(
    document: &ModelDocument,
    analysis: &ObservedSamAnalysis,
    annotation: &crate::engine::PrimerPayloadAnnotation,
    limits: Limits,
    control: &mut dyn FnMut() -> Result<(), EngineError>,
) -> Result<Projection, EngineError> {
    control()?;
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
    let mut meter = Meter {
        limits,
        work: 0,
        control: Some(control),
    };
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
    if annotation.reference_id != bindings[0].reference_id
        || annotation.assembly != bindings[0].assembly
        || annotation.reference_fasta_sha256 != bindings[0].reference_fasta_sha256
    {
        return Err(fail("model_annotation_binding"));
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
    let normalized = projection::annotation_index(document, annotation, &mut meter)?;
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
        units.push(projection::unit(
            unit,
            &observations,
            &models,
            &normalized,
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
