//! Annotation-checked hypothesis/candidate records and their conservative unit collapse.

use super::*;
use crate::engine::{
    PRIMER_PAYLOAD_ANNOTATION_SCHEMA, PrimerPayloadAnnotation, PrimerPayloadTranscript,
};

fn unknown_clip() -> EndClip {
    EndClip {
        unaligned_query_bases: None,
        kind: ClipKind::Unknown,
    }
}

pub(super) fn query_ends(
    path: &[&ObservedSamObservation],
    strand: char,
    meter: &mut Meter,
) -> Result<Option<([EndClip; 2], bool)>, EngineError> {
    let mut frame = None;
    let (mut start, mut end) = (usize::MAX, 0);
    for observation in path {
        meter.charge(1)?;
        if observation.query_axis_length == 0
            || observation.blocks.is_empty()
            || !matches!(strand, '+' | '-')
            || !matches!(observation.query_to_reference_strand, '+' | '-')
        {
            return Ok(None);
        }
        let current = (
            observation.query_axis_length,
            strand == observation.query_to_reference_strand,
        );
        if frame.is_some_and(|old| old != current) {
            return Ok(None);
        }
        frame = Some(current);
        for block in &observation.blocks {
            meter.charge(1)?;
            if block.query.0 >= block.query.1 || block.query.1 > current.0 {
                return Err(fail("observation_geometry"));
            }
            start = start.min(block.query.0);
            end = end.max(block.query.1);
        }
    }
    Ok(frame.map(|(axis, same)| (clips(path, strand, start, end, axis, same), same)))
}

pub(super) fn clips(
    path: &[&ObservedSamObservation],
    strand: char,
    qmin: usize,
    qmax: usize,
    axis: usize,
    same: bool,
) -> [EndClip; 2] {
    let end = |five: bool| {
        let query_left = five == same;
        let amount = if query_left { qmin } else { axis - qmax };
        if amount == 0 {
            return EndClip {
                unaligned_query_bases: Some(0),
                kind: ClipKind::None,
            };
        }
        let raw_side = usize::from(five != (strand == '+'));
        let mut kind = None;
        for observation in path {
            let terminal = observation.blocks.iter().any(|b| {
                if query_left {
                    b.query.0 == qmin
                } else {
                    b.query.1 == qmax
                }
            });
            if !terminal {
                continue;
            }
            let current = observation
                .clipping
                .as_ref()
                .and_then(|c| {
                    (c.soft[raw_side].checked_add(c.hard[raw_side]) == Some(amount)).then_some(
                        match (c.soft[raw_side] > 0, c.hard[raw_side] > 0) {
                            (true, true) => ClipKind::Mixed,
                            (true, false) => ClipKind::Soft,
                            (false, true) => ClipKind::Hard,
                            _ => ClipKind::Unknown,
                        },
                    )
                })
                .unwrap_or(ClipKind::Unknown);
            kind = Some(match kind {
                Some(previous) if previous != current => ClipKind::Unknown,
                _ => current,
            });
        }
        EndClip {
            unaligned_query_bases: Some(amount),
            kind: kind.unwrap_or(ClipKind::Unknown),
        }
    };
    [end(true), end(false)]
}

pub(super) fn annotation_index<'a>(
    document: &ModelDocument,
    annotation: &'a PrimerPayloadAnnotation,
    meter: &mut Meter,
) -> Result<BTreeMap<&'a str, &'a PrimerPayloadTranscript>, EngineError> {
    if annotation.schema != PRIMER_PAYLOAD_ANNOTATION_SCHEMA
        || annotation.assembly != document.assembly
        || annotation.reference_fasta_sha256 != document.reference_fasta_sha256
        || annotation.coordinate_convention != "zero_based_half_open"
    {
        return Err(fail("model_annotation_binding"));
    }
    if annotation.transcripts.len() > meter.limits.max_models {
        return Err(fail("model_annotation_budget"));
    }
    let mut out = BTreeMap::new();
    for transcript in &annotation.transcripts {
        meter.charge(1)?;
        if out
            .insert(transcript.transcript_id.as_str(), transcript)
            .is_some()
        {
            return Err(fail("model_annotation_identity"));
        }
    }
    Ok(out)
}

fn checked_annotation<'a>(
    id: &str,
    model: Option<&Prepared<'_>>,
    annotation: &'a BTreeMap<&str, &PrimerPayloadTranscript>,
    meter: &mut Meter,
) -> Result<&'a PrimerPayloadTranscript, EngineError> {
    meter.charge(1)?;
    let transcript = annotation
        .get(id)
        .copied()
        .ok_or_else(|| fail("model_annotation_geometry"))?;
    if let Some(model) = model {
        meter.charge(model.model.genomic_exons_0based_half_open.len())?;
        if transcript.genomic_contig.as_deref() != Some(model.model.contig.as_str())
            || transcript.target_name != model.model.contig
            || transcript.genomic_strand != Some(model.model.strand)
            || transcript.genomic_exons != model.model.genomic_exons_0based_half_open
            || !transcript.gene_ids.contains(&model.model.gene_id)
        {
            return Err(fail("model_annotation_geometry"));
        }
    }
    Ok(transcript)
}

pub(super) fn unit(
    unit: &ObservedSamUnitAssessment,
    observations: &BTreeMap<&str, &ObservedSamObservation>,
    models: &BTreeMap<&str, Prepared<'_>>,
    annotation: &BTreeMap<&str, &PrimerPayloadTranscript>,
    source: &str,
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
        candidate_transcript_ids: vec![],
        model_gene_ids: vec![],
        annotation_gene_ids: vec![],
        hypotheses: vec![],
    };
    for (index, hypothesis) in unit.hypotheses.iter().enumerate() {
        meter.charge(1)?;
        let mut record = HypothesisProjection {
            hypothesis_index: index,
            observation_ids: hypothesis.clone(),
            candidates: vec![],
            unavailable: None,
        };
        let mut path = vec![];
        let mut seen = BTreeSet::new();
        let mut common: Option<BTreeSet<&str>> = None;
        if hypothesis.is_empty() {
            record.unavailable = Some(Reason::InvalidOptimalHypothesis);
        }
        for id in hypothesis {
            meter.charge(1)?;
            if !seen.insert(id) {
                return Err(fail("hypothesis_identity"));
            }
            let Some(&observation) = observations.get(id.as_str()) else {
                record.unavailable = Some(Reason::UnboundOptimalHypothesis);
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
            let current = candidates(observation, source, meter)?;
            common = Some(match common {
                None => current,
                Some(previous) => {
                    meter.charge(previous.len())?;
                    previous.intersection(&current).copied().collect()
                }
            });
            path.push(observation);
        }
        if record.unavailable.is_none() {
            let common = common.unwrap_or_default();
            if common.is_empty() {
                record.unavailable = Some(Reason::NoCompatibleCdsModel);
            }
            for id in common {
                let model = models.get(id);
                let normalized = checked_annotation(id, model, annotation, meter)?;
                let mut relations = Vec::new();
                for observation in &path {
                    for comparison in &observation.comparisons {
                        meter.charge(1)?;
                        if comparison.source_id != source {
                            continue;
                        }
                        for model in comparison.comparison.iter().flat_map(|c| &c.models) {
                            meter.charge(1)?;
                            if model.transcript_id == id {
                                meter.charge(model.gene_ids.len())?;
                                if normalized.genomic_strand != Some(model.model_strand) {
                                    return Err(fail("model_annotation_geometry"));
                                }
                                relations.push(CandidateRelation {
                                    observation_id: observation.observation_id.clone(),
                                    relation: model.relation,
                                    gene_ids: model.gene_ids.clone(),
                                    model_strand: model.model_strand,
                                });
                            }
                        }
                    }
                }
                let genes: Vec<_> = relations
                    .iter()
                    .flat_map(|r| r.gene_ids.iter().cloned())
                    .collect::<BTreeSet<_>>()
                    .into_iter()
                    .collect();
                if genes.iter().any(|g| !normalized.gene_ids.contains(g)) {
                    return Err(fail("model_annotation_geometry"));
                }
                let mut candidate = CandidateProjection {
                    transcript_id: id.into(),
                    annotation_gene_ids: genes,
                    model_gene_ids: model
                        .map(|m| vec![m.model.gene_id.clone()])
                        .unwrap_or_default(),
                    relations,
                    transcript_intervals: vec![],
                    utr5: missing(Reason::CandidateModelUnavailable),
                    utr3: missing(Reason::CandidateModelUnavailable),
                    clip5: unknown_clip(),
                    clip3: unknown_clip(),
                    sense_to_annotation: None,
                    model_available: model.is_some(),
                    projection_available: false,
                };
                if model.is_none()
                    && let Some(strand) = normalized.genomic_strand
                    && let Some((ends, same)) = query_ends(&path, strand, meter)?
                {
                    [candidate.clip5, candidate.clip3] = ends;
                    candidate.sense_to_annotation = Some(same);
                }
                if let Some(model) = model {
                    if let Some(value) = one_path(&path, model, meter)? {
                        candidate.transcript_intervals = value.intervals;
                        [candidate.utr5, candidate.utr3] = value.metrics;
                        [candidate.clip5, candidate.clip3] = value.clips;
                        candidate.sense_to_annotation = Some(value.sense);
                        candidate.projection_available = value.available;
                    } else {
                        candidate.utr5 = missing(Reason::ModelProjectionIncompatible);
                        candidate.utr3 = missing(Reason::ModelProjectionIncompatible);
                    }
                }
                record.candidates.push(candidate);
            }
        }
        out.hypotheses.push(record);
    }
    let records: Vec<_> = out.hypotheses.iter().flat_map(|h| &h.candidates).collect();
    out.annotation_gene_ids = records
        .iter()
        .flat_map(|r| r.annotation_gene_ids.iter().cloned())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    if let Some(reason) = out.hypotheses.iter().find_map(|h| h.unavailable) {
        out.utr5 = missing(reason);
        out.utr3 = missing(reason);
        return Ok(out);
    }
    out.candidate_transcript_ids = records
        .iter()
        .map(|r| r.transcript_id.clone())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let unavailable = if records.iter().any(|r| !r.model_available) {
        Some(Reason::CandidateModelUnavailable)
    } else if records.iter().any(|r| !r.projection_available) {
        Some(Reason::ModelProjectionIncompatible)
    } else {
        None
    };
    if !records.is_empty() && records.iter().all(|r| r.model_available) {
        out.model_gene_ids = records
            .iter()
            .flat_map(|r| r.model_gene_ids.iter().cloned())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
    }
    if let Some(reason) = unavailable {
        out.utr5 = missing(reason);
        out.utr3 = missing(reason);
        return Ok(out);
    }
    if let Some(first) = records.first() {
        out.utr5 = first.utr5.clone();
        out.utr3 = first.utr3.clone();
        for candidate in &records[1..] {
            if out.utr5 != candidate.utr5 {
                out.utr5 = missing(Reason::AmbiguousTranscriptEndGeometry);
            }
            if out.utr3 != candidate.utr3 {
                out.utr3 = missing(Reason::AmbiguousTranscriptEndGeometry);
            }
        }
    }
    Ok(out)
}
