//! C3: all eligible same-unit hypotheses must carry the allele and compatible geometry.

use super::*;

fn compatible(
    input: &GeneInput,
    geometry: &AlignmentGeometry,
    model: &CodingModel,
    policy: &Policy,
    meter: &mut Meter<'_>,
) -> Result<bool, String> {
    meter
        .charge(geometry.blocks.len() + geometry.junctions.len() + model.transcript.exons.len())?;
    let introns = model
        .transcript
        .exons
        .windows(2)
        .map(|w| (w[0].1, w[1].0))
        .collect::<Vec<_>>();
    let mut chain = Vec::new();
    for (contig, a, z) in &geometry.junctions {
        meter.charge(1)?;
        if contig != &input.report.reference.contig {
            return Ok(false);
        }
        let Some(index) = introns.iter().position(|i| i == &(*a, *z)) else {
            return Ok(false);
        };
        chain.push(index);
    }
    chain.sort_unstable();
    chain.dedup();
    if chain.windows(2).any(|w| w[0] + 1 != w[1]) {
        return Ok(false);
    }
    let mut ranges = geometry.blocks.clone();
    ranges.sort();
    let mut merged = Vec::<(String, usize, usize)>::new();
    for (contig, a, z) in ranges {
        if a >= z {
            return Err(fail("linkage_geometry"));
        }
        if let Some(last) = merged.last_mut()
            && last.0 == contig
            && a <= last.2
        {
            last.2 = last.2.max(z);
        } else {
            merged.push((contig, a, z));
        }
    }
    if merged.is_empty() {
        return Ok(false);
    }
    let mut outside = 0usize;
    for (contig, a, z) in merged {
        meter.charge(model.transcript.exons.len() + 1)?;
        if contig != input.report.reference.contig {
            return Ok(false);
        }
        let covered = if contig == input.report.reference.contig {
            model
                .transcript
                .exons
                .iter()
                .map(|&(x, y)| z.min(y).saturating_sub(a.max(x)))
                .sum()
        } else {
            0
        };
        outside = outside
            .checked_add((z - a).saturating_sub(covered))
            .ok_or_else(|| fail("linkage_geometry"))?;
        if outside > policy.maximum_out_of_model_aligned_bases {
            return Ok(false);
        }
    }
    Ok(true)
}

pub(super) fn link(
    input: &GeneInput,
    evidence: &LinkageEvidence,
    sample: &evidence::CdnaSampleReport,
    model: &CodingModel,
    p: usize,
    alt: char,
    policy: &Policy,
    meter: &mut Meter<'_>,
) -> Result<IsoformLinkage, String> {
    let mut result = IsoformLinkage {
        status: "no_compatible_units".into(),
        compatible_units: Vec::new(),
        discriminating_units: Vec::new(),
        unique_model_units: Vec::new(),
        excluded_transcripts: BTreeMap::new(),
        reasons: Vec::new(),
    };
    if position(model, p).is_none() {
        result.status = "not_present".into();
        return Ok(result);
    }
    let source = evidence
        .samples
        .iter()
        .find(|s| s.sample_id == sample.sample_id)
        .ok_or_else(|| fail("linkage_sample"))?;
    let competitors = input
        .models
        .iter()
        .filter(|m| position(m, p).is_some())
        .collect::<Vec<_>>();
    for call in sample
        .site_calls
        .iter()
        .filter(|c| c.position_0based == p && c.base == Some(alt) && definite(c))
    {
        meter.charge(
            source.units.len()
                + input
                    .models
                    .iter()
                    .map(|m| m.transcript.exons.len())
                    .sum::<usize>()
                + 1,
        )?;
        let unit = source
            .units
            .iter()
            .find(|u| u.unit_id == call.unit_id)
            .ok_or_else(|| fail("linkage_unit"))?;
        if !unit.all_query_views_assessed {
            continue;
        }
        let hypotheses = unit
            .hypotheses
            .iter()
            .filter(|h| h.eligible)
            .collect::<Vec<_>>();
        if hypotheses.is_empty() {
            continue;
        }
        let mut admissible = Vec::new();
        for candidate in &competitors {
            let mut all = true;
            for hypothesis in &hypotheses {
                meter.charge(hypothesis.bases.len() + 1)?;
                let carries = hypothesis.bases.iter().any(|b| {
                    b.reference_position_0based == p
                        && b.observed_base == alt
                        && b.phred
                            .is_some_and(|q| q >= input.report.policy.minimum_base_quality)
                        && b.mapq.is_some_and(|q| {
                            q != 255 && q >= input.report.policy.minimum_mapping_quality
                        })
                });
                let geometry = evidence
                    .geometry
                    .get(&(
                        sample.sample_id.clone(),
                        unit.unit_id.clone(),
                        hypothesis.hypothesis_id.clone(),
                    ))
                    .ok_or_else(|| fail("linkage_hypothesis"))?;
                if !carries || !compatible(input, geometry, candidate, policy, meter)? {
                    all = false;
                    break;
                }
            }
            if all {
                admissible.push(candidate.transcript.transcript_id.as_str());
            }
        }
        if !admissible.contains(&model.transcript.transcript_id.as_str()) {
            continue;
        }
        result.compatible_units.push(unit.unit_id.clone());
        let excluded = competitors
            .iter()
            .filter(|c| !admissible.contains(&c.transcript.transcript_id.as_str()))
            .collect::<Vec<_>>();
        if !excluded.is_empty() {
            result.discriminating_units.push(unit.unit_id.clone());
            for c in excluded {
                result
                    .excluded_transcripts
                    .entry(c.transcript.transcript_id.clone())
                    .or_default()
                    .push(unit.unit_id.clone());
            }
        }
        if admissible.len() == 1 {
            result.unique_model_units.push(unit.unit_id.clone());
        }
    }
    result.compatible_units.sort();
    result.discriminating_units.sort();
    result.unique_model_units.sort();
    for units in result.excluded_transcripts.values_mut() {
        units.sort();
        units.dedup();
    }
    result.status = if !result.discriminating_units.is_empty() {
        "same_unit_discriminating_support"
    } else if !result.compatible_units.is_empty() {
        "same_unit_compatible_only"
    } else {
        "no_compatible_units"
    }
    .into();
    result.reasons.push("Compatibility is conditional on supplied eligible hypotheses, not identification of a full-length or expressed isoform.".into());
    Ok(result)
}
