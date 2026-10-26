//! Versioned exon-based source admission; record orientation is not transcript strand.

use super::*;

type Interval = (usize, usize); // 1-based inclusive, ascending genomic coordinates.
type Mapping = Result<MappedLocus, RnaReadPanelOverlapFallbackReason>;

struct MappedLocus {
    templates: Vec<MappedTemplate>,
    unions: [Vec<Interval>; 2],
}

struct MappedTemplate {
    member_id: String,
    strand: char,
    exons: Vec<Interval>,
}

fn anchor_span(anchor: &SequenceGenomeAnchorSummary) -> Option<usize> {
    if anchor.start_1based == 0 {
        return None;
    }
    anchor
        .end_1based
        .checked_sub(anchor.start_1based)?
        .checked_add(1)
}

fn union(mut intervals: Vec<Interval>) -> Vec<Interval> {
    intervals.sort_unstable();
    let mut merged: Vec<Interval> = Vec::new();
    for (start, end) in intervals {
        if let Some(last) = merged.last_mut()
            && (start <= last.1 || last.1.checked_add(1) == Some(start))
        {
            last.1 = last.1.max(end);
        } else {
            merged.push((start, end));
        }
    }
    merged
}

fn intersections(a: &[Interval], b: &[Interval]) -> Vec<Interval> {
    let (mut i, mut j) = (0, 0);
    let mut result = Vec::new();
    while i < a.len() && j < b.len() {
        let start = a[i].0.max(b[j].0);
        let end = a[i].1.min(b[j].1);
        if start <= end {
            result.push((start, end));
        }
        if a[i].1 < b[j].1 {
            i += 1;
        } else {
            j += 1;
        }
    }
    result
}

fn map_locus(locus: &RnaReadPanelLocus, templates: &[RnaReadPanelTemplate]) -> Mapping {
    use RnaReadPanelOverlapFallbackReason::*;
    let anchor = locus.genome_anchor.as_ref().ok_or(InvalidAnchorGeometry)?;
    let span = anchor_span(anchor).ok_or(InvalidAnchorGeometry)?;
    if anchor.seq_id != locus.source.seq_id || anchor.chromosome.trim().is_empty() {
        return Err(InvalidAnchorGeometry);
    }
    let orientation = match anchor.strand {
        Some(strand @ ('+' | '-')) => strand,
        _ => return Err(UnknownAnchorStrand),
    };
    let mut result = Vec::new();
    for template in templates
        .iter()
        .filter(|t| t.locus_id == locus.source.locus_id)
    {
        let strand = match (orientation, template.strand.as_str()) {
            ('+', "+") | ('-', "-") => '+',
            ('+', "-") | ('-', "+") => '-',
            _ => return Err(UnmappableExons),
        };
        let mut exons = Vec::new();
        for &(start, end) in &template.source_exons {
            if start >= end || end > span {
                return Err(UnmappableExons);
            }
            // The shared mapper defaults unknown strands to +; the explicit guard
            // above is required here so absence never becomes a known orientation.
            let interval = GentleEngine::genomic_interval_from_anchor(anchor, start, end)
                .filter(|&(s, e)| s >= anchor.start_1based && e <= anchor.end_1based)
                .ok_or(UnmappableExons)?;
            exons.push(interval);
        }
        if exons.is_empty() {
            return Err(UnmappableExons);
        }
        result.push(MappedTemplate {
            member_id: template.member_id.clone(),
            strand,
            exons: union(exons),
        });
    }
    if result.is_empty() {
        return Err(UnmappableExons);
    }
    Ok(MappedLocus {
        unions: [strand_union(&result, '+'), strand_union(&result, '-')],
        templates: result,
    })
}

fn strand_union(templates: &[MappedTemplate], strand: char) -> Vec<Interval> {
    union(
        templates
            .iter()
            .filter(|t| t.strand == strand)
            .flat_map(|t| t.exons.iter().copied())
            .collect(),
    )
}

fn overlapping_members(
    templates: &[MappedTemplate],
    strand: char,
    other: &[Interval],
) -> Vec<String> {
    templates
        .iter()
        .filter(|t| t.strand == strand && !intersections(&t.exons, other).is_empty())
        .map(|t| t.member_id.clone())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect()
}

pub(super) fn assess(
    panel: &RnaReadPanel,
) -> (Vec<(String, String)>, RnaReadPanelSourceOverlapAssessment) {
    let mut loci = panel.loci.iter().collect::<Vec<_>>();
    loci.sort_by(|a, b| a.source.locus_id.cmp(&b.source.locus_id));
    let mapped = loci
        .iter()
        .map(|l| map_locus(l, &panel.templates))
        .collect::<Vec<_>>();
    let mut pairs = BTreeSet::new();
    let mut assessment = RnaReadPanelSourceOverlapAssessment {
        rule: RnaReadPanelSourceOverlapRule::GenomicExonsV1,
        antisense_overlapping_loci: vec![],
        conservative_source_pairs: vec![],
    };
    for (i, a) in loci.iter().enumerate() {
        for (j, b) in loci.iter().enumerate().skip(i + 1) {
            let ids = (a.source.locus_id.clone(), b.source.locus_id.clone());
            if same_source_identity(a, b) {
                pairs.insert(ids.clone());
            }
            let (Some(x), Some(y)) = (&a.genome_anchor, &b.genome_anchor) else {
                continue;
            };
            if a.source.assembly != b.source.assembly || x.chromosome != y.chromosome {
                continue;
            }
            match (&mapped[i], &mapped[j]) {
                (Ok(left), Ok(right)) => {
                    let left_unions = &left.unions;
                    let right_unions = &right.unions;
                    for (index, strand) in ['+', '-'].into_iter().enumerate() {
                        if !intersections(&left_unions[index], &right_unions[index]).is_empty() {
                            pairs.insert(ids.clone());
                        }
                        let other = 1 - index;
                        let intervals = intersections(&left_unions[index], &right_unions[other]);
                        if !intervals.is_empty() {
                            let opposite = if strand == '+' { '-' } else { '+' };
                            assessment.antisense_overlapping_loci.push(
                                RnaReadPanelAntisenseOverlap {
                                    locus_ids: ids.clone(),
                                    assembly: a.source.assembly.clone(),
                                    chromosome: x.chromosome.clone(),
                                    genomic_strands: (strand, opposite),
                                    member_ids_a: overlapping_members(
                                        &left.templates,
                                        strand,
                                        &right_unions[other],
                                    ),
                                    member_ids_b: overlapping_members(
                                        &right.templates,
                                        opposite,
                                        &left_unions[index],
                                    ),
                                    overlapping_intervals_1based: intervals,
                                },
                            );
                        }
                    }
                }
                _ => {
                    // Invalid bounds cannot establish disjointness. Otherwise use
                    // whole-record overlap, regardless of extraction orientation.
                    if anchor_span(x).is_none()
                        || anchor_span(y).is_none()
                        || (x.start_1based <= y.end_1based && y.start_1based <= x.end_1based)
                    {
                        pairs.insert(ids.clone());
                        let reasons = [mapped[i].as_ref().err(), mapped[j].as_ref().err()]
                            .into_iter()
                            .flatten()
                            .cloned()
                            .collect::<BTreeSet<_>>()
                            .into_iter()
                            .collect();
                        assessment.conservative_source_pairs.push(
                            RnaReadPanelConservativeOverlap {
                                locus_ids: ids,
                                assembly: a.source.assembly.clone(),
                                chromosome: x.chromosome.clone(),
                                reasons,
                            },
                        );
                    }
                }
            }
        }
    }
    assessment.antisense_overlapping_loci.sort();
    assessment.conservative_source_pairs.sort();
    (pairs.into_iter().collect(), assessment)
}

/// Summaries retain a bound source observation, not the full geometry needed to recompute it.
/// Validate its portable structure without claiming fresh source verification.
pub(super) fn validate_summary(summary: &RnaReadPanelSummary) -> Result<(), EngineError> {
    let Some(assessment) = &summary.source_overlap_assessment else {
        return Ok(());
    };
    let valid_ids = |ids: &(String, String)| safe_id(&ids.0) && safe_id(&ids.1) && ids.0 < ids.1;
    let valid_members = |ids: &[String], locus: &str| {
        !ids.is_empty()
            && !ids.windows(2).any(|p| p[0] >= p[1])
            && ids.iter().all(|id| {
                id.rsplit_once(':')
                    .is_some_and(|(l, f)| l == locus && f.parse::<usize>().is_ok())
            })
    };
    if assessment
        .antisense_overlapping_loci
        .windows(2)
        .any(|p| p[0] >= p[1])
        || assessment
            .conservative_source_pairs
            .windows(2)
            .any(|p| p[0] >= p[1])
        || assessment.antisense_overlapping_loci.iter().any(|record| {
            !valid_ids(&record.locus_ids)
                || record.assembly.trim().is_empty()
                || record.chromosome.trim().is_empty()
                || !matches!(record.genomic_strands, ('+', '-') | ('-', '+'))
                || !valid_members(&record.member_ids_a, &record.locus_ids.0)
                || !valid_members(&record.member_ids_b, &record.locus_ids.1)
                || record.overlapping_intervals_1based.is_empty()
                || record
                    .overlapping_intervals_1based
                    .iter()
                    .any(|&(s, e)| s == 0 || s > e)
                || union(record.overlapping_intervals_1based.clone())
                    != record.overlapping_intervals_1based
        })
        || assessment.conservative_source_pairs.iter().any(|record| {
            !valid_ids(&record.locus_ids)
                || record.assembly.trim().is_empty()
                || record.chromosome.trim().is_empty()
                || record.reasons.is_empty()
                || record.reasons.windows(2).any(|p| p[0] >= p[1])
        })
    {
        return Err(bad("panel_summary_source_overlap_binding"));
    }
    Ok(())
}
