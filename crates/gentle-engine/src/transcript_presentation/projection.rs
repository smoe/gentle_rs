//! Binding-checked local geometry; clipping never changes full-chain identity.

use gentle_protocol::transcript_presentation::*;
use gentle_protocol::{GeneLocusGenomeAnchorBinding, GeneLocusSequenceBinding};

/// One complete CDS/phase alternative, with display coordinates and original records.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LocalAnnotationRow {
    pub comparison: TranscriptSourceComparisonRow,
    pub structure_id: String,
    pub records: Vec<TranscriptPresentationRecord>,
    /// Local 1-based inclusive intervals, in transcription order.
    pub exons: Vec<TranscriptInterval>,
    /// Display-clipped intervals. Phases refer to the original source boundaries
    /// in `records`, not a newly inferred reading frame at the clipping edge.
    pub cds: Option<Vec<TranscriptCds>>,
    pub annotated_start_1based: Option<u64>,
    pub local_strand: i8,
    pub genomic_strand: i8,
    /// A cropped chain must not match a complete loaded transcript by accident.
    pub complete_chain_in_locus: bool,
}

impl LocalAnnotationRow {
    /// Exact complete-chain matching, independent of labels and CDS equivalence.
    pub fn matches_local_chain(&self, exons: &[TranscriptInterval], strand: i8) -> bool {
        if !self.complete_chain_in_locus || self.local_strand != strand || exons.is_empty() {
            return false;
        }
        let mut expected = self.exons.clone();
        let mut observed = exons.to_vec();
        expected.sort();
        observed.sort();
        expected == observed
    }
}

/// Prepared once per report. Original source provenance remains portable and intact.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LocalAnnotationComparison {
    pub sources: Vec<TranscriptSourceBinding>,
    pub rows: Vec<LocalAnnotationRow>,
}

/// Project one interval onto a validated anchor, clipping only the visible geometry.
fn local_interval(
    i: TranscriptInterval,
    a: &GeneLocusGenomeAnchorBinding,
) -> Option<TranscriptInterval> {
    let start = i.start_1based.max(a.start_1based as u64);
    let end = i.end_1based.min(a.end_1based as u64);
    if start > end {
        return None;
    }
    let (start_1based, end_1based) = if a.strand == Some('-') {
        (
            a.end_1based as u64 - end + 1,
            a.end_1based as u64 - start + 1,
        )
    } else {
        (
            start - a.start_1based as u64 + 1,
            end - a.start_1based as u64 + 1,
        )
    };
    Some(TranscriptInterval {
        start_1based,
        end_1based,
    })
}

/// Use sequence orientation, not gene strand, for the genomic-to-local transform.
/// The caller must additionally verify this saved binding against the live document.
pub fn project(
    p: &TranscriptStructurePresentation,
    binding: &GeneLocusSequenceBinding,
) -> Result<LocalAnnotationComparison, String> {
    let comparisons = super::compare_sources(p)?;
    let a = binding
        .genome_anchor
        .as_ref()
        .ok_or("Annotation comparison requires a genomic anchor")?;
    if a.start_1based == 0
        || a.end_1based < a.start_1based
        || a.end_1based - a.start_1based + 1 != binding.sequence_length_bp
        || !matches!(a.strand, Some('+' | '-'))
        || p.chromosome.strip_prefix("chr").unwrap_or(&p.chromosome)
            != a.chromosome.strip_prefix("chr").unwrap_or(&a.chromosome)
        || p.locus_sequence_sha256 != binding.sequence_sha256.trim_start_matches("sha256:")
    {
        return Err(
            "Annotation comparison sequence/anchor binding is inconsistent or unoriented".into(),
        );
    }
    let mut rows = Vec::new();
    for comparison in comparisons {
        for id in &comparison.structure_ids {
            let group = p
                .structure_groups
                .iter()
                .find(|g| &g.structure_id == id)
                .ok_or("Missing structure")?;
            let records: Vec<_> = p
                .records
                .iter()
                .filter(|r| group.member_record_ids.contains(&r.record_id))
                .cloned()
                .collect();
            let r = &records.first().ok_or("Empty structure")?.structure;
            let start = if r.strand == -1 {
                r.exons[0].interval.end_1based
            } else {
                r.exons[0].interval.start_1based
            };
            rows.push(LocalAnnotationRow {
                comparison: comparison.clone(),
                structure_id: id.clone(),
                exons: r
                    .exons
                    .iter()
                    .filter_map(|e| local_interval(e.interval, a))
                    .collect(),
                cds: group.cds.as_ref().map(|cds| {
                    cds.iter()
                        .filter_map(|c| {
                            local_interval(c.interval, a).map(|interval| TranscriptCds {
                                interval,
                                phase: c.phase,
                            })
                        })
                        .collect()
                }),
                annotated_start_1based: local_interval(
                    TranscriptInterval {
                        start_1based: start,
                        end_1based: start,
                    },
                    a,
                )
                .map(|i| i.start_1based),
                local_strand: r.strand * if a.strand == Some('-') { -1 } else { 1 },
                genomic_strand: r.strand,
                complete_chain_in_locus: r.exons.iter().all(|e| {
                    e.interval.start_1based >= a.start_1based as u64
                        && e.interval.end_1based <= a.end_1based as u64
                }),
                records,
            });
        }
    }
    Ok(LocalAnnotationComparison {
        sources: p.sources.clone(),
        rows,
    })
}
