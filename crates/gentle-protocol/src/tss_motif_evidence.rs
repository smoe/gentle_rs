//! Pure coordinate and admission rules for imported TSS motif annotations.
//! No scoring, database access, contig aliasing or biological absence inference.

use crate::genomic_motif_evidence::*;
use crate::tss_profiles::*;
use std::collections::BTreeSet;

pub const MAX_EVIDENCE_REPORTS: usize = 64;
pub const MAX_EVIDENCE_HITS: usize = 100_000;

/// An original hit plus its clipped interval in the displayed sequence.
pub struct ProjectedMotifHit<'a> {
    pub hit: &'a GenomicMotifEvidenceHit,
    pub start: usize,
    pub end: usize,
    /// Original unclipped footprint in displayed-axis coordinates, possibly outside the window.
    pub full_start: i128,
    pub full_end: i128,
    pub local_strand: TssStrand,
    pub clipped: bool,
}

/// Check declared identity and geometry. Source-byte hashes are verified by the
/// engine before export; these checks alone do not authenticate a genome/package.
pub fn validate(report: &TssProfileReport) -> Result<(), String> {
    let fail = |why: &str| Err(format!("Imported motif evidence: {why}"));
    if report.imported_motif_evidence.len() > MAX_EVIDENCE_REPORTS {
        return fail("too many source reports");
    }
    let mut total = 0usize;
    let mut ids = BTreeSet::new();
    for source in &report.imported_motif_evidence {
        let evidence = &source.report;
        total = total.saturating_add(evidence.hits.len());
        if total > MAX_EVIDENCE_HITS {
            return fail("too many retained hits");
        }
        if evidence.schema != GENOMIC_MOTIF_EVIDENCE_SCHEMA
            || evidence.availability != GenomicMotifEvidenceAvailability::Available
            || evidence.report_id.is_empty()
            || !ids.insert(&evidence.report_id)
            || evidence.returned_hit_count != evidence.hits.len()
            || (evidence.query_complete && evidence.truncated)
            || evidence.regions.is_empty()
            || evidence.request.motif_ids.is_empty()
            || evidence.request.motif_ids.len() > MAX_GENOMIC_MOTIF_EVIDENCE_QUERY_MOTIFS
        {
            return fail("expected distinct available query reports with consistent hit counts");
        }
        let provider = evidence
            .provider
            .as_ref()
            .ok_or("Imported motif evidence: missing provider")?;
        if provider.genome_id != report.reference.genome_id
            || ![
                provider.assembly_name.as_deref(),
                provider.assembly_accession.as_deref(),
            ]
            .contains(&Some(report.reference.assembly.as_str()))
        {
            return fail("genome ID/assembly must match exactly; no aliases or liftover");
        }
        if provider.coordinate_mode != "bed_0based_half_open"
            || provider.score_mode.is_empty()
            || provider.run_id.is_empty()
            || provider.manifest_sha256.is_empty()
        {
            return fail("missing score/provenance or unsupported coordinates");
        }
        let mut regions = BTreeSet::new();
        for region in &evidence.regions {
            if region.start_0based >= region.end_0based_exclusive
                || !regions.insert(&region.interval_id)
                || region
                    .resolved_chromosome
                    .as_ref()
                    .is_none_or(|s| s.is_empty())
                || !matches!(
                    region.compatibility_status,
                    GenomicMotifEvidenceCompatibilityStatus::ContigGeometryMatchedOnly
                        | GenomicMotifEvidenceCompatibilityStatus::AssemblyAndContigGeometryMatched
                )
            {
                return fail("invalid or incompatible query region");
            }
        }
        for motif in &evidence.request.motif_ids {
            if evidence
                .motif_coverage
                .iter()
                .filter(|c| &c.motif_id == motif)
                .count()
                != 1
            {
                return fail("each requested accession requires one coverage record");
            }
        }
        if evidence
            .request
            .motif_ids
            .iter()
            .collect::<BTreeSet<_>>()
            .len()
            != evidence.request.motif_ids.len()
            || evidence.motif_coverage.len() != evidence.request.motif_ids.len()
            || evidence.motif_coverage.iter().any(|c| {
                c.returned_hit_count
                    != evidence
                        .hits
                        .iter()
                        .filter(|h| h.motif_id == c.motif_id)
                        .count()
            })
        {
            return fail("coverage counts or requested accessions are inconsistent");
        }
        for hit in &evidence.hits {
            let region = evidence
                .regions
                .iter()
                .find(|r| r.interval_id == hit.interval_id)
                .ok_or("Imported motif evidence: hit has unknown query region")?;
            if hit.start_0based >= hit.end_0based_exclusive
                || hit.end_0based_exclusive > i64::MAX as u64
                || Some(hit.chromosome.as_str()) != region.resolved_chromosome.as_deref()
                || hit.end_0based_exclusive <= region.start_0based
                || hit.start_0based >= region.end_0based_exclusive
                || !matches!(hit.strand.as_str(), "+" | "-")
                || !hit.score.is_finite()
                || hit.pwm_relative_score.is_some_and(|s| !s.is_finite())
                || hit.score_mode != provider.score_mode
                || !evidence.request.motif_ids.contains(&hit.motif_id)
            {
                return fail("invalid hit score, accession, direction or region binding");
            }
            if let Some(matrix) = report
                .panel_resolution
                .matrices
                .iter()
                .find(|m| m.specification.source_id == hit.motif_id)
                && hit.end_0based_exclusive - hit.start_0based != matrix.matrix_counts.len() as u64
            {
                return fail("hit width differs from its exact matrix accession");
            }
        }
    }
    Ok(())
}

/// Clip using the shared transcript-oriented window convention; genomic strand
/// remains available on the original hit. No dependence on query-local indices.
pub fn project<'a>(
    evidence: &'a GenomicMotifEvidenceReport,
    geometry: &TssGeometry,
    accession: &str,
) -> Vec<ProjectedMotifHit<'a>> {
    let low = geometry.start_1based - 1;
    let high = geometry.end_1based;
    evidence
        .hits
        .iter()
        .filter_map(|hit| {
            if hit.motif_id != accession
                || hit.chromosome != geometry.chromosome
                || hit.end_0based_exclusive <= low
                || hit.start_0based >= high
            {
                return None;
            }
            let a = hit.start_0based.max(low);
            let b = hit.end_0based_exclusive.min(high);
            let (start, end) = match geometry.strand {
                TssStrand::Plus => (a - low, b - low),
                TssStrand::Minus => (high - b, high - a),
            };
            let (full_start, full_end) = match geometry.strand {
                TssStrand::Plus => (
                    hit.start_0based as i128 - low as i128,
                    hit.end_0based_exclusive as i128 - low as i128,
                ),
                TssStrand::Minus => (
                    high as i128 - hit.end_0based_exclusive as i128,
                    high as i128 - hit.start_0based as i128,
                ),
            };
            Some(ProjectedMotifHit {
                hit,
                start: start as usize,
                end: end as usize,
                full_start,
                full_end,
                local_strand: if hit.strand == geometry.strand.as_str() {
                    TssStrand::Plus
                } else {
                    TssStrand::Minus
                },
                clipped: a != hit.start_0based || b != hit.end_0based_exclusive,
            })
        })
        .collect()
}

/// Covers only the actual query region union, not a bounding box spanning gaps.
pub fn covers_window(evidence: &GenomicMotifEvidenceReport, geometry: &TssGeometry) -> bool {
    let mut spans = evidence
        .regions
        .iter()
        .filter(|r| r.resolved_chromosome.as_deref() == Some(geometry.chromosome.as_str()))
        .map(|r| (r.start_0based, r.end_0based_exclusive))
        .collect::<Vec<_>>();
    spans.sort_unstable();
    let mut cursor = geometry.start_1based - 1;
    for (start, end) in spans {
        if start > cursor {
            break;
        }
        cursor = cursor.max(end);
        if cursor >= geometry.end_1based {
            return true;
        }
    }
    false
}

/// One provider/report/accession scale across all its TSS windows and both
/// strands. Signed raw scores are never confused with the triangle's direction.
pub fn score_range(evidence: &GenomicMotifEvidenceReport, accession: &str) -> (f64, f64) {
    let (mut min, mut max) = (0.0_f64, 0.0_f64);
    for hit in evidence.hits.iter().filter(|h| h.motif_id == accession) {
        min = min.min(hit.score);
        max = max.max(hit.score);
    }
    if min == max {
        max = 1.0;
    }
    (min, max)
}
