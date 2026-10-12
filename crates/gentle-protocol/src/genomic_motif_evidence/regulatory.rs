//! Annotation-selected motif evidence. Physical hits and n:m gene ownership
//! remain separate; completeness is only within the declared annotation subset.

use serde::{Deserialize, Serialize};
use serde_json::Value;

pub const REGULATORY_MOTIF_PROVIDER: &str = "genome_regulatory_tfbs_subset";

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum RegulatoryMotifCoverageState {
    Available,
    KnownEmptyIntersection,
    ChromosomeNotInPackage,
    MotifNotInPackage,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct RegulatoryMotifCoverage {
    pub chromosome: String,
    pub motif_id: String,
    pub state: RegulatoryMotifCoverageState,
}

/// Exact source-file bytes verified during this query, not a digital signature.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct RegulatoryMotifFileBinding {
    pub path: String,
    pub bytes: u64,
    pub sha256: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct RegulatoryMotifTssWindow {
    pub promoter_id: String,
    pub tss_id: String,
    pub chromosome: String,
    pub start_0based: u64,
    pub end_0based_exclusive: u64,
    pub tss_0based: u64,
    pub strand: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct RegulatoryMotifTranscriptOwner {
    pub tss_id: String,
    pub transcript_id: String,
    pub gene_id: String,
    pub gene_name: Option<String>,
}

/// Index into the report's physical hit list; ownership does not duplicate hits.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct RegulatoryMotifHitAnnotation {
    pub hit_index: usize,
    pub query_interval_ids: Vec<String>,
    pub regulation_tags: u16,
    pub overlaps_regulatory_tss_intersection: bool,
    /// IDs of overlapping physical windows, not predicted target genes.
    pub promoter_ids: Vec<String>,
    pub regulatory_feature_ids: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
pub struct RegulatoryMotifSubset {
    pub scope: String,
    pub score_selection: String,
    pub complete_genome_scan: bool,
    pub requires_tp73: bool,
    pub annotation_release: String,
    pub regulatory_release: String,
    pub promoter_definition_id: String,
    pub upstream_bp: u64,
    pub downstream_bp: u64,
    pub production_plan_sha256: String,
    pub source_commit: String,
    /// None means not supplied/assessed, not a zero pseudocount. In particular
    /// catalog inspection and empty payloads may not expose a scoring policy.
    pub score_pseudocount: Option<f64>,
    pub verified_files: Vec<RegulatoryMotifFileBinding>,
    pub coverage: Vec<RegulatoryMotifCoverage>,
    /// Original catalog rows retain taxa and provider metadata without human-only filtering.
    pub motif_metadata: Vec<Value>,
    pub catalog_total_motifs: usize,
    pub catalog_matched_motifs: usize,
    pub catalog_offset: usize,
    pub catalog_has_more: bool,
    pub tss_windows: Vec<RegulatoryMotifTssWindow>,
    pub transcript_owners: Vec<RegulatoryMotifTranscriptOwner>,
    pub regulatory_features: Vec<Value>,
    /// Source annotation records; these links are not nearest-gene inference.
    pub regulatory_gene_links: Vec<Value>,
    pub hit_annotations: Vec<RegulatoryMotifHitAnnotation>,
}

impl RegulatoryMotifSubset {
    /// Short portable non-claim for GUI and publication renderers.
    pub fn summary(&self) -> String {
        format!(
            "Regulatory/TSS intersection only; {} (-{}/+{} bp, TSS included); transcript annotation {}; regulatory release {}. Not a complete genome scan or experimental binding evidence.",
            self.promoter_definition_id,
            self.upstream_bp,
            self.downstream_bp,
            self.annotation_release,
            self.regulatory_release
        )
    }
}

/// Replay admission: never lose an annotation-selected source's scope or turn
/// its absent rows into an exhaustive full-genome claim. Byte verification is
/// performed by the reader/exporter; this validates portable internal bindings.
pub fn validate_regulatory_subset(
    report: &super::GenomicMotifEvidenceReport,
) -> Result<(), String> {
    use std::collections::BTreeSet;
    let kind = report.provider.as_ref().map(|p| p.provider_kind.as_str());
    let Some(subset) = &report.regulatory_subset else {
        return if kind == Some(REGULATORY_MOTIF_PROVIDER) {
            Err("Missing regulatory subset provenance".into())
        } else {
            Ok(())
        };
    };
    let digest = |s: &str| {
        s.len() == 64
            && s.bytes()
                .all(|b| b.is_ascii_hexdigit() && !b.is_ascii_uppercase())
    };
    if kind != Some(REGULATORY_MOTIF_PROVIDER)
        || subset.scope != "regulatory_and_tss"
        || subset.score_selection != "source_retention"
        || subset.complete_genome_scan
        || subset.requires_tp73
        || !digest(&subset.production_plan_sha256)
        || [
            &subset.annotation_release,
            &subset.regulatory_release,
            &subset.promoter_definition_id,
            &subset.source_commit,
        ]
        .iter()
        .any(|s| s.is_empty())
        || subset.hit_annotations.len() != report.hits.len()
        || subset
            .verified_files
            .iter()
            .any(|f| !digest(&f.sha256) || f.path.is_empty())
    {
        return Err("Invalid regulatory subset provenance or hit bindings".into());
    }
    for path in ["file_inventory.json", "annotation/manifest.json"] {
        if subset
            .verified_files
            .iter()
            .filter(|f| f.path == path)
            .count()
            != 1
        {
            return Err("Missing/duplicate regulatory metadata binding".into());
        }
    }
    let mut coverage = BTreeSet::new();
    for c in &subset.coverage {
        if !coverage.insert((&c.chromosome, &c.motif_id)) {
            return Err("Duplicate regulatory coverage row".into());
        }
        if c.state == RegulatoryMotifCoverageState::KnownEmptyIntersection
            && report
                .hits
                .iter()
                .any(|h| h.chromosome == c.chromosome && h.motif_id == c.motif_id)
        {
            return Err("Known-empty intersection contains hits".into());
        }
        if report.query_complete
            && matches!(
                c.state,
                RegulatoryMotifCoverageState::ChromosomeNotInPackage
                    | RegulatoryMotifCoverageState::MotifNotInPackage
            )
        {
            return Err("Unavailable subset coverage cannot be complete".into());
        }
    }
    let expected_coverage: BTreeSet<_> = report
        .regions
        .iter()
        .flat_map(|r| {
            report
                .request
                .motif_ids
                .iter()
                .map(move |m| (&r.requested_chromosome, m))
        })
        .collect();
    if coverage != expected_coverage {
        return Err("Missing or unexpected regulatory coverage rows".into());
    }
    if subset.score_pseudocount.is_some_and(|p| {
        !p.is_finite() || p < 0.0 || report.provider.as_ref().is_none_or(|v| v.pseudocount != p)
    }) {
        return Err("Invalid regulatory scoring policy".into());
    }
    let mut seen = BTreeSet::new();
    for a in &subset.hit_annotations {
        let h = report
            .hits
            .get(a.hit_index)
            .ok_or("Unknown regulatory hit index")?;
        if !seen.insert(a.hit_index)
            || !a.overlaps_regulatory_tss_intersection
            || a.regulation_tags > 255
            || a.regulation_tags & 127 == 0
            || a.regulation_tags & 128 == 0
            || !a.query_interval_ids.contains(&h.interval_id)
            || a.query_interval_ids.iter().any(|id| {
                !report.regions.iter().any(|r| {
                    r.interval_id == *id
                        && r.resolved_chromosome.as_deref() == Some(h.chromosome.as_str())
                        && r.start_0based < h.end_0based_exclusive
                        && r.end_0based_exclusive > h.start_0based
                })
            })
            || a.promoter_ids.iter().any(|id| {
                !subset.tss_windows.iter().any(|w| {
                    w.promoter_id == *id
                        && w.chromosome == h.chromosome
                        && w.start_0based < h.end_0based_exclusive
                        && w.end_0based_exclusive > h.start_0based
                })
            })
            || a.regulatory_feature_ids.iter().any(|id| {
                !subset
                    .regulatory_features
                    .iter()
                    .any(|f| f["regulatory_feature_id"] == *id)
            })
        {
            return Err("Invalid regulatory hit/tag/ownership relation".into());
        }
    }
    let mut windows = BTreeSet::new();
    for w in &subset.tss_windows {
        if !windows.insert(&w.promoter_id)
            || w.start_0based > w.tss_0based
            || w.tss_0based >= w.end_0based_exclusive
            || !matches!(w.strand.as_str(), "+" | "-")
        {
            return Err("Invalid/duplicate stored TSS window".into());
        }
    }
    if subset.transcript_owners.iter().any(|o| {
        o.transcript_id.is_empty()
            || o.gene_id.is_empty()
            || !subset.tss_windows.iter().any(|w| w.tss_id == o.tss_id)
    }) {
        return Err("Invalid transcript owner".into());
    }
    Ok(())
}
