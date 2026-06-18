//! Portable ortholog and cross-species promoter comparison contracts.
//!
//! Ortholog resources are local, reviewable mapping tables. The engine resolves
//! them into promoter windows through prepared genome catalogs, then emits
//! evidence-separated comparison reports for sequence, TFBS, expression, and
//! CUT&RUN/occupancy signals.

use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// Local offline ortholog mapping resource schema.
pub const ORTHOLOG_RESOURCE_SCHEMA: &str = "gentle.ortholog_resource.v1";
/// Resolved cross-species promoter cohort report schema.
pub const ORTHOLOG_PROMOTER_COHORT_SCHEMA: &str = "gentle.ortholog_promoter_cohort.v1";
/// Cross-species promoter comparison report schema.
pub const ORTHOLOG_PROMOTER_COMPARISON_SCHEMA: &str = "gentle.ortholog_promoter_comparison.v1";

/// How to handle multiple local ortholog rows for one anchor/target pair.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum OrthologAmbiguityPolicy {
    #[default]
    Reject,
    First,
}

impl OrthologAmbiguityPolicy {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Reject => "reject",
            Self::First => "first",
        }
    }
}

/// Anchor or target row role in a resolved ortholog promoter cohort.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum OrthologPromoterRole {
    #[default]
    Anchor,
    Target,
}

/// Cross-species CUT&RUN/occupancy support status.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum OrthologCutRunSupportStatus {
    Confirmed,
    Nearby,
    MotifOnly,
    OccupancyOnly,
    #[default]
    NoData,
    NotComparable,
}

impl OrthologCutRunSupportStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Confirmed => "confirmed",
            Self::Nearby => "nearby",
            Self::MotifOnly => "motif_only",
            Self::OccupancyOnly => "occupancy_only",
            Self::NoData => "no_data",
            Self::NotComparable => "not_comparable",
        }
    }
}

/// Species alias mapping used by a local ortholog resource.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct OrthologSpeciesAlias {
    pub species: String,
    pub aliases: Vec<String>,
}

/// One directional row in a local ortholog mapping resource.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct OrthologMappingRow {
    pub source_species: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_gene_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_gene_symbol: Option<String>,
    pub target_species: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub target_gene_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub target_gene_symbol: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub orthology_type: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub confidence: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source: Option<String>,
    #[serde(default)]
    pub evidence: Vec<String>,
}

/// Reviewable local ortholog mapping resource.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct OrthologResource {
    pub schema: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    #[serde(default)]
    pub species_aliases: Vec<OrthologSpeciesAlias>,
    #[serde(default)]
    pub rows: Vec<OrthologMappingRow>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

/// Request echoed into a resolved ortholog promoter cohort.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct OrthologPromoterCohortRequest {
    pub anchor_species: String,
    pub anchor_genome_id: String,
    pub anchor_gene_query: String,
    pub target_species: Vec<String>,
    #[serde(default)]
    pub target_genome_ids: BTreeMap<String, String>,
    #[serde(default)]
    pub transcript_ids: BTreeMap<String, String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ortholog_resource_path: Option<String>,
    pub upstream_bp: usize,
    pub downstream_bp: usize,
    pub ambiguity_policy: OrthologAmbiguityPolicy,
}

/// One unresolved species/gene row from ortholog promoter resolution.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct OrthologUnresolvedRow {
    pub species: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub genome_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_query: Option<String>,
    pub reason: String,
    #[serde(default)]
    pub candidates: Vec<String>,
}

/// One promoter window resolved for an anchor or target ortholog.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct OrthologPromoterRow {
    pub species: String,
    pub genome_id: String,
    pub role: OrthologPromoterRole,
    pub gene_query: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_symbol: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_id_requested: Option<String>,
    pub transcript_id: String,
    pub display_label: String,
    pub chromosome: String,
    pub strand: String,
    pub promoter_start_1based: usize,
    pub promoter_end_1based: usize,
    pub promoter_length_bp: usize,
    pub tss_1based: usize,
    pub tss_position_0based: usize,
    pub sequence_orientation: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub promoter_sequence: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub orthology_type: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub confidence: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub orthology_source: Option<String>,
    #[serde(default)]
    pub orthology_evidence: Vec<String>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

/// Portable resolved ortholog promoter cohort.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct OrthologPromoterCohortReport {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub request: OrthologPromoterCohortRequest,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ortholog_resource_label: Option<String>,
    pub resolved_promoter_count: usize,
    pub unresolved_count: usize,
    #[serde(default)]
    pub rows: Vec<OrthologPromoterRow>,
    #[serde(default)]
    pub unresolved_rows: Vec<OrthologUnresolvedRow>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

/// Per-species motif evidence summary in an ortholog promoter comparison.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct OrthologTfbsSummaryRow {
    pub species: String,
    pub gene_label: String,
    pub transcript_id: String,
    pub tf_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tf_name: Option<String>,
    pub max_score: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub peak_position_0based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub peak_position_promoter_relative_bp: Option<i64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub peak_genomic_position_1based: Option<usize>,
    pub positive_fraction: f64,
}

/// Pairwise TFBS score-track similarity across two ortholog promoters.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct OrthologPairwiseTfbsSimilarity {
    pub left_species: String,
    pub right_species: String,
    pub left_gene_label: String,
    pub right_gene_label: String,
    pub shared_motif_count: usize,
    pub mean_raw_pearson: f64,
    pub mean_smoothed_spearman: f64,
    #[serde(default)]
    pub motif_ids: Vec<String>,
}

/// Motif peak summary shared across, or specific to, species in a cohort.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct OrthologTfbsPeakSummary {
    pub tf_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tf_name: Option<String>,
    pub promoter_count: usize,
    #[serde(default)]
    pub species: Vec<String>,
    #[serde(default)]
    pub gene_labels: Vec<String>,
    pub max_score: f64,
    #[serde(default)]
    pub peak_positions_promoter_relative_bp: Vec<i64>,
}

/// Simple promoter-sequence similarity row; kept distinct from TFBS evidence.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct OrthologSequenceSimilarityRow {
    pub left_species: String,
    pub right_species: String,
    pub left_gene_label: String,
    pub right_gene_label: String,
    pub alignment_mode: String,
    pub compared_length_bp: usize,
    pub identical_bp: usize,
    pub identity_fraction: f64,
}

/// Cross-species CUT&RUN/occupancy assignment row.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct OrthologCutRunSupportRow {
    pub species: String,
    pub gene_label: String,
    pub status: OrthologCutRunSupportStatus,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub nearest_peak_distance_bp: Option<i64>,
    #[serde(default)]
    pub contributing_dataset_ids: Vec<String>,
    #[serde(default)]
    pub contributing_read_report_ids: Vec<String>,
    pub detail: String,
}

/// Expression evidence row carried into an ortholog promoter comparison.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct OrthologExpressionAssignment {
    pub species: String,
    pub gene_label: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sample_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub condition: Option<String>,
    pub value: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub unit: Option<String>,
    pub source: String,
    pub assignment_note: String,
}

/// Portable cross-species promoter comparison report.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct OrthologPromoterComparisonReport {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub cohort: OrthologPromoterCohortReport,
    #[serde(default)]
    pub motifs_requested: Vec<String>,
    pub score_kind: String,
    pub clip_negative: bool,
    #[serde(default)]
    pub promoter_summaries: Vec<OrthologTfbsSummaryRow>,
    #[serde(default)]
    pub pairwise_tfbs_similarity: Vec<OrthologPairwiseTfbsSimilarity>,
    #[serde(default)]
    pub conserved_tfbs_peaks: Vec<OrthologTfbsPeakSummary>,
    #[serde(default)]
    pub species_specific_tfbs_peaks: Vec<OrthologTfbsPeakSummary>,
    #[serde(default)]
    pub sequence_similarity: Vec<OrthologSequenceSimilarityRow>,
    #[serde(default)]
    pub cutrun_support: Vec<OrthologCutRunSupportRow>,
    #[serde(default)]
    pub expression_assignments: Vec<OrthologExpressionAssignment>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ortholog_resource_defaults_to_offline_schema_shape() {
        let resource: OrthologResource = serde_json::from_value(serde_json::json!({
            "schema": ORTHOLOG_RESOURCE_SCHEMA,
            "species_aliases": [{"species": "Homo sapiens", "aliases": ["human"]}],
            "rows": [{
                "source_species": "Homo sapiens",
                "source_gene_symbol": "TP73",
                "target_species": "Mus musculus",
                "target_gene_symbol": "Trp73",
                "orthology_type": "one_to_one",
                "confidence": "high"
            }]
        }))
        .expect("deserialize resource");
        assert_eq!(resource.rows.len(), 1);
        assert_eq!(resource.rows[0].confidence.as_deref(), Some("high"));
        assert_eq!(resource.species_aliases[0].aliases, vec!["human"]);
    }

    #[test]
    fn old_promoter_cohort_defaults_new_optional_rows() {
        let report: OrthologPromoterCohortReport = serde_json::from_value(serde_json::json!({
            "schema": ORTHOLOG_PROMOTER_COHORT_SCHEMA,
            "generated_at_unix_ms": 1,
            "request": {
                "anchor_species": "Homo sapiens",
                "anchor_genome_id": "HumanToy",
                "anchor_gene_query": "TP73"
            }
        }))
        .expect("deserialize old-shaped cohort");
        assert_eq!(
            report.request.ambiguity_policy,
            OrthologAmbiguityPolicy::Reject
        );
        assert!(report.rows.is_empty());
        assert!(report.warnings.is_empty());
    }

    #[test]
    fn cutrun_status_serializes_snake_case() {
        let row = OrthologCutRunSupportRow {
            species: "Homo sapiens".to_string(),
            gene_label: "TP73".to_string(),
            status: OrthologCutRunSupportStatus::MotifOnly,
            detail: "motif present; no occupancy data".to_string(),
            ..OrthologCutRunSupportRow::default()
        };
        let value = serde_json::to_value(row).expect("serialize row");
        assert_eq!(value["status"], "motif_only");
    }
}
