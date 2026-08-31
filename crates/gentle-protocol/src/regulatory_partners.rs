//! Portable regulatory-partner motif tuple and decision-tree contracts.
//!
//! These records keep exact promoter-oriented and genomic coordinates so a
//! client can replay the evidence trace without rescanning DNA. Motif and
//! promoter occupancy evidence remain associations, not causal assignments.

use crate::{GeneSetCutRunEvaluationState, GeneSetPromoterCohortReport, GeneSetUnresolvedMember};
use serde::{Deserialize, Serialize};

pub const REGULATORY_PARTNER_TUPLE_LEDGER_SCHEMA: &str =
    "gentle.regulatory_partner_tuple_ledger.v1";
pub const REGULATORY_PARTNER_DECISION_TREE_SCHEMA: &str =
    "gentle.regulatory_partner_decision_tree.v1";
pub const REGULATORY_PARTNER_SCREEN_SCHEMA: &str = "gentle.regulatory_partner_screen.v1";

/// How promoter occupancy participates in the evidence trace.
#[derive(
    Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Hash, Default,
)]
#[serde(rename_all = "snake_case")]
pub enum RegulatoryPartnerAnchorMode {
    /// Prefer promoter-level occupancy support, while preserving unknown and
    /// evaluated-with-zero-support as distinct outcomes.
    #[default]
    OccupancyPreferred,
    /// Evaluate motif tuples without treating occupancy as a required branch.
    MotifOnly,
}

/// Role played by one exact motif hit in the screen.
#[derive(
    Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Hash, Default,
)]
#[serde(rename_all = "snake_case")]
pub enum RegulatoryPartnerMotifRole {
    #[default]
    Anchor,
    Partner,
}

/// Exact relationship between the anchor and partner intervals.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RegulatoryPartnerIntervalRelation {
    Overlap,
    Abutting,
    Contained,
    #[default]
    Disjoint,
}

/// Motif-strand relationship in promoter-oriented coordinates.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RegulatoryPartnerOrientation {
    Same,
    Convergent,
    Divergent,
    #[default]
    Ambiguous,
}

/// Outcome of one evidence decision for one gene.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RegulatoryPartnerDecisionOutcome {
    Pass,
    Fail,
    #[default]
    Unknown,
}

/// Absolute threshold override for one motif token or resolved motif id.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct RegulatoryPartnerMotifThreshold {
    pub motif: String,
    pub min_llr_bits: f64,
}

/// Parameters that fully declare promoter motif enumeration and tuple rules.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct RegulatoryPartnerScreenRequest {
    pub genome_id: String,
    pub anchor_motifs: Vec<String>,
    pub partner_motifs: Vec<String>,
    pub default_min_llr_bits: f64,
    pub per_motif_thresholds: Vec<RegulatoryPartnerMotifThreshold>,
    pub max_anchor_partner_distance_bp: usize,
    pub upstream_bp: usize,
    pub downstream_bp: usize,
    pub anchor_mode: RegulatoryPartnerAnchorMode,
}

/// Content-bound pointer to an input report embedded or summarized by a run.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct RegulatoryPartnerSourceReportRef {
    pub report_kind: String,
    pub schema: String,
    pub content_sha256: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
}

/// One uncapped threshold-passing motif hit with both coordinate systems.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct RegulatoryPartnerMotifHit {
    pub hit_id: String,
    pub member_dedup_key: String,
    pub gene_symbol: String,
    pub role: RegulatoryPartnerMotifRole,
    pub tf_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tf_name: Option<String>,
    pub motif_consensus_iupac: String,
    pub motif_length_bp: usize,
    pub promoter_start_0based: usize,
    pub promoter_end_0based_exclusive: usize,
    pub genomic_start_1based: usize,
    pub genomic_end_1based: usize,
    pub promoter_forward_strand: bool,
    pub genomic_forward_strand: bool,
    pub matched_sequence: String,
    pub llr_bits: f64,
    pub llr_quantile_within_promoter: f64,
    pub signed_start_to_tss_bp: i64,
    pub signed_center_to_tss_bp: i64,
}

/// One exact anchor-plus-partner pair retained in the tuple ledger.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct RegulatoryPartnerTupleRow {
    pub tuple_id: String,
    pub member_dedup_key: String,
    pub gene_symbol: String,
    pub anchor_hit_id: String,
    pub partner_hit_id: String,
    pub signed_anchor_to_partner_center_distance_bp: f64,
    pub signed_anchor_to_partner_edge_distance_bp: i64,
    pub absolute_anchor_to_partner_distance_bp: usize,
    pub interval_relation: RegulatoryPartnerIntervalRelation,
    pub orientation: RegulatoryPartnerOrientation,
    pub within_requested_distance: bool,
    pub occupancy_evaluation_state: GeneSetCutRunEvaluationState,
    pub promoter_occupancy_supported: bool,
    #[serde(default)]
    pub occupancy_dataset_ids: Vec<String>,
    #[serde(default)]
    pub occupancy_read_report_ids: Vec<String>,
}

/// One step in a gene's replayable decision-tree path.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct RegulatoryPartnerDecisionTraceStep {
    pub node_id: String,
    pub outcome: RegulatoryPartnerDecisionOutcome,
    pub detail: String,
    #[serde(default)]
    pub evidence_ids: Vec<String>,
}

/// One requested member and the evidence needed by the linked GUI views.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct RegulatoryPartnerGeneRow {
    pub member_dedup_key: String,
    pub gene_symbol: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_id: Option<String>,
    pub promoter_resolved: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub unresolved_reason: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub chromosome: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub strand: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub promoter_start_1based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub promoter_end_1based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tss_1based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tss_position_0based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub promoter_sequence_5to3: Option<String>,
    #[serde(default)]
    pub motif_hit_ids: Vec<String>,
    #[serde(default)]
    pub tuple_ids: Vec<String>,
    pub occupancy_evaluation_state: GeneSetCutRunEvaluationState,
    pub promoter_occupancy_supported: bool,
    pub terminal_class: String,
    #[serde(default)]
    pub trace: Vec<RegulatoryPartnerDecisionTraceStep>,
}

/// Exact evidence rows from which the decision tree can be replayed.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct RegulatoryPartnerTupleLedger {
    pub schema: String,
    pub request: RegulatoryPartnerScreenRequest,
    pub promoter_cohort: Box<GeneSetPromoterCohortReport>,
    #[serde(default)]
    pub source_reports: Vec<RegulatoryPartnerSourceReportRef>,
    pub requested_member_count: usize,
    pub resolved_promoter_count: usize,
    pub motif_hit_count: usize,
    pub tuple_count: usize,
    #[serde(default)]
    pub genes: Vec<RegulatoryPartnerGeneRow>,
    #[serde(default)]
    pub motif_hits: Vec<RegulatoryPartnerMotifHit>,
    #[serde(default)]
    pub tuples: Vec<RegulatoryPartnerTupleRow>,
    #[serde(default)]
    pub unresolved_members: Vec<GeneSetUnresolvedMember>,
}

/// One node in the fixed, human-readable evidence decision tree.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct RegulatoryPartnerDecisionNode {
    pub node_id: String,
    pub label: String,
    pub question: String,
    pub terminal: bool,
}

/// Directed edge for one pass/fail/unknown branch.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct RegulatoryPartnerDecisionEdge {
    pub from_node_id: String,
    pub to_node_id: String,
    pub outcome: RegulatoryPartnerDecisionOutcome,
}

/// Fixed decision tree and gene membership at every traversed node.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct RegulatoryPartnerDecisionTree {
    pub schema: String,
    pub root_node_id: String,
    #[serde(default)]
    pub nodes: Vec<RegulatoryPartnerDecisionNode>,
    #[serde(default)]
    pub edges: Vec<RegulatoryPartnerDecisionEdge>,
}

/// Regulatory-partner screen output for machine and linked GUI inspection.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct RegulatoryPartnerScreenReport {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub ledger: RegulatoryPartnerTupleLedger,
    pub decision_tree: RegulatoryPartnerDecisionTree,
    #[serde(default)]
    pub warnings: Vec<String>,
    #[serde(default)]
    pub non_claims: Vec<String>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn old_screen_payload_defaults_additive_trace_fields() {
        let report: RegulatoryPartnerScreenReport = serde_json::from_value(serde_json::json!({
            "schema": REGULATORY_PARTNER_SCREEN_SCHEMA,
            "ledger": {
                "schema": REGULATORY_PARTNER_TUPLE_LEDGER_SCHEMA,
                "request": {"genome_id": "GRCh38"},
                "genes": [{"member_dedup_key": "gene:TP73", "gene_symbol": "TP73"}]
            },
            "decision_tree": {"schema": REGULATORY_PARTNER_DECISION_TREE_SCHEMA}
        }))
        .expect("deserialize additive screen payload");
        assert_eq!(
            report.ledger.request.anchor_mode,
            RegulatoryPartnerAnchorMode::OccupancyPreferred
        );
        assert_eq!(
            report.ledger.genes[0].occupancy_evaluation_state,
            GeneSetCutRunEvaluationState::Unevaluated
        );
        assert!(report.ledger.genes[0].trace.is_empty());
    }
}
