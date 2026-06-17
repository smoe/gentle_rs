//! Portable gene-set analysis contracts.
//!
//! Gene sets are resolved analysis operands: they expand curated catalog
//! groups, explicit user members, local ontology mappings, prepared-genome
//! neighborhoods, or deterministic random samples into auditable gene rows.

use serde::{Deserialize, Serialize};

/// Report schema for resolving a gene set into concrete member rows.
pub const GENE_SET_RESOLUTION_SCHEMA: &str = "gentle.gene_set_resolution.v1";
/// Report schema for promoter windows derived from a resolved gene set.
pub const GENE_SET_PROMOTER_COHORT_SCHEMA: &str = "gentle.gene_set_promoter_cohort.v1";
/// Report schema for CUT&RUN support summarized over a resolved gene set.
pub const GENE_SET_CUTRUN_REGULATORY_SUPPORT_SCHEMA: &str =
    "gentle.gene_set_cutrun_regulatory_support.v1";

/// Request source for resolving a gene set.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(tag = "source_kind", rename_all = "snake_case")]
pub enum GeneSetRequest {
    CatalogGroup {
        query: String,
    },
    ExplicitMembers {
        members: Vec<String>,
    },
    ExternalMapping {
        namespace: String,
        id: String,
    },
    GenomicNeighbors {
        anchor: String,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        flank_gene_count: Option<usize>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        flank_bp: Option<usize>,
        #[serde(default)]
        exclude_anchor: bool,
    },
    Random {
        count: usize,
        random_seed: u64,
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        exclude_members: Vec<String>,
    },
}

impl Default for GeneSetRequest {
    fn default() -> Self {
        Self::ExplicitMembers { members: vec![] }
    }
}

impl GeneSetRequest {
    pub fn source_kind_label(&self) -> &'static str {
        match self {
            Self::CatalogGroup { .. } => "catalog_group",
            Self::ExplicitMembers { .. } => "explicit_members",
            Self::ExternalMapping { .. } => "external_mapping",
            Self::GenomicNeighbors { .. } => "genomic_neighbors",
            Self::Random { .. } => "random",
        }
    }
}

/// One row of provenance attached to a resolved gene-set member or report.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneSetProvenanceRow {
    pub source_kind: String,
    pub source_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_path: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub note: Option<String>,
}

/// One gene resolved for a gene set.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneSetResolvedMember {
    pub dedup_key: String,
    pub symbol: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_id: Option<String>,
    #[serde(default)]
    pub aliases: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub chromosome: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub start_1based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub end_1based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub strand: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub biotype: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub member_status: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub confidence: Option<String>,
    #[serde(default)]
    pub contributing_group_ids: Vec<String>,
    #[serde(default)]
    pub provenance: Vec<GeneSetProvenanceRow>,
}

/// A requested member that could not be resolved.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneSetUnresolvedMember {
    pub query: String,
    pub reason: String,
    pub source_kind: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_id: Option<String>,
}

/// Random-source accounting that makes deterministic sampling auditable.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneSetRandomProvenance {
    pub genome_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub genome_build: Option<String>,
    pub gene_index_source: String,
    pub random_seed: u64,
    pub universe_size: usize,
    pub foreground_exclusion_count: usize,
}

/// Resolved, provenance-rich gene-set report.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneSetResolutionReport {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub request: GeneSetRequest,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub genome_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_group_catalog_label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub genome_catalog_label: Option<String>,
    #[serde(default)]
    pub contributing_group_ids: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub random: Option<GeneSetRandomProvenance>,
    pub requested_member_count: usize,
    pub resolved_member_count: usize,
    pub unresolved_member_count: usize,
    #[serde(default)]
    pub resolved_members: Vec<GeneSetResolvedMember>,
    #[serde(default)]
    pub unresolved_members: Vec<GeneSetUnresolvedMember>,
    #[serde(default)]
    pub warnings: Vec<String>,
    #[serde(default)]
    pub provenance: Vec<GeneSetProvenanceRow>,
}

/// One promoter window derived for a resolved gene-set member.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneSetPromoterWindow {
    pub member_dedup_key: String,
    pub symbol: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_id: Option<String>,
    pub gene_query: String,
    pub occurrence: usize,
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
    pub sequence_orientation: String,
    pub used_fuzzy_gene_match: bool,
}

/// Promoter-window cohort derived from a gene-set resolution.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneSetPromoterCohortReport {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub genome_id: String,
    pub upstream_bp: usize,
    pub downstream_bp: usize,
    pub gene_set_resolution: GeneSetResolutionReport,
    pub requested_member_count: usize,
    pub returned_window_count: usize,
    #[serde(default)]
    pub windows: Vec<GeneSetPromoterWindow>,
    #[serde(default)]
    pub unresolved_members: Vec<GeneSetUnresolvedMember>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

/// Evaluation state for one gene-set member in CUT&RUN support analysis.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum GeneSetCutRunEvaluationState {
    #[default]
    Unevaluated,
    Evaluated,
}

/// CUT&RUN support summary for one gene-set member.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GeneSetCutRunMemberSupport {
    pub member_dedup_key: String,
    pub symbol: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_id: Option<String>,
    pub evaluation_state: GeneSetCutRunEvaluationState,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub unevaluated_reason: Option<String>,
    pub promoter_start_1based: usize,
    pub promoter_end_1based: usize,
    pub chromosome: String,
    pub support_window_count: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub strongest_support_strength: Option<String>,
    pub overlapping_peak_count: usize,
    pub supporting_fragment_count: usize,
    pub cut_site_count: u32,
    #[serde(default)]
    pub contributing_dataset_ids: Vec<String>,
    #[serde(default)]
    pub contributing_read_report_ids: Vec<String>,
}

/// Aggregate statistics for gene-set CUT&RUN support.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GeneSetCutRunSupportAggregate {
    pub member_count: usize,
    pub evaluated_member_count: usize,
    pub unevaluated_member_count: usize,
    pub members_with_support_windows: usize,
    pub members_with_strong_support: usize,
    pub evaluated_fraction_with_support_windows: f64,
    pub evaluated_fraction_with_strong_support: f64,
    pub mean_support_window_count_evaluated: f64,
}

/// CUT&RUN regulatory-support report over a promoter gene-set cohort.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GeneSetCutRunRegulatorySupportReport {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub genome_id: String,
    pub promoter_cohort: GeneSetPromoterCohortReport,
    #[serde(default)]
    pub dataset_ids: Vec<String>,
    #[serde(default)]
    pub read_report_ids: Vec<String>,
    pub aggregate: GeneSetCutRunSupportAggregate,
    #[serde(default)]
    pub member_support: Vec<GeneSetCutRunMemberSupport>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gene_set_request_serializes_source_kind() {
        let request = GeneSetRequest::ExternalMapping {
            namespace: "GO".to_string(),
            id: "GO:0000381".to_string(),
        };
        let value = serde_json::to_value(&request).expect("serialize request");
        assert_eq!(value["source_kind"], "external_mapping");
        assert_eq!(value["id"], "GO:0000381");
    }

    #[test]
    fn resolution_report_defaults_to_schema_ready_shape() {
        let report = GeneSetResolutionReport {
            schema: GENE_SET_RESOLUTION_SCHEMA.to_string(),
            request: GeneSetRequest::CatalogGroup {
                query: "p53_family".to_string(),
            },
            ..GeneSetResolutionReport::default()
        };
        let json = serde_json::to_string(&report).expect("serialize report");
        assert!(json.contains("gentle.gene_set_resolution.v1"));
        assert!(json.contains("catalog_group"));
    }
}
