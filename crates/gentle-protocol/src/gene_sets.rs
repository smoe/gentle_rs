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
/// Local offline cache schema for direct gene-list producer inputs.
pub const GENE_SET_DIRECT_LIST_CACHE_SCHEMA: &str = "gentle.gene_set_direct_list_cache.v1";
/// Local offline cache schema for ontology-assignment producer inputs.
pub const GENE_SET_ONTOLOGY_ASSIGNMENT_CACHE_SCHEMA: &str =
    "gentle.gene_set_ontology_assignment_cache.v1";
/// Local offline cache schema for co-regulated cohort producer inputs.
pub const GENE_SET_CO_REGULATED_CACHE_SCHEMA: &str = "gentle.gene_set_co_regulated_cache.v1";

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

/// Retrieval producer family that supplied candidate members before resolution.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum GeneSetProducerKind {
    #[default]
    DirectGeneList,
    OntologyAssignment,
    CoRegulatedCohort,
}

/// Review state for a resolved set that came from a retrieval producer.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum GeneSetResolutionReviewStatus {
    #[default]
    Unreviewed,
    Reviewed,
    Included,
    Draft,
    Deprecated,
}

/// User-declared regulatory expectation among members of a resolved cohort.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum GeneSetCohortRelationship {
    #[default]
    Unspecified,
    Manual,
    CoRegulated,
    AntiCoRegulated,
}

/// Non-blocking relationship-expectation flag derived from available evidence.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneSetCohortRelationshipFlag {
    pub flag_kind: String,
    pub evidence_kind: String,
    #[serde(default)]
    pub member_symbols: Vec<String>,
    #[serde(default)]
    pub member_dedup_keys: Vec<String>,
    pub detail: String,
}

/// Report-level provenance for a retrieval producer or imported cache.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneSetProducerProvenance {
    pub producer_kind: GeneSetProducerKind,
    pub provider_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub provider_label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub provider_version: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cache_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cache_path: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cache_version: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cache_digest: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub import_op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub imported_at_unix_ms: Option<u128>,
}

/// Structured filter used by a retrieval producer.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneSetProducerFilter {
    pub field: String,
    pub operator: String,
    pub value: String,
}

/// Structured query metadata for a gene-set retrieval producer.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneSetProducerQueryMetadata {
    pub query_kind: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub query_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub query_label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub organism: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub taxon_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub symbol_namespace: Option<String>,
    #[serde(default)]
    pub filters: Vec<GeneSetProducerFilter>,
}

/// Retrieval metadata for an evidence-derived co-regulated cohort.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(default)]
pub struct GeneSetCoRegulatedProducerMetadata {
    #[serde(default)]
    pub dataset_ids: Vec<String>,
    #[serde(default)]
    pub contrast_labels: Vec<String>,
    #[serde(default)]
    pub condition_labels: Vec<String>,
    pub normalization_method: String,
    pub scoring_method: String,
    pub threshold_rule: String,
    pub sign_direction_rule: String,
    pub relationship: GeneSetCohortRelationship,
    pub interpretation_note: String,
}

impl Default for GeneSetCoRegulatedProducerMetadata {
    fn default() -> Self {
        Self {
            dataset_ids: Vec::new(),
            contrast_labels: Vec::new(),
            condition_labels: Vec::new(),
            normalization_method: String::new(),
            scoring_method: String::new(),
            threshold_rule: String::new(),
            sign_direction_rule: String::new(),
            relationship: GeneSetCohortRelationship::Unspecified,
            interpretation_note:
                "This evidence-derived cohort is a retrieval result and does not prove regulation."
                    .to_string(),
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
    pub review_status: GeneSetResolutionReviewStatus,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub genome_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub organism: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub taxon_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub symbol_namespace: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_group_catalog_label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub genome_catalog_label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub producer: Option<GeneSetProducerProvenance>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub query_metadata: Option<GeneSetProducerQueryMetadata>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub co_regulated_metadata: Option<GeneSetCoRegulatedProducerMetadata>,
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
    pub relationship: GeneSetCohortRelationship,
    pub gene_set_resolution: GeneSetResolutionReport,
    pub requested_member_count: usize,
    pub returned_window_count: usize,
    #[serde(default)]
    pub windows: Vec<GeneSetPromoterWindow>,
    #[serde(default)]
    pub relationship_flags: Vec<GeneSetCohortRelationshipFlag>,
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
    pub relationship: GeneSetCohortRelationship,
    #[serde(default)]
    pub dataset_ids: Vec<String>,
    #[serde(default)]
    pub read_report_ids: Vec<String>,
    pub aggregate: GeneSetCutRunSupportAggregate,
    #[serde(default)]
    pub member_support: Vec<GeneSetCutRunMemberSupport>,
    #[serde(default)]
    pub relationship_flags: Vec<GeneSetCohortRelationshipFlag>,
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

    #[test]
    fn old_resolution_payload_defaults_producer_metadata() {
        let report: GeneSetResolutionReport = serde_json::from_value(serde_json::json!({
            "schema": GENE_SET_RESOLUTION_SCHEMA,
            "generated_at_unix_ms": 1,
            "request": {"source_kind": "explicit_members", "members": ["TP73"]},
            "requested_member_count": 1,
            "resolved_member_count": 0,
            "unresolved_member_count": 1
        }))
        .expect("deserialize old gene-set resolution payload");

        assert_eq!(
            report.review_status,
            GeneSetResolutionReviewStatus::Unreviewed
        );
        assert!(report.organism.is_none());
        assert!(report.taxon_id.is_none());
        assert!(report.symbol_namespace.is_none());
        assert!(report.producer.is_none());
        assert!(report.query_metadata.is_none());
        assert!(report.co_regulated_metadata.is_none());
    }

    #[test]
    fn resolution_report_serializes_retrieval_producer_metadata() {
        let report = GeneSetResolutionReport {
            schema: GENE_SET_RESOLUTION_SCHEMA.to_string(),
            request: GeneSetRequest::ExplicitMembers {
                members: vec!["TP73".to_string(), "TP53".to_string()],
            },
            review_status: GeneSetResolutionReviewStatus::Reviewed,
            organism: Some("Homo sapiens".to_string()),
            taxon_id: Some("9606".to_string()),
            symbol_namespace: Some("HGNC".to_string()),
            producer: Some(GeneSetProducerProvenance {
                producer_kind: GeneSetProducerKind::OntologyAssignment,
                provider_id: "local_go_cache".to_string(),
                provider_label: Some("Local GO assignment cache".to_string()),
                provider_version: Some("2026-06".to_string()),
                cache_id: Some("go-cache-v1".to_string()),
                cache_path: Some("resources/go_assignments.json".to_string()),
                cache_version: Some("1".to_string()),
                cache_digest: Some("sha256:abc123".to_string()),
                import_op_id: Some("op_import_go_cache".to_string()),
                imported_at_unix_ms: Some(42),
            }),
            query_metadata: Some(GeneSetProducerQueryMetadata {
                query_kind: "ontology_assignment".to_string(),
                query_id: Some("GO:0000381".to_string()),
                query_label: Some("regulation of alternative mRNA splicing".to_string()),
                organism: Some("Homo sapiens".to_string()),
                taxon_id: Some("9606".to_string()),
                symbol_namespace: Some("HGNC".to_string()),
                filters: vec![GeneSetProducerFilter {
                    field: "evidence_code".to_string(),
                    operator: "equals".to_string(),
                    value: "IDA".to_string(),
                }],
            }),
            co_regulated_metadata: Some(GeneSetCoRegulatedProducerMetadata {
                dataset_ids: vec!["rnaseq_toy".to_string()],
                contrast_labels: vec!["treated_vs_control".to_string()],
                condition_labels: vec!["treated".to_string(), "control".to_string()],
                normalization_method: "log2_tpm".to_string(),
                scoring_method: "signed_delta".to_string(),
                threshold_rule: "abs(delta) >= 1.0".to_string(),
                sign_direction_rule: "same_sign".to_string(),
                relationship: GeneSetCohortRelationship::CoRegulated,
                ..GeneSetCoRegulatedProducerMetadata::default()
            }),
            ..GeneSetResolutionReport::default()
        };

        let value = serde_json::to_value(&report).expect("serialize gene-set producer metadata");
        assert_eq!(value["review_status"], "reviewed");
        assert_eq!(value["organism"], "Homo sapiens");
        assert_eq!(value["producer"]["producer_kind"], "ontology_assignment");
        assert_eq!(value["producer"]["provider_id"], "local_go_cache");
        assert_eq!(
            value["query_metadata"]["filters"][0]["field"],
            "evidence_code"
        );
        assert_eq!(
            value["co_regulated_metadata"]["relationship"],
            "co_regulated"
        );
        assert_eq!(
            value["co_regulated_metadata"]["interpretation_note"],
            "This evidence-derived cohort is a retrieval result and does not prove regulation."
        );

        let round_trip: GeneSetResolutionReport =
            serde_json::from_value(value).expect("round-trip gene-set producer metadata");
        assert_eq!(round_trip, report);
    }

    #[test]
    fn old_promoter_cohort_payload_defaults_relationship_to_unspecified() {
        let report: GeneSetPromoterCohortReport = serde_json::from_value(serde_json::json!({
            "schema": GENE_SET_PROMOTER_COHORT_SCHEMA,
            "generated_at_unix_ms": 1,
            "genome_id": "ToyGenome",
            "upstream_bp": 100,
            "downstream_bp": 20,
            "gene_set_resolution": {
                "schema": GENE_SET_RESOLUTION_SCHEMA,
                "generated_at_unix_ms": 1,
                "request": {"source_kind": "explicit_members", "members": ["A"]}
            },
            "requested_member_count": 1,
            "returned_window_count": 0
        }))
        .expect("deserialize old gene-set promoter cohort");
        assert_eq!(report.relationship, GeneSetCohortRelationship::Unspecified);
        assert!(report.relationship_flags.is_empty());
    }

    #[test]
    fn old_cutrun_payload_defaults_relationship_to_unspecified() {
        let report: GeneSetCutRunRegulatorySupportReport =
            serde_json::from_value(serde_json::json!({
                "schema": GENE_SET_CUTRUN_REGULATORY_SUPPORT_SCHEMA,
                "generated_at_unix_ms": 1,
                "genome_id": "ToyGenome",
                "promoter_cohort": {
                    "schema": GENE_SET_PROMOTER_COHORT_SCHEMA,
                    "generated_at_unix_ms": 1,
                    "genome_id": "ToyGenome",
                    "upstream_bp": 100,
                    "downstream_bp": 20,
                    "gene_set_resolution": {
                        "schema": GENE_SET_RESOLUTION_SCHEMA,
                        "generated_at_unix_ms": 1,
                        "request": {"source_kind": "explicit_members", "members": ["A"]}
                    },
                    "requested_member_count": 1,
                    "returned_window_count": 0
                }
            }))
            .expect("deserialize old gene-set CUT&RUN support report");
        assert_eq!(report.relationship, GeneSetCohortRelationship::Unspecified);
        assert!(report.relationship_flags.is_empty());
        assert_eq!(
            report.promoter_cohort.relationship,
            GeneSetCohortRelationship::Unspecified
        );
        assert!(report.promoter_cohort.relationship_flags.is_empty());
    }
}
