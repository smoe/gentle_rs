//! Portable gene-group catalog contracts.

use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// Catalog schema for local, ontology-mappable gene groups.
pub const GENE_GROUP_CATALOG_SCHEMA: &str = "gentle.gene_group_catalog.v1";
/// Report schema for listing gene-group catalog entries.
pub const GENE_GROUP_LIST_REPORT_SCHEMA: &str = "gentle.gene_group_list.v1";
/// Report schema for showing one resolved gene-group entry.
pub const GENE_GROUP_SHOW_REPORT_SCHEMA: &str = "gentle.gene_group_show.v1";
/// Report schema for resolving a user token against gene-group entries.
pub const GENE_GROUP_RESOLVE_REPORT_SCHEMA: &str = "gentle.gene_group_resolve.v1";
/// Report schema for validating gene-group catalog overlays.
pub const GENE_GROUP_DOCTOR_REPORT_SCHEMA: &str = "gentle.gene_group_doctor.v1";
/// Report schema for creating a review-gated draft gene-group fragment.
pub const GENE_GROUP_DRAFT_REPORT_SCHEMA: &str = "gentle.gene_group_draft.v1";

/// External ontology/resource namespace that gene groups may map to.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GeneGroupExternalResource {
    pub id: String,
    pub label: String,
    pub namespace: String,
    pub role: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub homepage_url: Option<String>,
    #[serde(default)]
    pub notes: Vec<String>,
}

/// Mapping from a GENtle-owned group to an external ontology/resource term.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GeneGroupExternalMapping {
    pub namespace: String,
    pub id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub relationship: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub url: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub note: Option<String>,
}

/// One gene membership row inside a group definition.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GeneGroupMember {
    pub symbol: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_id: Option<String>,
    #[serde(default)]
    pub aliases: Vec<String>,
    #[serde(default)]
    pub external_ids: BTreeMap<String, Vec<String>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub evidence_note: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub confidence: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub status: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub provenance: Option<String>,
}

/// One deterministic local gene-group term.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GeneGroupRecord {
    pub id: String,
    pub label: String,
    #[serde(default)]
    pub aliases: Vec<String>,
    pub short_description: String,
    pub long_definition: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub organism: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub taxon_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub symbol_namespace: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reference_context: Option<String>,
    #[serde(default)]
    pub usages: Vec<String>,
    #[serde(default)]
    pub tags: Vec<String>,
    #[serde(default)]
    pub members: Vec<GeneGroupMember>,
    #[serde(default)]
    pub external_mappings: Vec<GeneGroupExternalMapping>,
    #[serde(default)]
    pub curation_status: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_kind: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub provenance: Option<String>,
    #[serde(default)]
    pub notes: Vec<String>,
}

/// Gene-group catalog file or fragment.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GeneGroupCatalog {
    pub schema: String,
    #[serde(default)]
    pub external_resources: Vec<GeneGroupExternalResource>,
    #[serde(default)]
    pub groups: Vec<GeneGroupRecord>,
}

/// Compact list row with source provenance.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GeneGroupListEntry {
    pub id: String,
    pub label: String,
    #[serde(default)]
    pub aliases: Vec<String>,
    pub short_description: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub organism: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub taxon_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub symbol_namespace: Option<String>,
    #[serde(default)]
    pub usages: Vec<String>,
    #[serde(default)]
    pub tags: Vec<String>,
    pub member_count: usize,
    pub curation_status: String,
    #[serde(default)]
    pub external_mappings: Vec<GeneGroupExternalMapping>,
    pub source_scope: String,
    pub source_path: String,
}

/// Source file/fragment diagnostics emitted by the doctor route.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GeneGroupCatalogSourceReport {
    pub scope: String,
    pub path: String,
    pub status: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sha1: Option<String>,
    pub group_count: usize,
    pub external_resource_count: usize,
    #[serde(default)]
    pub warnings: Vec<String>,
    #[serde(default)]
    pub errors: Vec<String>,
}

/// List report for gene-group entries.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GeneGroupListReport {
    pub schema: String,
    pub catalog_label: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub filter: Option<String>,
    pub returned_group_count: usize,
    #[serde(default)]
    pub groups: Vec<GeneGroupListEntry>,
    #[serde(default)]
    pub external_resources: Vec<GeneGroupExternalResource>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

/// Show report for one gene-group entry.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GeneGroupShowReport {
    pub schema: String,
    pub catalog_label: String,
    pub group: GeneGroupRecord,
    pub source_scope: String,
    pub source_path: String,
    #[serde(default)]
    pub external_resources: Vec<GeneGroupExternalResource>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

/// Resolve report for one user-facing gene-group token.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GeneGroupResolveReport {
    pub schema: String,
    pub catalog_label: String,
    pub query: String,
    pub normalized_query: String,
    pub matched_group_count: usize,
    #[serde(default)]
    pub groups: Vec<GeneGroupListEntry>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

/// Doctor report for catalog overlay validation.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GeneGroupDoctorReport {
    pub schema: String,
    pub catalog_label: String,
    pub source_count: usize,
    pub parsed_source_count: usize,
    pub group_count: usize,
    pub external_resource_count: usize,
    pub warning_count: usize,
    pub error_count: usize,
    #[serde(default)]
    pub sources: Vec<GeneGroupCatalogSourceReport>,
    #[serde(default)]
    pub warnings: Vec<String>,
    #[serde(default)]
    pub errors: Vec<String>,
}

/// Draft report for one generated catalog fragment.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GeneGroupDraftReport {
    pub schema: String,
    pub generation_method: String,
    pub review_required: bool,
    pub input_description_sha1: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub agent_provider: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub agent_model: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub agent_generated_at_utc: Option<String>,
    #[serde(default)]
    pub user_member_count: usize,
    #[serde(default)]
    pub candidate_member_count: usize,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub candidate_members: Vec<GeneGroupMember>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub unresolved_candidates: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub output_path: Option<String>,
    pub group: GeneGroupRecord,
    pub catalog_fragment: GeneGroupCatalog,
    #[serde(default)]
    pub warnings: Vec<String>,
}
