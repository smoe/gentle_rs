//! Reference identities and independent primer-specificity selection receipts.

use super::*;

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default, deny_unknown_fields)]
pub struct PrimerSpecificitySourceReference {
    pub genome_id: String,
    pub assembly: Option<String>,
    pub release: Option<String>,
    pub sequence_sha1: Option<String>,
    pub annotation_sha1: Option<String>,
}

/// Prefixes and display labels are not biological identities.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default, deny_unknown_fields)]
pub struct PrimerSpecificityReferenceIdentity {
    pub genome_id: String,
    pub index_kind: crate::genomes::BlastDatabaseIndexKind,
    pub assembly: Option<String>,
    pub release: Option<String>,
    pub content_fingerprint: String,
    pub fingerprint_algorithm: String,
    pub subject_annotation_fingerprint: Option<String>,
    pub subject_annotation_fingerprint_algorithm: Option<String>,
}

/// Companion evidence; it does not widen the panel's single-reference gates.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct PrimerSpecificityReferenceSelection {
    pub schema: String,
    pub selection_id: String,
    pub panel_report_id: String,
    pub panel_digest: String,
    pub reference: PrimerSpecificityReferenceIdentity,
    pub policy: PrimerSpecificityPolicy,
    pub report_ids_by_assay: BTreeMap<String, String>,
    pub report_content_sha256_by_assay: BTreeMap<String, String>,
    pub pair_bindings_by_assay: BTreeMap<String, String>,
    pub acceptance: TranscriptAssayPanelSpecificityAcceptance,
    pub replaces_selection_id: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PrimerSpecificityPanelSource {
    pub panel_design_digest: String,
    pub source_reference: PrimerSpecificitySourceReference,
}
