//! Reference identities, selection receipts and standalone multi-reference handoffs.

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

pub const PRIMER_SPECIFICITY_MULTI_REQUEST_SCHEMA: &str =
    "gentle.primer_pair_multi_reference_request.v1";
pub const PRIMER_SPECIFICITY_MULTI_HANDOFF_SCHEMA: &str =
    "gentle.primer_pair_multi_reference_handoff.v1";
pub const PRIMER_SPECIFICITY_MULTI_MANIFEST_SCHEMA: &str =
    "gentle.primer_pair_multi_reference_execution_manifest.v1";
pub const PRIMER_SPECIFICITY_MULTI_SUMMARY_SCHEMA: &str =
    "gentle.primer_pair_multi_reference_summary.v1";

/// Exactly one saved pair or two explicit full oligos with declared tail boundaries.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum PrimerSpecificityMultiPair {
    SavedPair {
        primer_report_id: String,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        pair_rank: Option<usize>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        pair_index: Option<usize>,
    },
    ExplicitPair {
        forward: PrimerSpecificityInputPrimer,
        reverse: PrimerSpecificityInputPrimer,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrimerSpecificityMultiReference {
    pub target_genome_id: String,
    pub expected_index_kind: crate::genomes::BlastDatabaseIndexKind,
    pub required: bool,
    /// Explicit geometry is caller-provided, not an inferred liftover.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub intended_target: Option<PrimerSpecificityIntendedTarget>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrimerSpecificityMultiRequest {
    pub schema: String,
    pub pair: PrimerSpecificityMultiPair,
    #[serde(default)]
    pub policy: PrimerSpecificityPolicy,
    pub references: Vec<PrimerSpecificityMultiReference>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub catalog_path: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cache_dir: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum PrimerSpecificityMultiAvailability {
    Prepared,
    Unavailable,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrimerSpecificityMultiHandoffReference {
    pub requested_genome_id: String,
    pub resolved_genome_id: String,
    pub required: bool,
    pub expected_index_kind: crate::genomes::BlastDatabaseIndexKind,
    pub availability: PrimerSpecificityMultiAvailability,
    pub diagnostic: Option<String>,
    pub reference: Option<PrimerSpecificityReferenceIdentity>,
    pub child: Option<PrimerSpecificityHandoff>,
}

/// Files and commands are bound; this is not proof that a process executed.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrimerSpecificityMultiHandoff {
    pub schema: String,
    pub handoff_id: String,
    pub content_sha256: String,
    pub request: PrimerSpecificityMultiRequest,
    pub pair_binding_sha256: String,
    pub source_snapshot_sha256: Option<String>,
    pub primers: Vec<PrimerSpecificityInputPrimer>,
    pub references: Vec<PrimerSpecificityMultiHandoffReference>,
    pub manifest_path: String,
    pub nonclaims: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum PrimerSpecificityMultiExecutionState {
    Pending,
    Completed,
    Failed,
    Cancelled,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrimerSpecificityMultiExecutionCommand {
    pub command_id: String,
    pub command_sha256: String,
    pub output_path: String,
    pub state: PrimerSpecificityMultiExecutionState,
    pub exit_code: Option<i32>,
    pub output_size_bytes: Option<u64>,
    pub output_sha256: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrimerSpecificityMultiExecutionManifest {
    pub schema: String,
    pub handoff_id: String,
    pub handoff_content_sha256: String,
    pub pair_binding_sha256: String,
    pub commands: Vec<PrimerSpecificityMultiExecutionCommand>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum PrimerSpecificityMultiVerdict {
    Pass,
    Fail,
    Incomplete,
    NotRequested,
    NotRequired,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrimerSpecificityMultiSummaryReference {
    pub genome_id: String,
    pub required: bool,
    pub index_kind: crate::genomes::BlastDatabaseIndexKind,
    pub reference: Option<PrimerSpecificityReferenceIdentity>,
    pub verdict: PrimerSpecificityMultiVerdict,
    /// Applicability was checked at import; display does not probe resources.
    pub applicability: String,
    pub diagnostics: Vec<String>,
    pub report: Option<PrimerSpecificityReport>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrimerSpecificityMultiSummary {
    pub schema: String,
    pub summary_id: String,
    pub content_sha256: String,
    pub handoff_id: String,
    pub pair_binding_sha256: String,
    pub source_snapshot_sha256: Option<String>,
    pub request: PrimerSpecificityMultiRequest,
    pub execution_manifest: PrimerSpecificityMultiExecutionManifest,
    pub execution_complete: bool,
    pub references: Vec<PrimerSpecificityMultiSummaryReference>,
    pub genomic: PrimerSpecificityMultiVerdict,
    pub transcriptome: PrimerSpecificityMultiVerdict,
    pub nonclaims: Vec<String>,
}
