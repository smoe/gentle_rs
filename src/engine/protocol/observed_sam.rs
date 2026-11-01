//! S2.3 portable SAM evidence contracts, separate from S1 and original-read identity.

use super::*;
use crate::observed_transcript_structures as structures;

pub use structures::{
    ObservedSamAnalysis, ObservedSamAnnotationAvailability, ObservedSamCompanion,
    ObservedSamConversion, ObservedSamIntakeEvidence, ObservedSamMaskWitness,
    ObservedSamMissingCompanion, ObservedSamQueryUnit, ObservedSamQueryUnitKey,
    ObservedSamQueryUnitManifest, ObservedSamRowCounts, ObservedSamRowWitness,
    ObservedSamSearchDeclaration, ObservedSamTextBudgets, ObservedSamUnavailableReason,
    ObservedSamVersionPolicy,
};

pub const OBSERVED_SAM_REQUEST_SCHEMA: &str = "gentle.observed_sam_structures_request.v1";
pub const OBSERVED_SAM_REPORT_SCHEMA: &str = "gentle.observed_sam_structures_report.v1";
pub const OBSERVED_SAM_DATA_SCHEMA: &str = "gentle.observed_sam_structures_data.v1";
pub const OBSERVED_SAM_BINDING_SCHEMA: &str = "gentle.observed_sam_structures_result_binding.v1";
pub const OBSERVED_SAM_QUERY_UNIT_MANIFEST_SCHEMA: &str = "gentle.observed_sam_query_units.v1";

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum ObservedSamAnnotationBinding {
    Prepared {
        selector: structures::ReferenceGene,
        report: PrimerPayloadFile,
        annotation: Option<PrimerPayloadFile>,
    },
    /// Verified normalized geometry, not a fabricated preparation receipt.
    Normalized {
        selector: structures::ReferenceGene,
        provider: Option<String>,
        annotation: PrimerPayloadFile,
    },
}

impl ObservedSamAnnotationBinding {
    pub fn selector(&self) -> &structures::ReferenceGene {
        match self {
            Self::Prepared { selector, .. } | Self::Normalized { selector, .. } => selector,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedSamBudgets {
    pub max_input_bytes: u64,
    pub max_annotation_bytes: u64,
    pub max_output_bytes: u64,
    pub max_output_row_bytes: usize,
    pub max_units: usize,
    pub max_retained_bytes: usize,
    pub max_projection_work: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SummarizeObservedTranscriptStructuresFromSamRequest {
    pub schema: String,
    pub sam: PrimerPayloadFile,
    pub manifest: PrimerPayloadFile,
    pub dictionary: PrimerPayloadFile,
    pub conversion: ObservedSamConversion,
    pub reference_id: String,
    pub reference: structures::ReferenceBinding,
    pub reference_gene: structures::ReferenceGene,
    pub comparison_genes: Vec<structures::ReferenceGene>,
    pub annotation_sources: Vec<ObservedSamAnnotationBinding>,
    #[serde(default)]
    pub version_policy: ObservedSamVersionPolicy,
    pub unpadded_reference: bool,
    pub required_flank_bases: usize,
    pub text_budgets: ObservedSamTextBudgets,
    pub kernel_budgets: structures::Budgets,
    pub budgets: ObservedSamBudgets,
    pub output_dir: String,
    pub cancel_path: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedSamAnnotationEvidence {
    pub binding: ObservedSamAnnotationBinding,
    pub availability: ObservedSamAnnotationAvailability,
    pub preparation: Option<ObservedTranscriptAnnotationReport>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ObservedSamExecutionStatus {
    Completed,
    Cancelled,
    ControlError,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ObservedSamCoordinateAuthority {
    Declared,
    /// Reserved: this producer has no reference-byte verifier and cannot emit it.
    Verified,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedSamCounts {
    pub query_units: usize,
    pub assessed_units: usize,
    pub unavailable_units: usize,
    pub observations: usize,
    pub optimal_observations: usize,
    pub chain_families: usize,
    pub junctions: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedSamReport {
    pub schema: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub content_sha256: String,
    pub request_sha256: String,
    pub policy_sha256: String,
    pub semantic_context_sha256: Option<String>,
    pub request: SummarizeObservedTranscriptStructuresFromSamRequest,
    pub resolved_output_dir: String,
    pub status: ObservedSamExecutionStatus,
    pub control_error: Option<String>,
    pub source_eof_verified: bool,
    pub coordinate_authority: ObservedSamCoordinateAuthority,
    pub reference_verification: Option<PrimerPayloadFile>,
    pub identity_scope: String,
    pub counts: Option<ObservedSamCounts>,
    pub gene_resolutions: Vec<structures::ObservedSamSourceComparison>,
    pub artifacts: Vec<RnaReadPanelArtifact>,
    pub nonclaims: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedSamData {
    pub schema: String,
    pub request_sha256: String,
    pub intake: Option<ObservedSamIntakeEvidence>,
    pub annotations: Vec<ObservedSamAnnotationEvidence>,
    pub analysis: Option<ObservedSamAnalysis>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedSamResultBinding {
    pub schema: String,
    pub request_sha256: String,
    pub report_content_sha256: String,
    pub report_path: String,
    pub outcome: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedSamInspection {
    pub report: ObservedSamReport,
    pub artifacts_verified: bool,
    pub table: String,
    pub rows: Vec<serde_json::Value>,
}
