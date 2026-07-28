//! Portable contracts for previewing and applying feature-record creation and deletion.

use serde::{Deserialize, Serialize};
use serde_json::Value;

use crate::FeatureLocationEditStrand;

/// Schema emitted by feature-record curation preview and apply operations.
pub const FEATURE_RECORD_CURATION_SCHEMA: &str = "gentle.feature_record_curation.v1";

/// Algorithm identifier for sequence annotation-state fingerprints.
pub const FEATURE_ANNOTATION_STATE_FINGERPRINT_ALGORITHM: &str =
    "sha256_sequence_id_length_topology_ordered_gb_io_features_serde_json_v1";

/// One ordered GenBank qualifier, including an explicitly valueless qualifier.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FeatureRecordQualifier {
    pub key: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub value: Option<String>,
}

/// Create one exact simple feature at the end of the ordered feature table.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FeatureRecordCreateRequest {
    pub seq_id: String,
    pub feature_kind: String,
    pub start_0based: i64,
    pub end_0based_exclusive: i64,
    pub strand: FeatureLocationEditStrand,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub qualifiers: Vec<FeatureRecordQualifier>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub expected_annotation_state_fingerprint_sha256: Option<String>,
}

/// Delete one complete feature record without constraining its location shape.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FeatureRecordDeleteRequest {
    pub seq_id: String,
    pub feature_index: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub expected_feature_fingerprint_sha256: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub expected_annotation_state_fingerprint_sha256: Option<String>,
}

/// One tagged feature-record curation request.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(
    deny_unknown_fields,
    tag = "operation_kind",
    content = "request",
    rename_all = "snake_case"
)]
pub enum FeatureRecordCurationRequest {
    Create(FeatureRecordCreateRequest),
    Delete(FeatureRecordDeleteRequest),
}

impl FeatureRecordCurationRequest {
    /// Sequence whose ordered feature table is being curated.
    pub fn seq_id(&self) -> &str {
        match self {
            Self::Create(request) => &request.seq_id,
            Self::Delete(request) => &request.seq_id,
        }
    }
}

/// Kind of complete feature-record curation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FeatureRecordCurationKind {
    Create,
    Delete,
}

/// Portable snapshot of a complete feature record.
///
/// `location` is the lossless serde representation supplied by `gb-io`;
/// `location_display` is included for human inspection.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FeatureRecordSnapshot {
    pub feature_kind: String,
    pub location: Value,
    pub location_display: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub bounding_start_0based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub bounding_end_0based_exclusive: Option<usize>,
    pub qualifiers: Vec<FeatureRecordQualifier>,
    pub feature_fingerprint_sha256: String,
}

/// Why another annotation is worth reviewing alongside the create/delete action.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "evidence_kind", rename_all = "snake_case")]
pub enum FeatureRecordReviewEvidence {
    OverlappingLocation,
    SharedIdentifier {
        qualifier_key: String,
        qualifier_value: String,
        interpretation: String,
    },
}

/// Existing feature with informational spatial or identifier evidence.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FeatureRecordReviewCandidate {
    pub feature_index: usize,
    pub feature_kind: String,
    pub location_display: String,
    pub evidence: Vec<FeatureRecordReviewEvidence>,
    pub modified: bool,
}

/// Operation-specific result details.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "outcome_kind", rename_all = "snake_case")]
pub enum FeatureRecordCurationOutcome {
    Create {
        proposed_feature: FeatureRecordSnapshot,
        /// Present only after apply; preview deliberately does not promise an index.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        created_feature_index: Option<usize>,
    },
    Delete {
        deleted_feature_index: usize,
        deleted_feature: FeatureRecordSnapshot,
        /// Records the index renumbering caused by removal from the ordered table.
        shifted_feature_count: usize,
    },
}

/// Deterministic preview or apply report for one complete feature-record action.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FeatureRecordCurationReport {
    pub schema: String,
    pub seq_id: String,
    pub operation_kind: FeatureRecordCurationKind,
    pub dry_run: bool,
    pub applied: bool,
    pub fingerprint_algorithm: String,
    pub before_annotation_state_fingerprint_sha256: String,
    pub after_annotation_state_fingerprint_sha256: String,
    pub outcome: FeatureRecordCurationOutcome,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub review_candidates: Vec<FeatureRecordReviewCandidate>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn create_request_preserves_qualifier_order_duplicates_and_no_value() {
        let request = FeatureRecordCurationRequest::Create(FeatureRecordCreateRequest {
            seq_id: "seq".to_string(),
            feature_kind: "misc_feature".to_string(),
            start_0based: 10,
            end_0based_exclusive: 20,
            strand: FeatureLocationEditStrand::Forward,
            qualifiers: vec![
                FeatureRecordQualifier {
                    key: "note".to_string(),
                    value: Some("first".to_string()),
                },
                FeatureRecordQualifier {
                    key: "pseudo".to_string(),
                    value: None,
                },
                FeatureRecordQualifier {
                    key: "note".to_string(),
                    value: Some("second".to_string()),
                },
            ],
            expected_annotation_state_fingerprint_sha256: None,
        });
        let encoded = serde_json::to_string(&request).expect("request serializes");
        let decoded: FeatureRecordCurationRequest =
            serde_json::from_str(&encoded).expect("request deserializes");
        assert_eq!(decoded, request);
        assert_eq!(
            serde_json::from_str::<Value>(&encoded).expect("JSON")["operation_kind"],
            "create"
        );
    }

    #[test]
    fn delete_request_round_trips_without_preview_tokens() {
        let request = FeatureRecordCurationRequest::Delete(FeatureRecordDeleteRequest {
            seq_id: "seq".to_string(),
            feature_index: 4,
            expected_feature_fingerprint_sha256: None,
            expected_annotation_state_fingerprint_sha256: None,
        });
        let encoded = serde_json::to_string(&request).expect("request serializes");
        let decoded: FeatureRecordCurationRequest =
            serde_json::from_str(&encoded).expect("request deserializes");
        assert_eq!(decoded, request);
        assert_eq!(
            serde_json::from_str::<Value>(&encoded).expect("JSON")["operation_kind"],
            "delete"
        );
    }
}
