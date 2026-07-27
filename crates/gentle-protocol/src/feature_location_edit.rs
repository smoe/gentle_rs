//! Portable contracts for previewing and applying simple feature-location edits.

use serde::{Deserialize, Serialize};

/// Schema emitted by feature-location preview and apply operations.
pub const FEATURE_LOCATION_EDIT_SCHEMA: &str = "gentle.feature_location_edit.v1";

/// Algorithm identifier for complete-feature optimistic-concurrency fingerprints.
pub const FEATURE_LOCATION_FINGERPRINT_ALGORITHM: &str =
    "sha256_gb_io_feature_serde_json_ordered_qualifiers_v1";

/// Strand represented by an exact simple GenBank feature location.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FeatureLocationEditStrand {
    Forward,
    Reverse,
}

/// Coordinate snapshot expressed in both engine-native and user-facing forms.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FeatureLocationSnapshot {
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub start_1based: usize,
    pub end_1based_inclusive: usize,
    pub strand: FeatureLocationEditStrand,
    pub five_prime_position_1based: usize,
    pub three_prime_position_1based: usize,
}

/// Request shared by preview and apply operations.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FeatureLocationEditRequest {
    pub seq_id: String,
    pub feature_index: usize,
    pub new_start_0based: i64,
    pub new_end_0based_exclusive: i64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub expected_feature_fingerprint_sha256: Option<String>,
}

/// Why another feature is worth reviewing after a boundary edit.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum RelatedFeatureBoundaryReason {
    SharesOldStartBoundary,
    SharesOldEndBoundary,
}

/// A nearby annotation that shares one of the edited feature's old boundaries.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RelatedFeatureBoundaryCandidate {
    pub feature_index: usize,
    pub feature_kind: String,
    pub location_display: String,
    pub bounding_start_0based: usize,
    pub bounding_end_0based_exclusive: usize,
    pub match_reasons: Vec<RelatedFeatureBoundaryReason>,
    pub modified: bool,
}

/// Deterministic preview or apply report.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FeatureLocationEditReport {
    pub schema: String,
    pub seq_id: String,
    pub feature_index: usize,
    pub feature_kind: String,
    pub dry_run: bool,
    pub applied: bool,
    pub before: FeatureLocationSnapshot,
    pub after: FeatureLocationSnapshot,
    pub fingerprint_algorithm: String,
    pub before_feature_fingerprint_sha256: String,
    pub after_feature_fingerprint_sha256: String,
    pub related_features: Vec<RelatedFeatureBoundaryCandidate>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn feature_location_edit_request_round_trips() {
        let request = FeatureLocationEditRequest {
            seq_id: "seq".to_string(),
            feature_index: 3,
            new_start_0based: 10,
            new_end_0based_exclusive: 40,
            expected_feature_fingerprint_sha256: Some("sha256:abc".to_string()),
        };
        let encoded = serde_json::to_string(&request).expect("request serializes");
        let decoded: FeatureLocationEditRequest =
            serde_json::from_str(&encoded).expect("request deserializes");
        assert_eq!(decoded, request);
    }
}
