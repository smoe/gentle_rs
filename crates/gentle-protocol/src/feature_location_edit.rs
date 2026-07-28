//! Portable contracts for previewing and applying feature-location edits.

use serde::{Deserialize, Serialize};

/// Schema emitted by feature-location preview and apply operations.
pub const FEATURE_LOCATION_EDIT_SCHEMA: &str = "gentle.feature_location_edit.v1";

/// Schema emitted when one segment of an exact flat compound is edited.
pub const FEATURE_LOCATION_EDIT_SCHEMA_V2: &str = "gentle.feature_location_edit.v2";

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
    /// Zero-based child index for an exact flat `Join` or `Order`.
    ///
    /// `None` retains the v1 whole-feature behavior for one exact simple range.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub segment_index: Option<usize>,
}

/// The portion of a feature changed by a v2 location edit.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FeatureLocationEditTargetScope {
    Segment,
}

/// Flat compound operator retained by a segment edit.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FeatureLocationCompoundKind {
    Join,
    Order,
}

/// Established genomic direction of the compound's stored child list.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FeatureLocationStoredDirection {
    Ascending,
    Descending,
}

/// Context that identifies one edited child without duplicating its coordinates.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FeatureLocationCompoundContext {
    pub container_kind: FeatureLocationCompoundKind,
    pub stored_direction: FeatureLocationStoredDirection,
    pub segment_index: usize,
    pub biological_segment_number: usize,
    pub total_segments: usize,
}

/// Non-blocking biological or geometric review note for a segment edit.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FeatureLocationCompoundWarning {
    OverlappingSegments {
        edited_segment_index: usize,
        overlapping_segment_index: usize,
    },
    EstablishedDirectionBroken {
        first_segment_index: usize,
        second_segment_index: usize,
        established_direction: FeatureLocationStoredDirection,
    },
    CdsCodingLengthDeltaNotDivisibleByThree {
        signed_length_delta_nt: i64,
    },
}

/// One side of a half-open genomic interval.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FeatureLocationIntervalBoundaryRole {
    Start0Based,
    End0BasedExclusive,
}

/// Another exact feature segment sharing one old interbase boundary.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RelatedSegmentBoundaryCandidate {
    pub related_feature_index: usize,
    pub related_feature_kind: String,
    pub related_location_display: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub related_segment_index: Option<usize>,
    pub edited_boundary: FeatureLocationIntervalBoundaryRole,
    pub related_boundary: FeatureLocationIntervalBoundaryRole,
    pub boundary_coordinate_0based: usize,
    pub modified: bool,
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
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub target_scope: Option<FeatureLocationEditTargetScope>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub compound_context: Option<FeatureLocationCompoundContext>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub compound_validation_warnings: Vec<FeatureLocationCompoundWarning>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub related_segment_boundaries: Vec<RelatedSegmentBoundaryCandidate>,
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
            segment_index: Some(2),
        };
        let encoded = serde_json::to_string(&request).expect("request serializes");
        let decoded: FeatureLocationEditRequest =
            serde_json::from_str(&encoded).expect("request deserializes");
        assert_eq!(decoded, request);
    }

    #[test]
    fn simple_v1_request_omits_segment_field() {
        let request = FeatureLocationEditRequest {
            seq_id: "seq".to_string(),
            feature_index: 3,
            new_start_0based: 10,
            new_end_0based_exclusive: 40,
            expected_feature_fingerprint_sha256: None,
            segment_index: None,
        };
        let encoded = serde_json::to_value(request).expect("request serializes");
        assert_eq!(
            encoded
                .as_object()
                .expect("object")
                .keys()
                .cloned()
                .collect::<Vec<_>>(),
            vec![
                "feature_index",
                "new_end_0based_exclusive",
                "new_start_0based",
                "seq_id",
            ]
        );
        let decoded: FeatureLocationEditRequest =
            serde_json::from_value(encoded).expect("v1 request deserializes");
        assert_eq!(decoded.segment_index, None);
    }

    #[test]
    fn v1_report_deserializes_without_v2_fields() {
        let report: FeatureLocationEditReport = serde_json::from_value(serde_json::json!({
            "schema": FEATURE_LOCATION_EDIT_SCHEMA,
            "seq_id": "seq",
            "feature_index": 1,
            "feature_kind": "gene",
            "dry_run": true,
            "applied": false,
            "before": {
                "start_0based": 10,
                "end_0based_exclusive": 20,
                "start_1based": 11,
                "end_1based_inclusive": 20,
                "strand": "forward",
                "five_prime_position_1based": 11,
                "three_prime_position_1based": 20
            },
            "after": {
                "start_0based": 11,
                "end_0based_exclusive": 21,
                "start_1based": 12,
                "end_1based_inclusive": 21,
                "strand": "forward",
                "five_prime_position_1based": 12,
                "three_prime_position_1based": 21
            },
            "fingerprint_algorithm": FEATURE_LOCATION_FINGERPRINT_ALGORITHM,
            "before_feature_fingerprint_sha256": "sha256:before",
            "after_feature_fingerprint_sha256": "sha256:after",
            "related_features": []
        }))
        .expect("v1 report deserializes");
        assert_eq!(report.target_scope, None);
        assert_eq!(report.compound_context, None);
        assert!(report.compound_validation_warnings.is_empty());
        assert!(report.related_segment_boundaries.is_empty());
    }
}
