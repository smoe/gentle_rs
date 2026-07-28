//! Engine helpers for deterministic complete feature-record curation.

use gb_io::seq::{Feature, Location};
use gentle_protocol::{
    EngineError, FEATURE_ANNOTATION_STATE_FINGERPRINT_ALGORITHM, FeatureLocationEditStrand,
    FeatureRecordQualifier, FeatureRecordReviewCandidate, FeatureRecordReviewEvidence,
    FeatureRecordSnapshot,
};
use serde::Serialize;

use crate::{
    digest_utils::sha256_prefixed_str,
    feature_location::{
        feature_fingerprint_sha256, feature_ranges_sorted_i64, location_overlaps_usize,
    },
};

const SHARED_IDENTIFIER_QUALIFIERS: &[&str] = &["locus_tag", "gene", "protein_id", "transcript_id"];
const SHARED_IDENTIFIER_INTERPRETATION: &str = "potential_shared_identity_not_dependency";

#[derive(Serialize)]
struct AnnotationStateFingerprintInput<'a> {
    seq_id: &'a str,
    sequence_length: usize,
    circular: bool,
    features: &'a [Feature],
}

/// Fingerprint the complete ordered annotation table and its sequence context.
pub fn annotation_state_fingerprint_sha256(
    seq_id: &str,
    sequence_length: usize,
    circular: bool,
    features: &[Feature],
) -> Result<String, EngineError> {
    let input = AnnotationStateFingerprintInput {
        seq_id,
        sequence_length,
        circular,
        features,
    };
    let encoded = serde_json::to_string(&input).map_err(|err| {
        EngineError::internal(format!(
            "Could not serialize annotation state for {}",
            FEATURE_ANNOTATION_STATE_FINGERPRINT_ALGORITHM
        ))
        .with_cause(err)
    })?;
    Ok(sha256_prefixed_str(&encoded))
}

/// Build the one exact simple feature shape supported by create v1.
pub fn build_simple_feature(
    feature_kind: &str,
    start_0based: usize,
    end_0based_exclusive: usize,
    strand: FeatureLocationEditStrand,
    qualifiers: &[FeatureRecordQualifier],
) -> Result<Feature, EngineError> {
    if feature_kind.trim().is_empty() {
        return Err(EngineError::invalid_input("Feature kind must not be empty"));
    }
    if feature_kind.trim() != feature_kind {
        return Err(EngineError::invalid_input(
            "Feature kind must not have leading or trailing whitespace",
        ));
    }
    for qualifier in qualifiers {
        if qualifier.key.trim().is_empty() {
            return Err(EngineError::invalid_input(
                "Feature qualifier key must not be empty",
            ));
        }
        if qualifier.key.trim() != qualifier.key {
            return Err(EngineError::invalid_input(
                "Feature qualifier keys must not have leading or trailing whitespace",
            ));
        }
    }
    let location = Location::simple_range(
        i64::try_from(start_0based).map_err(|_| {
            EngineError::invalid_input(
                "Feature start exceeds the supported GenBank coordinate range",
            )
        })?,
        i64::try_from(end_0based_exclusive).map_err(|_| {
            EngineError::invalid_input("Feature end exceeds the supported GenBank coordinate range")
        })?,
    );
    Ok(Feature {
        kind: feature_kind.to_string().into(),
        location: match strand {
            FeatureLocationEditStrand::Forward => location,
            FeatureLocationEditStrand::Reverse => Location::Complement(Box::new(location)),
        },
        qualifiers: qualifiers
            .iter()
            .map(|qualifier| (qualifier.key.clone().into(), qualifier.value.clone()))
            .collect(),
    })
}

/// Convert one gb-io feature into a lossless portable report snapshot.
pub fn feature_record_snapshot(feature: &Feature) -> Result<FeatureRecordSnapshot, EngineError> {
    let location = serde_json::to_value(&feature.location).map_err(|err| {
        EngineError::internal("Could not serialize feature location for curation report")
            .with_cause(err)
    })?;
    let (bounding_start_0based, bounding_end_0based_exclusive) = feature
        .location
        .find_bounds()
        .ok()
        .and_then(|(start, end)| Some((usize::try_from(start).ok()?, usize::try_from(end).ok()?)))
        .map_or((None, None), |(start, end)| (Some(start), Some(end)));
    Ok(FeatureRecordSnapshot {
        feature_kind: feature.kind.to_string(),
        location,
        location_display: feature.location.to_string(),
        bounding_start_0based,
        bounding_end_0based_exclusive,
        qualifiers: feature
            .qualifiers
            .iter()
            .map(|(key, value)| FeatureRecordQualifier {
                key: key.to_string(),
                value: value.clone(),
            })
            .collect(),
        feature_fingerprint_sha256: feature_fingerprint_sha256(feature)?,
    })
}

fn features_overlap(left: &Feature, right: &Feature) -> bool {
    feature_ranges_sorted_i64(left)
        .into_iter()
        .filter_map(|(start, end)| Some((usize::try_from(start).ok()?, usize::try_from(end).ok()?)))
        .any(|(start, end)| {
            start < end && location_overlaps_usize(&right.location, start, end).unwrap_or(false)
        })
}

fn shared_identifier_evidence(
    target: &Feature,
    other: &Feature,
) -> Vec<FeatureRecordReviewEvidence> {
    let mut evidence = Vec::new();
    for supported_key in SHARED_IDENTIFIER_QUALIFIERS {
        let target_values = target
            .qualifiers
            .iter()
            .filter(|(key, _)| key.as_ref().eq_ignore_ascii_case(supported_key))
            .filter_map(|(_, value)| value.as_deref())
            .filter(|value| !value.is_empty());
        for target_value in target_values {
            let shared = other.qualifiers.iter().any(|(key, value)| {
                key.as_ref().eq_ignore_ascii_case(supported_key)
                    && value.as_deref() == Some(target_value)
            });
            if shared {
                let candidate = FeatureRecordReviewEvidence::SharedIdentifier {
                    qualifier_key: (*supported_key).to_string(),
                    qualifier_value: target_value.to_string(),
                    interpretation: SHARED_IDENTIFIER_INTERPRETATION.to_string(),
                };
                if !evidence.contains(&candidate) {
                    evidence.push(candidate);
                }
            }
        }
    }
    evidence
}

/// Find existing annotations worth inspecting before create/delete.
pub fn feature_record_review_candidates(
    features: &[Feature],
    target: &Feature,
    excluded_feature_index: Option<usize>,
) -> Vec<FeatureRecordReviewCandidate> {
    features
        .iter()
        .enumerate()
        .filter_map(|(feature_index, feature)| {
            if excluded_feature_index == Some(feature_index) {
                return None;
            }
            let mut evidence = Vec::new();
            if features_overlap(target, feature) {
                evidence.push(FeatureRecordReviewEvidence::OverlappingLocation);
            }
            evidence.extend(shared_identifier_evidence(target, feature));
            if evidence.is_empty() {
                return None;
            }
            Some(FeatureRecordReviewCandidate {
                feature_index,
                feature_kind: feature.kind.to_string(),
                location_display: feature.location.to_string(),
                evidence,
                modified: false,
            })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shared_identifier_is_informational_and_qualifiers_round_trip() {
        let target = Feature {
            kind: "gene".into(),
            location: Location::simple_range(10, 30),
            qualifiers: vec![
                ("gene".into(), Some("TEST".to_string())),
                ("pseudo".into(), None),
            ],
        };
        let other = Feature {
            kind: "mRNA".into(),
            location: Location::simple_range(20, 40),
            qualifiers: vec![("gene".into(), Some("TEST".to_string()))],
        };
        let candidates = feature_record_review_candidates(&[other], &target, None);
        assert_eq!(candidates.len(), 1);
        assert!(
            candidates[0]
                .evidence
                .contains(&FeatureRecordReviewEvidence::OverlappingLocation)
        );
        assert!(
            candidates[0]
                .evidence
                .contains(&FeatureRecordReviewEvidence::SharedIdentifier {
                    qualifier_key: "gene".to_string(),
                    qualifier_value: "TEST".to_string(),
                    interpretation: SHARED_IDENTIFIER_INTERPRETATION.to_string(),
                })
        );
        let snapshot = feature_record_snapshot(&target).expect("snapshot");
        assert_eq!(
            snapshot.qualifiers,
            vec![
                FeatureRecordQualifier {
                    key: "gene".to_string(),
                    value: Some("TEST".to_string()),
                },
                FeatureRecordQualifier {
                    key: "pseudo".to_string(),
                    value: None,
                },
            ]
        );
    }
}
