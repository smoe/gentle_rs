//! Feature location/strand utilities used across render and engine code.

use crate::digest_utils::sha256_prefixed_str;
use gb_io::seq::{Feature, Location};
use gentle_protocol::{
    EngineError, ErrorCode, FEATURE_LOCATION_FINGERPRINT_ALGORITHM, FeatureLocationEditStrand,
    FeatureLocationSnapshot, RelatedFeatureBoundaryCandidate, RelatedFeatureBoundaryReason,
};

/// Return the exact editable interval represented by a simple feature location.
///
/// The editor deliberately accepts only an exact `Range` or one
/// `Complement(Range)`. Compound, fuzzy, external, gap, and between-base
/// locations require a richer biological editing contract.
pub fn exact_feature_location_snapshot(
    feature: &Feature,
) -> Result<FeatureLocationSnapshot, EngineError> {
    let (start, end, strand) = match &feature.location {
        Location::Range((start, before), (end, after)) if !before.0 && !after.0 => {
            (*start, *end, FeatureLocationEditStrand::Forward)
        }
        Location::Complement(inner) => match inner.as_ref() {
            Location::Range((start, before), (end, after)) if !before.0 && !after.0 => {
                (*start, *end, FeatureLocationEditStrand::Reverse)
            }
            _ => {
                return Err(EngineError::new(
                    ErrorCode::Unsupported,
                    "Feature location editing supports only exact simple Range or Complement(Range) locations",
                ));
            }
        },
        _ => {
            return Err(EngineError::new(
                ErrorCode::Unsupported,
                "Feature location editing supports only exact simple Range or Complement(Range) locations",
            ));
        }
    };
    let start = usize::try_from(start).map_err(|_| {
        EngineError::invalid_input("Feature location start cannot be represented as usize")
    })?;
    let end = usize::try_from(end).map_err(|_| {
        EngineError::invalid_input("Feature location end cannot be represented as usize")
    })?;
    feature_location_snapshot(start, end, strand)
}

/// Build a coordinate snapshot after validating a non-empty half-open range.
pub fn feature_location_snapshot(
    start_0based: usize,
    end_0based_exclusive: usize,
    strand: FeatureLocationEditStrand,
) -> Result<FeatureLocationSnapshot, EngineError> {
    if start_0based >= end_0based_exclusive {
        return Err(EngineError::invalid_input(
            "Feature location must have start_0based < end_0based_exclusive",
        ));
    }
    let start_1based = start_0based.checked_add(1).ok_or_else(|| {
        EngineError::invalid_input("Feature location start exceeds supported coordinate range")
    })?;
    let end_1based_inclusive = end_0based_exclusive;
    let (five_prime_position_1based, three_prime_position_1based) = match strand {
        FeatureLocationEditStrand::Forward => (start_1based, end_1based_inclusive),
        FeatureLocationEditStrand::Reverse => (end_1based_inclusive, start_1based),
    };
    Ok(FeatureLocationSnapshot {
        start_0based,
        end_0based_exclusive,
        start_1based,
        end_1based_inclusive,
        strand,
        five_prime_position_1based,
        three_prime_position_1based,
    })
}

/// Construct a simple GenBank location while preserving the existing strand.
pub fn simple_location_from_snapshot(
    snapshot: &FeatureLocationSnapshot,
) -> Result<Location, EngineError> {
    let start = i64::try_from(snapshot.start_0based).map_err(|_| {
        EngineError::invalid_input("Feature location start exceeds the GenBank coordinate range")
    })?;
    let end = i64::try_from(snapshot.end_0based_exclusive).map_err(|_| {
        EngineError::invalid_input("Feature location end exceeds the GenBank coordinate range")
    })?;
    let location = Location::simple_range(start, end);
    Ok(match snapshot.strand {
        FeatureLocationEditStrand::Forward => location,
        FeatureLocationEditStrand::Reverse => Location::Complement(Box::new(location)),
    })
}

/// Fingerprint the complete feature record, including qualifier order.
pub fn feature_fingerprint_sha256(feature: &Feature) -> Result<String, EngineError> {
    let encoded = serde_json::to_string(feature).map_err(|err| {
        EngineError::internal(format!(
            "Could not serialize feature for {} fingerprint",
            FEATURE_LOCATION_FINGERPRINT_ALGORITHM
        ))
        .with_cause(err)
    })?;
    Ok(sha256_prefixed_str(&encoded))
}

/// Find annotations sharing either old boundary with the feature being edited.
pub fn related_feature_boundary_candidates(
    features: &[Feature],
    edited_feature_index: usize,
    before: &FeatureLocationSnapshot,
) -> Vec<RelatedFeatureBoundaryCandidate> {
    features
        .iter()
        .enumerate()
        .filter_map(|(feature_index, feature)| {
            if feature_index == edited_feature_index {
                return None;
            }
            let (start, end) = feature.location.find_bounds().ok()?;
            let (start, end) = (usize::try_from(start).ok()?, usize::try_from(end).ok()?);
            let mut match_reasons = Vec::new();
            if start == before.start_0based {
                match_reasons.push(RelatedFeatureBoundaryReason::SharesOldStartBoundary);
            }
            if end == before.end_0based_exclusive {
                match_reasons.push(RelatedFeatureBoundaryReason::SharesOldEndBoundary);
            }
            if match_reasons.is_empty() {
                return None;
            }
            Some(RelatedFeatureBoundaryCandidate {
                feature_index,
                feature_kind: feature.kind.to_string(),
                location_display: feature.location.to_string(),
                bounding_start_0based: start,
                bounding_end_0based_exclusive: end,
                match_reasons,
                modified: false,
            })
        })
        .collect()
}

pub fn collect_location_strands(location: &Location, reverse: bool, strands: &mut Vec<bool>) {
    match location {
        Location::Range(_, _) | Location::Between(_, _) => strands.push(reverse),
        Location::Complement(inner) => collect_location_strands(inner, !reverse, strands),
        Location::Join(parts)
        | Location::Order(parts)
        | Location::Bond(parts)
        | Location::OneOf(parts) => {
            for part in parts {
                collect_location_strands(part, reverse, strands);
            }
        }
        Location::External(_, maybe_loc) => {
            if let Some(loc) = maybe_loc {
                collect_location_strands(loc, reverse, strands);
            }
        }
        Location::Gap(_) => {}
    }
}

pub fn feature_is_reverse(feature: &Feature) -> bool {
    let mut strands = Vec::new();
    collect_location_strands(&feature.location, false, &mut strands);
    if strands.is_empty() {
        return feature_strand_qualifier_is_reverse(feature).unwrap_or(false);
    }
    if location_contains_complement(&feature.location) {
        strands.iter().filter(|is_reverse| **is_reverse).count() > strands.len() / 2
    } else {
        feature_strand_qualifier_is_reverse(feature).unwrap_or(false)
    }
}

fn location_contains_complement(location: &Location) -> bool {
    match location {
        Location::Complement(_) => true,
        Location::Join(parts)
        | Location::Order(parts)
        | Location::Bond(parts)
        | Location::OneOf(parts) => parts.iter().any(location_contains_complement),
        Location::External(_, maybe_loc) => maybe_loc
            .as_deref()
            .is_some_and(location_contains_complement),
        Location::Range(_, _) | Location::Between(_, _) | Location::Gap(_) => false,
    }
}

fn feature_strand_qualifier_is_reverse(feature: &Feature) -> Option<bool> {
    for value in feature.qualifier_values("strand") {
        let normalized = value.trim().to_ascii_lowercase();
        match normalized.as_str() {
            "-" | "-1" | "minus" | "reverse" => return Some(true),
            "+" | "+1" | "1" | "plus" | "forward" => return Some(false),
            _ => {}
        }
    }
    None
}

pub fn collect_location_ranges_i64(location: &Location, ranges: &mut Vec<(i64, i64)>) {
    match location {
        Location::Range((from, _), (to, _)) | Location::Between(from, to) => {
            if *from < 0 || *to < 0 {
                return;
            }
            let mut start = *from;
            let mut end = *to;
            if end < start {
                std::mem::swap(&mut start, &mut end);
            }
            ranges.push((start, end));
        }
        Location::Complement(inner) => collect_location_ranges_i64(inner, ranges),
        Location::Join(parts)
        | Location::Order(parts)
        | Location::Bond(parts)
        | Location::OneOf(parts) => {
            for part in parts {
                collect_location_ranges_i64(part, ranges);
            }
        }
        Location::External(_, maybe_loc) => {
            if let Some(loc) = maybe_loc {
                collect_location_ranges_i64(loc, ranges);
            }
        }
        Location::Gap(_) => {}
    }
}

pub fn collect_location_ranges_usize(location: &Location, ranges: &mut Vec<(usize, usize)>) {
    let mut ranges_i64 = Vec::new();
    collect_location_ranges_i64(location, &mut ranges_i64);
    ranges.extend(ranges_i64.into_iter().filter_map(|(start, end)| {
        Some((usize::try_from(start).ok()?, usize::try_from(end).ok()?))
    }));
}

/// Test whether any concrete location segment overlaps a half-open range.
///
/// `None` means the location contained no usable concrete segment, allowing a
/// caller to fall back to `find_bounds` without allocating a temporary vector.
pub fn location_overlaps_usize(
    location: &Location,
    range_start: usize,
    range_end_exclusive: usize,
) -> Option<bool> {
    match location {
        Location::Range((from, _), (to, _)) | Location::Between(from, to) => {
            let (Ok(mut start), Ok(mut end)) = (usize::try_from(*from), usize::try_from(*to))
            else {
                return None;
            };
            if end < start {
                std::mem::swap(&mut start, &mut end);
            }
            Some(start < range_end_exclusive && end > range_start)
        }
        Location::Complement(inner) => {
            location_overlaps_usize(inner, range_start, range_end_exclusive)
        }
        Location::Join(parts)
        | Location::Order(parts)
        | Location::Bond(parts)
        | Location::OneOf(parts) => {
            let mut found_segment = false;
            for part in parts {
                if let Some(overlaps) =
                    location_overlaps_usize(part, range_start, range_end_exclusive)
                {
                    found_segment = true;
                    if overlaps {
                        return Some(true);
                    }
                }
            }
            found_segment.then_some(false)
        }
        Location::External(_, maybe_loc) => maybe_loc
            .as_deref()
            .and_then(|inner| location_overlaps_usize(inner, range_start, range_end_exclusive)),
        Location::Gap(_) => None,
    }
}

pub fn feature_ranges_sorted_i64(feature: &Feature) -> Vec<(i64, i64)> {
    let mut ranges = Vec::new();
    collect_location_ranges_i64(&feature.location, &mut ranges);
    if ranges.is_empty()
        && let Ok((from, to)) = feature.location.find_bounds()
        && from >= 0
        && to >= 0
    {
        let mut start = from;
        let mut end = to;
        if end < start {
            std::mem::swap(&mut start, &mut end);
        }
        ranges.push((start, end));
    }
    ranges.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    ranges
}

pub fn normalize_range(seq_len: i64, from: i64, to: i64) -> Option<(i64, i64)> {
    if seq_len <= 0 {
        return None;
    }
    let mut start = from.rem_euclid(seq_len);
    let mut end = to.rem_euclid(seq_len);
    if end < start {
        std::mem::swap(&mut start, &mut end);
    }
    Some((start, end))
}

pub fn unwrap_ranges_monotonic(seq_len: i64, ranges: &[(i64, i64)]) -> Vec<(i64, i64)> {
    if seq_len <= 0 {
        return vec![];
    }
    let mut out = Vec::new();
    let mut previous_end: Option<i64> = None;
    for (raw_start, raw_end) in ranges {
        let Some((mut start, mut end)) = normalize_range(seq_len, *raw_start, *raw_end) else {
            continue;
        };
        if let Some(prev_end) = previous_end {
            while start <= prev_end {
                start += seq_len;
                end += seq_len;
            }
        }
        if end < start {
            continue;
        }
        out.push((start, end));
        previous_end = Some(end);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_feature(location: Location) -> Feature {
        Feature {
            kind: "mRNA".into(),
            location,
            qualifiers: vec![],
        }
    }

    #[test]
    fn parses_join_ranges_and_reverse_strand() {
        let feature = make_feature(Location::Complement(Box::new(Location::Join(vec![
            Location::simple_range(10, 20),
            Location::simple_range(40, 50),
        ]))));
        let ranges = feature_ranges_sorted_i64(&feature);
        assert_eq!(ranges, vec![(10, 20), (40, 50)]);
        assert!(feature_is_reverse(&feature));
    }

    #[test]
    fn overlap_query_handles_join_without_collecting_ranges() {
        let location = Location::Join(vec![
            Location::simple_range(10, 20),
            Location::simple_range(40, 50),
        ]);
        assert_eq!(location_overlaps_usize(&location, 15, 16), Some(true));
        assert_eq!(location_overlaps_usize(&location, 25, 30), Some(false));
    }

    #[test]
    fn complete_feature_fingerprint_survives_json_round_trip() {
        let feature = Feature {
            kind: "gene".into(),
            location: Location::Complement(Box::new(Location::simple_range(10, 30))),
            qualifiers: vec![
                ("gene".into(), Some("TEST".to_string())),
                ("note".into(), Some("ordered".to_string())),
            ],
        };
        let before = feature_fingerprint_sha256(&feature).expect("fingerprint");
        assert_eq!(
            before,
            "sha256:3f73e11799f0ae94d5ef9494781254bad52f2a8cdb96517b899571659a4c725b"
        );
        let encoded = serde_json::to_string(&feature).expect("serialize");
        let decoded: Feature = serde_json::from_str(&encoded).expect("deserialize");
        assert_eq!(
            feature_fingerprint_sha256(&decoded).expect("fingerprint"),
            before
        );
    }

    #[test]
    fn exact_feature_location_rejects_every_unsupported_location_shape() {
        let exact = Location::simple_range(10, 30);
        let unsupported = vec![
            Location::Between(10, 11),
            Location::Join(vec![exact.clone()]),
            Location::Order(vec![exact.clone()]),
            Location::Bond(vec![exact.clone()]),
            Location::OneOf(vec![exact.clone()]),
            Location::External("REMOTE".to_string(), Some(Box::new(exact.clone()))),
            Location::Gap(gb_io::seq::GapLength::Known(20)),
            Location::Complement(Box::new(Location::Join(vec![exact.clone()]))),
            Location::Range(
                (10, gb_io::seq::Before(true)),
                (30, gb_io::seq::After(false)),
            ),
            Location::Range(
                (10, gb_io::seq::Before(false)),
                (30, gb_io::seq::After(true)),
            ),
        ];
        for location in unsupported {
            let feature = Feature {
                kind: "misc_feature".into(),
                location,
                qualifiers: vec![],
            };
            let error =
                exact_feature_location_snapshot(&feature).expect_err("location must be rejected");
            assert_eq!(error.code, ErrorCode::Unsupported);
        }
    }

    #[test]
    fn unwraps_ranges_across_origin_in_order() {
        let wrapped = vec![(990, 999), (10, 30), (40, 50)];
        let unwrapped = unwrap_ranges_monotonic(1000, &wrapped);
        assert_eq!(unwrapped, vec![(990, 999), (1010, 1030), (1040, 1050)]);
    }

    #[test]
    fn falls_back_to_strand_qualifier_when_location_lacks_complement() {
        let feature = Feature {
            kind: "gene".into(),
            location: Location::simple_range(10, 20),
            qualifiers: vec![("strand".into(), Some("-".to_string()))],
        };
        assert!(feature_is_reverse(&feature));
    }

    #[test]
    fn explicit_complement_location_wins_over_conflicting_qualifier() {
        let feature = Feature {
            kind: "gene".into(),
            location: Location::Complement(Box::new(Location::simple_range(10, 20))),
            qualifiers: vec![("strand".into(), Some("+".to_string()))],
        };
        assert!(feature_is_reverse(&feature));
    }
}
