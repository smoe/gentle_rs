//! Feature location/strand utilities used across render and engine code.

use gb_io::seq::{Feature, Location};

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
        false
    } else {
        strands.iter().filter(|is_reverse| **is_reverse).count() > strands.len() / 2
    }
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

pub fn feature_ranges_sorted_i64(feature: &Feature) -> Vec<(i64, i64)> {
    let mut ranges = Vec::new();
    collect_location_ranges_i64(&feature.location, &mut ranges);
    if ranges.is_empty() {
        if let Ok((from, to)) = feature.location.find_bounds() {
            if from >= 0 && to >= 0 {
                let mut start = from;
                let mut end = to;
                if end < start {
                    std::mem::swap(&mut start, &mut end);
                }
                ranges.push((start, end));
            }
        }
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
    use gb_io::seq::FeatureKind;

    fn make_feature(location: Location) -> Feature {
        Feature {
            kind: FeatureKind::from("mRNA"),
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
    fn unwraps_ranges_across_origin_in_order() {
        let wrapped = vec![(990, 999), (10, 30), (40, 50)];
        let unwrapped = unwrap_ranges_monotonic(1000, &wrapped);
        assert_eq!(unwrapped, vec![(990, 999), (1010, 1030), (1040, 1050)]);
    }
}
