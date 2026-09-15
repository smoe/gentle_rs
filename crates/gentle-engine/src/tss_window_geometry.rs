//! Checked transcript-oriented windows and fixed-window unions/intersections.
//! No signal-driven boundary selection, sequence retrieval, or state mutation.

use gentle_protocol::{EngineError, tss_profiles::TssStrand, tss_window_geometry::*};
use std::collections::{BTreeMap, BTreeSet};

/// Positive, inclusive genomic bounds; minus-strand upstream is toward larger coordinates.
pub fn window_bounds(
    tss: u64,
    strand: TssStrand,
    upstream: usize,
    downstream: usize,
) -> Result<(u64, u64), EngineError> {
    let invalid = || {
        EngineError::invalid_input("TSS window crosses the contig start or overflows coordinates")
    };
    let (left, right) = match strand {
        TssStrand::Plus => (upstream, downstream),
        TssStrand::Minus => (downstream, upstream),
    };
    let start = tss
        .checked_sub(u64::try_from(left).map_err(|_| invalid())?)
        .filter(|&n| n > 0)
        .ok_or_else(invalid)?;
    let end = tss
        .checked_add(u64::try_from(right).map_err(|_| invalid())?)
        .ok_or_else(invalid)?;
    Ok((start, end))
}

fn require(condition: bool, message: &str) -> Result<(), EngineError> {
    if condition {
        Ok(())
    } else {
        Err(EngineError::invalid_input(message))
    }
}

fn identifiers<'a>(ids: impl IntoIterator<Item = &'a str>) -> Result<(), EngineError> {
    let mut seen = BTreeSet::new();
    for id in ids {
        require(
            !id.trim().is_empty() && id.len() <= 512 && seen.insert(id),
            "Empty, duplicate or oversized geometry identifier",
        )?;
    }
    Ok(())
}

/// Compute deterministic geometry only. All bounds are 1-based inclusive.
pub fn compute(
    mut request: TssWindowGeometryRequest,
) -> Result<TssWindowGeometryReport, EngineError> {
    require(
        request.schema == TSS_WINDOW_GEOMETRY_REQUEST_SCHEMA,
        "Unsupported TSS window geometry schema",
    )?;
    identifiers([request.assembly.as_str()])?;
    require(
        !request.groups.is_empty()
            && request.groups.len() <= 64
            && request.upstream_bp <= 100_000
            && request.downstream_bp <= 100_000,
        "Invalid TSS geometry request bounds",
    )?;
    identifiers(request.groups.iter().map(|g| g.group_id.as_str()))?;
    request.groups.sort_by(|a, b| a.group_id.cmp(&b.group_id));
    let mut groups = Vec::new();
    let mut work = 0usize;
    for group in &mut request.groups {
        identifiers([group.chromosome.as_str()])?;
        require(
            !group.anchors.is_empty()
                && group.anchors.len() <= 2048
                && group.features.len() <= 20_000,
            "Invalid TSS group size",
        )?;
        work = work.saturating_add(
            group
                .anchors
                .len()
                .saturating_mul(group.features.len().max(1)),
        );
        require(
            work <= 2_000_000,
            "TSS intersection evaluation budget exceeded",
        )?;
        identifiers(group.anchors.iter().map(|a| a.anchor_id.as_str()))?;
        identifiers(group.features.iter().map(|f| f.feature_id.as_str()))?;
        group.anchors.sort_by(|a, b| a.anchor_id.cmp(&b.anchor_id));
        group
            .features
            .sort_by(|a, b| a.feature_id.cmp(&b.feature_id));
        let mut windows = Vec::new();
        for anchor in &group.anchors {
            crate::tss_profiles::validate_geometry(&anchor.source)?;
            require(
                anchor.source.chromosome == group.chromosome
                    && anchor.source.strand == group.strand,
                "TSS group mixes chromosomes or strands",
            )?;
            let (start, end) = window_bounds(
                anchor.source.tss_1based,
                group.strand,
                request.upstream_bp,
                request.downstream_bp,
            )?;
            require(
                anchor.source.start_1based <= start && end <= anchor.source.end_1based,
                "Requested TSS window exceeds its prepared source",
            )?;
            windows.push(TssWindowGeometryWindow {
                anchor_id: anchor.anchor_id.clone(),
                start_1based: start,
                end_1based: end,
            });
        }
        windows.sort_by(|a, b| {
            (a.start_1based, a.end_1based, &a.anchor_id).cmp(&(
                b.start_1based,
                b.end_1based,
                &b.anchor_id,
            ))
        });
        let mut stretches: Vec<TssConnectedWindowStretch> = Vec::new();
        for window in &windows {
            if let Some(last) = stretches
                .last_mut()
                .filter(|last| window.start_1based <= last.end_1based.saturating_add(1))
            {
                last.end_1based = last.end_1based.max(window.end_1based);
                last.anchor_ids.push(window.anchor_id.clone());
            } else {
                stretches.push(TssConnectedWindowStretch {
                    start_1based: window.start_1based,
                    end_1based: window.end_1based,
                    anchor_ids: vec![window.anchor_id.clone()],
                });
            }
        }
        let sources = group
            .anchors
            .iter()
            .map(|a| (a.anchor_id.as_str(), &a.source))
            .collect::<BTreeMap<_, _>>();
        let mut intersections = Vec::new();
        for feature in &group.features {
            require(
                feature.start_1based > 0 && feature.end_1based >= feature.start_1based,
                "Invalid regulatory feature interval",
            )?;
            for (index, stretch) in stretches.iter().enumerate() {
                let start = feature.start_1based.max(stretch.start_1based);
                let end = feature.end_1based.min(stretch.end_1based);
                if start > end {
                    continue;
                }
                let mut containing_anchor_ids = stretch
                    .anchor_ids
                    .iter()
                    .filter(|id| {
                        sources[id.as_str()].start_1based <= start
                            && end <= sources[id.as_str()].end_1based
                    })
                    .cloned()
                    .collect::<Vec<_>>();
                containing_anchor_ids.sort();
                require(
                    !containing_anchor_ids.is_empty(),
                    "No single prepared source contains the feature intersection",
                )?;
                require(
                    intersections.len() < 100_000,
                    "TSS intersection output budget exceeded",
                )?;
                intersections.push(TssWindowFeatureIntersection {
                    stretch_index_1based: index + 1,
                    feature_id: feature.feature_id.clone(),
                    start_1based: start,
                    end_1based: end,
                    containing_anchor_ids,
                });
            }
        }
        intersections.sort_by(|a, b| {
            (a.stretch_index_1based, &a.feature_id).cmp(&(b.stretch_index_1based, &b.feature_id))
        });
        groups.push(TssWindowGeometryGroupResult {
            group_id: group.group_id.clone(),
            windows,
            stretches,
            intersections,
        });
    }
    let bytes =
        serde_json::to_vec(&request).map_err(|e| EngineError::invalid_input(e.to_string()))?;
    let hash = ring::digest::digest(&ring::digest::SHA256, &bytes);
    let request_sha256 = format!(
        "sha256:{}",
        hash.as_ref()
            .iter()
            .map(|b| format!("{b:02x}"))
            .collect::<String>()
    );
    Ok(TssWindowGeometryReport {
        schema: TSS_WINDOW_GEOMETRY_REPORT_SCHEMA.into(),
        request,
        request_sha256,
        groups,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use gentle_protocol::tss_profiles::TssGeometry;

    // Synthetic coordinates, no reference data or biological confidence implied.
    fn request(strand: TssStrand) -> TssWindowGeometryRequest {
        TssWindowGeometryRequest {
            schema: TSS_WINDOW_GEOMETRY_REQUEST_SCHEMA.into(),
            assembly: "synthetic".into(),
            upstream_bp: 500,
            downstream_bp: 200,
            groups: vec![TssWindowGeometryGroup {
                group_id: "TOY".into(),
                chromosome: "1".into(),
                strand,
                anchors: [1000, 1701, 2403]
                    .into_iter()
                    .map(|tss| TssWindowAnchor {
                        anchor_id: format!("tss{tss}"),
                        source: TssGeometry {
                            chromosome: "1".into(),
                            strand,
                            tss_1based: tss,
                            start_1based: tss - 800,
                            end_1based: tss + 800,
                            upstream_bp: 800,
                            downstream_bp: 800,
                        },
                    })
                    .collect(),
                features: vec![TssWindowFeature {
                    feature_id: "cross-gap".into(),
                    start_1based: 1901,
                    end_1based: 1903,
                }],
            }],
        }
    }

    #[test]
    fn tss_window_geometry_orients_checks_and_counts_the_tss_base() {
        assert_eq!(
            window_bounds(1000, TssStrand::Plus, 500, 200).unwrap(),
            (500, 1200)
        );
        assert_eq!(
            window_bounds(1000, TssStrand::Minus, 500, 200).unwrap(),
            (800, 1500)
        );
        let (a, b) = window_bounds(1000, TssStrand::Minus, 700, 300).unwrap();
        assert_eq!(b - a + 1, 1001);
        assert!(window_bounds(1, TssStrand::Plus, 1, 0).is_err());
        assert!(window_bounds(u64::MAX, TssStrand::Minus, 1, 0).is_err());
        assert_eq!(
            compute(request(TssStrand::Minus)).unwrap().groups[0].windows[0].start_1based,
            800
        );
    }

    #[test]
    fn tss_window_geometry_merges_touching_not_gaps_and_preserves_membership() {
        let r = compute(request(TssStrand::Plus)).unwrap();
        let g = &r.groups[0];
        assert_eq!(
            g.stretches
                .iter()
                .map(|s| (s.start_1based, s.end_1based))
                .collect::<Vec<_>>(),
            [(500, 1901), (1903, 2603)]
        );
        assert_eq!(g.stretches[0].anchor_ids, ["tss1000", "tss1701"]);
        assert_eq!(
            g.intersections
                .iter()
                .map(|i| (i.start_1based, i.end_1based))
                .collect::<Vec<_>>(),
            [(1901, 1901), (1903, 1903)]
        );
        assert_eq!(g.intersections[0].containing_anchor_ids, ["tss1701"]);
        let mut input = r.request.clone();
        input.groups[0].anchors.reverse();
        assert_eq!(r, compute(input).unwrap());
        assert_eq!(
            r,
            serde_json::from_slice(&serde_json::to_vec(&r).unwrap()).unwrap()
        );
    }

    #[test]
    fn tss_window_geometry_rejects_ambiguous_or_unavailable_sources() {
        for change in 0..6 {
            let mut r = request(TssStrand::Plus);
            match change {
                0 => r.groups[0].anchors[0].source.strand = TssStrand::Minus,
                1 => r.groups[0].anchors[0].source.chromosome = "2".into(),
                2 => r.groups[0].anchors[0].source.start_1based += 1,
                3 => r.groups[0].anchors[1].anchor_id = "tss1000".into(),
                4 => r.upstream_bp = 801,
                _ => r.groups[0].features[0].start_1based = 0,
            }
            assert!(compute(r).is_err(), "case {change}");
        }
        let mut r = request(TssStrand::Plus);
        r.groups[0].features = vec![TssWindowFeature {
            feature_id: "too-wide".into(),
            start_1based: 1,
            end_1based: 1901,
        }];
        assert!(compute(r).is_err());
    }

    #[test]
    fn tss_window_geometry_keeps_separate_groups_separate() {
        let mut r = request(TssStrand::Plus);
        let mut other = r.groups[0].clone();
        other.group_id = "OTHER".into();
        r.groups.push(other);
        assert_eq!(compute(r).unwrap().groups.len(), 2);
    }
}
