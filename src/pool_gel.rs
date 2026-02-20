use crate::DNA_LADDERS;
use std::collections::{BTreeMap, BTreeSet};
use svg::node::element::{Line, Rectangle, Text};
use svg::Document;

const SVG_WIDTH: f32 = 1040.0;
const SVG_HEIGHT: f32 = 760.0;
const GEL_LEFT: f32 = 90.0;
const GEL_RIGHT: f32 = SVG_WIDTH - 270.0;
const GEL_TOP: f32 = 90.0;
const GEL_BOTTOM: f32 = SVG_HEIGHT - 110.0;

#[derive(Clone, Debug)]
pub struct PoolGelBand {
    pub bp: usize,
    pub intensity: f32,
    pub count: usize,
    pub labels: Vec<String>,
}

#[derive(Clone, Debug)]
pub struct PoolGelLane {
    pub name: String,
    pub is_ladder: bool,
    pub bands: Vec<PoolGelBand>,
}

#[derive(Clone, Debug)]
pub struct PoolGelLayout {
    pub lanes: Vec<PoolGelLane>,
    pub selected_ladders: Vec<String>,
    pub sample_count: usize,
    pub pool_member_count: usize,
    pub range_min_bp: usize,
    pub range_max_bp: usize,
}

#[derive(Clone, Debug)]
pub struct GelSampleInput {
    pub name: String,
    pub members: Vec<(String, usize)>,
}

impl PoolGelLayout {
    pub fn y_for_bp(&self, bp: usize, top: f32, bottom: f32) -> f32 {
        let min_bp = self.range_min_bp.max(1) as f64;
        let max_bp = self.range_max_bp.max(self.range_min_bp + 1) as f64;
        let bp = bp.clamp(self.range_min_bp.max(1), self.range_max_bp.max(2)) as f64;
        let log_min = min_bp.log10();
        let log_max = max_bp.log10();
        let denom = (log_max - log_min).max(1e-6);
        let f = ((log_max - bp.log10()) / denom) as f32;
        top + f.clamp(0.0, 1.0) * (bottom - top)
    }
}

fn normalize_ladder_name(name: &str) -> String {
    name.trim().to_ascii_lowercase()
}

fn resolve_ladder_names(requested: &[String], min_bp: usize, max_bp: usize) -> Vec<String> {
    let names = DNA_LADDERS.names_sorted();
    if names.is_empty() {
        return vec![];
    }

    if !requested.is_empty() {
        let mut picked: Vec<String> = vec![];
        for requested_name in requested.iter().take(2) {
            let needle = normalize_ladder_name(requested_name);
            if needle.is_empty() {
                continue;
            }
            if let Some(found) = names
                .iter()
                .find(|name| normalize_ladder_name(name) == needle)
                .cloned()
            {
                if !picked.iter().any(|n| n == &found) {
                    picked.push(found);
                }
                continue;
            }
            if let Some(found) = names
                .iter()
                .find(|name| normalize_ladder_name(name).contains(&needle))
                .cloned()
            {
                if !picked.iter().any(|n| n == &found) {
                    picked.push(found);
                }
            }
        }
        if !picked.is_empty() {
            return picked;
        }
    }

    DNA_LADDERS.choose_for_range(min_bp, max_bp, 2)
}

fn sample_bands(members: &[(String, usize)]) -> Vec<PoolGelBand> {
    let mut by_bp: BTreeMap<usize, Vec<String>> = BTreeMap::new();
    for (id, bp) in members {
        by_bp.entry(*bp).or_default().push(id.clone());
    }
    by_bp
        .iter()
        .rev()
        .map(|(bp, ids)| {
            let count = ids.len();
            let intensity = (0.42 + (count as f32 * 0.2)).clamp(0.3, 1.0);
            PoolGelBand {
                bp: *bp,
                intensity,
                count,
                labels: ids.clone(),
            }
        })
        .collect::<Vec<_>>()
}

pub fn build_serial_gel_layout(
    samples: &[GelSampleInput],
    requested_ladders: &[String],
) -> Result<PoolGelLayout, String> {
    if samples.is_empty() {
        return Err("Serial gel needs at least one sample lane".to_string());
    }
    let mut normalized_samples: Vec<GelSampleInput> = vec![];
    let mut all_members: Vec<(String, usize)> = vec![];
    for (sample_idx, sample) in samples.iter().enumerate() {
        let mut valid_members = sample
            .members
            .iter()
            .filter(|(_, bp)| *bp > 0)
            .cloned()
            .collect::<Vec<_>>();
        if valid_members.is_empty() {
            return Err(format!(
                "Sample lane {} has no sequence lengths > 0 bp",
                sample_idx + 1
            ));
        }
        valid_members.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
        all_members.extend(valid_members.clone());
        let lane_name = if sample.name.trim().is_empty() {
            format!("Sample {} (n={})", sample_idx + 1, valid_members.len())
        } else {
            sample.name.clone()
        };
        normalized_samples.push(GelSampleInput {
            name: lane_name,
            members: valid_members,
        });
    }

    let pool_min = all_members.iter().map(|(_, bp)| *bp).min().unwrap_or(1);
    let pool_max = all_members.iter().map(|(_, bp)| *bp).max().unwrap_or(pool_min);
    let selected_ladders = resolve_ladder_names(requested_ladders, pool_min, pool_max);
    if selected_ladders.is_empty() {
        return Err("No DNA ladders available for pool-gel rendering".to_string());
    }

    let mut lanes: Vec<PoolGelLane> = vec![];
    let mut all_band_bps: Vec<usize> = all_members.iter().map(|(_, bp)| *bp).collect();

    for ladder_name in &selected_ladders {
        let Some(ladder) = DNA_LADDERS.get(ladder_name) else {
            continue;
        };
        let mut parts = ladder.bands().clone();
        parts.sort_by(|a, b| b.length_bp().total_cmp(&a.length_bp()));
        let max_strength = parts
            .iter()
            .filter_map(|p| p.relative_strength)
            .fold(0.0_f64, f64::max)
            .max(1.0);
        let bands = parts
            .into_iter()
            .map(|p| {
                let bp = p.length_bp().round().max(1.0) as usize;
                all_band_bps.push(bp);
                let raw = p.relative_strength.unwrap_or(1.0).max(0.1);
                let intensity = (raw / max_strength).clamp(0.18, 1.0) as f32;
                PoolGelBand {
                    bp,
                    intensity,
                    count: 1,
                    labels: vec![format!("{bp} bp")],
                }
            })
            .collect::<Vec<_>>();
        lanes.push(PoolGelLane {
            name: ladder.name().to_string(),
            is_ladder: true,
            bands,
        });
    }

    if lanes.is_empty() {
        return Err("None of the selected DNA ladders could be resolved".to_string());
    }

    for sample in &normalized_samples {
        lanes.push(PoolGelLane {
            name: sample.name.clone(),
            is_ladder: false,
            bands: sample_bands(&sample.members),
        });
    }

    let min_band = all_band_bps.iter().copied().min().unwrap_or(pool_min);
    let max_band = all_band_bps.iter().copied().max().unwrap_or(pool_max);
    let range_min_bp = ((min_band as f64) * 0.72).floor().max(1.0) as usize;
    let mut range_max_bp = ((max_band as f64) * 1.30).ceil().max(2.0) as usize;
    if range_max_bp <= range_min_bp {
        range_max_bp = range_min_bp + 1;
    }

    Ok(PoolGelLayout {
        lanes,
        selected_ladders,
        sample_count: normalized_samples.len(),
        pool_member_count: all_members.len(),
        range_min_bp,
        range_max_bp,
    })
}

pub fn build_pool_gel_layout(
    pool_members: &[(String, usize)],
    requested_ladders: &[String],
) -> Result<PoolGelLayout, String> {
    let members = pool_members
        .iter()
        .filter(|(_, bp)| *bp > 0)
        .cloned()
        .collect::<Vec<_>>();
    if members.is_empty() {
        return Err("Pool gel needs sequence lengths > 0 bp".to_string());
    }
    let sample = GelSampleInput {
        name: format!("Pool (n={})", members.len()),
        members,
    };
    build_serial_gel_layout(&[sample], requested_ladders)
}

pub fn export_pool_gel_svg(layout: &PoolGelLayout) -> String {
    let lane_count = layout.lanes.len().max(1);
    let lane_gap = (GEL_RIGHT - GEL_LEFT) / (lane_count as f32 + 1.0);
    let gel_width = GEL_RIGHT - GEL_LEFT;
    let gel_height = GEL_BOTTOM - GEL_TOP;
    let sample_lane_indices = layout
        .lanes
        .iter()
        .enumerate()
        .filter_map(|(idx, lane)| (!lane.is_ladder).then_some(idx))
        .collect::<Vec<_>>();

    let mut doc = Document::new()
        .set("viewBox", (0, 0, SVG_WIDTH, SVG_HEIGHT))
        .set("width", SVG_WIDTH)
        .set("height", SVG_HEIGHT)
        .add(
            Rectangle::new()
                .set("x", 0)
                .set("y", 0)
                .set("width", SVG_WIDTH)
                .set("height", SVG_HEIGHT)
                .set("fill", "#f9fafb"),
        )
        .add(
            Rectangle::new()
                .set("x", GEL_LEFT)
                .set("y", GEL_TOP)
                .set("width", gel_width)
                .set("height", gel_height)
                .set("rx", 10)
                .set("ry", 10)
                .set("fill", "#111315"),
        );

    let mut tick_bps = BTreeSet::new();
    for lane in layout.lanes.iter().filter(|l| l.is_ladder) {
        for band in &lane.bands {
            tick_bps.insert(band.bp);
        }
    }
    if tick_bps.is_empty() {
        tick_bps.insert(layout.range_min_bp);
        tick_bps.insert(layout.range_max_bp);
    }
    let mut accepted_ticks: Vec<usize> = vec![];
    let mut last_y: Option<f32> = None;
    for bp in tick_bps.iter().rev() {
        let y = layout.y_for_bp(*bp, GEL_TOP, GEL_BOTTOM);
        if last_y.map(|v| (v - y).abs() >= 16.0).unwrap_or(true) {
            accepted_ticks.push(*bp);
            last_y = Some(y);
        }
        if accepted_ticks.len() >= 20 {
            break;
        }
    }
    for bp in accepted_ticks {
        let y = layout.y_for_bp(bp, GEL_TOP, GEL_BOTTOM);
        doc = doc
            .add(
                Line::new()
                    .set("x1", GEL_LEFT)
                    .set("y1", y)
                    .set("x2", GEL_RIGHT)
                    .set("y2", y)
                    .set("stroke", "#2d3238")
                    .set("stroke-width", 1),
            )
            .add(
                Text::new(format!("{bp} bp"))
                    .set("x", GEL_RIGHT + 12.0)
                    .set("y", y + 4.0)
                    .set("font-family", "monospace")
                    .set("font-size", 12)
                    .set("fill", "#374151"),
            );
    }

    for (lane_idx, lane) in layout.lanes.iter().enumerate() {
        let x = GEL_LEFT + lane_gap * (lane_idx as f32 + 1.0);
        let lane_fill = if lane.is_ladder { "#1a2028" } else { "#1f252e" };
        doc = doc
            .add(
                Rectangle::new()
                    .set("x", x - 34.0)
                    .set("y", GEL_TOP + 10.0)
                    .set("width", 68.0)
                    .set("height", gel_height - 20.0)
                    .set("rx", 6)
                    .set("ry", 6)
                    .set("fill", lane_fill),
            )
            .add(
                Text::new(lane.name.clone())
                    .set("x", x)
                    .set("y", GEL_BOTTOM + 26.0)
                    .set("text-anchor", "middle")
                    .set("font-family", "monospace")
                    .set("font-size", 13)
                    .set("fill", "#0f172a"),
            );

        for band in &lane.bands {
            let y = layout.y_for_bp(band.bp, GEL_TOP + 14.0, GEL_BOTTOM - 14.0);
            let width = if lane.is_ladder {
                32.0 + 18.0 * band.intensity
            } else {
                38.0 + 24.0 * band.intensity
            };
            let height = if lane.is_ladder {
                2.5 + 2.5 * band.intensity
            } else {
                3.0 + 3.0 * band.intensity
            };
            let fill = if lane.is_ladder { "#e5e7eb" } else { "#f59e0b" };
            doc = doc.add(
                Rectangle::new()
                    .set("x", x - width * 0.5)
                    .set("y", y - height * 0.5)
                    .set("width", width)
                    .set("height", height)
                    .set("rx", 2)
                    .set("ry", 2)
                    .set("fill", fill)
                    .set("opacity", (0.42 + 0.58 * band.intensity).clamp(0.35, 1.0)),
            );

            if !lane.is_ladder && band.count > 0 {
                let mut label = format!("{} bp", band.bp);
                if band.count > 1 {
                    label.push_str(&format!(" (x{})", band.count));
                }
                doc = doc.add(
                    Text::new(label)
                        .set("x", x + 44.0)
                        .set("y", y + 4.0)
                        .set("font-family", "monospace")
                        .set("font-size", 11)
                        .set("fill", "#111827"),
                );
            }
        }
    }

    let ladder_caption = if layout.selected_ladders.is_empty() {
        "auto ladder".to_string()
    } else {
        layout.selected_ladders.join(" + ")
    };
    let title = format!(
        "Serial Gel Preview ({} sample lane(s), {} member(s)) | ladders: {}",
        layout.sample_count, layout.pool_member_count, ladder_caption
    );
    doc = doc
        .add(
            Text::new(title)
                .set("x", GEL_LEFT)
                .set("y", 42.0)
                .set("font-family", "monospace")
                .set("font-size", 16)
                .set("fill", "#0f172a"),
        )
        .add(
            Text::new(format!(
                "range {}..{} bp | lanes: {}",
                layout.range_min_bp, layout.range_max_bp, lane_count
            ))
            .set("x", GEL_LEFT)
            .set("y", 62.0)
            .set("font-family", "monospace")
            .set("font-size", 12)
            .set("fill", "#334155"),
        );

    for sample_idx in sample_lane_indices {
        let x = GEL_LEFT + lane_gap * (sample_idx as f32 + 1.0);
        doc = doc.add(
            Line::new()
                .set("x1", x)
                .set("y1", GEL_TOP)
                .set("x2", x)
                .set("y2", GEL_BOTTOM)
                .set("stroke", "#f59e0b")
                .set("stroke-width", 1.5)
                .set("opacity", 0.25),
        );
    }

    doc.to_string()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[cfg(feature = "snapshot-tests")]
    use std::fs;

    #[test]
    fn test_build_pool_gel_layout_auto_ladders() {
        let members = vec![
            ("frag_a".to_string(), 420),
            ("frag_b".to_string(), 950),
            ("frag_c".to_string(), 1210),
            ("frag_d".to_string(), 1210),
        ];
        let layout = build_pool_gel_layout(&members, &[]).unwrap();
        assert!(!layout.selected_ladders.is_empty());
        assert!(layout.lanes.iter().any(|l| l.is_ladder));
        assert!(layout.lanes.iter().any(|l| !l.is_ladder));
        assert_eq!(layout.sample_count, 1);
        assert!(layout.range_max_bp > layout.range_min_bp);
    }

    #[test]
    fn test_export_pool_gel_svg() {
        let members = vec![
            ("frag_a".to_string(), 350),
            ("frag_b".to_string(), 820),
            ("frag_c".to_string(), 1650),
        ];
        let layout = build_pool_gel_layout(&members, &[]).unwrap();
        let svg = export_pool_gel_svg(&layout);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("Serial Gel Preview"));
    }

    #[test]
    #[cfg(feature = "snapshot-tests")]
    fn snapshot_pool_gel_svg() {
        let members = vec![
            ("frag_a".to_string(), 350),
            ("frag_b".to_string(), 820),
            ("frag_c".to_string(), 1650),
        ];
        let ladders = vec![
            "NEB 100bp DNA Ladder".to_string(),
            "NEB 1kb DNA Ladder".to_string(),
        ];
        let layout = build_pool_gel_layout(&members, &ladders).unwrap();
        let svg = export_pool_gel_svg(&layout);
        let expected = include_str!("../tests/snapshots/pool_gel/minimal.svg");
        assert_eq!(svg, expected);
    }

    #[test]
    #[cfg(feature = "snapshot-tests")]
    #[ignore]
    fn write_pool_gel_snapshot() {
        let members = vec![
            ("frag_a".to_string(), 350),
            ("frag_b".to_string(), 820),
            ("frag_c".to_string(), 1650),
        ];
        let ladders = vec![
            "NEB 100bp DNA Ladder".to_string(),
            "NEB 1kb DNA Ladder".to_string(),
        ];
        let layout = build_pool_gel_layout(&members, &ladders).unwrap();
        let svg = export_pool_gel_svg(&layout);
        fs::create_dir_all("tests/snapshots/pool_gel").unwrap();
        fs::write("tests/snapshots/pool_gel/minimal.svg", svg).unwrap();
    }
}
